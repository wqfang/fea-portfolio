"""Single-element tests of the Neo-Hookean UMAT in CalculiX, against closed form.

One C3D8 unit cube, NLGEOM, four load cases:
  uniaxial    stretch to lambda = 2 on symmetry planes, lateral faces free
  biaxial     equibiaxial stretch to lambda = 1.5, thickness free
  shear       simple shear to gamma = 1 (every node prescribed): includes rotation
  rotation    stretch 1.5 along x, then the whole cube rotated 90 deg about z
              (every node prescribed): checks the stress turns with the body

Each case runs twice in CalculiX: with the UMAT (material ABAQUSNL...) and with the
built-in *HYPERELASTIC, NEO HOOKE. Both are compared with the closed-form Cauchy
stress at the final deformation gradient F, which is computed from the nodal
displacements (homogeneous deformation). For uniaxial and biaxial the free lateral
stretch is also solved independently from the zero-stress condition.

Abaqus Learning Edition decks with the built-in model are written to
runs/abaqus/ (it cannot run user subroutines).

Usage (from umat/):  uv run tests/single_element.py --ccx /path/to/ccx_umat
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "ref"))
from neohooke import cauchy  # noqa: E402

RUN = ROOT / "runs" / "single_element"
C10, D1 = 0.5, 0.02  # MPa, 1/MPa: mu0 = 1 MPa, K0 = 100 MPa
XYZ = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def rot_z(deg: float) -> np.ndarray:
    t = np.radians(deg)
    return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])


def all_nodes_bc(F: np.ndarray) -> list[str]:
    """Prescribe u = (F - I) X at every node."""
    u = XYZ @ (F - np.eye(3)).T
    return [f"{n + 1}, {d + 1}, {d + 1}, {u[n, d]:.12g}" for n in range(8) for d in range(3)]


CASES = {
    "uniaxial": ["1, 1, 3, 0.", "4, 1, 1, 0.", "5, 1, 1, 0.", "8, 1, 1, 0.",
                 "2, 2, 3, 0.", "6, 2, 2, 0.", "3, 3, 3, 0.", "4, 3, 3, 0.",
                 "2, 1, 1, 1.", "3, 1, 1, 1.", "6, 1, 1, 1.", "7, 1, 1, 1."],
    "biaxial": ["1, 1, 3, 0.", "4, 1, 1, 0.", "5, 1, 1, 0.", "8, 1, 1, 0.",
                "2, 2, 2, 0.", "5, 2, 2, 0.", "6, 2, 2, 0.", "2, 3, 3, 0.",
                "3, 3, 3, 0.", "4, 3, 3, 0.",
                "2, 1, 1, 0.5", "3, 1, 1, 0.5", "6, 1, 1, 0.5", "7, 1, 1, 0.5",
                "3, 2, 2, 0.5", "4, 2, 2, 0.5", "7, 2, 2, 0.5", "8, 2, 2, 0.5"],
    "shear": all_nodes_bc(np.array([[1, 1.0, 0], [0, 1, 0], [0, 0, 1]])),
    "rotation": all_nodes_bc(rot_z(90) @ np.diag([1.5, 1, 1])),
}


def deck(case: str, material: str) -> str:
    if material == "umat":
        mat = f"*MATERIAL, NAME=ABAQUSNL_NEOHOOKE\n*USER MATERIAL, CONSTANTS=2\n{C10}, {D1}"
        name = "ABAQUSNL_NEOHOOKE"
    else:
        mat = f"*MATERIAL, NAME=NEOHOOKE\n*HYPERELASTIC, NEO HOOKE\n{C10}, {D1}"
        name = "NEOHOOKE"
    nodes = "\n".join(f"{i + 1}, {x:g}, {y:g}, {z:g}" for i, (x, y, z) in enumerate(XYZ))
    return f"""*HEADING
Single element {case}, {material} Neo-Hookean, C10={C10} D1={D1}
*NODE, NSET=NALL
{nodes}
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
{mat}
*SOLID SECTION, ELSET=EALL, MATERIAL={name}
*STEP, NLGEOM, INC=200
*STATIC
0.05, 1., 1e-6, 0.05
*BOUNDARY
{chr(10).join(CASES[case])}
*NODE PRINT, NSET=NALL
U
*EL PRINT, ELSET=EALL
S
*END STEP
"""


_NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"


def last_block(text: str, header: str, ncols: int) -> list[list[float]]:
    blk = text[text.rindex(header):]
    rows = []
    for line in blk.splitlines()[1:]:
        vals = re.findall(_NUM, line)
        if len(vals) == ncols:
            rows.append([float(v) for v in vals])
        elif rows and line.strip() and not line.strip()[0].isdigit():
            break
    return rows


def run_ccx(exe: str, job: Path) -> dict:
    env = {**os.environ, "OMP_NUM_THREADS": "1"}  # multithreaded ccx 2.23 gave silent wrong answers
    res = subprocess.run([exe, "-i", job.stem], cwd=job.parent, capture_output=True, text=True, env=env)
    (job.parent / f"{job.stem}.ccx.log").write_text(res.stdout + res.stderr)
    dat = job.with_suffix(".dat").read_text()
    if res.returncode != 0 or "*ERROR" in res.stdout or "stresses" not in dat:
        raise RuntimeError(f"ccx failed on {job.name}; see {job.stem}.ccx.log")
    u = {int(r[0]): np.array(r[1:]) for r in last_block(dat, "displacements", 4)}
    s = np.array([r[2:] for r in last_block(dat, "stresses", 8)])  # elem, ip, 6 comps
    sta = job.with_suffix(".sta").read_text().splitlines()
    incs = [ln.split() for ln in sta if re.match(r"^\s*1\s+\d+\s+\d+\s+\d+", ln)]
    iters = [int(r[3]) for r in incs]
    return {"u": u, "s": s, "incs": len(iters), "iters": sum(iters)}


def read_abaqus_dat(path: Path) -> dict:
    """Final nodal displacements and integration-point stresses from an Abaqus .dat file."""
    text = path.read_text(errors="replace")
    e_mark, n_mark = "E L E M E N T   O U T P U T", "N O D E   O U T P U T"
    if e_mark not in text or n_mark not in text:
        raise ValueError(f"no element/node output in {path}")
    def rows(block: str, n: int) -> list[list[float]]:
        out = []
        for line in block.splitlines():
            vals = re.findall(_NUM, line)
            if len(vals) == n and re.match(r"^\s*\d+\s", line):
                out.append([float(v.replace("D", "E")) for v in vals])
        return out
    i_e, i_n = text.rindex(e_mark), text.rindex(n_mark)
    s = np.array([r[2:] for r in rows(text[i_e:i_n], 8)])
    u = {int(r[0]): np.array(r[1:]) for r in rows(text[i_n:], 4)}
    return {"u": u, "s": s}


def F_from_u(u: dict) -> np.ndarray:
    """Homogeneous F from nodal displacements: x = F X, least squares over 8 nodes."""
    X = XYZ
    x = X + np.array([u[n + 1] for n in range(8)])
    A = np.hstack([X, np.ones((8, 1))])
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    return coef[:3].T


def voigt(T: np.ndarray) -> np.ndarray:
    return np.array([T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[0, 2], T[1, 2]])


def lateral_stretch(case: str, lam: float) -> float:
    """Independent closed form: free lateral stretch from zero lateral stress."""
    if case == "uniaxial":
        f = lambda l2: cauchy(np.diag([lam, l2, l2]), C10, D1)[1, 1]
    else:
        f = lambda l3: cauchy(np.diag([lam, lam, l3]), C10, D1)[2, 2]
    return brentq(f, 0.2, 1.2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccx", default="ccx_umat", help="CalculiX built with the UMAT")
    ap.add_argument("--collect-abaqus", action="store_true",
                    help="also compare Abaqus results in runs/single_element/abaqus/*.dat")
    args = ap.parse_args()
    exe = str(Path(shutil.which(args.ccx) or args.ccx).resolve())
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / "abaqus").mkdir(exist_ok=True)

    rows = []
    for case in CASES:
        for material in ("umat", "builtin"):
            job = RUN / f"{case}_{material}.inp"
            job.write_text(deck(case, material))
            if material == "builtin":  # Abaqus spells the geometric-nonlinearity flag NLGEOM=YES
                (RUN / "abaqus" / f"{case}_builtin.inp").write_text(
                    deck(case, material).replace("*STEP, NLGEOM,", "*STEP, NLGEOM=YES,"))
            r = run_ccx(exe, job)
            F = F_from_u(r["u"])
            exact = voigt(cauchy(F, C10, D1))
            s_mean = r["s"].mean(axis=0)
            scale = np.abs(exact).max()
            err = np.abs(r["s"] - exact).max() / scale  # worst integration point
            row = {"case": case, "material": material, "increments": r["incs"],
                   "iterations": r["iters"], "max_rel_stress_error": err,
                   **{f"s{c}": v for c, v in zip(("11", "22", "33", "12", "13", "23"), s_mean)}}
            if case in ("uniaxial", "biaxial"):
                lam = 2.0 if case == "uniaxial" else 1.5
                free = F[1, 1] if case == "uniaxial" else F[2, 2]
                row["free_stretch_fe"] = free
                row["free_stretch_exact"] = lateral_stretch(case, lam)
            rows.append(row)
            extra = (f"  free stretch {row['free_stretch_fe']:.6f} vs exact "
                     f"{row['free_stretch_exact']:.6f}") if "free_stretch_fe" in row else ""
            print(f"{case:9s} {material:8s} incs={r['incs']:2d} iters={r['iters']:3d} "
                  f"max rel stress error={err:.1e}  s11={s_mean[0]:9.5f} s12={s_mean[3]:9.5f}{extra}")

    if args.collect_abaqus:
        for case in CASES:
            dat = RUN / "abaqus" / f"{case}_builtin.dat"
            if not dat.exists():
                print(f"(no Abaqus result for {case})")
                continue
            r = read_abaqus_dat(dat)
            F = F_from_u(r["u"])
            exact = voigt(cauchy(F, C10, D1))
            err = np.abs(r["s"] - exact).max() / np.abs(exact).max()
            umat = next(x for x in rows if x["case"] == case and x["material"] == "umat")
            s_mean = r["s"].mean(axis=0)
            row = {"case": case, "material": "abaqus_builtin", "max_rel_stress_error": err,
                   **{f"s{c}": v for c, v in zip(("11", "22", "33", "12", "13", "23"), s_mean)}}
            rows.append(row)
            d = np.abs(s_mean - np.array([umat[f"s{c}"] for c in ("11", "22", "33", "12", "13", "23")]))
            print(f"{case:9s} abaqus   max rel stress error={err:.1e}  "
                  f"|ccx UMAT - Abaqus built-in| max = {d.max():.2e} MPa")

    with open(RUN / "results.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    bat = ["@echo off"] + [f"call abaqus job={c}_builtin interactive" for c in CASES]
    (RUN / "abaqus" / "run_all.bat").write_text("\r\n".join(bat) + "\r\n")
    print(f"\nResults: {RUN / 'results.csv'}\nAbaqus decks (built-in model): {RUN / 'abaqus'}")


if __name__ == "__main__":
    main()
