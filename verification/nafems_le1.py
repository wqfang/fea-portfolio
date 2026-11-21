"""Benchmark 4: NAFEMS LE1, elliptic membrane under outward pressure.

NAFEMS standard benchmark LE1 (The Standard NAFEMS Benchmarks, 1990). Plane
stress, a quarter of an elliptic annulus:

  inner edge  (x / 2000)^2 + (y / 1000)^2 = 1     free
  outer edge  (x / 3250)^2 + (y / 2750)^2 = 1     uniform outward normal pressure 10 MPa
  x = 0 edge (A-B)                                 symmetry, u_x = 0
  y = 0 edge (D-C)                                 symmetry, u_y = 0

Thickness 100 mm, E = 210 GPa, nu = 0.3 (NAFEMS gives the geometry in metres;
here it is in mm to match the rest of the suite). Target: sigma_yy at point D
(2000, 0), the inner end of the y = 0 edge, = 92.7 MPa.

Mesh: a structured mapping between the two ellipses. A point at parametric
angle t lies at (A cos t, B sin t) on each ellipse; mesh lines join equal t on
the inner and outer ellipse, and are divided evenly between them. Level n has
2n elements across the membrane and 3n round it, so levels 1 and 2 have 6 and
24 elements, the element counts of the NAFEMS coarse and fine meshes (the
NAFEMS meshes place their interior nodes differently).

Sweep: CPS4, CPS4R, CPS8, CPS8R in CalculiX at levels 1-64, each at the
specified 100 mm section and at a 0.01 mm section. CalculiX expands plane-stress
elements into 3D bricks of the section thickness, so the thin section shows how
much the 100 mm brick moves its answer (in true plane stress, sigma_yy(D) does
not depend on the thickness). Abaqus decks (true plane stress, T = 100 mm) are
written for every case under the Learning Edition node cap, including CPS4I,
which CalculiX 2.23 does not support.

The stress at D is the averaged nodal sigma_yy, as a postprocessor would report
it. Errors are given against the NAFEMS target (92.7) and against a Richardson
extrapolation of the three finest thin-section CPS8 runs.

Usage
    uv run nafems_le1.py                      # CalculiX sweep + plot (~2 minutes)
    uv run nafems_le1.py --levels 1 2 4       # quicker
    uv run nafems_le1.py --collect-abaqus     # also read Abaqus .dat results
    uv run nafems_le1.py --check-loads        # mesh sizes and pressure-resultant check only
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from fea_common import parse_dat_nodal_stress_2d, parse_frd_nodal_stress, run_ccx, write_deck_2d
from plate_hole import Registry, block, edges_along

# Geometry (mm), material (MPa), load (MPa)
A_IN, B_IN = 2000.0, 1000.0       # inner ellipse semi-axes (x, y)
A_OUT, B_OUT = 3250.0, 2750.0     # outer ellipse semi-axes
T = 100.0                         # thickness
E, NU = 210_000.0, 0.3
P = 10.0                          # outward normal pressure on the outer edge
TARGET_SYY_D = 92.7               # NAFEMS target, sigma_yy at D (2000, 0)
ELEMENT_ORDER_2D = {"CPS4": 1, "CPS4R": 1, "CPS4I": 1, "CPS8": 2, "CPS8R": 2}
CCX_ELEMENTS = ["CPS4", "CPS4R", "CPS8", "CPS8R"]   # ccx 2.23 rejects CPS4I
LEVELS = [1, 2, 4, 8, 16, 32, 64]
THIN_T = 0.01                     # near-plane-stress section for the CalculiX check
ABAQUS_LE_NODE_CAP = 1000

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE / "runs" / "nafems_le1"

# Named points (NAFEMS labels)
POINTS = {"A": (0.0, B_IN), "B": (0.0, B_OUT), "C": (A_OUT, 0.0), "D": (A_IN, 0.0)}


def ellipse_point(a: float, b: float, t: float) -> tuple[float, float]:
    return a * math.cos(t), b * math.sin(t)


def outward_normal(a: float, b: float, x: float, y: float) -> np.ndarray:
    n = np.array([x / a**2, y / b**2])
    return n / np.linalg.norm(n)


def build_mesh(level: int, order: int):
    """Quarter annulus: i across (inner -> outer), j round (t = 0 at D-C -> 90 deg at A-B)."""
    n_r, n_t = 2 * level, 3 * level
    reg = Registry()

    def point(s, u):
        t = 0.5 * math.pi * u
        xi, yi = ellipse_point(A_IN, B_IN, t)
        xo, yo = ellipse_point(A_OUT, B_OUT, t)
        return xi + s * (xo - xi), yi + s * (yo - yi)

    elems, grid = block(reg, point, n_r, n_t, order)
    st = 2 if order == 2 else 1
    outer = edges_along(grid, st * n_r, None, n_t, order)
    inner = edges_along(grid, 0, None, n_t, order)
    return reg.xy, elems, outer, inner


def pressure_loads(xy, edges, order: int, thickness: float = T) -> dict[tuple[int, int], float]:
    """Consistent nodal forces for the outward pressure P on the outer edge (6-point Gauss)."""
    gp, gw = np.polynomial.legendre.leggauss(6)
    f: dict[tuple[int, int], float] = {}
    for e in edges:
        X = np.array([xy[n] for n in e])
        for xi, w in zip(gp, gw):
            if order == 1:
                N, dN = np.array([(1 - xi) / 2, (1 + xi) / 2]), np.array([-0.5, 0.5])
            else:
                N = np.array([xi * (xi - 1) / 2, xi * (xi + 1) / 2, 1 - xi * xi])
                dN = np.array([xi - 0.5, xi + 0.5, -2 * xi])
            tan = dN @ X
            ds = np.linalg.norm(tan)
            # traction along the outward normal of the element edge itself (as a solver's
            # pressure load does). The edge runs anticlockwise (t increasing), so the
            # outward normal is the tangent turned clockwise.
            t = P * np.array([tan[1], -tan[0]]) / ds
            for k, n in enumerate(e):
                f[(n, 1)] = f.get((n, 1), 0.0) + N[k] * t[0] * ds * w * thickness
                f[(n, 2)] = f.get((n, 2), 0.0) + N[k] * t[1] * ds * w * thickness
    return f


def build_case(etype: str, level: int, thickness: float = T) -> dict:
    order = ELEMENT_ORDER_2D[etype]
    xy, elems, outer, inner = build_mesh(level, order)
    loads = pressure_loads(xy, outer, order, thickness)
    cloads = [(n, d, v) for (n, d), v in sorted(loads.items()) if abs(v) > 1e-12]
    bcs = {"XSYM": (sorted(n for n, (x, _) in xy.items() if abs(x) < 1e-9), 1),
           "YSYM": (sorted(n for n, (_, y) in xy.items() if abs(y) < 1e-9), 2)}
    node_d = next(n for n, (x, y) in xy.items() if abs(x - A_IN) < 1e-6 and abs(y) < 1e-9)
    return {"xy": xy, "elems": elems, "outer": outer, "inner": inner, "bcs": bcs,
            "cloads": cloads, "node_D": node_d, "order": order, "thickness": thickness,
            "h_mm": (A_OUT - A_IN) / (2 * level)}   # element size across the membrane at D


def check_loads() -> None:
    # Resultant of the outward pressure on the quarter outer edge: P * T * (projected length)
    fx_exact, fy_exact = P * T * B_OUT, P * T * A_OUT
    print(f"Pressure resultant (exact): Fx = {fx_exact:.6g} N, Fy = {fy_exact:.6g} N")
    for etype in ("CPS4", "CPS8"):
        for lvl in LEVELS:
            c = build_case(etype, lvl)
            fx = sum(v for _, d, v in c["cloads"] if d == 1)
            fy = sum(v for _, d, v in c["cloads"] if d == 2)
            print(f"{etype} n={lvl:2d}: {len(c['elems']):5d} elements, {len(c['xy']):6d} nodes, "
                  f"load error Fx {fx / fx_exact - 1:+.1e}, Fy {fy / fy_exact - 1:+.1e}")


def tag(t: float) -> str:
    return f"T{t:g}".replace(".", "p")


def richardson(h, d):
    (h1, h2, h3), (d1, d2, d3) = h, d
    r = h1 / h2
    p = math.log(abs((d1 - d2) / (d2 - d3))) / math.log(r)
    return d3 + (d3 - d2) / (r**p - 1.0), p


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", nargs="+", default=list(ELEMENT_ORDER_2D))
    ap.add_argument("--levels", nargs="+", type=int, default=LEVELS)
    ap.add_argument("--ccx", default="ccx", help="CalculiX executable")
    ap.add_argument("--collect-abaqus", action="store_true",
                    help="read Abaqus results from runs/nafems_le1/abaqus/*.dat")
    ap.add_argument("--check-loads", action="store_true", help="only print mesh sizes and the load check")
    args = ap.parse_args()
    if args.check_loads:
        check_loads()
        return

    ccx_dir, abq_dir = RUN_DIR / "ccx", RUN_DIR / "abaqus"
    ccx_dir.mkdir(parents=True, exist_ok=True)
    abq_dir.mkdir(parents=True, exist_ok=True)
    rows, abq_jobs = [], []
    for etype in args.elements:
        for lvl in args.levels:
            for thick in (T, THIN_T):
                c = build_case(etype, lvl, thick)
                nn, d = len(c["xy"]), c["node_D"]
                base = {"element": etype, "level": lvl, "elements": len(c["elems"]), "nodes": nn,
                        "h_mm": c["h_mm"], "thickness_mm": thick}
                job = f"le1_{etype}_n{lvl}"
                head = f"NAFEMS LE1 elliptic membrane, {etype}, level {lvl}, T = {thick:g} mm, {nn} nodes"
                # Abaqus: true plane stress, so only the specified section is needed
                if thick == T and nn <= ABAQUS_LE_NODE_CAP:
                    write_deck_2d(abq_dir / f"{job}.inp", c["xy"], c["elems"], etype, E, NU, thick, c["bcs"],
                                  c["cloads"], head, "abaqus")
                    abq_jobs.append(job)
                    dat = abq_dir / f"{job}.dat"
                    if args.collect_abaqus and dat.exists():
                        syy = float(parse_dat_nodal_stress_2d(dat)[d][1])
                        rows.append({**base, "solver": "abaqus", "syy_D": syy})
                        print(f"abaqus {etype:5s} n={lvl:<2d} nodes={nn:5d}  syy(D) = {syy:8.3f}")
                if etype not in CCX_ELEMENTS:
                    continue
                deck = ccx_dir / f"{job}_{tag(thick)}.inp"
                write_deck_2d(deck, c["xy"], c["elems"], etype, E, NU, thick, c["bcs"], c["cloads"], head, "ccx")
                try:
                    run_ccx(deck, args.ccx)
                    syy = float(parse_frd_nodal_stress(deck.with_suffix(".frd"))[d][1])
                except (RuntimeError, ValueError, KeyError) as exc:
                    print(f"ccx    {etype:5s} n={lvl:<2d} T={thick:<6g} FAILED: {exc}", flush=True)
                    continue
                rows.append({**base, "solver": "ccx", "syy_D": syy})
                print(f"ccx    {etype:5s} n={lvl:<2d} nodes={nn:5d} T={thick:<6g} syy(D) = {syy:8.3f}"
                      f"  ({syy / TARGET_SYY_D - 1:+.3%} vs 92.7)", flush=True)

    bat = ["@echo off", "rem Run every NAFEMS LE1 case that fits the Learning Edition node cap"]
    bat += [f"call abaqus job={j} interactive" for j in abq_jobs]
    (abq_dir / "run_all.bat").write_text("\r\n".join(bat) + "\r\n")
    print(f"\nAbaqus decks under {ABAQUS_LE_NODE_CAP} nodes: {len(abq_jobs)} -> {abq_dir}")

    # Own converged reference: Richardson on the three finest thin-section CPS8 runs (plane stress)
    q = sorted((r for r in rows if r["solver"] == "ccx" and r["element"] == "CPS8" and r["thickness_mm"] == THIN_T),
               key=lambda r: r["h_mm"], reverse=True)
    ref = float("nan")
    if len(q) >= 3:
        ref, p_obs = richardson([r["h_mm"] for r in q[-3:]], [r["syy_D"] for r in q[-3:]])
        print(f"Richardson (CPS8, T = {THIN_T:g}, levels {[r['level'] for r in q[-3:]]}): syy(D) = {ref:.3f} MPa, "
              f"observed order {p_obs:.2f}; NAFEMS target {TARGET_SYY_D}")
    for r in rows:
        r["err_vs_target"] = r["syy_D"] / TARGET_SYY_D - 1
        r["err_vs_extrapolated"] = r["syy_D"] / ref - 1 if ref == ref else ""
    if rows:
        with open(RUN_DIR / "results.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"Results: {RUN_DIR / 'results.csv'}")
        make_plots(rows, ref)


def make_plots(rows, ref) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(ELEMENT_ORDER_2D, ["#c0392b", "#e67e22", "#8e44ad", "#2471a3", "#17a589"]))
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    series = (("ccx", T, "o", "-", "CalculiX, T = 100"), ("ccx", THIN_T, "s", ":", f"CalculiX, T = {THIN_T:g}"),
              ("abaqus", T, "x", "--", "Abaqus"))
    for ax, key, ylab in ((axes[0], "syy_D", "sigma_yy at D (MPa)"),
                          (axes[1], "err", "|error| vs extrapolated value")):
        for et in ELEMENT_ORDER_2D:
            for solver, thick, mk, ls, _ in series:
                rs = sorted((r for r in rows if r["element"] == et and r["solver"] == solver
                             and r["thickness_mm"] == thick), key=lambda r: r["h_mm"])
                if not rs:
                    continue
                x = [r["elements"] for r in rs]
                y = [r["syy_D"] for r in rs] if key == "syy_D" else [abs(r["err_vs_extrapolated"]) or 1e-6 for r in rs]
                lab = et if (solver, thick) == ("ccx", T) or (et == "CPS4I" and solver == "abaqus") else None
                (ax.semilogx if key == "syy_D" else ax.loglog)(x, y, marker=mk, ls=ls, color=colors[et], label=lab,
                                                                ms=5, mfc="none" if mk == "s" else None)
        ax.set_xlabel("number of elements (quarter model)")
        ax.set_ylabel(ylab)
    axes[0].axhline(TARGET_SYY_D, color="k", lw=1, label="NAFEMS 92.7")
    if ref == ref:
        axes[0].axhline(ref, color="k", lw=1, ls=":", label=f"extrapolated {ref:.2f}")
    axes[0].set_ylim(60, 100)
    axes[0].legend(fontsize=8, frameon=False, ncol=2)
    axes[0].set_title("sigma_yy at D")
    axes[1].set_title(f"error vs extrapolated {ref:.3f} MPa\n(o CalculiX T=100, squares T=0.01, x Abaqus)")
    axes[0].text(0.02, 0.03, "CPS4R at T = 0.01 is off scale (9-22 MPa)", transform=axes[0].transAxes,
                 fontsize=8, color=colors["CPS4R"])
    fig.tight_layout()
    fig.savefig(RUN_DIR / "le1_convergence.png", dpi=130)
    print(f"Plot: {RUN_DIR / 'le1_convergence.png'}")


if __name__ == "__main__":
    main()
