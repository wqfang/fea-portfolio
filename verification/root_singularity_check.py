"""Why the C3D20 cantilever converges at order ~1.6 instead of ~4.

For quadratic elements on a smooth problem, the load-weighted deflection (the
compliance, an energy quantity) should converge as h^4. The benchmark converges
much more slowly. Two features of the benchmark make the exact solution
singular:

  root   the fully clamped root blocks the Poisson contraction that the beam
         just beside it wants, so stresses are singular along the root edges
  tip    a uniform shear traction on the tip face meets the traction-free top
         and bottom faces with a jump (tau_xy = P/A on the face, 0 on the
         surface), a weaker singularity along the tip edges

This script runs C3D20 at n = 1, 2, 4, 8 in CalculiX for three variants:

  A clamped + uniform     the benchmark as in cantilever.py
  B clamped + parabolic   tip shear parabolic in y (zero at the top and bottom
                          faces): removes the tip-edge jump only
  C smooth                root face fixed in x only (free to contract), with
                          the root shear applied as an equal and opposite
                          parabolic traction. The loads balance exactly, so the
                          three point supports that stop rigid-body motion
                          carry no load: no singularity left

and reports the observed order from each pair of meshes. The compared quantity
is the compliance sum(f_i * v_i) / P over every loaded node (for A and B this is
the benchmark's load-weighted tip deflection).

Usage
    uv run root_singularity_check.py              # ~1-2 minutes
    uv run root_singularity_check.py --levels 1 2 4
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

import cantilever as cant
from fea_common import (ELEMENT_ORDER, box_hex_mesh, face_traction_loads,
                        parse_dat_displacements, run_ccx, write_deck)

OUT = cant.RUN_DIR / "root_check"
ETYPE = "C3D20"
L, H, B, P = cant.L, cant.H, cant.B, cant.P

# 8-node serendipity face: local node offsets on the (j, k) grid and natural coordinates
FACE_NODES = [((0, 0), (-1, -1)), ((2, 0), (1, -1)), ((2, 2), (1, 1)), ((0, 2), (-1, 1)),
              ((1, 0), (0, -1)), ((2, 1), (1, 0)), ((1, 2), (0, 1)), ((0, 1), (-1, 0))]
GAUSS3 = [(-np.sqrt(0.6), 5 / 9), (0.0, 8 / 9), (np.sqrt(0.6), 5 / 9)]


def serendipity8(eta: float, zeta: float) -> np.ndarray:
    N = []
    for _, (ei, zi) in FACE_NODES:
        if ei and zi:
            N.append(0.25 * (1 + eta * ei) * (1 + zeta * zi) * (eta * ei + zeta * zi - 1))
        elif ei:
            N.append(0.5 * (1 + eta * ei) * (1 - zeta**2))
        else:
            N.append(0.5 * (1 - eta**2) * (1 + zeta * zi))
    return np.array(N)


def parabolic_shear_loads(mesh, i_face: int, total: float) -> dict[int, float]:
    """Consistent nodal y-forces for tau(y) = 6 total y (H - y) / (B H^3) on grid plane i_face.

    The traction integrates to `total` and is zero at y = 0 and y = H. 3x3 Gauss
    is exact (integrand degree <= 5 in each direction).
    """
    nx, ny, nz = mesh.shape
    dy, dz = H / ny, B / nz
    loads: dict[int, float] = {}
    for k in range(nz):
        for j in range(ny):
            for ge, we in GAUSS3:
                for gz, wz in GAUSS3:
                    y = (j + 0.5 * (ge + 1)) * dy
                    tau = 6.0 * total * y * (H - y) / (B * H**3)
                    N = serendipity8(ge, gz)
                    for a, ((dj, dk), _) in enumerate(FACE_NODES):
                        nid = mesh.grid_id[(i_face, 2 * j + dj, 2 * k + dk)]
                        loads[nid] = loads.get(nid, 0.0) + N[a] * tau * we * wz * dy * dz / 4
    return loads


def build(variant: str, level: int):
    nx, ny, nz = int(round(L / H)) * level, level, level
    mesh = box_hex_mesh(L, H, B, nx, ny, nz, ELEMENT_ORDER[ETYPE])
    root = mesh.node_set(lambda x, y, z: x < 1e-9)
    if variant == "A":
        loads = face_traction_loads(mesh, 2 * nx, P)
    else:
        loads = parabolic_shear_loads(mesh, 2 * nx, P)
    if variant == "C":
        for n, f in parabolic_shear_loads(mesh, 0, -P).items():
            loads[n] = loads.get(n, 0.0) + f
        corner, top = mesh.grid_id[(0, 0, 0)], mesh.grid_id[(0, 2 * ny, 0)]
        # root face: u_x = 0 (free to contract); u_y at one corner, u_z at two
        bcs = [f"ROOT, 1, 1, 0.", f"{corner}, 2, 3, 0.", f"{top}, 3, 3, 0."]
    else:
        bcs = ["ROOT, 1, 3, 0."]
    assert abs(sum(f for n, f in loads.items() if mesh.nodes[n][0] > L - 1e-9) - P) < 1e-9 * P
    return mesh, root, loads, bcs


def run(variant: str, level: int) -> float:
    mesh, root, loads, bcs = build(variant, level)
    job = f"{variant}_n{level}"
    deck = OUT / f"{job}.inp"
    loaded = sorted(loads)
    write_deck(deck, mesh, ETYPE, cant.E, cant.NU, root, {n: (2, f) for n, f in loads.items()}, loaded,
               f"Root singularity check, variant {variant}, n={level}")
    text = deck.read_text().replace("*NSET, NSET=FIXED", "*NSET, NSET=ROOT")
    text = text.replace("*BOUNDARY\nFIXED, 1, 3, 0.", "*BOUNDARY\n" + "\n".join(bcs))
    deck.write_text(text)
    disp = parse_dat_displacements(run_ccx(deck))
    return sum(f * disp[n][1] for n, f in loads.items()) / P


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--levels", nargs="+", type=int, default=[1, 2, 4, 8])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    names = {"A": "clamped + uniform", "B": "clamped + parabolic", "C": "smooth"}
    rows = []
    for v in names:
        vals = []
        for lvl in args.levels:
            d = run(v, lvl)
            vals.append(d)
            row = {"variant": v, "description": names[v], "level": lvl, "h_mm": H / lvl, "compliance_mm": d,
                   "change_from_previous": "", "observed_order": ""}
            if len(vals) >= 2:
                row["change_from_previous"] = vals[-1] - vals[-2]
            if len(vals) >= 3:
                r = np.log(abs((vals[-3] - vals[-2]) / (vals[-2] - vals[-1]))) / np.log(2.0)
                row["observed_order"] = r
            rows.append(row)
            print(f"{v} {names[v]:20s} n={lvl}  compliance={d:.9f} mm"
                  + (f"  change={row['change_from_previous']:+.3e}" if row["change_from_previous"] != "" else "")
                  + (f"  order={row['observed_order']:.2f}" if row["observed_order"] != "" else ""), flush=True)
    with open(OUT.parent / "root_check_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nResults: {OUT.parent / 'root_check_results.csv'}")


if __name__ == "__main__":
    main()
