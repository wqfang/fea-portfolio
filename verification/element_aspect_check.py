"""Does adding layers through the depth cure C3D8 shear locking, or is it the element shape?

In the main cantilever sweep, refining C3D8 adds layers through the depth AND
shortens the elements along the beam, so the elements stay cubes. This script
separates the two effects with four CalculiX C3D8 mesh series (n = 1, 2, 4, 8):

  S1 cubes        nx = 20n, ny = nz = n   element 10/n x 10/n x 10/n  (main sweep)
  S2 layers only  nx = 20,  ny = nz = n   element 10 x 10/n x 10/n    (length/depth = n)
  S3 length only  nx = 20n, ny = nz = 1   element 10/n x 10 x 10      (length/depth = 1/n)
  S4 as S3 with nu = 0                                                 (Poisson effect removed)

S1-S3 use nu = 0.3 and are compared with the converged 3D value (Richardson
extrapolation of the three finest C3D20 meshes). S4 is compared with Timoshenko
at nu = 0: with no Poisson effect the clamped root blocks nothing, so the 3D
answer and beam theory coincide. S3 and S4 also run n = 16.

The compared quantity is the load-weighted tip deflection.

Usage
    uv run element_aspect_check.py        # ~30 s; writes runs/cantilever/aspect_check_results.csv
"""
from __future__ import annotations

import csv

import cantilever as cant
from fea_common import (ELEMENT_ORDER, box_hex_mesh, face_traction_loads,
                        parse_dat_displacements, run_ccx, write_deck)
from c3d8r_hourglass_check import converged_reference

OUT = cant.RUN_DIR / "aspect_check"
ETYPE = "C3D8"
NX0 = int(round(cant.L / cant.H))  # 20 elements along the beam for cube-shaped elements at n = 1
# name: (description, mesh divisions for level n, Poisson's ratio, levels)
SERIES = {
    "S1": ("cubes", lambda n: (NX0 * n, n, n), 0.3, [1, 2, 4, 8]),
    "S2": ("layers only", lambda n: (NX0, n, n), 0.3, [1, 2, 4, 8]),
    "S3": ("length only", lambda n: (NX0 * n, 1, 1), 0.3, [1, 2, 4, 8, 16]),
    "S4": ("length, nu=0", lambda n: (NX0 * n, 1, 1), 0.0, [1, 2, 4, 8, 16]),
}


def timoshenko(nu: float) -> float:
    I, A, G = cant.B * cant.H**3 / 12, cant.B * cant.H, cant.E / (2 * (1 + nu))
    kappa = 10 * (1 + nu) / (12 + 11 * nu)  # Cowper, rectangle
    return cant.P * cant.L**3 / (3 * cant.E * I) + cant.P * cant.L / (kappa * G * A)


def run(series: str, level: int) -> dict:
    _, divisions, nu, _ = SERIES[series]
    nx, ny, nz = divisions(level)
    mesh = box_hex_mesh(cant.L, cant.H, cant.B, nx, ny, nz, ELEMENT_ORDER[ETYPE])
    fixed = mesh.node_set(lambda x, y, z: x < 1e-9)
    loads = face_traction_loads(mesh, i_face=nx, total_force=cant.P)
    deck = OUT / f"{series}_n{level}.inp"
    write_deck(deck, mesh, ETYPE, cant.E, nu, fixed, {n: (2, f) for n, f in loads.items()},
               sorted(loads), f"Element aspect check, {series} ({SERIES[series][0]}), n={level}")
    d = cant.tip_deflection(parse_dat_displacements(run_ccx(deck)), loads)
    le, h = cant.L / nx, cant.H / ny
    return {"series": series, "description": SERIES[series][0], "nu": nu, "level": level, "nx": nx, "ny": ny, "nz": nz,
            "nodes": mesh.n_nodes, "elem_length_mm": le, "elem_depth_mm": h, "length_over_depth": le / h,
            "tip_mm": d}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ref = converged_reference()
    if ref is None:
        raise SystemExit("Run cantilever.py first: the converged reference comes from its results.csv")
    rows = []
    print(f"{'series':15s} {'n':>2s} {'nodes':>6s} {'L_e/h':>6s} {'tip (mm)':>10s} {'x reference':>12s}")
    for s, (_, _, nu, levels) in SERIES.items():
        reference = ref if nu == cant.NU else timoshenko(nu)
        for lvl in levels:
            r = run(s, lvl)
            r["reference_mm"] = reference
            r["ratio_to_reference"] = r["tip_mm"] / reference
            rows.append(r)
            print(f"{s} {r['description']:12s} {lvl:2d} {r['nodes']:6d} {r['length_over_depth']:6.3g} "
                  f"{r['tip_mm']:10.6f} {r['ratio_to_reference']:12.4f}", flush=True)
    with open(OUT.parent / "aspect_check_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nReferences: converged 3D {ref:.6f} mm (nu = 0.3); Timoshenko {timoshenko(0.0):.6f} mm (nu = 0)."
          f" Results: {OUT.parent / 'aspect_check_results.csv'}")


if __name__ == "__main__":
    main()
