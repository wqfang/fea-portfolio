"""Benchmark 4 modelling check: why CPS4 overshoots sigma_yy at D.

Linear elements have straight edges, so a CPS4 mesh models the elliptic hole as
an inscribed polygon, with a corner at D. This check isolates that geometry
error: it runs CPS8 twice on the same mapped mesh, once with the midside nodes
on the curves (the normal mesh) and once with every midside node moved to the
midpoint of its chord, which gives exactly the CPS4 polygon but keeps the
quadratic field. Both at the thin 0.01 mm section, so CalculiX's 3D expansion
does not enter. The CPS4 rows of runs/nafems_le1/results.csv are printed next
to them for comparison.

Usage
    uv run le1_straight_edge_check.py   # ~1 minute; writes runs/nafems_le1/straight_edge/results.csv
"""
from __future__ import annotations

import csv
import math

import nafems_le1 as le1
from fea_common import parse_frd_nodal_stress, run_ccx, write_deck_2d

OUT = le1.RUN_DIR / "straight_edge"
LEVELS = [1, 2, 4, 8, 16, 32]
REF = 92.658   # Richardson value from nafems_le1.py (CPS8, T = 0.01, levels 16-32-64)


def straighten(xy: dict, elems: list) -> dict:
    """Move each CPS8 midside node to the midpoint of its two corner nodes."""
    xy = dict(xy)
    for e in elems:
        for k in range(4):
            a, b, mid = e[k], e[(k + 1) % 4], e[4 + k]
            xy[mid] = ((xy[a][0] + xy[b][0]) / 2, (xy[a][1] + xy[b][1]) / 2)
    return xy


def solve(level: int, straight: bool) -> float:
    xy, elems, outer, _ = le1.build_mesh(level, 2)
    if straight:
        xy = straighten(xy, elems)
    loads = le1.pressure_loads(xy, outer, 2, le1.THIN_T)
    cloads = [(n, d, v) for (n, d), v in sorted(loads.items()) if abs(v) > 1e-12]
    bcs = {"XSYM": (sorted(n for n, (x, _) in xy.items() if abs(x) < 1e-9), 1),
           "YSYM": (sorted(n for n, (_, y) in xy.items() if abs(y) < 1e-9), 2)}
    node_d = next(n for n, (x, y) in xy.items() if abs(x - le1.A_IN) < 1e-6 and abs(y) < 1e-9)
    name = "straight" if straight else "curved"
    deck = OUT / f"cps8_{name}_n{level}.inp"
    write_deck_2d(deck, xy, elems, "CPS8", le1.E, le1.NU, le1.THIN_T, bcs, cloads,
                  f"NAFEMS LE1, CPS8 {name} edges, level {level}, T = {le1.THIN_T:g} mm", "ccx")
    run_ccx(deck)
    return float(parse_frd_nodal_stress(deck.with_suffix(".frd"))[node_d][1])


def corner_turn_deg(level: int) -> float:
    """Angle the polygon boundary turns through at D (two chords meeting at the symmetry line)."""
    t1 = 0.5 * math.pi / (3 * level)
    x1, y1 = le1.ellipse_point(le1.A_IN, le1.B_IN, t1)
    return 2 * math.degrees(math.atan2(le1.A_IN - x1, y1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cps4 = {}
    sweep = le1.RUN_DIR / "results.csv"
    if sweep.exists():
        for r in csv.DictReader(sweep.open()):
            if r["element"] == "CPS4" and (r["solver"] == "abaqus" or float(r["thickness_mm"]) < 1):
                cps4[(r["solver"], int(r["level"]))] = float(r["syy_D"])

    rows = []
    print(f"{'level':>5} {'turn':>6} {'CPS8 curved':>12} {'CPS8 straight':>14} {'CPS4 ccx':>9} {'CPS4 abq':>9}")
    for lvl in LEVELS:
        r = {"level": lvl, "elements": 6 * lvl * lvl, "corner_turn_deg": corner_turn_deg(lvl),
             "cps8_curved": solve(lvl, False), "cps8_straight": solve(lvl, True),
             "cps4_ccx_thin": cps4.get(("ccx", lvl)), "cps4_abaqus": cps4.get(("abaqus", lvl))}
        rows.append(r)
        err = lambda v: f"{v / REF - 1:+.2%}" if v is not None else "-"
        print(f"{lvl:5d} {r['corner_turn_deg']:5.1f}d {err(r['cps8_curved']):>12} {err(r['cps8_straight']):>14}"
              f" {err(r['cps4_ccx_thin']):>9} {err(r['cps4_abaqus']):>9}", flush=True)

    with (OUT / "results.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Results: {OUT / 'results.csv'}")


if __name__ == "__main__":
    main()
