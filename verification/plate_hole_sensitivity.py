"""Benchmark 2 modelling checks: section thickness in CalculiX, and plate length.

1. Thickness (Kirsch problem, CalculiX). CalculiX expands each plane-stress element
   into a 3D element of the section thickness, so it really solves a thin 3D plate.
   In true plane stress the thickness only scales the loads; here it changes the
   answer. Scanned for CPS8 (3D effect near the hole), CPS4, and CPS4R (CalculiX's
   expanded reduced-integration brick fails when it is much thinner than wide).
2. Plate length (finite plate, CPS8 n = 4): is L/W = 3 long enough to stand in for
   the infinitely long plate the handbook Kt assumes?

Usage
    uv run plate_hole_sensitivity.py     # ~1 minute; writes runs/plate_hole/sensitivity_results.csv
"""
from __future__ import annotations

import csv

import plate_hole as ph
from fea_common import parse_frd_nodal_stress, run_ccx, write_deck_2d

OUT = ph.RUN_DIR / "sensitivity"
THICKNESSES = [0.01, 0.1, 1.0, 5.0]
THICKNESS_CASES = [("CPS8", [2, 4, 8]), ("CPS4", [2, 4, 8]), ("CPS4R", [2, 4, 8])]
HALF_LENGTHS = [100.0, 150.0, 300.0, 600.0]


def peak_and_szz(problem: str, etype: str, level: int, tag: str) -> tuple[float, float]:
    xy, elems, bcs, cloads, _, peak = ph.build_case(problem, etype, level)
    deck = OUT / f"{tag}.inp"
    write_deck_2d(deck, xy, elems, etype, ph.E, ph.NU, ph.T, bcs, cloads, tag, "ccx")
    run_ccx(deck)
    s = parse_frd_nodal_stress(deck.with_suffix(".frd"))[peak]
    return s[0], s[2]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0, lh0 = ph.T, ph.LH
    rows = []
    print("1. Section thickness, Kirsch problem (exact peak 300 MPa)")
    for etype, levels in THICKNESS_CASES:
        for t in THICKNESSES:
            ph.T = t
            for lvl in levels:
                sxx, szz = peak_and_szz("kirsch", etype, lvl, f"kirsch_{etype}_n{lvl}_t{t:g}".replace(".", "p"))
                rows.append({"check": "thickness", "problem": "kirsch", "element": etype, "level": lvl,
                             "thickness_mm": t, "length_over_width": "", "sigma_peak": sxx, "szz_peak": szz,
                             "peak_rel_err": (sxx - 3 * ph.S) / (3 * ph.S)})
                print(f"   {etype:5s} T={t:<5g} n={lvl}  peak={sxx:8.3f}  szz={szz:8.4f}", flush=True)
    ph.T = t0

    print("2. Plate length, finite plate, CPS8 n=4")
    for lh in HALF_LENGTHS:
        ph.LH = lh
        sxx, szz = peak_and_szz("plate", "CPS8", 4, f"plate_CPS8_n4_L{int(2 * lh)}")
        rows.append({"check": "length", "problem": "plate", "element": "CPS8", "level": 4, "thickness_mm": ph.T,
                     "length_over_width": 2 * lh / (2 * ph.C), "sigma_peak": sxx, "szz_peak": szz,
                     "peak_rel_err": ""})
        print(f"   L/W={2 * lh / (2 * ph.C):4.1f}  peak={sxx:8.3f}  Kt_net={sxx * (1 - ph.A / ph.C) / ph.S:.4f}")
    ph.LH = lh0

    with open(ph.RUN_DIR / "sensitivity_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nResults: {ph.RUN_DIR / 'sensitivity_results.csv'}")


if __name__ == "__main__":
    main()
