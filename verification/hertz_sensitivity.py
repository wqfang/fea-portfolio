"""Benchmark 3 modelling checks: penalty stiffness, contact size a/R, body size.

1. Penalty stiffness K (CAX4 n=8, CAX8 n=4). CalculiX surface-to-surface contact
   is a penalty method: pressure = K * overclosure. Too soft and the bodies
   interpenetrate; too stiff and the nodal contact pressure oscillates.
2. Contact size a/R (CAX4 n=8). Hertz theory assumes a << R. The model is scaled
   so that only a/R changes: the core, the element size and the penalty scale
   with a, so a/h and overclosure/approach stay fixed.
3. Body size (CAX4 n=8): far boundary at R/2 and R from the contact.

Usage
    uv run hertz_sensitivity.py     # ~3 minutes; writes runs/hertz/sensitivity_results.csv
"""
from __future__ import annotations

import csv

import hertz_contact as hc
from fea_common import run_ccx

OUT = hc.RUN_DIR / "sensitivity"
PENALTIES = [1e7, 3e7, 1e8, 3e8, 1e9, 1e10]
PENALTY_CASES = [("CAX4", 8), ("CAX8", 4)]
A_OVER_R = [0.025, 0.05, 0.1]
BODY_SIZES = [0.5 * hc.R, hc.R]


def solve(tag: str, etype: str, level: int, k_pen: float, a_target: float, rho_max: float) -> dict:
    hc.A_TARGET, hc.RHO_MAX = a_target, rho_max
    hc.DELTA, hc.RHO_C = a_target**2 / hc.R, 1.25 * a_target
    case = hc.build_case(etype, level)
    deck = OUT / f"{tag}.inp"
    hc.write_deck(deck, case, etype, "ccx", k_pen=k_pen)
    run_ccx(deck)
    m = hc.metrics(case, *hc.ccx_results(deck))
    m.pop("profile")
    return m


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a0, rm0, k0 = hc.A_TARGET, hc.RHO_MAX, hc.K_PEN
    rows = []

    def record(check, etype, level, k, a_target, rho_max, tag):
        m = solve(tag, etype, level, k, a_target, rho_max)
        rows.append({"check": check, "element": etype, "level": level, "k_pen": k,
                     "a_over_R": a_target / hc.R, "rho_max_over_R": rho_max / hc.R, **m})
        print(f"   {etype:4s} n={level}  K={k:<6g} a/R={a_target / hc.R:<5g} L/R={rho_max / hc.R:<4g}"
              f"  F={m['F']:8.1f}  p0 err={m['p0_err']:+.3%}  a err={m['a_err']:+.3%}"
              f"  rms={m['p_rms_err']:.3%}  pmax/p0={m['p_max_ratio']:.2f}", flush=True)

    print("1. Penalty stiffness")
    for etype, lvl in PENALTY_CASES:
        for k in PENALTIES:
            record("penalty", etype, lvl, k, a0, rm0, f"penalty_{etype}_n{lvl}_K{k:.0e}")
    print("2. Contact size a/R (core, element size and penalty scale with a)")
    for ar in A_OVER_R:
        a = ar * hc.R
        record("a_over_R", "CAX4", 8, k0 * a0 / a, a, rm0, f"aR_{ar:g}".replace(".", "p"))
    print("3. Body size")
    for rm in BODY_SIZES:
        record("body_size", "CAX4", 8, k0, a0, rm, f"body_{rm:g}".replace(".", "p"))
    hc.A_TARGET, hc.RHO_MAX = a0, rm0
    hc.DELTA, hc.RHO_C = a0**2 / hc.R, 1.25 * a0

    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(hc.RUN_DIR / "sensitivity_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"\nResults: {hc.RUN_DIR / 'sensitivity_results.csv'}")


if __name__ == "__main__":
    main()
