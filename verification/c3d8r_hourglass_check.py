"""C3D8R in bending: how much the answer depends on Abaqus hourglass control.

With one element through the depth, bending of a reduced-integration brick is
an hourglass mode, so its only bending stiffness is the artificial hourglass
stiffness. This script writes the cantilever C3D8R decks (n = 1, 2: the meshes
under the Abaqus Learning Edition cap) in four hourglass-control variants:

  default     no *SECTION CONTROLS (what cantilever.py runs)
  stiff1      HOURGLASS=STIFFNESS, scale factor 1
  stiff10     HOURGLASS=STIFFNESS, scale factor 10
  enhanced    HOURGLASS=ENHANCED

CalculiX 2.23 ignores *SECTION CONTROLS (it warns "parameter not recognized" and
gives the same answer), so this is an Abaqus-only study.

Usage
    uv run c3d8r_hourglass_check.py      # write decks + run_hg.bat; tabulate any .dat present

Copy runs/cantilever/hourglass/abaqus/ to the Windows VM, run run_hg.bat, copy
the .dat files back, and run this script again.
"""
from __future__ import annotations

import csv

import cantilever as cant
from fea_common import parse_dat_displacements, write_deck

OUT = cant.RUN_DIR / "hourglass" / "abaqus"
LEVELS = [1, 2]
SECTION = "*SOLID SECTION, ELSET=EALL, MATERIAL=MAT"
VARIANTS = {
    "default": None,
    "stiff1": "*SECTION CONTROLS, NAME=HG, HOURGLASS=STIFFNESS\n1.",
    "stiff10": "*SECTION CONTROLS, NAME=HG, HOURGLASS=STIFFNESS\n10.",
    "enhanced": "*SECTION CONTROLS, NAME=HG, HOURGLASS=ENHANCED",
}


def converged_reference() -> float | None:
    """Richardson value from the three finest ccx C3D20 meshes in results.csv, if available."""
    path = cant.RUN_DIR / "results.csv"
    if not path.exists():
        return None
    with open(path) as fh:
        q = sorted((r for r in csv.DictReader(fh) if r["solver"] == "ccx" and r["element"] == "C3D20"),
                   key=lambda r: float(r["h_mm"]), reverse=True)
    if len(q) < 3:
        return None
    return cant.richardson([float(r["h_mm"]) for r in q[-3:]], [float(r["tip_mm"]) for r in q[-3:]])[0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    for level in LEVELS:
        mesh, fixed, cloads, tip_nodes, _ = cant.build_case("C3D8R", level)
        for name, controls in VARIANTS.items():
            job = f"c3d8r_n{level}_{name}"
            deck = OUT / f"{job}.inp"
            write_deck(deck, mesh, "C3D8R", cant.E, cant.NU, fixed, cloads, tip_nodes,
                       f"C3D8R hourglass check, n={level}, {name}")
            if controls:
                text = deck.read_text()
                assert text.count(SECTION) == 1
                deck.write_text(text.replace(SECTION, f"{SECTION}, CONTROLS=HG\n{controls}"))
            jobs.append(job)
    bat = ["@echo off", "rem C3D8R hourglass-control variants"]
    bat += [f"call abaqus job={j} interactive" for j in jobs]
    (OUT / "run_hg.bat").write_text("\r\n".join(bat) + "\r\n")
    print(f"Wrote {len(jobs)} decks and run_hg.bat -> {OUT}")

    timo = cant.timoshenko_tip_deflection()["timoshenko"]
    ref = converged_reference()
    rows = []
    for level in LEVELS:
        _, _, _, _, tip_loads = cant.build_case("C3D8R", level)
        for name in VARIANTS:
            dat = OUT / f"c3d8r_n{level}_{name}.dat"
            if dat.exists():
                d = cant.tip_deflection(parse_dat_displacements(dat), tip_loads)
                rows.append({"level": level, "variant": name, "tip_mm": d, "ratio_to_timoshenko": d / timo,
                             "ratio_to_extrapolated": d / ref if ref else ""})
    if not rows:
        print("No Abaqus .dat files yet; run run_hg.bat on the VM and copy the .dat files back.")
        return
    with open(OUT.parent / "hourglass_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\n{'level':>5s} {'variant':>9s} {'tip (mm)':>12s} {'x Timoshenko':>13s} {'x converged':>12s}")
    for r in rows:
        conv = f"{r['ratio_to_extrapolated']:12.4f}" if ref else f"{'-':>12s}"
        print(f"{r['level']:5d} {r['variant']:>9s} {r['tip_mm']:12.6f} {r['ratio_to_timoshenko']:13.4f} {conv}")
    print(f"\nTimoshenko {timo:.6f} mm" + (f"; converged 3D (C3D20 Richardson) {ref:.6f} mm" if ref else ""))


if __name__ == "__main__":
    main()
