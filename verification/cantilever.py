"""Benchmark 1: end-loaded cantilever beam, mesh convergence in CalculiX and Abaqus.

A steel beam (L x H x B = 200 x 10 x 10 mm) is clamped at x = 0 and carries a
uniform shear traction totalling P on the free end. The script sweeps element
types and mesh levels, runs CalculiX, and compares the tip deflection with
Timoshenko beam theory and with a Richardson-extrapolated converged value.

Usage
    python cantilever.py                       # run the CalculiX sweep + plots
    python cantilever.py --levels 1 2 4        # quicker sweep
    python cantilever.py --collect-abaqus      # also read Abaqus .dat results

Abaqus Learning Edition caps models at 1,000 nodes. Every case under the cap is
also written to runs/cantilever/abaqus/ with run_all.bat; copy that folder to the
Windows VM, run run_all.bat, copy the .dat files back, then use --collect-abaqus.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from fea_common import (ELEMENT_ORDER, box_hex_mesh, face_traction_loads,
                        parse_dat_displacements, run_ccx, write_deck)

# Geometry (mm), material (MPa), load (N)
L, H, B = 200.0, 10.0, 10.0
E, NU = 200_000.0, 0.3
P = 100.0
ABAQUS_LE_NODE_CAP = 1000

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE / "runs" / "cantilever"


def timoshenko_tip_deflection() -> dict[str, float]:
    """Euler-Bernoulli bending plus Timoshenko shear term (Cowper shear coefficient)."""
    I = B * H**3 / 12.0
    A = B * H
    G = E / (2.0 * (1.0 + NU))
    kappa = 10.0 * (1.0 + NU) / (12.0 + 11.0 * NU)  # Cowper (1966), rectangle
    bending = P * L**3 / (3.0 * E * I)
    shear = P * L / (kappa * G * A)
    return {"euler_bernoulli": bending, "timoshenko": bending + shear}


def build_case(etype: str, level: int):
    """Mesh, loads and supports for one element type at one refinement level."""
    ny = nz = level
    nx = int(round(L / H)) * level  # cube-shaped elements
    mesh = box_hex_mesh(L, H, B, nx, ny, nz, ELEMENT_ORDER[etype])
    s = 2 if mesh.order == 2 else 1
    fixed = mesh.node_set(lambda x, y, z: x < 1e-9)
    tip_loads = face_traction_loads(mesh, i_face=s * nx, total_force=P)
    cloads = {n: (2, f) for n, f in tip_loads.items()}  # DOF 2 = y direction
    tip_nodes = sorted(tip_loads)
    return mesh, fixed, cloads, tip_nodes, tip_loads


def tip_deflection(disp: dict[int, np.ndarray], tip_loads: dict[int, float]) -> float:
    """Load-weighted tip deflection: sum(f_i * v_i) / P.

    This is the displacement conjugate to the applied load (its work divided by P),
    so it is the same quantity beam theory predicts and is insensitive to the
    local distortion of the loaded face.
    """
    return sum(f * disp[n][1] for n, f in tip_loads.items()) / P


def richardson(h: list[float], d: list[float]) -> tuple[float, float]:
    """Extrapolated value and observed order from three solutions with a constant refinement ratio."""
    (h1, h2, h3), (d1, d2, d3) = h, d  # coarse -> fine
    r = h1 / h2
    p = np.log(abs((d1 - d2) / (d2 - d3))) / np.log(r)
    return d3 + (d3 - d2) / (r**p - 1.0), p


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", nargs="+", default=list(ELEMENT_ORDER))
    ap.add_argument("--levels", nargs="+", type=int, default=[1, 2, 4, 8])
    ap.add_argument("--ccx", default="ccx", help="CalculiX executable")
    ap.add_argument("--collect-abaqus", action="store_true",
                    help="read Abaqus results from runs/cantilever/abaqus/*.dat")
    args = ap.parse_args()

    ccx_dir, abq_dir = RUN_DIR / "ccx", RUN_DIR / "abaqus"
    ccx_dir.mkdir(parents=True, exist_ok=True)
    abq_dir.mkdir(parents=True, exist_ok=True)
    theory = timoshenko_tip_deflection()

    rows, abq_jobs = [], []
    for etype in args.elements:
        for lvl in args.levels:
            mesh, fixed, cloads, tip_nodes, tip_loads = build_case(etype, lvl)
            job = f"cant_{etype}_n{lvl}"
            heading = f"Cantilever benchmark, {etype}, level {lvl}, {mesh.n_nodes} nodes"
            deck = ccx_dir / f"{job}.inp"
            write_deck(deck, mesh, etype, E, NU, fixed, cloads, tip_nodes, heading)
            if mesh.n_nodes <= ABAQUS_LE_NODE_CAP:
                write_deck(abq_dir / f"{job}.inp", mesh, etype, E, NU, fixed, cloads, tip_nodes, heading)
                abq_jobs.append(job)

            delta = tip_deflection(parse_dat_displacements(run_ccx(deck, args.ccx)), tip_loads)
            row = {"solver": "ccx", "element": etype, "level": lvl, "nodes": mesh.n_nodes,
                   "dof": mesh.n_dof, "h_mm": H / lvl, "tip_mm": delta}
            rows.append(row)
            print(f"ccx    {etype:7s} n={lvl:<2d} nodes={mesh.n_nodes:6d}  tip={delta:.6f} mm"
                  f"  ({delta / theory['timoshenko']:.4f} x Timoshenko)")

            if args.collect_abaqus and (abq_dir / f"{job}.dat").exists():
                d_abq = tip_deflection(parse_dat_displacements(abq_dir / f"{job}.dat"), tip_loads)
                rows.append({**row, "solver": "abaqus", "tip_mm": d_abq})
                print(f"abaqus {etype:7s} n={lvl:<2d} nodes={mesh.n_nodes:6d}  tip={d_abq:.6f} mm"
                      f"  (ccx - abaqus = {(delta - d_abq) / d_abq:+.2e} relative)")

    # Batch file for the Windows VM
    bat = ["@echo off", "rem Run every cantilever case that fits the Learning Edition node cap"]
    bat += [f"call abaqus job={j} interactive" for j in abq_jobs]
    (abq_dir / "run_all.bat").write_text("\r\n".join(bat) + "\r\n")

    # Converged reference from the three finest quadratic (C3D20) meshes
    ref, p_obs = None, None
    q = sorted((r for r in rows if r["solver"] == "ccx" and r["element"] == "C3D20"),
               key=lambda r: r["h_mm"], reverse=True)
    if len(q) >= 3:
        ref, p_obs = richardson([r["h_mm"] for r in q[-3:]], [r["tip_mm"] for r in q[-3:]])

    for r in rows:
        r["ratio_to_timoshenko"] = r["tip_mm"] / theory["timoshenko"]
        r["rel_error_vs_extrapolated"] = abs(r["tip_mm"] - ref) / abs(ref) if ref else ""
    with open(RUN_DIR / "results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"\nEuler-Bernoulli: {theory['euler_bernoulli']:.6f} mm")
    print(f"Timoshenko:      {theory['timoshenko']:.6f} mm")
    if ref:
        print(f"Extrapolated FE (C3D20): {ref:.6f} mm, observed order {p_obs:.2f}, "
              f"{ref / theory['timoshenko']:.4f} x Timoshenko")
    print(f"Abaqus decks under {ABAQUS_LE_NODE_CAP} nodes: {len(abq_jobs)} -> {abq_dir}")
    make_plots(rows, theory, ref)


def make_plots(rows, theory, ref) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(ELEMENT_ORDER, ["#c0392b", "#e67e22", "#8e44ad", "#2471a3", "#17a589"]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    for et in ELEMENT_ORDER:
        rs = sorted((r for r in rows if r["solver"] == "ccx" and r["element"] == et), key=lambda r: r["dof"])
        if rs:
            ax.plot([r["dof"] for r in rs], [r["ratio_to_timoshenko"] for r in rs], "o-",
                    color=colors[et], label=et)
        ra = sorted((r for r in rows if r["solver"] == "abaqus" and r["element"] == et), key=lambda r: r["dof"])
        if ra:
            ax.plot([r["dof"] for r in ra], [r["ratio_to_timoshenko"] for r in ra], "x",
                    color=colors[et], ms=10, mew=2)
    ax.axhline(1.0, color="0.3", lw=1, ls="--")
    ax.text(0.01, 1.0, "Timoshenko beam theory", va="bottom", fontsize=8, color="0.3",
            transform=ax.get_yaxis_transform())
    if ref:
        ax.axhline(ref / theory["timoshenko"], color="0.3", lw=1, ls=":")
        ax.text(0.01, ref / theory["timoshenko"], "extrapolated FE", va="top", fontsize=8,
                color="0.3", transform=ax.get_yaxis_transform())
    # Keep the axis readable; label any point that falls off it (e.g. hourglassing)
    lo, hi = 0.6, 1.4
    for r in rows:
        if r["solver"] == "ccx" and not lo <= r["ratio_to_timoshenko"] <= hi:
            y = hi if r["ratio_to_timoshenko"] > hi else lo
            ax.annotate(f"{r['element']}, n={r['level']}: {r['ratio_to_timoshenko']:.1f}x (off scale)",
                        xy=(r["dof"], y), xytext=(r["dof"] * 9, y - 0.05 if y == hi else y + 0.05),
                        fontsize=8, color=colors[r["element"]],
                        arrowprops=dict(arrowstyle="->", color=colors[r["element"]], lw=1))
    ax.set_ylim(lo, hi)
    ax.set_xscale("log")
    ax.set_xlabel("degrees of freedom")
    ax.set_ylabel("tip deflection / Timoshenko")
    ax.set_title("Tip deflection (o = CalculiX, x = Abaqus)")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    if ref:
        for et in ELEMENT_ORDER:
            rs = sorted((r for r in rows if r["solver"] == "ccx" and r["element"] == et), key=lambda r: r["h_mm"])
            rs = [r for r in rs if r["rel_error_vs_extrapolated"] and r["rel_error_vs_extrapolated"] > 1e-9]
            if rs:
                ax.loglog([r["h_mm"] for r in rs], [r["rel_error_vs_extrapolated"] for r in rs], "o-",
                          color=colors[et], label=et)
        from matplotlib.ticker import FixedLocator, NullLocator, FormatStrFormatter
        hs = sorted({r["h_mm"] for r in rows})
        ax.xaxis.set_major_locator(FixedLocator(hs))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.set_xlabel("element size h (mm)")
        ax.set_ylabel("|error| vs extrapolated FE value")
        ax.set_title("Convergence")
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    out = RUN_DIR / "cantilever_convergence.png"
    fig.savefig(out, dpi=160)
    print(f"Plot: {out}")


if __name__ == "__main__":
    main()
