"""Why "C3D8" gives different answers in CalculiX and Abaqus.

An independent 8-node hex solver (numpy, ~50 lines) is run on the same mesh and
loads as cantilever.py in two variants:

  full   standard 2x2x2 Gauss integration         (what CalculiX calls C3D8)
  B-bar  same, but the volumetric strain at each  (what Abaqus/Standard calls C3D8:
         integration point is replaced by the      "selectively reduced integration")
         element-average volumetric strain

Each variant is compared against the matching solver's result for C3D8 at
n = 1 and n = 2 (nu = 0.3) and at n = 1 with nu = 0. Agreement to the printed
digits shows the solvers differ by formulation, not by a modelling error.

Usage
    uv run c3d8_formulation_check.py

CalculiX is run here. Abaqus results are read from .dat files when present:
nu = 0.3 from runs/cantilever/abaqus/ (written by cantilever.py), nu = 0 from
runs/cantilever/diag/abaqus/ (deck and run_diag.bat are written by this script).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import cantilever as cant
from fea_common import parse_dat_displacements, run_ccx, write_deck

DIAG = cant.RUN_DIR / "diag"
CASES = [(1, 0.3), (2, 0.3), (1, 0.0)]  # (mesh level, Poisson's ratio)

_G = 1.0 / np.sqrt(3.0)
GAUSS = [np.array((a, b, c)) for a in (-_G, _G) for b in (-_G, _G) for c in (-_G, _G)]
CORNERS = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
VOL = np.array([1, 1, 1, 0, 0, 0], float)  # Voigt vector picking the volumetric strain


def strain_matrix(X: np.ndarray, xi: np.ndarray) -> tuple[np.ndarray, float]:
    """B (6x24, Voigt order xx yy zz xy yz zx) and det J for a hex8 at natural point xi."""
    t = 1.0 + CORNERS * xi  # (8, 3): the three linear factors of each shape function
    dN = 0.125 * CORNERS * np.stack([t[:, 1] * t[:, 2], t[:, 0] * t[:, 2], t[:, 0] * t[:, 1]], 1)
    J = dN.T @ X
    dNdx = dN @ np.linalg.inv(J).T
    B = np.zeros((6, 24))
    for a, (x, y, z) in enumerate(dNdx):
        B[:, 3 * a:3 * a + 3] = [[x, 0, 0], [0, y, 0], [0, 0, z], [y, x, 0], [0, z, y], [z, 0, x]]
    return B, np.linalg.det(J)


def tip_deflection(level: int, nu: float, bbar: bool) -> float:
    E = cant.E
    lam, mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    D = np.diag([2 * mu] * 3 + [mu] * 3) + lam * np.outer(VOL, VOL)

    mesh, fixed, _, tip_nodes, tip_loads = cant.build_case("C3D8", level)
    index = {n: i for i, n in enumerate(sorted(mesh.nodes))}
    ndof = 3 * len(index)
    K = np.zeros((ndof, ndof))
    P_vol = np.outer(VOL, VOL) / 3.0  # projects a strain onto its volumetric part
    for conn in mesh.elements:
        X = np.array([mesh.nodes[n] for n in conn])
        pts = [strain_matrix(X, xi) for xi in GAUSS]
        if bbar:
            volume = sum(w for _, w in pts)
            B_vol_avg = sum(P_vol @ B * w for B, w in pts) / volume
            pts = [(B - P_vol @ B + B_vol_avg, w) for B, w in pts]
        dofs = np.array([3 * index[n] + k for n in conn for k in range(3)])
        K[np.ix_(dofs, dofs)] += sum(B.T @ D @ B * w for B, w in pts)

    F = np.zeros(ndof)
    for n, f in tip_loads.items():
        F[3 * index[n] + 1] = f
    free = np.setdiff1d(np.arange(ndof), [3 * index[n] + k for n in fixed for k in range(3)])
    u = np.zeros(ndof)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    disp = {n: u[3 * index[n]:3 * index[n] + 3] for n in tip_nodes}
    return cant.tip_deflection(disp, tip_loads)


def ccx_result(level: int, nu: float) -> float:
    mesh, fixed, cloads, tip_nodes, tip_loads = cant.build_case("C3D8", level)
    job = f"c3d8_n{level}_nu{nu:.1f}".replace(".", "p")
    deck = DIAG / "ccx" / f"{job}.inp"
    deck.parent.mkdir(parents=True, exist_ok=True)
    write_deck(deck, mesh, "C3D8", cant.E, nu, fixed, cloads, tip_nodes, f"C3D8 formulation check, {job}")
    return cant.tip_deflection(parse_dat_displacements(run_ccx(deck)), tip_loads)


def abaqus_result(level: int, nu: float) -> float | None:
    mesh, fixed, cloads, tip_nodes, tip_loads = cant.build_case("C3D8", level)
    if nu == cant.NU:
        dat = cant.RUN_DIR / "abaqus" / f"cant_C3D8_n{level}.dat"
    else:  # the main sweep only uses nu = 0.3, so write this deck here
        job = f"c3d8_n{level}_nu{nu:.1f}".replace(".", "p")
        abq = DIAG / "abaqus"
        abq.mkdir(parents=True, exist_ok=True)
        write_deck(abq / f"{job}.inp", mesh, "C3D8", cant.E, nu, fixed, cloads, tip_nodes,
                   f"C3D8 formulation check, {job}")
        (abq / "run_diag.bat").write_text(f"@echo off\r\ncall abaqus job={job} interactive\r\n")
        dat = abq / f"{job}.dat"
    return cant.tip_deflection(parse_dat_displacements(dat), tip_loads) if dat.exists() else None


def main() -> None:
    print(f"{'case':14s} {'numpy full':>11s} {'ccx C3D8':>11s} {'numpy B-bar':>12s} {'Abaqus C3D8':>12s}")
    for level, nu in CASES:
        full, bbar = tip_deflection(level, nu, False), tip_deflection(level, nu, True)
        ccx, abq = ccx_result(level, nu), abaqus_result(level, nu)
        abq_s = f"{abq:12.6f}" if abq is not None else f"{'(not run)':>12s}"
        print(f"n={level}, nu={nu:<4g}   {full:11.6f} {ccx:11.6f} {bbar:12.6f} {abq_s}")
    print("\nTip deflection in mm; Timoshenko reference "
          f"{cant.timoshenko_tip_deflection()['timoshenko']:.6f} mm at nu = 0.3.")


if __name__ == "__main__":
    main()
