"""Benchmark 2: plate with a circular hole in tension, stress concentration factor Kt.

Plane stress, quarter model with symmetry, remote tension S along x. Two problems:

  kirsch  A square region [0, C] x [0, C] around the hole (radius A) loaded on its
          outer edges with the exact Kirsch infinite-plate stresses. The exact
          solution is then known everywhere: peak stress 3S at the top of the hole,
          hoop stress S(1 - 2 cos 2theta) round the hole edge.
  plate   A finite plate, width W = 2C, length 2 * LH, uniform tension S on the
          ends. Kt is compared with Peterson's handbook fit and with a
          Richardson-extrapolated FE value.

Peak stress is the averaged nodal sigma_xx at (0, A), the top of the hole.

Usage
    uv run plate_hole.py                       # CalculiX sweep + plots
    uv run plate_hole.py --levels 1 2 4        # quicker
    uv run plate_hole.py --collect-abaqus      # also read Abaqus .dat results

Abaqus Learning Edition caps models at 1,000 nodes. Every case under the cap is
written to runs/plate_hole/abaqus/ with run_all.bat (including CPS4I, which ccx
does not support); copy that folder to the Windows VM, run run_all.bat, copy the
.dat files back, then use --collect-abaqus.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from fea_common import (ELEMENT_ORDER_2D, parse_dat_nodal_stress_2d, parse_frd_nodal_stress,
                        run_ccx, write_deck_2d)

# Geometry (mm), material (MPa), load (MPa)
A = 10.0          # hole radius (d = 20)
C = 50.0          # half-width (W = 100, d/W = 0.2); also the Kirsch square half-size
LH = 150.0        # half-length of the finite plate (L = 300, L/W = 3)
T = 0.01          # thickness: see note below
# CalculiX does not solve 2D elements directly: it expands each plane-stress element
# into a 3D element of the section thickness, so it really models a thin 3D plate.
# Near a hole a real plate is not in exact plane stress (sigma_zz != 0 inside), and
# the mid-plane peak rises with thickness/radius: +0.2% at T = 1 mm, +0.85% at
# T = 5 mm. A very thin section recovers plane stress. Abaqus solves true plane
# stress, where thickness only scales the loads, so the same T is used for both.
E, NU = 200_000.0, 0.3
S = 100.0         # remote tension
ABAQUS_LE_NODE_CAP = 1000
BETA = 3.3        # radial grading: first ring roughly square at the hole

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE / "runs" / "plate_hole"
CCX_ELEMENTS = ["CPS4", "CPS4R", "CPS8", "CPS8R"]


# --------------------------------------------------------------------------- reference solutions

def kirsch_stress(x, y):
    """Exact (sxx, syy, sxy) for an infinite plate with a hole of radius A under remote sxx = S."""
    r2 = x * x + y * y
    th = np.arctan2(y, x)
    a2, a4 = A * A / r2, A**4 / r2**2
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    srr = S / 2 * (1 - a2) + S / 2 * (1 - 4 * a2 + 3 * a4) * c2
    stt = S / 2 * (1 + a2) - S / 2 * (1 + 3 * a4) * c2
    srt = -S / 2 * (1 + 2 * a2 - 3 * a4) * s2
    c, s = np.cos(th), np.sin(th)
    sxx = srr * c * c + stt * s * s - 2 * srt * s * c
    syy = srr * s * s + stt * c * c + 2 * srt * s * c
    sxy = (srr - stt) * s * c + srt * (c * c - s * s)
    return sxx, syy, sxy


def peterson_kt_net(d_over_w: float) -> float:
    """Peterson's curve fit for a central hole in a finite-width plate in tension (net section)."""
    q = d_over_w
    return 3.00 - 3.13 * q + 3.66 * q**2 - 1.53 * q**3


# --------------------------------------------------------------------------- mesh

class Registry:
    """Node numbering with merging of coincident points (shared block edges)."""

    def __init__(self):
        self.ids: dict[tuple[float, float], int] = {}
        self.xy: dict[int, tuple[float, float]] = {}

    def node(self, x: float, y: float) -> int:
        key = (round(x, 8), round(y, 8))
        if key not in self.ids:
            nid = len(self.ids) + 1
            self.ids[key] = nid
            self.xy[nid] = tuple(0.0 if abs(v) < 1e-9 else float(v) for v in (x, y))
        return self.ids[key]


def grade(s: float) -> float:
    return (np.exp(BETA * s) - 1) / (np.exp(BETA) - 1)


def block(reg: Registry, point, ni: int, nj: int, order: int):
    """Structured quad block. point(s, t) maps the unit square; returns elements and the
    node grid (indexed in half-steps for quadratic elements)."""
    st = 2 if order == 2 else 1
    grid = {}
    for I in range(st * ni + 1):
        for J in range(st * nj + 1):
            if order == 2 and I % 2 and J % 2:
                continue  # serendipity: no centre node
            grid[(I, J)] = reg.node(*point(I / (st * ni), J / (st * nj)))
    elems = []
    for i in range(ni):
        for j in range(nj):
            I, J = st * i, st * j
            if order == 1:
                c = [(I, J), (I + 1, J), (I + 1, J + 1), (I, J + 1)]
            else:
                c = [(I, J), (I + 2, J), (I + 2, J + 2), (I, J + 2),
                     (I + 1, J), (I + 2, J + 1), (I + 1, J + 2), (I, J + 1)]
            elems.append([grid[p] for p in c])
    return elems, grid


def edges_along(grid, fixed_I: int | None, fixed_J: int | None, n: int, order: int):
    """Element edges (corner, corner[, mid]) along a grid line of a block."""
    st = 2 if order == 2 else 1
    out = []
    for k in range(n):
        K = st * k
        pts = [K, K + st] + ([K + 1] if order == 2 else [])
        out.append([grid[(fixed_I, p)] if fixed_I is not None else grid[(p, fixed_J)] for p in pts])
    return out


def build_mesh(problem: str, level: int, order: int):
    """Quarter model. O-grid around the hole: i radial, j round the hole (theta 0 -> 90 deg).
    The outer edge of the O-grid runs along x = C (theta < 45 deg) and y = C (theta > 45 deg).
    The finite plate adds a block from x = C to x = LH."""
    n_th, n_r = 8 * level, 4 * level
    reg = Registry()

    def ogrid(s, t):
        th = 0.5 * np.pi * t
        inner = np.array([A * np.cos(th), A * np.sin(th)])
        outer = np.array([C, 2 * C * t]) if t <= 0.5 else np.array([C * (2 - 2 * t), C])
        return inner + grade(s) * (outer - inner)

    st = 2 if order == 2 else 1
    elems, g = block(reg, ogrid, n_r, n_th, order)
    outer_edges = edges_along(g, st * n_r, None, n_th, order)
    half = len(outer_edges) // 2
    load_edges = []  # (edges, traction function)
    if problem == "kirsch":
        def on_x(x, y):  # outward normal (1, 0)
            sxx, _, sxy = kirsch_stress(x, y)
            return sxx, sxy

        def on_y(x, y):  # outward normal (0, 1)
            _, syy, sxy = kirsch_stress(x, y)
            return sxy, syy
        load_edges += [(outer_edges[:half], on_x), (outer_edges[half:], on_y)]
    else:
        n_x = 4 * level
        e2, g2 = block(reg, lambda s, t: (C + (LH - C) * s, C * t), n_x, n_th // 2, order)
        elems += e2
        load_edges.append((edges_along(g2, st * n_x, None, n_th // 2, order), lambda x, y: (S, 0.0)))
    return reg.xy, elems, load_edges


def edge_loads(xy, load_edges, order: int) -> dict[tuple[int, int], float]:
    """Consistent nodal forces for tractions on element edges (6-point Gauss per edge)."""
    gp, gw = np.polynomial.legendre.leggauss(6)
    f: dict[tuple[int, int], float] = {}
    for edges, traction in load_edges:
        for e in edges:
            X = np.array([xy[n] for n in e])
            for xi, w in zip(gp, gw):
                if order == 1:
                    N, dN = np.array([(1 - xi) / 2, (1 + xi) / 2]), np.array([-0.5, 0.5])
                else:
                    N = np.array([xi * (xi - 1) / 2, xi * (xi + 1) / 2, 1 - xi * xi])
                    dN = np.array([xi - 0.5, xi + 0.5, -2 * xi])
                x, y = N @ X
                ds = np.linalg.norm(dN @ X)
                tx, ty = traction(x, y)
                for a, n in enumerate(e):
                    f[(n, 1)] = f.get((n, 1), 0.0) + N[a] * tx * ds * w * T
                    f[(n, 2)] = f.get((n, 2), 0.0) + N[a] * ty * ds * w * T
    return f


def build_case(problem: str, etype: str, level: int):
    order = ELEMENT_ORDER_2D[etype]
    xy, elems, load_edges = build_mesh(problem, level, order)
    loads = edge_loads(xy, load_edges, order)
    cloads = [(n, d, v) for (n, d), v in sorted(loads.items()) if abs(v) > 1e-12]
    bcs = {"XSYM": (sorted(n for n, (x, _) in xy.items() if abs(x) < 1e-9), 1),
           "YSYM": (sorted(n for n, (_, y) in xy.items() if abs(y) < 1e-9), 2)}
    hole = sorted((n for n, (x, y) in xy.items() if abs(np.hypot(x, y) - A) < 1e-7),
                  key=lambda n: np.arctan2(xy[n][1], xy[n][0]))
    peak = next(n for n in hole if abs(xy[n][0]) < 1e-9)
    return xy, elems, bcs, cloads, hole, peak


# --------------------------------------------------------------------------- results

def hole_metrics(xy, stress, hole, peak) -> dict[str, float]:
    """Peak sigma_xx at the top of the hole, and the hoop-stress error round the hole edge
    (relative to the exact Kirsch profile; meaningful for the kirsch problem)."""
    sxx_peak = stress[peak][0]
    err = []
    for n in hole:
        x, y = xy[n]
        th = np.arctan2(y, x)
        sxx, syy, sxy = stress[n]
        hoop = sxx * np.sin(th)**2 + syy * np.cos(th)**2 - 2 * sxy * np.sin(th) * np.cos(th)
        err.append(hoop - S * (1 - 2 * np.cos(2 * th)))
    return {"sigma_peak": sxx_peak, "hoop_rms_err": float(np.sqrt(np.mean(np.square(err))))}


def metrics(problem, xy, stress, hole, peak) -> dict:
    m = hole_metrics(xy, stress, hole, peak)
    if problem != "kirsch":
        m["hoop_rms_err"] = ""  # the Kirsch profile is not the finite-plate answer
    return m


def richardson(h, d):
    (h1, h2, h3), (d1, d2, d3) = h, d
    r = h1 / h2
    p = np.log(abs((d1 - d2) / (d2 - d3))) / np.log(r)
    return d3 + (d3 - d2) / (r**p - 1.0), p


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", nargs="+", default=list(ELEMENT_ORDER_2D))
    ap.add_argument("--levels", nargs="+", type=int, default=[1, 2, 4, 8])
    ap.add_argument("--problems", nargs="+", default=["kirsch", "plate"])
    ap.add_argument("--ccx", default="ccx", help="CalculiX executable")
    ap.add_argument("--collect-abaqus", action="store_true",
                    help="read Abaqus results from runs/plate_hole/abaqus/*.dat")
    args = ap.parse_args()

    ccx_dir, abq_dir = RUN_DIR / "ccx", RUN_DIR / "abaqus"
    ccx_dir.mkdir(parents=True, exist_ok=True)
    abq_dir.mkdir(parents=True, exist_ok=True)
    kt_ref = peterson_kt_net(2 * A / (2 * C))
    print(f"Kirsch exact peak 3S = {3 * S:.1f} MPa; Peterson Kt_net (d/W = {A / C:.2f}) = {kt_ref:.4f}"
          f" -> peak {kt_ref * S / (1 - A / C):.2f} MPa\n")

    rows, abq_jobs = [], []
    for problem in args.problems:
        for etype in args.elements:
            for lvl in args.levels:
                xy, elems, bcs, cloads, hole, peak = build_case(problem, etype, lvl)
                job = f"{problem}_{etype}_n{lvl}"
                head = f"Plate with hole, {problem}, {etype}, level {lvl}, {len(xy)} nodes"
                base = {"problem": problem, "element": etype, "level": lvl, "nodes": len(xy),
                        "h_mm": A * np.pi / 2 / (8 * lvl)}
                if len(xy) <= ABAQUS_LE_NODE_CAP:
                    write_deck_2d(abq_dir / f"{job}.inp", xy, elems, etype, E, NU, T, bcs, cloads, head, "abaqus")
                    abq_jobs.append(job)
                    dat = abq_dir / f"{job}.dat"
                    if args.collect_abaqus and dat.exists():
                        m = metrics(problem, xy, parse_dat_nodal_stress_2d(dat), hole, peak)
                        rows.append({**base, "solver": "abaqus", **m})
                        print(f"abaqus {problem:6s} {etype:5s} n={lvl:<2d} nodes={len(xy):6d}"
                              f"  peak={m['sigma_peak']:9.3f}")
                if etype in CCX_ELEMENTS:
                    deck = ccx_dir / f"{job}.inp"
                    write_deck_2d(deck, xy, elems, etype, E, NU, T, bcs, cloads, head, "ccx")
                    run_ccx(deck, args.ccx)
                    frd = parse_frd_nodal_stress(deck.with_suffix(".frd"))
                    m = metrics(problem, xy, {n: v[[0, 1, 3]] for n, v in frd.items()}, hole, peak)  # sxx syy sxy
                    rows.append({**base, "solver": "ccx", **m})
                    print(f"ccx    {problem:6s} {etype:5s} n={lvl:<2d} nodes={len(xy):6d}"
                          f"  peak={m['sigma_peak']:9.3f}"
                          + (f"  hoop rms err={m['hoop_rms_err']:.3e}" if m["hoop_rms_err"] != "" else ""), flush=True)

    bat = ["@echo off", "rem Run every plate-with-hole case that fits the Learning Edition node cap"]
    bat += [f"call abaqus job={j} interactive" for j in abq_jobs]
    (abq_dir / "run_all.bat").write_text("\r\n".join(bat) + "\r\n")

    # References: exact for kirsch; Richardson on the three finest ccx CPS8 meshes for plate
    ref = {"kirsch": 3 * S}
    q = sorted((r for r in rows if r["solver"] == "ccx" and r["problem"] == "plate" and r["element"] == "CPS8"),
               key=lambda r: r["h_mm"], reverse=True)
    if len(q) >= 3:
        ref["plate"], p_obs = richardson([r["h_mm"] for r in q[-3:]], [r["sigma_peak"] for r in q[-3:]])
        print(f"\nPlate, Richardson (CPS8): peak {ref['plate']:.3f} MPa, observed order {p_obs:.2f},"
              f" Kt_net = {ref['plate'] * (1 - A / C) / S:.4f} (Peterson {kt_ref:.4f})")
    for r in rows:
        r["kt_gross"] = r["sigma_peak"] / S
        r["kt_net"] = r["sigma_peak"] * (1 - A / C) / S
        r["peak_rel_err"] = (r["sigma_peak"] - ref[r["problem"]]) / ref[r["problem"]] if r["problem"] in ref else ""
    if rows:
        with open(RUN_DIR / "results.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print(f"Abaqus decks under {ABAQUS_LE_NODE_CAP} nodes: {len(abq_jobs)} -> {abq_dir}")
    if rows:
        make_plots(rows, ref, kt_ref)


def make_plots(rows, ref, kt_ref) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(ELEMENT_ORDER_2D, ["#c0392b", "#e67e22", "#8e44ad", "#2471a3", "#17a589"]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    for et in ELEMENT_ORDER_2D:
        for solver, mk, ls in (("ccx", "o", "-"), ("abaqus", "x", "--")):
            rs = sorted((r for r in rows if r["problem"] == "kirsch" and r["element"] == et
                         and r["solver"] == solver), key=lambda r: r["h_mm"])
            if rs:
                ax.loglog([r["h_mm"] for r in rs], [abs(r["peak_rel_err"]) for r in rs], marker=mk, ls=ls,
                          color=colors[et], label=et if solver == "ccx" or et == "CPS4I" else None)
    ax.set_xlabel("element size round the hole (mm)")
    ax.set_ylabel("|peak stress error| vs exact 3S")
    ax.set_title("Kirsch problem (o = CalculiX, x = Abaqus)")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    for et in ELEMENT_ORDER_2D:
        for solver, mk, ls in (("ccx", "o", "-"), ("abaqus", "x", "--")):
            rs = sorted((r for r in rows if r["problem"] == "plate" and r["element"] == et
                         and r["solver"] == solver), key=lambda r: r["nodes"])
            if rs:
                ax.semilogx([r["nodes"] for r in rs], [r["kt_net"] for r in rs], marker=mk, ls=ls,
                            color=colors[et], label=et if solver == "ccx" or et == "CPS4I" else None)
    ax.axhline(kt_ref, color="0.3", ls=":", lw=1)
    ax.text(ax.get_xlim()[0] * 1.1, kt_ref - 0.004, "Peterson fit", fontsize=8, color="0.3", va="top")
    if "plate" in ref:
        kt_fe = ref["plate"] * (1 - A / C) / S
        ax.axhline(kt_fe, color="0.3", ls="--", lw=1)
        ax.text(ax.get_xlim()[0] * 1.1, kt_fe + 0.004, "extrapolated FE", fontsize=8, color="0.3")
    ax.set_ylim(kt_ref - 0.06, kt_ref + 0.06)
    ax.set_xlabel("nodes")
    ax.set_ylabel("Kt (net section)")
    ax.set_title("Finite plate, d/W = 0.2")
    ax.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    out = RUN_DIR / "plate_hole_convergence.png"
    fig.savefig(out, dpi=150)
    print(f"Plot: {out}")


if __name__ == "__main__":
    main()
