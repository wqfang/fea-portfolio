"""Benchmark 1 figure: geometry, boundary conditions and variables of the cantilever.

Draws the beam as cantilever.py builds it, with the variables labelled:

  (a) 3D sketch (not to scale): length L, depth H, width B, clamped root,
      uniform shear traction on the free end totalling P, axes, E and nu
  (b) side view at true scale with the level-2 mesh (2 elements through the
      depth, cube-shaped elements), and the deflected shape (exaggerated) with
      the tip deflection measured by the script
  (c) the loaded end face: uniform traction tau = P / (B H) and the second
      moment of area I = B H^3 / 12

No solver run is needed.

Usage
    uv run cantilever_geometry.py        # writes runs/cantilever/cantilever_geometry.png
                                         # and cantilever_schematic.png
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle

import cantilever as cb

BODY_FC, BODY_C, EDGE = "#cfe0f3", "#2b6cb0", "#7a8796"
BC_C, LOAD_C, DIM_C = "#b7312c", "#2f855a", "#222222"
LEVEL = 2


def dim_arrow(ax, p, q, text, offset=(0, 0), colour=DIM_C, fs=10):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>", mutation_scale=9, color=colour, lw=1, zorder=5))
    mid = ((p[0] + q[0]) / 2 + offset[0], (p[1] + q[1]) / 2 + offset[1])
    ax.text(*mid, text, color=colour, fontsize=fs, ha="center", va="center", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))


def arrow(ax, p, q, colour, ms=9, lw=1.2):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=colour, lw=lw, zorder=5))


def hatch_wall(ax, x, y0, y1, side=-1, n=8, d=None):
    """Hatched clamp line at x from y0 to y1, hatching on `side` (-1 = left)."""
    d = d or 0.08 * (y1 - y0)
    ax.plot([x, x], [y0, y1], color="k", lw=2, zorder=4)
    for y in np.linspace(y0, y1, n):
        ax.plot([x, x + side * d], [y, y - d], color="k", lw=0.8, zorder=4)


def panel_sketch(ax) -> None:
    """Oblique 3D sketch; drawn length is shortened so the section stays readable."""
    Ld, Hd, Bd = 100.0, 14.0, 14.0           # drawn sizes (not to scale)
    ox, oy = 0.55 * np.cos(np.radians(35)), 0.55 * np.sin(np.radians(35))

    def P(x, y, z):
        return (x + z * ox, y + z * oy)

    front = [P(0, 0, 0), P(Ld, 0, 0), P(Ld, Hd, 0), P(0, Hd, 0)]          # z = 0 face
    top = [P(0, Hd, 0), P(Ld, Hd, 0), P(Ld, Hd, Bd), P(0, Hd, Bd)]
    end = [P(Ld, 0, 0), P(Ld, 0, Bd), P(Ld, Hd, Bd), P(Ld, Hd, 0)]
    for poly, fc in ((front, BODY_FC), (top, "#e3edf8"), (end, "#b9d0ea")):
        ax.add_patch(Polygon(poly, closed=True, fc=fc, ec=BODY_C, lw=1.2, zorder=2))

    # clamp at x = 0 (wall)
    wall = [P(0, -5, -5), P(0, -5, Bd + 5), P(0, Hd + 5, Bd + 5), P(0, Hd + 5, -5)]
    ax.add_patch(Polygon(wall, closed=True, fc="#d9d9d9", ec="k", lw=1, hatch="///", zorder=1))
    ax.text(*P(-3, Hd + 7, Bd / 2), "clamped root, x = 0:\n" + r"$u_x = u_y = u_z = 0$ (all nodes)",
            color=BC_C, fontsize=9.5, ha="right", va="bottom")

    # shear traction on the end face (+y)
    for y in (2.0, 7.0, 12.0):
        for z in (2.0, 7.0, 12.0):
            x0, y0 = P(Ld, y - 1.5, z)
            arrow(ax, (x0, y0), (x0, y0 + 5.0), LOAD_C, ms=7, lw=1)
    xe, ye = P(Ld, Hd, Bd)
    ax.text(xe + 4, ye + 5, "uniform shear traction on the free end\n"
            r"$\tau = P / (BH)$ in $+y$, total $P$ = " + f"{cb.P:g} N", color=LOAD_C, fontsize=9.5, va="bottom", ha="right")

    # dimensions
    dim_arrow(ax, P(0, -3, 0), P(Ld, -3, 0), f"L = {cb.L:g} mm", offset=(0, -2.2))
    xa, ya = P(Ld + 3, 0, 0)
    xb, yb = P(Ld + 3, Hd, 0)
    dim_arrow(ax, (xa + 9, ya), (xb + 9, yb), f"H = {cb.H:g}", offset=(6, 0))
    dim_arrow(ax, P(Ld + 2, -2, 0), P(Ld + 2, -2, Bd), f"B = {cb.B:g}", offset=(5.5, -1))
    ax.text(*P(Ld / 2, Hd / 2, 0), r"steel: $E$ = " + f"{cb.E / 1000:g} GPa, " + r"$\nu$ = " + f"{cb.NU:g}",
            fontsize=10, ha="center", va="center", color=BODY_C, weight="bold")

    # axes
    o = P(-14, -12, 0)
    for d, lab in (((9, 0), "x"), ((0, 9), "y"), ((9 * ox, 9 * oy), "z")):
        arrow(ax, o, (o[0] + d[0], o[1] + d[1]), "k", ms=8, lw=1)
        ax.text(o[0] + 1.25 * d[0], o[1] + 1.25 * d[1], lab, fontsize=10, ha="center", va="center")

    ax.set_xlim(-50, 165)
    ax.set_ylim(-16, 40)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("(a) Model (3D sketch, not to scale)", fontsize=10.5)


def panel_side(ax, theory) -> None:
    mesh, *_ = cb.build_case("C3D8", LEVEL)
    nx, ny, _ = mesh.shape
    L, H = cb.L, cb.H
    # side view (z = 0 face) of the mesh
    polys = [[(L * i / nx, H * j / ny), (L * (i + 1) / nx, H * j / ny), (L * (i + 1) / nx, H * (j + 1) / ny),
              (L * i / nx, H * (j + 1) / ny)] for i in range(nx) for j in range(ny)]
    ax.add_collection(PolyCollection(polys, facecolors=BODY_FC, edgecolors=EDGE, linewidths=0.5, zorder=1))
    hatch_wall(ax, 0, -6, H + 6, d=3)
    ax.text(-4, H + 8, r"$u = 0$", color=BC_C, fontsize=10, ha="right")

    # deflected neutral axis, exaggerated (Euler-Bernoulli shape)
    scale = 8.0
    d_tip = theory["timoshenko"]
    x = np.linspace(0, L, 200)
    v = d_tip * (x / L) ** 2 * (3 - x / L) / 2
    ax.plot(x, H / 2 + scale * v, ls="--", color=BODY_C, lw=1.4, zorder=3)
    ax.text(0.55 * L, H / 2 + scale * v[110] + 2.2, f"deflected axis (x{scale:g})", color=BODY_C, fontsize=9,
            ha="center")
    for j in range(5):                       # end traction
        y = H * (j + 0.5) / 5
        arrow(ax, (L + 2.5, y - 2.4), (L + 2.5, y + 2.4), LOAD_C, ms=7, lw=1)
    ax.text(L + 6, H / 2, r"$P$", color=LOAD_C, fontsize=12, va="center")
    dim_arrow(ax, (L - 4, H / 2), (L - 4, H / 2 + scale * d_tip), r"$\delta_{tip}$", offset=(-7, 0))
    dim_arrow(ax, (0, -9.5), (L, -9.5), f"L = {L:g} mm")
    dim_arrow(ax, (-12, 0), (-12, H), "H", offset=(-4, 0))
    h = H / ny
    dim_arrow(ax, (10, -3), (10 + h, -3), "h", fs=9, offset=(-5, 0))
    ax.text(10 + h + 3, -3, f"level n = {LEVEL}: n elements through H, h = H/n (cube-shaped elements)",
            fontsize=8.5, va="center")
    ax.set_xlim(-30, L + 18)
    ax.set_ylim(-14, H + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"(b) Side view, true scale (level-{LEVEL} mesh: {nx} x {ny} x {ny} elements)", fontsize=10.5)


def panel_section(ax) -> None:
    B, H = cb.B, cb.H
    ax.add_patch(Rectangle((0, 0), B, H, fc=BODY_FC, ec=BODY_C, lw=1.4))
    for k in range(1, LEVEL):
        ax.plot([B * k / LEVEL] * 2, [0, H], color=EDGE, lw=0.6)
        ax.plot([0, B], [H * k / LEVEL] * 2, color=EDGE, lw=0.6)
    for y in np.linspace(1.5, H - 1.5, 4):
        for z in np.linspace(1.5, B - 1.5, 4):
            arrow(ax, (z, y - 1.1), (z, y + 1.1), LOAD_C, ms=6, lw=0.9)
    dim_arrow(ax, (0, -1.5), (B, -1.5), "B", fs=10)
    dim_arrow(ax, (-1.5, 0), (-1.5, H), "H", fs=10)
    ax.text(B / 2, H + 1.4, r"$\tau = P/(BH)$" + f" = {cb.P / (cb.B * cb.H):g} MPa", ha="center", fontsize=10,
            color=LOAD_C)
    ax.text(B + 1.5, H / 2,
            r"$I = BH^3/12,\ \ A = BH$" + "\n" + r"$G = E\,/\,[2(1+\nu)]$" + "\n"
            r"$\kappa = 10(1+\nu)\,/\,(12+11\nu)$" + "  (Cowper)\n\n" + "Timoshenko tip deflection:\n"
            r"$\delta = PL^3/(3EI) + PL/(\kappa G A)$" + f"\n= {cb.timoshenko_tip_deflection()['timoshenko']:.4f} mm",
            fontsize=9.5, va="center", linespacing=1.5)
    ax.set_xlim(-4, B + 24)
    ax.set_ylim(-4, H + 3.5)
    ax.set_aspect("equal")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("(c) Free-end face (x = L)", fontsize=10.5)


def schematic() -> None:
    """Simple problem schematic: the beam as built, no mesh, key variables only."""
    fig, ax = plt.subplots(figsize=(13, 4.2))
    Ld, Hd = 120.0, 16.0                      # drawn sizes (beam not to scale)
    ax.add_patch(Rectangle((0, 0), Ld, Hd, fc=BODY_FC, ec=BODY_C, lw=1.6, zorder=2))
    hatch_wall(ax, 0, -7, Hd + 7, d=3, n=10)
    ax.text(-5, Hd + 9, "clamped: " + r"$u_x = u_y = u_z = 0$", color=BC_C, fontsize=11, ha="center", va="bottom")
    # tip load
    arrow(ax, (Ld, Hd / 2 - 13), (Ld, Hd / 2 - 0.5), LOAD_C, ms=14, lw=2.2)
    ax.text(Ld + 2, Hd / 2 - 11, f"P = {cb.P:g} N\n(uniform shear traction\non the end face)", color=LOAD_C,
            fontsize=10.5, va="center")
    # dimensions and material
    dim_arrow(ax, (0, -11), (Ld, -11), f"L = {cb.L:g} mm", fs=11)
    dim_arrow(ax, (-9, 0), (-9, Hd), f"H = {cb.H:g} mm", offset=(-9, 0), fs=11)
    ax.text(Ld / 2, Hd / 2, r"steel: $E$ = " + f"{cb.E / 1000:g} GPa, " + r"$\nu$ = " + f"{cb.NU:g}",
            fontsize=11, ha="center", va="center", color=BODY_C, weight="bold")
    # cross-section
    x0, y0, s = Ld + 48, 0.0, 16.0
    ax.add_patch(Rectangle((x0, y0), s, s, fc=BODY_FC, ec=BODY_C, lw=1.6))
    dim_arrow(ax, (x0, y0 - 3), (x0 + s, y0 - 3), f"B = {cb.B:g}", offset=(0, -2.5), fs=10)
    dim_arrow(ax, (x0 + s + 3, y0), (x0 + s + 3, y0 + s), f"H = {cb.H:g}", offset=(6, 0), fs=10)
    ax.text(x0 + s / 2, s + 3, "cross-section", ha="center", fontsize=10)
    # axes
    o = (-22, -16)
    for d, lab in (((8, 0), "x"), ((0, 8), "y")):
        arrow(ax, o, (o[0] + d[0], o[1] + d[1]), "k", ms=8, lw=1)
        ax.text(o[0] + 1.3 * d[0], o[1] + 1.3 * d[1], lab, fontsize=10, ha="center", va="center")
    ax.set_xlim(-32, x0 + s + 18)
    ax.set_ylim(-20, Hd + 15)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Benchmark 1: end-loaded cantilever (schematic, not to scale)", fontsize=12.5)
    out = cb.RUN_DIR / "cantilever_schematic.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")


def main() -> None:
    theory = cb.timoshenko_tip_deflection()
    fig = plt.figure(figsize=(15, 7.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.75, 1], height_ratios=[1.45, 1], hspace=0.18, wspace=0.08)
    panel_sketch(fig.add_subplot(gs[0, 0]))
    panel_section(fig.add_subplot(gs[0, 1]))
    panel_side(fig.add_subplot(gs[1, :]), theory)
    fig.suptitle("Benchmark 1: end-loaded cantilever beam", fontsize=13, y=0.97)
    out = cb.RUN_DIR / "cantilever_geometry.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")
    schematic()


if __name__ == "__main__":
    main()
