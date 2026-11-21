"""Benchmark 2 figure: geometry, boundary conditions and variables of the plate with a hole.

Draws both quarter models as plate_hole.py builds them (level-1 CPS4 meshes),
with the variables labelled:

  (a) Kirsch problem: square [0, C] x [0, C] round the hole of radius a, symmetry
      supports, the exact Kirsch tractions on the outer edges (arrows drawn to
      scale), the peak point (0, a) and the hoop angle theta
  (b) exact hoop stress round the hole edge, sigma_tt = S (1 - 2 cos 2 theta),
      with the peak 3S at theta = 90 deg
  (c) finite plate: width W = 2C, length 2 LH, uniform tension S on the end,
      symmetry supports, Peterson's Kt_net

No solver run is needed.

Usage
    uv run plate_hole_geometry.py        # writes runs/plate_hole/plate_hole_geometry.png
                                         # and plate_hole_schematic.png
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Rectangle

import plate_hole as ph

BODY_FC, BODY_C, EDGE = "#cfe0f3", "#2b6cb0", "#7a8796"
BC_C, LOAD_C, DIM_C, PEAK_C = "#b7312c", "#2f855a", "#222222", "#6b46c1"


def dim_arrow(ax, p, q, text, offset=(0, 0), colour=DIM_C, fs=10):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>", mutation_scale=9, color=colour, lw=1, zorder=5))
    mid = ((p[0] + q[0]) / 2 + offset[0], (p[1] + q[1]) / 2 + offset[1])
    ax.text(*mid, text, color=colour, fontsize=fs, ha="center", va="center", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))


def arrow(ax, p, q, colour, ms=8, lw=1.1):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=colour, lw=lw, zorder=5))


def draw_mesh(ax, xy, elems) -> None:
    polys = [[xy[n] for n in c[:4]] for c in elems]
    ax.add_collection(PolyCollection(polys, facecolors=BODY_FC, edgecolors=EDGE, linewidths=0.4, zorder=1))


def rollers(ax, pts, normal, n=6) -> None:
    """Roller symbols (circles) just outside an edge: normal = outward direction."""
    for x, y in pts[:: max(1, len(pts) // n)]:
        ax.plot(x + 1.6 * normal[0], y + 1.6 * normal[1], marker="o", ms=5, mfc="white", mec=BC_C, zorder=4)


def symmetry_supports(ax, xy, x_max, y_max) -> None:
    on_x0 = sorted((p for p in xy.values() if abs(p[0]) < 1e-9), key=lambda p: p[1])
    on_y0 = sorted((p for p in xy.values() if abs(p[1]) < 1e-9), key=lambda p: p[0])
    rollers(ax, on_x0, (-1, 0))
    rollers(ax, on_y0, (0, -1), n=8)
    ax.text(-3.5, 0.62 * y_max, r"symmetry: $u_x = 0$", rotation=90, color=BC_C, fontsize=9.5, ha="right",
            va="center")
    ax.text(0.62 * x_max, -3.8, r"symmetry: $u_y = 0$", color=BC_C, fontsize=9.5, ha="center", va="top")


def panel_kirsch(ax) -> None:
    a, C, S = ph.A, ph.C, ph.S
    xy, elems, _ = ph.build_mesh("kirsch", 1, 1)
    # the rest of the plate (not modelled)
    for sx, sy in ((-1, 1), (1, -1), (-1, -1)):
        ax.add_patch(Rectangle((min(0, sx * C), min(0, sy * C)), C, C, fc="#f3f3f3", ec="#c5ccd3", lw=0.8,
                               ls="--", zorder=0))
    ax.add_patch(Circle((0, 0), a, fc="white", ec="#9aa5b1", ls="--", lw=0.8, zorder=0))
    ax.text(-C / 2, -C / 2, "rest of the plate\n(not modelled:\nquarter symmetry)", ha="center", va="center",
            fontsize=8.5, color="#7b8794", style="italic")
    draw_mesh(ax, xy, elems)
    symmetry_supports(ax, xy, C, C)

    # Kirsch tractions on the outer edges, to scale: 10 mm per 100 MPa
    k = 0.1
    for y in np.linspace(2.5, C - 2.5, 10):
        sxx, _, sxy = ph.kirsch_stress(C, y)
        arrow(ax, (C, y), (C + k * sxx, y + k * sxy), LOAD_C, ms=6, lw=0.9)
    for x in np.linspace(2.5, C - 2.5, 10):
        _, syy, sxy = ph.kirsch_stress(x, C)
        if np.hypot(sxy, syy) > 4:
            arrow(ax, (x, C), (x + k * sxy, C + k * syy), LOAD_C, ms=5, lw=0.9)
    ax.text(C + 28, -9, "outer edges: exact\nKirsch tractions\n" + r"$t = \sigma^{Kirsch} n$" + "\n(arrows to scale)",
            color=LOAD_C, fontsize=9, ha="right", va="top")
    arrow(ax, (C + 13, 0.5 * C), (C + 21, 0.5 * C), LOAD_C, ms=10, lw=2)
    ax.text(C + 17, 0.5 * C + 2.5, f"remote\nS = {S:g} MPa", color=LOAD_C, fontsize=9, ha="center", va="bottom")

    # hole radius, theta, peak point
    t = np.radians(35)
    dim_arrow(ax, (0, 0), (a * np.cos(t), a * np.sin(t)), "a", fs=10, offset=(-0.6, 1.4))
    ax.add_patch(Arc((0, 0), 2.9 * a, 2.9 * a, theta1=0, theta2=62, color=DIM_C, lw=0.9))
    ax.text(1.5 * a * np.cos(np.radians(25)), 1.5 * a * np.sin(np.radians(25)) + 0.5, r"$\theta$", fontsize=11)
    ax.plot(0, a, marker="o", ms=7, color=PEAK_C, zorder=6)
    ax.annotate("peak: " + r"$\sigma_{xx}(0, a) = 3S$", xy=(0, a), xytext=(-46, 30), fontsize=9.5, color=PEAK_C,
                arrowprops=dict(arrowstyle="->", color=PEAK_C), zorder=7,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    ax.plot(a, 0, marker="s", ms=6, color=PEAK_C, zorder=6)
    ax.annotate(r"$u_x(a, 0) = 3Sa/E$", xy=(a, 0), xytext=(-46, 20), fontsize=9, color=PEAK_C,
                arrowprops=dict(arrowstyle="->", color=PEAK_C), zorder=7,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    dim_arrow(ax, (0, C + 7), (C, C + 7), f"C = {C:g} mm", offset=(0, 0))
    ax.text(0.66 * C, 0.8 * C, "plane stress\n" + r"$E,\ \nu$, thickness $T$", ha="center", va="center", fontsize=9.5,
            color=BODY_C, weight="bold", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
    ax.set_xlim(-C - 4, C + 30)
    ax.set_ylim(-C - 4, C + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("(a) Kirsch problem: exact solution known everywhere\n(level-1 CPS4 mesh, true scale)",
                 fontsize=10.5)


def panel_hoop(ax) -> None:
    S = ph.S
    th = np.linspace(0, 90, 300)
    stt = S * (1 - 2 * np.cos(np.radians(2 * th)))
    ax.axhline(0, color="#999", lw=0.8)
    ax.plot(th, stt, color=BODY_C, lw=2)
    ax.plot(90, 3 * S, marker="o", ms=7, color=PEAK_C)
    ax.annotate(r"$3S$ = " + f"{3 * S:g} MPa at the top of the hole (0, a)", xy=(90, 3 * S), xytext=(18, 2.75 * S),
                fontsize=9.5, color=PEAK_C, arrowprops=dict(arrowstyle="->", color=PEAK_C))
    ax.annotate(r"$-S$ at $(a, 0)$", xy=(0, -S), xytext=(8, -0.55 * S), fontsize=9.5, color=DIM_C,
                arrowprops=dict(arrowstyle="->", color="#555"))
    ax.text(4, 1.9 * S, r"$\sigma_{\theta\theta}(a, \theta) = S\,(1 - 2\cos 2\theta)$", fontsize=11, color=BODY_C,
            ha="left")
    ax.text(88, -0.95 * S, "Kirsch: infinite plate, remote S along x.\n"
            r"$K_t = \sigma_{max}/S = 3$; displacements at the hole:" + "\n"
            r"$u_x(a,0) = 3Sa/E$,  $u_y(0,a) = -Sa/E$", fontsize=9, ha="right", va="bottom",
            bbox=dict(fc="white", ec="#bbb", pad=4))
    ax.set_xlim(0, 90)
    ax.set_ylim(-1.15 * S, 3.3 * S)
    ax.set_xticks(range(0, 91, 15))
    ax.set_xlabel(r"$\theta$ (deg), measured from the x axis")
    ax.set_ylabel("hoop stress at the hole edge (MPa)")
    ax.set_title("(b) Exact hoop stress round the hole edge", fontsize=10.5)


def panel_plate(ax) -> None:
    a, C, LH, S = ph.A, ph.C, ph.LH, ph.S
    xy, elems, _ = ph.build_mesh("plate", 1, 1)
    draw_mesh(ax, xy, elems)
    symmetry_supports(ax, xy, LH, C)
    for y in np.linspace(3, C - 3, 7):
        arrow(ax, (LH, y), (LH + 9, y), LOAD_C, ms=7)
    ax.text(LH + 11, C / 2, f"uniform\ntension\nS = {S:g} MPa", color=LOAD_C, fontsize=9.5, va="center")
    dim_arrow(ax, (0, C + 7), (LH, C + 7), f"L/2 = {LH:g} mm  (L = {2 * LH:g}, L/W = {LH / C:g})")
    dim_arrow(ax, (LH + 36, 0), (LH + 36, C), f"W/2 = {C:g}", offset=(0, 0))
    t = np.radians(40)
    dim_arrow(ax, (0, 0), (a * np.cos(t), a * np.sin(t)), "a", fs=10, offset=(-1, 1.6))
    ax.plot(0, a, marker="o", ms=7, color=PEAK_C, zorder=6)
    ax.annotate(r"peak $\sigma_{xx}(0, a)$", xy=(0, a), xytext=(14, 33), fontsize=9.5, color=PEAK_C,
                arrowprops=dict(arrowstyle="->", color=PEAK_C), bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    ax.text(0.55 * LH, 0.5 * C,
            r"$K_{t,net} = \sigma_{max}\,(1 - d/W)\,/\,S$" + "\n"
            r"Peterson: $K_{t,net} = 3.00 - 3.13q + 3.66q^2 - 1.53q^3$, $q = d/W$" + "\n"
            f"d = {2 * a:g} mm, W = {2 * C:g} mm: q = {a / C:g}, Peterson " + r"$K_{t,net}$" +
            f" = {ph.peterson_kt_net(a / C):.3f}",
            ha="center", va="center", fontsize=9.5, bbox=dict(fc="white", ec="#bbb", pad=4), zorder=6)
    ax.set_xlim(-12, LH + 45)
    ax.set_ylim(-9, C + 13)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("(c) Finite plate, quarter model (level-1 CPS4 mesh, true scale)", fontsize=10.5)


def schematic() -> None:
    """Simple problem schematic: the full plates, no mesh, key variables only."""
    a, C, LH, S = ph.A, ph.C, ph.LH, ph.S
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw=dict(width_ratios=[1, 1.55], wspace=0.08))

    # (a) infinite plate (Kirsch)
    Rp = 55.0
    ax1.add_patch(Rectangle((-Rp, -Rp), 2 * Rp, 2 * Rp, fc=BODY_FC, ec="none", zorder=0))
    for side in (-1, 1):                       # fade the edges: the plate is unbounded
        for k, al in enumerate((0.55, 0.35, 0.18)):
            ax1.add_patch(Rectangle((-Rp - 1, side * (Rp - 4 * k) - (4 if side > 0 else 0)), 2 * Rp + 2, 4,
                                    fc="white", ec="none", alpha=al, zorder=1))
    ax1.add_patch(Circle((0, 0), a, fc="white", ec=BODY_C, lw=1.6, zorder=2))
    for y in np.linspace(-Rp + 8, Rp - 8, 7):
        arrow(ax1, (Rp + 2, y), (Rp + 14, y), LOAD_C, ms=8)
        arrow(ax1, (-Rp - 2, y), (-Rp - 14, y), LOAD_C, ms=8)
    ax1.text(Rp + 8, Rp + 2, f"S = {S:g} MPa", color=LOAD_C, fontsize=11, ha="center", va="bottom")
    ax1.text(-Rp - 8, Rp + 2, "S", color=LOAD_C, fontsize=11, ha="center", va="bottom")
    t = np.radians(40)
    dim_arrow(ax1, (0, 0), (a * np.cos(t), a * np.sin(t)), "a", fs=11, offset=(-1.5, 2))
    ax1.text(0, a + 1.5, "", fontsize=1)
    ax1.add_patch(Arc((0, 0), 2.8 * a, 2.8 * a, theta1=0, theta2=60, color=DIM_C, lw=0.9))
    ax1.plot([a, 1.6 * a], [0, 0], color=DIM_C, lw=0.9)
    ax1.text(1.55 * a * np.cos(np.radians(28)), 1.55 * a * np.sin(np.radians(28)), r"$\theta$", fontsize=12)
    ax1.plot(0, a, marker="o", ms=7, color=PEAK_C, zorder=4)
    ax1.annotate(r"$\sigma_{max} = 3S$", xy=(0, a), xytext=(-40, 30), fontsize=11, color=PEAK_C,
                 arrowprops=dict(arrowstyle="->", color=PEAK_C))
    ax1.text(0, -36, "infinite plate, plane stress\n" + r"$E$ = " + f"{ph.E / 1000:g} GPa, " + r"$\nu$ = "
             + f"{ph.NU:g}", ha="center", fontsize=10.5, color=BODY_C, weight="bold")
    o = (-Rp + 4, -Rp + 4)
    for d, lab in (((9, 0), "x"), ((0, 9), "y")):
        arrow(ax1, o, (o[0] + d[0], o[1] + d[1]), "k", ms=8, lw=1)
        ax1.text(o[0] + 1.3 * d[0], o[1] + 1.3 * d[1], lab, fontsize=10, ha="center", va="center")
    ax1.set_xlim(-Rp - 18, Rp + 18)
    ax1.set_ylim(-Rp - 6, Rp + 10)
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax1.set_title(f"(a) Kirsch problem: hole of radius a = {a:g} mm\nin an infinite plate under remote tension S",
                  fontsize=11)

    # (b) finite plate, to scale
    ax2.add_patch(Rectangle((-LH, -C), 2 * LH, 2 * C, fc=BODY_FC, ec=BODY_C, lw=1.6, zorder=1))
    ax2.add_patch(Circle((0, 0), a, fc="white", ec=BODY_C, lw=1.6, zorder=2))
    for y in np.linspace(-C + 8, C - 8, 6):
        arrow(ax2, (LH + 2, y), (LH + 22, y), LOAD_C, ms=8)
        arrow(ax2, (-LH - 2, y), (-LH - 22, y), LOAD_C, ms=8)
    ax2.text(LH + 12, C + 4, f"S = {S:g} MPa", color=LOAD_C, fontsize=11, ha="center", va="bottom")
    ax2.text(-LH - 12, C + 4, "S", color=LOAD_C, fontsize=11, ha="center", va="bottom")
    dim_arrow(ax2, (-LH, -C - 14), (LH, -C - 14), f"L = {2 * LH:g} mm", fs=11)
    dim_arrow(ax2, (LH + 32, -C), (LH + 32, C), f"W = {2 * C:g} mm", offset=(0, 0), fs=11)
    dim_arrow(ax2, (-a, -a - 4), (a, -a - 4), f"d = {2 * a:g}", fs=10, offset=(0, -5))
    ax2.plot(0, a, marker="o", ms=7, color=PEAK_C, zorder=4)
    ax2.annotate(r"$\sigma_{max}$", xy=(0, a), xytext=(28, 30), fontsize=11, color=PEAK_C,
                 arrowprops=dict(arrowstyle="->", color=PEAK_C))
    ax2.text(-LH / 2, 0, "plane stress\n" + r"$E,\ \nu$, thickness $T$", ha="center", va="center", fontsize=10.5,
             color=BODY_C, weight="bold")
    ax2.text(LH / 2 + 5, -C / 2,
             r"$K_{t,net} = \sigma_{max}(1 - d/W)/S$" + "\n" + f"d/W = {a / C:g}, L/W = {LH / C:g}",
             ha="center", va="center", fontsize=10.5)
    ax2.set_xlim(-LH - 30, LH + 48)
    ax2.set_ylim(-C - 22, C + 14)
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.set_title("(b) Finite plate with a central hole under uniform end tension (to scale)", fontsize=11)
    fig.suptitle("Benchmark 2: plate with a circular hole in tension (schematic)", fontsize=12.5, y=1.0)
    out = ph.RUN_DIR / "plate_hole_schematic.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")


def main() -> None:
    fig = plt.figure(figsize=(15, 13))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1], height_ratios=[1.35, 1], hspace=0.2, wspace=0.18)
    panel_kirsch(fig.add_subplot(gs[0, 0]))
    panel_hoop(fig.add_subplot(gs[0, 1]))
    panel_plate(fig.add_subplot(gs[1, :]))
    fig.suptitle("Benchmark 2: plate with a circular hole in tension (plane stress)", fontsize=13, y=0.94)
    out = ph.RUN_DIR / "plate_hole_geometry.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")
    schematic()


if __name__ == "__main__":
    main()
