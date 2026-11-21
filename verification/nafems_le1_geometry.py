"""Benchmark 4 figures: NAFEMS LE1 elliptic membrane.

  schematic   the whole membrane (original geometry, no mesh): the two ellipses,
              outward pressure on the outer edge, free inner edge, the quarter
              that is modelled, point D and the target stress
  geometry    the quarter model as meshed: level-2 mesh (24 elements, as the
              NAFEMS fine mesh), points A-D, symmetry supports, the pressure as
              applied, how mesh lines are built from the parametric angle t, and
              the meshes at levels 1, 4 and 8

No solver run is needed.

Usage
    uv run nafems_le1_geometry.py   # writes runs/nafems_le1/nafems_le1_schematic.png
                                    # and nafems_le1_geometry.png
"""
from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch

import nafems_le1 as le

BODY_FC, BODY_C, EDGE = "#cfe0f3", "#2b6cb0", "#7a8796"
BC_C, LOAD_C, DIM_C, PEAK_C = "#b7312c", "#2f855a", "#222222", "#6b46c1"


def dim_arrow(ax, p, q, text, offset=(0, 0), colour=DIM_C, fs=10):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>", mutation_scale=9, color=colour, lw=1, zorder=5))
    mid = ((p[0] + q[0]) / 2 + offset[0], (p[1] + q[1]) / 2 + offset[1])
    ax.text(*mid, text, color=colour, fontsize=fs, ha="center", va="center", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))


def arrow(ax, p, q, colour, ms=8, lw=1.1):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=colour, lw=lw, zorder=5))


def ellipse(a, b, t0=0.0, t1=2 * math.pi, n=400):
    t = np.linspace(t0, t1, n)
    return a * np.cos(t), b * np.sin(t)


def pressure_arrows(ax, t_values, length=420.0, ms=8):
    for t in t_values:
        x, y = le.ellipse_point(le.A_OUT, le.B_OUT, t)
        n = le.outward_normal(le.A_OUT, le.B_OUT, x, y)
        arrow(ax, (x + 40 * n[0], y + 40 * n[1]), (x + length * n[0], y + length * n[1]), LOAD_C, ms=ms)


def draw_mesh(ax, xy, elems, lw=0.5) -> None:
    polys = []
    for c in elems:
        if len(c) == 8:   # follow the midside nodes so curved edges show
            order = [0, 4, 1, 5, 2, 6, 3, 7]
            polys.append([xy[c[k]] for k in order])
        else:
            polys.append([xy[n] for n in c])
    ax.add_collection(PolyCollection(polys, facecolors=BODY_FC, edgecolors=EDGE, linewidths=lw, zorder=1))


# --------------------------------------------------------------------------- schematic

def schematic() -> None:
    fig, ax = plt.subplots(figsize=(11, 9.4))
    xo, yo = ellipse(le.A_OUT, le.B_OUT)
    xi, yi = ellipse(le.A_IN, le.B_IN)
    ax.fill(np.r_[xo, xi[::-1]], np.r_[yo, yi[::-1]], fc=BODY_FC, ec="none", zorder=1)
    ax.plot(xo, yo, color=BODY_C, lw=1.6, zorder=2)
    ax.plot(xi, yi, color=BODY_C, lw=1.6, zorder=2)
    # the modelled quarter
    xq, yq = ellipse(le.A_OUT, le.B_OUT, 0, math.pi / 2, 100)
    xiq, yiq = ellipse(le.A_IN, le.B_IN, 0, math.pi / 2, 100)
    ax.fill(np.r_[xq, xiq[::-1]], np.r_[yq, yiq[::-1]], fc="#a9c6e8", ec="none", zorder=1)
    ax.text(2350, 1500, "modelled\nquarter", ha="center", va="center", fontsize=10, color=BODY_C, style="italic")
    for xs, ys in (([-3900, 3900], [0, 0]), ([0, 0], [-3400, 3400])):
        ax.plot(xs, ys, ls="-.", lw=0.9, color="#555", zorder=3)
    ax.text(3870, 90, "symmetry", fontsize=8.5, color="#555", ha="right", va="bottom")
    ax.text(70, 3380, "symmetry", fontsize=8.5, color="#555", ha="left", va="top")

    pressure_arrows(ax, np.linspace(0, 2 * math.pi, 28, endpoint=False), length=480)
    ax.text(-2900, 2850, f"uniform outward pressure\np = {le.P:g} MPa on the outer edge", color=LOAD_C,
            fontsize=11, ha="center", va="center")
    ax.text(-950, 430, "inner edge\ntraction free", fontsize=9.5, ha="center", va="center", color=DIM_C)

    # semi-axes
    dim_arrow(ax, (0, -150), (-le.A_IN, -150), f"{le.A_IN:g}", fs=10)
    dim_arrow(ax, (150, 0), (150, -le.B_IN), f"{le.B_IN:g}", fs=10, offset=(170, 0))
    dim_arrow(ax, (0, -le.B_OUT - 650), (-le.A_OUT, -le.B_OUT - 650), f"{le.A_OUT:g}", fs=10)
    dim_arrow(ax, (-le.A_OUT - 750, 0), (-le.A_OUT - 750, -le.B_OUT), f"{le.B_OUT:g}", fs=10, offset=(-260, 0))
    ax.plot([-le.A_OUT, -le.A_OUT], [-le.B_OUT - 750, -150], ls=":", color="#888", lw=0.8)
    ax.plot([-le.A_OUT - 850, -150], [-le.B_OUT, -le.B_OUT], ls=":", color="#888", lw=0.8)
    ax.text(-1000, -1700,
            r"inner ellipse  $(x/2000)^2 + (y/1000)^2 = 1$" + "\n" + r"outer ellipse  $(x/3250)^2 + (y/2750)^2 = 1$",
            fontsize=10, ha="center", va="center", bbox=dict(fc="white", ec="#bbb", pad=4), zorder=6)

    # point D and target
    ax.plot(le.A_IN, 0, marker="o", ms=8, color=PEAK_C, zorder=6)
    ax.annotate("D (2000, 0)\ntarget " + r"$\sigma_{yy}$" + f" = {le.TARGET_SYY_D:g} MPa",
                xy=(le.A_IN, 0), xytext=(1100, -1250), fontsize=11, color=PEAK_C, ha="center",
                arrowprops=dict(arrowstyle="->", color=PEAK_C), bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.text(2250, -3300, "plane stress, thickness T = 100 mm\n" + r"$E$ = " + f"{le.E / 1000:g} GPa, "
            + r"$\nu$ = " + f"{le.NU:g}", fontsize=10.5, ha="center", va="center", color=BODY_C, weight="bold")
    o = (-3900, 2000)
    for d, lab in (((450, 0), "x"), ((0, 450), "y")):
        arrow(ax, o, (o[0] + d[0], o[1] + d[1]), "k", ms=8, lw=1)
        ax.text(o[0] + 1.3 * d[0], o[1] + 1.3 * d[1], lab, fontsize=10, ha="center", va="center")
    ax.set_xlim(-4300, 4100)
    ax.set_ylim(-3700, 3500)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Benchmark 4: NAFEMS LE1, elliptic membrane (schematic, dimensions in mm)", fontsize=12.5)
    out = le.RUN_DIR / "nafems_le1_schematic.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")


# --------------------------------------------------------------------------- detailed figure

def panel_model(ax) -> None:
    c = le.build_case("CPS8", 2)
    xy = c["xy"]
    draw_mesh(ax, xy, c["elems"])
    ax.scatter(*zip(*xy.values()), s=7, color=BODY_C, zorder=3)

    # symmetry supports
    on_x0 = sorted((p for p in xy.values() if abs(p[0]) < 1e-9), key=lambda p: p[1])
    on_y0 = sorted((p for p in xy.values() if abs(p[1]) < 1e-9), key=lambda p: p[0])
    for x, y in on_x0:
        ax.plot(x - 70, y, marker="o", ms=5, mfc="white", mec=BC_C, zorder=4)
    for x, y in on_y0:
        ax.plot(x, y - 70, marker="o", ms=5, mfc="white", mec=BC_C, zorder=4)
    ax.text(-180, 1900, r"A-B: $u_x = 0$", rotation=90, color=BC_C, fontsize=10, ha="right", va="center")
    ax.text(2650, -200, r"D-C: $u_y = 0$", color=BC_C, fontsize=10, ha="center", va="top")

    # pressure, normal to the outer edge
    pressure_arrows(ax, np.linspace(0.02, math.pi / 2 - 0.02, 12), length=380, ms=7)
    ax.text(3250, 2650, f"p = {le.P:g} MPa\noutward, normal\nto the outer edge", color=LOAD_C, fontsize=10,
            ha="center", va="center")
    ax.text(640, 430, "inner edge free", fontsize=9, rotation=-28, ha="center", va="center", color=DIM_C)

    # named points
    offs = {"A": (-120, 0, "right"), "B": (-120, 60, "right"), "C": (60, -160, "left"), "D": (-60, -160, "right")}
    for k, (x, y) in le.POINTS.items():
        dx, dy, ha = offs[k]
        ax.plot(x, y, marker="s" if k != "D" else "o", ms=7 if k == "D" else 5,
                color=PEAK_C if k == "D" else "k", zorder=6)
        ax.text(x + dx, y + dy, f"{k} ({x:g}, {y:g})", fontsize=9.5, ha=ha, va="center",
                color=PEAK_C if k == "D" else "k", weight="bold" if k == "D" else None)
    ax.annotate(r"target $\sigma_{yy}(D)$" + f" = {le.TARGET_SYY_D:g} MPa", xy=(le.A_IN, 0), xytext=(1700, 900),
                fontsize=10.5, color=PEAK_C, arrowprops=dict(arrowstyle="->", color=PEAK_C),
                bbox=dict(fc="white", ec="none", alpha=0.85))

    # mesh construction at one parametric angle
    t = math.radians(50)
    pi, po = le.ellipse_point(le.A_IN, le.B_IN, t), le.ellipse_point(le.A_OUT, le.B_OUT, t)
    ax.plot([pi[0], po[0]], [pi[1], po[1]], color="#d69e2e", lw=2.2, zorder=4)
    for (xx, yy), lab, dx in ((pi, r"$(2000\cos t,\ 1000\sin t)$", -60), (po, r"$(3250\cos t,\ 2750\sin t)$", 120)):
        ax.plot(xx, yy, marker="o", ms=6, color="#d69e2e", zorder=6)
        ax.text(xx + dx, yy + (60 if dx > 0 else -60), lab, fontsize=9, color="#975a16",
                ha="left" if dx > 0 else "right", va="center", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    ax.text(-150, 3550, "mesh line: equal parametric angle t on both ellipses, divided evenly\n"
            "(level n: 2n elements across, 3n round; here level 2, CPS8)", fontsize=9, color="#975a16", va="top")

    ax.set_xlim(-650, 4000)
    ax.set_ylim(-500, 3650)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("(a) Quarter model as meshed (level 2: 24 elements, the NAFEMS fine-mesh count)", fontsize=10.5)


def panel_level(ax, level: int, letter: str) -> None:
    c = le.build_case("CPS4", level)
    draw_mesh(ax, c["xy"], c["elems"], lw=0.35 if level >= 8 else 0.6)
    ax.plot(le.A_IN, 0, marker="o", ms=5, color=PEAK_C, zorder=6)
    note = {1: " (NAFEMS coarse count)", 2: " (NAFEMS fine count)"}.get(level, "")
    ax.set_title(f"({letter}) level {level}: {len(c['elems'])} elements{note}", fontsize=9.5)
    ax.set_xlim(-100, 3350)
    ax.set_ylim(-100, 2850)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def detailed() -> None:
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 2, width_ratios=[2.4, 1], hspace=0.22, wspace=0.05)
    panel_model(fig.add_subplot(gs[:, 0]))
    for k, (lvl, letter) in enumerate(((1, "b"), (4, "c"), (8, "d"))):
        panel_level(fig.add_subplot(gs[k, 1]), lvl, letter)
    fig.suptitle("Benchmark 4: NAFEMS LE1 elliptic membrane, model and meshes", fontsize=13, y=0.95)
    out = le.RUN_DIR / "nafems_le1_geometry.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")


def main() -> None:
    schematic()
    detailed()


if __name__ == "__main__":
    main()
