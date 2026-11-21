"""Benchmark 3 figure: geometry, boundary conditions and variables of the Hertz model.

Draws the axisymmetric sphere-on-flat model exactly as hertz_contact.py builds it
(level-1 CAX4 mesh), with the variables labelled:

  (a) the whole model: sphere radius R, modelled region rho <= RHO_MAX = R round
      the contact point O, boundary conditions (axis u_r = 0, sphere's far
      boundary u_z = -delta, flat's far boundary u_z = 0), material E, nu
  (b) the contact zone: contact radius a, initial gap z = r^2 / 2R, the uniform
      mesh core RHO_C = 1.25 a and element size h
  (c) the Hertz pressure p(r) = p0 sqrt(1 - r^2/a^2) and the fit window used
      to measure p0 and a

No solver run is needed.

Usage
    uv run hertz_geometry.py        # writes runs/hertz/hertz_geometry.png and hertz_schematic.png
"""
from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Arc, FancyArrowPatch, Rectangle

import hertz_contact as hc

SPHERE_FC, FLAT_FC, EDGE = "#cfe0f3", "#f6dcc0", "#7a8796"
SPHERE_C, FLAT_C, BC_C, DIM_C = "#2b6cb0", "#c05621", "#b7312c", "#222222"


def hertz_design() -> tuple[float, float, float]:
    """Half-space Hertz values for the prescribed approach DELTA: a, F, p0."""
    a = math.sqrt(hc.DELTA * hc.R)
    F = 4 * hc.ESTAR * a**3 / (3 * hc.R)
    return a, F, 3 * F / (2 * math.pi * a * a)


def draw_mesh(ax, case) -> None:
    xy = case["xy"]
    polys, colours = [], []
    for c in case["elems"]:
        p = [xy[n] for n in c[:4]]
        polys.append(p)
        colours.append(SPHERE_FC if np.mean([q[1] for q in p]) > 0 else FLAT_FC)
    ax.add_collection(PolyCollection(polys, facecolors=colours, edgecolors=EDGE, linewidths=0.35, zorder=1))


def dim_arrow(ax, p, q, text, offset=(0, 0), colour=DIM_C, fs=10, **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>", mutation_scale=9, color=colour, lw=1, zorder=5))
    mid = ((p[0] + q[0]) / 2 + offset[0], (p[1] + q[1]) / 2 + offset[1])
    ax.text(*mid, text, color=colour, fontsize=fs, ha="center", va="center", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1), **kw)


def panel_model(ax, case) -> None:
    R, rho_max, delta = hc.R, hc.RHO_MAX, hc.DELTA
    xy = case["xy"]
    # context: the rest of the sphere and the half-space (not modelled)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(R * np.sin(th), R - R * np.cos(th), ls="--", lw=0.8, color="#9aa5b1", zorder=0)
    ax.plot([-1.25 * R, 1.25 * R], [0, 0], ls="--", lw=0.8, color="#9aa5b1", zorder=0)
    ax.fill_between([-1.25 * R, 0], -1.1 * R, 0, color="#f3f3f3", zorder=-1)
    ax.text(-0.62 * R, -0.55 * R, "mirror half\n(not modelled:\naxisymmetric)", ha="center", va="center",
            fontsize=8.5, color="#7b8794", style="italic")
    draw_mesh(ax, case)

    # axis of symmetry
    ax.plot([0, 0], [-1.08 * R, 2.05 * R], ls="-.", lw=1.1, color="#444", zorder=3)
    ax.text(0.15, 2.0 * R, "axis of symmetry", rotation=90, va="top", fontsize=8.5, color="#444")
    for z in np.linspace(-0.9 * R, 0.9 * R, 7):          # roller symbols: u_r = 0
        ax.plot([-0.22], [z], marker="o", ms=4.5, mfc="white", mec=BC_C, zorder=4)
    ax.text(-0.5, 0.42 * R, r"$u_r = 0$", color=BC_C, fontsize=10, ha="right", va="center")

    # sphere far boundary: prescribed u_z = -delta
    far_s = sorted((xy[n] for n in case["far_sphere"]), key=lambda p: p[0])
    for x, z in far_s[::2]:
        ax.add_patch(FancyArrowPatch((x, z + 1.7), (x, z + 0.25), arrowstyle="-|>", mutation_scale=8,
                                     color=BC_C, lw=1, zorder=4))
    ax.text(0.62 * R, 1.08 * R, r"$u_z = -\delta$" + f"\n$\\delta$ = {delta:g} mm", color=BC_C, fontsize=10,
            ha="left", va="bottom")

    # flat far boundary: u_z = 0 (hatched support)
    far_f = sorted((xy[n] for n in case["far_flat"]), key=lambda p: p[0])
    for x, z in far_f[::2]:
        ax.plot([x], [z - 0.3], marker="^", ms=6, mfc="white", mec=BC_C, zorder=4)
    ax.text(0.78 * R, -0.82 * R, r"$u_z = 0$", color=BC_C, fontsize=10, ha="left")

    # labels
    C = (0.0, R)
    ax.plot(*C, marker="+", ms=9, color=SPHERE_C, mew=1.5, zorder=5)
    ax.text(-0.45, R + 0.45, "C (sphere centre)", fontsize=8.5, color=SPHERE_C, ha="right")
    t = math.radians(-48)
    dim_arrow(ax, C, (R * math.sin(t), R - R * math.cos(t)), f"R = {R:g} mm", colour=SPHERE_C, offset=(0.2, 0.5))
    ax.plot(0, 0, marker="o", ms=4, color="k", zorder=6)
    ax.text(-0.5, -1.0, "O (initial contact point)", fontsize=8.5, ha="right", va="center")
    t = math.radians(62)
    dim_arrow(ax, (0, 0), (rho_max * math.cos(t), rho_max * math.sin(t)),
              r"$\rho_{max} = R$", offset=(0.75, -0.15), fs=9.5)
    t = math.radians(-35)
    dim_arrow(ax, (0, 0), (rho_max * math.cos(t), rho_max * math.sin(t)),
              r"$\rho_{max} = R$", offset=(0.5, 0.6), fs=9.5)
    lab = dict(ha="center", va="center", fontsize=10, weight="bold", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
    ax.text(0.24 * R, 0.8 * R, "sphere\n" + r"$E,\ \nu$", color=SPHERE_C, **lab)
    ax.text(0.33 * R, -0.72 * R, "flat (half-space)\n" + r"$E,\ \nu$", color=FLAT_C, **lab)
    ax.text(1.02 * R, 0.6, "frictionless contact\n(sphere = slave)", fontsize=8.5, ha="right", va="bottom")

    # coordinate triad
    o = (-1.1 * R, 1.6 * R)
    for d, lab in (((1.6, 0), "r"), ((0, 1.6), "z")):
        ax.add_patch(FancyArrowPatch(o, (o[0] + d[0], o[1] + d[1]), arrowstyle="-|>", mutation_scale=9, color="k"))
        ax.text(o[0] + 1.25 * d[0] + (0.15 if lab == "r" else 0), o[1] + 1.2 * d[1], lab, fontsize=10,
                ha="center", va="center")

    ax.set_xlim(-1.25 * R, 1.25 * R)
    ax.set_ylim(-1.12 * R, 2.08 * R)
    ax.set_aspect("equal")
    ax.set_xlabel("r (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("(a) Axisymmetric model and boundary conditions\n(level-1 CAX4 mesh, true scale; "
                 f"E = {hc.E / 1000:g} GPa, $\\nu$ = {hc.NU:g} for both bodies)", fontsize=10.5)
    ax.annotate("zoom (b)", xy=(-0.05, 0.3), xytext=(-5.5, 2.4), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", color="#555"))
    ax.add_patch(plt.Rectangle((-0.02, -0.5), 1.6, 1.0, fill=False, ec="#555", lw=0.8, ls=":", zorder=6))


def panel_contact(ax, case, a) -> None:
    R, rc, h = hc.R, hc.RHO_C, case["h"]
    draw_mesh(ax, case)
    r = np.linspace(0, 1.6, 300)
    ax.plot(r, R - np.sqrt(R * R - r * r), color=SPHERE_C, lw=1.6, zorder=3)
    ax.plot(r, 0 * r, color=FLAT_C, lw=1.6, zorder=3)
    ax.plot([0, 0], [-0.6, 0.6], ls="-.", lw=1.1, color="#444", zorder=3)

    # contact radius
    ax.axvline(a, ymin=0.1, ymax=0.9, ls=":", color="k", lw=1)
    dim_arrow(ax, (0, -0.07), (a, -0.07), "a", fs=11)
    # gap
    rg = 1.35
    gz = R - math.sqrt(R * R - rg * rg)
    ax.add_patch(FancyArrowPatch((rg, 0), (rg, gz), arrowstyle="<|-|>", mutation_scale=7, color=DIM_C, lw=1, zorder=5))
    ax.text(rg - 0.04, 0.17, r"initial gap" + "\n" + r"$z \approx r^2 / 2R$", fontsize=9, va="center", ha="right",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    # core and element size
    ax.add_patch(plt.Rectangle((0, -rc), rc, 2 * rc, fill=False, ec="#2f855a", lw=1.4, ls="--", zorder=4))
    ax.text(rc + 0.03, -0.3, r"uniform core $\rho_c = 1.25\,a$" + "\n(both bodies)", fontsize=9,
            color="#2f855a", va="center", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    dim_arrow(ax, (rc - h, -0.48), (rc, -0.48), "h", fs=10, offset=(0, -0.06))
    ax.text(0.04, 0.52, f"a/h = 4 x level (here level 1: a/h = {a / h:.0f})", fontsize=8.5, va="center",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    ax.text(-0.02, -0.03, "O", fontsize=9, ha="right", va="top")
    ax.set_xlim(-0.05, 1.6)
    ax.set_ylim(-0.6, 0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("r (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("(b) Contact zone (undeformed)", fontsize=10.5)


def panel_pressure(ax, a, p0) -> None:
    r = np.linspace(0, a, 300)
    p = p0 * np.sqrt(1 - (r / a) ** 2)
    lo, hi = hc.FIT_RANGE
    ax.axvspan(lo * a, hi * a, color="#e6f4ea", zorder=0)
    ax.text((lo + hi) / 2 * a, 0.12 * p0, f"fit window\n{lo:g} a < r < {hi:g} a", ha="center", fontsize=8.5,
            color="#2f855a")
    ax.fill_between(r, 0, p, color=SPHERE_FC, zorder=1)
    ax.plot(r, p, color=SPHERE_C, lw=2, zorder=2)
    ax.plot([a, 1.6 * a], [0, 0], color=SPHERE_C, lw=2)
    ax.annotate(r"$p_0$", xy=(0, p0), xytext=(0.12 * a, 1.12 * p0), fontsize=11,
                arrowprops=dict(arrowstyle="->", color="#555"))
    ax.annotate(r"$a$", xy=(a, 0), xytext=(1.08 * a, 0.3 * p0), fontsize=11,
                arrowprops=dict(arrowstyle="->", color="#555"))
    ax.text(0.04 * a, 0.45 * p0, r"$p(r) = p_0\sqrt{1 - r^2/a^2}$", fontsize=10.5, color=SPHERE_C)
    ax.text(1.56 * a, 1.36 * p0,
            "Hertz, identical materials:\n"
            r"$E^* = E\,/\,[2(1-\nu^2)]$" + "\n"
            r"$a = (3FR\,/\,4E^*)^{1/3}$" + "\n"
            r"$p_0 = 3F\,/\,(2\pi a^2)$" + "\n"
            "F = reaction on the\nsphere's far boundary",
            fontsize=9.5, ha="right", va="top", linespacing=1.5,
            bbox=dict(fc="white", ec="#bbb", pad=4))
    ax.set_xlim(0, 1.6 * a)
    ax.set_ylim(0, 1.42 * p0)
    ax.set_xlabel("r (mm)")
    ax.set_ylabel("contact pressure p (MPa)")
    ax.set_title(f"(c) Hertz pressure (design values: a = {a:.2f} mm, p0 = {p0:.0f} MPa)", fontsize=10.5)


def schematic() -> None:
    """Simple problem schematic: whole sphere on a flat, no mesh, key variables only.

    The contact radius is drawn much larger than it is (a/R = 0.05 in the model)."""
    fig, ax = plt.subplots(figsize=(9, 8.2))
    R, ad = 10.0, 2.6                         # drawn contact radius (exaggerated)
    zc = math.sqrt(R * R - ad * ad)           # sphere centre height: sphere flattened over r < a
    th = np.linspace(0, 2 * np.pi, 400)
    xs, ys = R * np.sin(th), zc - R * np.cos(th)
    ys = np.maximum(ys, 0.0)
    ax.fill(xs, ys, fc=SPHERE_FC, ec=SPHERE_C, lw=1.6, zorder=2)
    W, D = 15.0, 9.0
    ax.add_patch(Rectangle((-W, -D), 2 * W, D, fc=FLAT_FC, ec=FLAT_C, lw=1.6, zorder=1))
    # support under the flat
    ax.plot([-W, W], [-D, -D], color="k", lw=2, zorder=3)
    for x in np.linspace(-W, W, 24):
        ax.plot([x, x - 0.7], [-D, -D - 0.7], color="k", lw=0.8)
    ax.text(W - 0.5, -D - 1.1, "far boundary held: " + r"$u_z = 0$", color=BC_C, fontsize=10.5, ha="right",
            va="top")
    # load on the sphere
    arrow_kw = dict(arrowstyle="-|>", mutation_scale=16, color=BC_C, lw=2.4, zorder=5)
    ax.add_patch(FancyArrowPatch((0, zc + R + 5.5), (0, zc + R + 0.3), **arrow_kw))
    ax.text(0.6, zc + R + 3.2, "pressed down by approach " + r"$\delta$" + "\n(contact force F)", color=BC_C,
            fontsize=11, va="center")
    # axis of symmetry
    ax.plot([0, 0], [-D - 0.3, zc + R + 0.2], ls="-.", lw=1, color="#444", zorder=4)
    ax.text(-0.35, zc + R - 1.2, "axis of\nsymmetry", fontsize=9, color="#444", ha="right", va="top")
    # sphere radius and centre
    ax.plot(0, zc, marker="+", ms=10, mew=1.6, color=SPHERE_C, zorder=5)
    t = math.radians(130)
    dim_arrow(ax, (0, zc), (R * math.cos(t), zc + R * math.sin(t)), f"R = {hc.R:g} mm", colour=SPHERE_C,
              offset=(0.3, 0.5), fs=11)
    ax.text(3.2, zc - 1.0, "sphere\n" + r"$E,\ \nu$", color=SPHERE_C, fontsize=12, weight="bold", ha="center")
    ax.text(-9.0, -4.5, "flat (half-space)\n" + r"$E,\ \nu$", color=FLAT_C, fontsize=12, weight="bold",
            ha="center", va="center")
    # contact pressure on the flat
    r = np.linspace(-ad, ad, 200)
    ph_ = 3.2 * np.sqrt(np.maximum(0, 1 - (r / ad) ** 2))
    ax.fill_between(r, -ph_, 0, color="#f1c9a0", alpha=0.9, zorder=2)
    ax.plot(r, -ph_, color=FLAT_C, lw=1.4, zorder=3)
    for x in np.linspace(-0.85 * ad, 0.85 * ad, 7):
        h = 3.2 * math.sqrt(1 - (x / ad) ** 2)
        ax.add_patch(FancyArrowPatch((x, -0.05), (x, -h + 0.1), arrowstyle="-|>", mutation_scale=7, color=FLAT_C,
                                     lw=0.9, zorder=4))
    ax.annotate(r"$p_0$", xy=(0.05, -3.2), xytext=(4.5, -4.6), fontsize=12, color=FLAT_C,
                arrowprops=dict(arrowstyle="->", color=FLAT_C))
    ax.text(4.0, -6.6, r"$p(r) = p_0\sqrt{1 - r^2/a^2}$", fontsize=11, color=FLAT_C)
    dim_arrow(ax, (0, 0.6), (ad, 0.6), "a", fs=12, offset=(0, 0.9))
    ax.text(W - 0.3, 0.4, "frictionless contact", fontsize=10, ha="right", va="bottom")
    ax.text(0, -D - 3.6,
            r"Hertz: $E^* = E/[2(1-\nu^2)]$,  $a = (3FR/4E^*)^{1/3}$,  $p_0 = 3F/(2\pi a^2)$" + "\n"
            + f"E = {hc.E / 1000:g} GPa, " + r"$\nu$ = " + f"{hc.NU:g} (steel, both bodies);  "
            + f"model: a/R = {hc.A_TARGET / hc.R:g}, " + r"$\delta$" + f" = {hc.DELTA:g} mm",
            ha="center", va="top", fontsize=10.5, bbox=dict(fc="white", ec="#bbb", pad=5))
    o = (-W + 0.5, zc + R - 1)
    for d, lab in (((2, 0), "r"), ((0, 2), "z")):
        ax.add_patch(FancyArrowPatch(o, (o[0] + d[0], o[1] + d[1]), arrowstyle="-|>", mutation_scale=9, color="k"))
        ax.text(o[0] + 1.3 * d[0], o[1] + 1.3 * d[1], lab, fontsize=10, ha="center", va="center")
    ax.set_xlim(-W - 1, W + 1)
    ax.set_ylim(-D - 7.5, zc + R + 6.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Benchmark 3: Hertz contact, sphere on flat\n(schematic; contact radius exaggerated)", fontsize=12.5)
    out = hc.RUN_DIR / "hertz_schematic.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")


def main() -> None:
    case = hc.build_case("CAX4", 1)
    a, F, p0 = hertz_design()
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], hspace=0.32, wspace=0.18)
    panel_model(fig.add_subplot(gs[:, 0]), case)
    panel_contact(fig.add_subplot(gs[0, 1]), case, a)
    panel_pressure(fig.add_subplot(gs[1, 1]), a, p0)
    fig.suptitle("Benchmark 3: Hertz contact, elastic sphere on elastic flat", fontsize=13, y=0.98)
    out = hc.RUN_DIR / "hertz_geometry.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure: {out}")
    schematic()


if __name__ == "__main__":
    main()
