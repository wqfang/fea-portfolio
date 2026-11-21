"""Benchmark 3: Hertz contact, elastic sphere pressed onto an elastic flat.

Axisymmetric, frictionless, linear elastic, both bodies steel. The sphere (radius
R) touches the flat at one point; its far boundary is pushed down by DELTA and the
flat's far boundary is held. Each body is the part of the sphere / half-space
within RHO_MAX of the initial contact point, meshed as an O-grid: a uniform square
core round the contact (RHO_C = 1.25 a) and a graded ring out to RHO_MAX = R.

Hertz, for the contact force F measured in the model:
    E* = E / (2 (1 - nu^2))      (two identical materials)
    a  = (3 F R / (4 E*))^(1/3)  contact radius
    p0 = 3 F / (2 pi a^2)        peak contact pressure
    p(r) = p0 sqrt(1 - r^2 / a^2)

Measured quantities: F from the reaction on the sphere's far boundary (checked
against the integral of the contact pressure, F_int_ratio); the contact pressure
CPRESS at the sphere's surface nodes. p0 and a are found by fitting the
Hertz ellipse, p^2 = p0^2 (1 - r^2 / a^2), to the nodes with 0.25 a < r < 0.9 a.
The fit leaves out the nodes next to the axis, where CalculiX's contact pressure
is unreliable (see README), and the edge, where a node-based pressure cannot
follow the sharp drop to zero.

Usage
    uv run hertz_contact.py                     # CalculiX sweep + plots (~5 minutes)
    uv run hertz_contact.py --levels 1 2 4      # quicker
    uv run hertz_contact.py --collect-abaqus    # also read Abaqus .dat results

Abaqus Learning Edition caps models at 1,000 nodes. Every case under the cap is
written to runs/hertz/abaqus/ with run_all.bat (including CAX4I, which ccx does not
support); copy that folder to the Windows VM, run run_all.bat, copy the .dat files
back, then use --collect-abaqus.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np

from fea_common import _id_lines, run_ccx
from frd2vtu import read_frd
from plate_hole import Registry, block, edges_along

# Geometry (mm), material (MPa)
R = 10.0                      # sphere radius
E, NU = 200_000.0, 0.3
ESTAR = E / (2 * (1 - NU**2))
A_TARGET = 0.5                # design contact radius, a/R = 0.05
DELTA = A_TARGET**2 / R       # prescribed approach of the far boundaries (Hertz half-space value)
RHO_C = 1.25 * A_TARGET       # uniform core round the contact point
RHO_MAX = R                   # far boundary of both bodies, 20 a from the contact
K_PEN = 1e8                   # penalty contact stiffness (MPa/mm): overclosure ~ p0 / K = 0.14% of DELTA
CCX_WEDGE_DEG = 2.0           # CalculiX expands CAX elements into a 2 degree wedge; forces are for that wedge
FIT_RANGE = (0.25, 0.9)       # r / a range used for the Hertz fit
ABAQUS_LE_NODE_CAP = 1000

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE / "runs" / "hertz"
ELEMENT_ORDER = {"CAX4": 1, "CAX4R": 1, "CAX4I": 1, "CAX8": 2, "CAX8R": 2}
CCX_ELEMENTS = ["CAX4", "CAX4R", "CAX8", "CAX8R"]
ELEMS_PER_LEVEL = 5           # core is 5 * level elements across: a/h = 4 * level


# --------------------------------------------------------------------------- reference

def hertz(F: float) -> tuple[float, float]:
    """Contact radius and peak pressure for total force F (sphere radius R, identical materials)."""
    a = (3 * F * R / (4 * ESTAR)) ** (1 / 3)
    return a, 3 * F / (2 * math.pi * a * a)


# --------------------------------------------------------------------------- mesh

def _ring_grading(n: int, first: float, total: float):
    """Geometric spacing of n steps: first step `first`, total length `total`."""
    lo, hi = 1.0 + 1e-12, 5.0
    for _ in range(200):
        g = 0.5 * (lo + hi)
        if first * (g**n - 1) / (g - 1) > total:
            hi = g
        else:
            lo = g
    return lambda s: (g ** (n * s) - 1) / (g**n - 1)


def _body(reg: Registry, kind: str, m: int, order: int):
    """One body as an O-grid in reference coordinates (u along the surface, v into the
    body), mapped to the physical shape. Returns elements, surface edges, far-boundary nodes."""
    h = RHO_C / m
    w = _ring_grading(m, h, RHO_MAX - RHO_C)

    def phys(u, v):
        rho, psi = math.hypot(u, v), math.atan2(v, u)
        if kind == "flat":
            return rho * math.cos(psi), -rho * math.sin(psi)
        # sphere: a point at distance rho from the contact point lies at angle asin(rho/2R)
        # above the surface tangent; squeeze the reference angles into [that, 90 deg]
        phs = math.asin(min(rho / (2 * R), 1.0))
        phi = phs + psi * (math.pi / 2 - phs) / (math.pi / 2)
        return rho * math.cos(phi), rho * math.sin(phi)

    def ring(edge, angle):
        def pt(s, t):
            ex, ey = edge(t)
            ax, ay = RHO_MAX * math.cos(angle(t)), RHO_MAX * math.sin(angle(t))
            f = w(s)
            return phys(ex + f * (ax - ex), ey + f * (ay - ey))
        return pt

    blocks = [("core", lambda s, t: phys(s * RHO_C, t * RHO_C)),
              ("ring_a", ring(lambda t: (RHO_C, t * RHO_C), lambda t: t * math.pi / 4)),
              ("ring_b", ring(lambda t: (RHO_C * (1 - t), RHO_C), lambda t: math.pi / 4 * (1 + t)))]
    st = 2 if order == 2 else 1
    elems, surf, far = [], [], set()
    for name, pt in blocks:
        e, grid = block(reg, pt, m, m, order)
        elems += e
        if name in ("core", "ring_a"):
            surf += edges_along(grid, None, 0, m, order)
        if name != "core":
            far |= {nid for (I, _), nid in grid.items() if I == st * m}
    if kind == "flat":  # the mirror map makes elements clockwise: reverse them
        elems = [[c[0], c[3], c[2], c[1]] + ([c[7], c[6], c[5], c[4]] if order == 2 else []) for c in elems]
    return elems, surf, sorted(far)


def _faces(elems, edges, first_id: int) -> list[tuple[int, int]]:
    """(element id, face number) for each edge, matched on its corner nodes."""
    by_corners = {}
    for k, c in enumerate(elems):
        for f in range(4):
            by_corners[frozenset((c[f], c[(f + 1) % 4]))] = (k + first_id, f + 1)
    return [by_corners[frozenset(e[:2])] for e in edges]


def build_case(etype: str, level: int):
    m, order = ELEMS_PER_LEVEL * level, ELEMENT_ORDER[etype]
    rs, rf = Registry(), Registry()
    es, ss, far_s = _body(rs, "sphere", m, order)
    ef, sf, far_f = _body(rf, "flat", m, order)
    off = len(rs.xy)
    xy = dict(rs.xy) | {n + off: p for n, p in rf.xy.items()}
    ef = [[n + off for n in c] for c in ef]
    sf = [[n + off for n in e] for e in sf]
    return {
        "xy": xy, "elems": es + ef,
        "surf_sphere": _faces(es, ss, 1), "surf_flat": _faces(ef, sf, len(es) + 1),
        "sphere_surface": {n for e in ss for n in e},
        "far_sphere": far_s, "far_flat": [n + off for n in far_f],
        "axis": [n for n, (x, _) in xy.items() if abs(x) < 1e-9],
        "h": RHO_C / m,
    }


# --------------------------------------------------------------------------- deck

def write_deck(path: Path, case: dict, etype: str, solver: str, k_pen: float = K_PEN) -> None:
    xy = case["xy"]
    o = [f"*HEADING\n{path.stem}: Hertz sphere on flat, axisymmetric", "*NODE, NSET=NALL"]
    o += [f"{n}, {x:.12g}, {y:.12g}" for n, (x, y) in xy.items()]
    o.append(f"*ELEMENT, TYPE={etype}, ELSET=EALL")
    o += [", ".join(map(str, [k, *c])) for k, c in enumerate(case["elems"], start=1)]
    for name, ids in (("NTOP", case["far_sphere"]), ("NBOT", case["far_flat"]), ("NAXIS", case["axis"])):
        o += [f"*NSET, NSET={name}", _id_lines(ids)]
    o += ["*SURFACE, NAME=SSPHERE, TYPE=ELEMENT"] + [f"{e}, S{f}" for e, f in case["surf_sphere"]]
    o += ["*SURFACE, NAME=SFLAT, TYPE=ELEMENT"] + [f"{e}, S{f}" for e, f in case["surf_flat"]]
    o += ["*MATERIAL, NAME=STEEL", "*ELASTIC", f"{E:.10g}, {NU:.10g}",
          "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL",
          "*SURFACE INTERACTION, NAME=FRICTIONLESS",
          "*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR", f"{k_pen:.6g},",
          "*CONTACT PAIR, INTERACTION=FRICTIONLESS, TYPE=SURFACE TO SURFACE", "SSPHERE, SFLAT",
          "*BOUNDARY", "NAXIS, 1, 1, 0.", "NBOT, 2, 2, 0.",
          "*STEP", "*STATIC", "0.25, 1.",
          "*BOUNDARY", f"NTOP, 2, 2, {-DELTA:.10g}",
          # ccx accepts TOTALS=ONLY; Abaqus allows only YES/NO (SUMMARY=NO drops the max/min table)
          "*NODE PRINT, NSET=NTOP, " + ("TOTALS=ONLY" if solver == "ccx" else "TOTALS=YES, SUMMARY=NO"), "RF"]
    if solver == "ccx":
        o += ["*NODE FILE", "U, S", "*CONTACT FILE", "CDIS, CSTR"]
    elif solver == "abaqus":
        o += ["*CONTACT PRINT, SLAVE=SSPHERE", "CPRESS, COPEN"]
    else:
        raise ValueError(solver)
    o += ["*END STEP", ""]
    path.write_text("\n".join(o))


# --------------------------------------------------------------------------- results

_NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"


def ccx_results(deck: Path) -> tuple[float, dict[int, float]]:
    """Total force on the sphere's far boundary (whole 360 deg) and CPRESS at slave nodes."""
    dat = deck.with_suffix(".dat").read_text()
    blocks = re.findall(r"total force \(fx,fy,fz\) for set NTOP.*?\n\s*\n\s*(" + _NUM + r")\s+(" + _NUM + ")", dat)
    F = -float(blocks[-1][1]) * 360.0 / CCX_WEDGE_DEG
    contact = [d for name, _, d, *_ in read_frd(deck.with_suffix(".frd"))[2] if name == "CONTACT"][-1]
    return F, {n: v[3] for n, v in contact.items()}


def abaqus_results(dat_path: Path) -> tuple[float, dict[int, float]]:
    """Total RF2 on NTOP (Abaqus reports axisymmetric forces for the full 360 deg) and CPRESS."""
    text = dat_path.read_text(errors="replace")
    totals = re.findall(r"^\s*TOTAL\s+(" + _NUM + r")\s+(" + _NUM + ")", text, re.M)
    if not totals:
        raise ValueError(f"no RF total in {dat_path}")
    F = -float(totals[-1][1].replace("D", "E"))
    # last CPRESS table (final increment). Rows: slave node, status (CL closed / OP open),
    # CPRESS, COPEN; Abaqus prints only the slave nodes near the contact, and writes
    # numbers with a trailing point ("3427.")
    p = {}
    tail = text[text.rindex("CPRESS"):]
    num = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
    for line in tail.splitlines()[1:]:
        m = re.match(r"^\s*(\d+)\s+(?:[A-Z]{2}\s+)?(" + num + r")\s+(" + num + r")\s*$", line)
        if m:
            p[int(m.group(1))] = float(m.group(2).replace("D", "E"))
        elif p and not line.strip():
            break
    if not p:
        raise ValueError(f"no CPRESS table in {dat_path}")
    return F, p


def metrics(case: dict, F: float, cpress: dict[int, float]) -> dict:
    xy = case["xy"]
    aH, p0H = hertz(F)
    pts = sorted((xy[n][0], p) for n, p in cpress.items() if n in case["sphere_surface"])
    r = np.array([q[0] for q in pts])
    p = np.array([q[1] for q in pts])
    lo, hi = FIT_RANGE
    sel = (r > lo * aH) & (r < hi * aH)
    slope, icpt = np.polyfit(r[sel] ** 2, p[sel] ** 2, 1) if sel.sum() >= 3 else (1.0, -1.0)
    if icpt > 0 and slope < 0:
        p0_fit, a_fit = math.sqrt(icpt), math.sqrt(-icpt / slope)
    else:  # too few nodes in contact to fit an ellipse (coarsest meshes)
        p0_fit = a_fit = float("nan")
    p_h = p0H * np.sqrt(np.clip(1 - (r / aH) ** 2, 0, None))
    axis = p[np.argmin(r)]
    in_contact = r[p > 1e-6 * p0H]
    F_int = float(np.trapezoid(p * 2 * np.pi * r, r))  # checks the solver's axisymmetric force convention
    return {"F": F, "F_int_ratio": F_int / F, "a_hertz": aH, "p0_hertz": p0H, "p0_fit": p0_fit, "a_fit": a_fit,
            "p0_err": p0_fit / p0H - 1, "a_err": a_fit / aH - 1,
            "p_rms_err": float(np.sqrt(np.mean(((p - p_h)[sel] / p0H) ** 2))),
            "p_axis_ratio": axis / p0H, "p_max_ratio": p.max() / p0H,
            "r_last_contact_ratio": (in_contact.max() / aH) if in_contact.size else float("nan"),
            "profile": list(zip((r / aH).tolist(), (p / p0H).tolist()))}


# --------------------------------------------------------------------------- driver

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--levels", type=int, nargs="+", default=[1, 2, 4, 8])
    ap.add_argument("--collect-abaqus", action="store_true")
    args = ap.parse_args()

    ccx_dir, abq_dir = RUN_DIR / "ccx", RUN_DIR / "abaqus"
    ccx_dir.mkdir(parents=True, exist_ok=True)
    abq_dir.mkdir(parents=True, exist_ok=True)
    rows, profiles, abq_jobs = [], {}, []
    for etype in ELEMENT_ORDER:
        for lvl in args.levels:
            case = build_case(etype, lvl)
            job = f"hertz_{etype}_n{lvl}"
            nn = len(case["xy"])
            base = {"element": etype, "level": lvl, "nodes": nn, "h_mm": case["h"],
                    "a_over_h": A_TARGET / case["h"]}
            if nn < ABAQUS_LE_NODE_CAP:
                write_deck(abq_dir / f"{job}.inp", case, etype, "abaqus")
                abq_jobs.append(job)
                dat = abq_dir / f"{job}.dat"
                if args.collect_abaqus and dat.exists():
                    m = metrics(case, *abaqus_results(dat))
                    profiles[(etype, lvl, "abaqus")] = m.pop("profile")
                    rows.append({**base, "solver": "abaqus", **m})
            if etype not in CCX_ELEMENTS:
                continue
            deck = ccx_dir / f"{job}.inp"
            write_deck(deck, case, etype, "ccx")
            try:
                run_ccx(deck)
                m = metrics(case, *ccx_results(deck))
            except (RuntimeError, ValueError, IndexError) as exc:
                print(f"{etype:5s} n={lvl}  FAILED: {exc}", flush=True)
                rows.append({**base, "solver": "ccx", "F": "failed"})
                continue
            profiles[(etype, lvl, "ccx")] = m.pop("profile")
            rows.append({**base, "solver": "ccx", **m})
            print(f"{etype:5s} n={lvl}  nodes={nn:6d}  F={m['F']:8.1f}  p0 err={m['p0_err']:+.3%}"
                  f"  a err={m['a_err']:+.3%}  rms={m['p_rms_err']:.3%}  p(axis)/p0={m['p_axis_ratio']:.2f}  pmax/p0={m['p_max_ratio']:.2f}",
                  flush=True)

    bat = ["@echo off", "rem Run every Hertz contact case that fits the Learning Edition node cap"]
    bat += [f"call abaqus job={j} interactive" for j in abq_jobs]
    (abq_dir / "run_all.bat").write_text("\r\n".join(bat) + "\r\n")
    print(f"Abaqus decks under {ABAQUS_LE_NODE_CAP} nodes: {len(abq_jobs)} -> {abq_dir}")

    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(RUN_DIR / "results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    make_plots(rows, profiles)


def make_plots(rows, profiles) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(ELEMENT_ORDER, ["C0", "C1", "C4", "C2", "C3"]))
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    rr = np.linspace(0, 1, 200)
    ax.plot(rr, np.sqrt(1 - rr**2), "k-", lw=2, label="Hertz")
    for et in ELEMENT_ORDER:
        lv = [k[1] for k in profiles if k[0] == et and k[2] == "ccx"]
        if lv:
            r, p = zip(*profiles[(et, max(lv), "ccx")])
            ax.plot(r, p, ".", ms=4, color=colors[et], label=f"{et} ccx n={max(lv)}")
    ax.axvspan(0, FIT_RANGE[0], color="0.9", zorder=0)
    pmax = [max(q[1] for q in v) for k, v in profiles.items() if k[2] == "ccx"]
    if pmax:
        ax.text(0.01, 1.52, f"off scale: node next to the axis\nup to {max(pmax):.1f} p0; axis node 0",
                fontsize=8, va="top")
    ax.set_xlim(0, 1.4)
    ax.set_ylim(0, 1.6)
    ax.set_xlabel("r / a (Hertz)")
    ax.set_ylabel("contact pressure / p0 (Hertz)")
    ax.set_title("Contact pressure, finest CalculiX mesh\n(grey: near-axis zone left out of the fit)")
    ax.legend(fontsize=8, frameon=False)

    for ax, key, title in ((axes[1], "p0_err", "peak pressure p0 (fit)"), (axes[2], "a_err", "contact radius a (fit)")):
        for et in ELEMENT_ORDER:
            for solver, mk, ls in (("ccx", "o", "-"), ("abaqus", "x", "--")):
                rs = sorted((r for r in rows if r["element"] == et and r["solver"] == solver
                             and r.get("F") not in ("failed", None)), key=lambda r: r["a_over_h"])
                if rs:
                    ax.loglog([r["a_over_h"] for r in rs], [max(abs(r[key]), 1e-5) for r in rs], marker=mk,
                              ls=ls, color=colors[et], label=et if solver == "ccx" or et == "CAX4I" else None)
        ax.set_xlabel("elements across the contact radius, a / h")
        ax.set_ylabel(f"|{title.split(' (')[0]} error| vs Hertz")
        ax.set_title(f"{title} (o = CalculiX, x = Abaqus)")
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(RUN_DIR / "hertz_convergence.png", dpi=130)
    print(f"Plot: {RUN_DIR / 'hertz_convergence.png'}")


if __name__ == "__main__":
    main()
