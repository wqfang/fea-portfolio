# FEA verification suite

Benchmarks with known answers, each run in **CalculiX** (ccx, open source) and
**Abaqus/Standard** from the same Abaqus-format input deck, with mesh-convergence
studies. Units: mm, N, MPa.

| # | Benchmark | Checks | Status |
|---|-----------|--------|--------|
| 0 | Single-element uniaxial tension (`tension.inp`) | Material and boundary-condition setup; ccx vs Abaqus | Done: both match closed form |
| 1 | End-loaded cantilever (`cantilever.py`) | Bending, shear locking, hourglassing, convergence rate | Done: CalculiX and Abaqus |
| 2 | Plate with a hole (`plate_hole.py`) | Stress concentration factor Kt, stress convergence | Done |
| 3 | Hertz sphere-on-flat contact | Contact pressure and radius | Planned |
| 4 | NAFEMS LE1 elliptic membrane | Standard benchmark from NAFEMS (National Agency for Finite Element Methods and Standards) | Planned |

## Setup

Python dependencies are managed with [uv](https://docs.astral.sh/uv/) from
`pyproject.toml`; `uv run` creates `.venv/` and installs them on first use.

```bash
brew install costerwi/homebrew-calculix/calculix-ccx
uv sync          # optional: uv run does this automatically
```

Homebrew installs CalculiX as a versioned binary (e.g. `ccx_2.23`) with no plain
`ccx`. The scripts find the newest `ccx_<version>` on `PATH` automatically; a
shell alias is not enough because Python cannot see aliases. Pass `--ccx <path>`
to use a specific binary.

The scripts run CalculiX on **one thread**. With ccx 2.23 the multithreaded
SPOOLES solver intermittently returned wrong answers with no error (the same deck
gave 0.06 mm on one 8-thread run and 1.403 mm on others), and these models are
too small to benefit from threads.

## Benchmark 1: cantilever beam

Steel beam, 200 x 10 x 10 mm, clamped at x = 0, uniform shear traction totalling
P = 100 N on the free end. Five element types (C3D8, C3D8R, C3D8I, C3D20, C3D20R)
at four mesh levels (1, 2, 4, 8 elements through the depth; cube-shaped elements).

```bash
uv run cantilever.py                    # full sweep, ~1 minute on a laptop
uv run cantilever.py --levels 1 2 4     # quicker
```

Outputs go to `runs/cantilever/`: `results.csv`, `cantilever_convergence.png`,
and the decks.

**Tip deflection measure.** The script reports the load-weighted tip deflection,
sum(f_i * v_i) / P, which is the displacement conjugate to the applied load. It
is the quantity beam theory predicts, and it is not biased by local distortion of
the loaded face.

**References.**
- Timoshenko beam theory (bending + shear, Cowper shear coefficient): 1.6031 mm.
- Richardson extrapolation of the three finest C3D20 meshes (n = 2, 4, 8):
  **1.5976 mm**, the converged answer of the 3D model. Three meshes give three
  equations for the unknowns in d(h) = d_exact + C h^p, so the script solves for
  d_exact and the observed convergence order p (1.56 here). This separates mesh
  error from modelling error: the converged 3D value sits 0.34% below Timoshenko
  because the fully clamped 3D root also blocks Poisson contraction, which beam
  theory ignores. The observed order is well below the h^4 expected for quadratic
  elements; see *Why the convergence order is 1.6* below.

## Findings, Benchmark 2: Abaqus cross-check (Learning Edition 2026)

Every case under the 1,000-node cap was run in Abaqus: levels 1, 2 and 4 for
the linear elements, and 1 and 2 for the quadratic ones (26 decks). All ran
without errors. The only warning is a generic note about section thickness in
contact, which does not apply here. Abaqus prints the averaged nodal stresses
(`*EL PRINT, POSITION=AVERAGED AT NODES`) to 4 significant figures, so its
peak stresses are resolved to 0.1 MPa (0.03%).

To rerun: copy `runs/plate_hole/abaqus/` to the Windows VM, run `run_all.bat`
in an Abaqus Command window in that folder, copy the `.dat` files back, then
run `uv run plate_hole.py --collect-abaqus`. The Abaqus points appear as
crosses on `plate_hole_convergence.png`.

**Kirsch problem**: peak stress error vs the exact 300 MPa, Abaqus (CalculiX in
brackets).

| Element | n=1 | n=2 | n=4 |
|---|---|---|---|
| CPS4 | -2.2% (-1.7%) | +0.43% (+0.96%) | +0.57% (+0.86%) |
| CPS4R | -19% (-65%) | -8.7% (-61%) | -4.0% (-44%) |
| CPS4I | -4.2% (n/a) | -0.63% (n/a) | +0.067% (n/a) |
| CPS8 | -0.53% (-0.65%) | -0.20% (-0.21%) | |
| CPS8R | +1.7% (+1.7%) | +0.47% (+0.48%) | |

RMS hoop-stress error round the hole (MPa), Abaqus (CalculiX in brackets):

| Element | n=1 | n=2 | n=4 |
|---|---|---|---|
| CPS4 | 8.5 (9.2) | 1.36 (1.94) | 0.98 (1.49) |
| CPS4R | 38 (119) | 17.6 (108) | 8.3 (78) |
| CPS4I | 10.6 | 2.12 | 0.47 |
| CPS8 | 1.39 (1.65) | 0.39 (0.42) | |
| CPS8R | 2.43 (2.43) | 0.62 (0.62) | |

**Finite plate**: Kt_net, Abaqus (converged FE value 2.522; Peterson 2.508).

| Element | n=1 | n=2 | n=4 |
|---|---|---|---|
| CPS4 | 2.446 | 2.526 | 2.534 |
| CPS4R | 2.040 | 2.299 | 2.418 |
| CPS4I | 2.410 | 2.504 | 2.523 |
| CPS8 | 2.506 | 2.517 | |
| CPS8R | 2.566 | 2.534 | |

- **The quadratic elements agree between the two codes.** CPS8R matches
  CalculiX to the 4 printed digits (305.1 vs 305.13 MPa at n=1). CPS8 agrees
  within 0.12% at n=1 and 0.01% at n=2, and the gap shrinks with the mesh; it
  is consistent with the small difference between CalculiX's thin 3D brick and
  Abaqus's true plane stress, plus slightly different stress extrapolation to
  the nodes. This independently confirms the CalculiX results, the converged
  Kt_net = 2.522, and its 0.55% gap from Peterson's curve fit.
- **CPS4 differs by about 0.3-0.5% at the peak, and Abaqus is closer to the
  exact answer** (+0.57% vs +0.86% at n=4; hoop error 0.98 vs 1.49 MPa). This is
  most likely the formulation difference found for C3D8 in Benchmark 1 (not
  re-derived here with a numpy solver): Abaqus's
  "fully integrated" first-order elements use selectively reduced integration
  (B-bar) of the volumetric part, while CalculiX expands CPS4 into a plain
  fully integrated C3D8.
- **Abaqus CPS4R behaves like a normal reduced-integration element**: 19%,
  8.7% and 4.0% low, halving with each mesh refinement (first order). The error
  comes from sampling stress at one point per element, the centre, and
  extrapolating to the hole edge. This confirms that CalculiX's 44-65% errors
  come from its thin expanded brick, not from CPS4R. Even in Abaqus, CPS4R
  needs a fine mesh for peak stresses.
- **CPS4I is the best first-order element.** Its incompatible (bubble) modes let
  the element bend without parasitic shear. It is the least accurate on the
  coarsest mesh, but it converges fastest: 0.067% at n=4, with a hoop error
  (0.47 MPa) half that of CPS4 and close to CPS8 at n=2. On the plate, n=4 gives
  Kt_net = 2.523, within 0.04% of the converged value. CalculiX 2.23 has no
  CPS4I.
- **Practical ranking for peak stress at a notch:** CPS8 ≈ CPS8R > CPS4I >
  CPS4 >> CPS4R, in both codes. In CalculiX, avoid CPS4R for 2D work entirely.

