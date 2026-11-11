# FEA verification suite

Benchmarks with known answers, each run in **CalculiX** (ccx, open source) and
**Abaqus/Standard** from the same Abaqus-format input deck, with mesh-convergence
studies. Units: mm, N, MPa.

| # | Benchmark | Checks | Status |
|---|-----------|--------|--------|
| 0 | Single-element uniaxial tension (`tension.inp`) | Material and boundary-condition setup; ccx vs Abaqus | Done: both match closed form |
| 1 | End-loaded cantilever (`cantilever.py`) | Bending, shear locking, hourglassing, convergence rate | CalculiX done; Abaqus pending |
| 2 | Plate with a hole | Stress concentration factor Kt | Planned |
| 3 | Hertz sphere-on-flat contact | Contact pressure and radius | Planned |
| 4 | NAFEMS LE1 elliptic membrane | Standard benchmark from NAFEMS (National Agency for Finite Element Methods and Standards) | Planned |

## Setup

Python dependencies are managed with [uv](https://docs.astral.sh/uv/) from
`pyproject.toml`; `uv run` creates `.venv/` and installs them on first use.

```bash
brew install costerwi/homebrew-calculix/calculix-ccx
uv sync          # optional: uv run does this automatically
```

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
- Richardson extrapolation of the three finest C3D20 meshes: the converged 3D answer.
  The converged 3D value sits about 0.3% below Timoshenko because the fully clamped
  3D root also blocks Poisson contraction, which beam theory ignores.

### Abaqus cross-check (Learning Edition, 1,000-node cap)

Every case under 1,000 nodes is also written to `runs/cantilever/abaqus/` with
`run_all.bat`.

1. Copy `runs/cantilever/abaqus/` to the Windows VM (Git or Remote Desktop).
2. In an Abaqus Command window in that folder, run `run_all.bat`.
3. Copy the `.dat` files back into `runs/cantilever/abaqus/`.
4. Run `uv run cantilever.py --collect-abaqus`. Abaqus results appear as x markers
   on the plot and as `solver=abaqus` rows in `results.csv`.

## Findings (CalculiX 2.21)

<!-- Update after the Abaqus runs -->
- C3D8 (full integration) locks in bending: 65% of the answer with one element
  through the depth, still 1.2% low at eight.
- C3D8R with one element through the depth hourglasses: 27.7x too flexible. Two
  elements through the depth are the minimum for reduced-integration linear bricks
  in bending.
- C3D8I (incompatible modes) is within 0.7% on the coarsest mesh: the cheapest
  accurate element for bending.
- C3D20 and C3D20R are within 1% everywhere and converge fastest.
