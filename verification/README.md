# FEA verification suite

Benchmarks with known answers, each run in **CalculiX** (ccx, open source) and
**Abaqus/Standard** from the same Abaqus-format input deck, with mesh-convergence
studies. Units: mm, N, MPa.

| # | Benchmark | Checks | Status |
|---|-----------|--------|--------|
| 0 | Single-element uniaxial tension (`tension.inp`) | Material and boundary-condition setup; ccx vs Abaqus | Done: both match closed form |
| 1 | End-loaded cantilever (`cantilever.py`) | Bending, shear locking, hourglassing, convergence rate | Done: CalculiX and Abaqus |
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
  theory ignores.

### Abaqus cross-check (Learning Edition, 1,000-node cap)

Every case under 1,000 nodes is also written to `runs/cantilever/abaqus/` with
`run_all.bat`.

1. Copy `runs/cantilever/abaqus/` to the Windows VM (Git or Remote Desktop).
2. In an Abaqus Command window in that folder, run `run_all.bat`.
3. Copy the `.dat` files back into `runs/cantilever/abaqus/`.
4. Run `uv run cantilever.py --collect-abaqus`. Abaqus results appear as x markers
   on the plot and as `solver=abaqus` rows in `results.csv`. Without the flag the
   script ignores the Abaqus `.dat` files.

## Findings (CalculiX 2.23, Abaqus Learning Edition 2026)

Tip deflection as a fraction of Timoshenko (1.6031 mm). Abaqus runs only the
cases under its 1,000-node cap.

| Element | Solver | n=1 | n=2 | n=4 | n=8 |
|---|---|---|---|---|---|
| C3D8 | ccx | 0.646 | 0.875 | 0.962 | 0.988 |
| C3D8 | Abaqus | 1.207 | 0.980 | | |
| C3D8R | ccx | 2.77 | 1.175 | 1.035 | 1.006 |
| C3D8R | Abaqus | 91.7 | 1.325 | | |
| C3D8I | ccx | 0.9936 | 0.9937 | 0.9953 | 0.9961 |
| C3D8I | Abaqus | 0.9948 | 0.9939 | | |
| C3D20 | ccx | 0.9901 | 0.9948 | 0.9960 | 0.9964 |
| C3D20 | Abaqus | 0.9901 | | | |
| C3D20R | ccx | 0.9912 | 0.9953 | 0.9961 | 0.9965 |
| C3D20R | Abaqus | 0.9912 | | | |

**Element behaviour**
- **C3D8 in CalculiX** (full 2x2x2 integration) shear-locks in bending: 65% of
  the answer with one element through the depth, still 1.2% low at eight.
- **C3D8R** with one element through the depth hourglasses: bending is then
  exactly an hourglass mode, so the only bending stiffness is the solver's
  artificial hourglass control. CalculiX gives 2.8x, Abaqus 92x. Two elements
  through the depth are the minimum for reduced-integration linear bricks in
  bending, and even then they are 17-32% too flexible.
- **C3D8I** (incompatible modes) is within 0.7% on the coarsest mesh in both
  solvers: the cheapest accurate element for bending.
- **C3D20 and C3D20R** are within 1% everywhere and converge fastest.

**CalculiX vs Abaqus**
- **C3D20 and C3D20R agree to 1e-7.** The two solvers read the same deck with the
  same loads and supports, which validates the deck writer, the load
  distribution and the result parsing.
- **C3D8I agrees to 1e-3** (n=1) and 1e-4 (n=2). The small gap most likely comes from
  differences in how the two programs implement incompatible modes.
- **C3D8R disagrees by design**: hourglass control is solver-specific, so its
  coarse-mesh results are not comparable across programs.
- **C3D8 is not the same element in the two programs.** Abaqus's C3D8 at n=1 is
  21% *too flexible*, where a full-integration element should lock and be too
  stiff. The run was clean (no warnings, linear, residual 6e-10), and the
  deformed shape matched CalculiX's up to a scale factor.

### Why Abaqus and CalculiX disagree on C3D8 (`c3d8_formulation_check.py`)

Abaqus/Standard's C3D8 uses *selectively reduced integration*: the volumetric
strain at each integration point is replaced by the element average (the
mean-dilatation, or B-bar, method), which avoids volumetric locking in nearly
incompressible materials. CalculiX's C3D8 is plain full integration.

`c3d8_formulation_check.py` implements an independent 8-node hex solver in numpy
in both variants and runs it on the same mesh and loads:

```bash
uv run c3d8_formulation_check.py
```

| Case | numpy full | ccx C3D8 | numpy B-bar | Abaqus C3D8 |
|---|---|---|---|---|
| n=1, nu=0.3 | 1.036184 | 1.036184 | 1.934318 | 1.934318 |
| n=2, nu=0.3 | 1.403038 | 1.403038 | 1.571293 | 1.571293 |
| n=1, nu=0 | 1.068000 | 1.068000 | 1.487937 | 1.487937 |

Each solver is reproduced to every printed digit, so the difference is entirely
the element formulation.

**Why B-bar makes a coarse bending mesh too flexible.** In bending, the axial
strain changes sign through the depth, so its element average is close to zero.
Replacing the volumetric part of the strain with that average removes part of
the real bending strain, and the element loses bending stiffness. With one
element through the depth this loss outweighs shear locking, which pushes the
other way. Both effects shrink with refinement (Abaqus 0.98, CalculiX 0.875 at
n=2).

**The nu = 0 case** was a check on that explanation. Setting nu = 0 removes the
Lamé lambda term, but B-bar still modifies the strain itself, so Abaqus stays
too flexible and even shows sideways displacement (U3 = 1.4e-5 mm) where plain
full integration gives zero. Its Abaqus deck and `run_diag.bat` are in
`runs/cantilever/diag/abaqus/`.

**Takeaway.** The same element name can hide different formulations. Compare
solvers on elements with a common definition (C3D20, C3D20R, approximately
C3D8I), and check the theory manual before treating a "C3D8" result from one
code as a reference for another.
