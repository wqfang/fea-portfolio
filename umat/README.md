# User material subroutines (UMATs) verified in CalculiX and Abaqus

Project 2 of the portfolio: user material subroutines (UMATs) written to the
**Abaqus/Standard interface**, compiled into **CalculiX** (ccx, open source), and
verified at three levels: a material-point driver, single-element tests against
closed-form solutions, and cross-checks against the built-in models in CalculiX
and Abaqus.

| Step | Material | Status |
|---|---|---|
| 2.1 | Compressible Neo-Hookean hyperelasticity | Done (this page) |
| 2.2 | Finite-strain viscoelasticity (Prony series) on top of 2.1 | Next |
| 2.3 | Part-level model (e.g. button dome or O-ring) and parameter fit to test data | Planned |

Units: mm, N, MPa.

## Why CalculiX hosts the UMAT

The Abaqus Learning Edition refuses user subroutines ("User subroutines are not
supported in the Learning Edition"). CalculiX can call a UMAT with the Abaqus
argument list when it is compiled into the solver, so the same Fortran source
runs in CalculiX here and would run unchanged in a licensed Abaqus. The Learning
Edition still provides the independent reference: the same decks with the
built-in `*HYPERELASTIC, NEO HOOKE` model.

## Step 2.1: compressible Neo-Hookean

`src/umat_neohooke.f` implements the Abaqus form of the strain energy

    U = C10 (I1bar - 3) + (J - 1)^2 / D1,    Bbar = J^(-2/3) F F^T,  I1bar = tr(Bbar)

with `PROPS(1) = C10`, `PROPS(2) = D1` (small-strain shear modulus 2 C10, bulk
modulus 2/D1), so it is directly comparable with the built-in model using the
same constants. It returns

- **STRESS**: the Cauchy stress, sigma = (2 C10 / J)(Bbar - I1bar/3 I) + (2/D1)(J - 1) I
- **DDSDDE**: (1/J) times the Jaumann tangent of the Kirchhoff stress, the
  tangent Abaqus/Standard expects from a finite-strain UMAT:

      J C_ijkl = 2 C10 [ 1/2 (d_ik B_jl + B_ik d_jl + d_il B_jk + B_il d_jk)
                         - 2/3 (B_ij d_kl + d_ij B_kl) + 2/9 I1bar d_ij d_kl ]
                 + (2/D1) J (2J - 1) d_ij d_kl          (B = Bbar, d = Kronecker delta)

Test constants: C10 = 0.5 MPa, D1 = 0.02 /MPa (shear modulus 1 MPa, bulk
modulus 100 MPa, small-strain Poisson's ratio 0.495), typical of a filled rubber.

### The frame convention (the main finding)

CalculiX passes the full deformation gradient F to the UMAT, but its `ABAQUSNL`
wrapper (`umat_abaqusnl.f`) expects the **corotated** Cauchy stress
R^T sigma R back, not the global-frame stress Abaqus expects for elements without
a local orientation. It then converts that to its own second Piola-Kirchhoff
stress with the right stretch tensor U. For an isotropic material the corotated
stress is the same formula with Cbar = J^(-2/3) F^T F in place of Bbar, so the
UMAT has one preprocessor switch:

- default (Abaqus): Bbar, global frame
- `-DCCX_COROT` (CalculiX, set by the build script): Cbar, corotated frame

Built without the switch, the UMAT is still exact in uniaxial tension (no
rotation) but wrong as soon as the material rotates:

| Case (CalculiX) | Global-frame UMAT | Corotated UMAT | Exact |
|---|---|---|---|
| Uniaxial, lambda = 2: sigma_11 (MPa) | 3.429 | 3.429 | 3.429 |
| Simple shear, gamma = 1: sigma_12 (MPa) | **0.200** | 1.000 | 1.000 |
| Simple shear, gamma = 1: sigma_11 (MPa) | **1.267** | 0.667 | 0.667 |
| Stretch 1.5, rotated 90 deg: sigma_11 (MPa) | **50.424** | 49.788 | 49.788 |

In the rotated case the wrong build returns sigma_22 in the sigma_11 slot: the
stress does not turn with the body. **A UMAT that passes uniaxial tension can
still be wrong**, which is why the test set includes simple shear and a rigid
rotation.

## Verification

### 1. Material-point driver (`tests/test_material_point.py`)

`tests/driver.f90` calls the UMAT the way Abaqus does for one integration point.
The test compiles it in both frames, feeds both 25 random deformation gradients
(large stretch and rotation, J from 0.3 to 3), and compares with an independent
NumPy implementation (`ref/neohooke.py`):

| Check | Max relative difference |
|---|---|
| Abaqus build: stress vs closed form | 3.0e-15 |
| Abaqus build: DDSDDE vs analytic tangent | 9.2e-16 |
| Abaqus build: DDSDDE vs finite-difference Jaumann tangent | 7.2e-07 |
| CalculiX build: stress vs R^T sigma R | 3.2e-15 |
| CalculiX build: DDSDDE vs corotated tangent | 9.2e-16 |

The finite-difference tangent perturbs F with a pure symmetric rate,
F_hat = F + (eps/2)(e_k e_l^T + e_l e_k^T) F, so the Jaumann rate equals the
ordinary change of the Kirchhoff stress (Sun, Chaikof and Levenston, 2008). Its
7e-7 difference is the finite-difference truncation error.

### 2. Single-element tests in CalculiX (`tests/single_element.py`)

One C3D8 unit cube, geometrically nonlinear, 20 increments. Each case runs with
the UMAT and with CalculiX's built-in `*HYPERELASTIC, NEO HOOKE`; both are
compared with the closed-form Cauchy stress at the final F, computed from the
nodal displacements. For uniaxial and biaxial tension the free lateral stretch
is also solved independently from the zero-stress condition.

| Case | Final state | UMAT vs closed form | Free stretch: FE / exact |
|---|---|---|---|
| Uniaxial tension | lambda = 2, sigma_11 = 3.4287 MPa | 6.6e-07 | 0.711136 / 0.711136 |
| Equibiaxial tension | lambda = 1.5, sigma_11 = 2.0024 MPa | 1.6e-06 | 0.450377 / 0.450378 |
| Simple shear | gamma = 1, sigma_12 = 1.0000 MPa | 3.3e-08 | |
| Stretch then 90 deg rotation | sigma_22 = 50.424 MPa | 8.2e-08 | |

The differences are at the 7 significant digits CalculiX prints. The UMAT and the
built-in model agree to every printed digit, with the same number of Newton
iterations (40 to 41 over 20 increments).

**Control.** The same UMAT deck run in stock CalculiX, whose `umat.f` is an empty
template, gives sigma_12 = 0.012 MPa after the first increment and stops, so the
rebuilt solver really is calling this UMAT.

### 3. Abaqus cross-check (built-in model)

`tests/single_element.py` also writes the four cases with the built-in
Neo-Hookean model to `runs/single_element/abaqus/`, with `run_all.bat`.

1. Copy `runs/single_element/abaqus/` to the Windows VM.
2. Run `run_all.bat` in an Abaqus Command window.
3. Copy the `.dat` files back and run
   `uv run tests/single_element.py --ccx bin/ccx_2.23_umat --collect-abaqus`.

Abaqus's C3D8 uses selectively reduced (B-bar) integration, but every case here
is a homogeneous deformation, so it gives the same answer as full integration.

## Build and run

```bash
brew install gcc arpack pkg-config          # gfortran, ARPACK eigen-solver library
./build_ccx_umat.sh                          # ~2 min: downloads ccx 2.23 + SPOOLES, builds bin/ccx_2.23_umat
uv run tests/test_material_point.py          # seconds
uv run tests/single_element.py --ccx bin/ccx_2.23_umat
```

`build_ccx_umat.sh` mirrors the Homebrew formula (same sources, patches and
compiler flags), replaces ccx's empty `umat.f` with the chosen UMAT
(`--umat src/....f`), and adds `-DCCX_COROT`. In a deck, select the UMAT with a
material name starting with `ABAQUSNL` (finite strain) and
`*USER MATERIAL, CONSTANTS=2`. The tests run ccx on one thread (see
`verification/README.md` on multithreaded ccx 2.23).

Tested here with CalculiX 2.21 on Linux (gfortran 13); the macOS build of 2.23
follows the same steps.

## Notes on the CalculiX interface

From the ccx 2.23 manual (section 8.5.2) and `umat_abaqusnl.f`:

- `ABAQUS` + name calls the routine with small-strain kinematics; `ABAQUSNL` +
  name with finite strain (corotated logarithmic strain in, corotated Cauchy
  stress out). The prefix is removed before the name reaches the routine.
- Not supplied: `SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, PREDEF, DPRED, DROT,
  PNEWDT, CELENT, LAYER, KSPT`. This UMAT computes stress from `DFGRD1` alone,
  so it needs none of them.
- CalculiX symmetrizes DDSDDE (21 constants instead of 36) and maps it to its
  own tangent with U^-1, which is exact for the stress but only approximate for
  the tangent; convergence was still two Newton iterations per increment here.
