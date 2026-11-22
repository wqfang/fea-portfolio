# fea-portfolio

Hands-on finite element analysis (FEA) projects, run in **CalculiX** (ccx, open
source) on a MacBook and cross-checked in **Abaqus/Standard** (Learning Edition,
on a Windows VM) from the same Abaqus-format input deck.

Three principles shape the work:

- **Methods, not just models.** Every project ends in a reusable script or tool
  that reruns from one command: mesh, deck, solve, parse, plot.
- **Two solvers, one input deck.** Each model runs in both CalculiX and Abaqus,
  and every difference between them is traced to its cause and logged.
- **Product-relevant physics.** Contact, impact on glass, hyperelastic and
  viscoelastic polymers, adhesives and thermal mismatch.

## Projects

| # | Project | Physics | Solvers | Status |
|---|---|---|---|---|
| 1 | [Verification and convergence suite](verification/README.md) | Linear static, contact | ccx + Abaqus | **Done**: 5 benchmarks |
| 2 | User material (UMAT) for a device polymer | Finite-strain hyperelastic + viscoelastic | ccx (UMAT host) + Abaqus (built-in reference) | In progress |
| 3 | Ball drop on cover glass | Impact, contact, explicit dynamics | Abaqus/Explicit + ccx | Planned |
| 4 | Automated design of experiments (DOE) + machine-learning surrogate | Parametric contact or bending | ccx + Python | Planned |
| 5 | Board warpage under thermal load (optional) | Thermo-mechanical, thermal-expansion mismatch | ccx | Planned |

## Project 1: verification and convergence suite

How do you know a model is right before anyone uses it? Five benchmarks with
known answers, each solved over a mesh-refinement sweep with up to five element
types, in both solvers. Every gap, against the reference or between the
solvers, is traced to its cause: element formulation, mesh, boundary
conditions, or the limits of the reference itself.

![NAFEMS LE1 convergence, CalculiX and Abaqus](verification/runs/nafems_le1/le1_convergence.png)

| Benchmark | Reference | Key result |
|---|---|---|
| 0. Single-element tension | Closed form | Exact in both solvers |
| 1. End-loaded cantilever | Timoshenko beam theory | Converged 3D tip deflection 1.5976 mm, 0.34% below beam theory, entirely from the clamped root |
| 2. Plate with a hole | Kirsch (exact); Peterson's handbook fit | Kirsch peak stress to 0.013%; converged Kt_net = 2.522, 0.55% above Peterson's fit |
| 3. Hertz sphere-on-flat contact | Hertz theory | Contact radius within 0.14%; peak pressure +0.7%, from Hertz's small-contact assumption |
| 4. NAFEMS LE1 elliptic membrane | NAFEMS target 92.7 MPa | Converged 92.658 MPa, the target rounded to three figures |

Full models, results tables and derivations:
[`verification/README.md`](verification/README.md).

### Selected findings

- **The same element name can hide different formulations.** Abaqus and
  CalculiX disagree by up to 87% on a coarse `C3D8` cantilever. An independent
  8-node hexahedron solver written in NumPy reproduces both programs to every
  printed digit and shows the cause: Abaqus uses selectively reduced
  integration (the B-bar method), CalculiX plain full integration.
- **A clamped root, not the mesh, limits convergence.** Quadratic elements
  converge at order 1.6 instead of 4 on the cantilever. Removing the stress
  singularity at the clamped edge restores order 3.8 and brings the 3D model
  within 0.003% of Timoshenko beam theory.
- **Shear locking depends on element shape, not layer count.** Adding layers
  through the depth stalls at 72% of the answer; shortening the elements
  removes the locking.
- **Coarse reduced-integration results can be pure artefact.** With one `C3D8R`
  element through the depth, scaling the hourglass stiffness by 10 changes the
  deflection by exactly 10. Enhanced hourglass control brings it within 0.03%.
- **Reference values carry their own error.** Peterson's curve fit is 0.55%
  from a mesh that reproduces the exact Kirsch solution to 0.013%; Hertz theory
  is 0.7% high on peak pressure at a/R = 0.05, an error that scales with a/R;
  the NAFEMS LE1 target is a rounding of 92.658 MPa.
- **CalculiX solves 2D elements in 3D, and it shows.** At the specified 100 mm
  section, CalculiX's LE1 answer is 0.38% high; Abaqus, solving true plane
  stress, matches CalculiX's thin-section result to 0.01%. The same expansion
  makes thin `CPS4R` elements unusable and zeroes the contact pressure on the
  axis of axisymmetric models.
- **Straight-edged elements put a corner on a curved boundary.** Linear `CPS4`
  overshoots the LE1 stress by up to 2% and converges from above. Quadratic
  elements with their edges forced straight overshoot too (+2.8% at 96
  elements, +1.0% at 384), which isolates the polygonal geometry as the main
  cause of the Abaqus overshoot.
- **Solver pitfall:** multithreaded CalculiX 2.23 returned wrong answers with no
  error message, so every run is single-threaded.

### Solver differences

What changes between the CalculiX and Abaqus decks, and what behaves
differently with the same deck.

| Topic | CalculiX 2.23 | Abaqus/Standard Learning Edition 2026 | Effect seen |
|---|---|---|---|
| `C3D8` | Plain full integration | Selectively reduced integration (B-bar) | Up to 87% apart on a coarse bending mesh |
| `C3D8R` hourglass control | Own default; ignores `*SECTION CONTROLS` (warns) | Stiffness method by default; `*SECTION CONTROLS` honoured | 2.8x vs 92x too flexible with one element through the depth |
| Plane-stress `CPS*` elements | Expanded to 3D bricks of the section thickness | True 2D | Thickness changes the stress (+0.38% on LE1); thin `CPS4R` far too stiff |
| Axisymmetric `CAX*` elements | Expanded to a 2-degree wedge; forces reported for the wedge (x180) | True 2D; forces for the full 360 degrees | Zero contact pressure on the axis in CalculiX |
| Incompatible-mode `CPS4I`, `CAX4I` | Not available | Available | Best first-order element for notch stress (Abaqus only) |
| Contact | Node-to-surface did not converge; surface-to-surface used | Surface-to-surface | Same force; quadratic contact pressure 70% off in CalculiX at 4 elements across the contact, 2% in Abaqus |
| Results output | `*NODE FILE`, `*CONTACT FILE` to `.frd` | `*EL PRINT, POSITION=AVERAGED AT NODES`, `*CONTACT PRINT` to `.dat` (4 significant figures) | Separate parsers; Abaqus resolution 0.03% |
| Reaction totals | `*NODE PRINT, TOTALS=ONLY` | `TOTALS=YES` (ONLY is rejected) | Keyword differs |
| User subroutines (UMAT) | Compiled into the solver | Refused by the Learning Edition | UMATs run in CalculiX; Abaqus built-in models are the reference |
| Size limit | None (laptop memory) | 1,000 nodes | Fine meshes in CalculiX only |

## Quick start

Requirements: [uv](https://docs.astral.sh/uv/) and CalculiX (`ccx`). Abaqus is
needed only for the cross-checks; its results are already committed.

```bash
brew install costerwi/homebrew-calculix/calculix-ccx   # macOS
cd verification
uv run cantilever.py        # about 1 minute
uv run plate_hole.py        # about 10 seconds
uv run hertz_contact.py     # about 3 minutes
uv run nafems_le1.py        # about 2 minutes
```

Each script writes its input decks, a `results.csv` and a convergence plot to
`verification/runs/<benchmark>/`, plus the Abaqus decks that fit under the
1,000-node cap. Results are also converted to VTK files for viewing in
ParaView.

## Repository layout

```
verification/                   Project 1
  README.md                      full write-up of every benchmark
  cantilever.py                  benchmark 1
  plate_hole.py                  benchmark 2
  hertz_contact.py               benchmark 3
  nafems_le1.py                  benchmark 4
  *_geometry.py                  problem schematics and mesh figures
  *_check.py, *_sensitivity.py   investigations
  fea_common.py                  deck writing, solver runs, result parsing
  frd2vtu.py                     CalculiX results to VTK
  tension*.inp                   benchmark 0
  runs/                          committed results, plots and Abaqus decks
```

## Tools

CalculiX 2.23 on macOS; Abaqus/Standard Learning Edition 2026 on a Windows VM;
Python (NumPy, Matplotlib) for meshing, deck generation, parsing and plots;
ParaView for viewing results.

The Abaqus Learning Edition is for learning; these are personal,
non-commercial projects.

## License

[GPL-3.0](LICENSE) (GNU General Public License, version 3).
