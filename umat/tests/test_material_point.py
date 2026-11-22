"""Material-point check of the UMAT: stress and tangent against the Python reference.

Compiles tests/driver.f90 with the UMAT twice (Abaqus frame and CalculiX corotated
frame), feeds both the same random deformation gradients, and compares:
  - Abaqus build:   STRESS vs Cauchy stress sigma; DDSDDE vs analytic and
                    finite-difference Jaumann tangents
  - CalculiX build: STRESS vs R^T sigma R; DDSDDE vs the tangent rotated the same way

Usage (from the umat/ folder):  uv run tests/test_material_point.py
Needs gfortran on PATH (brew install gcc).
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ref"))
from neohooke import VOIGT, cauchy, ddsdde_analytic, ddsdde_fd  # noqa: E402

C10, D1 = 0.5, 0.02  # MPa, 1/MPa  (mu0 = 1 MPa, K0 = 100 MPa)
N_CASES = 25


def gfortran() -> str:
    for name in ("gfortran", *(f"gfortran-{v}" for v in range(16, 10, -1))):
        if shutil.which(name):
            return name
    sys.exit("gfortran not found: brew install gcc")


def build(tmp: Path, corot: bool) -> Path:
    (tmp / "ABA_PARAM.INC").write_text("      IMPLICIT REAL*8(A-H,O-Z)\n")
    exe = tmp / ("driver_ccx" if corot else "driver_abaqus")
    cmd = [gfortran(), "-cpp", "-O2", f"-I{tmp}", "-o", str(exe),
           str(HERE / "driver.f90"), str(HERE.parent / "src" / "umat_neohooke.f")]
    if corot:
        cmd.insert(2, "-DCCX_COROT")
    subprocess.run(cmd, check=True)
    return exe


def run(exe: Path, Fs: list[np.ndarray]) -> list[tuple[np.ndarray, np.ndarray]]:
    inp = f"{C10} {D1}\n" + "\n".join(" ".join(f"{v:.17g}" for v in F.ravel()) for F in Fs)
    out = np.loadtxt(subprocess.run([str(exe)], input=inp, capture_output=True,
                                    text=True, check=True).stdout.splitlines())
    return [(out[7 * c], out[7 * c + 1:7 * c + 7]) for c in range(len(Fs))]


def voigt(T: np.ndarray) -> np.ndarray:
    return np.array([T[a, b] for (a, b) in VOIGT])


def polar_R(F: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(F)
    return U @ Vt


def rel(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.abs(a - b).max() / max(np.abs(b).max(), 1e-300))


def main() -> None:
    rng = np.random.default_rng(42)
    Fs = []
    while len(Fs) < N_CASES:  # random large deformations with rotation, J in a sane range
        F = np.eye(3) + 0.5 * rng.standard_normal((3, 3))
        if 0.3 < np.linalg.det(F) < 3.0:
            Fs.append(F)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        abq = run(build(tmp, corot=False), Fs)
        ccx = run(build(tmp, corot=True), Fs)

    worst = dict(stress_abq=0.0, tan_abq_analytic=0.0, tan_abq_fd=0.0, stress_ccx=0.0, tan_ccx=0.0)
    for F, (s_a, D_a), (s_c, D_c) in zip(Fs, abq, ccx):
        sig = cauchy(F, C10, D1)
        worst["stress_abq"] = max(worst["stress_abq"], rel(s_a, voigt(sig)))
        worst["tan_abq_analytic"] = max(worst["tan_abq_analytic"], rel(D_a, ddsdde_analytic(F, C10, D1)))
        worst["tan_abq_fd"] = max(worst["tan_abq_fd"], rel(D_a, ddsdde_fd(F, C10, D1)))
        R = polar_R(F)
        J = np.linalg.det(F)
        worst["stress_ccx"] = max(worst["stress_ccx"], rel(s_c, voigt(R.T @ sig @ R)))
        Cb = J ** (-2.0 / 3.0) * F.T @ F  # = R^T Bbar R
        worst["tan_ccx"] = max(worst["tan_ccx"], rel(D_c, ddsdde_analytic(F, C10, D1, Bb=Cb)))

    labels = {
        "stress_abq": "Abaqus build, stress vs closed form",
        "tan_abq_analytic": "Abaqus build, DDSDDE vs analytic tangent",
        "tan_abq_fd": "Abaqus build, DDSDDE vs finite-difference tangent",
        "stress_ccx": "CalculiX build, stress vs R^T sigma R",
        "tan_ccx": "CalculiX build, DDSDDE vs corotated tangent",
    }
    tol = {"tan_abq_fd": 1e-5}
    ok = True
    print(f"{N_CASES} random deformation gradients, C10 = {C10} MPa, D1 = {D1} /MPa")
    for k, v in worst.items():
        passed = v < tol.get(k, 1e-10)
        ok &= passed
        print(f"  {'PASS' if passed else 'FAIL'}  {labels[k]:52s} max rel. diff {v:.1e}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
