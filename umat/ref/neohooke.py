"""Reference implementation of Abaqus-form compressible Neo-Hookean hyperelasticity.

Strain energy (Abaqus *HYPERELASTIC, NEO HOOKE):
    U = C10 (I1bar - 3) + (J - 1)^2 / D1
with Bbar = J^(-2/3) F F^T and I1bar = tr(Bbar).
Small-strain equivalents: shear modulus mu0 = 2 C10, bulk modulus K0 = 2 / D1.

Used to check the Fortran UMAT and the CalculiX runs. Voigt order follows Abaqus:
11, 22, 33, 12, 13, 23 (shear strains are engineering strains in DDSDDE).
"""
from __future__ import annotations

import numpy as np

I3 = np.eye(3)
VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def kirchhoff(F: np.ndarray, C10: float, D1: float) -> np.ndarray:
    """Kirchhoff stress tau = J sigma."""
    J = np.linalg.det(F)
    Bb = J ** (-2.0 / 3.0) * F @ F.T
    I1b = np.trace(Bb)
    return 2.0 * C10 * (Bb - I1b / 3.0 * I3) + 2.0 / D1 * J * (J - 1.0) * I3


def cauchy(F: np.ndarray, C10: float, D1: float) -> np.ndarray:
    return kirchhoff(F, C10, D1) / np.linalg.det(F)


def ddsdde_analytic(F: np.ndarray, C10: float, D1: float, Bb: np.ndarray | None = None) -> np.ndarray:
    """Abaqus DDSDDE for finite-strain hyperelasticity: (1/J) x Jaumann tangent of Kirchhoff stress.

    Bb can be replaced by the corotated Cbar = R^T Bbar R (CalculiX frame).
    """
    J = np.linalg.det(F)
    if Bb is None:
        Bb = J ** (-2.0 / 3.0) * F @ F.T
    I1b = np.trace(Bb)
    d = I3
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    C[i, j, k, l] = 2.0 * C10 * (
                        0.5 * (d[i, k] * Bb[j, l] + Bb[i, k] * d[j, l] + d[i, l] * Bb[j, k] + Bb[i, l] * d[j, k])
                        - 2.0 / 3.0 * (Bb[i, j] * d[k, l] + d[i, j] * Bb[k, l])
                        + 2.0 / 9.0 * I1b * d[i, j] * d[k, l]
                    ) + 2.0 / D1 * J * (2.0 * J - 1.0) * d[i, j] * d[k, l]
    C /= J
    return np.array([[C[a, b, c, e] for (c, e) in VOIGT] for (a, b) in VOIGT])


def ddsdde_fd(F: np.ndarray, C10: float, D1: float, eps: float = 1e-7) -> np.ndarray:
    """Finite-difference Jaumann tangent (Sun, Chaikof & Levenston 2008 perturbation).

    Perturb F_hat = F + (eps/2)(e_k e_l^T + e_l e_k^T) F, a pure symmetric rate with no spin,
    so the Jaumann rate of tau equals its ordinary change. Column for (k,l):
        (tau(F_hat) - tau(F)) / (J eps)
    """
    J = np.linalg.det(F)
    tau0 = kirchhoff(F, C10, D1)
    D = np.zeros((6, 6))
    for col, (k, l) in enumerate(VOIGT):
        E = np.zeros((3, 3))
        E[k, l] += 0.5
        E[l, k] += 0.5
        Fh = F + eps * E @ F
        dtau = (kirchhoff(Fh, C10, D1) - tau0) / (J * eps)
        D[:, col] = [dtau[a, b] for (a, b) in VOIGT]
    return D


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    C10, D1 = 0.5, 0.02  # MPa, 1/MPa: mu0 = 1 MPa, K0 = 100 MPa (nu0 ~ 0.495)
    worst = 0.0
    for _ in range(20):
        F = I3 + 0.4 * rng.standard_normal((3, 3))
        if np.linalg.det(F) < 0.3:
            continue
        A, N = ddsdde_analytic(F, C10, D1), ddsdde_fd(F, C10, D1)
        worst = max(worst, np.abs(A - N).max() / np.abs(A).max())
    print(f"max relative difference analytic vs finite-difference tangent: {worst:.2e}")
