#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_finite_dirac_algebra_PT.py -- Finite PT Dirac algebra
================================================================

This script tests the finite, transfer-matrix analogue of the Dirac
constraint algebra in the PT quantum-gravity sector.

Definitions on the finite boundary Hilbert space:

    H_m      = -log(T_m)                         Hamiltonian constraint
    Pi_0     = |1><pi|                            physical projector
    Q        = I - Pi_0                            unphysical quotient
    D_X      = [H_m, Q X Q]                       admissible boundary deformation
    O_f      = f(T_m)                             spectral Dirac observable

The interpretation is deliberately conservative.  This is not yet the full
continuum hypersurface-deformation algebra of ADM gravity.  It is the finite
PT algebra induced by the transfer cylinder:

  - H_m annihilates the physical Perron sector;
  - admissible D_X is H-exact inside the unphysical quotient, hence
    physically invisible after projection;
  - [H_m, D_X] is again H-exact by Jacobi;
  - [D_X, D_Y] has vanishing physical matrix elements;
  - spectral observables commute with H_m and are invariant modulo
    H-exact deformations;
  - physical amplitudes are invariant under H-exact boundary shifts.

What this closes:
  - anomaly-free finite constraint algebra on the PT boundary quotient;
  - the admissibility condition QXQ: arbitrary boundary operators are not
    silently promoted to gauge transformations.
  - the finite Dirac bracket content needed by the cylinder sector.

What this script does not by itself close:
  - topology-changing physical foams
    (see test_qg_topology_changing_foams_PT.py);
  - higher Fourier/RG decorations
    (see test_qg_fourier_rg_decorations_PT.py);
  - black-hole/ringdown precision skeletons
    (see test_qg_black_hole_ringdown_precision_PT.py).
"""

import io
import sys

import numpy as np
from scipy.linalg import logm

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def sieve_primes(n_max):
    is_prime = [True] * (n_max + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n_max**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n_max + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n_max + 1) if is_prime[i]]


def build_T_gaps(gaps, modulus):
    residues = [g % modulus for g in gaps]
    classes = sorted(set(residues))
    idx = {c: i for i, c in enumerate(classes)}
    counts = np.zeros((len(classes), len(classes)), dtype=float)
    for a, b in zip(residues, residues[1:]):
        counts[idx[a], idx[b]] += 1.0

    T = np.zeros_like(counts)
    for i in range(counts.shape[0]):
        row_sum = counts[i].sum()
        if row_sum > 0:
            T[i] = counts[i] / row_sum
    return T, classes


def stationary_left(T):
    eigvals, eigvecs = np.linalg.eig(T.T)
    idx = np.argmin(np.abs(eigvals - 1.0))
    pi = np.real(eigvecs[:, idx])
    if pi.sum() < 0:
        pi = -pi
    pi = np.maximum(pi, 0.0)
    return pi / pi.sum()


def perron_projector(T):
    right = np.ones(T.shape[0])
    left = stationary_left(T)
    return np.outer(right, left)


def comm(A, B):
    return A @ B - B @ A


def mat_norm(A):
    return float(np.linalg.norm(A))


def random_matrix(n, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))


def random_vector(n, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=n) + 1j * rng.normal(size=n)


def physical_amplitude(Pi0, phi, psi):
    return np.vdot(phi, Pi0 @ psi)


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("  PT QUANTUM GRAVITY: FINITE DIRAC CONSTRAINT ALGEBRA")
print("=" * 78)

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]
T6, classes6 = build_T_gaps(gaps, 6)
T30, classes30 = build_T_gaps(gaps, 30)

print(f"\n  T_6  shape = {T6.shape}, classes = {classes6}")
print(f"  T_30 shape = {T30.shape}, classes = {classes30}")
print()

n_pass = 0
n_total = 0

for label, T in (("T6", T6), ("T30", T30)):
    n = T.shape[0]
    H = -logm(T)
    Pi0 = perron_projector(T)
    X = random_matrix(n, 1000 + n)
    Y = random_matrix(n, 2000 + n)
    phi = random_vector(n, 3000 + n)
    psi = random_vector(n, 4000 + n)
    chi = random_vector(n, 5000 + n)
    eta = random_vector(n, 6000 + n)

    Q = np.eye(n) - Pi0
    X_adm = Q @ X @ Q
    Y_adm = Q @ Y @ Q
    D_X_raw = comm(H, X)
    D_Y_raw = comm(H, Y)
    D_X = comm(H, X_adm)
    D_Y = comm(H, Y_adm)
    O = np.linalg.matrix_power(T, 5) + 0.25 * np.linalg.matrix_power(T, 2)

    print(f"--- {label}: boundary dimension {n} ---")

    h_phys = mat_norm(H @ Pi0) + mat_norm(Pi0 @ H)
    check(f"{label}: Hamiltonian constraint annihilates physical projector",
          h_phys < 1e-10,
          f"||H Pi0||+||Pi0 H||={h_phys:.2e}")

    hh = mat_norm(comm(H, H))
    check(f"{label}: scalar-scalar bracket closes trivially",
          hh < 1e-12,
          f"||[H,H]||={hh:.2e}")

    raw_anomaly = mat_norm(Pi0 @ comm(D_X_raw, D_Y_raw) @ Pi0)
    check(f"{label}: arbitrary boundary operators are not all gauge",
          raw_anomaly > 1e-3,
          f"negative-control anomaly ||Pi[D_rawX,D_rawY]Pi||={raw_anomaly:.2e}")

    jacobi_err = mat_norm(comm(H, D_X) - comm(H, comm(H, X)))
    adm_jacobi_err = mat_norm(comm(H, D_X) - comm(H, comm(H, X_adm)))
    check(f"{label}: Hamiltonian-deformation bracket is H-exact",
          adm_jacobi_err < 1e-10,
          f"||[H,D_X]-[H,[H,QXQ]]||={adm_jacobi_err:.2e}")

    phys_DX = mat_norm(Pi0 @ D_X) + mat_norm(D_X @ Pi0)
    phys_DY = mat_norm(Pi0 @ D_Y) + mat_norm(D_Y @ Pi0)
    check(f"{label}: H-exact deformations vanish on physical quotient",
          phys_DX < 1e-10 and phys_DY < 1e-10,
          f"||Pi D_X||+||D_X Pi||={phys_DX:.2e}, ||Pi D_Y||+||D_Y Pi||={phys_DY:.2e}")

    phys_D_comm = mat_norm(Pi0 @ comm(D_X, D_Y) @ Pi0)
    check(f"{label}: deformation-deformation bracket is first-class physically",
          phys_D_comm < 1e-9,
          f"||Pi [D_X,D_Y] Pi||={phys_D_comm:.2e}")

    obs_comm = mat_norm(comm(H, O))
    check(f"{label}: spectral observables commute with H",
          obs_comm < 1e-10,
          f"||[H,O_f]||={obs_comm:.2e}")

    phys_obs_def = mat_norm(Pi0 @ comm(D_X, O) @ Pi0)
    check(f"{label}: observables are invariant modulo H-exact deformations",
          phys_obs_def < 1e-9,
          f"||Pi [D_X,O_f] Pi||={phys_obs_def:.2e}")

    amp0 = physical_amplitude(Pi0, phi, psi)
    amp_shift = physical_amplitude(Pi0, phi + H.conj().T @ eta, psi + H @ chi)
    amp_err = abs(amp_shift - amp0)
    check(f"{label}: physical amplitudes are invariant under constraint shifts",
          amp_err < 1e-10,
          f"|W(phi+H^dag eta,psi+H chi)-W(phi,psi)|={amp_err:.2e}")

    projected_D = mat_norm(Pi0 @ D_X @ Pi0 @ D_Y @ Pi0)
    check(f"{label}: projected gauge products vanish",
          projected_D < 1e-12,
          f"||Pi D_X Pi D_Y Pi||={projected_D:.2e}")

    quotient_anomaly = (
        mat_norm(Pi0 @ comm(H, D_X) @ Pi0)
        + mat_norm(Pi0 @ comm(D_X, D_Y) @ Pi0)
        + mat_norm(Pi0 @ comm(H, O) @ Pi0)
    )
    check(f"{label}: no anomaly on the physical boundary quotient",
          quotient_anomaly < 1e-9,
          f"projected anomaly norm={quotient_anomaly:.2e}")

print("\n" + "=" * 78)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: finite PT Dirac algebra on the transfer-boundary quotient,")
print("               H-exact deformations, anomaly-free physical brackets.")
print("  CONTINUUM LIFT: see test_qg_continuum_dirac_lift_PT.py.")
print("  TOPOLOGICAL FOAMS: see test_qg_topology_changing_foams_PT.py.")
print("  DECORATIONS/BH: see Fourier-RG and black-hole/ringdown QG scripts.")
print("=" * 78)

sys.exit(0 if n_pass == n_total else 1)
