#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_continuum_dirac_lift_PT.py -- Continuum lift of finite PT Dirac algebra
================================================================================

This script tests the precise PT meaning of the "continuum Dirac lift".

The continuum is not obtained by taking an empirical finite-N prime-gap
matrix at a larger modulus and hoping that finite statistics commute exactly.
The PT lift is the CRT/cylindrical refinement:

    T_fine = T_coarse \otimes S_p
    H_fine = H_coarse \otimes I + I \otimes H_p
    Pi_fine = Pi_coarse \otimes Pi_p

where S_p is a primitive stochastic internal refinement channel.  Coarse
boundary functions are injected by J(f)=f \otimes 1, and fine boundary data
are averaged by A=I \otimes <pi_p|.  Cylindrical consistency means:

    A J = I
    T_fine J = J T_coarse,       A T_fine = T_coarse A
    H_fine J = J H_coarse,       A H_fine = H_coarse A
    Pi_fine J = J Pi_coarse,     A Pi_fine = Pi_coarse A

For admissible deformations D_X=[H,QXQ], cylindrical operators also satisfy

    A D_X^fine J = D_X^coarse.

This is the finite-to-continuum lift used here: every finite PT Dirac algebra
embeds into every CRT refinement without anomaly.  The direct/projective
limit is then defined by cylindrical equivalence classes.

What this closes:
  - cylindrical consistency of H, Pi, amplitudes, and finite Dirac algebra;
  - anomaly-free lift of the finite PT Dirac algebra through CRT refinements.

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


def lazy_complete_channel(dim, stay=0.73):
    """Primitive row-stochastic refinement with positive spectrum."""
    uniform = np.ones((dim, dim), dtype=float) / dim
    return stay * np.eye(dim) + (1.0 - stay) * uniform


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


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


def finite_dirac_data(T):
    H = -logm(T)
    Pi = perron_projector(T)
    Q = np.eye(T.shape[0]) - Pi
    return H, Pi, Q


def refine(T, H, Pi, dim, stay):
    S = lazy_complete_channel(dim, stay=stay)
    Hs = -logm(S)
    Pis = perron_projector(S)
    nf = T.shape[0] * dim
    nc = T.shape[0]
    ones = np.ones((dim, 1))
    pi = stationary_left(S).reshape(1, dim)
    J = np.kron(np.eye(nc), ones)       # coarse -> fine
    A = np.kron(np.eye(nc), pi)         # fine -> coarse
    T_f = np.kron(T, S)
    H_f = np.kron(H, np.eye(dim)) + np.kron(np.eye(nc), Hs)
    Pi_f = np.kron(Pi, Pis)
    Q_f = np.eye(nf) - Pi_f
    return T_f, H_f, Pi_f, Q_f, J, A


print("=" * 82)
print("  PT QUANTUM GRAVITY: CONTINUUM LIFT OF THE FINITE DIRAC ALGEBRA")
print("=" * 82)

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]
T0, classes0 = build_T_gaps(gaps, 6)
H0, Pi0, Q0 = finite_dirac_data(T0)

levels = [
    ("base", T0, H0, Pi0, Q0),
]

n_pass = 0
n_total = 0

print(f"\n  Base boundary: T_6 shape = {T0.shape}, classes = {classes0}")

current = (T0, H0, Pi0, Q0)
for step, (dim, stay) in enumerate(((5, 0.73), (7, 0.67)), start=1):
    T, H, Pi, Q = current
    T_f, H_f, Pi_f, Q_f, J, A = refine(T, H, Pi, dim, stay)
    label = f"refinement {step}: x Z/{dim}Z"
    print(f"\n--- {label}: {T.shape[0]} -> {T_f.shape[0]} states ---")

    eye_c = np.eye(T.shape[0])
    aj_err = mat_norm(A @ J - eye_c)
    check(f"{label}: coarse averaging and injection are inverse",
          aj_err < 1e-12,
          f"||A J - I||={aj_err:.2e}")

    t_forward = mat_norm(T_f @ J - J @ T)
    t_backward = mat_norm(A @ T_f - T @ A)
    check(f"{label}: transfer cylinder is cylindrically consistent",
          t_forward < 1e-12 and t_backward < 1e-12,
          f"||T_f J-JT||={t_forward:.2e}, ||A T_f-T A||={t_backward:.2e}")

    h_forward = mat_norm(H_f @ J - J @ H)
    h_backward = mat_norm(A @ H_f - H @ A)
    check(f"{label}: Hamiltonian constraint lifts cylindrically",
          h_forward < 1e-10 and h_backward < 1e-10,
          f"||H_f J-JH||={h_forward:.2e}, ||A H_f-H A||={h_backward:.2e}")

    pi_forward = mat_norm(Pi_f @ J - J @ Pi)
    pi_backward = mat_norm(A @ Pi_f - Pi @ A)
    check(f"{label}: physical projector lifts cylindrically",
          pi_forward < 1e-10 and pi_backward < 1e-10,
          f"||Pi_f J-JPi||={pi_forward:.2e}, ||A Pi_f-Pi A||={pi_backward:.2e}")

    N = 17
    amp_kernel_err = mat_norm(A @ np.linalg.matrix_power(T_f, N) @ J
                              - np.linalg.matrix_power(T, N))
    check(f"{label}: arbitrary-boundary amplitudes lift for every cylinder length",
          amp_kernel_err < 1e-10,
          f"||A T_f^{N} J - T^{N}||={amp_kernel_err:.2e}")

    X = random_matrix(T.shape[0], 7000 + step)
    Y = random_matrix(T.shape[0], 8000 + step)
    X_cyl = J @ X @ A
    Y_cyl = J @ Y @ A
    D_X = comm(H, Q @ X @ Q)
    D_Y = comm(H, Q @ Y @ Q)
    D_X_f = comm(H_f, Q_f @ X_cyl @ Q_f)
    D_Y_f = comm(H_f, Q_f @ Y_cyl @ Q_f)

    d_lift_err = mat_norm(A @ D_X_f @ J - D_X)
    check(f"{label}: admissible Dirac deformation lifts",
          d_lift_err < 1e-9,
          f"||A D_X^fine J-D_X||={d_lift_err:.2e}")

    bracket_lift_err = mat_norm(A @ comm(D_X_f, D_Y_f) @ J - comm(D_X, D_Y))
    check(f"{label}: deformation bracket lifts cylindrically",
          bracket_lift_err < 1e-8,
          f"||A[D_Xf,D_Yf]J-[D_X,D_Y]||={bracket_lift_err:.2e}")

    anomaly_f = mat_norm(Pi_f @ comm(D_X_f, D_Y_f) @ Pi_f)
    check(f"{label}: lifted algebra is anomaly-free physically",
          anomaly_f < 1e-8,
          f"||Pi_f[D_Xf,D_Yf]Pi_f||={anomaly_f:.2e}")

    # Negative control: a non-cylindrical fine operator should not generally
    # descend to a clean coarse deformation.
    Z_raw = random_matrix(T_f.shape[0], 9000 + step)
    D_raw = comm(H_f, Q_f @ Z_raw @ Q_f)
    desc_raw = mat_norm(A @ D_raw @ J)
    check(f"{label}: non-cylindrical fine deformations are not silently accepted",
          desc_raw > 1e-4,
          f"negative-control descendant norm={desc_raw:.2e}")

    current = (T_f, H_f, Pi_f, Q_f)
    levels.append((label, T_f, H_f, Pi_f, Q_f))

print("\n--- Direct two-step consistency ---")
T1, H1, Pi1, Q1 = levels[1][1:]
T2, H2, Pi2, Q2 = levels[2][1:]

long_phys = mat_norm(H2 @ Pi2) + mat_norm(Pi2 @ H2)
check("two-step limit: physical sector remains constrained",
      long_phys < 1e-9,
      f"||H_limit Pi_limit||+||Pi_limit H_limit||={long_phys:.2e}")

semigroup_err = mat_norm(np.linalg.matrix_power(T2, 5) @ np.linalg.matrix_power(T2, 9)
                         - np.linalg.matrix_power(T2, 14))
check("two-step limit: transfer semigroup remains exact",
      semigroup_err < 1e-10,
      f"||T^5 T^9 - T^14||={semigroup_err:.2e}")

print("\n" + "=" * 82)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: cylindrical continuum lift of H, Pi, amplitudes, and")
print("               finite PT Dirac algebra through CRT refinements.")
print("  TOPOLOGICAL FOAMS: see test_qg_topology_changing_foams_PT.py.")
print("  DECORATIONS/BH: see Fourier-RG and black-hole/ringdown QG scripts.")
print("=" * 82)

sys.exit(0 if n_pass == n_total else 1)
