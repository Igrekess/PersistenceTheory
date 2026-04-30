#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_canonical_constraints_PT.py -- Minimal canonical QG sector in PT
========================================================================

This script tests the part of the quantum-gravity programme that can be
closed without importing an external canonical formalism.

The PT input is the stochastic transfer matrix T_m of sieve gaps.  From it
we define the canonical Hamiltonian constraint

    H_m = -log(T_m).

Because T_m is stochastic, the Perron-Frobenius mode has lambda_0 = 1, hence
H_m annihilates the physical state.  The physical projector is the rank-one
Perron projector

    Pi_0 = |1><pi|,

where |1> is the right PF vector and <pi| is the stationary left measure.

What this closes:
  - a minimal physical Hilbert sector ker(H_m);
  - the WDW constraint H_m |Psi_0> = 0;
  - a canonical physical projector Pi_0;
  - spectral Dirac observables O=f(T_m), commuting with H_m;
  - Z_N = Tr(T_m^N) dominated by dim ker(H_m)=1.

What this does NOT close:
  - the continuum lift of the finite PT Dirac algebra
    (see test_qg_continuum_dirac_lift_PT.py);
  - topology-changing physical foams
    (see test_qg_topology_changing_foams_PT.py);
  - higher Fourier/RG decorations
    (see test_qg_fourier_rg_decorations_PT.py);
  - black-hole/ringdown precision skeletons
    (see test_qg_black_hole_ringdown_precision_PT.py).

Status: canonical minimal sector, not yet full quantum gravity.
"""

import io
import math
import sys

import numpy as np
from scipy.linalg import logm

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def sieve_primes(n_max):
    """Return primes up to n_max by Eratosthenes sieve."""
    is_prime = [True] * (n_max + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n_max**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n_max + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n_max + 1) if is_prime[i]]


def build_T_gaps(gaps, modulus):
    """Build the row-stochastic transition matrix of prime gaps mod modulus."""
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
    """Return the stationary left probability vector pi with pi T = pi."""
    eigvals, eigvecs = np.linalg.eig(T.T)
    idx = np.argmin(np.abs(eigvals - 1.0))
    pi = np.real(eigvecs[:, idx])
    if pi.sum() < 0:
        pi = -pi
    pi = np.maximum(pi, 0.0)
    total = pi.sum()
    if total == 0:
        raise RuntimeError("stationary vector is numerically degenerate")
    return pi / total


def canonical_data(T):
    """Return H, right PF vector, left PF vector, and Perron projector."""
    n = T.shape[0]
    right = np.ones(n)
    left = stationary_left(T)
    projector = np.outer(right, left)  # left @ right = 1
    H = -logm(T)
    return H, right, left, projector


def mat_norm(A):
    return float(np.linalg.norm(A))


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


print("=" * 76)
print("  PT QUANTUM GRAVITY: MINIMAL CANONICAL CONSTRAINT SECTOR")
print("=" * 76)

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]

T6, classes6 = build_T_gaps(gaps, 6)
T30, classes30 = build_T_gaps(gaps, 30)

print(f"\n  T_6  shape = {T6.shape}, classes = {classes6}")
print(f"  T_30 shape = {T30.shape}, classes = {classes30}")
print()

n_pass = 0
n_total = 0

# ---------------------------------------------------------------------
# Stochastic input and PF uniqueness
# ---------------------------------------------------------------------
row_err6 = np.max(np.abs(T6.sum(axis=1) - 1.0))
row_err30 = np.max(np.abs(T30.sum(axis=1) - 1.0))
check("T_m is row-stochastic",
      row_err6 < 1e-12 and row_err30 < 1e-12,
      f"max row error: T6={row_err6:.2e}, T30={row_err30:.2e}")

for label, T in (("T6", T6), ("T30", T30)):
    eigvals = np.linalg.eigvals(T)
    eig_abs = np.sort(np.abs(eigvals))[::-1]
    pf_err = abs(eig_abs[0] - 1.0)
    gap = 1.0 - eig_abs[1]
    check(f"{label}: unique Perron mode and spectral gap",
          pf_err < 1e-12 and gap > 0,
          f"|lambda0|-1={pf_err:.2e}, Delta=1-|lambda1|={gap:.4f}")

# ---------------------------------------------------------------------
# Hamiltonian constraint and physical projector
# ---------------------------------------------------------------------
for label, T in (("T6", T6), ("T30", T30)):
    H, right, left, Pi0 = canonical_data(T)
    h_right = mat_norm(np.real(H @ right))
    left_h = mat_norm(np.real(left @ H))
    proj_err = mat_norm(Pi0 @ Pi0 - Pi0)
    hp_err = mat_norm(np.real(H @ Pi0))
    ph_err = mat_norm(np.real(Pi0 @ H))
    check(f"{label}: WDW constraint H|Psi0>=0 and <pi|H=0",
          h_right < 1e-10 and left_h < 1e-10,
          f"||H1||={h_right:.2e}, ||pi H||={left_h:.2e}")
    check(f"{label}: Perron projector is physical",
          proj_err < 1e-12 and hp_err < 1e-10 and ph_err < 1e-10,
          f"||Pi0^2-Pi0||={proj_err:.2e}, ||H Pi0||={hp_err:.2e}, ||Pi0 H||={ph_err:.2e}")

# ---------------------------------------------------------------------
# First-class minimal algebra and Dirac observables
# ---------------------------------------------------------------------
for label, T in (("T6", T6), ("T30", T30)):
    H, right, left, Pi0 = canonical_data(T)
    comm_HT = mat_norm(np.real(H @ T - T @ H))
    O = np.linalg.matrix_power(T, 7)
    comm_HO = mat_norm(np.real(H @ O - O @ H))
    dirac_err = mat_norm(np.real(Pi0 @ O @ Pi0 - Pi0))
    check(f"{label}: spectral observables commute with the constraint",
          comm_HT < 1e-10 and comm_HO < 1e-10 and dirac_err < 1e-10,
          f"||[H,T]||={comm_HT:.2e}, ||[H,T^7]||={comm_HO:.2e}, ||Pi O Pi-Pi||={dirac_err:.2e}")

# ---------------------------------------------------------------------
# Partition function / path integral sector
# ---------------------------------------------------------------------
for label, T in (("T6", T6), ("T30", T30)):
    eigvals = np.linalg.eigvals(T)
    for N in (10, 50, 200):
        trace_power = np.trace(np.linalg.matrix_power(T, N))
        spectral_sum = np.sum(eigvals**N)
        if N == 50:
            trace_err = abs(trace_power - spectral_sum)
            check(f"{label}: Z_N=Tr(T^N)=sum lambda^N",
                  trace_err < 1e-8,
                  f"N={N}, |trace-spectral|={trace_err:.2e}")
    Z200 = float(np.real(np.trace(np.linalg.matrix_power(T, 200))))
    check(f"{label}: path integral projects onto dim ker(H)=1",
          abs(Z200 - 1.0) < 1e-8,
          f"Tr(T^200)={Z200:.12f}")

# ---------------------------------------------------------------------
# Excitation spectrum: zero mode plus damped modes
# ---------------------------------------------------------------------
for label, T in (("T6", T6), ("T30", T30)):
    eigvals = np.linalg.eigvals(T)
    energies = []
    for ev in eigvals:
        if abs(ev - 1.0) > 1e-10:
            energies.append(-np.log(complex(ev)))
    min_re = min(e.real for e in energies)
    has_complex = any(abs(e.imag) > 1e-10 for e in energies)
    check(f"{label}: excited modes have positive damping",
          min_re > 0,
          f"min Re(E_excited)={min_re:.4f}, complex modes={has_complex}")

print("\n" + "=" * 76)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: minimal canonical sector ker(H), Pi0, Dirac spectral observables.")
print("  NEXT CLOSED SECTORS: continuum Dirac lift and topology-changing physical foams.")
print("  LATER CLOSED SECTORS: Fourier/RG decorations and BH/ringdown skeleton.")
print("=" * 76)

sys.exit(0 if n_pass == n_total else 1)
