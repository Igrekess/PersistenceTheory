#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_covariant_boundary_amplitudes_PT.py -- PT covariant amplitudes
========================================================================

This script closes the next minimal quantum-gravity step after the canonical
Perron sector:

    W_N(phi, psi) = <phi| T_m^N |psi>

for arbitrary boundary data phi, psi on the finite sieve boundary space.

Interpretation:
  - T_m is the finite PT foam cylinder with N layers;
  - psi is incoming boundary data;
  - phi is outgoing dual boundary data;
  - gluing is matrix multiplication / resolution of identity;
  - the closed foam is Tr(T_m^N);
  - the physical infinite-cylinder amplitude is
        W_phys(phi, psi) = <phi|Pi_0|psi>,
    where Pi_0 = |1><pi| is the Perron physical projector.

What this closes:
  - arbitrary finite boundary amplitudes on the PT transfer boundary;
  - gluing / composition for arbitrary boundaries;
  - closed boundary trace as the periodic path integral;
  - projection of arbitrary boundary data onto the physical Perron sector;
  - decoupling of H-exact boundary components in the physical amplitude.

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


def amplitude(T, N, phi, psi):
    """Finite covariant cylinder amplitude <phi|T^N|psi>."""
    return np.vdot(phi, np.linalg.matrix_power(T, N) @ psi)


def physical_amplitude(Pi0, phi, psi):
    """Infinite-cylinder physical amplitude <phi|Pi_0|psi>."""
    return np.vdot(phi, Pi0 @ psi)


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


def deterministic_boundaries(n):
    rng = np.random.default_rng(20260429 + n)
    phi1 = rng.normal(size=n) + 1j * rng.normal(size=n)
    phi2 = rng.normal(size=n) + 1j * rng.normal(size=n)
    psi1 = rng.normal(size=n) + 1j * rng.normal(size=n)
    psi2 = rng.normal(size=n) + 1j * rng.normal(size=n)
    chi = rng.normal(size=n) + 1j * rng.normal(size=n)
    eta = rng.normal(size=n) + 1j * rng.normal(size=n)
    rho = rng.random(size=n)
    rho = rho / rho.sum()
    return phi1, phi2, psi1, psi2, chi, eta, rho


print("=" * 78)
print("  PT QUANTUM GRAVITY: COVARIANT AMPLITUDES WITH ARBITRARY BOUNDARIES")
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
    Pi0 = perron_projector(T)
    H = -logm(T)
    phi1, phi2, psi1, psi2, chi, eta, rho = deterministic_boundaries(n)
    ones = np.ones(n)
    basis = np.eye(n, dtype=complex)
    N = 7
    M = 11
    a = 0.37 - 0.22j
    b = -0.41 + 0.13j

    print(f"--- {label}: boundary dimension {n} ---")

    row_err = np.max(np.abs(T.sum(axis=1) - 1.0))
    check(f"{label}: stochastic transfer cylinder",
          row_err < 1e-12,
          f"max row error={row_err:.2e}")

    lhs_bra = amplitude(T, N, a * phi1 + b * phi2, psi1)
    rhs_bra = np.conj(a) * amplitude(T, N, phi1, psi1) + np.conj(b) * amplitude(T, N, phi2, psi1)
    lhs_ket = amplitude(T, N, phi1, a * psi1 + b * psi2)
    rhs_ket = a * amplitude(T, N, phi1, psi1) + b * amplitude(T, N, phi1, psi2)
    sesq_err = abs(lhs_bra - rhs_bra) + abs(lhs_ket - rhs_ket)
    check(f"{label}: arbitrary-boundary sesquilinearity",
          sesq_err < 1e-10,
          f"sesquilinearity error={sesq_err:.2e}")

    glued = sum(amplitude(T, N, phi1, basis[i]) * amplitude(T, M, basis[i], psi1)
                for i in range(n))
    direct = amplitude(T, N + M, phi1, psi1)
    glue_err = abs(glued - direct)
    check(f"{label}: gluing over complete boundary basis",
          glue_err < 1e-10,
          f"|sum_i W_N(phi,e_i)W_M(e_i,psi)-W_N+M|={glue_err:.2e}")

    closed = sum(amplitude(T, N, basis[i], basis[i]) for i in range(n))
    trace = np.trace(np.linalg.matrix_power(T, N))
    trace_err = abs(closed - trace)
    check(f"{label}: closed foam equals Tr(T^N)",
          trace_err < 1e-10,
          f"|closed-trace|={trace_err:.2e}")

    norm_amp = amplitude(T, N, rho, ones)
    check(f"{label}: normalized probabilistic boundary is conserved",
          abs(norm_amp - 1.0) < 1e-12,
          f"<rho|T^N|1>={norm_amp:.12f}")

    fixed_err = mat_norm(T @ Pi0 - Pi0) + mat_norm(Pi0 @ T - Pi0) + mat_norm(Pi0 @ Pi0 - Pi0)
    check(f"{label}: Perron boundary projector is a fixed cylinder",
          fixed_err < 1e-10,
          f"||T Pi-Pi||+||Pi T-Pi||+||Pi^2-Pi||={fixed_err:.2e}")

    largeN = 80
    limit_err = mat_norm(np.linalg.matrix_power(T, largeN) - Pi0)
    amp_limit_err = abs(amplitude(T, largeN, phi1, psi1) - physical_amplitude(Pi0, phi1, psi1))
    check(f"{label}: long cylinder converges to physical amplitude",
          limit_err < 1e-10 and amp_limit_err < 1e-10,
          f"||T^{largeN}-Pi0||={limit_err:.2e}, |W_N-W_phys|={amp_limit_err:.2e}")

    ket_exact = H @ chi
    bra_exact = H.conj().T @ eta
    exact_ket_amp = physical_amplitude(Pi0, phi1, ket_exact)
    exact_bra_amp = physical_amplitude(Pi0, bra_exact, psi1)
    exact_err = abs(exact_ket_amp) + abs(exact_bra_amp)
    check(f"{label}: H-exact boundary components decouple physically",
          exact_err < 1e-10,
          f"|W_phys(phi,Hchi)|+|W_phys(H^dag eta,psi)|={exact_err:.2e}")

    projected_finite = amplitude(T, N, Pi0.conj().T @ phi1, Pi0 @ psi1)
    projected_physical = physical_amplitude(Pi0, phi1, psi1)
    proj_amp_err = abs(projected_finite - projected_physical)
    check(f"{label}: projected boundaries are N-independent",
          proj_amp_err < 1e-10,
          f"|W_N(Pi phi,Pi psi)-W_phys|={proj_amp_err:.2e}")

print("\n" + "=" * 78)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: arbitrary finite boundary amplitudes, gluing, closed trace,")
print("               physical Perron projection, H-exact boundary decoupling.")
print("  DIRAC LIFT: see test_qg_continuum_dirac_lift_PT.py.")
print("  TOPOLOGICAL FOAMS: see test_qg_topology_changing_foams_PT.py.")
print("  DECORATIONS/BH: see Fourier-RG and black-hole/ringdown QG scripts.")
print("=" * 78)

sys.exit(0 if n_pass == n_total else 1)
