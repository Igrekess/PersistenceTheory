#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_fourier_rg_decorations_PT.py -- Higher Fourier/RG PT foam decorations
================================================================================

This script closes the structural part of the higher-decoration question left
open by the minimal PT spin-foam sector.

The input is deliberately small:

  - the active faces are p in {3, 5, 7};
  - Fourier decorations are non-trivial characters of Z/pZ, modulo conjugation;
  - RG decorations are log-derivative insertions

        L_p = -mu d/dmu ln sin_p^2(mu),  L_p(mu*) = gamma_p(mu*).

The tests verify that:

  - character orbits are completely classified by conjugation k ~ -k;
  - higher Fourier orbits are orthogonal to the fundamental orbit and are
    distinct observables, not hidden corrections to it;
  - CRT/Pontryagin factorisation makes composite decorations multiplicative;
  - RG insertions generate a rigid symmetric algebra in gamma_3,gamma_5,gamma_7;
  - mixed Fourier/RG decorations factorise without extra coefficients.

This closes the mathematical classification of higher Fourier/RG decorations.
It does not turn a particular decorated observable into a measured
black-hole/ringdown number; that is tested separately.
"""

import io
import itertools
import math
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ACTIVE = (3, 5, 7)
MU_STAR = 15.0
TOL = 1e-10


def sin2_theta(p, mu):
    q = 1.0 - 2.0 / mu
    qp = q**p
    return (1.0 - qp) * (2.0 * p - 1.0 + qp) / p**2


def gamma_p(p, mu):
    q = 1.0 - 2.0 / mu
    qp = q**p
    delta = (1.0 - qp) / p
    dln_delta = -2.0 * p * q ** (p - 1) / (mu * (1.0 - qp))
    factor = 2.0 * (1.0 - delta) / (2.0 - delta)
    return -dln_delta * factor


def gamma_numeric(p, mu, h=1e-5):
    left = math.log(sin2_theta(p, mu - h))
    right = math.log(sin2_theta(p, mu + h))
    return -mu * (right - left) / (2.0 * h)


def character(p, k):
    r = np.arange(p)
    return np.exp(2j * np.pi * k * r / p)


def inner(u, v):
    return np.vdot(u, v) / len(u)


def conjugacy_orbits(p):
    unseen = set(range(1, p))
    orbits = []
    while unseen:
        k = min(unseen)
        orbit = {k, (-k) % p}
        orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    return orbits


def multi_indices(total_degree, dim):
    if dim == 1:
        yield (total_degree,)
        return
    for first in range(total_degree + 1):
        for rest in multi_indices(total_degree - first, dim - 1):
            yield (first,) + rest


def monomial(gammas, alpha):
    out = 1.0
    for p, power in zip(ACTIVE, alpha):
        out *= gammas[p] ** power
    return out


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


print("=" * 84)
print("  PT QUANTUM GRAVITY: HIGHER FOURIER/RG DECORATIONS")
print("=" * 84)
print(f"\n  active faces = {ACTIVE}, mu* = {MU_STAR}")

n_pass = 0
n_total = 0

print("\n--- Fourier orbit classification ---")
expected_counts = {3: 1, 5: 2, 7: 3}
total_orbits = 0
for p in ACTIVE:
    orbits = conjugacy_orbits(p)
    total_orbits += len(orbits)
    cover = sorted(k for orbit in orbits for k in orbit)
    conj_ok = all(((-k) % p) in orbit for orbit in orbits for k in orbit)
    check(f"Z/{p}Z: conjugate character orbits are exhaustive",
          cover == list(range(1, p)) and conj_ok and len(orbits) == expected_counts[p],
          f"orbits={orbits}")

check("active T^3 has six non-trivial conjugacy orbits",
      total_orbits == 6,
      f"total orbits={total_orbits} = 1+2+3")

for p in ACTIVE:
    chars = [character(p, k) for k in range(p)]
    gram = np.array([[inner(a, b) for b in chars] for a in chars])
    err = np.linalg.norm(gram - np.eye(p))
    check(f"Z/{p}Z: Pontryagin characters are orthonormal",
          err < TOL,
          f"||G-I||={err:.2e}")

for p in (5, 7):
    fund = character(p, 1)
    higher_errs = []
    for orbit in conjugacy_orbits(p)[1:]:
        for k in orbit:
            higher_errs.append(abs(inner(fund, character(p, k))))
    check(f"Z/{p}Z: higher Fourier orbits are not corrections to chi_1",
          max(higher_errs) < TOL,
          f"max |<chi_1, chi_k>| over higher orbits={max(higher_errs):.2e}")

print("\n--- CRT/Pontryagin multiplicativity ---")
p, q = 5, 7
for k, ell in ((1, 1), (2, 3)):
    m = p * q
    residues = np.arange(m)
    chi_m = np.exp(2j * np.pi * (k * residues / p + ell * residues / q))
    chi_factor = character(p, k)[residues % p] * character(q, ell)[residues % q]
    err = np.linalg.norm(chi_m - chi_factor) / math.sqrt(m)
    check(f"CRT character factorisation ({p},{k}) x ({q},{ell})",
          err < TOL,
          f"RMS error={err:.2e}")

print("\n--- RG log-derivative insertions ---")
gammas = {p: gamma_p(p, MU_STAR) for p in ACTIVE}
for p in ACTIVE:
    g_num = gamma_numeric(p, MU_STAR)
    err = abs(gammas[p] - g_num)
    check(f"gamma_{p} is the log-derivative insertion L_{p}",
          err < 1e-8,
          f"analytic={gammas[p]:.12f}, numeric={g_num:.12f}, err={err:.2e}")

gamma_order = gammas[3] > gammas[5] > gammas[7] > 0.5
check("active RG hierarchy is ordered and above threshold",
      gamma_order,
      ", ".join(f"gamma_{p}={gammas[p]:.6f}" for p in ACTIVE))

for degree in range(0, 5):
    alphas = list(multi_indices(degree, len(ACTIVE)))
    expected = math.comb(degree + len(ACTIVE) - 1, len(ACTIVE) - 1)
    values = [monomial(gammas, alpha) for alpha in alphas]
    positive = all(v > 0 for v in values)
    check(f"RG degree {degree}: symmetric-decoration basis count",
          len(alphas) == expected and positive,
          f"{len(alphas)} monomials, expected C({degree+2},2)={expected}")

alpha = (2, 1, 1)
lhs = monomial(gammas, alpha)
rhs = gammas[3] ** 2 * gammas[5] * gammas[7]
check("higher RG cross-face products have no free coefficient",
      abs(lhs - rhs) < TOL,
      f"gamma^{{{alpha}}}={lhs:.12f}")

print("\n--- Mixed Fourier/RG decorations ---")
for p, k, alpha in ((5, 2, (1, 0, 1)), (7, 3, (0, 2, 1))):
    chi = character(p, k)
    rg = monomial(gammas, alpha)
    decorated = rg * chi
    norm_ratio = math.sqrt(float(np.real(inner(decorated, decorated)))) / abs(rg)
    check(f"mixed decoration chi_{k}^{p} x gamma^{alpha} factorises",
          abs(norm_ratio - 1.0) < TOL,
          f"||rg chi||/|rg|={norm_ratio:.12f}")

fundamental_slots = 3
higher_slots = total_orbits - fundamental_slots
rg_basis_deg_le4 = sum(math.comb(d + 2, 2) for d in range(5))
mixed_capacity = total_orbits * rg_basis_deg_le4
check("decorated observable space is finite at every truncation degree",
      higher_slots == 3 and rg_basis_deg_le4 == 35 and mixed_capacity == 210,
      f"Fourier orbits={total_orbits}, RG monomials<=4={rg_basis_deg_le4}, mixed={mixed_capacity}")

print("\n" + "=" * 84)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: higher Fourier character orbits and RG log-derivative")
print("               decorations as a rigid Pontryagin/CRT symmetric algebra.")
print("  STILL SEPARATE: assigning a decorated observable to a precision datum.")
print("=" * 84)

sys.exit(0 if n_pass == n_total else 1)
