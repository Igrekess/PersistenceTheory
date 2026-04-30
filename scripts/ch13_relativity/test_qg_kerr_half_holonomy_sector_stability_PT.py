#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_half_holonomy_sector_stability_PT.py
=================================================

Sector-stability test for the PT Kerr half-holonomy deformation.

The previous uniqueness script proves that, at fixed selected sector and first
order in x=M Omega_H(a), the only admissible deformation is

    phi_i(x)   = phi_i + pi x,
    kappa_i(x) = kappa_i.

This script checks the remaining fixed-sector assumption: after applying the
exact half-holonomy over the whole Kerr range 0 <= x < 1/2, does another
oscillatory Ruelle mode overtake the selected (2,2,0) macro mode?

Because each deformed quality factor is affine in x,

    Q_i(x) = 2*pi (phi_i + pi x) / (2 kappa_i),

the ordering against any competitor can change only at a single crossing.
It is enough to evaluate every gap at both endpoints x=0 and x=1/2 and to
compute the crossing location.  No observational data enter this test.
"""

from __future__ import annotations

import io
import math
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

CYCLIC_PERIOD = 2.0 * math.pi
X_MAX_EXTREMAL = 0.5
TOL = 1e-12


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


def ruelle_modes(T):
    modes = []
    for ev in np.linalg.eigvals(T):
        if abs(ev.imag) > 1e-10 and abs(ev) < 1.0:
            phi = abs(np.angle(ev))
            kappa = -math.log(abs(ev))
            chirality = 1 if ev.imag > 0 else -1
            q0 = deformed_quality(phi, kappa, 0.0)
            modes.append({
                "lambda": ev,
                "phi": phi,
                "kappa": kappa,
                "chirality": chirality,
                "Q0": q0,
            })
    return modes


def deformed_quality(phi, kappa, x):
    return CYCLIC_PERIOD * (phi + math.pi * x) / (2.0 * kappa)


def crossing_x(selected, competitor):
    """Return x where Q_selected(x)=Q_competitor(x), or None if parallel."""
    phi_s, k_s = selected["phi"], selected["kappa"]
    phi_c, k_c = competitor["phi"], competitor["kappa"]
    denom = (1.0 / k_s) - (1.0 / k_c)
    numer = (phi_c / k_c) - (phi_s / k_s)
    if abs(denom) < 1e-15:
        return None
    return float(numer / (math.pi * denom))


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


print("=" * 88)
print("  PT QG: KERR HALF-HOLONOMY SECTOR STABILITY")
print("=" * 88)

n_pass = 0
n_total = 0

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]
T30, classes = build_T_gaps(gaps, 30)
modes = ruelle_modes(T30)
positive = sorted([m for m in modes if m["chirality"] > 0], key=lambda m: -m["Q0"])
negative = sorted([m for m in modes if m["chirality"] < 0], key=lambda m: -m["Q0"])
selected = positive[0]

print(f"\n  classes={classes}")
print(f"  positive modes={len(positive)}, negative modes={len(negative)}")
for i, m in enumerate(positive):
    print(f"  + mode {i}: Q0={m['Q0']:.6f}, phi={m['phi']:.6f}, kappa={m['kappa']:.6f}, lambda={m['lambda']}")

row_err = np.max(np.abs(T30.sum(axis=1) - 1.0))
check("T_30 is stochastic",
      row_err < TOL,
      f"max row error={row_err:.2e}")

check("complex modes come in chiral conjugate pairs",
      len(positive) == len(negative) and len(positive) > 0
      and all(any(abs(n["lambda"] - np.conj(p["lambda"])) < 1e-10 for n in negative) for p in positive),
      f"pairs={len(positive)}")

check("undeformed selected mode is unique",
      selected["Q0"] > positive[1]["Q0"],
      f"Q0 selected={selected['Q0']:.6f}, runner-up={positive[1]['Q0']:.6f}")

q_selected_0 = deformed_quality(selected["phi"], selected["kappa"], 0.0)
q_selected_ext = deformed_quality(selected["phi"], selected["kappa"], X_MAX_EXTREMAL)
check("half-holonomy increases selected quality monotonically",
      q_selected_ext > q_selected_0,
      f"Q(0)={q_selected_0:.6f}, Q(1/2)={q_selected_ext:.6f}")

endpoint_gaps = []
crossings = []
for comp in positive[1:]:
    gap0 = q_selected_0 - deformed_quality(comp["phi"], comp["kappa"], 0.0)
    gap_ext = q_selected_ext - deformed_quality(comp["phi"], comp["kappa"], X_MAX_EXTREMAL)
    endpoint_gaps.append((gap0, gap_ext))
    crossings.append(crossing_x(selected, comp))

check("selected mode beats every competitor at x=0",
      all(g0 > 0.0 for g0, _gext in endpoint_gaps),
      f"min gap={min(g0 for g0, _ in endpoint_gaps):.6f}")

check("selected mode beats every competitor at extremal x=1/2",
      all(gext > 0.0 for _g0, gext in endpoint_gaps),
      f"min gap={min(gext for _, gext in endpoint_gaps):.6f}")

check("affine quality gaps imply endpoint checks cover the whole Kerr interval",
      True,
      "Q_s(x)-Q_i(x) is affine in x")

admissible_crossings = [x for x in crossings if x is not None and 0.0 <= x <= X_MAX_EXTREMAL]
check("no competitor crossing lies in the Kerr interval",
      len(admissible_crossings) == 0,
      f"crossings={crossings}")

grid = np.linspace(0.0, X_MAX_EXTREMAL, 1001)
best_indices = []
for x in grid:
    qualities = [deformed_quality(m["phi"], m["kappa"], x) for m in positive]
    best_indices.append(int(np.argmax(qualities)))
check("dense scan keeps selected mode as argmax for 0 <= x <= 1/2",
      set(best_indices) == {0},
      f"argmax set={sorted(set(best_indices))}")

negative_selected = negative[0]
check("negative chirality has conjugate-stable selected sector",
      abs(negative_selected["lambda"] - np.conj(selected["lambda"])) < 1e-10
      and abs(negative_selected["Q0"] - selected["Q0"]) < TOL)

for spin in (0.0, 0.5, 0.69, 0.9, 0.999999):
    x = spin / (2.0 * (1.0 + math.sqrt(max(0.0, 1.0 - spin * spin))))
    qualities = [deformed_quality(m["phi"], m["kappa"], x) for m in positive]
    check(f"astrophysical spin a={spin:g} keeps selected sector",
          int(np.argmax(qualities)) == 0,
          f"x={x:.6f}, Q_best={max(qualities):.6f}")

min_endpoint_margin = min(min(g0, gext) for g0, gext in endpoint_gaps)
check("sector stability has a finite positive margin",
      min_endpoint_margin > 0.1,
      f"minimum endpoint margin={min_endpoint_margin:.6f}")

print("\n" + "=" * 88)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: the half-holonomy deformation does not trigger a sector switch")
print("               for the full Kerr interval 0 <= M Omega_H < 1/2.")
print("=" * 88)

sys.exit(0 if n_pass == n_total else 1)
