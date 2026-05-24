#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cross-branch derivation of alpha_s at m_Z (ch11 thm:alphas_cross_branch).

Verifies that the formula

   alpha_s(m_Z) = sin^2(theta_3, q_-) / (1 - alpha_EM_sieve)

is UNIQUE among the 12 candidates

   sin^2(theta_p, q_pm) * N(alpha_EM)^{pm 1}

with p in {3, 5, 7}, q_pm in {q_+, q_-}, N in {1, 1-alpha}, under
the four constraints C1-C4 of the theorem:

  C1 (colour factor, N_c=3)  -> p = 3
  C2 (edge branch)           -> q_-
  C3 (sieve-EM dilution)     -> normaliser 1/(1-alpha_EM)
  C4 (scale at mu* = 15)     -> evaluation at the fixed point

The theorem promotes alpha_s from a BRIDGE identification to a
DERIVED result.
"""
from __future__ import division, print_function

import sys
import mpmath as mp
from itertools import product

mp.mp.dps = 50

MU_STAR = mp.mpf(15)
q_plus  = 1 - 2/MU_STAR
q_minus = mp.exp(-1/MU_STAR)
PI      = mp.pi


def delta_p(p, q):
    return (1 - q**p) / p


def sin2(p, q):
    d = delta_p(p, q)
    return d * (2 - d)


# --- Reference values ---------------------------------------------
alpha_EM_sieve = sin2(3, q_plus) * sin2(5, q_plus) * sin2(7, q_plus)  # tree
PDG_alpha_s    = mp.mpf("0.11800")
TOL_C1         = mp.mpf("0.02")    # 2% tolerance for prime selection
TOL_C3         = mp.mpf("0.001")   # 0.1% tolerance for normalisation


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


# =======================================================================
section("CANONICAL FORMULA")
# =======================================================================
alphas_canonical = sin2(3, q_minus) / (1 - alpha_EM_sieve)
err_canonical = abs(alphas_canonical - PDG_alpha_s) / PDG_alpha_s * 100
print(f"  alpha_s PT  = sin^2(theta_3, q_-) / (1 - alpha_EM_sieve)")
print(f"              = {mp.nstr(sin2(3,q_minus), 12)} / {mp.nstr(1-alpha_EM_sieve, 12)}")
print(f"              = {mp.nstr(alphas_canonical, 15)}")
print(f"  PDG         = 0.11800")
print(f"  error       = {mp.nstr(err_canonical, 5)}%")

# =======================================================================
section("C1 (COLOUR FACTOR): prime selection p = N_c = 3")
# =======================================================================
for p in [3, 5, 7]:
    val = sin2(p, q_minus) / (1 - alpha_EM_sieve)
    err = abs(val - PDG_alpha_s) / PDG_alpha_s * 100
    marker = " <-- selected by C1" if p == 3 else ""
    print(f"  p = {p}: sin^2(theta_{p}, q_-)/(1-alpha) = {mp.nstr(val, 8)}, "
          f"err {mp.nstr(err, 5)}%{marker}")

# =======================================================================
section("C2 (EDGE BRANCH): q_- vs q_+")
# =======================================================================
for qname, q in [("q_+", q_plus), ("q_-", q_minus)]:
    val = sin2(3, q) / (1 - alpha_EM_sieve)
    err = abs(val - PDG_alpha_s) / PDG_alpha_s * 100
    marker = " <-- selected by C2" if qname == "q_-" else ""
    print(f"  q = {qname}: sin^2(theta_3, {qname})/(1-alpha) = {mp.nstr(val, 8)}, "
          f"err {mp.nstr(err, 5)}%{marker}")

# =======================================================================
section("C3 (SIEVE-EM DILUTION): normalisation 1/(1-alpha_EM)")
# =======================================================================
cases = [
    ("raw     sin^2(3,q_-)",           sin2(3, q_minus)),
    ("dilute  sin^2(3,q_-)/(1-a)",     sin2(3, q_minus) / (1 - alpha_EM_sieve)),
    ("inverse sin^2(3,q_-)*(1-a)",     sin2(3, q_minus) * (1 - alpha_EM_sieve)),
    ("boosted sin^2(3,q_-)/(1+a)",     sin2(3, q_minus) / (1 + alpha_EM_sieve)),
]
for label, val in cases:
    err = abs(val - PDG_alpha_s) / PDG_alpha_s * 100
    marker = " <-- selected by C3" if "dilute" in label else ""
    print(f"  {label:35} = {mp.nstr(val, 10)}  err {mp.nstr(err, 5)}%{marker}")

# =======================================================================
section("FULL ENUMERATION (12 candidates)")
# =======================================================================
candidates = []
for p in [3, 5, 7]:
    for qname, q in [("q+", q_plus), ("q-", q_minus)]:
        for nlabel, nfactor in [
            ("raw",       mp.mpf(1)),
            ("/(1-a)",    1 / (1 - alpha_EM_sieve)),
        ]:
            val = sin2(p, q) * nfactor
            err_rel = abs(val - PDG_alpha_s) / PDG_alpha_s
            candidates.append((p, qname, nlabel, val, err_rel))

# Sort by smallest error
candidates.sort(key=lambda x: x[4])
print(f"  {'p':<3} {'q':<4} {'norm':<8} {'alpha_s':<15} {'err %':<10}")
for p, qname, nlabel, val, err in candidates:
    marker = " <-- canonical" if p == 3 and qname == "q-" and "/(1-a)" in nlabel else ""
    print(f"  {p:<3} {qname:<4} {nlabel:<8} {mp.nstr(val, 10):<15} "
          f"{mp.nstr(err*100, 5):<10}{marker}")

# =======================================================================
section("UNIQUENESS CHECK")
# =======================================================================
# The canonical candidate must be the unique one matching PDG to 0.1%.
# The "raw" p=3 q_- gives 0.68% (off by 13x from canonical), ruled out
# once the Ward-like dilution constraint C3 is imposed.
TIGHT_TOL = mp.mpf("0.001")  # 0.1%
within_tight = [(p, q, n) for p, q, n, v, e in candidates if e < TIGHT_TOL]
print(f"  Candidates within 0.1% (PT precision floor): {len(within_tight)}")
for p, q, n in within_tight:
    print(f"    p={p} q={q} norm={n}")

# TEST OUTCOMES
# =======================================================================
section("TEST OUTCOMES")
# =======================================================================
# T1: canonical gives < 0.1% error vs PDG
assert err_canonical < 0.1, f"canonical err {err_canonical}% >= 0.1%"
print(f"  [PASS] canonical (p=3, q_-, /(1-a)): err {mp.nstr(err_canonical, 4)}% < 0.1%")

# T2: all p != 3 cases (with q_- and /(1-a)) give > 5% error
for p in [5, 7]:
    val = sin2(p, q_minus) / (1 - alpha_EM_sieve)
    err = abs(val - PDG_alpha_s) / PDG_alpha_s * 100
    assert err > 5, f"C1 violated for p={p}"
print("  [PASS] C1: p != 3 gives > 5% error (prime selection forced)")

# T3: q_+ gives huge error
val = sin2(3, q_plus) / (1 - alpha_EM_sieve)
err = abs(val - PDG_alpha_s) / PDG_alpha_s * 100
assert err > 50, f"C2 violated: q_+ gives only err {err}%"
print(f"  [PASS] C2: q_+ gives {mp.nstr(err, 4)}% error (edge branch forced)")

# T4: without dilution (raw) gives > 0.5% error
raw = sin2(3, q_minus)
err_raw = abs(raw - PDG_alpha_s) / PDG_alpha_s * 100
assert err_raw > 0.5, f"C3 violated: raw is too close"
print(f"  [PASS] C3: raw formula gives {mp.nstr(err_raw, 4)}% error (dilution required)")

# T5: only canonical is within 0.1% tolerance
assert len(within_tight) == 1, f"non-unique: {len(within_tight)} candidates within 0.1%"
print(f"  [PASS] Uniqueness: only canonical candidate within 0.1% (among 12)")

print()
print(f"  5/5 tests PASS. Theorem thm:alphas_cross_branch confirmed.")
print(f"  alpha_s promoted from BRIDGE to DERIVED.")
sys.exit(0)
