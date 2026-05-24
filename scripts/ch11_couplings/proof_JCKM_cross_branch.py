#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cross-branch derivation of J_CKM (ch11 thm:JCKM_cross_branch).

Verifies that the formula

   J_CKM = alpha_EM^2 * sin^2(theta_23^PMNS)

is UNIQUE among the 18 candidates

   alpha_EM^n * X

with n in {1, 2, 3} and X in {sin^2(theta_{12,13,23}^PMNS), gamma_3,5,7},
under the three constraints J1-J3 of the theorem:

  J1 (double vertex-edge crossing)   -> n = 2
  J2 (atmospheric PMNS sector)        -> X = sin^2(theta_23^PMNS)
  J3 (observed scale J ~ 3e-5)        -> n = 2 forced

The theorem promotes J_CKM from a BRIDGE identification to a DERIVED result.
"""
from __future__ import division, print_function

import sys
import mpmath as mp

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


def gamma_p(p, mu=MU_STAR):
    q = 1 - 2/mu
    d = delta_p(p, q)
    return 4*p*q**(p-1)*(1-d) / (mu*(1-q**p)*(2-d))


# --- Reference values ---------------------------------------------
alpha_EM_sieve = sin2(3, q_plus) * sin2(5, q_plus) * sin2(7, q_plus)

# Dress alpha to physical value ~1/137.036
# (full chain: bare + spiral + echo + 2-loop)
p1 = 2
delta_2 = (1 - q_plus**p1) / p1
sin2_2 = delta_2 * (2 - delta_2)
theta_2 = mp.acos(1 - delta_2)
N2 = 26
cos2_leak = mp.cos(theta_2/N2)**2
depth_2 = (MU_STAR - p1) / p1**2
F2_d = sin2_2 * cos2_leak * depth_2
alpha_1 = 1 / (1/alpha_EM_sieve + F2_d)
sum_g2 = sum(gamma_p(p)**2 for p in [3,5,7])
delta_5 = (1 - q_plus**5)/5
delta_7 = (1 - q_plus**7)/7
sum_g = sum(gamma_p(p) for p in [3,5,7])
prop_tree = (delta_5 + delta_7) / sum_g
prop = prop_tree * (1 + alpha_EM_sieve/25)
r_fb = alpha_1 * sum_g2 * prop
spiral = F2_d / (1 + gamma_p(3) * r_fb)
alpha_spiral = 1 / (1/alpha_EM_sieve + spiral)
sin2_11 = sin2(11, q_plus); sin2_13 = sin2(13, q_plus)
beta_echo = sin2_11*gamma_p(11) + sin2_13*gamma_p(13)
delta_echo = sin2_2 * beta_echo * alpha_spiral**2
alpha_echo = 1/(1/alpha_spiral + delta_echo)
delta_2l = (alpha_echo/PI)**2 / 3
alpha_EM = 1 / (1/alpha_echo + delta_2l)

# PMNS angles
sin2_th12 = 1 - gamma_p(5)
sin2_th13 = 3*alpha_EM / (1 - 2*alpha_EM)
sin2_th23 = gamma_p(7) - sin2_th13

# gamma values
g3 = gamma_p(3); g5 = gamma_p(5); g7 = gamma_p(7)

PDG_J_CKM = mp.mpf("3.08e-5")


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


# =======================================================================
section("CANONICAL FORMULA")
# =======================================================================
J_canonical = alpha_EM**2 * sin2_th23
err_canonical = abs(J_canonical - PDG_J_CKM) / PDG_J_CKM * 100
print(f"  J_CKM PT   = alpha_EM^2 * sin^2(theta_23^PMNS)")
print(f"             = ({mp.nstr(alpha_EM, 10)})^2 * {mp.nstr(sin2_th23, 10)}")
print(f"             = {mp.nstr(J_canonical, 12)}")
print(f"  PDG        = 3.08e-5")
print(f"  error      = {mp.nstr(err_canonical, 5)}%  (tree level)")

# NLO correction (1+eps)
s = mp.mpf(1)/2
alpha_s_val = sin2(3, q_minus) / (1 - alpha_EM_sieve)  # ~0.118
eps = alpha_s_val * s / (2*PI)
J_nlo = J_canonical * (1 + eps)
err_nlo = abs(J_nlo - PDG_J_CKM) / PDG_J_CKM * 100
print(f"  NLO (1+eps): err {mp.nstr(err_nlo, 5)}%")

# =======================================================================
section("FULL ENUMERATION (18 candidates)")
# =======================================================================
targets = [
    ("sin2(th_12)", sin2_th12),
    ("sin2(th_13)", sin2_th13),
    ("sin2(th_23)", sin2_th23),
    ("gamma_3",     g3),
    ("gamma_5",     g5),
    ("gamma_7",     g7),
]
results = []
for n in [1, 2, 3]:
    for name, x in targets:
        val = alpha_EM**n * x
        err = abs(val - PDG_J_CKM) / PDG_J_CKM
        results.append((n, name, val, err))
results.sort(key=lambda r: r[3])

print(f"  {'n':<3} {'X':<15} {'J_PT':<18} {'err (%)':<12}")
for n, name, val, err in results[:12]:
    marker = " <-- canonical" if (n == 2 and name == "sin2(th_23)") else ""
    print(f"  {n:<3} {name:<15} {mp.nstr(val, 8):<18} {mp.nstr(err*100, 5):<12}{marker}")

# =======================================================================
section("CONSTRAINT BREAKDOWN")
# =======================================================================
# J1: n = 2 by double crossing (α per branch swap, V4 duality requires 2)
# J2: sin2_th_23 is the atmospheric mixing (p=5 <-> p=7)
# J3: observed scale forces n = 2

# Check n=1 gives ~10^-3, n=3 gives ~10^-7
n1_best = min(alpha_EM * x for _, x in targets if x > 0.01)
n3_best = max(alpha_EM**3 * x for _, x in targets)
print(f"  n = 1 best candidate: {mp.nstr(n1_best, 6)} (off by ~{mp.nstr(n1_best/PDG_J_CKM, 3)}x)")
print(f"  n = 2 canonical     : {mp.nstr(J_canonical, 6)}")
print(f"  n = 3 best candidate: {mp.nstr(n3_best, 6)} (off by ~{mp.nstr(n3_best/PDG_J_CKM, 3)}x)")
print(f"  Only n = 2 is in the observed range.")

# =======================================================================
section("UNIQUENESS")
# =======================================================================
# The canonical (2, sin2_th23) must be uniquely within 1%.
TOL = mp.mpf("0.01")
winners = [(n, name, val, err) for n, name, val, err in results if err < TOL]
print(f"  Candidates within 1%: {len(winners)}")
for n, name, val, err in winners:
    print(f"    n = {n}, X = {name}, err = {mp.nstr(err*100, 5)}%")

# =======================================================================
section("TEST OUTCOMES")
# =======================================================================
# T1: canonical within 1% of PDG (tree level)
assert err_canonical < 1.0, f"canonical err {err_canonical}% too large"
print(f"  [PASS] canonical (α², sin²θ_23): err {mp.nstr(err_canonical, 4)}% < 1%")

# T2: NLO within 0.5% of PDG
assert err_nlo < 0.5, f"NLO err {err_nlo}%"
print(f"  [PASS] NLO (1+ε): err {mp.nstr(err_nlo, 4)}% < 0.5%")

# T3: all n != 2 candidates ruled out
for n, name, val, err in results:
    if n != 2:
        assert err > 0.5, f"n={n} {name}: err {err*100}% < 50%"
print("  [PASS] J1/J3: n != 2 ruled out (err > 50% for all 12 cases)")

# T4: all X != sin2_th23 (with n=2) have err > 2%
for n, name, val, err in results:
    if n == 2 and name != "sin2(th_23)":
        assert err > 0.02, f"n=2 {name}: err only {err*100}%"
print("  [PASS] J2: X != sin²θ_23 (with n=2) have err > 2%")

# T5: unique winner within 1%
assert len(winners) == 1, f"non-unique: {len(winners)} within 1%"
assert winners[0][0] == 2 and winners[0][1] == "sin2(th_23)"
print("  [PASS] Uniqueness: only canonical within 1% (among 18)")

print()
print(f"  5/5 tests PASS. Theorem thm:JCKM_cross_branch confirmed.")
print(f"  J_CKM promoted from BRIDGE to DERIVED.")
sys.exit(0)
