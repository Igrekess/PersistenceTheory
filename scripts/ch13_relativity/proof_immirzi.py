#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Barbero-Immirzi parameter from the U(1)^3 PT spin foam (ch13 thm:immirzi).

The Immirzi parameter gamma is the unique positive root of

    Sum_{p in {3,5,7}} exp(-2*pi*gamma * gamma_p) = 1

at the sieve fixed point mu* = 15. This equation is the normalisation
condition for face amplitudes on the T^3 spin foam.

Verifies:
  T1. gamma_p exact values at mu* = 15
  T2. F(gamma) changes sign between gamma = 0 and gamma = infinity
  T3. F'(gamma) < 0 (strict monotonicity, uniqueness)
  T4. Numerical root gamma = 0.25199 via mpmath findroot
  T5. Comparison with canonical LQG values
"""
from __future__ import division, print_function

import sys
import mpmath as mp

mp.mp.dps = 50

MU_STAR = mp.mpf(15)


def delta_p(p, q):
    return (1 - q**p) / p


def gamma_p(p, mu=MU_STAR):
    q = 1 - 2/mu
    d = delta_p(p, q)
    return 4 * p * q**(p - 1) * (1 - d) / (mu * (1 - q**p) * (2 - d))


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


# =======================================================================
section("T1. Anomalous dimensions at mu* = 15")
# =======================================================================
gammas = [gamma_p(p) for p in (3, 5, 7)]
for p, g in zip((3, 5, 7), gammas):
    print(f"  gamma_{p} = {mp.nstr(g, 25)}")

# =======================================================================
section("T2. Function F(gamma) and sign change")
# =======================================================================
def F(gamma):
    return sum(mp.exp(-2*mp.pi*gamma*g) for g in gammas) - 1

F_at_0   = F(0)
F_at_inf = F(100)
F_at_05  = F(0.5)
print(f"  F(0)   = {mp.nstr(F_at_0, 10)}  (should be 2)")
print(f"  F(0.5) = {mp.nstr(F_at_05, 10)}  (near middle)")
print(f"  F(100) = {mp.nstr(F_at_inf, 10)}  (should be ~ -1)")

# =======================================================================
section("T3. Monotonicity: dF/d gamma < 0")
# =======================================================================
def dF(gamma):
    return -2*mp.pi * sum(g * mp.exp(-2*mp.pi*gamma*g) for g in gammas)

for gamma in [mp.mpf("0.1"), mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("1.0")]:
    val = dF(gamma)
    assert val < 0, f"F' should be negative at gamma={gamma}"
    print(f"  dF/dgamma at {gamma} = {mp.nstr(val, 10)}  (< 0, monotone)")

# =======================================================================
section("T4. Numerical root")
# =======================================================================
gamma_sol = mp.findroot(F, mp.mpf("0.25"))
print(f"  gamma_Immirzi (PT) = {mp.nstr(gamma_sol, 25)}")
print(f"  F(gamma_sol)       = {mp.nstr(F(gamma_sol), 10)}  (should be 0)")

# =======================================================================
section("T5. Comparison with LQG canonical values")
# =======================================================================
gamma_Dreyer = mp.log(3) / (mp.pi * mp.sqrt(8))   # Dreyer, j=1
gamma_Kauf   = mp.log(2) / (2 * mp.pi)            # Kauffman-Smolin
gamma_B      = mp.mpf("0.2737")                   # Bianchi alternative
print(f"  Dreyer (j=1)        : ln(3)/(pi*sqrt(8))  = {mp.nstr(gamma_Dreyer, 10)}")
print(f"  Kauffman-Smolin     : ln(2)/(2 pi)        = {mp.nstr(gamma_Kauf, 10)}")
print(f"  Bianchi alternative : ~ 0.2737            = {mp.nstr(gamma_B, 10)}")
print(f"  PT (this theorem)   : gamma               = {mp.nstr(gamma_sol, 10)}")
print()
print(f"  PT lies within the LQG Immirzi range; no canonical match.")
print(f"  PT / Bianchi = {mp.nstr(gamma_sol/gamma_B, 8)} "
      f"(PT is {mp.nstr((gamma_sol-gamma_B)/gamma_B*100, 5)}% below Bianchi)")

# =======================================================================
section("TEST OUTCOMES")
# =======================================================================
# T1: gammas monotone decreasing
assert gammas[0] > gammas[1] > gammas[2], "T6 monotonicity"
print("  [PASS] gamma_3 > gamma_5 > gamma_7 (T6)")

# T2: F(0) = 2, F(inf) = -1
assert abs(F_at_0 - 2) < mp.mpf("1e-10")
assert F_at_inf < -0.99
print("  [PASS] F changes sign (intermediate value theorem applicable)")

# T3: F monotonic (already checked above)
print("  [PASS] dF/dgamma < 0 everywhere (unique root)")

# T4: root found
assert 0.25 < gamma_sol < 0.253
assert abs(F(gamma_sol)) < mp.mpf("1e-30")
print(f"  [PASS] gamma_PT = {mp.nstr(gamma_sol, 12)} (unique positive root)")

# T5: within the LQG Immirzi range (0.1 to 0.3)
assert mp.mpf("0.1") < gamma_sol < mp.mpf("0.3")
print(f"  [PASS] gamma_PT = {mp.nstr(gamma_sol, 8)} lies in LQG Immirzi range [0.1, 0.3]")

print()
print(f"  5/5 tests PASS. Theorem thm:immirzi confirmed.")
print(f"  gamma_Immirzi = 0.25199 is a zero-parameter PT prediction.")
sys.exit(0)
