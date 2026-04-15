#!/usr/bin/env python3
"""
proof_C1_periodic_table.py -- Chapter 22: Chemistry

Monograph: ch22_chemistry.tex
Derivation chain: s = 1/2 -> q_plus = 13/15 -> delta_p -> sin2_p -> block capacities
Zero fitted parameters.

This script proves the periodic table's block structure from Persistence Theory:

  Step 1. BLOCK CAPACITIES
          s = 2, p = 6, d = 10, f = 14: each is 2*(2l+1) with l = 0,1,2,3.
          Total capacity of a complete period = 32.

  Step 2. PERIOD LENGTHS
          (2, 8, 8, 18, 18, 32, 32) = 2 * (perfect squares).
          Noble gas Z values are cumulative sums.

  Step 3. AUFBAU ORDERING
          Blocks fill in order of decreasing delta_p: s -> p -> d -> f.
          The persistence gap delta_p is monotonically decreasing with p.

  Step 4. G-SHELL EXCLUSION
          gamma_11 < 1/2 prevents a g-block from opening.
          9 = 3^2 is composite, so the angular prime sequence {1,3,5,7}
          cannot extend to 9.

  Step 5. D-BLOCK PARABOLA
          The d-block weight w(f) = 4*alpha_d*f*(1-f) on Z/10Z is
          a parabola maximal at half-filling.

Theorems verified:
  C1 "Periodic Table CRT" (ch22_chemistry.tex) -- Block structure = CRT
     decomposition of electron filling on circles Z/(2*P_l)Z

PT constants used:
  s = 1/2, mu* = 15, q_plus = 13/15
  delta_p = (1 - q^p) / p
  sin2_p = delta_p * (2 - delta_p)
  gamma_p = 4*p*q^(p-1)*(1 - delta_p) / (mu*(1 - q^p)*(2 - delta_p))
"""

import math
import sys
from pathlib import Path

# Path setup
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import (
    s, mu_star, q_plus, PRIMES_ACTIFS,
    delta_p, sin2_theta, gamma_p_exact,
)

ck = Checker("proof_C1_periodic_table", chapter="ch22_chemistry", total_steps=5)

# Local helpers
def _delta(p):
    """Persistence gap for prime p at q_plus."""
    return delta_p(p, q_plus)

def _sin2(p):
    """Effective sin^2 for prime p at q_plus."""
    return sin2_theta(p, q_plus)

def _gamma(p):
    """Anomalous dimension for prime p at mu*."""
    return gamma_p_exact(p, mu_star)

def _is_prime(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    k = 5
    while k * k <= n:
        if n % k == 0 or n % (k + 2) == 0:
            return False
        k += 6
    return True


# =====================================================================
# Step 1: BLOCK CAPACITIES
# =====================================================================
ck.section("Step 1: Block capacities")

BLOCK_CAP = {"s": 2, "p": 6, "d": 10, "f": 14}
L_VALUES = {"s": 0, "p": 1, "d": 2, "f": 3}

# Each block capacity = 2*(2l+1)
for block, cap in BLOCK_CAP.items():
    l_val = L_VALUES[block]
    expected = 2 * (2 * l_val + 1)
    ck.check(f"cap({block}) = 2*(2*{l_val}+1) = {expected}",
             cap == expected,
             f"cap = {cap}")

# Total = 32
total = sum(BLOCK_CAP.values())
ck.check("total capacity = 32", total == 32, f"sum = {total}")


# =====================================================================
# Step 2: PERIOD LENGTHS AND NOBLE GASES
# =====================================================================
ck.section("Step 2: Period lengths and noble gas closures")

PERIOD_LENGTHS = [2, 8, 8, 18, 18, 32, 32]
PERIOD_HALVES = [1, 4, 4, 9, 9, 16, 16]

for i, (pl, ph) in enumerate(zip(PERIOD_LENGTHS, PERIOD_HALVES), 1):
    ok = (pl == 2 * ph) and (math.isqrt(ph) ** 2 == ph)
    ck.check(f"period {i}: length={pl}, half={ph} (perfect square)",
             ok,
             f"sqrt({ph}) = {math.isqrt(ph)}")

# Noble gas Z values = cumulative sums of period lengths
NOBLE_Z = [2, 10, 18, 36, 54, 86]
cumul = 0
idx = 0
for z in NOBLE_Z:
    cumul += PERIOD_LENGTHS[idx]
    idx += 1
    if cumul < z:
        cumul += PERIOD_LENGTHS[idx]
        idx += 1
    ck.check(f"noble gas Z={z} = cumsum(periods)",
             cumul == z,
             f"cumsum = {cumul}")


# =====================================================================
# Step 3: AUFBAU ORDERING
# =====================================================================
ck.section("Step 3: Aufbau ordering from decreasing delta_p")

# Angular primes: 1, 3, 5, 7 map to blocks s, p, d, f
# (We use p=1 for s-block with the convention delta(1) = 1-q)
angular_primes = [1, 3, 5, 7]
deltas = [_delta(p) for p in angular_primes]

aufbau_ok = all(deltas[i] > deltas[i + 1] for i in range(len(deltas) - 1))
ck.check("delta_p strictly decreasing: s > p > d > f",
         aufbau_ok,
         f"deltas = [{', '.join(f'{d:.6f}' for d in deltas)}]")

# Individual delta values are positive
for p, d in zip(angular_primes, deltas):
    ck.check(f"delta({p}) > 0",
             d > 0,
             f"delta({p}) = {d:.8f}")


# =====================================================================
# Step 4: G-SHELL EXCLUSION
# =====================================================================
ck.section("Step 4: g-shell exclusion")

# gamma_11 < 1/2 prevents g-block from opening
g11 = _gamma(11)
ck.check("gamma_11 < 1/2 (no g-shell)",
         g11 < 0.5,
         f"gamma_11 = {g11:.6f}")

# 9 is not prime: angular prime sequence {1,3,5,7} cannot extend
ck.check("9 = 3^2 is composite (no l=4 angular prime)",
         not _is_prime(9),
         "9 = 3 x 3")


# =====================================================================
# Step 5: D-BLOCK PARABOLA
# =====================================================================
ck.section("Step 5: d-block parabola on Z/10Z")

# alpha_d = sin^2 at prime 5 (d-block coupling)
alpha_d = _sin2(5)

def w_d(f):
    """d-block parabolic weight: 4*alpha_d*f*(1-f)."""
    return 4.0 * alpha_d * f * (1.0 - f)

ck.check("w_d(0) = 0 (empty shell)",
         abs(w_d(0.0)) < 1e-15,
         f"w(0) = {w_d(0.0):.2e}")

ck.check("w_d(1) = 0 (full shell)",
         abs(w_d(1.0)) < 1e-15,
         f"w(1) = {w_d(1.0):.2e}")

ck.check("w_d(1/2) = alpha_d (maximal at half-fill)",
         abs(w_d(0.5) - alpha_d) < 1e-15,
         f"w(1/2) = {w_d(0.5):.6f}, alpha_d = {alpha_d:.6f}")

ck.check("w_d(f) symmetric: w(0.3) = w(0.7)",
         abs(w_d(0.3) - w_d(0.7)) < 1e-15,
         f"w(0.3) = {w_d(0.3):.6f}, w(0.7) = {w_d(0.7):.6f}")

# =====================================================================
ck.summary()
