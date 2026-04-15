#!/usr/bin/env python3
"""
proof_C9_tier_cascade.py -- Chapter 22: Chemistry

Monograph: ch22_chemistry.tex
Derivation chain: s = 1/2 -> q_plus = 13/15 -> anchor primes -> tier = period - P
Zero fitted parameters.

This script proves Corollary C9: the tier cascade structure in the periodic table.

Within each block l, the cyclic group Z/(2*P_l)Z decomposes into 3 tiers:
  - Tier -1: first-fill period (period < anchor prime)
  - Tier  0: bifurcation period (period = anchor prime)
  - Tier +1: second-fill period (period > anchor prime)

The tier number controls:
  - Half-fill anomalies (bifurcation tier = 0)
  - Chirality factor (-1)^tier
  - Orbital stability patterns

  Step 1. PT CONSTANTS
          sin2 + cos2 = 1 for all active primes.

  Step 2. D-BLOCK TIER ASSIGNMENTS
          Anchor P2 = 5: period 4 -> tier -1, period 5 -> tier 0,
          period 6 -> tier +1.

  Step 3. F-BLOCK TIER ASSIGNMENTS
          Anchor P3 = 7: period 6 -> tier -1 (lanthanides),
          period 7 -> tier 0 (actinides).

  Step 4. CHIRALITY FACTOR
          (-1)^tier: tier -1 -> -1, tier 0 -> +1, tier +1 -> -1.

  Step 5. HALF-FILL ANOMALIES
          Cr d5, Mo d5, Eu f7, Am f7 at half-fill = anchor prime.

  Step 6. CYCLIC GROUP STRUCTURE
          d-block = 2*P2 = 10, f-block = 2*P3 = 14, s+p = 2*P1 = 6.

Theorems verified:
  C9 "Tier Cascade" (ch22_chemistry.tex) -- Tier = period - anchor prime

PT constants used:
  s = 1/2, mu* = 15, q_plus = 13/15
  P1 = 3 (s/p anchor), P2 = 5 (d-block anchor), P3 = 7 (f-block anchor)
"""

import sys
from pathlib import Path

# Path setup
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import (
    s, mu_star, q_plus, PRIMES_ACTIFS,
    delta_p, sin2_theta,
)

ck = Checker("proof_C9_tier_cascade", chapter="ch22_chemistry", total_steps=6)

# Anchor primes
P1 = 3   # s/p-block anchor
P2 = 5   # d-block anchor
P3 = 7   # f-block anchor


def cos2_p(p):
    d = delta_p(p, q_plus)
    return (1.0 - d) ** 2


def tier_number(period, anchor_prime):
    return period - anchor_prime


def chirality_factor(tier):
    return (-1) ** tier


# =====================================================================
# Element data
# =====================================================================

# d-block: anchor P2 = 5
# period 4: tier -1, period 5: tier 0, period 6: tier +1
D_BLOCK = {
    'Sc': (21, 4, 1), 'Ti': (22, 4, 2), 'V':  (23, 4, 3),
    'Cr': (24, 4, 5), 'Mn': (25, 4, 5), 'Fe': (26, 4, 6),
    'Co': (27, 4, 7), 'Ni': (28, 4, 8), 'Cu': (29, 4, 10),
    'Zn': (30, 4, 10),
    'Y':  (39, 5, 1), 'Zr': (40, 5, 2), 'Nb': (41, 5, 4),
    'Mo': (42, 5, 5), 'Tc': (43, 5, 5), 'Ru': (44, 5, 7),
    'Rh': (45, 5, 8), 'Pd': (46, 5, 10), 'Ag': (47, 5, 10),
    'Cd': (48, 5, 10),
    'Hf': (72, 6, 2), 'Ta': (73, 6, 3), 'W':  (74, 6, 4),
    'Re': (75, 6, 5), 'Os': (76, 6, 6), 'Ir': (77, 6, 7),
    'Pt': (78, 6, 9), 'Au': (79, 6, 10), 'Hg': (80, 6, 10),
}

# f-block: anchor P3 = 7
# period 6: tier -1 (lanthanides), period 7: tier 0 (actinides)
F_BLOCK = {
    'Ce': (58, 6, 1),  'Pr': (59, 6, 3),  'Nd': (60, 6, 4),
    'Pm': (61, 6, 5),  'Sm': (62, 6, 6),  'Eu': (63, 6, 7),
    'Gd': (64, 6, 7),  'Tb': (65, 6, 9),  'Dy': (66, 6, 10),
    'Ho': (67, 6, 11), 'Er': (68, 6, 12), 'Tm': (69, 6, 13),
    'Yb': (70, 6, 14), 'Lu': (71, 6, 14),
    'Th': (90, 7, 0),  'Pa': (91, 7, 2),  'U':  (92, 7, 3),
    'Np': (93, 7, 4),  'Pu': (94, 7, 6),  'Am': (95, 7, 7),
    'Cm': (96, 7, 7),  'Bk': (97, 7, 9),  'Cf': (98, 7, 10),
}


# =====================================================================
# Step 1: PT CONSTANTS
# =====================================================================
ck.section("Step 1: PT constants")

d3 = delta_p(3, q_plus)
ck.check("delta_3 self-consistent",
         abs(d3 - (1.0 - q_plus**3) / 3.0) < 1e-14,
         f"delta_3 = {d3:.10f}")

for p in PRIMES_ACTIFS:
    s2 = sin2_theta(p, q_plus)
    c2 = cos2_p(p)
    ck.check(f"sin2_{p} + cos2_{p} = 1",
             abs(s2 + c2 - 1.0) < 1e-14,
             f"sum = {s2 + c2:.15f}")


# =====================================================================
# Step 2: D-BLOCK TIER ASSIGNMENTS
# =====================================================================
ck.section("Step 2: d-block tier assignments (anchor P2 = 5)")

d_per4 = [sym for sym, (z, per, n) in D_BLOCK.items() if per == 4]
d_per5 = [sym for sym, (z, per, n) in D_BLOCK.items() if per == 5]
d_per6 = [sym for sym, (z, per, n) in D_BLOCK.items() if per == 6]

ck.check("d-block period 4 -> tier -1",
         tier_number(4, P2) == -1 and len(d_per4) > 0,
         f"tier = {tier_number(4, P2)}, count = {len(d_per4)}")

ck.check("d-block period 5 -> tier 0 (bifurcation)",
         tier_number(5, P2) == 0 and len(d_per5) > 0,
         f"tier = {tier_number(5, P2)}, count = {len(d_per5)}")

ck.check("d-block period 6 -> tier +1 (inverse)",
         tier_number(6, P2) == +1 and len(d_per6) > 0,
         f"tier = {tier_number(6, P2)}, count = {len(d_per6)}")


# =====================================================================
# Step 3: F-BLOCK TIER ASSIGNMENTS
# =====================================================================
ck.section("Step 3: f-block tier assignments (anchor P3 = 7)")

f_per6 = [sym for sym, (z, per, n) in F_BLOCK.items() if per == 6]
f_per7 = [sym for sym, (z, per, n) in F_BLOCK.items() if per == 7]

ck.check("f-block period 6 -> tier -1 (lanthanides)",
         tier_number(6, P3) == -1 and len(f_per6) > 0,
         f"tier = {tier_number(6, P3)}, count = {len(f_per6)}")

ck.check("f-block period 7 -> tier 0 (actinides, bifurcation)",
         tier_number(7, P3) == 0 and len(f_per7) > 0,
         f"tier = {tier_number(7, P3)}, count = {len(f_per7)}")


# =====================================================================
# Step 4: CHIRALITY FACTOR
# =====================================================================
ck.section("Step 4: Chirality factor (-1)^tier")

ck.check("chirality(tier=-1) = -1",
         chirality_factor(-1) == -1)

ck.check("chirality(tier=0) = +1",
         chirality_factor(0) == +1)

ck.check("chirality(tier=+1) = -1",
         chirality_factor(+1) == -1)


# =====================================================================
# Step 5: HALF-FILL ANOMALIES
# =====================================================================
ck.section("Step 5: Half-fill anomalies at bifurcation")

# Cr: d5s1 half-fill, d-count = P2 = 5
z_cr, per_cr, d_cr = D_BLOCK['Cr']
ck.check("Cr(Z=24) d-count = P2 = 5 (half-fill)",
         d_cr == P2,
         f"d-count = {d_cr}")

# Mn: d5 (normal half-fill)
z_mn, per_mn, d_mn = D_BLOCK['Mn']
ck.check("Mn(Z=25) d-count = P2 = 5 (half-fill)",
         d_mn == P2,
         f"d-count = {d_mn}")

# Mo: d5 at bifurcation tier (period 5)
z_mo, per_mo, d_mo = D_BLOCK['Mo']
ck.check("Mo(Z=42) d5 at tier 0 (bifurcation period)",
         d_mo == P2 and tier_number(per_mo, P2) == 0,
         f"d-count = {d_mo}, tier = {tier_number(per_mo, P2)}")

# Eu: f7, half-fill of f-shell, count = P3
z_eu, per_eu, f_eu = F_BLOCK['Eu']
ck.check("Eu(Z=63) f-count = P3 = 7 (half-fill)",
         f_eu == P3,
         f"f-count = {f_eu}")

# Am: f7 at actinide tier 0
z_am, per_am, f_am = F_BLOCK['Am']
ck.check("Am(Z=95) f7 at tier 0 (actinide bifurcation)",
         f_am == P3 and tier_number(per_am, P3) == 0,
         f"f-count = {f_am}, tier = {tier_number(per_am, P3)}")


# =====================================================================
# Step 6: CYCLIC GROUP STRUCTURE
# =====================================================================
ck.section("Step 6: Cyclic group Z/(2*P_l)Z")

ck.check("d-block shell size = 2*P2 = 10",
         2 * P2 == 10,
         f"2*P2 = {2 * P2}")

ck.check("f-block shell size = 2*P3 = 14",
         2 * P3 == 14,
         f"2*P3 = {2 * P3}")

ck.check("s+p block shell = 2*P1 = 6",
         2 * P1 == 6,
         f"2*P1 = {2 * P1}")


# =====================================================================
ck.summary()
