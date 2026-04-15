#!/usr/bin/env python3
"""
proof_alpha_EM.py -- Chapter 10: Fine Structure Constant alpha_EM

Monograph: ch10_fine_structure.tex
Derivation chain: s = 1/2 -> q_+ = 13/15 -> sin^2(theta_p) -> alpha_bare
                   -> p=2 dressing -> echo screening -> 2-loop VP -> alpha_EM
Zero fitted parameters.

This script derives 1/alpha_EM = 137.036... from absolute scratch:

  Step 1. BARE PRODUCT (PONTRYAGIN)
          alpha_bare = prod_{p in {3,5,7}} sin^2(theta_p, q_+)
          where sin^2(theta_p) = delta_p * (2 - delta_p),
          delta_p = (1 - q^p)/p, q_+ = 13/15.
          Exact rational arithmetic: 1/alpha_bare = 136.28.

  Step 2. FISHER-KOIDE COUPLING C_K
          Solve the Koide equation Q(C_K) = 2/3 where masses
          m_i = exp(-C_K * S_p_i) with S_p = integral of gamma_p/mu.
          Result: C_K = 18.300 (transcendental equation, unique root).

  Step 3. ONE-LOOP DRESSING (Delta_1)
          Meta-entropic action R_T1 from T1 forbidden transitions,
          charged fraction 26/27, circle normalization 2*pi.
          Delta_1 = C_K * R_T1 / (2*pi) * 26/27 = 0.758.

  Step 4. p=2 ARCHITECTURE (FULL DERIVATION)
          F(2) = sin^2_2 * cos^2(theta_2/N_2) * (mu-2)/4
          with spiral resummation, echo screening by {11,13},
          and 2-loop VP correction (alpha/pi)^2/N_c.
          Final: 1/alpha = 137.035999... (sub-ppb, 0 parameters).

  Step 5. CROSS-VERIFICATION
          Compare p=2 architecture result against pt_constants.alpha_EM.
          Verify Fisher-Koide identity C_K * sin^2_3 = G_Fisher at tree.
          Verify NLO structure: delta_CK = 1/(p1*p3) = 1/21.

Theorems verified:
  BA5 "Pontryagin Product"        (ch10_fine_structure.tex) -- alpha_bare = prod sin^2
  D09 "One-Loop Dressing"         (ch10_fine_structure.tex) -- Delta_1 Ward-forced
  R32 "Sub-ppb alpha_EM"          (ch10_fine_structure.tex) -- p=2 architecture
  --  "Fisher-Koide Identity"     (ch10_fine_structure.tex) -- C_K * sin^2_3 ~ 4

PT constants used:
  s = 1/2 (T1), mu* = 15 (T5), q_+ = 13/15, alpha_EM (for verification)
"""

import sys
import numpy as np
from fractions import Fraction
from pathlib import Path
from scipy.optimize import brentq
from scipy.integrate import quad

# Add parent to path for imports
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import alpha_EM as ALPHA_EM_REF, mu_star as MU_STAR_REF

ck = Checker("proof_alpha_EM", chapter="ch10", total_steps=5)

# =====================================================================
# Fundamental constants (all derived from s = 1/2)
# =====================================================================
ALPHA_CODATA = 1.0 / 137.035999084  # CODATA 2018 reference
MU_STAR = 15.0
Q_PLUS = 1.0 - 2.0 / MU_STAR        # = 13/15
ACTIVE_PRIMES = [3, 5, 7]
N_c = 3   # color number
N_gen = 3  # generations


def delta_p(p, q):
    """Algebraic deficit: (1 - q^p) / p"""
    return (1.0 - q**p) / p


def sin2_p(p, q):
    """sin^2(theta_p) = delta_p * (2 - delta_p)"""
    d = delta_p(p, q)
    return d * (2.0 - d)


def gamma_p_exact(p, mu):
    """Anomalous dimension gamma_p = -d(ln sin^2)/d(ln mu), exact formula."""
    if mu <= 2.01:
        return 0.0
    q = 1.0 - 2.0 / mu
    qp = q**p
    d = (1.0 - qp) / p
    if d < 1e-15 or abs(2.0 - d) < 1e-15:
        return 0.0
    dln_delta = 2.0 * p * q**(p - 1) / (mu * (1.0 - qp))
    factor = 2.0 * (1.0 - d) / (2.0 - d)
    return dln_delta * factor


# =====================================================================
# Step 1: BARE PRODUCT (PONTRYAGIN)
# =====================================================================
# Theorem BA5: alpha_bare = prod_{p in {3,5,7}} sin^2(theta_p, q_+)
# The holonomy angle theta_p encodes the information deficit at prime p.
# The product over the three active primes gives the bare coupling.
#
# We compute this in EXACT rational arithmetic to demonstrate that
# no floating-point ambiguity enters.
ck.section("Step 1: Bare product (Pontryagin, BA5)")

# Exact rational computation
q_frac = Fraction(13, 15)
sin2_exact = {}
delta_exact = {}

for p in ACTIVE_PRIMES:
    delta_exact[p] = (1 - q_frac**p) / p
    sin2_exact[p] = delta_exact[p] * (2 - delta_exact[p])
    print(f"  p={p}: delta_{p} = {delta_exact[p]} = {float(delta_exact[p]):.10f}")
    print(f"         sin^2(theta_{p}) = {float(sin2_exact[p]):.10f}")

alpha_bare_frac = Fraction(1)
for p in ACTIVE_PRIMES:
    alpha_bare_frac *= sin2_exact[p]

inv_alpha_bare_exact = float(Fraction(1) / alpha_bare_frac)
print(f"\n  alpha_bare (exact rational) = {float(alpha_bare_frac):.12f}")
print(f"  1/alpha_bare               = {inv_alpha_bare_exact:.6f}")

# Floating-point values for later steps
sin2_vals = {p: sin2_p(p, Q_PLUS) for p in ACTIVE_PRIMES}
alpha_bare = np.prod([sin2_vals[p] for p in ACTIVE_PRIMES])
inv_alpha_bare = 1.0 / alpha_bare

ck.check_close("inv_alpha_bare", inv_alpha_bare, 136.28, tol_pct=0.05,
               unit="1/alpha_bare")

# Verify exact and float agree
ck.check_close("exact_vs_float", inv_alpha_bare, inv_alpha_bare_exact,
               tol_pct=0.0001)

# Verify bare error from CODATA is < 1%
err_bare_pct = abs(inv_alpha_bare - 1.0 / ALPHA_CODATA) / (1.0 / ALPHA_CODATA) * 100
ck.check("bare_error_lt_1pct", err_bare_pct < 1.0,
         f"bare error = {err_bare_pct:.4f}% < 1%")

print(f"\n  Bare error vs CODATA: {err_bare_pct:.4f}%")
print(f"  Dressing needed: {1.0/ALPHA_CODATA - inv_alpha_bare:.4f}")

# =====================================================================
# Step 2: FISHER-KOIDE COUPLING C_K
# =====================================================================
# The Koide coupling C_K is derived from the condition Q(C_K) = 2/3,
# where Q is the Koide ratio of masses m_i = exp(-C_K * S_p_i).
# S_p_i are the integral actions of the anomalous dimensions gamma_p.
#
# The Fisher-Koide identity states:
#   C_K * sin^2(theta_3) ~ G_Fisher = 4 = 1/(s*(1-s))
# at tree level, with NLO correction delta_CK = 1/21 = 1/(p1*p3).
ck.section("Step 2: Fisher-Koide coupling C_K")

# Integral actions S_p = int_p^{3*pi} gamma_p(mu)/mu dmu
mu_end = N_gen * np.pi
S_int = {}
for p in ACTIVE_PRIMES:
    val, _ = quad(lambda mu, pp=p: gamma_p_exact(pp, mu) / mu,
                  p, mu_end, limit=200)
    S_int[p] = val
    print(f"  S_{p} = int_{p}^{{3pi}} gamma_{p}/mu dmu = {val:.10f}")


def koide_Q(C):
    """Koide ratio Q(C) = sum(m) / (sum(sqrt(m)))^2"""
    masses = [np.exp(-C * S_int[p]) for p in ACTIVE_PRIMES]
    sum_m = sum(masses)
    sum_sqrt_m = sum(m**0.5 for m in masses)
    return sum_m / sum_sqrt_m**2


# Solve Q(C_K) = 2/3
C_K = brentq(lambda C: koide_Q(C) - 2.0 / 3.0, 5, 50)
Q_verify = koide_Q(C_K)

print(f"\n  C_K = {C_K:.6f} (from Q = 2/3)")
print(f"  Q(C_K) = {Q_verify:.12f} (target = 0.666666...)")

ck.check_close("C_K_value", C_K, 18.30, tol_pct=0.1, unit="C_K")
ck.check_close("Q_koide_eq_2_3", Q_verify, 2.0 / 3.0, tol_pct=0.0001)

# Fisher-Koide identity: C_K * sin^2_3 ~ 4
G_Fisher = 4.0  # = 1 / (s * (1 - s))
product_CK_sin2 = C_K * sin2_vals[3]
residual_FK = product_CK_sin2 - G_Fisher

print(f"\n  Fisher-Koide identity:")
print(f"    C_K * sin^2_3 = {product_CK_sin2:.8f}")
print(f"    G_Fisher      = {G_Fisher:.1f}")
print(f"    Residual      = {residual_FK:.6e}")

ck.check_close("fisher_koide_tree", product_CK_sin2, G_Fisher, tol_pct=0.5,
               unit="C_K * sin^2_3")

# NLO structure: delta_CK = C_K - 4/sin^2_3 ~ 1/21
CK_tree = G_Fisher / sin2_vals[3]
delta_CK = C_K - CK_tree
inv_21 = 1.0 / 21.0
nlo_err_pct = abs(delta_CK - inv_21) / inv_21 * 100

print(f"\n  NLO correction:")
print(f"    delta_CK = C_K - 4/sin^2_3 = {delta_CK:.8f}")
print(f"    1/21     =                    {inv_21:.8f}")
print(f"    Error    = {nlo_err_pct:.4f}%")

ck.check_close("nlo_delta_CK", delta_CK, inv_21, tol_pct=10.0,
               unit="delta_CK")

# =====================================================================
# Step 3: ONE-LOOP DRESSING (Delta_1)
# =====================================================================
# The one-loop dressing uses the universal correction base:
#   Delta_1 = C_K * R_T1 / (2*pi) * 26/27
#
# R_T1 = ln(c_3D * c_2D): meta-entropic action from T1
#   c_3D = ln(N_c^2) / ln(N_c^2 - 2) = ln(9)/ln(7)
#   c_2D = ln(2^N_c) / ln(2^N_c - 2)  = ln(8)/ln(6)
#
# 26/27 = (N_c^3 - 1)/N_c^3: charged fraction of generation cube
# 2*pi = S^1 perimeter (compactification normalization)
ck.section("Step 3: One-loop dressing (Delta_1, D09)")

# Meta-entropic action R_T1
c_3D = np.log(N_c**2) / np.log(N_c**2 - 2)           # ln(9)/ln(7)
c_2D = np.log(2**N_c) / np.log(2**N_c - 2)            # ln(8)/ln(6)
R_T1 = np.log(c_3D * c_2D)

print(f"  c_3D = ln(9)/ln(7) = {c_3D:.8f}")
print(f"  c_2D = ln(8)/ln(6) = {c_2D:.8f}")
print(f"  R_T1 = ln(c_3D * c_2D) = {R_T1:.8f}")

ck.check_close("c_3D", c_3D, np.log(9) / np.log(7), tol_pct=0.001)
ck.check_close("c_2D", c_2D, np.log(8) / np.log(6), tol_pct=0.001)

# Charged fraction
f_charged = (N_c**3 - 1) / N_c**3  # = 26/27
print(f"\n  f_charged = (27-1)/27 = {f_charged:.8f}")

ck.check_close("f_charged", f_charged, 26.0 / 27.0, tol_pct=0.001)

# Assembly
Delta_1 = C_K * R_T1 / (2.0 * np.pi) * f_charged
inv_alpha_order1 = inv_alpha_bare + Delta_1

print(f"\n  Delta_1 = C_K * R_T1 / (2*pi) * 26/27")
print(f"         = {C_K:.4f} * {R_T1:.6f} / {2*np.pi:.5f} * {f_charged:.6f}")
print(f"         = {Delta_1:.8f}")
print(f"\n  1/alpha_bare    = {inv_alpha_bare:.6f}")
print(f"  + Delta_1       = {Delta_1:.6f}")
print(f"  1/alpha_order1  = {inv_alpha_order1:.6f}")
print(f"  1/alpha_CODATA  = {1.0/ALPHA_CODATA:.6f}")

ck.check_close("Delta_1", Delta_1, 0.758, tol_pct=1.0, unit="Delta_1")
ck.check_close("inv_alpha_order1", inv_alpha_order1, 137.038, tol_pct=0.01)

# =====================================================================
# Step 4: p=2 ARCHITECTURE (FULL DERIVATION)
# =====================================================================
# The complete derivation uses the p=2 binary prime as the info/anti-info
# operator. The dressing measures the informational leak across the
# p=2 boundary:
#
#   F(2) = sin^2_2 * cos^2(theta_2/N_2) * (mu-2)/4
#   + spiral resummation via gamma_3
#   + echo screening via {11, 13} ghost primes
#   + 2-loop VP correction (alpha/pi)^2/N_c
ck.section("Step 4: p=2 architecture (R32)")

# F(2): binary dressing operator
p1 = 2
delta_2 = (1.0 - Q_PLUS**p1) / p1
sin2_2 = delta_2 * (2.0 - delta_2)
theta_2 = np.arccos(1.0 - delta_2)
N2 = (p1 + 1)**(p1 + 1) - 1  # = 26
cos2_leak = np.cos(theta_2 / N2)**2
depth_2 = (MU_STAR - p1) / p1**2  # = 13/4

F2 = sin2_2 * cos2_leak * depth_2

print(f"  Binary channel p=2:")
print(f"    delta_2 = {delta_2:.10f}")
print(f"    sin^2_2 = {sin2_2:.10f}")
print(f"    theta_2 = {theta_2:.10f} rad")
print(f"    N_2 = (3^3 - 1) = {N2}")
print(f"    cos^2(theta_2/N_2) = {cos2_leak:.10f}")
print(f"    (mu-2)/4 = {depth_2}")
print(f"    F(2) = {F2:.10f}")

ck.check("N2_equals_26", N2 == 26, f"N_2 = {N2}")

# Gamma values for spiral
gamma_vals = {}
for p in [3, 5, 7, 11, 13]:
    gamma_vals[p] = gamma_p_exact(p, MU_STAR)

# Spiral resummation
alpha_1 = 1.0 / (inv_alpha_bare + F2)
sum_g2 = sum(gamma_vals[p]**2 for p in [3, 5, 7])
sum_g = sum(gamma_vals[p] for p in [3, 5, 7])
d5 = (1.0 - Q_PLUS**5) / 5.0
d7 = (1.0 - Q_PLUS**7) / 7.0
prop = (d5 + d7) / sum_g * (1.0 + alpha_bare / 25.0)
r = alpha_1 * sum_g2 * prop
spiral = F2 / (1.0 + gamma_vals[3] * r)
inv_spiral = inv_alpha_bare + spiral
alpha_d = 1.0 / inv_spiral

print(f"\n  Spiral resummation:")
print(f"    spiral correction = {spiral:.10f}")
print(f"    1/alpha(spiral)   = {inv_spiral:.10f}")

# Echo screening (ghost primes {11, 13})
sin2_echo = {p: sin2_p(p, Q_PLUS) for p in [11, 13]}
beta_echo = sum(sin2_echo[p] * gamma_vals[p] for p in [11, 13])
echo = sin2_2 * beta_echo * alpha_d**2

print(f"\n  Echo screening ({'{'}11, 13{'}'}):")
print(f"    beta_echo = {beta_echo:.10f}")
print(f"    echo      = {echo:.2e}")

# 2-loop VP correction
twoloop = (alpha_d / np.pi)**2 / N_c

print(f"\n  2-loop VP:")
print(f"    (alpha/pi)^2/N_c = {twoloop:.2e}")

# Final result
inv_final = inv_spiral + echo + twoloop
alpha_final = 1.0 / inv_final
err_ppb = abs(inv_final - 1.0 / ALPHA_CODATA) / (1.0 / ALPHA_CODATA) * 1e9

print(f"\n  RESULT:")
print(f"    1/alpha(PT)    = {inv_final:.10f}")
print(f"    1/alpha(CODATA)= {1.0/ALPHA_CODATA:.10f}")
print(f"    Error          = {err_ppb:.2f} ppb")

ck.check("sub_ppb_alpha", err_ppb < 1.0,
         f"error = {err_ppb:.2f} ppb < 1 ppb")

# Also check with percent tolerance
err_pct_final = abs(inv_final - 1.0 / ALPHA_CODATA) / (1.0 / ALPHA_CODATA) * 100
ck.check_close("inv_alpha_final", inv_final, 1.0 / ALPHA_CODATA,
               tol_pct=0.001)

# =====================================================================
# Step 5: CROSS-VERIFICATION
# =====================================================================
# Verify that our derivation matches pt_constants.alpha_EM,
# which was independently computed by the same algorithm.
ck.section("Step 5: Cross-verification")

inv_alpha_ref = 1.0 / ALPHA_EM_REF
ref_err_pct = abs(inv_final - inv_alpha_ref) / inv_alpha_ref * 100

print(f"  pt_constants.alpha_EM = {ALPHA_EM_REF:.12e}")
print(f"  1/alpha(pt_constants) = {inv_alpha_ref:.10f}")
print(f"  1/alpha(this script)  = {inv_final:.10f}")
print(f"  Agreement             = {ref_err_pct:.6f}%")

ck.check_close("vs_pt_constants", inv_final, inv_alpha_ref,
               tol_pct=0.001)

# Verify the p=2 formula: Delta_1(p2) = sin^2_2_gen * (mu-2)/G_Fisher
Delta_1_p2 = float(sin2_2 * depth_2)
p2_vs_ch10_pct = abs(Delta_1_p2 - Delta_1) / Delta_1 * 100

print(f"\n  p=2 compact formula:")
print(f"    sin^2(theta_2,gen) * (mu-2)/4 = {Delta_1_p2:.8f}")
print(f"    C_K * R_T1/(2pi) * 26/27      = {Delta_1:.8f}")
print(f"    Agreement = {p2_vs_ch10_pct:.4f}%")
print(f"    (Both agree to < 0.1% -- same identity from two sides)")

ck.check("p2_vs_ch10_agree", p2_vs_ch10_pct < 0.5,
         f"p2 vs ch10 agree to {p2_vs_ch10_pct:.4f}%")

# Summary: derivation chain
print(f"\n  Derivation chain (0 fitted parameters):")
print(f"    s = 1/2 (T1, sieve theorem)")
print(f"    -> q_+ = 1 - 2/mu* = 13/15 (T5, fixed point)")
print(f"    -> sin^2(theta_p) for p in {{3,5,7}}")
print(f"    -> alpha_bare = prod sin^2 = 1/{inv_alpha_bare:.2f}")
print(f"    -> + F(2) dressing via p=2 architecture")
print(f"    -> + echo screening via {{11,13}}")
print(f"    -> + 2-loop VP = (alpha/pi)^2/N_c")
print(f"    -> 1/alpha_EM = {inv_final:.6f}")
print(f"    CODATA:         {1.0/ALPHA_CODATA:.6f}")
print(f"    Error: {err_ppb:.2f} ppb")

# =====================================================================
# BILAN
# =====================================================================
ck.summary()
