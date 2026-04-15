#!/usr/bin/env python3
"""
proof_couplings.py -- Chapter 11: Coupling Constants (Weinberg Angle, Strong Coupling)

Monograph: ch11_couplings.tex
Derivation chain: s = 1/2 -> gamma_p -> sin^2(theta_W) tree -> dressing
                   s = 1/2 -> q_- = exp(-1/mu*) -> alpha_s
Zero fitted parameters.

This script derives the electroweak mixing angle and the strong coupling
constant from the sieve structure:

  Step 1. WEINBERG ANGLE (TREE LEVEL)
          sin^2(theta_W) = gamma_7^2 / sum(gamma_p^2) for p in {3,5,7}.
          The anomalous dimension gamma_7 of the heaviest active prime
          selects the weak mixing angle via the spectral ratio.
          Tree value: ~0.2380.

  Step 2. UNIVERSAL CORRECTION BASE
          The same correction base used for alpha_EM (ch10):
            C_base = C_Koide * ln(c_3D * c_2D) / (2*pi)
          Geometric rule:
            alpha_EM: 1/alpha += C_base * 26/27  (vertex, charged, order 0)
            sin^2_W:  sin^2   -= C_base * alpha  (metric, neutral, order 1)
          26/27 absent for Z (neutral current couples to ALL 27 states).
          Factor alpha: one level deeper in the hierarchy.

  Step 3. sin^2(theta_W) DRESSED
          sin^2_dressed = sin^2_tree - C_base * alpha_EM.
          Result: 0.23119 (vs PDG 0.23121, error < 0.05%).

  Step 4. STRONG COUPLING alpha_s
          alpha_s = sin^2(theta_3, q_-) / (1 - alpha_EM)
          where q_- = exp(-1/mu*) is the thermal branch.
          The strong coupling lives on the q_- branch (thermal, edge),
          while alpha_EM lives on the q_+ branch (statistical, vertex).
          Result: 0.1181 (vs PDG 0.1180, error < 0.1%).

  Step 5. CROSS-VERIFICATION
          Compare all results against pt_constants values.
          Verify the geometric rule: same C_base for both alpha and sin^2.
          Verify the q_+/q_- duality.

Theorems verified:
  D09  "Weinberg Tree Level"        (ch11_couplings.tex) -- gamma_7^2 / sum(gamma_p^2)
  D20  "sin^2(theta_W) Dressing"    (ch11_couplings.tex) -- C_base * alpha_EM correction
  D13  "Strong Coupling"            (ch11_couplings.tex) -- alpha_s from q_- branch

PT constants used:
  s = 1/2 (T1), mu* = 15 (T5), q_+ = 13/15, q_- = exp(-1/15),
  alpha_EM (ch10), C_Koide (ch10)
"""

import sys
import numpy as np
from pathlib import Path
from scipy.optimize import brentq
from scipy.integrate import quad

# Add parent to path for imports
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import (
    alpha_EM as ALPHA_EM_CONST,
    sin2_thetaW as SIN2_TW_CONST,
    sin2_thetaW_tree as SIN2_TW_TREE_CONST,
    alpha_s as ALPHA_S_CONST,
)

ck = Checker("proof_couplings", chapter="ch11", total_steps=5)

# =====================================================================
# Fundamental constants
# =====================================================================
MU_STAR = 15.0
Q_PLUS = 1.0 - 2.0 / MU_STAR     # = 13/15 (statistical branch)
Q_MINUS = np.exp(-1.0 / MU_STAR)  # thermal branch
ACTIVE_PRIMES = [3, 5, 7]
N_c = 3

# Experimental references (PDG 2024)
SIN2_PDG = 0.23121     # sin^2(theta_W) MS-bar at Z-pole
ALPHA_S_PDG = 0.1180   # alpha_s(M_Z)
ALPHA_EM_CODATA = 1.0 / 137.035999084


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
# Step 1: WEINBERG ANGLE (TREE LEVEL)
# =====================================================================
# The tree-level Weinberg angle is a spectral ratio:
#   sin^2(theta_W) = gamma_7^2 / sum_{p in {3,5,7}} gamma_p^2
#
# gamma_p(mu*) measures the RG scaling of sin^2(theta_p) at the fixed
# point. The heaviest active prime p=7 has the smallest gamma_p,
# and its fractional weight in the gamma^2 sum gives sin^2(theta_W).
ck.section("Step 1: Weinberg angle (tree level)")

# Compute anomalous dimensions at mu* = 15
gamma = {}
for p in ACTIVE_PRIMES:
    gamma[p] = gamma_p_exact(p, MU_STAR)
    qp = Q_PLUS**p
    dp = (1.0 - qp) / p
    sp = dp * (2.0 - dp)
    print(f"  p={p}: delta_{p}={dp:.8f}, sin^2={sp:.8f}, gamma_{p}={gamma[p]:.8f}")

# Tree-level Weinberg angle
sum_gamma2 = sum(g**2 for g in gamma.values())
sin2_tree = gamma[7]**2 / sum_gamma2

print(f"\n  sum(gamma_p^2) = {sum_gamma2:.8f}")
print(f"  sin^2(theta_W) tree = gamma_7^2 / sum(gamma_p^2)")
print(f"                      = {gamma[7]:.6f}^2 / {sum_gamma2:.6f}")
print(f"                      = {sin2_tree:.6f}")
print(f"  PDG value           = {SIN2_PDG}")

# Tree-level error
err_tree_pct = abs(sin2_tree - SIN2_PDG) / SIN2_PDG * 100
print(f"  Tree error          = {err_tree_pct:.4f}%")

ck.check_close("sin2_tree", sin2_tree, 0.2380, tol_pct=0.5,
               unit="sin^2_tree")
ck.check("tree_error_gt_2pct", err_tree_pct > 2.0,
         f"tree error = {err_tree_pct:.2f}% > 2% (dressing needed)")

# =====================================================================
# Step 2: UNIVERSAL CORRECTION BASE
# =====================================================================
# The correction base is shared with alpha_EM (ch10):
#   C_base = C_Koide * ln(c_3D * c_2D) / (2*pi)
#
# c_3D = ln(N_c^2) / ln(N_c^2 - 2) = ln(9)/ln(7)
# c_2D = ln(2^N_c) / ln(2^N_c - 2) = ln(8)/ln(6)
#
# Geometric rule distinguishes vertex (alpha) from metric (sin^2):
#   alpha_EM: x 26/27 (charged sector, order 0)
#   sin^2_W:  x alpha (neutral sector, order 1)
ck.section("Step 2: Universal correction base")

# Derive C_Koide
mu_end = len(ACTIVE_PRIMES) * np.pi  # = 3*pi
S_int = {}
for p in ACTIVE_PRIMES:
    val, _ = quad(lambda mu, pp=p: gamma_p_exact(pp, mu) / mu,
                  p, mu_end, limit=200)
    S_int[p] = val


def koide_Q(m1, m2, m3):
    return (m1 + m2 + m3) / (m1**0.5 + m2**0.5 + m3**0.5)**2


C_Koide = brentq(lambda C: koide_Q(np.exp(-C * S_int[3]),
                                     np.exp(-C * S_int[5]),
                                     np.exp(-C * S_int[7])) - 2.0 / 3.0,
                  5, 50)

print(f"  C_Koide = {C_Koide:.4f} (Q = 2/3 self-consistency)")

# Entropic costs from T1 (Catalan: 3^2 - 2^3 = 1)
cost_3D = np.log(N_c**2) / np.log(N_c**2 - 2)     # ln(9)/ln(7)
cost_2D = np.log(2**N_c) / np.log(2**N_c - 2)      # ln(8)/ln(6)

C_base = C_Koide * np.log(cost_3D * cost_2D) / (2.0 * np.pi)

print(f"  cost_3D = ln(9)/ln(7) = {cost_3D:.6f}")
print(f"  cost_2D = ln(8)/ln(6) = {cost_2D:.6f}")
print(f"  C_base  = C_Koide * ln(cost_3D * cost_2D) / (2*pi) = {C_base:.6f}")

ck.check_close("C_Koide", C_Koide, 18.30, tol_pct=0.1)
ck.check_close("C_base", C_base, 0.788, tol_pct=1.0, unit="C_base")

# Verify geometric rule for alpha_EM
alpha_bare = np.prod([sin2_p(p, Q_PLUS) for p in ACTIVE_PRIMES])
f_charged = (N_c**3 - 1) / N_c**3  # = 26/27
dressing_alpha = C_base * f_charged
inv_alpha_ch10 = 1.0 / alpha_bare + dressing_alpha

print(f"\n  Geometric rule verification:")
print(f"    alpha_EM: 1/alpha += C_base * 26/27 = {dressing_alpha:.6f}")
print(f"    1/alpha(tree+Delta_1) = {inv_alpha_ch10:.4f}")
print(f"    sin^2_W: sin^2 -= C_base * alpha (next step)")

ck.check_close("alpha_dressing_check", inv_alpha_ch10, 137.038, tol_pct=0.01)

# =====================================================================
# Step 3: sin^2(theta_W) DRESSED
# =====================================================================
# For the Weinberg angle, the correction uses alpha_EM instead of 26/27:
#   sin^2_dressed = sin^2_tree - C_base * alpha_EM
#
# Physical interpretation:
# - No 26/27 factor: Z boson (neutral current) couples to ALL 27 states
# - Factor alpha_EM: one level deeper in the hierarchy (metric, not vertex)
ck.section("Step 3: sin^2(theta_W) dressed (D20)")

# Use physical alpha_EM for self-consistent dressing
alpha_EM = ALPHA_EM_CODATA  # self-consistent value
delta_sin2 = C_base * alpha_EM

sin2_dressed = sin2_tree - delta_sin2

print(f"  sin^2_tree    = {sin2_tree:.6f}")
print(f"  delta_sin2    = C_base * alpha_EM")
print(f"                = {C_base:.6f} * {alpha_EM:.8f}")
print(f"                = {delta_sin2:.6f}")
print(f"  sin^2_dressed = {sin2_dressed:.6f}")
print(f"  sin^2_PDG     = {SIN2_PDG}")

err_dressed_pct = abs(sin2_dressed - SIN2_PDG) / SIN2_PDG * 100
improvement = err_tree_pct / err_dressed_pct

print(f"\n  Error tree    = {err_tree_pct:.4f}%")
print(f"  Error dressed = {err_dressed_pct:.4f}%")
print(f"  Improvement   = {improvement:.1f}x")

ck.check_close("sin2_dressed", sin2_dressed, SIN2_PDG, tol_pct=0.5,
               unit="sin^2_W")
ck.check("improvement_gt_5x", improvement > 5.0,
         f"improvement = {improvement:.1f}x > 5x")

# =====================================================================
# Step 4: STRONG COUPLING alpha_s
# =====================================================================
# alpha_s = sin^2(theta_3, q_-) / (1 - alpha_EM)
#
# where q_- = exp(-1/mu*) is the thermal branch (eigenvalue -1 of T_3).
# The q_+ branch (statistical) gives the EM coupling;
# the q_- branch (thermal) gives the strong coupling.
#
# Physical: color confinement lives on the thermal branch (edge of the
# sieve), while electromagnetism lives on the statistical branch (vertex).
ck.section("Step 4: Strong coupling alpha_s (D13)")

# q_- branch
sin2_therm = {p: sin2_p(p, Q_MINUS) for p in ACTIVE_PRIMES}

print(f"  q_- = exp(-1/mu*) = exp(-1/15) = {Q_MINUS:.10f}")
print(f"  q_+ = 1 - 2/mu*  = 13/15      = {Q_PLUS:.10f}")
print(f"\n  sin^2(theta_p, q_-) at mu* = {MU_STAR}:")
for p in ACTIVE_PRIMES:
    print(f"    p={p}: sin^2 = {sin2_therm[p]:.10f}")

# alpha_s derivation
alpha_s_derived = sin2_therm[3] / (1.0 - alpha_EM)

print(f"\n  alpha_s = sin^2(theta_3, q_-) / (1 - alpha_EM)")
print(f"         = {sin2_therm[3]:.8f} / {1.0 - alpha_EM:.8f}")
print(f"         = {alpha_s_derived:.6f}")
print(f"  PDG    = {ALPHA_S_PDG}")

err_alpha_s_pct = abs(alpha_s_derived - ALPHA_S_PDG) / ALPHA_S_PDG * 100
print(f"  Error  = {err_alpha_s_pct:.4f}%")

ck.check_close("alpha_s", alpha_s_derived, ALPHA_S_PDG, tol_pct=0.5,
               unit="alpha_s")

# Verify q_+/q_- duality
print(f"\n  q_+/q_- duality:")
print(f"    q_+ = 1 - 2/mu* = {Q_PLUS:.10f} (algebraic, statistical)")
print(f"    q_- = exp(-1/mu*)= {Q_MINUS:.10f} (exponential, thermal)")
print(f"    Ratio q_+/q_- = {Q_PLUS/Q_MINUS:.8f}")
print(f"    q_+ - q_- = {Q_PLUS - Q_MINUS:.8f}")

ck.check("q_plus_ne_q_minus", abs(Q_PLUS - Q_MINUS) > 0.01,
         f"q_+ = {Q_PLUS:.6f} != q_- = {Q_MINUS:.6f} (distinct branches)")

# =====================================================================
# Step 5: CROSS-VERIFICATION
# =====================================================================
# Compare all derived values against pt_constants.
ck.section("Step 5: Cross-verification")

# sin^2(theta_W) tree
ref_err_tree = abs(sin2_tree - SIN2_TW_TREE_CONST) / SIN2_TW_TREE_CONST * 100
print(f"  sin^2_tree (this)         = {sin2_tree:.8f}")
print(f"  sin^2_tree (pt_constants) = {SIN2_TW_TREE_CONST:.8f}")
print(f"  Agreement                 = {ref_err_tree:.6f}%")

ck.check_close("sin2_tree_vs_ref", sin2_tree, SIN2_TW_TREE_CONST,
               tol_pct=0.001)

# alpha_s
ref_err_as = abs(alpha_s_derived - ALPHA_S_CONST) / ALPHA_S_CONST * 100
print(f"\n  alpha_s (this)            = {alpha_s_derived:.8f}")
print(f"  alpha_s (pt_constants)    = {ALPHA_S_CONST:.8f}")
print(f"  Agreement                 = {ref_err_as:.6f}%")

ck.check_close("alpha_s_vs_ref", alpha_s_derived, ALPHA_S_CONST,
               tol_pct=0.01)

# Summary of all three couplings
print(f"\n  Summary of coupling constants (0 fitted parameters):")
print(f"  {'Observable':<25} {'PT derived':<15} {'Experimental':<15} {'Error'}")
print(f"  {'-'*70}")
print(f"  {'1/alpha_EM':<25} {1.0/ALPHA_EM_CONST:<15.6f} {'137.035999':<15} "
      f"{abs(1.0/ALPHA_EM_CONST - 137.035999084)/137.035999084*100:.6f}%")
print(f"  {'sin^2(theta_W)':<25} {sin2_dressed:<15.6f} {SIN2_PDG:<15} "
      f"{err_dressed_pct:.4f}%")
print(f"  {'alpha_s':<25} {alpha_s_derived:<15.6f} {ALPHA_S_PDG:<15} "
      f"{err_alpha_s_pct:.4f}%")

# =====================================================================
# BILAN
# =====================================================================
ck.summary()
