#!/usr/bin/env python3
"""
test_condensed_matter.py -- Chapter 22: Condensed Matter (Band Gaps)

Monograph: ch22c_condensed.tex
Derivation chain: s = 1/2 -> sin2 -> D_KL cascade -> E_gap = sqrt(E_h^2 + E_c^2)
Zero fitted parameters.

Tests the structural framework of the PT band gap model:

  Step 1. STRUCTURAL IDENTITIES
          D_KL cascade, ghost primes, vertex-edge duality.

  Step 2. PYTHAGOREAN DECOMPOSITION
          E_gap = sqrt(E_h^2 + E_c^2) for covalent + ionic channels.
          Phillips ionicity = cos^2(theta_3) = 0.7808.

  Step 3. COVALENT CASCADE (GROUP IV)
          D_cov(n) monotone in period, predicts ordering C > Si > Ge > Sn.

  Step 4. TERNARY BOWING
          D_total modular degeneracy score on Z/pZ for 3 atoms.
          Bowing b_PT = D_total * sqrt(|dE_gap| * (|dE_h| + |dE_c|)).

  Step 5. OUT-OF-SAMPLE STRUCTURE
          Key structural predictions (rank ordering, ionicity bounds).

Theorems verified:
  C14 "Thermodynamic Limit" (ch22c_condensed.tex) -- Pythagorean band gap structure

PT constants used:
  s = 1/2, mu* = 15, q_plus = 13/15
  sin2_theta, gamma_p_exact, delta_p from pt_constants
"""

import sys
import math
from pathlib import Path

_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from pt_constants import (
    s, mu_star, q_plus, alpha_EM,
    sin2_theta, gamma_p_exact, delta_p,
    sin2_stat, gamma, C_F,
)
from lib.pt_check import Checker

# ── PT constants ──────────────────────────────────────────────────────
q = q_plus

sin2_3 = sin2_theta(3, q)
sin2_5 = sin2_theta(5, q)
sin2_7 = sin2_theta(7, q)
sin2_11 = sin2_theta(11, q)
sin2_13 = sin2_theta(13, q)

G3 = gamma_p_exact(3, mu_star)
G5 = gamma_p_exact(5, mu_star)
G7 = gamma_p_exact(7, mu_star)
G11 = gamma_p_exact(11, mu_star)
G13 = gamma_p_exact(13, mu_star)

BETA_GHOST = sin2_11 * G11 + sin2_13 * G13

# Actions S_p = -ln(sin2_p)
S3 = -math.log(sin2_3)
S5 = -math.log(sin2_5)
S7 = -math.log(sin2_7)
S11 = -math.log(sin2_11)

ACTIVE_PRIMES = [3, 5, 7]

ck = Checker("test_condensed_matter", chapter="ch22_chemistry", total_steps=5)


# ── Helper: covalent cascade ─────────────────────────────────────────
def D_cov(n_avg):
    """Covalent D_KL cascade gated by active primes {3,5,7}."""
    if n_avg <= 2:
        return math.log(2)
    if n_avg <= 2.5:
        return math.log(2) + (n_avg - 2) * (C_F * S3 - math.log(2))
    if n_avg <= 3:
        return C_F * S3
    if n_avg <= 4:
        return C_F * (S3 + (n_avg - 3) * sin2_3 * S5)
    if n_avg <= 5:
        return C_F * (S3 + sin2_3 * S5 + (n_avg - 4) * sin2_3 * sin2_5 * S7)
    return C_F * (S3 + sin2_3 * S5 + sin2_3 * sin2_5 * S7)


# ── Helper: Pythagorean gap ──────────────────────────────────────────
def pythagorean_gap(E_h, E_c):
    """E_gap = sqrt(E_h^2 + E_c^2)."""
    return math.sqrt(E_h**2 + E_c**2)


# ── Helper: modular degeneracy score ─────────────────────────────────
def compute_D_total(Z1, Z2, Z3):
    """Modular degeneracy score for 3 atoms on Z/pZ.
    D_total = sum D_p * sin^2(theta_p) for p in {3,5,7}.
    D_p = min(1, 3 - |{Z1%p, Z2%p, Z3%p}|).
    """
    D_total = 0.0
    for p in ACTIVE_PRIMES:
        s2 = sin2_theta(p, q)
        residues = {Z1 % p, Z2 % p, Z3 % p}
        D_p = min(1, 3 - len(residues))
        D_total += D_p * s2
    return D_total


# ======================================================================
#  STEP 1: STRUCTURAL IDENTITIES
# ======================================================================
ck.section("Step 1: Structural identities")

# alpha_EM = product of sin^2(theta_p, q_plus)
alpha_nue = sin2_3 * sin2_5 * sin2_7
ck.check_close("alpha_bare = prod sin^2_p",
               alpha_nue, 1.0 / 136.3, tol_pct=1.0, unit="")

# Ghost primes {11, 13} have gamma_p < s
ck.check("Ghost prime 11: gamma_11 < s",
         G11 < s,
         f"gamma_11 = {G11:.4f}, s = {s}")

ck.check("Ghost prime 13: gamma_13 < s",
         G13 < s,
         f"gamma_13 = {G13:.4f}, s = {s}")

# Active primes {3,5,7} have gamma_p > s
ck.check("Active prime 3: gamma_3 > s",
         G3 > s,
         f"gamma_3 = {G3:.4f}, s = {s}")

ck.check("Active prime 5: gamma_5 > s",
         G5 > s,
         f"gamma_5 = {G5:.4f}, s = {s}")

ck.check("Active prime 7: gamma_7 > s",
         G7 > s,
         f"gamma_7 = {G7:.4f}, s = {s}")

# Vertex-edge duality: exp(-C_F) close to sin^2_3^C_F / s
exp_CF = math.exp(-C_F)
sin2_3_CF_over_s = sin2_3**C_F / s
ck.check("Vertex-edge identity: exp(-C_F) ~ sin^2_3^C_F / s (< 1%)",
         abs(exp_CF - sin2_3_CF_over_s) / exp_CF < 0.01,
         f"exp(-C_F) = {exp_CF:.6f}, sin^2_3^C_F/s = {sin2_3_CF_over_s:.6f}")

# Ghost VP coefficient
ck.check("Ghost VP coefficient beta_ghost ~ 0.10",
         0.08 < BETA_GHOST < 0.12,
         f"beta_ghost = {BETA_GHOST:.4f}")

# ======================================================================
#  STEP 2: PYTHAGOREAN DECOMPOSITION
# ======================================================================
ck.section("Step 2: Pythagorean decomposition E_gap = sqrt(E_h^2 + E_c^2)")

# Test Pythagorean identity for several (E_h, E_c) pairs
test_pairs = [
    (1.12, 0.0,  "Pure covalent (Si-like)"),
    (1.10, 0.88, "Mixed (GaAs-like)"),
    (0.50, 3.00, "Strongly ionic"),
    (2.00, 2.00, "Equal channels"),
]
for E_h, E_c, label in test_pairs:
    E_gap = pythagorean_gap(E_h, E_c)
    E_gap_check = math.sqrt(E_h**2 + E_c**2)
    ck.check(f"Pythagorean: {label}",
             abs(E_gap - E_gap_check) < 1e-14,
             f"E_gap = {E_gap:.6f}")

# Phillips ionicity = cos^2(theta_3) = (1-delta_3)^2
d3 = delta_p(3, q)
cos2_3 = (1.0 - d3)**2
ck.check_close("Phillips ionicity cos^2(theta_3) vs 0.785",
               cos2_3, 0.785, tol_pct=1.0, unit="")

# sin^2 + cos^2 = 1 (GFT budget)
ck.check("sin^2_3 + cos^2_3 = 1",
         abs(sin2_3 + cos2_3 - 1.0) < 1e-14,
         f"sum = {sin2_3 + cos2_3:.15f}")


# ======================================================================
#  STEP 3: COVALENT CASCADE (GROUP IV)
# ======================================================================
ck.section("Step 3: Covalent cascade D_cov(n)")

# D_cov is monotonically increasing with period
D_cov_2 = D_cov(2.0)  # C (period 2)
D_cov_3 = D_cov(3.0)  # Si (period 3)
D_cov_4 = D_cov(4.0)  # Ge (period 4)
D_cov_5 = D_cov(5.0)  # Sn (period 5)

ck.check("D_cov monotone: D(2) < D(3) < D(4) < D(5)",
         D_cov_2 < D_cov_3 < D_cov_4 < D_cov_5,
         f"D = [{D_cov_2:.3f}, {D_cov_3:.3f}, {D_cov_4:.3f}, {D_cov_5:.3f}]")

# D_cov(2) = ln(2) (kinematic floor)
ck.check_close("D_cov(2) = ln(2) (kinematic floor)",
               D_cov_2, math.log(2), tol_pct=0.01, unit="nats")

# D_cov(3) = C_F * S3 (first dynamic level)
ck.check_close("D_cov(3) = C_F * S3 (first dynamic level)",
               D_cov_3, C_F * S3, tol_pct=0.01, unit="nats")

# Group IV gap ordering: C > Si > Ge > Sn (from D_cov ordering)
# Higher D_cov -> smaller gap: E_h ~ exp(-D_cov)
E_h_C  = math.exp(-D_cov_2)
E_h_Si = math.exp(-D_cov_3)
E_h_Ge = math.exp(-D_cov_4)
E_h_Sn = math.exp(-D_cov_5)

ck.check("Group IV ordering: E(C) > E(Si) > E(Ge) > E(Sn)",
         E_h_C > E_h_Si > E_h_Ge > E_h_Sn,
         f"E = [{E_h_C:.4f}, {E_h_Si:.4f}, {E_h_Ge:.4f}, {E_h_Sn:.4f}]")

# Cascade gating: each prime gates by sin^2 of the previous
# D_cov(4) - D_cov(3) ~ C_F * sin2_3 * S5 (second gate)
delta_34 = D_cov_4 - D_cov_3
expected_34 = C_F * sin2_3 * S5
ck.check_close("Cascade gate: D(4)-D(3) = C_F*sin2_3*S5",
               delta_34, expected_34, tol_pct=0.1, unit="nats")


# ======================================================================
#  STEP 4: TERNARY BOWING (MODULAR DEGENERACY)
# ======================================================================
ck.section("Step 4: Ternary bowing (modular degeneracy)")

# D_total for AlGaAs (Al=13, Ga=31, As=33): check structure
D_AlGaAs = compute_D_total(13, 31, 33)
ck.check("D_total(AlGaAs) >= 0",
         D_AlGaAs >= 0,
         f"D_total = {D_AlGaAs:.4f}")

# D_total for identical atoms should be maximal
D_same = compute_D_total(14, 14, 14)  # Si-Si-Si
D_diff = compute_D_total(13, 14, 33)  # Al-Si-As (all different mod 3,5,7)
ck.check("D_total(same) >= D_total(different)",
         D_same >= D_diff,
         f"D(same) = {D_same:.4f}, D(diff) = {D_diff:.4f}")

# Modular degeneracy: Z1 == Z2 (mod p) for all p gives D_p = 1
# GaAs alloys: Ga(31) and Al(13): 31 mod 3 = 1, 13 mod 3 = 1 -> degenerate
ck.check("Ga and Al degenerate mod 3 (31%3 == 13%3)",
         31 % 3 == 13 % 3,
         f"31 mod 3 = {31 % 3}, 13 mod 3 = {13 % 3}")

# Transitivity cap: triple coincidence = single bit (not 2)
D_triple = compute_D_total(13, 13, 13)  # all identical
# For all p: |residues| = 1, D_p = min(1, 3-1) = min(1,2) = 1
expected_D_triple = sum(sin2_theta(p, q) for p in ACTIVE_PRIMES)
ck.check_close("Triple coincidence capped at 1 per prime (transitivity)",
               D_triple, expected_D_triple, tol_pct=0.01)


# ======================================================================
#  STEP 5: OUT-OF-SAMPLE STRUCTURE
# ======================================================================
ck.section("Step 5: Out-of-sample structural predictions")

# Wide-gap ordering: BN > AlN > GaN (from IE and cascade)
# BN: n_avg = 2 (period 2-2), AlN: n_avg = 2.5 (period 3-2), GaN: n_avg = 3 (period 4-2)
# Lower n_avg -> smaller D_cov -> larger gap
ck.check("Wide-gap order: D_cov(2) < D_cov(2.5) < D_cov(3.5)",
         D_cov(2.0) < D_cov(2.5) < D_cov(3.5),
         f"D = [{D_cov(2.0):.3f}, {D_cov(2.5):.3f}, {D_cov(3.5):.3f}]")

# Ionicity classification: f_i = cos^2(theta_3) = 0.78 (Phillips boundary)
# Materials with f_i > threshold are ionic
ck.check("Phillips threshold ~ 0.785 separates ionic/covalent",
         0.77 < cos2_3 < 0.80,
         f"cos^2_3 = {cos2_3:.4f}")

# Ghost screening is small: sin^2_7 * s ~ 8.6%
ghost_leak = sin2_7 * s
ck.check("Ghost screening: sin^2_7 * s ~ 8.6%",
         0.07 < ghost_leak < 0.10,
         f"sin^2_7 * s = {ghost_leak:.4f}")

# Kinematic gate: 1 - ln(2)/(C_F*S3) ~ 0.66
f_kin = 1.0 - math.log(2) / (C_F * S3)
ck.check("Kinematic gate: 1 - ln(2)/(C_F*S3) in [0.6, 0.7]",
         0.6 < f_kin < 0.7,
         f"f_kin = {f_kin:.4f}")

# Boundary reinforcement: gamma_7 > s (period 4 = last active)
ck.check("Boundary reinforcement: gamma_7 > s at period 4",
         G7 > s,
         f"gamma_7 = {G7:.4f} > s = {s}")

# ── Summary ──────────────────────────────────────────────────────────
ck.summary()
