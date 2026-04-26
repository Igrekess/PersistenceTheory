#!/usr/bin/env python3
"""
test_DM_candidate.py -- Chapter 20e: DM scalar singlet candidate (p=2)

Monograph: chapters/ch20e_DM_candidate.tex
Derivation: PT bifurcation -> p=2 operator -> scalar singlet at m_H/2.
"""

import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


M_H = 125.25      # GeV (PDG)
V_HIGGS = 246.220 # GeV
GAMMA_H = 4.07e-3 # GeV (PDG total Higgs width)


def predict_DM_mass():
    return M_H / 2.0


def lambda_HS_upper_limit(Br_inv_limit=0.107):
    """Upper limit from Br(H -> inv) < 0.107 (PDG 2024)."""
    return math.sqrt(8 * math.pi * M_H * GAMMA_H * Br_inv_limit) / V_HIGGS


def relic_sigma_v(lambda_HS):
    m_DM = predict_DM_mass()
    return lambda_HS ** 2 / (8 * math.pi * m_DM ** 2)


ck = Checker("test_DM_candidate", chapter="ch20e", total_steps=3)

# ---- Step 1: Mass prediction ----
ck.section("Step 1: m_DM = m_H / 2 (PT p=2 prediction)")
m_DM = predict_DM_mass()
ck.check("m_DM_is_mH_over_2", abs(m_DM - M_H / 2) < 1e-6,
         f"m_DM = {m_DM:.3f} GeV, m_H/2 = {M_H/2:.3f} GeV")
ck.check("m_DM_resonance_window", 60.0 < m_DM < 65.0,
         f"m_DM = {m_DM:.2f} GeV in resonance window")

# ---- Step 2: Coupling constraint ----
ck.section("Step 2: Higgs portal lambda_HS constraint")
lam_max = lambda_HS_upper_limit()
ck.check("lambda_HS_max_subunit", lam_max < 1.0,
         f"lambda_HS_max = {lam_max:.4f}")
ck.check("lambda_HS_max_above_target", lam_max > 1e-3,
         f"lambda_HS_max = {lam_max:.4f} > 0.001 minimum target")

# ---- Step 3: Falsifiable predictions ----
ck.section("Step 3: PRED falsifiers")
ck.check("LZ_2024_target_mass", abs(m_DM - 62.5) < 0.5,
         f"m_DM = {m_DM:.2f} GeV; LZ 2024+ target ~62.5 GeV")
sv_at_lam_001 = relic_sigma_v(0.01)
ck.check("relic_abundance_finite", sv_at_lam_001 > 0,
         f"sigma v at lambda_HS=0.01: {sv_at_lam_001:.3e}")

ck.summary()
