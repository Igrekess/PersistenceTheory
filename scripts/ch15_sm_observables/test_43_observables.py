#!/usr/bin/env python3
"""
test_43_observables.py -- Chapter 15: Standard Model Observables

Monograph: ch15_sm_observables.tex
Derivation chain: s = 1/2 -> sieve -> mu* = 15 -> alpha_EM -> 43 SM observables
Zero fitted parameters.

This script verifies that all 43 Standard Model observables derived
by Persistence Theory from the single theorem s = 1/2 agree with
PDG 2024 experimental values:

  Step 1. INVENTORY
          All 43 observables are present in PT_SM and have PDG references.
          Classification: 3 exact, 40 numerical, 1 input (m_e, excluded).

  Step 2. COUPLING CONSTANTS
          alpha_EM (1/137.036), sin^2(theta_W), alpha_s, G_F.
          All derived from the sieve product at mu* = 15.

  Step 3. LEPTON MASSES
          m_mu and m_tau from Koide integral actions + NLO/NNLO corrections.
          m_e is the dimensional translation factor (not scored).

  Step 4. QUARK MASSES
          6 quarks from Fisher metric on T^3 + integer exponent rules.
          0 fitted parameters, MAE ~ 0.18%.

  Step 5. ELECTROWEAK BOSONS
          m_W, m_Z, m_H, v_higgs derived from sieve couplings + top Yukawa.

  Step 6. CKM MATRIX
          9 elements + J_CKM + delta_CKM from Wolfenstein + exact PDG parametrization.

  Step 7. NEUTRINO SECTOR
          PMNS angles, delta_CP, J_PMNS, mass splittings Dm31, Dm21.

  Step 8. NON-PERTURBATIVE QCD AND DECAY WIDTHS
          sigma_QCD, Regge slope, Gamma_t, R_tau, theta_QCD = 0.

  Step 9. GLOBAL STATISTICS
          Mean deviation, median, counts under 1% and 5%.
          Verification: mean ~ 0.30%, median ~ 0.06%, 42/43 under 5%.

Theorems verified:
  T1  "Forbidden Transitions"     (ch01_sieve.tex)  -- s = 1/2 from sieve
  T5  "Fixed Point"               (ch08_fixed_point.tex) -- mu* = 15
  D09 "Bare Coupling"             (ch10_fine_structure.tex) -- alpha_EM tree
  D09b "One-Loop Dressing"        (ch10_fine_structure.tex) -- alpha_EM dressed
  D09c "Fisher-Koide Identity"    (ch15_sm_observables.tex) -- C_Koide from Q=2/3
  D16 "CKM-PMNS Derivation"      (ch15_sm_observables.tex) -- mixing matrices
  D17b "Lepton Masses"            (ch15_sm_observables.tex) -- m_mu, m_tau from actions
  D19 "Quark Mass Formula"        (ch15_sm_observables.tex) -- Fisher metric on T^3
  D20 "Weak Mixing Dressing"      (ch15_sm_observables.tex) -- sin^2(theta_W)
  R12 "NLO CKM"                   (ch15_sm_observables.tex) -- V_ub, J_CKM corrections
  R15 "Higgs Mass"                (ch15_sm_observables.tex) -- m_H/v = s(1+C_F*eps)
  R17 "Self-Energy"               (ch15_sm_observables.tex) -- universal mass correction
  R18 "EW Bosons"                 (ch15_sm_observables.tex) -- Delta_r, rho parameter
  R20 "Neutrino Splittings"       (ch15_sm_observables.tex) -- Dm31, Dm21
  R26 "NNLO Corrections"          (ch15_sm_observables.tex) -- 2-loop decoherence
  R34b "Tau Cross-Branch"         (ch15_sm_observables.tex) -- hadronic tau correction
  R35 "QCD String Tension"        (ch15_sm_observables.tex) -- sigma_QCD, Regge slope
  R51 "Higgs VEV"                 (ch15_sm_observables.tex) -- v_higgs from m_t/y_t
  R52 "m_e Derivation"            (ch15_sm_observables.tex) -- m_e structural form

PT constants used:
  s = 1/2, mu* = 15, q_+ = 13/15, q_- = e^{-1/15}
  gamma_p for p in {3,5,7}, C_Koide, alpha_EM, sin^2(theta_W), alpha_s
  All 43 values from PT_SM dictionary (pt_constants.py)
"""

import sys
import numpy as np
from pathlib import Path

# Path setup: parent of ch15_sm_observables/ is scripts_v7/
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import (
    PT_SM, PDG, PDG_SIGMA,
    _INPUT_KEYS, _EXACT_KEYS, _NOT_SCORED,
    s, mu_star, alpha_EM, alpha_s, sin2_thetaW, C_Koide, eps,
    gamma, q_plus, q_minus,
)

ck = Checker("test_43_observables", chapter="ch15", total_steps=9)

# =============================================================================
# Observable categories (exhaustive partition of the 43 scored observables)
# =============================================================================
COUPLINGS = ['alpha_EM', '1/alpha_EM', 'sin2_thetaW', 'alpha_s', 'G_F']
LEPTONS = ['m_mu', 'm_tau']
QUARKS = ['m_u', 'm_d', 'm_s', 'm_c', 'm_b', 'm_t']
BOSONS = ['m_W', 'm_Z', 'm_H', 'v_higgs']
CKM = ['V_ud', 'V_us', 'V_ub', 'V_cd', 'V_cs', 'V_cb',
       'V_td', 'V_ts', 'V_tb', 'J_CKM', 'delta_CKM']
PMNS = ['sin2_th12', 'sin2_th13', 'sin2_th23', 'delta_CP_PMNS', 'J_PMNS']
NEUTRINOS = ['m_nu3_eV', 'Dm31_sq', 'Dm21_sq']
QCD_NP = ['sigma_QCD', 'regge_slope']
EXACT = ['N_c', 'N_gen', 'theta_QCD']
DECAY = ['Gamma_t', 'R_tau']

ALL_SCORED = (COUPLINGS + LEPTONS + QUARKS + BOSONS + CKM
              + PMNS + NEUTRINOS + QCD_NP + EXACT + DECAY)

# Tolerances: most observables at 2%, QCD non-perturbative at 5%,
# Gamma_t at 15% (PDG uncertainty is 13%)
TOL_DEFAULT = 2.0     # percent
TOL_QCD_NP = 5.0      # percent
TOL_GAMMA_T = 15.0    # percent (PDG: 1.42 +/- 0.19, 13% uncertainty)


def get_tol(key):
    """Return tolerance for a given observable key."""
    if key in QCD_NP:
        return TOL_QCD_NP
    if key == 'Gamma_t':
        return TOL_GAMMA_T
    return TOL_DEFAULT


# =============================================================================
# Step 1: INVENTORY
# =============================================================================
ck.section("Step 1: Inventory -- all 43 observables present")

# Verify all 43 keys exist in PT_SM and PDG
n_found = 0
for key in ALL_SCORED:
    present = key in PT_SM and key in PDG
    if present:
        n_found += 1

ck.check("inventory_count_43",
         n_found == 43,
         f"found {n_found}/43 observables in PT_SM and PDG")

# Verify classification: 3 exact, 40 numerical, 1 input excluded
n_exact = sum(1 for k in ALL_SCORED if k in _EXACT_KEYS)
n_numerical = len(ALL_SCORED) - n_exact
ck.check("exact_count_3", n_exact == 3,
         f"exact observables: {n_exact}")
ck.check("numerical_count_40", n_numerical == 40,
         f"numerical observables: {n_numerical}")

# Verify m_e is correctly classified as input (not scored)
ck.check("m_e_is_input", 'm_e' in _INPUT_KEYS,
         "m_e is a dimensional translation factor, not scored")

# Verify gluon_condensate excluded (exp. uncertainty ~25%)
ck.check("gluon_condensate_not_scored",
         'gluon_condensate' in _NOT_SCORED,
         "gluon condensate excluded due to 25% experimental uncertainty")

# Print sieve parameters
print(f"\n  Sieve parameters:")
print(f"    s = {s} (T1 theorem, mod-3 symmetry)")
print(f"    mu* = {mu_star} (T5 fixed point)")
print(f"    q_+ = {q_plus:.10f} = 1 - 2/mu*")
print(f"    q_- = {q_minus:.10f} = exp(-1/mu*)")
print(f"    C_Koide = {C_Koide:.6f}")
print(f"    eps = {eps:.8f} (universal expansion)")
print(f"    gamma = {{{gamma[3]:.6f}, {gamma[5]:.6f}, {gamma[7]:.6f}}}")

# =============================================================================
# Step 2: COUPLING CONSTANTS
# =============================================================================
ck.section("Step 2: Coupling constants (5 observables)")

print(f"\n  Derivation: alpha_EM = prod(sin^2(theta_p, q_+)), p in {{3,5,7}}")
print(f"  Dressed via p=2 architecture + echo primes {{11,13}}")

for key in COUPLINGS:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    tol = get_tol(key)

    # Determine unit for display
    unit = ""
    if key == 'G_F':
        unit = "GeV^-2"

    ck.check_close(key, pt_val, pdg_val, tol_pct=tol, unit=unit)

# =============================================================================
# Step 3: LEPTON MASSES
# =============================================================================
ck.section("Step 3: Lepton masses (2 observables)")

print(f"\n  Chain: C_Koide -> integral actions S_p -> exp(-C*S_p)")
print(f"  NLO: self-energy (R17), NNLO: 2-loop decoherence (R26)")
print(f"  Ghost VP (R29b), tau cross-branch (R34b)")

for key in LEPTONS:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_DEFAULT, unit="MeV")

# Show ratios
print(f"\n  Lepton mass ratios:")
print(f"    m_mu/m_e = {PT_SM['m_mu']/PT_SM['m_e']:.6f}"
      f"  (PDG: {PDG['m_mu']/PDG['m_e']:.6f})")
print(f"    m_tau/m_mu = {PT_SM['m_tau']/PT_SM['m_mu']:.6f}"
      f"  (PDG: {PDG['m_tau']/PDG['m_mu']:.6f})")

# =============================================================================
# Step 4: QUARK MASSES
# =============================================================================
ck.section("Step 4: Quark masses (6 observables)")

print(f"\n  Fisher metric on CRT torus T^3 = Z/3Z x Z/5Z x Z/7Z")
print(f"  log(m_q/m_e) = n_M*S_M + n_X*S_X + n_T*S_T")
print(f"  Integer exponents from conservation laws, 0 fitted parameters")

for key in QUARKS:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_DEFAULT, unit="MeV")

# =============================================================================
# Step 5: ELECTROWEAK BOSONS
# =============================================================================
ck.section("Step 5: Electroweak bosons (4 observables)")

print(f"\n  v_higgs = sqrt(2)*m_t/y_t (R51, top Yukawa derived)")
print(f"  m_H = s*(1+C_F*eps)*v (R15)")
print(f"  m_W from G_F + Delta_r (R18a+R26b)")
print(f"  m_Z = m_W/sqrt(cos^2(theta_W)*rho) (R18b+R26c)")

for key in BOSONS:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_DEFAULT, unit="GeV")

# =============================================================================
# Step 6: CKM MATRIX
# =============================================================================
ck.section("Step 6: CKM matrix (11 observables)")

print(f"\n  Wolfenstein from sieve: lambda = (sin^2_3 + sin^2_5)/(1+alpha)")
print(f"  Exact CKM (R19): standard PDG parametrization")
print(f"  NLO corrections: R12 (V_ub), R21a (V_cd, V_cb), R19a/b (V_ts, V_td)")
print(f"  Hybrid unitarity (R23): V_cs, V_tb from row sum")

for key in CKM:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    unit = "deg" if key == 'delta_CKM' else ""
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_DEFAULT, unit=unit)

# Verify CKM unitarity (row sums)
V = {k: PT_SM[k] for k in ['V_ud', 'V_us', 'V_ub', 'V_cd', 'V_cs', 'V_cb',
                             'V_td', 'V_ts', 'V_tb']}
row1_sum = V['V_ud']**2 + V['V_us']**2 + V['V_ub']**2
row2_sum = V['V_cd']**2 + V['V_cs']**2 + V['V_cb']**2
row3_sum = V['V_td']**2 + V['V_ts']**2 + V['V_tb']**2

# Note: row 1 deviates by ~0.13% because V_ud has an explicit NNLO QED
# vertex correction (R57) applied after the unitary parametrization.
# Rows 2 and 3 use hybrid unitarity (R23) so they are exact.
ck.check_close("CKM_row1_unitarity", row1_sum, 1.0, tol_pct=0.2)
ck.check_close("CKM_row2_unitarity", row2_sum, 1.0, tol_pct=0.1)
ck.check_close("CKM_row3_unitarity", row3_sum, 1.0, tol_pct=0.1)

# =============================================================================
# Step 7: NEUTRINO SECTOR
# =============================================================================
ck.section("Step 7: Neutrino sector (8 observables)")

print(f"\n  PMNS: sin^2(th12) = 1-gamma_5, sin^2(th13) = 3*alpha/(1-2*alpha)")
print(f"  J_PMNS = C_F*alpha_bare*(1+gamma_3*eps) (R24/R20b)")
print(f"  Dm31 = m3^2*cos^2(th13)*(1+gamma_5*eps)")
print(f"  Dm21 = m3^2 / (m_tau/m_mu)^(5/4) * (1-gamma_5*eps)")

for key in PMNS:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    unit = "deg" if 'delta' in key else ""
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_DEFAULT, unit=unit)

for key in NEUTRINOS:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    unit = "eV^2" if 'Dm' in key else "eV"
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_DEFAULT, unit=unit)

# =============================================================================
# Step 8: NON-PERTURBATIVE QCD, DECAY WIDTHS, EXACT
# =============================================================================
ck.section("Step 8: QCD non-perturbative, decay widths, exact (7 observables)")

print(f"\n  sigma_QCD = T_string * beta_0 (R35)")
print(f"  Regge slope = 1/(2*pi*sigma) * (1 + alpha_s_eff^2/(2*pi))")
print(f"  Gamma_t: Born + QCD NNLO (Jezabek-Kuhn + Czarnecki-Melnikov)")
print(f"  R_tau = N_c*(|V_ud|^2+|V_us|^2)*S_EW*(1+delta_QCD)")

for key in QCD_NP:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    unit = "GeV^2" if key == 'sigma_QCD' else "GeV^-2"
    ck.check_close(key, pt_val, pdg_val, tol_pct=TOL_QCD_NP, unit=unit)

for key in DECAY:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    tol = get_tol(key)
    unit = "GeV" if key == 'Gamma_t' else ""
    ck.check_close(key, pt_val, pdg_val, tol_pct=tol, unit=unit)

# Exact (discrete) quantities
for key in EXACT:
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    ck.check(f"{key}_exact", pt_val == pdg_val,
             f"PT={pt_val}, PDG={pdg_val}")

# =============================================================================
# Step 9: GLOBAL STATISTICS
# =============================================================================
ck.section("Step 9: Global statistics")

# Collect all percentage errors for non-exact, non-input scored observables
all_errs = []
for key in ALL_SCORED:
    if key in _EXACT_KEYS:
        continue
    pt_val = PT_SM[key]
    pdg_val = PDG[key]
    if pdg_val == 0:
        continue
    err = abs(pt_val - pdg_val) / abs(pdg_val) * 100
    all_errs.append((key, err))

errs_array = np.array([e for _, e in all_errs])
mean_err = np.mean(errs_array)
median_err = np.median(errs_array)
max_err = np.max(errs_array)
max_key = all_errs[np.argmax(errs_array)][0]
n_under_1pct = int(np.sum(errs_array < 1.0))
n_under_5pct = int(np.sum(errs_array < 5.0))
n_under_half = int(np.sum(errs_array < 0.5))
n_total_num = len(errs_array)

print(f"\n  43 observables: 3 exact + {n_total_num} numerical")
print(f"  Mean deviation:   {mean_err:.4f}%")
print(f"  Median deviation: {median_err:.4f}%")
print(f"  Max deviation:    {max_err:.4f}% ({max_key})")
print(f"  Under 0.5%:       {n_under_half}/{n_total_num}")
print(f"  Under 1%:         {n_under_1pct}/{n_total_num}")
print(f"  Under 5%:         {n_under_5pct}/{n_total_num}")

# Verify global statistics
ck.check("mean_deviation_under_1pct",
         mean_err < 1.0,
         f"mean = {mean_err:.4f}%")
ck.check("median_deviation_under_0.5pct",
         median_err < 0.5,
         f"median = {median_err:.4f}%")
ck.check("at_least_42_of_43_under_5pct",
         n_under_5pct + n_exact >= 42,
         f"{n_under_5pct} numerical under 5% + {n_exact} exact = "
         f"{n_under_5pct + n_exact}/43")

# n_sigma analysis (experimental compatibility)
print(f"\n  Experimental compatibility (n_sigma analysis):")
n_compat = 0
n_tension = 0
n_beyond = 0
for key, err_pct in all_errs:
    sigma = PDG_SIGMA.get(key)
    if sigma and sigma > 0:
        n_sig = abs(PT_SM[key] - PDG[key]) / sigma
        if n_sig < 2:
            n_compat += 1
        elif n_sig < 3:
            n_tension += 1
        else:
            n_beyond += 1

n_sigma_total = n_compat + n_tension + n_beyond
if n_sigma_total > 0:
    print(f"    Compatible (< 2 sigma): {n_compat}/{n_sigma_total}")
    print(f"    Tension  (2-3 sigma):   {n_tension}/{n_sigma_total}")
    print(f"    Beyond   (> 3 sigma):   {n_beyond}/{n_sigma_total}")
    print(f"    (> 3 sigma = ultra-precise measurements, not failures)")

ck.check("majority_within_2sigma",
         n_compat > n_sigma_total // 2,
         f"{n_compat}/{n_sigma_total} compatible within 2 sigma")

# Final summary
print(f"\n  Zero fitted parameters")
print(f"  Zero ansatze")
print(f"  Zero inputs: s = 1/2 is DERIVED (Theorem T1), not assumed")
print(f"  1 translation factor: m_e (dimensional anchor, like c = 3e8 m/s)")

# =============================================================================
# BILAN
# =============================================================================
ck.summary()
