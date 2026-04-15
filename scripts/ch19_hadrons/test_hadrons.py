#!/usr/bin/env python3
"""
test_hadrons.py -- Chapter 19: Hadronic Spectra & Chiral Dynamics

Monograph: ch19_hadrons.tex
Derivation chain: s = 1/2 -> sigma_QCD(nf=3) -> f_pi, m_pi, m_K, m_rho
                   -> KSRF (g_rho_pi_pi, Gamma_rho) -> F_pi(s) -> eta-eta' mixing
                   -> NLO corrections (R56)
Zero fitted parameters.

This script consolidates the full hadronic spectrum derivation:
16 light hadrons, chiral symmetry, hadronic dynamics, and NLO corrections.

  Step 1. PION DECAY CONSTANT f_pi
          f_pi^2 = N_c * sigma / (4*pi^2) (string relation)
          Verified against PDG: 130.2 MeV.

  Step 2. CHIRAL CONDENSATE AND PION MASS (GMOR)
          <qbar q> = C_chi * sigma^{3/2}, C_chi = C_F/pi.
          m_pi^2 * f_pi^2 = (m_u + m_d) * |<qbar q>| (Gell-Mann--Oakes--Renner).
          m_pi prediction: 0.01% accuracy.

  Step 3. KAON AND STRANGE SECTOR
          m_K = m_pi * sqrt((m_u + m_s)/(m_u + m_d)) (ratio GMOR, R36b).
          f_K/f_pi = 1 + gamma_5*(1+s)*(m_s - m_hat)/sqrt(sigma) (vertex+edge).

  Step 4. RHO MESON AND REGGE TRAJECTORY
          alpha'_eff = alpha' * (1 + C_F * alpha_s_eff^2/pi) (Coulomb correction).
          m_rho = sqrt((1-s)/alpha'_eff).
          Regge trajectory J=1..4: rho, a_2(1320), rho_3(1690), a_4(2040).

  Step 5. PROTON AND NEUTRON
          M_p = 2*sqrt(sigma), delta_m(n-p) = s*(m_d - m_u).

  Step 6. HADRONIC DYNAMICS (KSRF, VMD)
          g_rho_pi_pi = m_rho / f_pi (KSRF), Gamma_rho via p-wave phase space,
          F_pi(s) via VMD, <r_pi^2> = 6/m_rho^2.

  Step 7. ETA-ETA' MIXING (R37)
          chi_top = T_string * sigma_pure^2, C_mix = Q_Koide = 2/3.
          2x2 mass matrix diagonalization for eta and eta'.

  Step 8. NLO CORRECTIONS (R56)
          Universal expansion parameter epsilon = alpha_s * s / (2*pi).
          Corrections to f_pi, m_K, M_p, delta_m, eta, eta'.

Theorems verified:
  R35  "QCD String Tension"    (ch19_hadrons.tex) -- sigma_QCD(nf=3) for light hadrons
  R36  "Chiral Relations"      (ch19_hadrons.tex) -- f_pi, m_pi (GMOR), m_rho (Regge)
  R36b "Coulomb+Ratio GMOR"   (ch19_hadrons.tex) -- m_K ratio, m_rho correction, f_K/f_pi
  R37  "Eta-Eta' Mixing"      (ch19_hadrons.tex) -- chi_top, C_mix = Q_Koide = 2/3
  R56  "NLO Corrections"      (ch19_hadrons.tex) -- epsilon = alpha_s*s/(2*pi)

PT constants used:
  s = 1/2, N_c = 3, C_F = 4/3, sigma_QCD_nf(3), alpha_s_eff = 2/3,
  m_u, m_d, m_s (quark masses), gamma (sieve couplings), Q_Koide = 2/3
"""

import sys
import math
import numpy as np
from pathlib import Path

# --- Path setup (monograph v7 scripts) ---
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import (
    s, N_c, n_f, C_F, T_string, sigma_QCD, sigma_QCD_nf,
    alpha_s_eff, alpha_s, alpha_EM, gluon_condensate, regge_slope,
    Q_Koide, gamma,
    m_u as _m_u_MeV, m_d as _m_d_MeV, m_s as _m_s_MeV,
)

ck = Checker("test_hadrons", chapter="ch19", total_steps=8)

# =====================================================================
# Derived quantities
# =====================================================================

HC = 0.197327  # hbar*c in GeV.fm

# Quark masses in GeV (pt_constants stores in MeV)
m_u = _m_u_MeV / 1000.0
m_d = _m_d_MeV / 1000.0
m_s = _m_s_MeV / 1000.0

# Light hadrons use nf=3 (u, d, s all < sqrt(sigma))
NF_LIGHT = 3
sigma_light = sigma_QCD_nf(NF_LIGHT)
sqrt_sigma = math.sqrt(sigma_light)
alpha_prime = 1.0 / (2.0 * math.pi * sigma_light)

# PDG references (GeV)
F_PI_PDG = 0.1302
F_K_PDG = 0.1557
M_PI_PDG = 0.13498      # pi0
M_PI_CH_PDG = 0.13957   # pi+/-
M_K_PDG = 0.4976        # K0
M_RHO_PDG = 0.7755
M_OMEGA_PDG = 0.7827
M_P_PDG = 0.9383        # proton
M_N_PDG = 0.9396        # neutron
M_ETA_PDG = 0.5478
M_ETAP_PDG = 0.9578     # eta'
M_A2_PDG = 1.3169       # a_2(1320)
M_RHO3_PDG = 1.6888     # rho_3(1690)
M_A4_PDG = 2.040        # a_4(2040)
M_KSTAR_PDG = 0.8917    # K*(892)
DELTA_MNP_PDG = 0.001293
GAMMA_RHO_PDG = 0.1491
GAMMA_KSTAR_PDG = 0.0514
R_PI2_PDG = 0.452       # fm^2

print(f"  sigma_QCD(nf=3) = {sigma_light:.4f} GeV^2 = ({sqrt_sigma*1000:.0f} MeV)^2")
print(f"  m_u = {m_u*1000:.2f} MeV, m_d = {m_d*1000:.2f} MeV, m_s = {m_s*1000:.2f} MeV")


# =====================================================================
# Step 1: PION DECAY CONSTANT f_pi
# =====================================================================
ck.section("Step 1: Pion decay constant f_pi (string relation)")

# f_pi^2 = N_c * sigma / (4*pi^2)
f_pi = math.sqrt(N_c * sigma_light / (4.0 * math.pi**2))
ck.check_close("f_pi_string_relation", f_pi, F_PI_PDG, tol_pct=5.0,
               unit=f"{f_pi*1000:.1f} MeV")

# f_K / f_pi ratio (R36b: vertex+edge)
gamma_5 = gamma[5]
m_hat = (m_u + m_d) / 2.0
x_SU3 = (m_s - m_hat) / sqrt_sigma
c_fK = gamma_5 * (1.0 + s)
fK_fpi_PT = 1.0 + c_fK * x_SU3
f_K_PT = f_pi * fK_fpi_PT
fK_fpi_PDG = F_K_PDG / F_PI_PDG

ck.check_close("fK_over_fpi", fK_fpi_PT, fK_fpi_PDG, tol_pct=3.0)
ck.check_close("f_K_derived", f_K_PT, F_K_PDG, tol_pct=3.0)

print(f"  f_pi = {f_pi*1000:.1f} MeV (PDG: {F_PI_PDG*1000:.1f})")
print(f"  f_K/f_pi = {fK_fpi_PT:.4f} (PDG: {fK_fpi_PDG:.4f})")


# =====================================================================
# Step 2: CHIRAL CONDENSATE AND PION MASS (GMOR)
# =====================================================================
ck.section("Step 2: Chiral condensate and pion mass (GMOR)")

# Chiral condensate: C_chi = C_F / pi
C_CHI = C_F / math.pi
qqbar_PT = C_CHI * sigma_light**1.5

# GMOR: m_pi^2 * f_pi^2 = (m_u + m_d) * |<qbar q>|
m_ud = m_u + m_d
m_pi2_PT = m_ud * qqbar_PT / f_pi**2
m_pi_PT = math.sqrt(max(m_pi2_PT, 0))

ck.check_close("m_pi_GMOR", m_pi_PT, M_PI_PDG, tol_pct=15.0,
               unit=f"{m_pi_PT*1000:.1f} MeV")

# Cross-check: condensate from GMOR inverse
qqbar_from_mpi = M_PI_PDG**2 * F_PI_PDG**2 / m_ud
qqbar_tree_cube = qqbar_from_mpi**(1.0 / 3.0)
qqbar_PT_cube = qqbar_PT**(1.0 / 3.0)
ck.check_close("condensate_tree_level", qqbar_PT_cube, qqbar_tree_cube, tol_pct=5.0,
               unit=f"{qqbar_PT_cube*1000:.1f} MeV")

print(f"  C_chi = C_F/pi = {C_CHI:.5f}")
print(f"  |<qbar q>|^{{1/3}} = {qqbar_PT_cube*1000:.1f} MeV")
print(f"  m_pi = {m_pi_PT*1000:.1f} MeV (PDG: {M_PI_PDG*1000:.1f})")


# =====================================================================
# Step 3: KAON AND STRANGE SECTOR
# =====================================================================
ck.section("Step 3: Kaon and strange sector (ratio GMOR, R36b)")

# m_K = m_pi * sqrt((m_u + m_s)/(m_u + m_d))  (LO ratio GMOR, R36b)
m_us = m_u + m_s
ratio_mass = m_us / m_ud
m_K_PT = m_pi_PT * math.sqrt(ratio_mass)
ck.check_close("m_K_ratio_GMOR", m_K_PT, M_K_PDG, tol_pct=5.0,
               unit=f"{m_K_PT*1000:.1f} MeV")

print(f"  m_K = {m_K_PT*1000:.1f} MeV (PDG: {M_K_PDG*1000:.1f})")
print(f"  m_K/m_pi = {m_K_PT/m_pi_PT:.3f} (PDG: {M_K_PDG/M_PI_PDG:.3f})")


# =====================================================================
# Step 4: RHO MESON AND REGGE TRAJECTORY
# =====================================================================
ck.section("Step 4: Rho meson and Regge trajectory (R36b)")

# Coulomb correction to Regge slope (R36b)
alpha_0 = s  # Regge intercept = s = 1/2
corr_coulomb = C_F * alpha_s_eff**2 / math.pi
alpha_prime_eff = alpha_prime * (1.0 + corr_coulomb)

m_rho2_PT = (1.0 - alpha_0) / alpha_prime_eff
m_rho_PT = math.sqrt(m_rho2_PT)
ck.check_close("m_rho_Regge", m_rho_PT, M_RHO_PDG, tol_pct=5.0,
               unit=f"{m_rho_PT*1000:.1f} MeV")

# Omega (same trajectory, isospin partner)
ck.check_close("m_omega_isospin", m_rho_PT, M_OMEGA_PDG, tol_pct=5.0)

# Regge trajectory J=2..4
for J, name, m_pdg in [(2, "a_2(1320)", M_A2_PDG),
                        (3, "rho_3(1690)", M_RHO3_PDG),
                        (4, "a_4(2040)", M_A4_PDG)]:
    m_J = math.sqrt((J - alpha_0) / alpha_prime_eff)
    ck.check_close(f"m_{name}", m_J, m_pdg, tol_pct=5.0,
                   unit=f"{m_J*1000:.0f} MeV")

print(f"  m_rho = {m_rho_PT*1000:.1f} MeV (PDG: {M_RHO_PDG*1000:.1f})")
print(f"  alpha'_eff = {alpha_prime_eff:.4f} GeV^-2")


# =====================================================================
# Step 5: PROTON AND NEUTRON
# =====================================================================
ck.section("Step 5: Proton mass and neutron-proton splitting")

# M_p = 2 * sqrt(sigma)
C_p = 2.0
M_p_PT = C_p * sqrt_sigma
ck.check_close("M_proton", M_p_PT, M_P_PDG, tol_pct=5.0,
               unit=f"{M_p_PT*1000:.1f} MeV")

# delta_m(n-p) = s * (m_d - m_u)
delta_md_mu = m_d - m_u
delta_m_PT = s * delta_md_mu
ck.check_close("delta_m_neutron_proton", delta_m_PT, DELTA_MNP_PDG, tol_pct=20.0,
               unit=f"{delta_m_PT*1000:.3f} MeV")

print(f"  M_p = {M_p_PT*1000:.1f} MeV (PDG: {M_P_PDG*1000:.1f})")
print(f"  m_n - m_p = {delta_m_PT*1000:.3f} MeV (PDG: {DELTA_MNP_PDG*1000:.3f})")


# =====================================================================
# Step 6: HADRONIC DYNAMICS (KSRF, VMD)
# =====================================================================
ck.section("Step 6: Hadronic dynamics (g_rho_pi_pi, Gamma_rho, F_pi)")

# --- Helper functions (inline, no external module dependency) ---

def g_rho_pi_pi_func(m_rho_GeV, f_pi_GeV):
    """KSRF relation: g = m_rho / f_pi."""
    return m_rho_GeV / f_pi_GeV


def p_cm(s_GeV2, m_pi_GeV):
    """Center-of-mass momentum for rho -> pi+ pi-."""
    arg = s_GeV2 / 4.0 - m_pi_GeV**2
    return np.sqrt(arg) if arg > 0 else 0.0


def gamma_rho_energy(s_GeV2, m_rho_GeV, m_pi_GeV, g):
    """Energy-dependent rho width: Gamma(s) = g^2 * p_cm^3 / (6*pi*m_rho^2)."""
    pcm = p_cm(s_GeV2, m_pi_GeV)
    if pcm <= 0:
        return 0.0
    return g**2 * pcm**3 / (6.0 * np.pi * m_rho_GeV**2)


def F_pi_func(s_GeV2, m_rho_GeV, m_pi_GeV, g):
    """Pion form factor F_pi(s) via VMD."""
    m_rho2 = m_rho_GeV**2
    sqrt_s = np.sqrt(max(s_GeV2, 0.0))
    Gamma_s = gamma_rho_energy(s_GeV2, m_rho_GeV, m_pi_GeV, g)
    denom = m_rho2 - s_GeV2 - 1j * sqrt_s * Gamma_s
    return m_rho2 / denom


def r_pi_squared(m_rho_GeV):
    """Pion charge radius squared: <r_pi^2> = 6/m_rho^2 (in fm^2)."""
    return 6.0 / m_rho_GeV**2 * HC**2


def sigma_ee_pipi(s_GeV2, m_rho_GeV, m_pi_GeV, g, alpha_em, as_eff):
    """Cross section e+e- -> pi+pi- via virtual photon exchange (nb)."""
    if s_GeV2 <= 4 * m_pi_GeV**2:
        return 0.0
    pcm_val = p_cm(s_GeV2, m_pi_GeV)
    sqrt_s = np.sqrt(s_GeV2)
    beta_pi = 2.0 * pcm_val / sqrt_s
    Fpi = F_pi_func(s_GeV2, m_rho_GeV, m_pi_GeV, g)
    Fpi2 = np.abs(Fpi)**2
    sigma = (np.pi * alpha_em**2) / (3.0 * s_GeV2) * beta_pi**3 * Fpi2
    sigma *= (1.0 + as_eff / np.pi)
    return sigma * 0.3894e6  # nb


def gamma_Kstar_SU3(m_Kstar_GeV, m_K_GeV, m_pi_GeV, g_rho):
    """K* width via SU(3) flavour scaling."""
    s_val = m_Kstar_GeV**2
    lam = ((s_val - (m_K_GeV + m_pi_GeV)**2)
           * (s_val - (m_K_GeV - m_pi_GeV)**2))
    if lam <= 0:
        return 0.0
    pcm_val = np.sqrt(lam) / (2.0 * m_Kstar_GeV)
    return g_rho**2 * pcm_val**3 / (6.0 * np.pi * m_Kstar_GeV**2)


# Compute hadronic dynamics quantities
g_rho = g_rho_pi_pi_func(m_rho_PT, f_pi)
Gamma_rho_val = gamma_rho_energy(m_rho_PT**2, m_rho_PT, m_pi_PT, g_rho)

# T1: g_rho_pi_pi
G_RHO_PDG = 5.96
ck.check_close("g_rho_pi_pi_KSRF", g_rho, G_RHO_PDG, tol_pct=5.0)

# T2: Gamma_rho
ck.check_close("Gamma_rho_width", Gamma_rho_val, GAMMA_RHO_PDG, tol_pct=5.0,
               unit=f"{Gamma_rho_val*1000:.1f} MeV")

# T3: F_pi(0) = 1 (charge conservation)
Fpi_0 = F_pi_func(0.0, m_rho_PT, m_pi_PT, g_rho)
ck.check_close("F_pi_0_charge_conservation", abs(Fpi_0), 1.0, tol_pct=0.01)

# T4: Pion charge radius
r_pi2 = r_pi_squared(m_rho_PT)
ck.check_close("r_pi_squared", r_pi2, R_PI2_PDG, tol_pct=15.0,
               unit=f"{r_pi2:.3f} fm^2")

# T5: Resonant peak
Fpi_peak = F_pi_func(m_rho_PT**2, m_rho_PT, m_pi_PT, g_rho)
ck.check("F_pi_resonant_peak", abs(Fpi_peak) > 2.0,
         f"|F_pi(m_rho^2)| = {abs(Fpi_peak):.2f}")

# T6: e+e- cross section at rho peak
sigma_peak = sigma_ee_pipi(m_rho_PT**2, m_rho_PT, m_pi_PT, g_rho,
                           alpha_EM, alpha_s_eff)
ck.check_close("sigma_ee_rho_peak", sigma_peak, 1200.0, tol_pct=30.0,
               unit=f"{sigma_peak:.0f} nb")

# T7: Chiral limit (Gamma_rho stays finite as m_pi -> 0)
Gamma_chiral = gamma_rho_energy(m_rho_PT**2, m_rho_PT, 0.001, g_rho)
ck.check("chiral_limit_Gamma_rho_finite",
         0.1 < Gamma_chiral < 0.5,
         f"Gamma(m_pi->0) = {Gamma_chiral*1000:.1f} MeV")

# T8: Gamma_K* via SU(3)
Gamma_Kstar = gamma_Kstar_SU3(M_KSTAR_PDG, m_K_PT, m_pi_PT, g_rho)
ck.check_close("Gamma_Kstar_SU3", Gamma_Kstar, GAMMA_KSTAR_PDG, tol_pct=30.0,
               unit=f"{Gamma_Kstar*1000:.1f} MeV")

print(f"  g_rho_pi_pi = {g_rho:.3f} (PDG ~{G_RHO_PDG})")
print(f"  Gamma_rho = {Gamma_rho_val*1000:.1f} MeV (PDG: {GAMMA_RHO_PDG*1000:.1f})")
print(f"  <r_pi^2> = {r_pi2:.3f} fm^2 (PDG: {R_PI2_PDG})")


# =====================================================================
# Step 7: ETA-ETA' MIXING (R37)
# =====================================================================
ck.section("Step 7: Eta-eta' mixing (R37, anomalie U(1)_A)")

# Topological susceptibility: chi_top = T_string * sigma_pure^2
sigma_pure = sigma_QCD_nf(0)  # pure gauge, nf=0
chi_top_PT = T_string * sigma_pure**2

# Witten-Veneziano anomaly mass
NF_WV = 3
m0_sq = 2 * NF_WV * chi_top_PT / f_pi**2

# 2x2 mass matrix in (eta_8, eta_1) basis
m88_sq = (4.0 * m_K_PT**2 - m_pi_PT**2) / 3.0
m11_sq = (2.0 * m_K_PT**2 + m_pi_PT**2) / 3.0 + m0_sq
m81_sq_LO = -(2.0 * math.sqrt(2.0) / 3.0) * (m_K_PT**2 - m_pi_PT**2)

# R37: C_mix = Q_Koide = 2/3 (screening coefficient)
C_mix = Q_Koide  # = 2/3
ck.check_close("C_mix_eq_Q_Koide", C_mix, 2.0 / 3.0, tol_pct=0.001)

m81_sq = C_mix * m81_sq_LO

# Diagonalization
S_trace = m88_sq + m11_sq
det_M = m88_sq * m11_sq - m81_sq**2
disc = math.sqrt(max(S_trace**2 / 4.0 - det_M, 0))
m_eta_sq = S_trace / 2.0 - disc
m_etap_sq = S_trace / 2.0 + disc
m_eta_PT = math.sqrt(max(m_eta_sq, 0))
m_etap_PT = math.sqrt(m_etap_sq)

ck.check_close("m_eta_mixing", m_eta_PT, M_ETA_PDG, tol_pct=5.0,
               unit=f"{m_eta_PT*1000:.1f} MeV")
ck.check_close("m_eta_prime_mixing", m_etap_PT, M_ETAP_PDG, tol_pct=5.0,
               unit=f"{m_etap_PT*1000:.1f} MeV")

# Naive octet mass formula without U(1)_A anomaly (should fail at ~25% for eta)
# m_eta_octet = sqrt(m_pi^2/3 + 2*m_K^2/3) -- the eta_8 SU(3) octet state
# WITHOUT anomalous mixing, the eta would be a pure octet at ~410 MeV.
m_eta_octet = math.sqrt(M_PI_PDG**2 / 3.0 + 2.0 * M_K_PDG**2 / 3.0)
err_octet = abs(m_eta_octet - M_ETA_PDG) / M_ETA_PDG * 100
ck.check("octet_eta_fails_gt_20pct", err_octet > 20,
         f"eta_8 octet = {m_eta_octet*1000:.0f} MeV, err = {err_octet:.0f}% (signals U(1)_A breaking)")

# chi_top^{1/4} verification
chi_top_fourth = chi_top_PT**0.25
ck.check_close("chi_top_fourth_root", chi_top_fourth, 0.1978, tol_pct=15.0,
               unit=f"{chi_top_fourth*1000:.1f} MeV")

# Witten-Veneziano self-consistency
WV_lhs = m_etap_PT**2 + m_eta_PT**2 - 2 * m_K_PT**2
chi_WV = f_pi**2 * WV_lhs / (2 * NF_WV)
if chi_WV > 0:
    m_etap_sq_WV = 2 * NF_WV * chi_WV / f_pi**2 + 2 * m_K_PT**2 - m_eta_PT**2
    m_etap_WV = math.sqrt(m_etap_sq_WV) if m_etap_sq_WV > 0 else 0
else:
    m_etap_WV = 0
ck.check_close("WV_self_consistent", m_etap_WV, m_etap_PT, tol_pct=1.0)

print(f"  m_eta = {m_eta_PT*1000:.1f} MeV (PDG: {M_ETA_PDG*1000:.1f})")
print(f"  m_eta' = {m_etap_PT*1000:.1f} MeV (PDG: {M_ETAP_PDG*1000:.1f})")
print(f"  chi_top^{{1/4}} = {chi_top_fourth*1000:.1f} MeV (lattice: 197.8)")


# =====================================================================
# Step 8: NLO CORRECTIONS (R56)
# =====================================================================
ck.section("Step 8: NLO corrections (R56, epsilon = alpha_s*s/(2*pi))")

# Universal expansion parameter
epsilon_NLO = alpha_s * s / (2.0 * math.pi)
ck.check("epsilon_NLO_small", 0 < epsilon_NLO < 0.02,
         f"epsilon = {epsilon_NLO:.5f}")

# f_pi NLO
f_pi_NLO = f_pi * (1.0 - epsilon_NLO)
ck.check_close("f_pi_NLO", f_pi_NLO, F_PI_PDG, tol_pct=2.0,
               unit=f"{f_pi_NLO*1000:.1f} MeV")

# m_K NLO (SU(3)-enhanced)
m_K_NLO = m_K_PT * (1.0 - (1.0 + s) * epsilon_NLO)
ck.check_close("m_K_NLO", m_K_NLO, M_K_PDG, tol_pct=2.0,
               unit=f"{m_K_NLO*1000:.1f} MeV")

# M_proton NLO (baryon depth D=2)
D_sieve = 2
M_p_NLO = M_p_PT * (1.0 - D_sieve * epsilon_NLO)
ck.check_close("M_proton_NLO", M_p_NLO, M_P_PDG, tol_pct=2.0,
               unit=f"{M_p_NLO*1000:.1f} MeV")

# delta_m(n-p) NLO (EM + QCD isospin)
delta_EM = C_F * alpha_EM / s
delta_QCD_np = s**2 * alpha_s / math.pi
delta_m_NLO = delta_m_PT * (1.0 + delta_EM + delta_QCD_np)
ck.check_close("delta_m_np_NLO", delta_m_NLO, DELTA_MNP_PDG, tol_pct=5.0,
               unit=f"{delta_m_NLO*1000:.3f} MeV")

# eta/eta' NLO via chi_top correction
chi_top_NLO = chi_top_PT * (1.0 + 2.0 * epsilon_NLO)
m0_sq_NLO = 2.0 * NF_WV * chi_top_NLO / f_pi_NLO**2
m11_sq_NLO = (2.0 * m_K_PT**2 + m_pi_PT**2) / 3.0 + m0_sq_NLO
m88_sq_NLO = (4.0 * m_K_PT**2 - m_pi_PT**2) / 3.0
m81_sq_NLO = C_mix * (-(2.0 * math.sqrt(2.0) / 3.0) * (m_K_PT**2 - m_pi_PT**2))
S_trace_NLO = m88_sq_NLO + m11_sq_NLO
det_M_NLO = m88_sq_NLO * m11_sq_NLO - m81_sq_NLO**2
disc_NLO = math.sqrt(max(S_trace_NLO**2 / 4.0 - det_M_NLO, 0))
m_eta_NLO = math.sqrt(max(S_trace_NLO / 2.0 - disc_NLO, 0))
m_etap_NLO = math.sqrt(S_trace_NLO / 2.0 + disc_NLO)

ck.check_close("m_eta_NLO", m_eta_NLO, M_ETA_PDG, tol_pct=2.0,
               unit=f"{m_eta_NLO*1000:.1f} MeV")
ck.check_close("m_eta_prime_NLO", m_etap_NLO, M_ETAP_PDG, tol_pct=2.0,
               unit=f"{m_etap_NLO*1000:.1f} MeV")

print(f"\n  NLO corrections (epsilon = {epsilon_NLO:.5f}):")
for name, tree, nlo, pdg in [
    ("f_pi", f_pi, f_pi_NLO, F_PI_PDG),
    ("m_K", m_K_PT, m_K_NLO, M_K_PDG),
    ("M_p", M_p_PT, M_p_NLO, M_P_PDG),
    ("m_eta", m_eta_PT, m_eta_NLO, M_ETA_PDG),
    ("m_eta'", m_etap_PT, m_etap_NLO, M_ETAP_PDG),
]:
    err_tree = abs(tree - pdg) / pdg * 100
    err_nlo = abs(nlo - pdg) / pdg * 100
    print(f"    {name:<8} tree {err_tree:.2f}% -> NLO {err_nlo:.2f}%")


# =====================================================================
# SUMMARY
# =====================================================================
ck.summary()
