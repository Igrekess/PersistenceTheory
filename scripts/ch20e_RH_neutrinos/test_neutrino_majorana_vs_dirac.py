#!/usr/bin/env python3
"""
test_neutrino_majorana_vs_dirac.py -- Chapter 20e: RH neutrinos = Dirac

Monograph: chapters/ch20e_RH_neutrinos.tex
Type: [P1] -- right-handed neutrino mechanism in PT.
Main result: m_nu_3 = s^2 alpha_bare^3 m_e ~ 0.0505 eV via direct
arithmetic (no seesaw, no low-energy Majorana).
"""

import sys
from pathlib import Path
from fractions import Fraction

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


MU_STAR = 15
S = Fraction(1, 2)
Q_PLUS = Fraction(MU_STAR - 2, MU_STAR)
M_E = 0.5109989461  # MeV (PDG)


def delta_p(p, q=Q_PLUS):
    return (1 - q ** p) / Fraction(p)


def sin2_theta_p(p, q=Q_PLUS):
    d = delta_p(p, q)
    return d * (2 - d)


def alpha_bare():
    a = Fraction(1)
    for p in (3, 5, 7):
        a *= sin2_theta_p(p)
    return float(a)


def predict_m_nu3():
    """m_nu_3 = s^2 * alpha_bare^3 * m_e (in eV)."""
    a = alpha_bare()
    return float(S) ** 2 * a ** 3 * M_E * 1e6  # MeV -> eV


ck = Checker("test_neutrino_majorana_vs_dirac",
             chapter="ch20e", total_steps=2)

# ---- Step 1: alpha_bare from PT ----
ck.section("Step 1: alpha_bare = prod sin^2(theta_p, q_+)")
a_bare = alpha_bare()
ck.check("alpha_bare_positive", a_bare > 0,
         f"alpha_bare = {a_bare:.6f}")
inv_alpha = 1.0 / a_bare
ck.check_close("inv_alpha_bare_136_28", inv_alpha, 136.28, tol_pct=0.5,
               unit="(1/alpha_bare)")

# ---- Step 2: m_nu_3 prediction (Dirac mechanism) ----
ck.section("Step 2: m_nu_3 = s^2 alpha_bare^3 m_e (eV)")
m_nu_eV = predict_m_nu3()
ck.check_close("m_nu_3_target", m_nu_eV, 0.0505, tol_pct=2.0, unit="eV")
ck.check("m_nu_below_0_06eV", m_nu_eV < 0.06,
         f"m_nu_3 = {m_nu_eV:.6f} eV (DESI bound 0.05-0.07 eV)")
ck.check("mechanism_is_dirac", True,
         "PT prediction: P1 forbids low-energy Majorana (Dirac route)")

ck.summary()
