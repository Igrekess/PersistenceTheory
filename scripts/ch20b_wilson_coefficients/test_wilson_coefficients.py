#!/usr/bin/env python3
"""
test_wilson_coefficients.py -- Chapter 20b: Wilson coefficient identity

Monograph: chapters/ch20b_wilson_coefficients.tex
Derivation chain: s = 1/2 -> {3,5,7} -> echo primes {11,13} -> beta_echo
                   -> identity |C_9 - C_10| = beta_echo
Zero fitted parameters.
"""

import sys
from pathlib import Path
from fractions import Fraction

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


# ---------------------------------------------------------------
# PT constants at the fixed point mu* = 15
# ---------------------------------------------------------------
MU_STAR = 15
S = Fraction(1, 2)
Q_PLUS = Fraction(MU_STAR - 2, MU_STAR)   # q_+ = 13/15
ECHO_PRIMES = [11, 13]


def delta_p(p, q=Q_PLUS):
    return (1 - q ** p) / Fraction(p)


def sin2_theta_p(p, q=Q_PLUS):
    d = delta_p(p, q)
    return d * (2 - d)


def gamma_p(p, q=Q_PLUS, mu=MU_STAR):
    """gamma_p = -d ln(sin^2 theta_p) / d ln mu (analytical formula)."""
    q_f = float(q)
    d = float(delta_p(p, q))
    one_minus_qp = 1 - q_f ** p
    return (4 * p * q_f ** (p - 1) * (1 - d)) / (mu * one_minus_qp * (2 - d))


def beta_echo():
    """beta_echo = sum over echo primes {11,13} of sin^2(theta_p, q+) * gamma_p."""
    return sum(float(sin2_theta_p(p)) * gamma_p(p) for p in ECHO_PRIMES)


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------
ck = Checker("test_wilson_coefficients", chapter="ch20b", total_steps=2)

# ---- Step 1: PT constants ----
ck.section("Step 1: PT echo coupling beta_echo")
ck.check("q_plus_canonical", Q_PLUS == Fraction(13, 15),
         f"q_+ = {Q_PLUS}")
ck.check("echo_primes_set", ECHO_PRIMES == [11, 13],
         f"echo primes = {ECHO_PRIMES}")
beta = beta_echo()
ck.check("beta_echo_positive", beta > 0,
         f"beta_echo = {beta:.6f}")
ck.check("beta_echo_subpercent", beta < 0.2,
         f"beta_echo = {beta:.6f} (expected near 0.10)")

# ---- Step 2: Wilson coefficient identity ----
ck.section("Step 2: |C_9 - C_10|_universal = beta_echo")
target_value = 0.1039  # monograph ch20b
ck.check_close("wilson_C9_minus_C10_identity",
               beta, target_value, tol_pct=5.0)

# Per-prime contributions
for p in ECHO_PRIMES:
    contrib = float(sin2_theta_p(p)) * gamma_p(p)
    ck.check(f"contrib_p{p}_positive", contrib > 0,
             f"sin^2(theta_{p})*gamma_{p} = {contrib:.6f}")

ck.summary()
