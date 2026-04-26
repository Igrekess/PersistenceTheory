#!/usr/bin/env python3
"""
test_wilson_coefficients.py -- Chapter 20b: Wilson coefficient identity

Monograph: chapters/ch20b_wilson_coefficients.tex
Derivation chain: s = 1/2 -> {3,5,7} -> echo primes {11,13} -> beta_echo
                   -> identity |C_9 - C_10| = beta_echo
Zero fitted parameters.

Main result of ch20b:

    |C_9 - C_10|_{universal} = beta_echo = 0.1039 (3.92%)

where beta_echo = sum over echo primes p in {11,13} of
    sin^2(theta_p, q_+) * gamma_p

This is a [PRED] identity falsifiable via b -> s ell+ ell- decay
spectroscopy at HL-LHC and FCC-hh.
"""

from fractions import Fraction
import math


# ---------------------------------------------------------------
# 1. PT constants at the fixed point mu* = 15
# ---------------------------------------------------------------
MU_STAR = 15
S = Fraction(1, 2)
Q_PLUS = Fraction(MU_STAR - 2, MU_STAR)   # q_+ = 13/15

ECHO_PRIMES = [11, 13]


def delta_p(p, q=Q_PLUS):
    """delta_p = (1 - q^p) / p  (rational at q = 13/15)."""
    return (1 - q ** p) / Fraction(p)


def sin2_theta_p(p, q=Q_PLUS):
    """sin^2(theta_p) = delta_p (2 - delta_p)."""
    d = delta_p(p, q)
    return d * (2 - d)


def gamma_p(p, q=Q_PLUS, mu=MU_STAR):
    """
    gamma_p (anomalous dimension) at mu*=15 by analytical formula:
        gamma_p = -d ln(sin^2 theta_p) / d ln mu
                = 4 p q^(p-1) (1 - delta_p) / (mu (1 - q^p) (2 - delta_p))
    where q = (mu - 2)/mu.
    """
    q_f = float(q)
    d = float(delta_p(p, q))
    s2 = float(sin2_theta_p(p, q))
    # Use the analytical derivative
    one_minus_qp = 1 - q_f ** p
    return (4 * p * q_f ** (p - 1) * (1 - d)) / (mu * one_minus_qp * (2 - d))


# ---------------------------------------------------------------
# 2. Echo coupling beta_echo on echo primes {11, 13}
# ---------------------------------------------------------------
def beta_echo():
    """
    beta_echo = sum_{p in {11,13}} sin^2(theta_p, q_+) * gamma_p

    This is the universal echo coupling at the fixed point.
    """
    total = 0.0
    for p in ECHO_PRIMES:
        s2p = float(sin2_theta_p(p))
        gp = gamma_p(p)
        contrib = s2p * gp
        total += contrib
    return total


# ---------------------------------------------------------------
# 3. PT prediction for |C_9 - C_10|_{universal}
# ---------------------------------------------------------------
def predict_C9_minus_C10():
    """
    PT identity: |C_9 - C_10|_{universal} = beta_echo
    """
    return beta_echo()


# ---------------------------------------------------------------
# 4. Main verification
# ---------------------------------------------------------------
def main():
    print("=" * 70)
    print("Chapter 20b: Wilson coefficient identity |C_9 - C_10| = beta_echo")
    print("=" * 70)
    print(f"\nPT constants at mu* = {MU_STAR}:")
    print(f"  q_+ = {Q_PLUS} = {float(Q_PLUS):.6f}")
    print()

    print("Echo prime contributions:")
    print(f"  {'p':>4} {'sin^2(theta_p, q+)':>20} {'gamma_p':>12} {'product':>12}")
    print(f"  {'-'*4} {'-'*20} {'-'*12} {'-'*12}")
    for p in ECHO_PRIMES:
        s2p = float(sin2_theta_p(p))
        gp = gamma_p(p)
        prod = s2p * gp
        print(f"  {p:>4} {s2p:>20.8f} {gp:>12.8f} {prod:>12.8f}")

    beta = beta_echo()
    print(f"\nbeta_echo = sum = {beta:.6f}")
    print()

    # Target from monograph (3.92% for the |C9-C10| identity)
    target_value = 0.1039
    error = abs(beta - target_value) / target_value * 100
    print(f"PT prediction: |C_9 - C_10|_universal = {beta:.5f}")
    print(f"Monograph value (ch20b):                {target_value:.5f}")
    print(f"Discrepancy:                            {error:.2f}%")

    if error < 5.0:
        status = "PASS"
    else:
        status = "REVIEW"
    print(f"\nStatus: {status}")
    print()
    print("This identity is falsifiable via b -> s ell+ ell- decay")
    print("spectroscopy at HL-LHC and FCC-hh (per ch20b [PRED]).")


if __name__ == "__main__":
    main()
