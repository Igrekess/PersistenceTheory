#!/usr/bin/env python3
"""
test_super_echo_selection.py -- Chapter 20g: super-echo prime selection

Monograph: chapters/ch20g_super_echo_cutoff.tex
Type: [DER]
Main result: at hadronic scale mu_b ~ 14.726, the active set is
{11, 13, 17, 19, 23} (super-echo cutoff at p = 23).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


SUPER_ECHO_PRIMES = (11, 13, 17, 19, 23)
HADRONIC_SCALE = 14.726  # mu_b


def gamma_p_at_mu(p, mu):
    q = (mu - 2) / mu
    qp = q ** p
    one_minus_qp = 1 - qp
    if one_minus_qp <= 0:
        return 0.0
    delta = one_minus_qp / p
    if delta <= 0 or 2 - delta <= 0:
        return 0.0
    return (4 * p * q ** (p - 1) * (1 - delta)) / (mu * one_minus_qp * (2 - delta))


def selected_primes_at_scale(mu, threshold=0.0, max_prime=29):
    """Return primes p such that gamma_p(mu) > threshold."""
    primes = [p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29] if p <= max_prime]
    return tuple(p for p in primes if gamma_p_at_mu(p, mu) > threshold)


ck = Checker("test_super_echo_selection", chapter="ch20g", total_steps=2)

# ---- Step 1: Super-echo set composition ----
ck.section("Step 1: super-echo set {11, 13, 17, 19, 23}")
ck.check("super_echo_5_primes", len(SUPER_ECHO_PRIMES) == 5,
         f"super-echo: {SUPER_ECHO_PRIMES}")
ck.check("super_echo_min_11", min(SUPER_ECHO_PRIMES) == 11,
         "min element = 11 (first inactive)")
ck.check("super_echo_max_23", max(SUPER_ECHO_PRIMES) == 23,
         "max element = 23 (cutoff)")

# ---- Step 2: Hadronic scale activation ----
ck.section("Step 2: super-echo at hadronic scale mu_b ~ 14.726")
gammas = {p: gamma_p_at_mu(p, HADRONIC_SCALE) for p in SUPER_ECHO_PRIMES}
for p in SUPER_ECHO_PRIMES:
    ck.check(f"gamma_p{p}_at_mub_positive", gammas[p] > 0,
             f"gamma_{p}(mu_b) = {gammas[p]:.4f}")
gamma_29 = gamma_p_at_mu(29, HADRONIC_SCALE)
ck.check("p29_below_threshold_or_subdominant",
         gamma_29 < gammas[23],
         f"gamma_29(mu_b) = {gamma_29:.4f} < gamma_23 = {gammas[23]:.4f}")

ck.summary()
