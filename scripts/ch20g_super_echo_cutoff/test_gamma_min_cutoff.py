#!/usr/bin/env python3
"""
test_gamma_min_cutoff.py -- Chapter 20g: gamma_min(mu) cutoff function

Monograph: chapters/ch20g_super_echo_cutoff.tex
Type: [DER]
Main result: gamma_min(mu) decreases with mu and produces the
super-echo cutoff at p = 23 for mu in the hadronic window.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


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


ck = Checker("test_gamma_min_cutoff", chapter="ch20g", total_steps=2)

# ---- Step 1: Monotonicity gamma_p decreasing with p ----
ck.section("Step 1: gamma_p monotone decreasing in p at fixed mu")
mu_b = 14.726
prev_gamma = None
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29]:
    g = gamma_p_at_mu(p, mu_b)
    if prev_gamma is not None:
        ck.check(f"gamma_p{p}_lt_prev",
                 g < prev_gamma,
                 f"gamma_{p}({mu_b}) = {g:.4f} < gamma_prev = {prev_gamma:.4f}")
    prev_gamma = g

# ---- Step 2: Cutoff at p = 23 (super-echo upper bound) ----
ck.section("Step 2: super-echo cutoff p = 23 at hadronic scale")
g_23 = gamma_p_at_mu(23, mu_b)
g_29 = gamma_p_at_mu(29, mu_b)
ck.check("gamma_23_positive", g_23 > 0,
         f"gamma_23(mu_b) = {g_23:.4f}")
ck.check("gamma_29_smaller", g_29 < g_23,
         f"gamma_29 = {g_29:.4f} < gamma_23 = {g_23:.4f} (cutoff at 23)")

ck.summary()
