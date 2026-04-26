#!/usr/bin/env python3
"""
test_bottom_loop_topology.py -- Chapter 20g: bottom-loop topology b->s gamma

Monograph: chapters/ch20g_bottom_loop_topology.tex
Type: [DER]
Delta C_7^bottom = -(1/N_c) (sin^2(13) gamma_13) / sin^2(3) * gamma_7
"""

import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


MU_STAR = 15
N_C = 3
Q_MINUS = math.exp(-1.0 / MU_STAR)


def delta_p(p, q):
    return (1 - q ** p) / p


def sin2_theta_p(p, q):
    d = delta_p(p, q)
    return d * (2 - d)


def gamma_p(p, q, mu=MU_STAR):
    q_f = float(q)
    d = float(delta_p(p, q_f))
    one_minus_qp = 1 - q_f ** p
    return (4 * p * q_f ** (p - 1) * (1 - d)) / (mu * one_minus_qp * (2 - d))


def delta_C7_bottom():
    q = Q_MINUS
    s2_13 = sin2_theta_p(13, q)
    g_13 = gamma_p(13, q)
    s2_3 = sin2_theta_p(3, q)
    g_7 = gamma_p(7, q)
    return -(1.0 / N_C) * (s2_13 * g_13) / s2_3 * g_7


ck = Checker("test_bottom_loop_topology", chapter="ch20g", total_steps=2)

# ---- Step 1: Component values on q_- branch ----
ck.section("Step 1: q_- branch components (hadronic)")
q = Q_MINUS
s2_13 = sin2_theta_p(13, q)
g_13 = gamma_p(13, q)
s2_3 = sin2_theta_p(3, q)
g_7 = gamma_p(7, q)
ck.check("sin2_p13_positive", s2_13 > 0,
         f"sin^2(theta_13, q_-) = {s2_13:.6f}")
ck.check("gamma_13_positive", g_13 > 0,
         f"gamma_13 = {g_13:.6f}")
ck.check("sin2_p3_positive", s2_3 > 0,
         f"sin^2(theta_3, q_-) = {s2_3:.6f}")
ck.check("gamma_7_positive", g_7 > 0,
         f"gamma_7 = {g_7:.6f}")

# ---- Step 2: Delta C_7^bottom amplitude ----
ck.section("Step 2: Delta C_7^bottom amplitude")
dC7 = delta_C7_bottom()
ck.check("dC7_negative", dC7 < 0,
         f"Delta C_7^bottom = {dC7:.6f} (negative as expected)")
ck.check("dC7_subunit", abs(dC7) < 1.0,
         f"|Delta C_7^bottom| = {abs(dC7):.6f} (perturbative)")
ck.check("dC7_above_threshold", abs(dC7) > 1e-4,
         f"|Delta C_7^bottom| measurable")

ck.summary()
