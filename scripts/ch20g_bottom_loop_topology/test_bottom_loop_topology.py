#!/usr/bin/env python3
"""
test_bottom_loop_topology.py -- Chapter 20g: Bottom-loop topology

Monograph: chapters/ch20g_bottom_loop_topology.tex
Type: [DER]

Main result:
  The bottom-loop contribution to b -> s gamma is the unique structure:
    Delta C_7^{bottom} = -(1/N_c) * (sin^2 theta_13(q_-) gamma_13)
                                / sin^2 theta_3(q_-) * gamma_7

  with N_c = 3, p = 13 (echo prime for bottom-loop), p_3 = 7 (last active).
"""

from fractions import Fraction


MU_STAR = 15
N_C = 3
Q_MINUS_FLOAT = 2.71828182845904523536 ** (-1.0 / MU_STAR)


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
    """
    PT prediction: Delta C_7^{bottom} = -(1/N_c) (sin^2(13)gamma_13)
                                        / sin^2(3) * gamma_7
    Evaluated on q_- branch.
    """
    q = Q_MINUS_FLOAT
    s2_13 = sin2_theta_p(13, q)
    g_13 = gamma_p(13, q)
    s2_3 = sin2_theta_p(3, q)
    g_7 = gamma_p(7, q)
    return -(1.0 / N_C) * (s2_13 * g_13) / s2_3 * g_7


def main():
    print("=" * 70)
    print("Chapter 20g: Bottom-loop topology for b -> s gamma")
    print("=" * 70)
    print()
    print(f"PT structure: Delta C_7^{{bottom}} = -(1/N_c) (sin^2(13) gamma_13)")
    print(f"                                  / sin^2(3) * gamma_7")
    print()
    print(f"Evaluated on q_- branch (hadronic, edge-channel).")
    print()

    q = Q_MINUS_FLOAT
    s2_13 = sin2_theta_p(13, q)
    g_13 = gamma_p(13, q)
    s2_3 = sin2_theta_p(3, q)
    g_7 = gamma_p(7, q)

    print(f"Component values at mu* = {MU_STAR}, q_- = {q:.5f}:")
    print(f"  sin^2(theta_13, q_-) = {s2_13:.6f}")
    print(f"  gamma_13            = {g_13:.6f}")
    print(f"  sin^2(theta_3, q_-)  = {s2_3:.6f}")
    print(f"  gamma_7             = {g_7:.6f}")
    print()

    dC7 = delta_C7_bottom()
    print(f"Delta C_7^{{bottom}} = {dC7:.6f}")
    print()
    print(f"This contribution combines with the dominant W-loop and")
    print(f"top-loop diagrams to give the observed b -> s gamma branching.")


if __name__ == "__main__":
    main()
