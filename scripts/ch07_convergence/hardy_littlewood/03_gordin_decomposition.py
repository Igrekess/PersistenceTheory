#!/usr/bin/env python3
"""
PT_HL -- Section 4 : Decomposition de Gordin
==============================================

LEMME 2 (Gordin 1969): La marche S_n se decompose en:
  S_n = M_n + R_n
ou:
  d_n = [(2*T_12 - 1)*chi_n + chi_{n+1}] / (2*T_12)
  M_n = sum d_i (quasi-martingale)
  |R_n| <= 1/T_12 <= 2
  Q_N = sum d_i^2 = sigma^2 * (N-1)  avec sigma^2 = (1-T_12)/T_12

Il suffit de montrer r_K = max|M|^2 / Q < 1.

Author: yan senez
Date: March 2026
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from infrastructure_00 import *


def main():
    print_header("SECTION 4: DECOMPOSITION DE GORDIN")
    print("""
  S_n = M_n + R_n (Gordin 1969)
  d_n = [(2T_12-1)*chi_n + chi_{n+1}] / (2T_12)
  Il suffit de montrer: r_K = max|M|^2 / Q < 1.
    """)

    tests = {}

    levels = {}
    for k in range(3, 9):
        levels[k] = compute_level(k)

    # 1. Verifier |R_n| <= 2
    print("  1. Reste borne: |R_n| <= 2\n")
    print(f"  {'K':>3} {'max|R|':>10} {'<= 2?':>8} {'1/T_12':>10}")

    for k in range(3, 9):
        lk = levels[k]
        bound = 1.0 / lk['T_12']
        ok = lk['max_R'] <= 2.0 + 1e-12
        tests[f'T3.{k}.1 |R|<=2'] = ok
        print(f"  {k:>3} {lk['max_R']:>10.6f} {'OUI' if ok else 'NON':>8} "
              f"{bound:>10.6f}")

    # 2. Verifier Q = sigma^2 * (N-1) avec sigma^2 = (1-T_12)/T_12
    print(f"\n  2. Identite energetique: Q = sigma^2 * (N-1)\n")
    print(f"  {'K':>3} {'Q':>12} {'sigma^2*(N-1)':>16} {'erreur rel.':>12}")

    for k in range(3, 9):
        lk = levels[k]
        sigma_sq = (1 - lk['T_12']) / lk['T_12']
        Q_pred = sigma_sq * (lk['N_K'] - 1)
        err = abs(lk['Q'] - Q_pred) / lk['Q'] if lk['Q'] > 0 else 0
        ok = err < 1e-10
        tests[f'T3.{k}.2 Q=sigma^2*(N-1)'] = ok
        print(f"  {k:>3} {lk['Q']:>12.2f} {Q_pred:>16.2f} {err:>12.2e}")

    # 3. Ratio de Gordin r_K = max|M|^2 / Q
    print(f"\n  3. Ratio de Gordin r_K = max|M|^2 / Q\n")
    print(f"  {'K':>3} {'max|M|':>10} {'Q':>12} {'r_K':>14} {'< 1?':>6}")

    gordin_r = []
    for k in range(3, 9):
        lk = levels[k]
        r = lk['r_K']
        gordin_r.append(r)
        ok = r < 1.0
        tests[f'T3.{k}.3 r_K<1'] = ok
        print(f"  {k:>3} {lk['max_M']:>10.3f} {lk['Q']:>12.1f} "
              f"{r:>14.6f} {'OUI' if ok else 'NON':>6}")

    # 4. Decroissance monotone
    decr = all(gordin_r[i+1] < gordin_r[i] for i in range(len(gordin_r)-1))
    tests['T3.dec r_K decroissant'] = decr
    contraction = gordin_r[0] / gordin_r[-1] if gordin_r[-1] > 0 else float('inf')
    print(f"\n  4. Decroissance monotone: {'OUI' if decr else 'NON'}")
    print(f"     Contraction totale: {contraction:.0f}x "
          f"(de {gordin_r[0]:.4f} a {gordin_r[-1]:.6f})")

    # 5. sigma^2 -> 1/2
    print(f"\n  5. Convergence sigma^2 -> 1/2\n")
    print(f"  {'K':>3} {'sigma^2':>10} {'|sigma^2-1/2|':>14}")

    for k in range(3, 9):
        lk = levels[k]
        sigma_sq = (1 - lk['T_12']) / lk['T_12']
        print(f"  {k:>3} {sigma_sq:>10.6f} {abs(sigma_sq-0.5):>14.6f}")

    return print_summary(tests, "SECTION 4: GORDIN")


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
