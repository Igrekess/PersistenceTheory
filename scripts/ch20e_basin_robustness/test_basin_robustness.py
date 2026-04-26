#!/usr/bin/env python3
"""
test_basin_robustness.py -- Chapter 20e: basin robustness of mu* = 15

Monograph: chapters/ch20e_basin_robustness.tex
Type: [THM + VAL]
Main result: basin radius Delta_mu_0 <= 2 from SM baseline mu* = 15.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


MU_STAR_SM = 15
ACTIVE_PRIMES_SM = (3, 5, 7)
THRESHOLD = 0.5


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


def primes_below(n):
    sieve = [True] * n
    sieve[0:2] = [False, False]
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, n, i):
                sieve[j] = False
    return [i for i in range(n) if sieve[i]]


def iterate_fixed_point(mu_0, max_prime=251, max_iter=20):
    """
    PT iteration: mu_{k+1} = sum {p in P>=3 : gamma_p(mu_k) > 1/2}.
    The prime p=2 is excluded by construction (info/anti-info channel,
    not part of the residue lattice; cf. monograph ch01_sieve.tex
    rem:p2_self_critical).
    """
    primes_list = [p for p in primes_below(max_prime + 1) if p >= 3]
    mu = mu_0
    for k in range(max_iter):
        active = [p for p in primes_list if gamma_p_at_mu(p, mu) > THRESHOLD]
        mu_new = sum(active) if active else 0
        if mu_new == mu:
            return mu, k, active
        mu = mu_new
    return None, max_iter, []


ck = Checker("test_basin_robustness", chapter="ch20e", total_steps=3)

# ---- Step 1: SM baseline check ----
ck.section("Step 1: SM baseline mu* = 15")
ck.check("active_primes_sum", sum(ACTIVE_PRIMES_SM) == MU_STAR_SM,
         f"3+5+7 = {sum(ACTIVE_PRIMES_SM)} = mu*")

# ---- Step 2: Iterate from mu_0 in basin ----
ck.section("Step 2: basin convergence Delta_mu_0 in [-2, +2]")
in_basin_results = []
for delta in [-2, -1, 0, 1, 2]:
    mu_0 = MU_STAR_SM + delta
    mu_inf, k, active = iterate_fixed_point(mu_0)
    in_basin_results.append((delta, mu_inf, k))
    print(f"  Delta_mu_0={delta:+d}, mu_0={mu_0}, mu_inf={mu_inf}, iter={k}")

all_converge = all(r[1] == MU_STAR_SM for r in in_basin_results)
ck.check("basin_radius_2_converges_to_15", all_converge,
         f"all Delta_mu_0 in [-2, +2] return to mu* = 15")

# ---- Step 3: Outside basin (Delta = +3) ----
ck.section("Step 3: outside basin (Delta_mu_0 = +3)")
mu_inf_3, k_3, active_3 = iterate_fixed_point(MU_STAR_SM + 3)
print(f"  Delta_mu_0=+3, mu_0=18, mu_inf={mu_inf_3}, iter={k_3}")
ck.check("delta_3_does_not_return_to_15",
         mu_inf_3 != MU_STAR_SM,
         f"mu_0=18 -> mu_inf={mu_inf_3} (not 15, basin exit confirmed)")

ck.summary()
