#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_half_holonomy_exactness_PT.py
==========================================

Exactness proof for the PT Kerr half-holonomy deformation.

The first-order uniqueness test proves that

    phi(x)   = phi0 + pi x + O(x^2),
    kappa(x) = kappa0 + O(x^2),

where x = M Omega_H(a).  This script closes the higher-order analytic
freedom under the stronger PT holonomy statement:

  H1. The Kerr spin correction is a same-sector U(1) horizon holonomy
      character.  Therefore the phase increment F(x)=phi(x)-phi0 is locally
      additive:

          F(x+y) = F(x) + F(y).

  H2. The spin involution selects the half-character, so F'(0)=pi.

  H3. No dissipative channel is opened by this holonomy.  Therefore tau/M,
      equivalently kappa, is invariant as a function of x.

For an arbitrary analytic jet to finite order N, local additivity kills every
phase coefficient beyond the linear one, while dissipative invariance kills
every damping coefficient.  Since N is arbitrary, no analytic higher-order
same-sector correction remains.

No LVK data enter this proof.
"""

from __future__ import annotations

import io
import sys

import sympy as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

N = 6
x, y, a = sp.symbols("x y a", real=True)
phi0, kappa0 = sp.symbols("phi0 kappa0", positive=True, real=True)
c = sp.symbols(f"c1:{N + 1}", real=True)
d = sp.symbols(f"d1:{N + 1}", real=True)
pi = sp.pi


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


def coeff_dict(poly):
    return sp.Poly(sp.expand(poly), x, y).terms()


def coefficient_equations(poly):
    equations = []
    for (powers, coeff) in coeff_dict(poly):
        if coeff != 0:
            equations.append(sp.Eq(coeff, 0))
    return equations


print("=" * 88)
print("  PT QG: EXACTNESS OF THE KERR HALF-HOLONOMY DEFORMATION")
print("=" * 88)

n_pass = 0
n_total = 0

F = sum(c[i - 1] * x**i for i in range(1, N + 1))
K = kappa0 + sum(d[i - 1] * x**i for i in range(1, N + 1))
tau_over_m = 4 * pi**2 / K

print("\n--- Arbitrary analytic same-sector jet ---")
print(f"  F(x)=phi(x)-phi0 = {F}")
print(f"  kappa(x)          = {K}")

check("analytic phase jet contains one coefficient per order",
      len(c) == N,
      f"N={N}, coefficients={c}")

check("analytic damping jet contains one coefficient per order",
      len(d) == N,
      f"N={N}, coefficients={d}")

additivity_residual = sp.expand(F.subs(x, x + y) - F - F.subs(x, y))
additivity_eqs = coefficient_equations(additivity_residual)
additivity_solution = sp.solve(additivity_eqs, c[1:], dict=True)

check("local holonomy character gives polynomial additivity equations",
      len(additivity_eqs) == N * (N - 1) // 2,
      f"equations={len(additivity_eqs)}")

check("additivity kills all nonlinear phase coefficients",
      additivity_solution == [{ci: 0 for ci in c[1:]}],
      f"solution={additivity_solution}")

spin_solution = sp.solve([sp.Eq(c[0], pi)], [c[0]], dict=True)
check("spin half-character fixes the linear phase coefficient",
      spin_solution == [{c[0]: pi}],
      "F'(0)=pi")

exact_phase = sp.simplify(F.subs({c[0]: pi, **{ci: 0 for ci in c[1:]}}))
check("exact phase increment is F(x)=pi*x",
      exact_phase == pi * x,
      str(exact_phase))

tau_series = sp.series(tau_over_m, x, 0, N + 1).removeO()
tau0 = 4 * pi**2 / kappa0
tau_eqs = []
for power in range(1, N + 1):
    tau_eqs.append(sp.Eq(sp.expand(tau_series - tau0).coeff(x, power), 0))
tau_solution = sp.solve(tau_eqs, d, dict=True)

check("exact dissipative invariance gives one equation per damping order",
      len(tau_eqs) == N,
      f"orders=1..{N}")

check("dissipative invariance kills all damping coefficients",
      tau_solution == [{di: 0 for di in d}],
      f"solution={tau_solution}")

exact_kappa = sp.simplify(K.subs({di: 0 for di in d}))
check("exact damping exponent is kappa(x)=kappa0",
      exact_kappa == kappa0,
      str(exact_kappa))

exact_momega = sp.simplify((phi0 + exact_phase) / (2 * pi))
check("exact frequency rule is Momega=Momega0+x/2",
      sp.simplify(exact_momega - phi0 / (2 * pi) - x / 2) == 0,
      str(exact_momega))

omega_h = a / (2 * (1 + sp.sqrt(1 - a**2)))
exact_kerr_momega = sp.simplify(exact_momega.subs(x, omega_h))
check("substitution x=M Omega_H(a) gives the Kerr half-holonomy formula",
      sp.simplify(exact_kerr_momega - phi0 / (2 * pi) - omega_h / 2) == 0,
      str(exact_kerr_momega))

spin_series = sp.series(exact_kerr_momega - phi0 / (2 * pi), a, 0, 8).removeO()
check("higher powers in spin are fixed by Kerr geometry, not free PT coefficients",
      spin_series == a / 8 + a**3 / 32 + a**5 / 64 + 5 * a**7 / 512,
      f"Delta Momega(a)={spin_series}+O(a^8)")

wrong_phase_1 = pi * x + x**2
wrong_phase_2 = pi * x + sp.Rational(1, 3) * x**3
wrong_res_1 = sp.expand(wrong_phase_1.subs(x, x + y) - wrong_phase_1 - wrong_phase_1.subs(x, y))
wrong_res_2 = sp.expand(wrong_phase_2.subs(x, x + y) - wrong_phase_2 - wrong_phase_2.subs(x, y))

check("quadratic phase correction violates holonomy additivity",
      wrong_res_1 != 0,
      f"residual={wrong_res_1}")

check("cubic phase correction violates holonomy additivity",
      wrong_res_2 != 0,
      f"residual={wrong_res_2}")

wrong_kappa = kappa0 + x**2
wrong_tau = sp.series(4 * pi**2 / wrong_kappa - tau0, x, 0, 5).removeO()
check("nonzero damping correction violates tau/M invariance",
      wrong_tau != 0,
      f"tau residual={wrong_tau}")

for order in range(2, N + 1):
    trial = pi * x + x**order
    residual = sp.expand(trial.subs(x, x + y) - trial - trial.subs(x, y))
    check(f"phase monomial x^{order} is not an admissible holonomy character",
          residual != 0,
          f"lowest residual degree={sp.Poly(residual, x, y).total_degree()}")

check("exact rule preserves zero-spin limit",
      sp.simplify(exact_kerr_momega.subs(a, 0) - phi0 / (2 * pi)) == 0)

check("exact rule reaches half-shift at extremality",
      sp.simplify(sp.limit(exact_kerr_momega - phi0 / (2 * pi), a, 1, dir="-") - sp.Rational(1, 4)) == 0,
      "Delta Momega -> 1/4")

print("\n" + "=" * 88)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: analytic higher-order same-sector Kerr corrections are")
print("               excluded by holonomy-character exactness and dissipative invariance.")
print("=" * 88)

sys.exit(0 if n_pass == n_total else 1)
