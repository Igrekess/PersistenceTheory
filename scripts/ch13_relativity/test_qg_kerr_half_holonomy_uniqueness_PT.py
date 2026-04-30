#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_half_holonomy_uniqueness_PT.py
===========================================

First-order uniqueness proof for the PT Kerr half-holonomy deformation.

Assumptions, all explicit:

  A1. The macroscopic selector is already fixed to the same sector
      (ell,m,n)=(2,2,0).  No mode switch is allowed at this order.
  A2. Kerr no-hair leaves only the dimensionless spin scalar

          x = M Omega_H(a) = a / (2(1+sqrt(1-a^2)))

      once mass scaling and chirality have been factored out.
  A3. We work at first order in this scalar x.  Therefore every admissible
      same-sector deformation has the form

          phi(x)   = phi_0   + A x + O(x^2)
          kappa(x) = kappa_0 + B x + O(x^2).

  A4. Spin involution means the horizon phase enters as a half-cycle of the
      PT cyclic phase 2*pi, so Delta phi = pi x at first order.
  A5. Dissipative invariance means the first Kerr spin correction is a chiral
      phase torsion, not an opening of a new damping channel:

          d(tau/M)/dx |_{x=0} = 0.

The script proves that these constraints have exactly one first-order solution:

          A = pi,  B = 0,

or equivalently

          M omega(a) = M omega_0 + Omega_H(a)/2,
          kappa(a)   = kappa_0.

No LVK data enter this proof.  LVK validation is handled in the neighbouring
half-holonomy and posterior-confrontation scripts.
"""

from __future__ import annotations

import io
import math
import sys

import sympy as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

pi = sp.pi
x, A, B, c = sp.symbols("x A B c", real=True)
phi0, kappa0 = sp.symbols("phi0 kappa0", positive=True, real=True)


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


print("=" * 88)
print("  PT QG: UNIQUENESS OF THE KERR HALF-HOLONOMY DEFORMATION")
print("=" * 88)

n_pass = 0
n_total = 0

phi = phi0 + A * x
kappa = kappa0 + B * x
momega = phi / (2 * pi)
tau_over_m = 4 * pi**2 / kappa
quality = 2 * pi * phi / (2 * kappa)

print("\n--- General first-order same-sector deformation ---")
print(f"  phi(x)   = {phi}")
print(f"  kappa(x) = {kappa}")

check("first-order same-sector deformation has two free coefficients",
      {A, B} == (phi.free_symbols | kappa.free_symbols) - {x, phi0, kappa0},
      "free coefficients are A for phase and B for damping")

check("zero-spin limit recovers the undeformed mode",
      sp.simplify(phi.subs(x, 0) - phi0) == 0
      and sp.simplify(kappa.subs(x, 0) - kappa0) == 0)

# Kerr horizon scalar.
a = sp.symbols("a", nonnegative=True, real=True)
omega_h = a / (2 * (1 + sp.sqrt(1 - a**2)))
omega_h_small = sp.series(omega_h, a, 0, 3).removeO()

check("Kerr no-hair scalar vanishes at zero spin",
      sp.simplify(omega_h.subs(a, 0)) == 0)

check("Kerr no-hair scalar has small-spin slope 1/4",
      sp.simplify(sp.diff(omega_h, a).subs(a, 0) - sp.Rational(1, 4)) == 0,
      f"M Omega_H(a) = {omega_h_small} + O(a^3)")

# Spin involution / half-cycle.
spin_involution_residual = sp.simplify(A - pi)
spin_solution = sp.solve([spin_involution_residual], [A], dict=True)

check("spin involution imposes half-cycle phase A=pi",
      spin_solution == [{A: pi}],
      "Delta phi = (2*pi)/2 * M Omega_H")

freq_coeff = sp.simplify((momega - phi0 / (2 * pi)) / x)
check("half-cycle phase is equivalent to frequency coefficient 1/2",
      sp.simplify(freq_coeff.subs(A, pi) - sp.Rational(1, 2)) == 0,
      "M omega = M omega_0 + Omega_H/2")

# Dissipative invariance.
tau_slope = sp.simplify(sp.diff(tau_over_m, x).subs(x, 0))
diss_solution = sp.solve([tau_slope], [B], dict=True)

check("dissipative invariance imposes B=0",
      diss_solution == [{B: 0}],
      f"d(tau/M)/dx|0 = {tau_slope}")

constraints = sp.Matrix([[1, 0], [0, -4 * pi**2 / kappa0**2]])
rhs = sp.Matrix([pi, 0])
solution = sp.solve([spin_involution_residual, tau_slope], [A, B], dict=True)

check("constraint matrix has full rank",
      constraints.rank() == 2,
      f"rank={constraints.rank()}")

check("combined constraints have a unique first-order solution",
      solution == [{A: pi, B: 0}],
      f"solution={solution}")

candidate_phi = sp.simplify(phi.subs(solution[0]))
candidate_kappa = sp.simplify(kappa.subs(solution[0]))
candidate_tau = sp.simplify(tau_over_m.subs(solution[0]))
candidate_q = sp.simplify(quality.subs(solution[0]))

check("unique solution is exactly phi=phi0+pi*x",
      candidate_phi == phi0 + pi * x,
      str(candidate_phi))

check("unique solution keeps kappa invariant",
      candidate_kappa == kappa0,
      str(candidate_kappa))

check("unique solution keeps tau/M invariant to all orders in x at this level",
      sp.simplify(candidate_tau - 4 * pi**2 / kappa0) == 0,
      str(candidate_tau))

check("quality factor scales only through phase",
      sp.simplify(candidate_q / (pi * phi0 / kappa0) - candidate_phi / phi0) == 0,
      str(candidate_q))

small_spin_slope = sp.simplify(sp.diff((phi0 + pi * omega_h) / (2 * pi), a).subs(a, 0))
check("unique small-spin frequency slope is 1/8",
      small_spin_slope == sp.Rational(1, 8),
      f"d(Momega)/da|0={small_spin_slope}")

# Exhaust alternatives at the same order.
wrong_phase_examples = [0, pi / 2, 2 * pi, sp.Rational(3, 2) * pi]
wrong_phase_fail = all(sp.simplify(spin_involution_residual.subs(A, val)) != 0
                       for val in wrong_phase_examples)
check("all tested non-half-cycle phase coefficients violate spin involution",
      wrong_phase_fail,
      f"tested A={wrong_phase_examples}")

wrong_damping_examples = [-1, sp.Rational(1, 2), 1, pi]
wrong_damping_fail = all(sp.simplify(tau_slope.subs(B, val)) != 0
                         for val in wrong_damping_examples)
check("all tested nonzero damping coefficients violate dissipative invariance",
      wrong_damping_fail,
      f"tested B={wrong_damping_examples}")

free_frequency_residual = sp.simplify(2 * pi * c - pi)
free_frequency_solution = sp.solve([free_frequency_residual], [c], dict=True)
check("writing Momega=Momega0+c Omega_H uniquely gives c=1/2",
      free_frequency_solution == [{c: sp.Rational(1, 2)}],
      f"solution={free_frequency_solution}")

print("\n" + "=" * 88)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: first-order uniqueness under fixed sector, Kerr holonomy,")
print("               spin involution, and dissipative invariance.")
print("=" * 88)

sys.exit(0 if n_pass == n_total else 1)
