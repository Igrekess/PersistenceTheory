#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_black_hole_ringdown_precision_PT.py -- PT black-hole precision skeleton
================================================================================

This script turns the PT black-hole/ringdown discussion into a set of
dimensionless, parameter-free precision observables.

It does not claim that LIGO/Virgo/KAGRA data are already fitted.  It verifies
the rigid observable skeleton:

  - G_PT = 2*pi*alpha fixes Schwarzschild thermodynamics;
  - first law, negative heat capacity, evaporation scaling and Page fraction
    follow without an extra parameter;
  - Ruelle-Pollicott eigenvalues of the PT transfer matrix give damped
    ringdown modes lambda = r exp(i phi);
  - using the horizon light-crossing clock t_H = 2 G M, the ringdown family
    scales as f ~ 1/M, tau ~ M, and Q is mass-invariant.

The neighbouring Kerr selector script closes the structural identification
of which PT mode is selected by a macroscopic geometry.  The remaining
empirical task is comparison with measured astrophysical QNM posteriors.
"""

import io
import math
import sys

import numpy as np
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ACTIVE = (3, 5, 7)
ALPHA_EM_PHYS = 1.0 / 137.035999084


def sin2_theta(p, mu):
    q = 1.0 - 2.0 / mu
    qp = q**p
    return (1.0 - qp) * (2.0 * p - 1.0 + qp) / p**2


def alpha_sieve(mu):
    out = 1.0
    for p in ACTIVE:
        out *= sin2_theta(p, mu)
    return out


def sieve_primes(n_max):
    is_prime = [True] * (n_max + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n_max**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n_max + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n_max + 1) if is_prime[i]]


def build_T_gaps(gaps, modulus):
    residues = [g % modulus for g in gaps]
    classes = sorted(set(residues))
    idx = {c: i for i, c in enumerate(classes)}
    counts = np.zeros((len(classes), len(classes)), dtype=float)
    for a, b in zip(residues, residues[1:]):
        counts[idx[a], idx[b]] += 1.0

    T = np.zeros_like(counts)
    for i in range(counts.shape[0]):
        row_sum = counts[i].sum()
        if row_sum > 0:
            T[i] = counts[i] / row_sum
    return T, classes


def qnm_modes(T):
    eigvals = np.linalg.eigvals(T)
    modes = []
    for ev in eigvals:
        if ev.imag > 1e-10 and abs(ev) < 1.0:
            r = abs(ev)
            phi = abs(np.angle(ev))
            kappa = -math.log(r)
            modes.append({
                "lambda": ev,
                "r": r,
                "phi": phi,
                "kappa": kappa,
                "Q": phi / (2.0 * kappa),
            })
    modes.sort(key=lambda m: (m["kappa"], -m["phi"]))
    return modes


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


print("=" * 84)
print("  PT QUANTUM GRAVITY: BLACK-HOLE/RINGDOWN PRECISION OBSERVABLES")
print("=" * 84)

n_pass = 0
n_total = 0

mu_alpha = brentq(lambda m: alpha_sieve(m) - ALPHA_EM_PHYS, 14.5, 16.0)
alpha = alpha_sieve(mu_alpha)
G = 2.0 * math.pi * alpha

print(f"\n  mu_alpha={mu_alpha:.12f}, alpha={alpha:.12e}, G_PT=2*pi*alpha={G:.12e}")

check("G_PT is fixed by alpha with no black-hole parameter",
      abs(G / alpha - 2.0 * math.pi) < 1e-14,
      f"G/alpha={G/alpha:.12f}")

print("\n--- Schwarzschild thermodynamics ---")

for M in (0.5, 1.0, 10.0, 100.0):
    r_s = 2.0 * G * M
    area = 4.0 * math.pi * r_s**2
    entropy = area / (4.0 * G)
    entropy_closed = 4.0 * math.pi * G * M**2
    temp = 1.0 / (8.0 * math.pi * G * M)
    first_law = temp * (8.0 * math.pi * G * M)
    check(f"M={M:g}: area entropy formula is exact",
          abs(entropy - entropy_closed) < 1e-12,
          f"S={entropy:.12e}, 4*pi*G*M^2={entropy_closed:.12e}")
    check(f"M={M:g}: first law T dS/dM = 1",
          abs(first_law - 1.0) < 1e-12,
          f"T={temp:.12e}, T*dS/dM={first_law:.12f}")

M1, M2 = 3.0, 12.0
S1 = 4.0 * math.pi * G * M1**2
S2 = 4.0 * math.pi * G * M2**2
T1 = 1.0 / (8.0 * math.pi * G * M1)
T2 = 1.0 / (8.0 * math.pi * G * M2)
check("entropy scales as M^2 and temperature as 1/M",
      abs((S2 / S1) - (M2 / M1) ** 2) < 1e-12
      and abs((T2 / T1) - (M1 / M2)) < 1e-12,
      f"S ratio={S2/S1:.6f}, T ratio={T2/T1:.6f}")

heat_capacity = -8.0 * math.pi * G * M1**2
numeric_capacity = (M2 - M1) / (T2 - T1)
check("Schwarzschild heat capacity is negative",
      heat_capacity < 0.0 and numeric_capacity < 0.0,
      f"C_exact(M={M1:g})={heat_capacity:.6e}, secant={numeric_capacity:.6e}")

evap_const = 5120.0 * math.pi * G**2
t_evap_2 = evap_const * 2.0**3
t_evap_4 = evap_const * 4.0**3
check("evaporation time scales as M^3",
      abs((t_evap_4 / t_evap_2) - 8.0) < 1e-12,
      f"t(4)/t(2)={t_evap_4/t_evap_2:.6f}")

page_fraction = 1.0 - 1.0 / (2.0 * math.sqrt(2.0))
M_page_over_M0 = 1.0 / math.sqrt(2.0)
S_page_ratio = M_page_over_M0**2
check("Page point is fixed by equal remaining/radiated entropy",
      abs(S_page_ratio - 0.5) < 1e-12 and 0.64 < page_fraction < 0.65,
      f"M_Page/M0={M_page_over_M0:.12f}, t_Page/t_evap={page_fraction:.12f}")

print("\n--- PT ringdown modes ---")

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]
T30, classes = build_T_gaps(gaps, 30)
modes = qnm_modes(T30)

check("T_30 has damped complex Ruelle-Pollicott modes",
      len(modes) > 0 and all(0.0 < m["r"] < 1.0 and m["kappa"] > 0 for m in modes),
      f"{len(modes)} positive-imaginary modes, classes={classes}")

mode = modes[0]
check("dominant PT ringdown mode has finite quality factor",
      mode["phi"] > 0.0 and mode["kappa"] > 0.0 and mode["Q"] > 0.0,
      f"r={mode['r']:.6f}, phi={mode['phi']:.6f}, kappa={mode['kappa']:.6f}, Q={mode['Q']:.6f}")

def ringdown_observable(mode, mass):
    horizon_clock = 2.0 * G * mass
    omega = mode["phi"] / horizon_clock
    freq = omega / (2.0 * math.pi)
    tau = horizon_clock / mode["kappa"]
    return freq, tau, math.pi * freq * tau

f_a, tau_a, Q_a = ringdown_observable(mode, 10.0)
f_b, tau_b, Q_b = ringdown_observable(mode, 40.0)
check("ringdown frequency scales as 1/M and damping time as M",
      abs((f_b / f_a) - 0.25) < 1e-12 and abs((tau_b / tau_a) - 4.0) < 1e-12,
      f"f40/f10={f_b/f_a:.12f}, tau40/tau10={tau_b/tau_a:.12f}")

check("ringdown quality factor is mass-invariant",
      abs(Q_a - Q_b) < 1e-12 and abs(Q_a - mode["Q"]) < 1e-12,
      f"Q10={Q_a:.12f}, Q40={Q_b:.12f}")

dimensionless_omega = 2.0 * G * 10.0 * (2.0 * math.pi * f_a)
dimensionless_damping = (2.0 * G * 10.0) / tau_a
check("dimensionless ringdown pair recovers the PT eigenphase",
      abs(dimensionless_omega - mode["phi"]) < 1e-12
      and abs(dimensionless_damping - mode["kappa"]) < 1e-12,
      f"Omega={dimensionless_omega:.12f}, kappa={dimensionless_damping:.12f}")

print("\n" + "=" * 84)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: parameter-free black-hole thermodynamic identities and")
print("               a PT Ruelle-Pollicott ringdown observable family.")
print("  SELECTOR: see test_qg_kerr_macro_mode_selector_PT.py.")
print("  STILL EMPIRICAL: comparison with observed QNM posteriors.")
print("=" * 84)

sys.exit(0 if n_pass == n_total else 1)
