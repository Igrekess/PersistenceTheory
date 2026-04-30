#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_spin_half_holonomy_PT.py
=====================================

Parameter-free PT candidate for the spin-dependent Kerr deformation of the
selected macroscopic ringdown mode.

Starting point:
  - the structural PT selector fixes the macroscopic sector (ell,m,n)=(2,2,0);
  - the undeformed T_30 Ruelle mode gives lambda=r exp(i phi_0);
  - LVK pyRing/pSEOBNR posteriors show that tau/M is already close, while
    M omega is systematically low.

PT half-holonomy rule:

    M Omega_H(a) = a / (2 (1 + sqrt(1-a^2)))
    phi(a)      = phi_0 + pi M Omega_H(a)
    kappa(a)    = kappa_0

Equivalently:

    M omega(a) = M omega_0 + Omega_H(a)/2.

The coefficient 1/2 is not fitted here: it is the spin involution/half-cycle
factor.  The deformation is a chiral phase torsion, not a dissipative-channel
opening, so kappa and tau/M remain invariant at this order.

This file is a structural test plus an optional LVK validation if the external
GWTC-4.0 remnants release has been extracted under PT_LVK_REMNANTS_DIR or
/tmp/pt_lvk_remnants.
"""

from __future__ import annotations

import io
import math
import os
import pathlib
import sys

import h5py
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

CYCLIC_PERIOD = 2.0 * math.pi
TSUN_SECONDS = 4.925490947e-6
TOL = 1e-12


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


def selected_pt_mode():
    primes = [p for p in sieve_primes(100000) if p >= 5]
    gaps = [b - a for a, b in zip(primes, primes[1:])]
    T30, _classes = build_T_gaps(gaps, 30)
    candidates = []
    for ev in np.linalg.eigvals(T30):
        if ev.imag > 1e-10 and abs(ev) < 1.0:
            phi = abs(np.angle(ev))
            kappa = -math.log(abs(ev))
            Q = CYCLIC_PERIOD * phi / (2.0 * kappa)
            candidates.append((Q, -kappa, ev, phi, kappa))
    Q, _minus_kappa, ev, phi, kappa = max(candidates)
    return {
        "lambda": ev,
        "phi": phi,
        "kappa": kappa,
        "Momega": phi / CYCLIC_PERIOD,
        "Q": Q,
        "tau_over_M": 4.0 * math.pi**2 / kappa,
    }


def horizon_omega(spin):
    a = abs(float(spin))
    return a / (2.0 * (1.0 + math.sqrt(max(0.0, 1.0 - a * a))))


def deform(mode, spin):
    omega_h = horizon_omega(spin)
    phi = mode["phi"] + math.pi * omega_h
    kappa = mode["kappa"]
    Momega = phi / CYCLIC_PERIOD
    Q = CYCLIC_PERIOD * phi / (2.0 * kappa)
    return {
        "phi": phi,
        "kappa": kappa,
        "Momega": Momega,
        "Q": Q,
        "tau_over_M": 4.0 * math.pi**2 / kappa,
    }


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


def quantiles(arr):
    return np.quantile(np.asarray(arr, dtype=float), [0.05, 0.16, 0.5, 0.84, 0.95])


def lvk_validation(mode):
    root = pathlib.Path(os.environ.get("PT_LVK_REMNANTS_DIR", "/tmp/pt_lvk_remnants"))
    summary = root / "pyring" / "plotting_scripts_and_data" / "event_summary_samples" / "events_summary_file.h5"
    if not summary.exists():
        print("\n--- LVK validation skipped ---")
        print(f"  Missing external pyRing summary at {summary}")
        return

    undeformed_domega = []
    deformed_domega = []
    deformed_dtau = []
    in90_undeformed = 0
    in90_deformed = 0
    events = 0

    with h5py.File(summary, "r") as f:
        for event in sorted(f.keys()):
            imr = f[event]["IMR"]
            af = np.asarray(imr["af"])
            mf = np.asarray(imr["Mf"])
            freq = np.asarray(imr["f_22"])
            tau_ms = np.asarray(imr["tau_22"])
            n = min(len(mf), len(freq), len(tau_ms), len(af))
            spin = float(np.median(af[:n]))
            momega_imr = np.median(2.0 * math.pi * freq[:n] * mf[:n] * TSUN_SECONDS)
            tau_over_m_imr = np.median((tau_ms[:n] / 1000.0) / (mf[:n] * TSUN_SECONDS))

            post_path = (
                root / "pyring" / "samples" / event
                / f"{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M"
                / f"pyRing_{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M.dat"
            )
            if not post_path.exists():
                continue
            samples = pd.read_csv(post_path, sep="\t", usecols=["domega_220", "dtau_220"])
            q_omega = quantiles(samples["domega_220"])
            q_tau = quantiles(samples["dtau_220"])

            d0_omega = mode["Momega"] / momega_imr - 1.0
            dm = deform(mode, spin)
            d1_omega = dm["Momega"] / momega_imr - 1.0
            d1_tau = dm["tau_over_M"] / tau_over_m_imr - 1.0

            undeformed_domega.append(d0_omega)
            deformed_domega.append(d1_omega)
            deformed_dtau.append(d1_tau)
            in90_undeformed += q_omega[0] <= d0_omega <= q_omega[4] and q_tau[0] <= d1_tau <= q_tau[4]
            in90_deformed += q_omega[0] <= d1_omega <= q_omega[4] and q_tau[0] <= d1_tau <= q_tau[4]
            events += 1

    print("\n--- Optional LVK pyRing validation ---")
    print(f"  events={events}")
    print(f"  undeformed domega q05/q16/q50/q84/q95 = {quantiles(undeformed_domega)}")
    print(f"  deformed   domega q05/q16/q50/q84/q95 = {quantiles(deformed_domega)}")
    print(f"  deformed   dtau   q05/q16/q50/q84/q95 = {quantiles(deformed_dtau)}")
    print(f"  inside pyRing 90% intervals, both coordinates:")
    print(f"    undeformed={in90_undeformed}/{events}")
    print(f"    deformed  ={in90_deformed}/{events}")

    check("LVK: half-holonomy reduces median frequency bias by >10x",
          abs(np.median(deformed_domega)) < abs(np.median(undeformed_domega)) / 10.0,
          f"median undeformed={np.median(undeformed_domega):+.3f}, deformed={np.median(deformed_domega):+.3f}")
    check("LVK: half-holonomy improves 90% event-level inclusion",
          in90_deformed > in90_undeformed,
          f"{in90_undeformed}/{events} -> {in90_deformed}/{events}")
    check("LVK: deformed tau/M remains near zero-bias",
          abs(np.median(deformed_dtau)) < 0.03,
          f"median dtau={np.median(deformed_dtau):+.3f}")


print("=" * 86)
print("  PT QG: KERR SPIN HALF-HOLONOMY DEFORMATION")
print("=" * 86)

n_pass = 0
n_total = 0

mode0 = selected_pt_mode()
print(f"\n  lambda={mode0['lambda']}")
print(f"  phi0={mode0['phi']:.12f}, kappa0={mode0['kappa']:.12f}")
print(f"  Momega0={mode0['Momega']:.6f}, Q0={mode0['Q']:.6f}, tau/M={mode0['tau_over_M']:.6f}")

check("undeformed selected mode is oscillatory and damped",
      mode0["phi"] > 0.0 and mode0["kappa"] > 0.0 and 0.0 < abs(mode0["lambda"]) < 1.0)

check("Kerr horizon angular velocity vanishes at zero spin",
      abs(horizon_omega(0.0)) < TOL)

check("Kerr horizon angular velocity tends to one half at extremality",
      abs(horizon_omega(1.0 - 1e-12) - 0.5) < 1e-6,
      f"Omega_H(1-eps)={horizon_omega(1.0 - 1e-12):.12f}")

check("zero spin recovers the undeformed PT mode",
      abs(deform(mode0, 0.0)["Momega"] - mode0["Momega"]) < TOL
      and abs(deform(mode0, 0.0)["tau_over_M"] - mode0["tau_over_M"]) < TOL)

spins = np.linspace(0.0, 0.99, 100)
deformed = [deform(mode0, a) for a in spins]
check("phase/frequency deformation is monotone in |a|",
      all(b["Momega"] >= a["Momega"] for a, b in zip(deformed, deformed[1:])))

check("damping exponent kappa is unchanged by phase holonomy",
      all(abs(d["kappa"] - mode0["kappa"]) < TOL for d in deformed))

check("tau/M is invariant under pure phase torsion",
      all(abs(d["tau_over_M"] - mode0["tau_over_M"]) < TOL for d in deformed),
      f"tau/M={mode0['tau_over_M']:.6f}")

check("quality factor scales with the deformed phase",
      all(abs(d["Q"] / mode0["Q"] - d["Momega"] / mode0["Momega"]) < TOL for d in deformed))

small = 1e-6
small_derivative = (deform(mode0, small)["Momega"] - mode0["Momega"]) / small
check("small-spin slope is the spin-involution value 1/8",
      abs(small_derivative - 0.125) < 1e-6,
      f"d(Momega)/da|0={small_derivative:.12f}")

gw150914_like = deform(mode0, 0.69)
check("GW150914-like spin lands in the observed Kerr-ringdown band",
      0.51 < gw150914_like["Momega"] < 0.54
      and 12.0 < gw150914_like["tau_over_M"] < 12.8,
      f"a=0.69: Momega={gw150914_like['Momega']:.6f}, tau/M={gw150914_like['tau_over_M']:.6f}")

lvk_validation(mode0)

print("\n" + "=" * 86)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  STATUS: derived zero-parameter Kerr spin-deformation candidate.")
print("  CLAIM: structural + LVK-supported diagnostic; not yet promoted to theorem.")
print("=" * 86)

sys.exit(0 if n_pass == n_total else 1)
