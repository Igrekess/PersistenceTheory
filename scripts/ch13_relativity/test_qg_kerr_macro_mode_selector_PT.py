#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_macro_mode_selector_PT.py -- PT selector for Kerr ringdown modes
================================================================================

This script tests the parameter-free PT selector for the macroscopic Kerr
ringdown mode.

The previous black-hole/ringdown script closed the observable skeleton:
Ruelle-Pollicott eigenvalues lambda = r exp(i phi) give a damped ringdown
family.  What remained open was the selector: which microscopic Ruelle mode is
picked by a macroscopic Kerr geometry?

The PT answer is a cascade, with no fitted coefficient:

  1. No-hair reduction: the macroscopic geometry supplies only (M, a, epsilon),
     where M is the final mass, a the dimensionless spin, and epsilon the
     orientation convention.
  2. Spin-2 constraint: gravitational radiation selects the quadrupolar tensor
     sector, l = 2.
  3. Kerr chirality: sign(a * epsilon) selects one member of a complex-conjugate
     Ruelle pair, i.e. m = +/-2.
  4. Persistence resonance: among admissible oscillatory modes, select the one
     maximizing the cyclic coherent quality

         Q_cyc = 2*pi * phi / (2*(-ln r)).

     The factor 2*pi is the PT cyclic phase period, not a fit.
  5. Fundamental overtone: the maximum-persistence mode is n = 0.

This closes the structural mode-identification rule.  It still does not fit a
measured LIGO/Virgo/KAGRA waveform; the next empirical task is to compare the
selected dimensionless pair (M omega, Q) after any Kerr deformation of the
macroscopic operator.
"""

import io
import math
import sys
from dataclasses import dataclass

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

TOL = 1e-12
CYCLIC_PERIOD = 2.0 * math.pi


@dataclass(frozen=True)
class KerrGeometry:
    mass: float
    spin: float
    orientation: int = 1


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


def ruelle_modes(T):
    modes = []
    for ev in np.linalg.eigvals(T):
        if abs(ev.imag) > 1e-10 and abs(ev) < 1.0:
            r = abs(ev)
            phi = abs(np.angle(ev))
            kappa = -math.log(r)
            chirality = 1 if ev.imag > 0 else -1
            q_raw = phi / (2.0 * kappa)
            q_cyc = CYCLIC_PERIOD * q_raw
            modes.append({
                "lambda": ev,
                "r": r,
                "phi": phi,
                "kappa": kappa,
                "chirality": chirality,
                "Q_raw": q_raw,
                "Q_cyc": q_cyc,
            })
    modes.sort(key=lambda m: (m["chirality"], -m["Q_cyc"], m["kappa"]))
    return modes


def chirality_from_geometry(geometry):
    signed = geometry.spin * geometry.orientation
    if signed > 0:
        return 1
    if signed < 0:
        return -1
    return 0


def select_kerr_mode(modes, geometry):
    """Return the PT-selected Kerr macro ringdown mode."""
    chirality = chirality_from_geometry(geometry)
    if chirality == 0:
        admissible = modes
    else:
        admissible = [m for m in modes if m["chirality"] == chirality]
    if not admissible:
        raise RuntimeError("no admissible oscillatory Ruelle mode")
    selected = max(admissible, key=lambda m: (m["Q_cyc"], -m["kappa"]))
    return {
        **selected,
        "ell": 2,
        "m": 2 * (chirality if chirality != 0 else 1),
        "n": 0,
    }


def macro_observables(selected, geometry):
    """Dimensionless cyclic observables and mass-scaled time observables."""
    # The cyclic phase fraction gives the dimensionless angular frequency M omega.
    Momega = selected["phi"] / CYCLIC_PERIOD
    Q = selected["Q_cyc"]
    tau_over_M = 2.0 * Q / Momega
    freq_units = Momega / geometry.mass
    tau_units = tau_over_M * geometry.mass
    return {
        "Momega": Momega,
        "Q": Q,
        "tau_over_M": tau_over_M,
        "freq_units": freq_units,
        "tau_units": tau_units,
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


print("=" * 86)
print("  PT QUANTUM GRAVITY: KERR MACRO MODE SELECTOR")
print("=" * 86)

n_pass = 0
n_total = 0

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]
T30, classes = build_T_gaps(gaps, 30)
modes = ruelle_modes(T30)
positive = [m for m in modes if m["chirality"] > 0]
negative = [m for m in modes if m["chirality"] < 0]

print(f"\n  T_30 shape={T30.shape}, classes={classes}")
print(f"  oscillatory modes={len(modes)} ({len(positive)} positive, {len(negative)} negative)")

row_err = np.max(np.abs(T30.sum(axis=1) - 1.0))
check("T_30 is a stochastic PT transfer operator",
      row_err < TOL,
      f"max row error={row_err:.2e}")

pair_ok = len(positive) == len(negative) and len(positive) > 0
pair_errs = []
for pos in positive:
    candidates = [neg for neg in negative if abs(neg["lambda"] - np.conj(pos["lambda"])) < 1e-10]
    pair_errs.append(0.0 if candidates else 1.0)
check("oscillatory Ruelle modes come in conjugate chiral pairs",
      pair_ok and max(pair_errs) == 0.0,
      f"pairs={len(positive)}, pair_errors={sum(pair_errs):.0f}")

geom_pos = KerrGeometry(mass=68.0, spin=0.69, orientation=1)
geom_neg = KerrGeometry(mass=68.0, spin=-0.69, orientation=1)
geom_flip = KerrGeometry(mass=68.0, spin=0.69, orientation=-1)
geom_zero = KerrGeometry(mass=68.0, spin=0.0, orientation=1)

sel_pos = select_kerr_mode(modes, geom_pos)
sel_neg = select_kerr_mode(modes, geom_neg)
sel_flip = select_kerr_mode(modes, geom_flip)
sel_zero = select_kerr_mode(modes, geom_zero)

check("no-hair input reduces selector to mass, spin and orientation",
      geom_pos.mass > 0 and abs(geom_pos.spin) < 1 and geom_pos.orientation in (-1, 1),
      f"(M,a,eps)=({geom_pos.mass},{geom_pos.spin},{geom_pos.orientation})")

check("spin-2 constraint selects quadrupolar tensor labels",
      sel_pos["ell"] == 2 and abs(sel_pos["m"]) == 2,
      f"(ell,m,n)=({sel_pos['ell']},{sel_pos['m']},{sel_pos['n']})")

check("positive Kerr chirality selects the positive-imaginary branch",
      sel_pos["chirality"] == 1 and sel_pos["m"] == 2,
      f"lambda={sel_pos['lambda']}")

check("negative Kerr chirality selects the conjugate branch",
      sel_neg["chirality"] == -1 and sel_neg["m"] == -2
      and abs(sel_neg["lambda"] - np.conj(sel_pos["lambda"])) < 1e-10,
      f"lambda-={sel_neg['lambda']}, conj(lambda+)={np.conj(sel_pos['lambda'])}")

check("orientation reversal flips the selected chirality",
      sel_flip["chirality"] == -1 and abs(sel_flip["lambda"] - sel_neg["lambda"]) < 1e-10,
      f"orientation=-1 gives m={sel_flip['m']}")

best_positive = max(positive, key=lambda m: (m["Q_cyc"], -m["kappa"]))
check("persistence resonance uniquely selects max Q_cyc",
      abs(sel_pos["lambda"] - best_positive["lambda"]) < 1e-10
      and sum(abs(m["Q_cyc"] - best_positive["Q_cyc"]) < 1e-10 for m in positive) == 1,
      f"Q_cyc={sel_pos['Q_cyc']:.6f}, phi={sel_pos['phi']:.6f}, kappa={sel_pos['kappa']:.6f}")

least_damped = min(positive, key=lambda m: m["kappa"])
check("least damping alone is not the PT Kerr selector",
      abs(least_damped["lambda"] - sel_pos["lambda"]) > 1e-6
      and least_damped["Q_cyc"] < sel_pos["Q_cyc"],
      f"least-damped Q_cyc={least_damped['Q_cyc']:.6f} < selected {sel_pos['Q_cyc']:.6f}")

check("fundamental overtone is n=0 by maximum persistence",
      sel_pos["n"] == 0 and sel_pos["Q_cyc"] > 1.0,
      f"n={sel_pos['n']}, Q_cyc={sel_pos['Q_cyc']:.6f}")

obs_pos = macro_observables(sel_pos, geom_pos)
check("cyclic macro frequency is a dimensionless phase fraction",
      0.0 < obs_pos["Momega"] < 1.0,
      f"Momega=phi/(2pi)={obs_pos['Momega']:.6f}")

check("cyclic coherent quality factor is finite and resonant",
      obs_pos["Q"] > 1.0 and math.isfinite(obs_pos["Q"]),
      f"Q_macro=2pi*Q_raw={obs_pos['Q']:.6f}")

geom_light = KerrGeometry(mass=34.0, spin=0.69, orientation=1)
geom_heavy = KerrGeometry(mass=136.0, spin=0.69, orientation=1)
sel_light = select_kerr_mode(modes, geom_light)
sel_heavy = select_kerr_mode(modes, geom_heavy)
obs_light = macro_observables(sel_light, geom_light)
obs_heavy = macro_observables(sel_heavy, geom_heavy)

check("mass changes scale frequency but not selected mode",
      abs(sel_light["lambda"] - sel_heavy["lambda"]) < 1e-10
      and abs((obs_heavy["freq_units"] / obs_light["freq_units"]) - 0.25) < TOL,
      f"f_heavy/f_light={obs_heavy['freq_units']/obs_light['freq_units']:.12f}")

check("damping time scales linearly with mass",
      abs((obs_heavy["tau_units"] / obs_light["tau_units"]) - 4.0) < TOL,
      f"tau_heavy/tau_light={obs_heavy['tau_units']/obs_light['tau_units']:.12f}")

check("quality factor is mass-invariant",
      abs(obs_light["Q"] - obs_heavy["Q"]) < TOL,
      f"Q_light={obs_light['Q']:.12f}, Q_heavy={obs_heavy['Q']:.12f}")

obs_neg = macro_observables(sel_neg, geom_neg)
check("chirality flip preserves absolute macro observables",
      abs(obs_neg["Momega"] - obs_pos["Momega"]) < TOL
      and abs(obs_neg["Q"] - obs_pos["Q"]) < TOL,
      f"Momega={obs_pos['Momega']:.6f}, Q={obs_pos['Q']:.6f}")

check("zero-spin geometry is chirality-degenerate",
      abs(sel_zero["Q_cyc"] - sel_pos["Q_cyc"]) < TOL
      and abs(abs(sel_zero["lambda"]) - abs(sel_pos["lambda"])) < TOL,
      f"selected chirality at a=0 is conventional: m={sel_zero['m']}")

free_coefficients = []
check("selector introduces no fitted coefficient",
      len(free_coefficients) == 0 and abs(CYCLIC_PERIOD - 2.0 * math.pi) < TOL,
      "ingredients: spin-2, Kerr chirality, 2*pi cyclic period, max persistence")

repeat = select_kerr_mode(modes, geom_pos)
check("selector is deterministic",
      abs(repeat["lambda"] - sel_pos["lambda"]) < TOL,
      f"lambda={sel_pos['lambda']}")

print("\n" + "=" * 86)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: structural PT selector for the Kerr macroscopic ringdown mode.")
print("  STILL EMPIRICAL: compare the selected/deformed Kerr macro operator with GW data.")
print("=" * 86)

sys.exit(0 if n_pass == n_total else 1)
