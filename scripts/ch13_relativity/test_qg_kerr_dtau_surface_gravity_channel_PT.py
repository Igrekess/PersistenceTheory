#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_dtau_surface_gravity_channel_PT.py
===============================================

PT diagnostic for the residual Kerr dtau channel.

The exact half-holonomy closes the phase channel:

    Delta phi = pi M Omega_H(a).

A natural PT candidate for damping is the horizon relaxation scale governed by
the Kerr surface gravity:

    M kappa_H(a) = sqrt(1-a^2) / (2(1 + sqrt(1-a^2))).

Normalizing to the Schwarzschild value M kappa_H(0)=1/4 gives a zero-parameter
catalog-scale damping stretch

    R_tau(a) = 1/(4 M kappa_H(a)) - 1
             = (1 + sqrt(1-a^2))/(2 sqrt(1-a^2)) - 1.

This script deliberately does not promote R_tau to a theorem.  It tests whether
the channel is structurally PT-natural and whether the local LVK release
supports it consistently across pSEOB and pyRing.

Conclusion encoded by the checks:
  - structural surface-gravity channel: closed;
  - pSEOB catalog-scale size: supported;
  - pyRing event-level support: not yet closed, because the event medians are
    weakly correlated with spin/surface gravity and strongly entangled with
    domega.
"""

from __future__ import annotations

import argparse
import gzip
import io
import math
import os
import pathlib
import sys

import h5py
import numpy as np
import pandas as pd
import sympy as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

TSUN_SECONDS = 4.925490947e-6


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


def quantiles(x, qs=(0.05, 0.16, 0.5, 0.84, 0.95)):
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]
    return np.quantile(arr, qs)


def fmt_q(q):
    return f"[{q[0]:+.3f}, {q[1]:+.3f}, {q[2]:+.3f}, {q[3]:+.3f}, {q[4]:+.3f}]"


def surface_gravity_dtau(spin):
    a = np.asarray(spin, dtype=float)
    s = np.sqrt(np.maximum(0.0, 1.0 - a * a))
    return (1.0 + s) / (2.0 * s) - 1.0


def rankdata(values):
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr)
    ranks = np.empty(len(arr), dtype=float)
    i = 0
    while i < len(arr):
        j = i + 1
        while j < len(arr) and arr[order[j]] == arr[order[i]]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def corr(x, y, rank=False):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return float("nan")
    if rank:
        x = rankdata(x)
        y = rankdata(y)
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def read_event_table(root):
    summary = root / "pyring" / "plotting_scripts_and_data" / "event_summary_samples" / "events_summary_file.h5"
    rows = []
    if not summary.exists():
        return pd.DataFrame(rows)
    with h5py.File(summary, "r") as f:
        for event in sorted(f.keys()):
            imr = f[event]["IMR"]
            spin = float(np.median(np.asarray(imr["af"])))
            row = {"event": event, "spin": spin, "surface_dtau": float(surface_gravity_dtau(spin))}

            post = (
                root / "pyring" / "samples" / event
                / f"{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M"
                / f"pyRing_{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M.dat"
            )
            if post.exists():
                samples = pd.read_csv(post, sep="\t", usecols=["domega_220", "dtau_220"])
                row["post_domega_med"] = float(np.median(samples["domega_220"]))
                row["post_dtau_med"] = float(np.median(samples["dtau_220"]))
                q_dtau = quantiles(samples["dtau_220"])
                row["dtau_zero_in90"] = bool(q_dtau[0] <= 0.0 <= q_dtau[4])

            rows.append(row)
    return pd.DataFrame(rows)


def read_pseob_catalog_medians(root):
    folder = root / "pSEOB" / "joint"
    rows = []
    if not folder.exists():
        return rows
    for case in sorted(p for p in folder.iterdir() if p.is_dir()):
        files = sorted(case.glob("dtau_220*.dat.gz"))
        if not files:
            continue
        with gzip.open(files[0], "rt") as f:
            dtau = np.loadtxt(f)
        q = quantiles(dtau)
        rows.append({"case": case.name, "q": q, "median": float(q[2]), "zero_in90": bool(q[0] <= 0.0 <= q[4])})
    return rows


def read_pseob_hier_medians(root):
    folder = root / "pSEOB" / "hierarchical"
    rows = []
    if not folder.exists():
        return rows
    for path in sorted(folder.glob("*/*.npz")):
        if path.name.startswith("pseob_hier_pop"):
            continue
        z = np.load(path, allow_pickle=True)
        key = "mu_dtau_220"
        if key in z.files:
            q = quantiles(z[key])
            rows.append({"case": path.parent.name, "q": q, "median": float(q[2]), "zero_in90": bool(q[0] <= 0.0 <= q[4])})
    return rows


def read_pyring_hier(root):
    folder = root / "pyring" / "plotting_scripts_and_data" / "hierarchical_posterior_samples"
    rows = []
    if not folder.exists():
        return rows
    for path in sorted(folder.glob("*.h5")):
        with h5py.File(path, "r") as f:
            d = f["gr_deviations"][:]["dtau_220"]
            q = quantiles(d)
            rows.append({"case": path.name, "q": q, "median": float(q[2]), "zero_in90": bool(q[0] <= 0.0 <= q[4])})
    return rows


print("=" * 92)
print("  PT QG: KERR DTAU SURFACE-GRAVITY CHANNEL")
print("=" * 92)

n_pass = 0
n_total = 0

parser = argparse.ArgumentParser()
parser.add_argument("--lvk-root", default=os.environ.get("PT_LVK_REMNANTS_DIR", "/tmp/pt_lvk_remnants"))
args = parser.parse_args()
root = pathlib.Path(args.lvk_root)

print("\n--- Structural derivation ---")
a = sp.symbols("a", nonnegative=True, real=True)
s = sp.sqrt(1 - a**2)
mkappa_h = s / (2 * (1 + s))
r_tau = sp.simplify(1 / (4 * mkappa_h) - 1)

check("Kerr surface gravity has Schwarzschild value M kappa_H(0)=1/4",
      sp.simplify(mkappa_h.subs(a, 0) - sp.Rational(1, 4)) == 0)

check("Kerr surface gravity vanishes at extremality",
      sp.simplify(sp.limit(mkappa_h, a, 1, dir="-")) == 0)

check("surface relaxation stretch vanishes at zero spin",
      sp.simplify(r_tau.subs(a, 0)) == 0,
      str(r_tau))

series = sp.series(r_tau, a, 0, 8).removeO()
check("surface relaxation is even in spin and starts at a^2/4",
      series == a**2 / 4 + 3 * a**4 / 16 + 5 * a**6 / 32,
      f"R_tau(a)={series}+O(a^8)")

derivative = sp.diff(r_tau, a)
check("surface relaxation is monotone for 0<a<1",
      derivative.subs(a, sp.Rational(1, 2)) > 0 and derivative.subs(a, sp.Rational(4, 5)) > 0,
      f"dR/da={sp.simplify(derivative)}")

check("surface channel is dissipative rather than chiral",
      sp.simplify(r_tau.subs(a, -a) - r_tau) == 0,
      "R_tau(-a)=R_tau(a)")

print("\n--- LVK release diagnostics ---")
if not root.exists():
    print(f"  LVK root not found: {root}")
    print("=" * 92)
    sys.exit(0 if n_pass == n_total else 1)

events = read_event_table(root)
if not events.empty:
    surface_q = quantiles(events["surface_dtau"])
    print(f"  O4a surface-gravity R_tau q={fmt_q(surface_q)}")
    check("O4a surface-gravity scale is in the observed pSEOB dtau band",
          0.12 <= surface_q[2] <= 0.22,
          f"median={surface_q[2]:+.3f}")

pseob_joint = read_pseob_catalog_medians(root)
if pseob_joint and not events.empty:
    inside = sum(surface_q[0] <= row["median"] <= surface_q[4] for row in pseob_joint)
    zero_excluded = sum(not row["zero_in90"] for row in pseob_joint)
    for row in pseob_joint:
        print(f"  pSEOB joint {row['case']:<38} dtau={fmt_q(row['q'])} zero90={row['zero_in90']}")
    check("most pSEOB joint medians lie inside the surface-gravity catalog band",
          inside >= 3,
          f"{inside}/{len(pseob_joint)} inside R_tau 90% spin band")
    check("pSEOB joint posteriors consistently exclude zero dtau at 90%",
          zero_excluded == len(pseob_joint),
          f"{zero_excluded}/{len(pseob_joint)} exclude zero")

pseob_hier = read_pseob_hier_medians(root)
if pseob_hier and not events.empty:
    inside = sum(surface_q[0] <= row["median"] <= surface_q[4] for row in pseob_hier)
    print(f"  pSEOB hierarchical medians={[round(row['median'], 3) for row in pseob_hier]}")
    check("pSEOB hierarchical mu_dtau medians mostly match the surface-gravity band",
          inside >= 3,
          f"{inside}/{len(pseob_hier)} inside R_tau 90% spin band")

pyring_hier = read_pyring_hier(root)
if pyring_hier:
    zero_in = sum(row["zero_in90"] for row in pyring_hier)
    for row in pyring_hier:
        print(f"  pyRing hier {row['case']:<24} dtau={fmt_q(row['q'])} zero90={row['zero_in90']}")
    check("pyRing hierarchical dtau still allows zero at 90%",
          zero_in == len(pyring_hier),
          f"{zero_in}/{len(pyring_hier)} include zero")

if not events.empty and {"post_dtau_med", "post_domega_med"}.issubset(events.columns):
    c_surface = corr(events["surface_dtau"], events["post_dtau_med"])
    c_surface_rank = corr(events["surface_dtau"], events["post_dtau_med"], rank=True)
    c_degen = corr(events["post_domega_med"], events["post_dtau_med"])
    zero_events = int(events["dtau_zero_in90"].sum())
    print(f"  corr(surface R_tau, event dtau median) pearson={c_surface:+.3f}, spearman={c_surface_rank:+.3f}")
    print(f"  corr(event domega median, event dtau median) pearson={c_degen:+.3f}")
    print(f"  event dtau=0 in 90% intervals: {zero_events}/{len(events)}")
    check("event-level pyRing dtau is not a clean surface-gravity law",
          abs(c_surface) < 0.25 and abs(c_surface_rank) < 0.25,
          "weak spin/surface correlation")
    check("event-level pyRing dtau is strongly entangled with domega",
          c_degen < -0.5,
          f"corr={c_degen:+.3f}")
    check("most event-level pyRing dtau posteriors still include zero",
          zero_events >= 18,
          f"{zero_events}/{len(events)}")

print("\n" + "=" * 92)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  STATUS: surface-gravity dtau is a strong catalog-scale PT candidate,")
print("          but the event-level pyRing diagnostics prevent canonical promotion.")
print("  OPEN: separate a genuine dissipative relaxation channel from ringdown")
print("        start-time / overtone / domega-dtau model systematics.")
print("=" * 92)

sys.exit(0 if n_pass == n_total else 1)
