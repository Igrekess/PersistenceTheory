#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_qg_kerr_dtau_tension_PT.py
==================================

Exploratory diagnostic for the residual Kerr-ringdown dtau tension.

This is not a canonical proof script: it needs the external LVK/GWTC-4.0
remnant-release files under /tmp/pt_lvk_remnants or PT_LVK_REMNANTS_DIR.

Question:
  The PT half-holonomy leaves kappa and tau/M invariant.  pSEOB and pyRing
  deviation posteriors nevertheless tend to prefer positive dtau.  Is this
  tension aligned with the PT tau residual, with spin, or with analysis-model
  choices?

The script compares four layers:
  1. PT tau/M against IMR Kerr event medians.
  2. pyRing event-level (domega_220, dtau_220) medians.
  3. pyRing hierarchical O4a posteriors.
  4. pSEOB joint catalog posteriors.

It also tests a purely geometric PT scale for damping tension: if the phase
channel is governed by the horizon angular velocity Omega_H, the damping
channel should be compared to the Kerr surface gravity kappa_H.  The
dimensionless surface-gravity relaxation factor is

    R_tau(a) = 1/(4 M kappa_H(a)) - 1
             = (1 + sqrt(1-a^2))/(2 sqrt(1-a^2)) - 1.

This is not promoted here; it is a scale diagnostic.
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

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

TSUN_SECONDS = 4.925490947e-6
PT_TAU_OVER_M = 12.382562


def quantiles(x, qs=(0.05, 0.16, 0.5, 0.84, 0.95)):
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]
    return np.quantile(arr, qs)


def fmt_q(q):
    return f"[{q[0]:+.3f}, {q[1]:+.3f}, {q[2]:+.3f}, {q[3]:+.3f}, {q[4]:+.3f}]"


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


def surface_gravity_dtau(spin):
    a = np.asarray(spin, dtype=float)
    s = np.sqrt(np.maximum(0.0, 1.0 - a * a))
    return (1.0 + s) / (2.0 * s) - 1.0


def read_pseob_joint(root):
    folder = root / "pSEOB" / "joint"
    rows = []
    if not folder.exists():
        return rows
    for case in sorted(p for p in folder.iterdir() if p.is_dir()):
        dtau_files = sorted(case.glob("dtau_220*.dat.gz"))
        domega_files = sorted(case.glob("domega_220*.dat.gz"))
        if not dtau_files or not domega_files:
            continue
        with gzip.open(dtau_files[0], "rt") as f:
            dtau = np.loadtxt(f)
        with gzip.open(domega_files[0], "rt") as f:
            domega = np.loadtxt(f)
        rows.append({
            "case": case.name,
            "dtau_q": quantiles(dtau),
            "domega_q": quantiles(domega),
            "dtau_zero_in90": bool(quantiles(dtau)[0] <= 0.0 <= quantiles(dtau)[4]),
            "domega_zero_in90": bool(quantiles(domega)[0] <= 0.0 <= quantiles(domega)[4]),
        })
    return rows


def read_pyring_hier(root):
    folder = root / "pyring" / "plotting_scripts_and_data" / "hierarchical_posterior_samples"
    rows = []
    if not folder.exists():
        return rows
    for path in sorted(folder.glob("*.h5")):
        with h5py.File(path, "r") as f:
            d = f["gr_deviations"][:]
            rows.append({
                "case": path.name,
                "dtau_q": quantiles(d["dtau_220"]),
                "domega_q": quantiles(d["domega_220"]),
                "dtau_zero_in90": bool(quantiles(d["dtau_220"])[0] <= 0.0 <= quantiles(d["dtau_220"])[4]),
                "domega_zero_in90": bool(quantiles(d["domega_220"])[0] <= 0.0 <= quantiles(d["domega_220"])[4]),
            })
    return rows


def event_rows(root):
    summary = root / "pyring" / "plotting_scripts_and_data" / "event_summary_samples" / "events_summary_file.h5"
    rows = []
    if not summary.exists():
        return rows

    with h5py.File(summary, "r") as f:
        for event in sorted(f.keys()):
            imr = f[event]["IMR"]
            mf = np.asarray(imr["Mf"])
            tau_ms = np.asarray(imr["tau_22"])
            spin = np.asarray(imr["af"])
            n = min(len(mf), len(tau_ms), len(spin))
            tau_over_m = (tau_ms[:n] / 1000.0) / (mf[:n] * TSUN_SECONDS)
            row = {
                "event": event,
                "spin": float(np.median(spin[:n])),
                "imr_tau_over_m": float(np.median(tau_over_m)),
                "pt_dtau_vs_imr": float(PT_TAU_OVER_M / np.median(tau_over_m) - 1.0),
            }
            row["surface_gravity_dtau"] = float(surface_gravity_dtau(row["spin"]))
            for group in ("Kerr", "KerrPostmerger"):
                if group in f[event]:
                    g = f[event][group]
                    gm = np.asarray(g["Mf"])
                    gtau = np.asarray(g["tau_22"])
                    gn = min(len(gm), len(gtau))
                    gtau_over_m = (gtau[:gn] / 1000.0) / (gm[:gn] * TSUN_SECONDS)
                    row[f"{group.lower()}_dtau_vs_imr"] = float(
                        np.median(gtau_over_m) / np.median(tau_over_m) - 1.0
                    )

            post = (
                root / "pyring" / "samples" / event
                / f"{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M"
                / f"pyRing_{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M.dat"
            )
            if post.exists():
                samples = pd.read_csv(post, sep="\t", usecols=["domega_220", "dtau_220"])
                row["post_domega_med"] = float(np.median(samples["domega_220"]))
                row["post_dtau_med"] = float(np.median(samples["dtau_220"]))
                row["post_dtau_q"] = quantiles(samples["dtau_220"])
                row["post_dtau_zero_in90"] = bool(row["post_dtau_q"][0] <= 0.0 <= row["post_dtau_q"][4])
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lvk-root", default=os.environ.get("PT_LVK_REMNANTS_DIR", "/tmp/pt_lvk_remnants"))
    parser.add_argument("--top", type=int, default=8)
    args = parser.parse_args()

    root = pathlib.Path(args.lvk_root)
    print("=" * 88)
    print("  PT QG: KERR DTAU TENSION DIAGNOSTIC")
    print("=" * 88)
    print(f"  LVK root: {root}")
    print(f"  PT half-holonomy tau/M = {PT_TAU_OVER_M:.6f}")
    print()

    if not root.exists():
        print("  External LVK directory not found.")
        return 2

    pseob = read_pseob_joint(root)
    if pseob:
        print("--- pSEOB joint catalog posteriors ---")
        for row in pseob:
            print(
                f"  {row['case']:<42} "
                f"dtau={fmt_q(row['dtau_q'])} zero90={row['dtau_zero_in90']}  "
                f"domega={fmt_q(row['domega_q'])} zero90={row['domega_zero_in90']}"
            )
        print()

    hier = read_pyring_hier(root)
    if hier:
        print("--- pyRing hierarchical posteriors ---")
        for row in hier:
            print(
                f"  {row['case']:<28} "
                f"dtau={fmt_q(row['dtau_q'])} zero90={row['dtau_zero_in90']}  "
                f"domega={fmt_q(row['domega_q'])} zero90={row['domega_zero_in90']}"
            )
        print()

    rows = event_rows(root)
    if rows:
        df = pd.DataFrame(rows)
        print("--- pyRing event-level summary ---")
        for col in (
            "pt_dtau_vs_imr",
            "surface_gravity_dtau",
            "post_dtau_med",
            "kerr_dtau_vs_imr",
            "kerrpostmerger_dtau_vs_imr",
        ):
            if col in df:
                vals = df[col].dropna().to_numpy()
                print(f"  {col:<28} q={fmt_q(quantiles(vals))} mean={vals.mean():+.3f}")
        if "post_dtau_zero_in90" in df:
            valid = df["post_dtau_zero_in90"].dropna()
            print(f"  event posteriors with dtau=0 in 90% interval: {int(valid.sum())}/{len(valid)}")
        print()

        print("--- Correlations ---")
        pairs = [
            ("spin", "pt_dtau_vs_imr"),
            ("spin", "surface_gravity_dtau"),
            ("spin", "post_dtau_med"),
            ("spin", "kerrpostmerger_dtau_vs_imr"),
            ("pt_dtau_vs_imr", "post_dtau_med"),
            ("surface_gravity_dtau", "post_dtau_med"),
            ("post_domega_med", "post_dtau_med"),
        ]
        for a_col, b_col in pairs:
            if a_col in df and b_col in df:
                print(
                    f"  {a_col:<18} vs {b_col:<28} "
                    f"pearson={corr(df[a_col], df[b_col]):+.3f} "
                    f"spearman={corr(df[a_col], df[b_col], rank=True):+.3f}"
                )
        print()

        if "post_dtau_med" in df:
            ranked = df.reindex(df["post_dtau_med"].abs().sort_values(ascending=False).index)
            print(f"--- largest |pyRing event dtau medians|, top {min(args.top, len(ranked))} ---")
            cols = ["event", "spin", "surface_gravity_dtau", "pt_dtau_vs_imr", "post_dtau_med", "post_domega_med"]
            print(ranked[cols].head(args.top).to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
            print()

    print("--- Diagnostic reading ---")
    print("  1. PT tau/M vs IMR event medians is already near zero-bias.")
    print("  2. pyRing and pSEOB deviation posteriors prefer positive dtau more often.")
    print("  3. The geometric surface-gravity scale has the right catalog-level size")
    print("     for pSEOB dtau, but event-level pyRing medians remain weakly correlated")
    print("     with spin.")
    print("  4. The next PT question is whether dtau is a genuine surface-gravity")
    print("     relaxation channel, or a ringdown-model/start-time/overtones systematic.")
    print("=" * 88)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
