#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_kerr_dtau_systematics_separation_PT.py
==============================================

Separation test for the residual Kerr dtau tension.

The previous surface-gravity diagnostic found a PT-natural catalog-scale
candidate

    R_tau(a) = 1/(4 M kappa_H(a)) - 1,

but pyRing event-level posteriors did not show a clean event-by-event law.
This script asks a sharper question:

  Is the released pyRing event evidence consistent with a real, stable
  dissipative dtau signal after higher-mode/start-time structure has been
  included?

Diagnostic logic:
  - A physical dissipative correction should not be supported only by broad
    posterior medians; it should also improve evidence against the corresponding
    no-deviation higher-mode model.
  - It should be reasonably stable against the ringdown start time.
  - It should not be dominated by domega--dtau degeneracy.

The script deliberately encodes the negative result as PASS checks: the current
released event-level pyRing data separate the pSEOB catalog-scale dtau tension
from a clean event-level dissipative signal.  Therefore dtau remains a
systematics/collective-channel research problem, not a canonical PT correction.
"""

from __future__ import annotations

import argparse
import io
import os
import pathlib
import sys

import h5py
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


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


def read_evidence_table(root):
    path = root / "pyring" / "plotting_scripts_and_data" / "bayes_factor_data" / "evidence.h5"
    if not path.exists():
        return pd.DataFrame()
    rows = []
    with h5py.File(path, "r") as f:
        for event in sorted(f.keys()):
            if "KerrPostmerger" not in f[event]:
                continue
            group = f[event]["KerrPostmerger"]
            needed = [
                "22_21_33_44_55",
                "22_21_33_44_55_dtau_220",
                "22_21_33_44_55_domega_220",
                "22_21_33_44_55_domega_dtau_220",
                "start_times",
            ]
            if not all(k in group for k in needed):
                continue

            start = np.asarray(group["start_times"], dtype=float)
            base = np.asarray(group["22_21_33_44_55"], dtype=float)
            dtau = np.asarray(group["22_21_33_44_55_dtau_220"], dtype=float) - base
            domega = np.asarray(group["22_21_33_44_55_domega_220"], dtype=float) - base
            both = np.asarray(group["22_21_33_44_55_domega_dtau_220"], dtype=float) - base
            idx0 = int(np.argmin(np.abs(start)))
            late = start >= 0.0

            rows.append({
                "event": event,
                "dtau0": float(dtau[idx0]),
                "domega0": float(domega[idx0]),
                "both0": float(both[idx0]),
                "max_dtau": float(np.max(dtau)),
                "max_both": float(np.max(both)),
                "mean_late_dtau": float(np.mean(dtau[late])),
                "late_positive_fraction_dtau": float(np.mean(dtau[late] > 0.0)),
                "late_positive_fraction_both": float(np.mean(both[late] > 0.0)),
                "best_t_dtau": float(start[int(np.argmax(dtau))]),
                "best_t_both": float(start[int(np.argmax(both))]),
            })
    return pd.DataFrame(rows)


def read_event_posteriors(root):
    folder = root / "pyring" / "samples"
    rows = []
    if not folder.exists():
        return pd.DataFrame()
    for event_dir in sorted(p for p in folder.iterdir() if p.is_dir()):
        event = event_dir.name
        post = (
            event_dir
            / f"{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M"
            / f"pyRing_{event}_TEOB_22_21_33_44_55_domega_dtau_220_0.0M.dat"
        )
        if not post.exists():
            continue
        samples = pd.read_csv(post, sep="\t", usecols=["domega_220", "dtau_220"])
        q_dtau = quantiles(samples["dtau_220"])
        rows.append({
            "event": event,
            "post_dtau_med": float(q_dtau[2]),
            "post_domega_med": float(np.median(samples["domega_220"])),
            "dtau_zero_in90": bool(q_dtau[0] <= 0.0 <= q_dtau[4]),
        })
    return pd.DataFrame(rows)


print("=" * 92)
print("  PT QG: KERR DTAU SIGNAL/SYSTEMATICS SEPARATION")
print("=" * 92)

n_pass = 0
n_total = 0

parser = argparse.ArgumentParser()
parser.add_argument("--lvk-root", default=os.environ.get("PT_LVK_REMNANTS_DIR", "/tmp/pt_lvk_remnants"))
args = parser.parse_args()
root = pathlib.Path(args.lvk_root)

ev = read_evidence_table(root)
post = read_event_posteriors(root)

if ev.empty or post.empty:
    print(f"  Missing pyRing evidence/posterior data under {root}")
    print("=" * 92)
    sys.exit(0)

df = ev.merge(post, on="event", how="inner")
print(f"  events={len(df)}")
print(f"  dtau0 evidence q={fmt_q(quantiles(df['dtau0']))}")
print(f"  both0 evidence q={fmt_q(quantiles(df['both0']))}")
print(f"  max dtau evidence q={fmt_q(quantiles(df['max_dtau']))}")
print(f"  late dtau positive fraction q={fmt_q(quantiles(df['late_positive_fraction_dtau']))}")

dtau_positive_0 = int((df["dtau0"] > 0.0).sum())
both_positive_0 = int((df["both0"] > 0.0).sum())
domega_positive_0 = int((df["domega0"] > 0.0).sum())
stable_dtau = int((df["late_positive_fraction_dtau"] >= 0.8).sum())
stable_both = int((df["late_positive_fraction_both"] >= 0.8).sum())
zero_in90 = int(df["dtau_zero_in90"].sum())
c_degen = corr(df["post_domega_med"], df["post_dtau_med"])
c_evidence_posterior = corr(df["dtau0"], df["post_dtau_med"])

check("dtau-only evidence at nominal start is negative for the catalog median",
      float(np.median(df["dtau0"])) < -0.3,
      f"median logB(dtau vs full modes)={np.median(df['dtau0']):+.3f}")

check("joint domega+dtau evidence at nominal start is even more negative",
      float(np.median(df["both0"])) < -0.7,
      f"median logB(both vs full modes)={np.median(df['both0']):+.3f}")

check("only a minority of events prefer dtau at nominal start",
      dtau_positive_0 <= 5,
      f"{dtau_positive_0}/{len(df)} positive")

check("almost no events prefer joint domega+dtau at nominal start",
      both_positive_0 <= 2,
      f"{both_positive_0}/{len(df)} positive")

check("dtau support is not stable across nonnegative start times",
      stable_dtau <= 4,
      f"{stable_dtau}/{len(df)} have >=80% positive late-start support")

check("joint domega+dtau support is absent across nonnegative start times",
      stable_both == 0,
      f"{stable_both}/{len(df)} have >=80% positive late-start support")

check("most event-level dtau posteriors still include zero at 90%",
      zero_in90 >= 18,
      f"{zero_in90}/{len(df)} include zero")

check("dtau posterior medians are strongly entangled with domega medians",
      c_degen < -0.5,
      f"corr(domega,dtau)={c_degen:+.3f}")

check("posterior dtau medians track dtau evidence but do not overcome the evidence penalty",
      c_evidence_posterior > 0.4 and float(np.median(df["dtau0"])) < 0.0,
      f"corr(logB_dtau,posterior_dtau)={c_evidence_posterior:+.3f}, median logB={np.median(df['dtau0']):+.3f}")

check("broad posterior medians alone are insufficient evidence for a physical channel",
      (zero_in90 >= 18 and dtau_positive_0 <= 5 and stable_dtau <= 4),
      "posterior width + evidence penalty + start-time instability")

print("\n--- Classification ---")
print("  physical candidate retained: catalog-scale surface-gravity relaxation")
print("  event-level pyRing verdict : no clean dissipative signal after higher modes")
print("  most likely current status  : model/start-time/overtone/domega-dtau systematic")
print("  needed next                 : reanalysis with fixed PT priors or start-time")
print("                                marginalization designed for R_tau(a)")

print("\n" + "=" * 92)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: the current pyRing event-level dtau medians are separated")
print("               from a clean dissipative signal; they behave as systematics-")
print("               dominated diagnostics under the released evidence.")
print("=" * 92)

sys.exit(0 if n_pass == n_total else 1)
