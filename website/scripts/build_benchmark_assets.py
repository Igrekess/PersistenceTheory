#!/usr/bin/env python3
"""build_benchmark_assets.py — génère les assets pour /benchmark-ptc.

Lit les fichiers de benchmark canoniques dans PT_PROJECTS/PTC/, produit :
- SVG charts (matplotlib) → website/public/benchmark/charts/
- JSON compact pour la table 1000 mol → website/src/data/benchmark1000.json
- JSON métriques agrégées (3 benchmarks) → website/src/data/benchmarkMetrics.json
- CSV téléchargeables → website/public/benchmark/data/

Inputs (chemins absolus) :
- /Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT_PROJECTS/PTC/
  - ptc/data/benchmark_1000.json
  - ptc/data/benchmark_1000_verified.json     (phase 0 output)
  - benchmarkb3lyp/ptc_fresh_2026-05-01.json
  - benchmarkb3lyp/260503_comparison_per_molecule.csv
  - benchmarkb3lyp/260503_comparison_metrics.json
  - benchmark_b3lyp_def2tzvp/comparison_3way_n12_summary.csv
  - benchmark_b3lyp_def2tzvp/comparison_3way_n12_per_molecule.csv
  - benchmark_b3lyp_def2tzvp/comparison_3way_n12_winrates.csv

Run: `python website/scripts/build_benchmark_assets.py`
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# PUBLIC/PT_CHEMISTRY is the canonical, audit-cleaned source (54 wrong CAS
# purged, dual-convention CCCBDB cross-validation done).
PTC_ROOT = Path(
    "/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PUBLIC/PT_CHEMISTRY"
)
AUDIT_FILE = PTC_ROOT / "benchmarkb3lyp" / "audit_dexp_2026-05-03.json"
COMBINED_AUDIT_FILE = PTC_ROOT / "benchmarkb3lyp" / "audit_combined_2026-05-04.json"
WEB_ROOT = Path(__file__).resolve().parent.parent  # website/
CHARTS_OUT = WEB_ROOT / "public" / "benchmark" / "charts"
DATA_OUT_PUB = WEB_ROOT / "public" / "benchmark" / "data"
DATA_OUT_SRC = WEB_ROOT / "src" / "data"

CHARTS_OUT.mkdir(parents=True, exist_ok=True)
DATA_OUT_PUB.mkdir(parents=True, exist_ok=True)
DATA_OUT_SRC.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style constants — match the website's palette
# ---------------------------------------------------------------------------
COLOR_PTC = "#1f4e79"      # thmblue
COLOR_B3LYP = "#c7661f"    # derorange
COLOR_DEF2 = "#7c3aed"     # violet
COLOR_NEUTRAL = "#64748b"
COLOR_GRID = "#e2e8f0"

plt.rcParams.update({
    "font.family": "ui-sans-serif, system-ui, -apple-system, sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": COLOR_GRID,
    "grid.linewidth": 0.6,
    "axes.edgecolor": COLOR_NEUTRAL,
    "axes.labelcolor": "#1e293b",
    "xtick.color": COLOR_NEUTRAL,
    "ytick.color": COLOR_NEUTRAL,
})


def save_svg(name: str) -> None:
    plt.tight_layout()
    out = CHARTS_OUT / f"{name}.svg"
    plt.savefig(out, format="svg", bbox_inches="tight")
    plt.close()
    print(f"  ✓ {out.relative_to(WEB_ROOT)}")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# Source quality tiers (mirroring PT_CHEMISTRY/ptc_app/components/benchmark_panel.py)
_SOURCE_TIER_A = {"ATcT_0K", "ATcT"}
_SOURCE_TIER_B = {"NIST", "JANAF", "HH", "CRC", "Gurvich", "CCCBDB",
                  "NIST-derived", "CRC/Luo"}
_SOURCE_TIER_C = {"3dMLBE20", "Luo2007", "Morse2019", "Bao2021",
                  "Gong2014", "Casey1993", "Peterson1997"}


def _source_tier(source):
    if not source:
        return "?"
    if source in _SOURCE_TIER_A: return "A"
    if source in _SOURCE_TIER_B: return "B"
    if source in _SOURCE_TIER_C: return "C"
    return "?"


def _reliability_label(audit_status, audit_diff_pct, source,
                       combined_status=None, combined_diff=None):
    """Reliability badge for a row.

    Prefers the combined CCCBDB+Burcat status when available :
      🟢🟢 ok_consensus    both sources <2 %
      🟢   ok_one_source   one source <2 %, other n/a
      🟡   one_source_fail one <2 %, other >2 %
      🟠   single_fail     one source >2 %, other n/a
      🔴   conflict        both sources >2 %
      ⚪   n/a             neither source has data — labelled by tier ABC

    Falls back to legacy CCCBDB-only logic when combined audit is missing.
    """
    if combined_status is not None:
        diff = combined_diff
        if combined_status == "ok_consensus":
            return f"🟢🟢 {abs(diff):.1f}%" if diff is not None else "🟢🟢"
        if combined_status == "ok_one_source":
            return f"🟢 {abs(diff):.1f}%" if diff is not None else "🟢"
        if combined_status == "one_source_fail":
            return f"🟡 {abs(diff):.1f}%" if diff is not None else "🟡"
        if combined_status == "single_fail":
            a = abs(diff) if diff is not None else 0
            if a < 5:   return f"🟠 {a:.1f}%"
            if a < 15:  return f"🟠 {a:.1f}%"
            return f"🔴 {a:.0f}%"
        if combined_status == "conflict":
            a = abs(diff) if diff is not None else 0
            return f"🔴 {a:.1f}%"
        if combined_status == "n/a":
            return f"⚪ {_source_tier(source)}"

    s = audit_status or "n/a"
    diff = audit_diff_pct
    if s in ("ok_<1%", "warn_1-2%"):
        return f"🟢 {abs(diff):.1f}%" if diff is not None else "🟢"
    if s == "fail_>2%":
        if diff is None:
            return "🟠 fail"
        a = abs(diff)
        if a < 5:    return f"🟡 {a:.1f}%"
        if a < 15:   return f"🟠 {a:.1f}%"
        return f"🔴 {a:.0f}%"
    return f"⚪ {_source_tier(source)}"


def _n_atoms_from_smiles(smiles):
    """Count total atoms (incl. H) for URL Mask logic. None if rdkit fails."""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        return Chem.AddHs(m).GetNumAtoms()
    except Exception:
        return None


def _best_source_url(cas, name, n_atoms, url_cccbdb):
    """Pick the most authoritative source URL.

    Priority chain (mirrors PT_CHEMISTRY benchmark_panel._best_source_url):
      1. CCCBDB exp atomization page (when CAS validated)
      2. NIST WebBook by CAS with Mask matching compound type
      3. NIST WebBook by name with Mask matching compound type

    Mask=400 (hex) = Constants of diatomic molecules (D₀ visible directly)
    Mask=1   (hex) = Thermochemistry (heats of formation; D_at derived)

    PubChem SMILES URLs avoided (break on 27 % of bracketed SMILES).
    """
    import urllib.parse
    if url_cccbdb:
        return url_cccbdb
    mask = "400" if n_atoms == 2 else "1"
    if cas:
        return (
            f"https://webbook.nist.gov/cgi/cbook.cgi?"
            f"ID={urllib.parse.quote(str(cas))}&Mask={mask}&Units=SI"
        )
    if name:
        return (
            f"https://webbook.nist.gov/cgi/cbook.cgi?"
            f"Name={urllib.parse.quote(str(name))}&Mask={mask}&Units=SI"
        )
    return ""


def load_bench1():
    """Benchmark 1 — PTC vs ATcT/expérimental, post-trim verified dataset.

    Iterates from benchmark_1000_verified.json (the canonical post-cleanup
    dataset, 994 mol after dedup + 5 conflict removals). Joined with :
    - ptc_fresh_2026-05-01.json    : PTC computed values per molecule
    - audit_dexp_2026-05-03.json   : CCCBDB cross-validation
    - audit_combined_2026-05-04.json (optional) : combined CCCBDB+Burcat
    """
    verified = json.loads((PTC_ROOT / "ptc/data/benchmark_1000_verified.json").read_text())
    fresh = json.loads((PTC_ROOT / "benchmarkb3lyp/ptc_fresh_2026-05-01.json").read_text())

    audit = json.loads(AUDIT_FILE.read_text()) if AUDIT_FILE.exists() else {"results": []}
    audit_by_key = {(r["name"], r["smiles"]): r for r in audit.get("results", [])}

    combined = (
        json.loads(COMBINED_AUDIT_FILE.read_text())
        if COMBINED_AUDIT_FILE.exists() else {"results": []}
    )
    combined_by_key = {(r["name"], r["smiles"]): r for r in combined.get("results", [])}

    fresh_by_name = {r["name"].lower(): r for r in fresh["results"]}

    rows = []
    for v in verified:
        name = v["name"]
        smi = v["smiles"]
        d_exp = v["D_exp"]
        cat = v.get("category", "?")
        f = fresh_by_name.get(name.lower(), {})
        a = audit_by_key.get((name, smi), {})
        c = combined_by_key.get((name, smi), {})

        cas = v.get("cas")
        source = v.get("source", "?")
        n_at = _n_atoms_from_smiles(smi)
        audit_status = a.get("status")
        audit_diff = a.get("best_diff_pct")
        audit_conv = a.get("best_convention")
        combined_status = c.get("combined_status")
        combined_diff = c.get("best_diff_pct")

        rows.append({
            "name": name,
            "smiles": smi,
            "category": cat,
            "n_atoms": n_at,
            "d_exp": round(d_exp, 4),
            "d_ptc": round(f.get("d_pt", 0), 4) if f else None,
            "err_eV": round(f.get("abs_err_eV", 0), 4) if f else None,
            "rel_pct": round(f.get("rel_pct", 0), 3) if f else None,
            "source": source,
            "tier": _source_tier(source),
            "cas": cas,
            "unc_eV": v.get("unc_eV"),
            "audit_status": audit_status,
            "audit_diff_pct": round(audit_diff, 2) if audit_diff is not None else None,
            "audit_convention": audit_conv,
            "combined_status": combined_status,
            "combined_diff_pct": round(combined_diff, 2) if combined_diff is not None else None,
            "reliability": _reliability_label(
                audit_status, audit_diff, source,
                combined_status=combined_status, combined_diff=combined_diff,
            ),
            "url": _best_source_url(cas, name, n_at, v.get("url_cccbdb")),
        })

    return rows, fresh, audit, combined


def load_bench2():
    """Benchmark 2 — PTC vs B3LYP/6-31G*, 860 mol."""
    csv_path = PTC_ROOT / "benchmarkb3lyp/260503_comparison_per_molecule.csv"
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "name": r["name"],
                "smiles": r["smiles"],
                "category": r["category"],
                "n_atoms": int(r["n_atoms"]),
                "d_exp": float(r["D_exp_eV"]),
                "d_b3lyp": float(r["D_b3lyp_eV"]),
                "err_b3lyp": float(r["err_b3lyp_eV"]),
                "d_ptc": float(r["D_ptc_eV"]),
                "err_ptc": float(r["err_ptc_eV"]),
                "winner": r["winner"],
            })
    metrics = json.loads(
        (PTC_ROOT / "benchmarkb3lyp/260503_comparison_metrics.json").read_text()
    )
    return rows, metrics


def load_bench3():
    """Benchmark 3 — 3-way comparison, 536 mol n≤12."""
    summary = []
    csv_summary = PTC_ROOT / "benchmark_b3lyp_def2tzvp/comparison_3way_n12_summary.csv"
    with open(csv_summary) as f:
        reader = csv.DictReader(f)
        for r in reader:
            summary.append({
                "bin": r["bin"],
                "n": int(r["n"]),
                "method": r["method"],
                "mae_eV": float(r["MAE_eV"]),
                "mae_kcal": float(r["MAE_kcal"]),
                "mean_signed_eV": float(r["mean_signed_eV"]),
                "within_1kcal_pct": float(r["within_1kcal_pct"]),
                "params": int(r["params"]),
            })
    winrates = []
    csv_wr = PTC_ROOT / "benchmark_b3lyp_def2tzvp/comparison_3way_n12_winrates.csv"
    with open(csv_wr) as f:
        reader = csv.DictReader(f)
        for r in reader:
            winrates.append({
                "bin": r["bin"],
                "n": int(r["n"]),
                "ptc_vs_631gs": float(r["PTC_vs_6-31Gs_pct"]),
                "ptc_vs_def2tzvp": float(r["PTC_vs_def2TZVP_pct"]),
                "def2tzvp_vs_631gs": float(r["def2TZVP_vs_6-31Gs_pct"]),
                "ptc_3way_best": float(r["PTC_3way_best_pct"]),
                "def2tzvp_3way_best": float(r["def2TZVP_3way_best_pct"]),
                "p631gs_3way_best": float(r["6-31Gs_3way_best_pct"]),
            })
    return summary, winrates


# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------
def chart_bench1_histogram(rows):
    """Bench 1 — Histogramme des erreurs PTC vs expérience (1000 mol)."""
    errs = [r["err_eV"] for r in rows if r["err_eV"] is not None]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(errs, bins=60, color=COLOR_PTC, alpha=0.85, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color=COLOR_NEUTRAL, linewidth=0.8, linestyle="--")
    ax.set_xlabel("|erreur| sur D_at (eV)")
    ax.set_ylabel("Nombre de molécules")
    ax.set_title(f"PTC vs expérience — distribution des erreurs absolues ({len(errs)} mol)",
                 fontsize=12, color="#1e293b", pad=10)
    save_svg("bench1_histogram_errors")


def chart_bench1_scatter(rows):
    """Bench 1 — Scatter D_at_PTC vs D_at_exp, coloured by category, with y=x line.

    Inspired by the PTC benchmark panel (ptc_app/components/benchmark_panel.py).
    """
    cat_colors = {
        "Organique":   "#3b82f6",
        "Inorganique": "#f97316",
        "d-block":     "#22c55e",
    }
    fig, ax = plt.subplots(figsize=(8, 7.5))

    # y=x line first so it sits behind the points
    all_d = [r["d_exp"] for r in rows if r["d_ptc"] is not None]
    all_d_ptc = [r["d_ptc"] for r in rows if r["d_ptc"] is not None]
    lo = min(min(all_d), min(all_d_ptc)) * 0.95
    hi = max(max(all_d), max(all_d_ptc)) * 1.05
    ax.plot([lo, hi], [lo, hi], color="#94a3b8", linestyle="--", linewidth=1, zorder=1)

    for cat, color in cat_colors.items():
        xs = [r["d_exp"] for r in rows if r["category"] == cat and r["d_ptc"] is not None]
        ys = [r["d_ptc"] for r in rows if r["category"] == cat and r["d_ptc"] is not None]
        ax.scatter(xs, ys, s=14, color=color, alpha=0.7,
                   label=f"{cat} (n={len(xs)})",
                   edgecolors="white", linewidths=0.3, zorder=2)

    ax.set_xlabel("D_at expérimental (eV)")
    ax.set_ylabel("D_at PTC (eV)")
    ax.set_title("PTC vs expérience — 1000 molécules (ligne y = x = parfaite)",
                 fontsize=12, color="#1e293b", pad=10)
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    save_svg("bench1_scatter_dat")


def chart_bench1_error_bins(rows):
    """Bench 1 — Histogram d'erreurs par bin de %, coloration par sévérité."""
    bins_def = [(0, 1, "#22c55e"),     # green: chemical accuracy
                (1, 2, "#84cc16"),
                (2, 3, "#eab308"),     # yellow
                (3, 5, "#f59e0b"),
                (5, 10, "#f97316"),    # orange
                (10, 1000, "#dc2626")] # red: >10%
    labels = ["0–1 %", "1–2 %", "2–3 %", "3–5 %", "5–10 %", "> 10 %"]
    counts = [0] * 6
    for r in rows:
        if r["rel_pct"] is None:
            continue
        e = r["rel_pct"]
        for i, (lo, hi, _) in enumerate(bins_def):
            if lo <= e < hi:
                counts[i] += 1
                break

    colors = [b[2] for b in bins_def]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(labels, counts, color=colors, alpha=0.9, edgecolor="white", linewidth=0.5)
    for bar, cnt in zip(bars, counts):
        if cnt > 0:
            pct = 100 * cnt / sum(counts)
            ax.text(bar.get_x() + bar.get_width()/2, cnt + max(counts)*0.01,
                    f"{cnt}\n({pct:.0f} %)",
                    ha="center", va="bottom", fontsize=9, color="#334155")
    ax.set_ylabel("Nombre de molécules")
    ax.set_xlabel("Erreur relative |D_PTC − D_exp| / D_exp")
    ax.set_title(f"Distribution des erreurs relatives PTC vs expérience ({sum(counts)} mol)",
                 fontsize=12, color="#1e293b", pad=10)
    ax.set_ylim(0, max(counts) * 1.18)
    save_svg("bench1_error_bins")


def chart_bench1_by_category(fresh):
    """Bench 1 — MAE par catégorie chimique."""
    pc = fresh["per_category"]
    cats = list(pc.keys())
    vals = [pc[c]["mae"] for c in cats]
    counts = [pc[c]["n"] for c in cats]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(cats, vals, color=[COLOR_PTC, COLOR_B3LYP, COLOR_DEF2], alpha=0.85, width=0.55)
    for bar, val, cnt in zip(bars, vals, counts):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.3,
                f"{val:.2f} %\n(n={cnt})",
                ha="center", va="bottom", fontsize=10, color="#334155")
    ax.set_ylabel("MAE relative (%)")
    ax.set_title("PTC — MAE relative par catégorie chimique (1000 mol)",
                 fontsize=12, color="#1e293b", pad=10)
    ax.set_ylim(0, max(vals) * 1.25)
    save_svg("bench1_mae_by_category")


def chart_bench2_overlay(rows):
    """Bench 2 — Histogramme superposé PTC vs B3LYP (860 mol)."""
    errs_ptc = [r["err_ptc"] for r in rows]
    errs_b3 = [r["err_b3lyp"] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = np.linspace(-5, 5, 60)
    ax.hist(errs_b3, bins=bins, color=COLOR_B3LYP, alpha=0.55,
            label=f"B3LYP/6-31G* (MAE 0.93 eV)",
            edgecolor="white", linewidth=0.3)
    ax.hist(errs_ptc, bins=bins, color=COLOR_PTC, alpha=0.75,
            label=f"PTC (MAE 0.60 eV)",
            edgecolor="white", linewidth=0.3)
    ax.axvline(0, color=COLOR_NEUTRAL, linewidth=0.8, linestyle="--")
    ax.set_xlabel("Erreur signée sur D_at (eV)")
    ax.set_ylabel("Nombre de molécules")
    ax.set_title(f"PTC vs B3LYP/6-31G* — distribution des erreurs signées ({len(rows)} mol)",
                 fontsize=12, color="#1e293b", pad=10)
    ax.legend(loc="upper left", frameon=False)
    save_svg("bench2_overlay_errors")


def chart_bench2_mae_by_size(metrics):
    """Bench 2 — MAE par bin de taille (PTC vs B3LYP)."""
    bins_data = metrics["per_size_bin"]
    bins = list(bins_data.keys())
    mae_ptc = [bins_data[b]["ptc"]["mae_eV"] for b in bins]
    mae_b3 = [bins_data[b]["b3lyp"]["mae_eV"] for b in bins]

    x = np.arange(len(bins))
    width = 0.38

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(x - width/2, mae_b3, width, label="B3LYP/6-31G*", color=COLOR_B3LYP, alpha=0.85)
    ax.bar(x + width/2, mae_ptc, width, label="PTC", color=COLOR_PTC, alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(bins, rotation=0)
    ax.set_ylabel("MAE (eV)")
    ax.set_title("MAE par bin de taille (n_atoms) — PTC vs B3LYP/6-31G*",
                 fontsize=12, color="#1e293b", pad=10)
    ax.legend(loc="upper left", frameon=False)
    save_svg("bench2_mae_by_size")


def chart_bench3_3way(summary):
    """Bench 3 — MAE 3-méthodes × 4 bins."""
    bins = ["A (n=2-4)", "B (n=5-8)", "C (n=9-12)", "A+B+C (n≤12)"]
    methods = ["B3LYP/6-31G*", "B3LYP/def2-TZVP", "PTC"]
    colors = {"B3LYP/6-31G*": COLOR_B3LYP, "B3LYP/def2-TZVP": COLOR_DEF2, "PTC": COLOR_PTC}

    data = {m: [] for m in methods}
    for b in bins:
        for m in methods:
            row = next((r for r in summary if r["bin"] == b and r["method"] == m), None)
            data[m].append(row["mae_eV"] if row else 0)

    x = np.arange(len(bins))
    width = 0.27

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(x - width, data["B3LYP/6-31G*"], width, label="B3LYP/6-31G*",
           color=COLOR_B3LYP, alpha=0.85)
    ax.bar(x, data["B3LYP/def2-TZVP"], width, label="B3LYP/def2-TZVP",
           color=COLOR_DEF2, alpha=0.85)
    ax.bar(x + width, data["PTC"], width, label="PTC", color=COLOR_PTC, alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(bins, rotation=0, fontsize=10)
    ax.set_ylabel("MAE (eV)")
    ax.set_title("MAE par bin de taille — 3-way (PTC vs deux bases B3LYP, 536 mol n≤12)",
                 fontsize=12, color="#1e293b", pad=10)
    ax.legend(loc="upper left", frameon=False)
    save_svg("bench3_mae_3way")


def chart_timing(timing):
    """Bench timing — PTC vs B3LYP/def2-TZVP, per bin.

    Log-scale bar chart (the gap is too large for linear scale).
    """
    bins = ["A+B (n≤8)", "C (n=9-12)", "Total (n≤12)"]
    ptc_t = [
        timing["ptc"]["per_bin"]["A+B (n≤8)"]["total_s"],
        timing["ptc"]["per_bin"]["C (n=9-12)"]["total_s"],
        timing["ptc"]["total_s"],
    ]
    def2_t = [
        timing["b3lyp_def2tzvp"]["per_bin"]["A+B (n≤8)"]["total_s"],
        timing["b3lyp_def2tzvp"]["per_bin"]["C (n=9-12)"]["total_s"],
        timing["b3lyp_def2tzvp"]["total_s"],
    ]
    speedup = [d / p for d, p in zip(def2_t, ptc_t)]

    x = np.arange(len(bins))
    width = 0.38

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(x - width/2, def2_t, width, label="B3LYP/def2-TZVP",
           color=COLOR_DEF2, alpha=0.9)
    ax.bar(x + width/2, ptc_t, width, label="PTC", color=COLOR_PTC, alpha=0.9)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_ylabel("Temps total de calcul (s, échelle log)")
    ax.set_title("Coût de calcul — PTC vs B3LYP/def2-TZVP",
                 fontsize=12, color="#1e293b", pad=10)
    ax.legend(loc="upper left", frameon=False)

    # Annotate speedup factors
    for i, (xv, sp, dv) in enumerate(zip(x, speedup, def2_t)):
        ax.text(xv, dv * 1.6, f"× {sp:,.0f}",
                ha="center", va="bottom", fontsize=10,
                color="#dc2626", fontweight="bold")

    save_svg("timing_log")


def chart_bench3_winrates(winrates):
    """Bench 3 — Win rate stacked bar (3-way best)."""
    bins = [w["bin"] for w in winrates]
    ptc = [w["ptc_3way_best"] for w in winrates]
    def2 = [w["def2tzvp_3way_best"] for w in winrates]
    b631 = [w["p631gs_3way_best"] for w in winrates]

    x = np.arange(len(bins))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(x, ptc, label="PTC le meilleur", color=COLOR_PTC, alpha=0.9)
    ax.bar(x, def2, bottom=ptc, label="def2-TZVP le meilleur", color=COLOR_DEF2, alpha=0.85)
    bot2 = [p + d for p, d in zip(ptc, def2)]
    ax.bar(x, b631, bottom=bot2, label="6-31G* le meilleur", color=COLOR_B3LYP, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(bins, rotation=0, fontsize=10)
    ax.set_ylabel("% de molécules où la méthode est la meilleure")
    ax.set_title("Win rate 3-way par bin de taille",
                 fontsize=12, color="#1e293b", pad=10)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_ylim(0, 100)
    save_svg("bench3_winrates")


# ---------------------------------------------------------------------------
# Compact JSON outputs
# ---------------------------------------------------------------------------
def write_bench1_json(rows, fresh, audit, combined):
    """Write benchmark1000.json. Stats are recomputed from `rows` so that
    headline numbers match the post-trim dataset (994 mol, not 1000).
    """
    from collections import Counter

    # Recompute headline stats from rows actually shipped
    valid = [r for r in rows if r["d_ptc"] is not None]
    n_total = len(rows)
    n_fails = n_total - len(valid)
    if valid:
        errs_eV  = [r["err_eV"]  for r in valid if r["err_eV"]  is not None]
        errs_pct = [r["rel_pct"] for r in valid if r["rel_pct"] is not None]
        mae_eV  = sum(errs_eV) / len(errs_eV) if errs_eV else 0.0
        mae_pct = sum(errs_pct) / len(errs_pct) if errs_pct else 0.0
        # Chemical accuracy : |err| < 1 kcal/mol = 0.0434 eV
        chem_acc = sum(1 for e in errs_eV if e < 1.0/23.0605) / len(errs_eV) * 100
    else:
        mae_eV = mae_pct = chem_acc = 0.0
    mae_kcal = mae_eV * 23.0605

    # Per-category stats
    by_cat = {}
    for r in valid:
        c = r["category"]
        by_cat.setdefault(c, []).append(r["rel_pct"])
    per_cat = {
        c: {"n": len(arr), "mae_pct": round(sum(arr)/len(arr), 3) if arr else 0}
        for c, arr in by_cat.items()
    }

    audit_counts = Counter(r.get("status") for r in audit.get("results", []))
    combined_counts = Counter(r.get("combined_status") for r in combined.get("results", []))

    tier_counts = Counter(r["tier"] for r in rows)
    rel_buckets = {"green": 0, "yellow": 0, "orange": 0, "red": 0, "white": 0}
    for r in rows:
        b = r["reliability"][:1]
        if b == "🟢": rel_buckets["green"] += 1
        elif b == "🟡": rel_buckets["yellow"] += 1
        elif b == "🟠": rel_buckets["orange"] += 1
        elif b == "🔴": rel_buckets["red"] += 1
        else: rel_buckets["white"] += 1

    out = {
        "n_total": n_total,
        "n_fails": n_fails,
        "mae_eV": round(mae_eV, 4),
        "mae_kcal": round(mae_kcal, 2),
        "mae_pct": round(mae_pct, 3),
        "chemical_acc_pct": round(chem_acc, 2),
        "per_category": per_cat,
        "audit": {
            "n_audited":      sum(audit_counts.values()),
            "ok_lt_1pct":     audit_counts.get("ok_<1%", 0),
            "warn_1_2pct":    audit_counts.get("warn_1-2%", 0),
            "fail_gt_2pct":   audit_counts.get("fail_>2%", 0),
            "n_a":            audit_counts.get("n/a", 0),
        },
        "audit_combined": {
            "n_audited":       sum(combined_counts.values()),
            "ok_consensus":    combined_counts.get("ok_consensus", 0),
            "ok_one_source":   combined_counts.get("ok_one_source", 0),
            "one_source_fail": combined_counts.get("one_source_fail", 0),
            "single_fail":     combined_counts.get("single_fail", 0),
            "conflict":        combined_counts.get("conflict", 0),
            "n_a":             combined_counts.get("n/a", 0),
        },
        "tiers": dict(tier_counts),
        "reliability_buckets": rel_buckets,
        "rows": rows,
    }
    out_path = DATA_OUT_SRC / "benchmark1000.json"
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"  ✓ {out_path.relative_to(WEB_ROOT)}  ({len(rows)} rows)")


def compute_timing() -> dict:
    """Compute PTC timing on the 536-mol 3-way subset, attach B3LYP/def2-TZVP
    timings provided by the user (logs from 2026-05-03 run on M-series Mac).

    Same machine class for both methods (Mac M-series), same molecule subset.
    """
    import csv as _csv
    import time as _time
    import statistics as _stat
    import sys as _sys
    _sys.path.insert(0, str(PTC_ROOT))
    from ptc.api import Molecule

    csv_path = PTC_ROOT / "benchmark_b3lyp_def2tzvp/comparison_3way_n12_per_molecule.csv"
    mols = []
    with open(csv_path) as f:
        for row in _csv.DictReader(f):
            mols.append({"name": row.get("name"), "smiles": row.get("smiles"),
                         "n_atoms": int(row["n_atoms"])})

    bins = {"A (n=2-4)": [], "B (n=5-8)": [], "C (n=9-12)": []}
    all_t = []
    fails = 0
    for m in mols:
        try:
            t0 = _time.perf_counter()
            mol = Molecule(m["smiles"])
            _ = mol.D_at
            t = _time.perf_counter() - t0
            all_t.append(t)
            n = m["n_atoms"]
            if n <= 4: bins["A (n=2-4)"].append(t)
            elif n <= 8: bins["B (n=5-8)"].append(t)
            elif n <= 12: bins["C (n=9-12)"].append(t)
        except Exception:
            fails += 1

    bin_a_total = sum(bins["A (n=2-4)"])
    bin_b_total = sum(bins["B (n=5-8)"])
    bin_c_total = sum(bins["C (n=9-12)"])
    total = sum(all_t)

    return {
        "subset": "536 mol, n ≤ 12 (Bin A+B+C from 3-way benchmark)",
        "machine_note": "Apple M1 Max (single-thread for both methods)",
        "ptc": {
            "n": len(all_t),
            "fails": fails,
            "total_s": round(total, 3),
            "median_ms": round(1000 * _stat.median(all_t), 2),
            "mean_ms": round(1000 * _stat.mean(all_t), 2),
            "max_ms": round(1000 * max(all_t), 1),
            "per_bin": {
                "A (n=2-4)":   {"n": len(bins["A (n=2-4)"]),  "total_s": round(bin_a_total, 3)},
                "B (n=5-8)":   {"n": len(bins["B (n=5-8)"]),  "total_s": round(bin_b_total, 3)},
                "A+B (n≤8)":   {"n": len(bins["A (n=2-4)"]) + len(bins["B (n=5-8)"]),
                                "total_s": round(bin_a_total + bin_b_total, 3)},
                "C (n=9-12)":  {"n": len(bins["C (n=9-12)"]), "total_s": round(bin_c_total, 3)},
            },
        },
        # B3LYP/def2-TZVP timings — from user's CLI logs 2026-05-03
        "b3lyp_def2tzvp": {
            "n": 536,
            "fails": 0,
            "total_s": 6970.8,         # 1h56 on 536 mol
            "per_bin": {
                "A+B (n≤8)":   {"n": 369, "total_s": 1853.4},   # 30min53
                "C (n=9-12)":  {"n": 167, "total_s": 5117.4},   # 6970.8 - 1853.4
            },
            "source": "User CLI logs from `benchmark_b3lyp_def2tzvp_n4` runs, 2026-05-03",
        },
    }


def write_metrics_json(fresh, b2_metrics, b3_summary, b3_winrates, timing,
                       bench1_rows=None):
    """Aggregated headline metrics for the page hero. Bench1 numbers are
    recomputed from the trimmed dataset rows when provided."""
    g = b2_metrics["global"]

    # Bench1 — recompute from rows shipped (post-trim)
    if bench1_rows:
        valid = [r for r in bench1_rows if r["d_ptc"] is not None]
        n_b1 = len(valid)
        errs_eV = [r["err_eV"] for r in valid if r["err_eV"] is not None]
        errs_pct = [r["rel_pct"] for r in valid if r["rel_pct"] is not None]
        mae_eV = sum(errs_eV) / len(errs_eV) if errs_eV else 0.0
        mae_pct = sum(errs_pct) / len(errs_pct) if errs_pct else 0.0
        mae_kcal = mae_eV * 23.0605
        chem_acc = sum(1 for e in errs_eV if e < 1.0/23.0605) / len(errs_eV) * 100 if errs_eV else 0
        by_cat = {}
        for r in valid:
            by_cat.setdefault(r["category"], []).append(r["rel_pct"])
        per_cat = {
            c: {"n": len(arr), "mae_pct": round(sum(arr)/len(arr), 2) if arr else 0}
            for c, arr in by_cat.items()
        }
    else:
        n_b1 = fresh["n_mol"]
        mae_eV = fresh["mae_eV"]
        mae_kcal = fresh["mae_kcal"]
        mae_pct = fresh["mae_pct"]
        chem_acc = fresh["chemical_acc_pct"]
        per_cat = {cat: {"n": d["n"], "mae_pct": round(d["mae"], 2)}
                   for cat, d in fresh["per_category"].items()}

    out = {
        "bench1": {
            "n": n_b1,
            "mae_eV": round(mae_eV, 3),
            "mae_kcal": round(mae_kcal, 2),
            "mae_pct": round(mae_pct, 2),
            "chem_acc_pct": round(chem_acc, 1),
            "per_category": per_cat,
        },
        "bench2": {
            "n": g["b3lyp"]["n"],
            "ptc": {
                "mae_eV": round(g["ptc"]["mae_eV"], 3),
                "median_eV": round(g["ptc"]["median_eV"], 3),
                "chem_acc_pct": round(g["ptc"]["chem_acc_pct"], 1),
                "within_5kcal": round(g["ptc"]["within_5kcal"], 1),
                "max_eV": round(g["ptc"]["max_eV"], 2),
            },
            "b3lyp": {
                "mae_eV": round(g["b3lyp"]["mae_eV"], 3),
                "median_eV": round(g["b3lyp"]["median_eV"], 3),
                "chem_acc_pct": round(g["b3lyp"]["chem_acc_pct"], 1),
                "within_5kcal": round(g["b3lyp"]["within_5kcal"], 1),
                "max_eV": round(g["b3lyp"]["max_eV"], 2),
            },
            "ptc_win_rate_pct": round(g["head_to_head"]["ptc_win_rate_pct"], 1),
        },
        "bench3": {
            "summary": b3_summary,
            "winrates": b3_winrates,
        },
        "timing": timing,
    }
    out_path = DATA_OUT_SRC / "benchmarkMetrics.json"
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"  ✓ {out_path.relative_to(WEB_ROOT)}")


def copy_csvs():
    """Copy the canonical CSVs to public/ for direct download."""
    import shutil
    pairs = [
        (PTC_ROOT / "benchmarkb3lyp/260503_comparison_per_molecule.csv",
         DATA_OUT_PUB / "benchmark2_per_molecule.csv"),
        (PTC_ROOT / "benchmark_b3lyp_def2tzvp/comparison_3way_n12_per_molecule.csv",
         DATA_OUT_PUB / "benchmark3_3way_per_molecule.csv"),
        (PTC_ROOT / "benchmark_b3lyp_def2tzvp/comparison_3way_n12_summary.csv",
         DATA_OUT_PUB / "benchmark3_3way_summary.csv"),
    ]
    for src, dst in pairs:
        shutil.copy(src, dst)
        print(f"  ✓ {dst.relative_to(WEB_ROOT)}")

    # Bench 1 export: per-mol CSV with sources + URLs + reliability
    bench1_rows, _, _, _ = load_bench1()
    csv_out = DATA_OUT_PUB / "benchmark1_per_molecule.csv"
    with open(csv_out, "w", newline="") as f:
        if not bench1_rows:
            return
        writer = csv.DictWriter(f, fieldnames=list(bench1_rows[0].keys()))
        writer.writeheader()
        writer.writerows(bench1_rows)
    print(f"  ✓ {csv_out.relative_to(WEB_ROOT)}  ({len(bench1_rows)} rows)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("Loading data...")
    bench1_rows, fresh, audit, combined = load_bench1()
    bench2_rows, b2_metrics = load_bench2()
    b3_summary, b3_winrates = load_bench3()
    print(f"  bench1: {len(bench1_rows)} rows, fresh stats loaded")
    print(f"  bench2: {len(bench2_rows)} rows, metrics loaded")
    print(f"  bench3: {len(b3_summary)} summary lines, {len(b3_winrates)} winrate lines")

    print("\nMeasuring PTC timing on the 536-mol 3-way subset...")
    timing = compute_timing()
    print(f"  PTC: {timing['ptc']['total_s']} s on {timing['ptc']['n']} mol "
          f"(median {timing['ptc']['median_ms']} ms)")
    print(f"  B3LYP/def2-TZVP: {timing['b3lyp_def2tzvp']['total_s']} s "
          f"(from user logs)")

    print("\nGenerating charts...")
    chart_bench1_histogram(bench1_rows)
    chart_bench1_scatter(bench1_rows)
    chart_bench1_error_bins(bench1_rows)
    chart_bench1_by_category(fresh)
    chart_bench2_overlay(bench2_rows)
    chart_bench2_mae_by_size(b2_metrics)
    chart_bench3_3way(b3_summary)
    chart_bench3_winrates(b3_winrates)
    chart_timing(timing)

    print("\nWriting JSON data...")
    write_bench1_json(bench1_rows, fresh, audit, combined)
    write_metrics_json(fresh, b2_metrics, b3_summary, b3_winrates, timing,
                       bench1_rows=bench1_rows)

    print("\nCopying canonical CSVs to public/ for download...")
    copy_csvs()

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
