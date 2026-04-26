#!/usr/bin/env python3
"""
run_all.py -- Master test runner for Persistence Theory monograph scripts (v7).

Discovers and runs all proof/test scripts in chapter subdirectories.
Reports per-domain results, generates aggregate summary.csv from JSON reports.

Usage:
    python run_all.py                  # run all domains
    python run_all.py --tree           # show full arborescence (no run)
    python run_all.py --list           # flat list of all scripts (no run)
    python run_all.py ch01             # run scripts matching 'ch01'
    python run_all.py --summary        # aggregate JSON reports into summary.csv
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# Domains in monograph order
DOMAIN_LABELS = {
    "ch01_sieve":          "Part I   Sieve (ch01) -- T1, T3, s=1/2",
    "ch02_uniqueness":     "Part I   Uniqueness (ch02) -- T6",
    "ch03_conservation":   "Part I   Conservation (ch03) -- q+/q-",
    "ch04_gft":            "Part I   GFT (ch04) -- log2(m)=D_KL+H",
    "ch05_geometry":       "Part I   Geometry (ch05) -- G1, G5 Fisher",
    "ch06_holonomy":       "Part I   Holonomy (ch06) -- sin2, gamma_p",
    "ch07_convergence":    "Part II  Convergence (ch07) -- T4",
    "ch08_fixed_point":    "Part II  Fixed point (ch08) -- T5: mu*=15",
    "ch_math_structures":  "Part II  Math structures -- M05-M35",
    "ch_PM":               "Part II  Predictive Mathematics",
    "ch25_BA0_closing":    "Part II  BA0 Closing (ch25) -- T0",
    "ch09_bridge":         "Part III Bridge (ch09) -- BA0-BA5",
    "ch10_fine_structure":  "Part III Fine structure (ch10) -- alpha_EM",
    "ch11_couplings":      "Part III Couplings (ch11) -- sin2W, alpha_s",
    "ch12_quantum":        "Part IV  Quantum (ch12) -- Born, Schrodinger",
    "ch13_relativity":     "Part IV  Relativity (ch13) -- SO(3,1), GR",
    "ch14_thermodynamics": "Part IV  Thermodynamics (ch14) -- KMS",
    "ch15_sm_observables": "Part IV  SM observables (ch15) -- 43 obs",
    "ch16_perturbative":   "Part V   Perturbative (ch16) -- NLO/NNLO",
    "ch17_feynman":        "Part V   Feynman (ch17) -- F1-F7",
    "ch18_quarkonium":     "Part V   Quarkonium (ch18)",
    "ch19_hadrons":        "Part V   Hadrons (ch19)",
    "ch20_collider":       "Part V   Collider (ch20) -- 72 obs",
    "ch21_predictions":    "Part V   Predictions (ch21) -- P1-P28",
    "ch22_chemistry":      "Part VI  Chemistry (ch22) -- C1-C14",
    # Part VII: Cosmological Scale (existing chapters)
    "ch20f_cosmological_DM":  "Part VII Cosmological DM (ch20f) -- F_inactive",
    "ch20f_dark_energy":      "Part VII Dark energy (ch20f) -- P20 routes A/B/D",
    "ch20f_neutrino_mass_bound": "Part VII Neutrino mass bound (ch20f) -- DESI DR3",
    # Part VIII: Class B Extensions (added 2026-04-26, 13 chapters)
    "ch20b_wilson_coefficients":      "Part VIII Wilson |C9-C10|=beta_echo (ch20b) [PRED]",
    "ch20c_hadronic_margin":          "Part VIII Hadronic margin b->sll (ch20c) [DER-PHYS]",
    "ch20d_BSM_taxonomy":             "Part VIII BSM taxonomy A/B/C (ch20d) [META]",
    "ch20e_DM_candidate":             "Part VIII DM scalar singlet p=2 (ch20e) [PRED]",
    "ch20e_RH_neutrinos":             "Part VIII Right-handed neutrinos Dirac (ch20e)",
    "ch20e_basin_robustness":         "Part VIII Basin robustness mu*=15 (ch20e) [THM+VAL]",
    "ch20f_hubble_tension":           "Part VIII Hubble tension H_0=67.43 (ch20f) [PRED]",
    "ch20f_meerkat_filament":         "Part VIII MeerKAT filament alignment (ch20f) [VAL]",
    "ch20g_bobeth_181_structure":     "Part VIII Bobeth 181-prefactor (ch20g) [DETAIL]",
    "ch20g_bottom_loop_topology":     "Part VIII Bottom-loop topology b->sgamma (ch20g) [DER]",
    "ch20g_counting_convention":      "Part VIII Counting convention mu*=4Nc+3 (ch20g)",
    "ch20g_higgs_portal_derivation":  "Part VIII Higgs portal lambda_HS (ch20g) [DER]",
    "ch20g_super_echo_cutoff":        "Part VIII Super-echo cutoff {11..23} (ch20g) [DER]",
    # Part IX: Audit & Appendices
    "ch23_audit":          "Part IX  Audit (ch23)",
    "verify_sota":         "Part IX  SOTA verification (requires ptc)",
}

# Prefixes that indicate library/utility modules (not runnable scripts)
_SKIP_PREFIXES = ("pt_", "_", "__", "conftest")
_SKIP_FILES = {"run_all.py", "conftest.py"}


def discover_domains():
    """Auto-discover all subdirectories containing .py scripts."""
    domains = []
    for entry in sorted(os.listdir(SCRIPT_DIR)):
        path = SCRIPT_DIR / entry
        if not path.is_dir():
            continue
        if entry.startswith((".", "__", "archive", "lib", "reports")):
            continue
        # Check it has at least one runnable script
        has_script = any(
            f.endswith(".py") and not any(f.startswith(p) for p in _SKIP_PREFIXES)
            for f in os.listdir(path) if os.path.isfile(path / f)
        )
        if has_script:
            domains.append(entry)
    return domains


def find_scripts(domain_dir):
    """Find all runnable scripts in a domain directory."""
    scripts = []
    for f in sorted(os.listdir(domain_dir)):
        path = domain_dir / f
        if not path.is_file() or not f.endswith(".py"):
            continue
        if any(f.startswith(p) for p in _SKIP_PREFIXES):
            continue
        if f in _SKIP_FILES:
            continue
        scripts.append(path)
    return scripts


def run_script(path, timeout=300):
    """Run a Python script, return (success, output, duration)."""
    t0 = time.time()
    env = {**os.environ, "PYTHONPATH": str(SCRIPT_DIR)}
    try:
        result = subprocess.run(
            [sys.executable, str(path)],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(SCRIPT_DIR),
            env=env,
            encoding="utf-8",
            errors="replace",
        )
        dt = time.time() - t0
        ok = result.returncode == 0
        out = result.stdout + result.stderr
        return ok, out, dt
    except subprocess.TimeoutExpired:
        return False, f"TIMEOUT after {timeout}s", time.time() - t0
    except Exception as e:
        return False, str(e), time.time() - t0


def show_tree(domains, filter_str=None):
    """Print the full arborescence of scripts."""
    total = 0
    for d in domains:
        label = DOMAIN_LABELS.get(d, d)
        domain_dir = SCRIPT_DIR / d
        scripts = find_scripts(domain_dir)
        if not scripts:
            continue
        if filter_str and filter_str.lower() not in d.lower():
            continue
        print(f"\n{d}/  ({len(scripts)} scripts) -- {label}")
        for i, s in enumerate(scripts):
            is_last = (i == len(scripts) - 1)
            connector = "\u2514\u2500\u2500" if is_last else "\u251c\u2500\u2500"
            print(f"  {connector} {s.name}")
            total += 1
    print(f"\n  Total: {total} scripts")


def _load_all_reports():
    """Load all JSON reports, return list of (chapter, script, report_dict)."""
    reports_dir = SCRIPT_DIR / "reports"
    reports = []
    for chapter_dir in sorted(reports_dir.iterdir()):
        if not chapter_dir.is_dir():
            continue
        for json_file in sorted(chapter_dir.glob("*.json")):
            try:
                with open(json_file, encoding="utf-8") as f:
                    report = json.load(f)
                reports.append((
                    report.get("chapter", ""),
                    report.get("script", ""),
                    report,
                ))
            except (json.JSONDecodeError, OSError):
                continue
    return reports


def aggregate_reports():
    """Generate complete report: summary.csv, values.csv, and readable text."""
    reports_dir = SCRIPT_DIR / "reports"
    all_reports = _load_all_reports()
    if not all_reports:
        print("No reports found.")
        return

    # ── 1. summary.csv: one row per script ──
    summary_rows = []
    for chapter, script, rpt in all_reports:
        summary_rows.append({
            "chapter": chapter,
            "script": script,
            "n_pass": rpt.get("n_pass", 0),
            "n_fail": rpt.get("n_fail", 0),
            "n_total": rpt.get("n_total", 0),
            "success": rpt.get("success", False),
            "duration_s": rpt.get("duration_s", 0),
            "timestamp": rpt.get("timestamp", ""),
        })

    csv_path = reports_dir / "summary.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    # ── 2. values.csv: one row per numerical result ──
    val_rows = []
    for chapter, script, rpt in all_reports:
        for r in rpt.get("results", []):
            if "value" not in r:
                continue
            val_rows.append({
                "chapter": chapter,
                "script": script,
                "label": r.get("label", ""),
                "value": r.get("value", ""),
                "expected": r.get("expected", ""),
                "err_pct": r.get("err_pct", ""),
                "unit": r.get("unit", ""),
                "status": r.get("status", ""),
            })

    val_path = reports_dir / "values.csv"
    if val_rows:
        with open(val_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=val_rows[0].keys())
            writer.writeheader()
            writer.writerows(val_rows)

    # ── 3. Readable text report ──
    total_pass = sum(r["n_pass"] for r in summary_rows)
    total_fail = sum(r["n_fail"] for r in summary_rows)
    total_ok = sum(1 for r in summary_rows if r["success"])
    total_dur = sum(r["duration_s"] for r in summary_rows)
    n_scripts = len(summary_rows)

    print(f"\n{'=' * 78}")
    print(f"  PERSISTENCE THEORY v7 -- COMPLETE REPORT")
    print(f"  Zero fitted parameters. All derived from s = 1/2 (Theorem T1).")
    print(f"{'=' * 78}")

    # Per-chapter table
    print(f"\n  {'Chapter':<20} {'Script':<35} {'Checks':>7} {'Status':>8} {'Time':>6}")
    print(f"  {'-'*20} {'-'*35} {'-'*7} {'-'*8} {'-'*6}")

    current_part = ""
    for sr in summary_rows:
        ch = sr["chapter"]
        # Detect part changes
        part = ""
        if ch.startswith("ch0") and ch <= "ch06":
            part = "Part I: Arithmetic Core"
        elif ch in ("ch07", "ch08", "ch25", "ch_PM", "ch_math"):
            part = "Part II: PT Closure"
        elif ch.startswith("ch09") or ch.startswith("ch10") or ch.startswith("ch11"):
            part = "Part III: Bridge"
        elif ch.startswith("ch1") and ch >= "ch12" and ch <= "ch15":
            part = "Part IV: Reconstruction"
        elif ch.startswith("ch1") and ch >= "ch16" or ch.startswith("ch2") and ch <= "ch21":
            part = "Part V: Validation"
        elif ch.startswith("ch22"):
            part = "Part VI: Chemistry"
        elif ch.startswith("ch23") or ch.startswith("verify"):
            part = "Part VII: Audit"

        if part and part != current_part:
            current_part = part
            print(f"\n  -- {part} --")

        tag = "OK" if sr["success"] else "FAIL"
        print(f"  {ch:<20} {sr['script']:<35} {sr['n_total']:>5}   {tag:>6}  {sr['duration_s']:>5.1f}s")

    # Totals
    print(f"\n  {'-'*78}")
    print(f"  {'TOTAL':<20} {n_scripts:>2} scripts{'':<24} {total_pass+total_fail:>5}   "
          f"{'ALL OK' if total_fail == 0 else f'{total_fail} FAIL':>6}  {total_dur:>5.1f}s")

    # Numerical values by chapter (compact)
    if val_rows:
        print(f"\n{'=' * 78}")
        print(f"  NUMERICAL RESULTS ({len(val_rows)} values)")
        print(f"{'=' * 78}")

        current_ch = ""
        for vr in val_rows:
            if vr["chapter"] != current_ch:
                current_ch = vr["chapter"]
                print(f"\n  [{current_ch}] {vr['script']}")
                print(f"  {'Label':<35} {'PT Value':>15} {'Expected':>15} {'Err%':>10} {'Unit':<8}")
                print(f"  {'-'*35} {'-'*15} {'-'*15} {'-'*10} {'-'*8}")

            val = vr["value"]
            exp = vr["expected"]
            err = vr["err_pct"]
            unit = vr["unit"] or ""
            label = vr["label"][:35]

            val_s = f"{val:>15.8g}" if isinstance(val, float) else f"{val!s:>15}"
            exp_s = f"{exp:>15.8g}" if isinstance(exp, (int, float)) and exp != "" else f"{exp!s:>15}"
            err_s = f"{err:>9.4f}%" if isinstance(err, (int, float)) and err != "" else f"{'':>10}"

            print(f"  {label:<35} {val_s} {exp_s} {err_s} {unit:<8}")

    # Files generated
    print(f"\n{'=' * 78}")
    print(f"  FILES GENERATED")
    print(f"{'=' * 78}")
    print(f"  {csv_path}  ({n_scripts} scripts)")
    if val_rows:
        print(f"  {val_path}  ({len(val_rows)} numerical results)")
    for chapter, script, _ in all_reports:
        ch_dir = reports_dir / chapter if chapter else reports_dir
        print(f"  {ch_dir / (script + '.json')}")

    print(f"\n  Zero fitted parameters. Zero inputs.")
    print(f"  s = 1/2 is DERIVED (Theorem T1), not assumed.")
    print(f"  m_e = translation factor (like c = 3e8 m/s).")
    print()


def main():
    parser = argparse.ArgumentParser(description="PT v7 monograph test runner")
    parser.add_argument("filter", nargs="?", default=None,
                        help="Filter domains by name fragment")
    parser.add_argument("--tree", action="store_true",
                        help="Show arborescence only (no run)")
    parser.add_argument("--list", action="store_true",
                        help="Flat list of scripts (no run)")
    parser.add_argument("--summary", action="store_true",
                        help="Aggregate JSON reports into summary.csv")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-script timeout in seconds (default: 300)")
    args = parser.parse_args()

    domains = discover_domains()

    if args.summary:
        aggregate_reports()
        return

    if args.tree:
        show_tree(domains, args.filter)
        return

    if args.list:
        for d in domains:
            domain_dir = SCRIPT_DIR / d
            for s in find_scripts(domain_dir):
                rel = s.relative_to(SCRIPT_DIR)
                if args.filter and args.filter.lower() not in str(rel).lower():
                    continue
                print(rel)
        return

    # Run mode
    t0_global = time.time()
    global_pass = 0
    global_fail = 0
    results_table = []

    for d in domains:
        label = DOMAIN_LABELS.get(d, d)
        domain_dir = SCRIPT_DIR / d
        scripts = find_scripts(domain_dir)
        if not scripts:
            continue
        if args.filter and args.filter.lower() not in d.lower():
            continue

        print(f"\n{'=' * 72}")
        print(f"  {label}")
        print(f"{'=' * 72}")

        domain_pass = 0
        domain_fail = 0

        for script_path in scripts:
            ok, output, dt = run_script(script_path, timeout=args.timeout)
            tag = "PASS" if ok else "FAIL"
            print(f"  [{tag}] {script_path.name}  ({dt:.1f}s)")

            if ok:
                domain_pass += 1
                global_pass += 1
            else:
                domain_fail += 1
                global_fail += 1
                # Show last few lines of output on failure
                lines = output.strip().split("\n")
                for line in lines[-5:]:
                    print(f"         {line}")

            results_table.append((d, script_path.name, tag, dt))

        total_d = domain_pass + domain_fail
        print(f"  --- {d}: {domain_pass}/{total_d} PASS")

    # Global summary
    dt_global = time.time() - t0_global
    total = global_pass + global_fail
    print(f"\n{'=' * 72}")
    print(f"  GLOBAL: {global_pass}/{total} PASS, {global_fail} FAIL")
    print(f"  Duration: {dt_global:.1f}s")
    print(f"{'=' * 72}")

    # Aggregate reports
    aggregate_reports()

    sys.exit(0 if global_fail == 0 else 1)


if __name__ == "__main__":
    main()
