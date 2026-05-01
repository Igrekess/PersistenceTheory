#!/usr/bin/env python3
"""
pt_check.py — Unified test framework for Persistence Theory monograph scripts.

Provides the Checker class that replaces 116 duplicate check() definitions
across the v6 scripts. Features:
  - Boolean and numerical tolerance checks
  - Section headers for step-by-step demonstrations
  - Progress bar for long computations
  - Structured JSON report generation
  - Cross-platform UTF-8 output

Usage:
    from lib.pt_check import Checker

    ck = Checker("proof_T1_sieve", chapter="ch01", total_steps=5)
    ck.section("Step 1: Sieve construction")
    ck.check("gaps_count", len(gaps) == 42, f"got {len(gaps)}")
    ck.check_close("alpha_EM", computed, 137.036, tol_pct=0.01, unit="")
    ck.progress(k, k_max, "Computing sieve level")
    ck.summary()
"""

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Cross-platform UTF-8 output
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass  # Python < 3.7


class Checker:
    """Structured test runner with progress bar and JSON report output."""

    def __init__(self, script_name, chapter="", total_steps=0):
        """
        Parameters
        ----------
        script_name : str
            Identifier for this script (used in report filename).
        chapter : str
            Monograph chapter code, e.g. "ch01".
        total_steps : int
            Number of sections (for global progress display). 0 = no bar.
        """
        self.name = script_name
        self.chapter = chapter
        self.total_steps = total_steps
        self.current_step = 0
        self.n_pass = 0
        self.n_fail = 0
        self.results = []
        self._t0 = time.time()
        self._step_t0 = self._t0

    # ------------------------------------------------------------------
    # Section headers
    # ------------------------------------------------------------------
    def section(self, title):
        """Print a section header and update step progress."""
        self.current_step += 1
        elapsed = time.time() - self._t0
        if self.total_steps > 0:
            pct = self.current_step / self.total_steps
            print(
                f"\n{'=' * 72}"
                f"\n  [{self.current_step}/{self.total_steps}] {title}"
                f"  ({elapsed:.1f}s elapsed)"
                f"\n{'=' * 72}"
            )
        else:
            print(f"\n{'=' * 72}\n  {title}\n{'=' * 72}")
        self._step_t0 = time.time()

    # ------------------------------------------------------------------
    # Assertions
    # ------------------------------------------------------------------
    def check(self, label, condition, detail=""):
        """Boolean assertion.

        Parameters
        ----------
        label : str
            Short description of the check.
        condition : bool
            True = PASS, False = FAIL.
        detail : str
            Additional info printed on FAIL.
        """
        status = "PASS" if condition else "FAIL"
        self.results.append({
            "label": label,
            "status": status,
            "detail": str(detail) if detail else "",
        })
        if condition:
            self.n_pass += 1
            print(f"  [PASS] {label}")
        else:
            self.n_fail += 1
            msg = f"  [FAIL] {label}"
            if detail:
                msg += f"  -- {detail}"
            print(msg)
        return condition

    def check_close(self, label, computed, expected, tol_pct=1.0, unit=""):
        """Numerical tolerance check: |computed - expected| / |expected| < tol_pct/100.

        Parameters
        ----------
        label : str
            Short description.
        computed : float
            Value obtained by the script.
        expected : float
            Reference / experimental value.
        tol_pct : float
            Tolerance in percent (e.g. 0.01 means 0.01%).
        unit : str
            Physical unit for display (e.g. "MeV", "eV").
        """
        if expected == 0:
            err_pct = abs(computed) * 100
        else:
            err_pct = abs(computed - expected) / abs(expected) * 100
        ok = err_pct <= tol_pct
        status = "PASS" if ok else "FAIL"

        detail = (
            f"computed={computed:.8g}, expected={expected:.8g}, "
            f"err={err_pct:.4f}%, tol={tol_pct}%"
        )
        if unit:
            detail += f" [{unit}]"

        self.results.append({
            "label": label,
            "status": status,
            "value": float(computed),
            "expected": float(expected),
            "err_pct": round(err_pct, 6),
            "tol_pct": tol_pct,
            "unit": unit,
            "detail": detail,
        })
        if ok:
            self.n_pass += 1
            print(f"  [PASS] {label}  (err={err_pct:.4f}%)")
        else:
            self.n_fail += 1
            print(f"  [FAIL] {label}  -- {detail}")
        return ok

    # ------------------------------------------------------------------
    # Record values (no assertion — just store for the JSON report)
    # ------------------------------------------------------------------
    def record_value(self, label, value, unit="", expected=None, detail=""):
        """Record a computed value in the JSON report without asserting.

        Use this to capture every numerical result from the monograph
        so the JSON report is a complete reproducibility record.

        Parameters
        ----------
        label : str
            Identifier for this value (e.g. "m_W_GeV", "IE_Na_eV").
        value : float
            The computed value.
        unit : str
            Physical unit (e.g. "GeV", "eV", "MeV").
        expected : float or None
            Reference value (PDG, experiment) if available.
        detail : str
            Additional context.
        """
        entry = {
            "label": label,
            "status": "RECORD",
            "value": float(value),
            "unit": unit,
        }
        if expected is not None:
            entry["expected"] = float(expected)
            if expected != 0:
                entry["err_pct"] = round(
                    abs(value - expected) / abs(expected) * 100, 6
                )
        if detail:
            entry["detail"] = detail
        self.results.append(entry)

    # ------------------------------------------------------------------
    # Progress bar (for long computations within a step)
    # ------------------------------------------------------------------
    def progress(self, current, total, label=""):
        """In-step progress bar for long computations.

        Parameters
        ----------
        current : int
            Current iteration (0-based).
        total : int
            Total iterations.
        label : str
            Description of current work.
        """
        if total <= 0:
            return
        pct = (current + 1) / total
        bar_len = 40
        filled = int(bar_len * pct)
        bar = "\u2588" * filled + "\u2591" * (bar_len - filled)
        elapsed = time.time() - self._step_t0
        if pct > 0.01:
            eta = elapsed / pct * (1 - pct)
            eta_str = f"ETA {eta:.0f}s"
        else:
            eta_str = "..."
        sys.stdout.write(f"\r  [{bar}] {pct:5.1%} \u2014 {label}  {eta_str}   ")
        sys.stdout.flush()
        if current + 1 >= total:
            sys.stdout.write("\n")

    # ------------------------------------------------------------------
    # Summary and JSON report
    # ------------------------------------------------------------------
    def summary(self):
        """Print final summary, write JSON report, exit with appropriate code."""
        duration = time.time() - self._t0
        total = self.n_pass + self.n_fail

        print(f"\n{'=' * 72}")
        print(f"  BILAN: {self.name}")
        print(f"  {self.n_pass}/{total} PASS, {self.n_fail} FAIL")
        print(f"  Duration: {duration:.2f}s")
        if self.n_fail == 0:
            print("  STATUS: ALL CHECKS PASSED")
        else:
            print(f"  STATUS: {self.n_fail} FAILURE(S)")
            for r in self.results:
                if r["status"] == "FAIL":
                    print(f"    - {r['label']}: {r.get('detail', '')}")
        print(f"{'=' * 72}")

        self._write_report(duration)
        sys.exit(0 if self.n_fail == 0 else 1)

    def _write_report(self, duration):
        """Write structured JSON report to the companion reports directory."""
        report = {
            "script": self.name,
            "chapter": self.chapter,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "duration_s": round(duration, 3),
            "n_pass": self.n_pass,
            "n_fail": self.n_fail,
            "n_total": self.n_pass + self.n_fail,
            "success": self.n_fail == 0,
            "results": self.results,
            "pt_framework": "persistence-theory",
            "zero_fitted_parameters": True,
        }

        # Determine report directory relative to this file
        lib_dir = Path(__file__).resolve().parent
        scripts_dir = lib_dir.parent
        report_dir = scripts_dir / "reports"
        if self.chapter:
            report_dir = report_dir / self.chapter
        report_dir.mkdir(parents=True, exist_ok=True)

        report_path = report_dir / f"{self.name}.json"
        try:
            with open(report_path, "w", encoding="utf-8") as f:
                json.dump(report, f, indent=2, ensure_ascii=False)
            print(f"  Report: {report_path}")
        except OSError as e:
            print(f"  Warning: could not write report: {e}")
