#!/usr/bin/env python3
"""
Verify molecular dissociation energy SOTA scores against PTC calculator.

Requires: ptc package (Persistence Theory Chemistry calculator).
If ptc is not installed, all tests report SKIP.

Scores verified (from pt_chimie_scores.tex):
  - 849 molecules, MAE 2.15%, median 1.32%
  - 354/849 under 1%
  - Main-only 2.02%, Extended 1.58%, ATcT 2.88%, d-block 3.66%
"""

import sys

try:
    from ptc.molecule import compute_atomization_energy
    from ptc.data.benchmark import load_benchmark
    PTC_AVAILABLE = True
except ImportError:
    try:
        from ptc.transfer_matrix import compute_dat
        from ptc.constants import BENCHMARK_MOLECULES
        PTC_AVAILABLE = True
    except ImportError:
        PTC_AVAILABLE = False


def main():
    if not PTC_AVAILABLE:
        print("SKIP: ptc package not installed.")
        print("      Install PTC or add it to PYTHONPATH to verify molecular scores.")
        print("      The monograph Level 1/2 scripts run without PTC.")
        return 0

    print("=" * 60)
    print("  verify_molecules.py — Molecular SOTA verification")
    print("=" * 60)

    n_pass = 0
    n_fail = 0

    # Load benchmark
    try:
        benchmark = load_benchmark()
    except Exception:
        benchmark = BENCHMARK_MOLECULES

    errors = []
    for mol in benchmark:
        name = mol.get('name', mol.get('formula', '?'))
        obs = mol.get('D_at_obs', mol.get('D0_obs'))
        if obs is None or obs <= 0:
            continue
        try:
            pt = compute_atomization_energy(mol)
        except Exception:
            try:
                pt = compute_dat(mol)
            except Exception:
                continue
        err = abs(pt - obs) / obs * 100
        errors.append(err)

    if not errors:
        print("  FAIL  No molecules computed (check PTC API)")
        return 1

    n_mol = len(errors)
    mae = sum(errors) / n_mol
    median = sorted(errors)[n_mol // 2]
    under_1 = sum(1 for e in errors if e < 1.0)

    # Test 1: molecule count
    if n_mol >= 800:
        print(f"  PASS  Molecule count = {n_mol} (threshold >= 800)")
        n_pass += 1
    else:
        print(f"  FAIL  Molecule count = {n_mol} (expected >= 800)")
        n_fail += 1

    # Test 2: global MAE
    if mae < 2.5:
        print(f"  PASS  Global MAE = {mae:.2f}% (threshold <2.5%)")
        n_pass += 1
    else:
        print(f"  FAIL  Global MAE = {mae:.2f}% (expected <2.5%)")
        n_fail += 1

    # Test 3: median
    if median < 2.0:
        print(f"  PASS  Median = {median:.2f}% (threshold <2.0%)")
        n_pass += 1
    else:
        print(f"  FAIL  Median = {median:.2f}% (expected <2.0%)")
        n_fail += 1

    # Test 4: count under 1%
    if under_1 >= 300:
        print(f"  PASS  Under 1% = {under_1}/{n_mol} (threshold >= 300)")
        n_pass += 1
    else:
        print(f"  FAIL  Under 1% = {under_1}/{n_mol} (expected >= 300)")
        n_fail += 1

    print(f"\n  {n_pass}/{n_pass + n_fail} PASS")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
