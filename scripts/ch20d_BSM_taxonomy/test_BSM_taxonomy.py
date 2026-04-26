#!/usr/bin/env python3
"""
test_BSM_taxonomy.py -- Chapter 20d: BSM scenario taxonomy A/B/C

Monograph: chapters/ch20d_BSM_taxonomy.tex
Type: META — taxonomy of BSM extensions relative to PT canonical core.

This script enumerates the three classes:
  Class A: scenarios INCOMPATIBLE with PT (excluded by canonical structure)
  Class B: scenarios PARTIALLY COMPATIBLE (within basin Delta_mu0 <= 2)
  Class C: scenarios that REQUIRE PT REVISION (outside basin)

The classification is structural, not numerical: each scenario is
placed in A/B/C based on whether it is compatible with the canonical
{3,5,7} structure and the mu* = 15 fixed point.
"""


# ---------------------------------------------------------------
# BSM scenarios catalogued in ch20d
# ---------------------------------------------------------------
SCENARIOS = [
    # (name, class, reason)
    ("4th generation",                   "A", "breaks q_+/q_- bifurcation; leptons no longer exclusively in vertex branch"),
    ("Majorana mass at low energy",      "A", "violates P1 (no Majorana phase); breaks bifurcation symmetry"),
    ("Real DM scalar singlet (p=2)",     "B", "compatible with basin; 1 real scalar = +1/2 unit (Delta_mu0 = 1/2)"),
    ("Axion (axial U(1) scalar)",        "B", "compatible; 1 real scalar (Delta_mu0 = 1/2)"),
    ("2 Higgs doublet model (2HDM)",     "B", "compatible; 4 real scalars + 4 complex (Delta_mu0 = 2 limit)"),
    ("Right-handed Dirac neutrinos",     "B", "Dirac mass mechanism (no Majorana); P1 preserved"),
    ("Full SUSY (squarks, gauginos, ..)", "C", "Delta_mu0 ~ 12 fermion + scalar partners; exits basin (Delta_mu0 > 2)"),
    ("3HDM",                             "C", "Delta_mu0 = 3 (3 doublets); exits basin"),
    ("Three heavy Majorana partners",    "C", "Delta_mu0 = 3 fermions; exits basin AND breaks bifurcation"),
    ("Light higgsino isolated",          "C", "Downgrade from B in this monograph (ch20e_basin_robustness.tex)"),
]


def main():
    print("=" * 70)
    print("Chapter 20d: BSM scenario taxonomy A/B/C")
    print("=" * 70)
    print()
    print("Class A: incompatible with PT canonical structure")
    print("Class B: compatible within basin Delta_mu0 <= 2 (preserves mu* = 15)")
    print("Class C: requires PT revision (outside basin or breaks bifurcation)")
    print()

    for cls in ["A", "B", "C"]:
        items = [s for s in SCENARIOS if s[1] == cls]
        print(f"--- Class {cls} ({len(items)} scenarios) ---")
        for name, _, reason in items:
            print(f"  - {name}")
            print(f"      reason: {reason}")
        print()

    # Summary
    print(f"Total scenarios catalogued: {len(SCENARIOS)}")
    print(f"  Class A (excluded):        {sum(1 for s in SCENARIOS if s[1] == 'A')}")
    print(f"  Class B (compatible):       {sum(1 for s in SCENARIOS if s[1] == 'B')}")
    print(f"  Class C (requires revision):{sum(1 for s in SCENARIOS if s[1] == 'C')}")
    print()
    print("This taxonomy is qualitative; the basin radius Delta_mu0 <= 2 ")
    print("is verified numerically in ch20e_basin_robustness.")


if __name__ == "__main__":
    main()
