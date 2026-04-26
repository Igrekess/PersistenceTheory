#!/usr/bin/env python3
"""
test_BSM_taxonomy.py -- Chapter 20d: BSM scenario taxonomy A/B/C

Monograph: chapters/ch20d_BSM_taxonomy.tex
Type: META — taxonomy of BSM extensions relative to PT canonical core.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


# (name, class, reason)
SCENARIOS = [
    ("4th generation", "A", "breaks q_+/q_- bifurcation"),
    ("Majorana mass at low energy", "A", "violates P1; breaks bifurcation"),
    ("Real DM scalar singlet (p=2)", "B", "compatible; +1/2 unit"),
    ("Axion (axial U(1) scalar)", "B", "+1/2 unit, compatible"),
    ("2HDM", "B", "+2 units, basin limit"),
    ("Right-handed Dirac neutrinos", "B", "Dirac mass; P1 preserved"),
    ("Full SUSY", "C", "Delta_mu0 ~ 12; exits basin"),
    ("3HDM", "C", "Delta_mu0 = 3; exits basin"),
    ("Three heavy Majorana partners", "C", "Delta_mu0 = 3; breaks bifurcation"),
    ("Light higgsino isolated", "C", "downgraded from B"),
]


ck = Checker("test_BSM_taxonomy", chapter="ch20d", total_steps=2)

# ---- Step 1: Taxonomy structure ----
ck.section("Step 1: A/B/C taxonomy of BSM scenarios")
n_A = sum(1 for s in SCENARIOS if s[1] == "A")
n_B = sum(1 for s in SCENARIOS if s[1] == "B")
n_C = sum(1 for s in SCENARIOS if s[1] == "C")
ck.check("class_A_nonempty", n_A >= 1, f"|Class A| = {n_A}")
ck.check("class_B_nonempty", n_B >= 1, f"|Class B| = {n_B}")
ck.check("class_C_nonempty", n_C >= 1, f"|Class C| = {n_C}")
ck.check("classes_partition", n_A + n_B + n_C == len(SCENARIOS),
         f"A+B+C = {n_A + n_B + n_C} = total {len(SCENARIOS)}")

# ---- Step 2: Specific exclusions / inclusions ----
ck.section("Step 2: structural classification consistency")
for cls in ["A", "B", "C"]:
    items = [s[0] for s in SCENARIOS if s[1] == cls]
    print(f"  Class {cls}: {items}")

# Specific structural claims from ch20d
dm_scalar = next(s for s in SCENARIOS if "DM scalar singlet" in s[0])
ck.check("DM_scalar_in_class_B", dm_scalar[1] == "B",
         "p=2 scalar singlet is canonical Class B")

susy = next(s for s in SCENARIOS if s[0] == "Full SUSY")
ck.check("full_SUSY_in_class_C", susy[1] == "C",
         "full SUSY exits the basin")

majorana = next(s for s in SCENARIOS if "Majorana" in s[0] and "low" in s[0])
ck.check("low_E_Majorana_in_class_A", majorana[1] == "A",
         "low-energy Majorana excluded by PT (P1)")

ck.summary()
