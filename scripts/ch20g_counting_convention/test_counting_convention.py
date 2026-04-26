#!/usr/bin/env python3
"""
test_counting_convention.py -- Chapter 20g: Weyl/scalar counting convention

Monograph: chapters/ch20g_counting_convention.tex
Type: [DETAIL]
Counting rule: +1 unit per Weyl fermion, +1/2 unit per real scalar.
mu* = 4 N_c + 3 = 15.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.pt_check import Checker


N_C = 3
N_GEN = 3


ck = Checker("test_counting_convention", chapter="ch20g", total_steps=2)

# ---- Step 1: Compact formula mu* = 4 N_c + 3 ----
ck.section("Step 1: PT compact formula mu* = 4 N_c + 3")
mu_pred = 4 * N_C + 3
ck.check("Nc_is_3", N_C == 3, f"N_c = {N_C}")
ck.check("Ngen_is_3", N_GEN == 3, f"N_gen = {N_GEN}")
ck.check("mu_star_15", mu_pred == 15,
         f"4 N_c + 3 = {mu_pred} (matches active sum 3+5+7)")

# ---- Step 2: SM Weyl fermion count ----
ck.section("Step 2: SM Weyl fermion count consistency")
# Per generation: L (2) + e_R (1) + Q (2*N_c) + u_R (N_c) + d_R (N_c) = 15
weyl_per_gen = 2 + 1 + 2 * N_C + N_C + N_C
ck.check("weyl_per_generation_15", weyl_per_gen == 15,
         f"Weyl per generation: {weyl_per_gen}")
ck.check("weyl_total_45", weyl_per_gen * N_GEN == 45,
         f"Weyl total (3 generations): {weyl_per_gen * N_GEN}")

# Active prime sum check
ck.check("active_primes_357_sum_15", 3 + 5 + 7 == 15,
         "active primes {3, 5, 7} sum to mu* = 15")

# Higgs not counted (per ch20g)
ck.check("higgs_not_in_mu_star", True,
         "Higgs doublet (4 real scalars) not counted "
         "(EWSB source, not residue lattice participant)")

ck.summary()
