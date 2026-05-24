#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
V_td: tree vs NLO vs NNLO, compared to PDG unitarity fit and direct measurement.

Context (ch15 footnote on V_td, April 2026):
  PDG fit (from unitarity constraint)  : V_td^fit    = 0.00800
  Direct measurements (B->X_d gamma, B-Bbar oscillation) : V_td^direct = 0.00860(20)
  Pull between fit and direct         : ~3 sigma

Goal: quantify to what extent the PT NNLO correction
  V_td^NNLO = V_td^NLO * (1 + gamma_11 * alpha_EM)
closes this gap.

Finding: the NNLO ghost correction gamma_11 * alpha_EM = 0.31% is ~24x
too small to close the 7.5% gap between fit and direct. The tension
is inter-experimental, not a missing PT correction. Notably, the
PT *tree-level* Wolfenstein value V_td^tree = 0.00854 already agrees
with the direct measurement to 0.65%, while the NLO correction
(1 - n_f * eps) = 0.933 brings it in line with the unitarity fit.
"""
from __future__ import division, print_function

import os
import sys
import pathlib

SCRIPTS_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPTS_ROOT))

from pt_constants import (
    V_td, _V_td_exact, n_f, eps, _gamma_ghost, alpha_EM, alpha_s,
)


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


# =======================================================================
section("V_td: tree vs NLO vs NNLO")
# =======================================================================
V_td_tree = _V_td_exact
nlo_factor = 1 - n_f * eps
V_td_NLO   = V_td_tree * nlo_factor
nnlo_factor = 1 + _gamma_ghost[11] * alpha_EM
V_td_NNLO  = V_td_NLO * nnlo_factor

print(f"  V_td tree (Wolfenstein)         = {V_td_tree:.8f}")
print(f"  NLO factor (1 - n_f*eps)        = {nlo_factor:.8f}")
print(f"  V_td NLO                        = {V_td_NLO:.8f}")
print(f"  NNLO factor (1 + g11*alpha)     = {nnlo_factor:.8f}")
print(f"  V_td NNLO (current PT)          = {V_td_NNLO:.8f}")
print(f"  V_td (from pt_constants export) = {V_td:.8f}")

# =======================================================================
section("Experimental values and errors")
# =======================================================================
V_td_fit    = 0.00800   # PDG 2024 unitarity fit
V_td_direct = 0.00860   # B->X_d gamma + B-Bbar oscillation

def pct(a, b):
    return abs(a - b) / b * 100

print(f"  PDG unitarity fit            V_td^fit    = {V_td_fit}")
print(f"  Direct B->X_d gamma          V_td^direct = {V_td_direct}")
print()
print(f"  Relative errors (PT vs ...):")
print(f"    tree vs fit    = {pct(V_td_tree, V_td_fit):.3f}%")
print(f"    tree vs direct = {pct(V_td_tree, V_td_direct):.3f}%  <- EXCELLENT agreement")
print(f"    NLO  vs fit    = {pct(V_td_NLO,  V_td_fit):.3f}%")
print(f"    NLO  vs direct = {pct(V_td_NLO,  V_td_direct):.3f}%")
print(f"    NNLO vs fit    = {pct(V_td_NNLO, V_td_fit):.3f}%  <- current PT-vs-fit")
print(f"    NNLO vs direct = {pct(V_td_NNLO, V_td_direct):.3f}%  <- current PT-vs-direct")

# =======================================================================
section("NNLO correction size vs required gap")
# =======================================================================
required_shift = (V_td_direct - V_td_fit) / V_td_fit
nnlo_shift     = nnlo_factor - 1

print(f"  Required shift fit -> direct : {required_shift*100:.3f}% (7.5%)")
print(f"  NNLO shift gamma_11*alpha    : {nnlo_shift*100:.3f}% (0.31%)")
print(f"  Ratio required/available     : {required_shift/nnlo_shift:.1f}x too small")

# =======================================================================
section("CONCLUSION")
# =======================================================================
print("""
  1. PT tree-level Wolfenstein V_td = 0.00854 agrees with the direct
     measurement 0.0086(2) at 0.65%.

  2. The NLO correction (1 - n_f*eps) = 0.933 shifts V_td to 0.00797,
     aligning with the PDG unitarity fit 0.008 but moving away from
     the direct measurement (new error = 7%).

  3. The NNLO echo correction gamma_11 * alpha_EM = +0.31% is too small
     (by a factor ~24) to close the direct-vs-fit gap.

  4. Verdict: the ~3 sigma V_td tension is inter-experimental
     (direct vs unitarity fit), NOT a missing PT correction. If future
     direct measurements (Belle II, LHCb 2026-2028) confirm 0.0086,
     the PT NLO coefficient -n_f for V_td should be re-examined,
     since tree-level already agrees with direct.
""")

sys.exit(0)
