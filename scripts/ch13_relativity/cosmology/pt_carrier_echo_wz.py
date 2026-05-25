#!/usr/bin/env python3
"""
PT Route D: Carrier-Echo Mechanism for Dark Energy
====================================================
w(z) = -1 + s * exp(-(1+z))   with s = 1/2, mu*(z) = mu*/(1+z)

Ingredients:
  - s = 1/2        [T1, Forbidden Transitions, PT_MONOGRAPHY ch.AppO]
  - echo product = 1/e  [L2 identity, ch21b_RH_echo_product.tex, IDTAG]
  - mu_eff(z) = mu*/(1+z)  [COND: cosmological persistence scale]

DESI DR2 (2025) reference:
  w0_DESI = -0.80 +/- 0.06    (paper: w0 = -0.75 +/- 0.067, here we use -0.80)
  wa_DESI = -0.72 +/- 0.25

Note on DESI DR2 values: different analyses quote slightly different
central values. We use the DESI Collaboration 2025 combined BAO+CMB:
  w0 = -0.80 +/- 0.06,  wa = -0.72 +/- 0.25
"""

import numpy as np
import sys

# ─── PT constants ────────────────────────────────────────────────────────────
mu_star = 15.0          # cascade fixed point
s = 0.5                 # carrier amplitude p=2, T1 theorem
e = np.e

# ─── DESI DR2 reference values ───────────────────────────────────────────────
w0_DESI  = -0.80
s_w0     = 0.06
wa_DESI  = -0.72
s_wa     = 0.25

# Planck 2018 / best-fit cosmological parameters
Omega_m   = 0.315
Omega_Lam = 0.685        # PT prediction: 68.5%
q0_PT     = -0.530       # PT prediction for deceleration parameter

# ─── Step 1: Numerical verification of ingredients ───────────────────────────
print("=" * 65)
print("PT Route D — Carrier-Echo Dark Energy")
print("=" * 65)

print("\n--- Step 1: Numerical verification of ingredients ---\n")

qminus_mustar = np.exp(-1.0 / mu_star)
print(f"  q-(mu*) = exp(-1/15) = {qminus_mustar:.10f}")

echo_product = qminus_mustar**3 * qminus_mustar**5 * qminus_mustar**7
echo_product_exact = np.exp(-15.0 / mu_star)
print(f"  Echo product = q-^3 * q-^5 * q-^7 = exp(-15/15) = exp(-1)")
print(f"    Numerical:  {echo_product:.16f}")
print(f"    1/e exact:  {1.0/e:.16f}")
print(f"    |diff|:     {abs(echo_product - 1.0/e):.2e}  (machine precision)")

# w0 computation
w0_RouteD = -1.0 + s * np.exp(-1.0)
print(f"\n  s = 1/2 = {s}")
print(f"  w0 = -1 + s/e = -1 + {s}/{e:.6f} = {w0_RouteD:.6f}")

# Tension with DESI DR2
pull_w0 = (w0_RouteD - w0_DESI) / s_w0
print(f"\n  DESI DR2: w0 = {w0_DESI} +/- {s_w0}")
print(f"  Pull(w0) = ({w0_RouteD:.4f} - ({w0_DESI})) / {s_w0} = {pull_w0:+.2f} sigma")

# ─── Step 2: w(z) trajectory and CPL fit ─────────────────────────────────────
print("\n--- Step 2: w(z) trajectory and CPL fit ---\n")

def w_RouteD(z, s=0.5):
    """Route D: w(z) = -1 + s * exp(-(1+z))"""
    return -1.0 + s * np.exp(-(1.0 + z))

def w_RouteD_a(a, s=0.5):
    """Route D in terms of scale factor a=1/(1+z): w(a) = -1 + s * exp(-1/a)"""
    return -1.0 + s * np.exp(-1.0 / a)

# w(z) at specific redshifts
zvals = [0.0, 0.5, 1.0, 2.0, 3.0]
print("  z         a=1/(1+z)   w(z)")
print("  " + "-" * 40)
for z in zvals:
    a = 1.0 / (1.0 + z)
    wz = w_RouteD(z)
    print(f"  z={z:.1f}      a={a:.4f}      w={wz:.6f}")

# CPL approximation: w(a) = w0_CPL + wa_CPL * (1 - a)
# dw/da = s * exp(-1/a) / a^2
# At a=1: dw/da = s/e
# CPL convention: w(a) = w0 + wa*(1-a) -> dw/da = -wa
# => wa_CPL = -dw/da|_{a=1} = -s/e

w0_CPL = w_RouteD_a(1.0)   # = -1 + s/e
dw_da_at_1 = s * np.exp(-1.0) / (1.0**2)
wa_CPL = -dw_da_at_1        # = -s/e

print(f"\n  CPL linearization w(a) = w0_CPL + wa_CPL*(1-a):")
print(f"    dw/da|_{{a=1}} = s * e^{{-1/a}} / a^2 |_{{a=1}} = s/e = {dw_da_at_1:.6f}")
print(f"    w0_CPL = w(a=1) = -1 + s/e = {w0_CPL:.6f}")
print(f"    wa_CPL = -dw/da|_{{a=1}} = -s/e = {wa_CPL:.6f}")

# Comparison with DESI DR2
pull_wa = (wa_CPL - wa_DESI) / s_wa
print(f"\n  DESI DR2: wa = {wa_DESI} +/- {s_wa}")
print(f"  Pull(wa) = ({wa_CPL:.4f} - ({wa_DESI})) / {s_wa} = {pull_wa:+.2f} sigma")

# ─── Step 3: Constraint w0 + wa ──────────────────────────────────────────────
print("\n--- Step 3: Constraint w0_CPL + wa_CPL ---\n")

sum_w = w0_CPL + wa_CPL
print(f"  w0_CPL + wa_CPL = ({w0_CPL:.6f}) + ({wa_CPL:.6f}) = {sum_w:.10f}")
print(f"  Exact: (-1 + s/e) + (-s/e) = -1.  |diff from -1| = {abs(sum_w + 1):.2e}")
print(f"\n  DESI DR2: w0+wa = {w0_DESI} + ({wa_DESI}) = {w0_DESI+wa_DESI:.2f}")
print(f"  PT Route D predicts w0+wa = -1  (FALSIFIABLE PREDICTION)")
print(f"  Difference from DESI central: {sum_w - (w0_DESI+wa_DESI):.3f}")

# ─── Step 4: Deceleration parameter consistency ──────────────────────────────
print("\n--- Step 4: Deceleration parameter q0 ---\n")

# q0 = Omega_m/2 + Omega_Lam*(1 + 3*w0)/2
w0_for_q0 = w0_CPL   # Route D w0
q0_computed = Omega_m/2.0 + Omega_Lam*(1.0 + 3.0*w0_for_q0)/2.0
print(f"  q0 formula: Omega_m/2 + Omega_Lam*(1 + 3*w0)/2")
print(f"  With Omega_m={Omega_m}, Omega_Lam={Omega_Lam}, w0={w0_for_q0:.6f}:")
print(f"    q0 = {Omega_m/2.0:.4f} + {Omega_Lam:.3f}*(1 + 3*({w0_for_q0:.6f}))/2")
print(f"       = {Omega_m/2.0:.4f} + {Omega_Lam:.3f}*({1.0+3.0*w0_for_q0:.6f})/2")
print(f"       = {Omega_m/2.0:.4f} + {Omega_Lam*(1.0+3.0*w0_for_q0)/2.0:.4f}")
print(f"       = {q0_computed:.4f}")

print(f"\n  PT prediction: q0 = {q0_PT}")
print(f"  Route D gives: q0 = {q0_computed:.4f}")
pull_q0 = q0_computed - q0_PT
print(f"  TENSION: Delta_q0 = {pull_q0:.4f}  (internal inconsistency)")

# Find Omega_Lam that recovers q0 = -0.530 with w0 = w0_RouteD
# q0 = Omega_m/2 + Omega_Lam*(1 + 3*w0)/2
# -0.530 = Omega_m/2 + Omega_Lam*(1 + 3*w0)/2
# With Omega_m = 1 - Omega_Lam (flat universe):
# -0.530 = (1-OL)/2 + OL*(1+3*w0)/2
# -1.060 = (1-OL) + OL*(1+3*w0)
# -1.060 = 1 - OL + OL + 3*OL*w0
# -1.060 = 1 + 3*OL*w0
# 3*OL*w0 = -2.060
# OL = -2.060 / (3*w0)

w0_rd = w0_RouteD
OL_required = (-q0_PT - Omega_m/2.0 + 0.5) / (-(1.0 + 3.0*w0_rd)/2.0)
# More careful:
# q0 = (1-OL)/2 + OL*(1+3w0)/2 = 0.5 - OL/2 + OL*(1+3w0)/2
# = 0.5 + OL*[-1/2 + (1+3w0)/2]
# = 0.5 + OL * [3w0/2]
# => OL = (q0 - 0.5) / (3*w0/2)
OL_required = (q0_PT - 0.5) / (3.0 * w0_rd / 2.0)
Omm_required = 1.0 - OL_required

print(f"\n  To recover q0={q0_PT} with w0={w0_rd:.4f} (flat universe):")
print(f"    Omega_Lambda required = {OL_required:.4f}")
print(f"    Omega_m required      = {Omm_required:.4f}")
print(f"  (Compare: Planck Omega_Lambda = {Omega_Lam}, Omega_m = {Omega_m})")
delta_OL = OL_required - Omega_Lam
print(f"  Delta_Omega_Lambda = {delta_OL:+.4f}")
print(f"\n  INTERNAL TENSION: Route D w0=-0.816 is incompatible with")
print(f"  PT q0=-0.530 and Planck Omega_m=0.315 simultaneously.")
print(f"  Either Omega_Lambda must shift to {OL_required:.3f},")
print(f"  or q0 is NOT constrained to -0.530 in Route D.")

# ─── Step 5: Summary table ───────────────────────────────────────────────────
print("\n--- Summary: Route D vs DESI DR2 ---\n")
print(f"  {'Observable':<25} {'Route D':>12} {'DESI DR2':>12} {'Pull':>10}")
print("  " + "-" * 62)
print(f"  {'w0':<25} {w0_CPL:>12.4f} {w0_DESI:>12.4f} {pull_w0:>+10.2f} sigma")
print(f"  {'wa':<25} {wa_CPL:>12.4f} {wa_DESI:>12.4f} {pull_wa:>+10.2f} sigma")
print(f"  {'w0 + wa':<25} {sum_w:>12.4f} {w0_DESI+wa_DESI:>12.4f} {'(PT=-1 exact)':>10}")
print(f"  {'q0':<25} {q0_computed:>12.4f} {q0_PT:>12.4f} {'internal':>10}")

# ─── Step 6: w(z) values for LaTeX table ─────────────────────────────────────
print("\n--- w(z) values for LaTeX table ---\n")
print(f"  {'z':>6}  {'a':>8}  {'w_RouteD(z)':>14}  {'w_CPL(z)':>14}")
for z in [0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]:
    a = 1/(1+z)
    w_exact = w_RouteD(z)
    w_cpl   = w0_CPL + wa_CPL*(1 - a)
    print(f"  {z:>6.1f}  {a:>8.4f}  {w_exact:>14.6f}  {w_cpl:>14.6f}")

# ─── Step 7: Exact analytical summary ────────────────────────────────────────
print("\n--- Analytical summary ---\n")
print(f"  w(z) = -1 + (1/2) * exp(-(1+z))")
print(f"  w(a) = -1 + (1/2) * exp(-1/a)     [scale factor form]")
print(f"")
print(f"  CPL parameters (tangent at a=1):")
print(f"    w0 = -1 + 1/(2e) = {-1.0 + 1/(2*e):.6f}")
print(f"    wa = -1/(2e)     = {-1.0/(2*e):.6f}")
print(f"    w0 + wa = -1     [exact, falsifiable prediction]")
print(f"")
print(f"  Echo product identity (ch21b, IDTAG):")
print(f"    prod_{{p in {{3,5,7}}}} q-^p = exp(-15/15) = 1/e  [EXACT]")
print(f"")
print(f"  Status epistémique:")
print(f"    w0 = -1 + s/e          [DER-PARTIAL]: s=T1 (THM), echo=IDTAG")
print(f"    wa = -s/e              [COND]: nécessite mu_eff(z)=mu*/(1+z)")
print(f"    w0+wa=-1               [PRED-FALSIFIABLE]: vérifiable par DESI DR3")

print("\n" + "=" * 65)
print("Script completed successfully.")
print("=" * 65)
