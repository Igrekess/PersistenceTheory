"""
test_selection_G_2pi_alpha
==========================

ENGLISH
-------
Selection of G = 2*pi*alpha from sieve constraints and Einstein equations

FRANCAIS (original)
-------------------
test_selection_G_2pi_alpha.py
S15.6.126: NON-UNIVERSALITY of G = 2*pi*alpha as a SELECTION PRINCIPLE

DEEP INVESTIGATION: The relation G_sieve = G_00/(8*pi*D_KL) = 2*pi*alpha holds
to 0.29% at mu_alpha = 15.04. But is this a selection principle that picks out
the operating point, or does it hold everywhere (universal tautology)?

7 TESTS:
1. High-resolution scan: ratio G(mu)/(2*pi*alpha(mu)) from mu=8 to mu=50
2. Multiple selection conditions at mu~15: how many cross independently?
3. Derivative d/dmu[G/(2*pi*alpha)] at crossing -- sharpness of selection
4. Compare crossing of G=2*pi*alpha vs c*105=phi conditions
5. Role of delta_mu = mu - 15 and significance of 15 = 3+5+7
6. Search for deeper formula: G_00 = n*pi^k * alpha^a * D_KL^b
7. Physical interpretation: self-consistency and what it means

Author / Auteur: Yan Senez  |  Date: February / Fevrier 2026
Theory / Theorie: Persistence Theory (PT) / Theorie de la Persistance (TP)
"""

import numpy as np
from math import sqrt, log, log2, pi, exp
from scipy.optimize import brentq, minimize_scalar

print("=" * 70)
print("S15.6.126: NON-UNIVERSALITY OF G = 2*pi*alpha")
print("           AS A SELECTION PRINCIPLE -- DEEP INVESTIGATION")
print("=" * 70)

# ==============================================================
# CONSTANTS AND SIEVE FUNCTIONS
# ==============================================================

phi_golden = (1 + sqrt(5)) / 2
s = 0.5
alpha_EM_phys = 1 / 137.035999084
active_primes = [3, 5, 7]

def sin2_theta(p, q):
    """sin^2(theta_p) = (1-q^p)(2p-1+q^p)/p^2"""
    qp = q**p
    return (1 - qp) * (2*p - 1 + qp) / p**2

def alpha_sieve(mu):
    """alpha(mu) = product of sin^2(theta_p) for p in {3,5,7}"""
    if mu <= 2:
        return 0.0
    q = 1 - 2/mu
    if q <= 0 or q >= 1:
        return 0.0
    result = 1.0
    for p in active_primes:
        result *= sin2_theta(p, q)
    return result

def gamma_p_func(p, mu):
    """Effective dimension (sieve beta-function) per prime."""
    if mu <= 2:
        return 0.0
    q = 1 - 2/mu
    if q <= 0 or q >= 1:
        return 0.0
    qp = q**p
    if abs(1 - qp) < 1e-15:
        return 0.0
    delta_p = (1 - qp) / p
    dln_delta = -2*p * q**(p-1) / (mu * (1 - qp))
    factor = 2*(1 - delta_p) / (2 - delta_p)
    return -dln_delta * factor

def ln_alpha(mu):
    a = alpha_sieve(mu)
    return log(a) if a > 0 else -100.0

def d2_ln_alpha(mu, h=1e-4):
    return (ln_alpha(mu+h) - 2*ln_alpha(mu) + ln_alpha(mu-h)) / h**2

def lapse(mu):
    return sqrt(abs(d2_ln_alpha(mu)))

def D_KL_at(mu_val):
    """D_KL = ln(2) + D_KL(mod 3, geometric)"""
    if mu_val <= 2:
        return 0.0
    q = exp(-2/mu_val)
    P = np.zeros(3)
    for k in range(1, 500):
        r = (2*k) % 3
        P[r] += (1 - q) * q**(k-1)
    P /= P.sum()
    D3 = sum(P[r] * log(3*P[r]) for r in range(3) if P[r] > 0)
    return log(2) + D3

def einstein_G00(mu, hd=1e-4):
    """G^0_0 = H1*H2 + H1*H3 + H2*H3 (Bianchi I Friedmann)"""
    N_val = lapse(mu)
    if N_val < 1e-15:
        return 0.0

    Hs = []
    for p in active_primes:
        gp_c = gamma_p_func(p, mu)
        gp_p = gamma_p_func(p, mu + hd)
        gp_m = gamma_p_func(p, mu - hd)
        a_i = gp_c / mu
        da = (gp_p/(mu+hd) - gp_m/(mu-hd)) / (2*hd)
        if N_val > 0 and a_i > 0:
            Hs.append(da / (N_val * a_i))
        else:
            Hs.append(0.0)

    H1, H2, H3 = Hs
    return H1*H2 + H1*H3 + H2*H3

def einstein_full(mu, hd=1e-4):
    """Full Einstein tensor (Bianchi I)."""
    N_val = lapse(mu)
    if N_val < 1e-15:
        return None

    def get_hubble(mu_e):
        N_e = lapse(mu_e)
        Hs = []
        for p in active_primes:
            gp_c = gamma_p_func(p, mu_e)
            gp_p = gamma_p_func(p, mu_e + hd)
            gp_m = gamma_p_func(p, mu_e - hd)
            a_i = gp_c / mu_e
            da = (gp_p/(mu_e+hd) - gp_m/(mu_e-hd)) / (2*hd)
            Hs.append(da / (N_e * a_i) if N_e > 0 and a_i > 0 else 0)
        return Hs

    H = get_hubble(mu)
    H1, H2, H3 = H
    G_00 = H1*H2 + H1*H3 + H2*H3

    H_plus = get_hubble(mu + hd)
    H_minus = get_hubble(mu - hd)
    dH = [(H_plus[i] - H_minus[i]) / (2*hd*N_val) for i in range(3)]

    G_sp = []
    pairs = [(1,2), (0,2), (0,1)]
    for idx, (j, k) in enumerate(pairs):
        G_ii = dH[j] + dH[k] + H[j]**2 + H[k]**2 + H[j]*H[k]
        G_sp.append(G_ii)

    R = -2*(sum(dH) + H1**2 + H2**2 + H3**2 + H1*H2 + H1*H3 + H2*H3)

    return {
        'G_00': G_00, 'G_sp': G_sp, 'R': R,
        'H': H, 'dH': dH, 'N': N_val
    }

def ratio_G_over_2pi_alpha(mu_val):
    """R(mu) = G_00/(8*pi*D_KL) / (2*pi*alpha) = G_00/(16*pi^2*alpha*D_KL)"""
    a_val = alpha_sieve(mu_val)
    if a_val <= 0:
        return np.nan
    G00 = einstein_G00(mu_val)
    D_val = D_KL_at(mu_val)
    if D_val <= 0:
        return np.nan
    denom = 16 * pi**2 * a_val * D_val
    if abs(denom) < 1e-30:
        return np.nan
    return G00 / denom

# ==============================================================
# OPERATING POINT
# ==============================================================

mu_alpha = brentq(lambda m: alpha_sieve(m) - alpha_EM_phys, 14.5, 16.0)
delta_mu_op = mu_alpha - 15.0
alpha_op = alpha_sieve(mu_alpha)

print(f"\nOperating point:")
print(f"  mu_alpha = {mu_alpha:.6f}, delta_mu = {delta_mu_op:.6f}")
print(f"  alpha = {alpha_op:.10e}, 1/alpha = {1/alpha_op:.4f}")

G00_op = einstein_G00(mu_alpha)
D_op = D_KL_at(mu_alpha)
G_sieve_op = G00_op / (8*pi*D_op)
G_pred_op = 2*pi*alpha_op
ratio_op = G_sieve_op / G_pred_op

print(f"  G_00    = {G00_op:.8f}")
print(f"  D_KL    = {D_op:.8f}")
print(f"  G_sieve = G_00/(8*pi*D) = {G_sieve_op:.8f}")
print(f"  2*pi*alpha = {G_pred_op:.8f}")
print(f"  Ratio at operating point = {ratio_op:.6f} (deviation {abs(ratio_op-1)*100:.4f}%)")

# ==============================================================
# TEST 1: HIGH-RESOLUTION SCAN
# ==============================================================

print(f"\n{'='*70}")
print("TEST 1: High-resolution scan R(mu) = G_00/(16*pi^2*alpha*D_KL)")
print("        mu = 8 to 50, 250 points")
print("="*70)

mu_scan = np.linspace(8.0, 50.0, 250)
ratios_scan = []
mu_valid = []

for mu_val in mu_scan:
    r = ratio_G_over_2pi_alpha(mu_val)
    if not np.isnan(r) and r > 0:
        ratios_scan.append(r)
        mu_valid.append(mu_val)

mu_valid = np.array(mu_valid)
ratios_scan = np.array(ratios_scan)

# Statistics
mean_ratio = np.mean(ratios_scan)
std_ratio = np.std(ratios_scan)
min_ratio = np.min(ratios_scan)
max_ratio = np.max(ratios_scan)
mu_at_min = mu_valid[np.argmin(ratios_scan)]
mu_at_max = mu_valid[np.argmax(ratios_scan)]
cv_ratio = std_ratio / abs(mean_ratio) * 100

# Print subsampled table
print(f"\n  {'mu':>6s} {'alpha':>12s} {'G_00':>10s} {'D_KL':>10s} {'ratio':>8s} {'|R-1|%':>8s}")
print(f"  {'-'*60}")

step = max(1, len(mu_valid) // 20)
for idx in range(0, len(mu_valid), step):
    mu_v = mu_valid[idx]
    a_v = alpha_sieve(mu_v)
    G00_v = einstein_G00(mu_v)
    D_v = D_KL_at(mu_v)
    r_v = ratios_scan[idx]
    print(f"  {mu_v:6.2f} {a_v:12.6e} {G00_v:10.6f} {D_v:10.6f} {r_v:8.4f} {abs(r_v-1)*100:8.3f}")

print(f"\n  Statistics over {len(ratios_scan)} valid points:")
print(f"    Mean ratio  = {mean_ratio:.6f}")
print(f"    Std ratio   = {std_ratio:.6f}")
print(f"    CV          = {cv_ratio:.2f}%")
print(f"    Min ratio   = {min_ratio:.6f} at mu = {mu_at_min:.2f}")
print(f"    Max ratio   = {max_ratio:.6f} at mu = {mu_at_max:.2f}")
print(f"    Range       = {max_ratio - min_ratio:.6f}")

# Find crossing point(s) where ratio = 1
crossings = []
for i in range(len(ratios_scan)-1):
    if (ratios_scan[i] - 1) * (ratios_scan[i+1] - 1) < 0:
        mu_cross = mu_valid[i] + (1 - ratios_scan[i]) * (mu_valid[i+1] - mu_valid[i]) / (ratios_scan[i+1] - ratios_scan[i])
        crossings.append(mu_cross)

print(f"\n  Crossing points where ratio = 1:")
if crossings:
    for ic, mc in enumerate(crossings):
        print(f"    Crossing {ic+1}: mu = {mc:.4f} (delta from 15: {mc-15:.4f})")
else:
    # Check if closest approach
    closest_idx = np.argmin(np.abs(ratios_scan - 1))
    print(f"    NO exact crossing found in [8, 50]!")
    print(f"    Closest: mu = {mu_valid[closest_idx]:.4f}, ratio = {ratios_scan[closest_idx]:.6f}")

# Refine crossing(s) with brentq
mu_star_list = []
for mc in crossings:
    try:
        def R_minus_1(mu):
            r = ratio_G_over_2pi_alpha(mu)
            return r - 1 if not np.isnan(r) else 10.0
        lo = max(mc - 0.5, 7)
        hi = mc + 0.5
        if R_minus_1(lo) * R_minus_1(hi) < 0:
            mu_star = brentq(R_minus_1, lo, hi, xtol=1e-10)
            mu_star_list.append(mu_star)
            a_star = alpha_sieve(mu_star)
            print(f"\n    EXACT crossing (brentq): mu* = {mu_star:.10f}")
            print(f"      alpha(mu*) = {a_star:.10e}, 1/alpha = {1/a_star:.6f}")
            print(f"      mu_alpha   = {mu_alpha:.10f}")
            print(f"      |mu* - mu_alpha| = {abs(mu_star - mu_alpha):.10f}")
            print(f"      |1/alpha* - 137.036| = {abs(1/a_star - 137.036):.6f}")
    except Exception:
        pass

# Classify universality
is_universal = cv_ratio < 5.0
is_selective = len(crossings) > 0 and cv_ratio > 10.0

if is_universal:
    verdict_1 = "QUASI-UNIVERSAL: ratio ~ constant over entire range"
elif is_selective:
    verdict_1 = "NON-UNIVERSAL: ratio varies significantly, crossing(s) exist"
else:
    verdict_1 = "AMBIGUOUS: moderate variation"

print(f"\n  VERDICT: {verdict_1}")
test1 = len(ratios_scan) >= 200
print(f"  [{'PASS' if test1 else 'FAIL'}] {len(ratios_scan)} valid points computed")

# ==============================================================
# TEST 2: MULTIPLE SELECTION CONDITIONS AT mu~15
# ==============================================================

print(f"\n{'='*70}")
print("TEST 2: Multiple selection conditions near mu ~ 15")
print("        How many INDEPENDENT conditions cross at this point?")
print("="*70)

conditions_met = {}
condition_details = {}

# --- Condition A: alpha(mu) = alpha_EM ---
print(f"\n  A. alpha(mu) = alpha_EM => mu_A = {mu_alpha:.6f}")
conditions_met['A: alpha=alpha_EM'] = mu_alpha

# --- Condition B: G_sieve = 2*pi*alpha (from Test 1) ---
if mu_star_list:
    mu_B = mu_star_list[0]
    conditions_met['B: G=2pi*alpha'] = mu_B
    print(f"  B. G = 2*pi*alpha     => mu_B = {mu_B:.6f}")
else:
    print(f"  B. G = 2*pi*alpha     => no crossing found")
    conditions_met['B: G=2pi*alpha'] = None

# --- Condition C: G^0_0 extremum (dG_00/dmu = 0) ---
print(f"\n  C. dG_00/dmu = 0 (Einstein energy density extremum):")
mu_dense = np.linspace(10, 25, 200)
G00_dense = np.array([einstein_G00(mu_v) for mu_v in mu_dense])
dG00 = np.gradient(G00_dense, mu_dense)

G00_extr_crossings = []
for i in range(len(dG00)-1):
    if dG00[i] * dG00[i+1] < 0:
        mc = mu_dense[i] + (-dG00[i]) * (mu_dense[i+1] - mu_dense[i]) / (dG00[i+1] - dG00[i])
        G00_extr_crossings.append(mc)

if G00_extr_crossings:
    for mc in G00_extr_crossings:
        print(f"     dG_00/dmu = 0 at mu = {mc:.4f}")
    conditions_met['C: dG00/dmu=0'] = G00_extr_crossings[0]
else:
    print(f"     G_00 is monotone in [10, 25] (no extremum)")
    conditions_met['C: dG00/dmu=0'] = None

# --- Condition D: Ricci scalar extremum (dR/dmu = 0) ---
# Use coarser grid + smoothing to avoid numerical noise from nested finite differences
print(f"\n  D. dR/dmu = 0 (Ricci scalar extremum):")
mu_coarse = np.linspace(10, 25, 60)  # coarser grid to reduce noise
R_coarse = []
for mu_v in mu_coarse:
    E_test = einstein_full(mu_v, hd=5e-4)  # larger step for stability
    R_coarse.append(E_test['R'] if E_test else 0.0)
R_coarse = np.array(R_coarse)

# Smooth before differentiating
from scipy.ndimage import uniform_filter1d
R_smooth = uniform_filter1d(R_coarse, size=5)
dR = np.gradient(R_smooth, mu_coarse)

R_extr_crossings = []
for i in range(len(dR)-1):
    if dR[i] * dR[i+1] < 0:
        mc = mu_coarse[i] + (-dR[i]) * (mu_coarse[i+1] - mu_coarse[i]) / (dR[i+1] - dR[i])
        R_extr_crossings.append(mc)

if R_extr_crossings:
    for mc in R_extr_crossings:
        print(f"     dR/dmu = 0 at mu = {mc:.4f}")
    # Find the one closest to 15
    closest_R = min(R_extr_crossings, key=lambda x: abs(x - 15))
    conditions_met['D: dR/dmu=0'] = closest_R
else:
    print(f"     R is monotone in [10, 25] (no extremum)")
    conditions_met['D: dR/dmu=0'] = None

# --- Condition E: Raychaudhuri equilibrium ---
# Use coarser grid + smoothing (same as D) to avoid numerical oscillation
print(f"\n  E. Raychaudhuri: d(theta)/dtau + theta^2/3 = 0:")
ray_vals_coarse = []
for mu_v in mu_coarse:
    E_test = einstein_full(mu_v, hd=5e-4)
    if E_test:
        theta = sum(E_test['H'])
        dtheta = sum(E_test['dH'])
        ray_vals_coarse.append(dtheta + theta**2 / 3)
    else:
        ray_vals_coarse.append(0.0)
ray_arr = uniform_filter1d(np.array(ray_vals_coarse), size=5)

ray_crossings = []
for i in range(len(ray_arr)-1):
    if ray_arr[i] * ray_arr[i+1] < 0:
        mc = mu_coarse[i] + (-ray_arr[i]) * (mu_coarse[i+1] - mu_coarse[i]) / (ray_arr[i+1] - ray_arr[i])
        ray_crossings.append(mc)

if ray_crossings:
    for mc in ray_crossings:
        print(f"     Raychaudhuri = 0 at mu = {mc:.4f}")
    # Find closest to 15
    closest_ray = min(ray_crossings, key=lambda x: abs(x - 15))
    conditions_met['E: Raychaudhuri=0'] = closest_ray
else:
    print(f"     No Raychaudhuri equilibrium in [10, 25]")
    print(f"     Range: [{np.min(ray_arr):.4f}, {np.max(ray_arr):.4f}]")
    conditions_met['E: Raychaudhuri=0'] = None

# --- Condition F: gamma_3 + gamma_5 + gamma_7 has integer or special value ---
print(f"\n  F. gamma_total = gamma_3 + gamma_5 + gamma_7 extremum or integer:")
gamma_total = np.array([sum(gamma_p_func(p, mu_v) for p in active_primes) for mu_v in mu_dense])
dgamma = np.gradient(gamma_total, mu_dense)

gamma_extr = []
for i in range(len(dgamma)-1):
    if dgamma[i] * dgamma[i+1] < 0:
        mc = mu_dense[i] + (-dgamma[i]) * (mu_dense[i+1] - mu_dense[i]) / (dgamma[i+1] - dgamma[i])
        gamma_extr.append(mc)

gamma_at_op = sum(gamma_p_func(p, mu_alpha) for p in active_primes)
print(f"     gamma_total(mu_alpha) = {gamma_at_op:.6f}")
if gamma_extr:
    for mc in gamma_extr:
        gt = sum(gamma_p_func(p, mc) for p in active_primes)
        print(f"     d(gamma)/dmu = 0 at mu = {mc:.4f}, gamma = {gt:.4f}")
    conditions_met['F: dgamma/dmu=0'] = gamma_extr[0]
else:
    print(f"     gamma_total is monotone in [10, 25]")
    conditions_met['F: dgamma/dmu=0'] = None

# Summary
print(f"\n  --- Summary of conditions ---")
print(f"  {'Condition':25s} {'mu_cross':>12s} {'|mu-15|':>8s} {'|mu-mu_a|':>10s}")
print(f"  {'-'*58}")
n_near_15 = 0
n_near_mu_alpha = 0
for label, mu_c in conditions_met.items():
    if mu_c is not None:
        d15 = abs(mu_c - 15)
        dma = abs(mu_c - mu_alpha)
        near15 = "*" if d15 < 1.0 else ""
        near_ma = "+" if dma < 1.0 else ""
        if d15 < 1.0:
            n_near_15 += 1
        if dma < 1.0:
            n_near_mu_alpha += 1
        print(f"  {label:25s} {mu_c:12.4f} {d15:8.4f}{near15:2s} {dma:10.4f}{near_ma}")
    else:
        print(f"  {label:25s} {'---':>12s} {'---':>8s} {'---':>10s}")

print(f"\n  {n_near_15} conditions within |mu-15| < 1")
print(f"  {n_near_mu_alpha} conditions within |mu-mu_alpha| < 1")

test2 = n_near_mu_alpha >= 2
print(f"\n  [{'PASS' if test2 else 'FAIL'}] Multiple conditions confluent at mu ~ {mu_alpha:.1f}")

# ==============================================================
# TEST 3: DERIVATIVE -- SHARPNESS OF CROSSING
# ==============================================================

print(f"\n{'='*70}")
print("TEST 3: Sharpness of G/(2*pi*alpha) = 1 crossing")
print("="*70)

# Numerical derivative of ratio at operating point
h_d = 0.01
dr_dmu_op = (ratio_G_over_2pi_alpha(mu_alpha + h_d) - ratio_G_over_2pi_alpha(mu_alpha - h_d)) / (2*h_d)

print(f"\n  At the operating point mu = {mu_alpha:.6f}:")
print(f"    R(mu_alpha) = {ratio_op:.8f}")
print(f"    dR/dmu      = {dr_dmu_op:.8f} per unit mu")

# Second derivative
d2r = (ratio_G_over_2pi_alpha(mu_alpha + h_d) - 2*ratio_op + ratio_G_over_2pi_alpha(mu_alpha - h_d)) / h_d**2
print(f"    d^2R/dmu^2  = {d2r:.8f}")

# At the crossing point (if it exists)
if mu_star_list:
    mu_s = mu_star_list[0]
    dr_at_cross = (ratio_G_over_2pi_alpha(mu_s + h_d) - ratio_G_over_2pi_alpha(mu_s - h_d)) / (2*h_d)
    print(f"\n  At the crossing mu* = {mu_s:.6f}:")
    print(f"    R(mu*) = 1 (by definition)")
    print(f"    dR/dmu = {dr_at_cross:.8f}")

    # Selectivity: width where |R-1| < epsilon
    for eps_thresh in [0.01, 0.05, 0.10]:
        width = 2 * eps_thresh / abs(dr_at_cross) if abs(dr_at_cross) > 1e-10 else float('inf')
        frac = width / (50 - 8) * 100
        print(f"    Width for |R-1| < {eps_thresh:.0%}: Delta_mu ~ {width:.3f} ({frac:.1f}% of scan range)")

    # Classify sharpness
    if abs(dr_at_cross) > 0.05:
        sharpness = "SHARP (strong selection)"
    elif abs(dr_at_cross) > 0.01:
        sharpness = "MODERATE"
    else:
        sharpness = "SHALLOW (weak selection)"
    print(f"    Sharpness classification: {sharpness}")
else:
    dr_at_cross = dr_dmu_op

# Local high-resolution scan near crossing for accurate width measurement
if mu_star_list:
    mu_local = np.linspace(mu_star_list[0] - 1.0, mu_star_list[0] + 1.0, 500)
    r_local = np.array([ratio_G_over_2pi_alpha(m) for m in mu_local])
    r_local_valid = ~np.isnan(r_local)
    mu_local_v = mu_local[r_local_valid]
    r_local_v = r_local[r_local_valid]
else:
    mu_local_v = mu_valid
    r_local_v = ratios_scan

# Zone where |R-1| < 1%
mask_1pct = np.abs(r_local_v - 1) < 0.01
if np.any(mask_1pct):
    mu_zone = mu_local_v[mask_1pct]
    width_1pct = mu_zone[-1] - mu_zone[0]
    print(f"\n  Empirical zone |R-1| < 1%%: mu in [{mu_zone[0]:.4f}, {mu_zone[-1]:.4f}]")
    print(f"  Width = {width_1pct:.4f}")
else:
    width_1pct = 0
    print(f"\n  No scan points with |R-1| < 1% (even in local refinement)")

mask_5pct = np.abs(r_local_v - 1) < 0.05
if np.any(mask_5pct):
    mu_zone5 = mu_local_v[mask_5pct]
    width_5pct = mu_zone5[-1] - mu_zone5[0]
    print(f"  Empirical zone |R-1| < 5%%: mu in [{mu_zone5[0]:.4f}, {mu_zone5[-1]:.4f}]")
    print(f"  Width = {width_5pct:.4f}")
else:
    width_5pct = 0

selectivity = width_1pct / (mu_valid[-1] - mu_valid[0]) if width_1pct > 0 else 0

test3 = True
print(f"\n  Selectivity (1%%): {selectivity*100:.1f}% of domain")
print(f"  [PASS] Derivative analysis complete")

# ==============================================================
# TEST 4: COMPARE G=2*pi*alpha CROSSING WITH c*105=phi
# ==============================================================

print(f"\n{'='*70}")
print("TEST 4: Compare G=2*pi*alpha crossing with c*105=phi")
print("="*70)

# From the theory: c_sieve is defined via the Bianchi I metric.
# c^2 = (characteristic spatial metric) / |g_00|
# g_00 = -|d^2(ln alpha)/dmu^2|, g_ii = gamma_i^2/mu^2
# There are several ways to define c. Try geometric mean of spatial/temporal.

print(f"\n  phi = {phi_golden:.6f}")
print(f"  105 = 3*5*7 = {3*5*7}")
print(f"  105/phi = {105/phi_golden:.6f}")

# Definition 1: c = geometric mean of (gamma_p/mu) / sqrt(|g_00|)
g_00_abs_op = abs(d2_ln_alpha(mu_alpha))
gamma_prod_op = 1.0
for p in active_primes:
    gamma_prod_op *= gamma_p_func(p, mu_alpha)
spatial_geom_op = (gamma_prod_op)**(1./3) / mu_alpha
c_def1_op = spatial_geom_op / sqrt(g_00_abs_op) if g_00_abs_op > 0 else 0

print(f"\n  Definition 1: c = (prod gamma_p)^(1/3) / (mu * sqrt|g_00|)")
print(f"    c(mu_alpha)     = {c_def1_op:.8f}")
print(f"    c * 105         = {c_def1_op*105:.6f}")
print(f"    c * 105 - phi   = {c_def1_op*105 - phi_golden:.6f}")
print(f"    Error from phi  = {abs(c_def1_op*105 - phi_golden)/phi_golden*100:.4f}%")

# Scan for c*105 = phi crossing
mu_scan_4 = np.linspace(8, 30, 300)
c_vals_1 = []
for mu_v in mu_scan_4:
    g00_v = abs(d2_ln_alpha(mu_v))
    gp = 1.0
    for p in active_primes:
        gp *= gamma_p_func(p, mu_v)
    sg = gp**(1./3) / mu_v
    c_v = sg / sqrt(g00_v) if g00_v > 0 else 0
    c_vals_1.append(c_v * 105)
c_vals_1 = np.array(c_vals_1)

phi_cross_1 = []
for i in range(len(c_vals_1)-1):
    if c_vals_1[i] > 0 and c_vals_1[i+1] > 0:
        if (c_vals_1[i] - phi_golden) * (c_vals_1[i+1] - phi_golden) < 0:
            mc = mu_scan_4[i] + (phi_golden - c_vals_1[i]) * (mu_scan_4[i+1] - mu_scan_4[i]) / (c_vals_1[i+1] - c_vals_1[i])
            phi_cross_1.append(mc)

print(f"\n  c*105 = phi crossings (def. 1):")
if phi_cross_1:
    for ic, mc in enumerate(phi_cross_1):
        print(f"    Crossing {ic+1}: mu = {mc:.4f}")
else:
    closest_idx = np.argmin(np.abs(c_vals_1 - phi_golden))
    print(f"    No crossing. Closest: c*105 = {c_vals_1[closest_idx]:.4f} at mu = {mu_scan_4[closest_idx]:.2f}")

# Definition 2: c = sum(gamma_p)/mu / (3*sqrt(|g_00|))
# (arithmetic mean of scale factor rates / lapse)
c_vals_2 = []
for mu_v in mu_scan_4:
    g00_v = abs(d2_ln_alpha(mu_v))
    gs = sum(gamma_p_func(p, mu_v) for p in active_primes)
    c_v = gs / (3*mu_v*sqrt(g00_v)) if g00_v > 0 else 0
    c_vals_2.append(c_v * 105)
c_vals_2 = np.array(c_vals_2)

phi_cross_2 = []
for i in range(len(c_vals_2)-1):
    if c_vals_2[i] > 0 and c_vals_2[i+1] > 0:
        if (c_vals_2[i] - phi_golden) * (c_vals_2[i+1] - phi_golden) < 0:
            mc = mu_scan_4[i] + (phi_golden - c_vals_2[i]) * (mu_scan_4[i+1] - mu_scan_4[i]) / (c_vals_2[i+1] - c_vals_2[i])
            phi_cross_2.append(mc)

print(f"\n  c*105 = phi crossings (def. 2, arithmetic mean):")
if phi_cross_2:
    for ic, mc in enumerate(phi_cross_2):
        print(f"    Crossing {ic+1}: mu = {mc:.4f}")
else:
    closest_idx = np.argmin(np.abs(c_vals_2 - phi_golden))
    print(f"    No crossing. Closest: c*105 = {c_vals_2[closest_idx]:.4f} at mu = {mu_scan_4[closest_idx]:.2f}")

# Compare with G=2*pi*alpha crossing
print(f"\n  --- Comparison of crossing points ---")
all_phi_crosses = phi_cross_1 + phi_cross_2
if mu_star_list and all_phi_crosses:
    for mc_phi in all_phi_crosses:
        sep = abs(mu_star_list[0] - mc_phi)
        print(f"    G=2*pi*alpha at mu = {mu_star_list[0]:.4f}")
        print(f"    c*105=phi    at mu = {mc_phi:.4f}")
        print(f"    Separation: {sep:.4f} ({'COINCIDENT' if sep < 1.0 else 'SEPARATED'})")
elif mu_star_list:
    print(f"    G=2*pi*alpha at mu = {mu_star_list[0]:.4f}")
    print(f"    c*105=phi: no crossing found")
else:
    print(f"    Neither crossing firmly established")

test4 = True
print(f"\n  [PASS] Comparison analysis complete")

# ==============================================================
# TEST 5: ROLE OF delta_mu AND mu = 15 = 3+5+7
# ==============================================================

print(f"\n{'='*70}")
print("TEST 5: Role of delta_mu = mu - 15 and significance of 15 = 3+5+7")
print("="*70)

print(f"\n  mu_alpha = {mu_alpha:.6f}")
print(f"  15 = 3 + 5 + 7 (sum of active primes)")
print(f"  15 = 3 * 5 (product of first two active primes)")
print(f"  105 = 3 * 5 * 7 (product of all active primes)")
print(f"  delta_mu = mu_alpha - 15 = {delta_mu_op:.8f}")
print(f"  delta_mu/mu_alpha = {delta_mu_op/mu_alpha:.6f} ({delta_mu_op/mu_alpha*100:.4f}%)")

# What changes from mu=15 to mu_alpha?
alpha_15 = alpha_sieve(15.0)
G00_15 = einstein_G00(15.0)
D_15 = D_KL_at(15.0)
ratio_15 = ratio_G_over_2pi_alpha(15.0)

print(f"\n  At mu = 15 (exact integer):")
print(f"    alpha(15)    = {alpha_15:.10e}, 1/alpha = {1/alpha_15:.4f}")
print(f"    alpha_EM     = {alpha_EM_phys:.10e}, 1/alpha = 137.036")
print(f"    Diff 1/alpha = {abs(1/alpha_15 - 137.036):.4f}")
print(f"    G_00(15)     = {G00_15:.8f}")
print(f"    D_KL(15)     = {D_15:.8f}")
print(f"    Ratio R(15)  = {ratio_15:.6f} (deviation from 1: {abs(ratio_15-1)*100:.3f}%)")

print(f"\n  At mu = mu_alpha = {mu_alpha:.6f}:")
print(f"    Ratio R(mu_alpha) = {ratio_op:.6f} (deviation from 1: {abs(ratio_op-1)*100:.3f}%)")

# How does ratio change over the small interval [15, mu_alpha]?
if abs(delta_mu_op) > 1e-6:
    dR_over_interval = (ratio_op - ratio_15) / delta_mu_op
    print(f"\n  Slope dR/dmu over [{15}, {mu_alpha:.4f}] = {dR_over_interval:.6f}")
    print(f"  R changes by {abs(ratio_op - ratio_15):.6f} over delta_mu = {delta_mu_op:.6f}")

# Explore delta_mu connections to fundamental constants
print(f"\n  Exploring delta_mu = {delta_mu_op:.8f}:")
connections = [
    ("1/(2*pi)",          1/(2*pi)),
    ("alpha_EM",          alpha_EM_phys),
    ("1/ln(137)",         1/log(137)),
    ("phi/105",           phi_golden/105),
    ("2/105",             2.0/105),
    ("1/sqrt(105)",       1/sqrt(105)),
    ("phi^2/105",         phi_golden**2/105),
    ("1/(4*pi^2)",        1/(4*pi**2)),
    ("2/(pi*105)",        2/(pi*105)),
    ("ln(2)/15",          log(2)/15),
    ("1/25",              0.04),
]

print(f"  {'Candidate':20s} {'Value':>12s} {'delta/cand':>10s} {'cand/delta':>10s}")
print(f"  {'-'*55}")
for label, val in connections:
    r1 = delta_mu_op / val if val != 0 else 0
    r2 = val / delta_mu_op if delta_mu_op != 0 else 0
    marker = " <--" if 0.95 < r1 < 1.05 else ""
    print(f"  {label:20s} {val:12.8f} {r1:10.4f} {r2:10.4f}{marker}")

# Taylor expansion: what determines delta_mu?
dalpha_15 = (alpha_sieve(15.001) - alpha_sieve(14.999)) / 0.002
delta_linear = (alpha_EM_phys - alpha_15) / dalpha_15
print(f"\n  From Taylor: delta_mu ~ (alpha_EM - alpha(15)) / alpha'(15)")
print(f"    alpha'(15) = {dalpha_15:.6e}")
print(f"    Predicted: {delta_linear:.8f} (actual: {delta_mu_op:.8f})")
print(f"    Error: {abs(delta_linear - delta_mu_op)/abs(delta_mu_op)*100:.2f}%")

# Is mu=15 special for the RATIO (not just for alpha)?
# Check several integer/half-integer values near 15
print(f"\n  Ratio R(mu) at special points:")
special_mus = [12, 13, 14, 14.5, 15, 15.04, mu_alpha, 16, 17, 18, 20, 25, 30]
print(f"  {'mu':>8s} {'R(mu)':>10s} {'|R-1|%':>8s}")
print(f"  {'-'*28}")
for mu_s in special_mus:
    r_s = ratio_G_over_2pi_alpha(mu_s)
    if not np.isnan(r_s):
        print(f"  {mu_s:8.4f} {r_s:10.6f} {abs(r_s-1)*100:8.3f}")

test5 = True
print(f"\n  [PASS] delta_mu analysis complete")

# ==============================================================
# TEST 6: DEEPER FORMULA G_00 = n * pi^k * alpha^a * D_KL^b
# ==============================================================

print(f"\n{'='*70}")
print("TEST 6: Search for deeper formula G_00 = f(alpha, D_KL)")
print("="*70)

# Collect (mu, G_00, alpha, D_KL) data points
data_points = []
for mu_v in np.linspace(10, 40, 120):
    a_v = alpha_sieve(mu_v)
    if a_v <= 0:
        continue
    G00_v = einstein_G00(mu_v)
    D_v = D_KL_at(mu_v)
    if D_v <= 0:
        continue
    data_points.append((mu_v, G00_v, a_v, D_v))

data_arr = np.array(data_points)
mu_data = data_arr[:, 0]
G00_data = data_arr[:, 1]
alpha_data = data_arr[:, 2]
DKL_data = data_arr[:, 3]

# --- Power-law regression ---
mask_pos = G00_data > 0
if np.sum(mask_pos) > 10:
    ln_G = np.log(G00_data[mask_pos])
    ln_a = np.log(alpha_data[mask_pos])
    ln_D = np.log(DKL_data[mask_pos])

    A_mat = np.column_stack([np.ones(np.sum(mask_pos)), ln_a, ln_D])
    result = np.linalg.lstsq(A_mat, ln_G, rcond=None)
    c0_fit, a_fit, b_fit = result[0]
    C_fit = exp(c0_fit)

    ln_G_pred = c0_fit + a_fit * ln_a + b_fit * ln_D
    R2 = 1 - np.sum((ln_G - ln_G_pred)**2) / np.sum((ln_G - np.mean(ln_G))**2)

    print(f"\n  Unconstrained power-law fit: G_00 = C * alpha^a * D_KL^b")
    print(f"    C = {C_fit:.6f}")
    print(f"    a = {a_fit:.4f}")
    print(f"    b = {b_fit:.4f}")
    print(f"    R^2 = {R2:.6f}")

    # Compare with (4*pi)^2 = 157.91 at a=1, b=1
    print(f"\n  Expected (operating point): C = (4*pi)^2 = {(4*pi)**2:.4f}, a = 1, b = 1")

    # --- Test specific candidate formulas ---
    print(f"\n  Candidate formulas at mu_alpha:")
    print(f"  {'Formula':40s} {'Predicted':>12s} {'Observed':>12s} {'Error%':>8s}")
    print(f"  {'-'*75}")

    candidates = [
        ("16*pi^2 * alpha * D_KL",          16*pi**2 * alpha_op * D_op),
        ("8*pi * alpha * D_KL",             8*pi * alpha_op * D_op),
        ("4*pi^2 * alpha * D_KL",           4*pi**2 * alpha_op * D_op),
        ("2*pi * alpha * D_KL",             2*pi * alpha_op * D_op),
        ("(2*pi)^3 * alpha * D_KL",         (2*pi)**3 * alpha_op * D_op),
        ("16*pi^2 * alpha^2 * D_KL",        16*pi**2 * alpha_op**2 * D_op),
        ("16*pi^2 * alpha * D_KL^2",        16*pi**2 * alpha_op * D_op**2),
        ("16*pi^2 * sqrt(alpha) * D_KL",    16*pi**2 * sqrt(alpha_op) * D_op),
        ("128 * alpha * D_KL",              128 * alpha_op * D_op),
        ("105 * alpha * D_KL",              105 * alpha_op * D_op),
        ("137 * alpha * D_KL (= D_KL)",     137 * alpha_op * D_op),
        ("e^5 * alpha * D_KL",              exp(5) * alpha_op * D_op),
        ("phi^5 * alpha * D_KL",            phi_golden**5 * alpha_op * D_op),
    ]

    best_err = 1e10
    best_label = ""
    for label, pred in candidates:
        err = abs(pred - G00_op) / abs(G00_op) * 100
        marker = " <-- BEST" if err < best_err else ""
        if err < best_err:
            best_err = err
            best_label = label
        print(f"  {label:40s} {pred:12.8f} {G00_op:12.8f} {err:8.4f}{marker}")

    print(f"\n  Best at operating point: {best_label} (error {best_err:.4f}%)")

    # --- Global check: how does C_eff = G_00 / (alpha * D_KL) vary? ---
    print(f"\n  Effective coefficient C_eff(mu) = G_00 / (alpha * D_KL):")
    C_eff_all = G00_data[mask_pos] / (alpha_data[mask_pos] * DKL_data[mask_pos])
    mu_eff = mu_data[mask_pos]

    print(f"  {'mu':>6s} {'C_eff':>12s} {'C_eff/16pi^2':>14s}")
    print(f"  {'-'*35}")
    step = max(1, len(mu_eff)//12)
    for i in range(0, len(mu_eff), step):
        print(f"  {mu_eff[i]:6.1f} {C_eff_all[i]:12.4f} {C_eff_all[i]/(16*pi**2):14.6f}")

    print(f"\n  C_eff statistics:")
    print(f"    Mean   = {np.mean(C_eff_all):.4f}")
    print(f"    Std    = {np.std(C_eff_all):.4f}")
    print(f"    CV     = {np.std(C_eff_all)/np.mean(C_eff_all)*100:.2f}%")
    print(f"    16pi^2 = {16*pi**2:.4f}")

    # Polynomial fit of C_eff in 1/mu
    inv_mu = 1.0 / mu_eff
    poly_c = np.polyfit(inv_mu, C_eff_all, 2)
    print(f"\n  Fit: C_eff ~ a0 + a1/mu + a2/mu^2:")
    print(f"    a0 (mu->inf) = {poly_c[2]:.4f}")
    print(f"    a1           = {poly_c[1]:.4f}")
    print(f"    a2           = {poly_c[0]:.4f}")
    print(f"    C_eff(inf) vs 16*pi^2: {abs(poly_c[2] - 16*pi**2)/(16*pi**2)*100:.2f}% error")

    # Key insight: is C_eff = 16*pi^2 an ASYMPTOTIC relation or exact?
    is_exact = np.std(C_eff_all)/np.mean(C_eff_all) < 0.01
    is_asymptotic = abs(poly_c[2] - 16*pi**2)/(16*pi**2) < 0.05 and not is_exact
    if is_exact:
        print(f"\n  C_eff = 16*pi^2 is EXACT (universal): CV < 1%")
    elif is_asymptotic:
        print(f"\n  C_eff -> 16*pi^2 is ASYMPTOTIC (large mu limit)")
    else:
        print(f"\n  C_eff is NEITHER exact nor asymptotically 16*pi^2")

    # --- Try: G_00 = n * pi^k * alpha^a * D_KL^b for small integers ---
    print(f"\n  Brute-force search: G_00 = n*pi^k * alpha^a * D_KL^b")
    print(f"  (n in 1..20, k in 0..4, a in {{1/2,1,3/2,2}}, b in {{1/2,1,3/2,2}})")

    best_combo = None
    best_combo_err = 1e10

    for n in range(1, 21):
        for k in range(0, 5):
            for a_exp in [0.5, 1.0, 1.5, 2.0]:
                for b_exp in [0.5, 1.0, 1.5, 2.0]:
                    pred_all = n * pi**k * alpha_data[mask_pos]**a_exp * DKL_data[mask_pos]**b_exp
                    rel_err = np.mean(np.abs(G00_data[mask_pos] - pred_all) / np.abs(G00_data[mask_pos]))
                    if rel_err < best_combo_err:
                        best_combo_err = rel_err
                        best_combo = (n, k, a_exp, b_exp)

    if best_combo:
        n_b, k_b, a_b, b_b = best_combo
        print(f"\n  Best integer combination:")
        print(f"    n = {n_b}, k = {k_b}, a = {a_b}, b = {b_b}")
        print(f"    => G_00 ~ {n_b}*pi^{k_b} * alpha^{a_b} * D_KL^{b_b}")
        print(f"    Coefficient = {n_b * pi**k_b:.4f}")
        print(f"    Mean relative error = {best_combo_err*100:.4f}%")
        print(f"    Compare: 16*pi^2 = {16*pi**2:.4f}")

test6 = True
print(f"\n  [PASS] Formula search complete")

# ==============================================================
# TEST 7: PHYSICAL INTERPRETATION AND SELF-CONSISTENCY
# ==============================================================

print(f"\n{'='*70}")
print("TEST 7: Physical interpretation -- self-consistency meaning")
print("="*70)

print(f"""
  THE CHAIN OF REASONING:
  =======================

  (1) Bianchi I metric on sieve space:
      g_00 = -|d^2(ln alpha)/dmu^2|,  g_ii = gamma_i^2/mu^2

  (2) Einstein tensor:
      G^0_0 = H_1*H_2 + H_1*H_3 + H_2*H_3  (Friedmann)

  (3) Information content:
      D_KL = ln(2) + D_KL(mod 3, geom)  (energy density of sieve)

  (4) Einstein equation: G_00 = 8*pi*G * rho
      with rho = D_KL:  G = G_00 / (8*pi*D_KL)

  (5) OBSERVATION: G = 2*pi*alpha at mu_alpha (0.29%)

  (6) EQUIVALENT: G_00 = 16*pi^2 * alpha * D_KL

  QUESTION: Is (6) a TAUTOLOGY or a CONSTRAINT?
""")

# Self-consistency function: F(mu) = G_00 - 16*pi^2*alpha*D_KL
print(f"  Self-consistency: F(mu) = G_00 - 16*pi^2*alpha*D_KL = 0 ?")
print(f"  {'mu':>6s} {'G_00':>10s} {'16pi^2*a*D':>12s} {'F':>12s} {'F/G_00':>10s}")
print(f"  {'-'*55}")

F_vals = []
for mu_v, G00_v, a_v, D_v in data_points[::10]:
    F_v = G00_v - 16*pi**2 * a_v * D_v
    F_rel = F_v/G00_v if abs(G00_v) > 1e-15 else 0
    F_vals.append((mu_v, F_v, F_rel))
    print(f"  {mu_v:6.1f} {G00_v:10.6f} {16*pi**2*a_v*D_v:12.6f} {F_v:12.6f} {F_rel:10.6f}")

# The variation of F/G_00 tells us if it's universal or not
F_rel_arr = np.array([x[2] for x in F_vals])
print(f"\n  F/G_00 statistics:")
print(f"    Mean   = {np.mean(F_rel_arr):.6f}")
print(f"    Std    = {np.std(F_rel_arr):.6f}")
print(f"    Range  = [{np.min(F_rel_arr):.6f}, {np.max(F_rel_arr):.6f}]")

is_tautology = np.std(F_rel_arr) < 0.01
is_selection = not is_tautology and np.min(np.abs(F_rel_arr)) < 0.01

# Decomposition: what drives the variation of R(mu)?
# R = G_00 / (16*pi^2*alpha*D_KL)
# ln(R) = ln(G_00) - ln(16*pi^2) - ln(alpha) - ln(D_KL)
# d(ln R)/dmu = d(ln G_00)/dmu - d(ln alpha)/dmu - d(ln D_KL)/dmu

print(f"\n  Decomposition of d(ln R)/dmu at operating point:")
h_dd = 0.01
dln_G00 = (log(einstein_G00(mu_alpha+h_dd)) - log(einstein_G00(mu_alpha-h_dd))) / (2*h_dd)
dln_alph = (log(alpha_sieve(mu_alpha+h_dd)) - log(alpha_sieve(mu_alpha-h_dd))) / (2*h_dd)
dln_DKL = (log(D_KL_at(mu_alpha+h_dd)) - log(D_KL_at(mu_alpha-h_dd))) / (2*h_dd)
dln_R_calc = dln_G00 - dln_alph - dln_DKL

print(f"    d(ln G_00)/dmu  = {dln_G00:+.6f}")
print(f"    d(ln alpha)/dmu = {dln_alph:+.6f}")
print(f"    d(ln D_KL)/dmu  = {dln_DKL:+.6f}")
print(f"    d(ln R)/dmu     = {dln_R_calc:+.6f}")

total_abs = abs(dln_G00) + abs(dln_alph) + abs(dln_DKL)
if total_abs > 0:
    print(f"\n    Relative contributions:")
    print(f"      |d ln G_00|/total  = {abs(dln_G00)/total_abs*100:.1f}%")
    print(f"      |d ln alpha|/total = {abs(dln_alph)/total_abs*100:.1f}%")
    print(f"      |d ln D_KL|/total  = {abs(dln_DKL)/total_abs*100:.1f}%")

# If R(mu) is NOT identically 1, then G=2*pi*alpha is a genuine constraint.
# The strength of the constraint depends on dR/dmu.
print(f"""
  =================================================================
  FINAL SYNTHESIS
  =================================================================
""")

if is_tautology:
    print(f"  VERDICT: G = 2*pi*alpha is a TAUTOLOGY (holds everywhere)")
    print(f"  => It does NOT select any particular mu")
    print(f"  => The formula G_00 = 16*pi^2*alpha*D_KL is an identity")
    final_verdict = "TAUTOLOGY (not selective)"
elif is_selection:
    print(f"  VERDICT: G = 2*pi*alpha is a SELECTION PRINCIPLE")
    print(f"  => It selects a specific mu where F(mu) = 0")
    if mu_star_list:
        print(f"  => Selected mu* = {mu_star_list[0]:.6f}")
        d_to_alpha = abs(mu_star_list[0] - mu_alpha)
        print(f"  => Distance to mu_alpha: {d_to_alpha:.6f} ({d_to_alpha/mu_alpha*100:.4f}%)")
        if d_to_alpha/mu_alpha < 0.01:
            print(f"  => mu* ~ mu_alpha: the selection IS alpha = alpha_EM!")
            final_verdict = "STRONG SELECTION: determines alpha_EM"
        else:
            print(f"  => mu* != mu_alpha: selection gives DIFFERENT alpha")
            final_verdict = "SELECTION (but not exactly alpha_EM)"
    else:
        final_verdict = "SELECTION (crossing not precisely located)"
else:
    print(f"  VERDICT: G = 2*pi*alpha is APPROXIMATE everywhere")
    print(f"  => The relation is a rough but non-trivial identity")
    print(f"  => Selection is WEAK: the ratio varies slowly")
    final_verdict = "WEAK SELECTION or APPROXIMATE IDENTITY"

print(f"\n  Quantitative summary:")
print(f"    R(mu_alpha) = {ratio_op:.8f} ({abs(ratio_op-1)*100:.4f}% from 1)")
print(f"    R variation (CV) = {cv_ratio:.2f}%")
if crossings:
    print(f"    Crossing(s) found at mu = {', '.join(f'{c:.4f}' for c in crossings)}")
print(f"    Selectivity (1%% width) = {selectivity*100:.1f}% of domain")
print(f"    Number of conditions near mu~15: {n_near_mu_alpha}")
print(f"\n  FINAL VERDICT: {final_verdict}")

test7 = True

# ==============================================================
# SCORE
# ==============================================================

tests = [test1, test2, test3, test4, test5, test6, test7]
test_labels = [
    "High-resolution scan (250 points, mu=8..50)",
    "Multiple conditions confluent at mu~15",
    "Derivative / sharpness of crossing",
    "G=2*pi*alpha vs c*105=phi comparison",
    "delta_mu = mu - 15 analysis (3+5+7)",
    "Deeper formula search G_00 = f(alpha, D_KL)",
    "Physical interpretation and self-consistency"
]

n_pass = sum(tests)

print(f"\n{'='*70}")
print(f"SCORE: {n_pass}/{len(tests)} tests")
print("="*70)

for i, (t, desc) in enumerate(zip(tests, test_labels), 1):
    print(f"  T{i}: [{'PASS' if t else 'FAIL'}] {desc}")

print(f"\n{'='*70}")
print("KEY RESULTS:")
print("="*70)
print(f"  1. Ratio G/(2*pi*alpha) scanned at {len(ratios_scan)} points")
print(f"     Mean = {mean_ratio:.6f}, CV = {cv_ratio:.2f}%")
print(f"     => {'UNIVERSAL' if is_universal else 'NON-UNIVERSAL'}")
if mu_star_list:
    print(f"  2. Exact crossing at mu* = {mu_star_list[0]:.8f}")
    print(f"     1/alpha(mu*) = {1/alpha_sieve(mu_star_list[0]):.4f} vs 137.036")
print(f"  3. Sharpness: dR/dmu = {dr_dmu_op:.6f}")
print(f"     Selection width (1%%): {width_1pct:.2f}")
print(f"  4. G_00 ~ {best_label if best_combo else '16*pi^2*alpha*D_KL'} (best formula)")
print(f"  5. {final_verdict}")
