"""
q1q2q3_extended.py

Extensions des Q1, Q2, Q3 :
  Q1+ : test complexe (caractères, parties Re/Im, modules) sur les constantes math
  Q2+ : running multi-énergies (1 GeV, 5 GeV, 10 GeV, M_Z, 200 GeV)
  Q3+ : a_e via ghost VP analogue à a_μ ch21 (PT-derivation)
"""
from __future__ import annotations
import math
import cmath
from itertools import product as iproduct
from pathlib import Path

# ============================================================
# Setup PT exact (mpmath précision 30 décimales)
# ============================================================

import mpmath as mp
mp.mp.dps = 30

MU_STAR = mp.mpf(15)
Q_STAT = mp.mpf(13) / mp.mpf(15)
Q_THERM = mp.exp(mp.mpf(-1) / mp.mpf(15))

ACTIVE = [3, 5, 7]
GHOST = [11, 13]

def delta_p(p, q): return (1 - q**p) / p
def sin2_p(p, q): d = delta_p(p, q); return d * (2 - d)

S2_STAT = {p: sin2_p(p, Q_STAT) for p in ACTIVE + GHOST}
S2_THERM = {p: sin2_p(p, Q_THERM) for p in ACTIVE + GHOST}

# gamma_p formula: 4p q^(p-1) (1-delta) / (mu (1-q^p) (2-delta))
def gamma_p(p, mu, q):
    qp = q**p
    d = (1 - qp) / p
    return 4*p*(q**(p-1))*(1-d)/(mu*(1-qp)*(2-d))

GAMMA = {p: gamma_p(p, MU_STAR, Q_STAT) for p in ACTIVE + GHOST + [17, 19, 23]}

# Caractères complexes (theta_p = arccos(1-delta), then complex chi)
# theta_p (réel from PT)
def theta_p(p, q):
    d = delta_p(p, q)
    cos_t = 1 - d
    if abs(cos_t) > 1: return 0
    return mp.acos(cos_t)

THETA_STAT = {p: theta_p(p, Q_STAT) for p in ACTIVE + GHOST}
THETA_THERM = {p: theta_p(p, Q_THERM) for p in ACTIVE + GHOST}

# Caractères PT typiques (issus du spin foam):
#   chi_4(p) = exp(2πi × p / 5) (quartic mod 5)
#   chi_3(p) = exp(2πi × p / 7 mod 3) (cubic mod 7)
# Mais surtout les phases naturelles à mu*=15 = exp(2πi × theta_p / 2pi) = exp(i theta_p)
PSI_STAT = {p: mp.exp(1j * THETA_STAT[p]) for p in ACTIVE + GHOST}

# ============================================================
# Q1+ : Constantes math classiques avec complexes
# ============================================================

print("=" * 90)
print("Q1+ : Constantes math classiques — base étendue avec complexes")
print("=" * 90)

# Build extended candidate list (real + complex projections)
candidates = {
    "q_stat": float(Q_STAT),
    "q_therm": float(Q_THERM),
    "1/q_stat": float(1/Q_STAT),
    "1/q_therm": float(1/Q_THERM),
    "ln(q_stat)": float(mp.log(Q_STAT)),  # negative
    "ln(q_therm)": float(mp.log(Q_THERM)),
    "ln(15)": math.log(15),
    "ln(8/30)": math.log(8/30),  # negative
    "log_2(15)": math.log2(15),
    "ln(2pi)": math.log(2*math.pi),
    "2/pi": 2/math.pi,
    "1/pi": 1/math.pi,
}

# Add real parts
for p in ACTIVE + GHOST:
    candidates[f"sin²_{p}(stat)"] = float(S2_STAT[p])
    candidates[f"sin²_{p}(therm)"] = float(S2_THERM[p])
    candidates[f"theta_{p}(stat)"] = float(THETA_STAT[p])
    candidates[f"theta_{p}(stat)/pi"] = float(THETA_STAT[p] / mp.pi)
    candidates[f"|psi_{p}(stat)|"] = float(abs(PSI_STAT[p]))  # = 1
    candidates[f"Re psi_{p}"] = float(PSI_STAT[p].real)
    candidates[f"Im psi_{p}"] = float(PSI_STAT[p].imag)
for p in ACTIVE + GHOST:
    candidates[f"gamma_{p}"] = float(GAMMA[p])

# Sums and products
candidates["Sum sin²(stat)"] = sum(float(S2_STAT[p]) for p in ACTIVE)
candidates["Sum sin²(therm)"] = sum(float(S2_THERM[p]) for p in ACTIVE)
candidates["Prod sin²(stat)"] = float(mp.fprod(S2_STAT[p] for p in ACTIVE))
candidates["Prod sin²(therm)"] = float(mp.fprod(S2_THERM[p] for p in ACTIVE))
candidates["Sum gamma_active"] = sum(float(GAMMA[p]) for p in ACTIVE)
candidates["Sum gamma_ghost"] = sum(float(GAMMA[p]) for p in GHOST)
candidates["Sum theta(stat)"] = sum(float(THETA_STAT[p]) for p in ACTIVE)
candidates["Sum Re psi"] = sum(float(PSI_STAT[p].real) for p in ACTIVE)
candidates["Sum Im psi"] = sum(float(PSI_STAT[p].imag) for p in ACTIVE)
candidates["|Sum psi|"] = float(abs(sum(PSI_STAT[p] for p in ACTIVE)))
candidates["Sum |psi|"] = sum(float(abs(PSI_STAT[p])) for p in ACTIVE)  # = 3
candidates["arg(Sum psi)"] = float(mp.arg(sum(PSI_STAT[p] for p in ACTIVE)))
# Cosine sums (related to PT geometry)
candidates["cos(theta_3)"] = float(mp.cos(THETA_STAT[3]))
candidates["cos(theta_5)"] = float(mp.cos(THETA_STAT[5]))
candidates["cos(theta_7)"] = float(mp.cos(THETA_STAT[7]))
candidates["sin(theta_3)"] = float(mp.sin(THETA_STAT[3]))
candidates["sin(theta_5)"] = float(mp.sin(THETA_STAT[5]))
candidates["sin(theta_7)"] = float(mp.sin(THETA_STAT[7]))

print(f"  Total candidates: {len(candidates)}")

# Targets
TARGETS = {
    "Mertens M": 0.2614972128476428,
    "Brun B_2 (twin)": 1.9021605823,
    "Brun B_4 (cousin)": 1.197044,
    "Khintchine K_0": 2.6854520010653,
    "Euler-Mascheroni γ": 0.5772156649015329,
    "Catalan G": 0.9159655941772190,
    "Apéry ζ(3)": 1.2020569031595943,
    "Twin-prime const HL": 0.6601618158,
    "e^γ": math.exp(0.5772156649),
}

def search(target, cands, tol_rel=1e-6):
    found = []
    for n, v in cands.items():
        if abs(v) < 1e-18: continue
        rel = abs(v - target) / abs(target)
        if rel < tol_rel: found.append((rel, f"{n} = {v}"))
    for n1, v1 in cands.items():
        for n2, v2 in cands.items():
            if n1 == n2: continue
            for op_name, op in [('+', lambda a,b:a+b), ('-', lambda a,b:a-b),
                                 ('*', lambda a,b:a*b), ('/', lambda a,b:a/b if abs(b)>1e-15 else float('inf'))]:
                try:
                    val = op(v1, v2)
                    if abs(val) < 1e-15: continue
                    rel = abs(val - target) / abs(target)
                    if rel < tol_rel:
                        found.append((rel, f"({n1}) {op_name} ({n2}) = {val:.10f}"))
                except (ZeroDivisionError, OverflowError):
                    pass
    found.sort()
    return found[:3]

print(f"\n{'Constante':<25} | {'Best (sub-ppm)':<60} | reluctance")
print("-" * 100)
for tname, tval in TARGETS.items():
    best = search(tval, candidates, tol_rel=1e-5)  # tol = 10 ppm
    if best:
        rel_ppm = best[0][0] * 1e6
        print(f"{tname:<25} | {best[0][1][:55]:<55} | {rel_ppm:.2f} ppm")
    else:
        # Bonferroni: with len(cands)^2 * 4 ops ~ 100k tests, 0.1% noise expected
        # Try at 100 ppm, then 1000 ppm
        best100 = search(tval, candidates, tol_rel=1e-4)
        if best100:
            rel_ppm = best100[0][0] * 1e6
            print(f"{tname:<25} | {best100[0][1][:55]:<55} | {rel_ppm:.0f} ppm (no sub-ppm)")
        else:
            print(f"{tname:<25} | {'no match below 100 ppm':<55} | n/a")

# ============================================================
# Q2+ : Running multi-énergies
# ============================================================

print("\n" + "=" * 90)
print("Q2+ : Running α_EM à 5 énergies — chercher μ_eff(E) tq α_PT(μ_eff) = α(E)")
print("=" * 90)

# alpha_EM running (PDG values, leptons + hadrons)
# These are 1/alpha at various scales:
RUNNING_DATA = {
    1.0: 137.0359,    # Q ~ 1 GeV (just above α(0))
    2.0: 134.4,       # Q = 2 GeV (below charm threshold)
    5.0: 132.5,       # Q = 5 GeV (below bottom)
    10.0: 131.0,      # Q = 10 GeV (above bottom)
    91.1876: 127.952,  # M_Z
    200.0: 127.5,      # 200 GeV
}

def alpha_inv_naive(mu_pt, primes=ACTIVE):
    if mu_pt <= 2: return float('inf')
    q = 1 - 2/mu_pt
    if q <= 0 or q >= 1: return float('inf')
    prod = 1.0
    for p in primes:
        prod *= float(sin2_p(p, mp.mpf(q)))
    return 1/prod if prod > 0 else float('inf')

def find_mu_eff(target_alpha_inv, dressing=0.758):
    """Find mu such that alpha_inv_naive(mu) + dressing = target."""
    target_naive = target_alpha_inv - dressing
    # binary search in [3, 100]
    from scipy.optimize import brentq
    def f(mu):
        return alpha_inv_naive(mu) - target_naive
    try:
        return brentq(f, 3.01, 99.9, xtol=1e-6)
    except ValueError:
        return None

print(f"  {'Q (GeV)':>10} | {'α^-1(Q) PDG':>14} | {'μ_eff (PT)':>12} | {'check α^-1':>12} | {'log10(Q/m_e)':>14}")
print("-" * 80)
m_e = 0.000511  # GeV
for Q, ainv in RUNNING_DATA.items():
    mu_eff = find_mu_eff(ainv)
    log_q = math.log10(Q / m_e)
    if mu_eff is not None:
        check = alpha_inv_naive(mu_eff) + 0.758
        print(f"  {Q:>10.4f} | {ainv:>14.4f} | {mu_eff:>12.4f} | {check:>12.4f} | {log_q:>14.4f}")
    else:
        print(f"  {Q:>10.4f} | {ainv:>14.4f} | {'no match':>12} | --- | {log_q:>14.4f}")

# Look for linear correlation: μ_eff vs log(Q)?
print(f"\n  Fit linéaire μ_eff = a + b·log10(Q/m_e) ?")
data_pts = [(Q, find_mu_eff(ainv)) for Q, ainv in RUNNING_DATA.items()]
data_pts = [(Q, m) for Q, m in data_pts if m is not None]
if len(data_pts) >= 3:
    xs = [math.log10(Q/m_e) for Q, _ in data_pts]
    ys = [m for _, m in data_pts]
    n = len(xs)
    mx = sum(xs)/n; my = sum(ys)/n
    var_x = sum((x-mx)**2 for x in xs)
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    if var_x > 1e-15:
        b = cov/var_x; a = my - b*mx
        ss_tot = sum((y-my)**2 for y in ys)
        ss_res = sum((y - (a+b*x))**2 for x,y in zip(xs,ys))
        r2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 0
        print(f"    Fit: μ_eff = {a:.3f} + {b:.4f} · log10(Q/m_e), R² = {r2:.4f}")
    else:
        print("    Cannot fit (var_x = 0)")

# ============================================================
# Q3+ : a_e via ghost VP (analogue à a_μ ch21)
# ============================================================

print("\n" + "=" * 90)
print("Q3+ : a_e via ghost VP analogue à a_μ ch21")
print("=" * 90)

# From ch21 formula:
#   a_μ = a_QED + a_EW + a_HVP + a_LBL + a_ghost
#   a_ghost = a_HVP^LO × β_ghost × (1 - γ_7)
#   β_ghost = sin²θ_11 · γ_11 + sin²θ_13 · γ_13
#   = 0.1039 (ch21)
#   (1 - γ_7) = 0.4045
# For muon: a_HVP^LO ~ 6.8e-8, a_ghost = 286e-11

# For electron, all hadronic contributions scale as (m_e/m_μ)^2
m_e_GeV = 0.00051100
m_mu_GeV = 0.10566
mass_ratio_sq = (m_e_GeV / m_mu_GeV)**2

# a_HVP for muon (Aoyama et al. 2020 White Paper):
A_HVP_LO_MU = 6938e-11  # leading order HVP
A_HVP_NLO_MU = -98e-11
A_HVP_NNLO_MU = 12e-11

# Scale to electron
A_HVP_LO_E = A_HVP_LO_MU * mass_ratio_sq
print(f"  m_e²/m_μ² = {mass_ratio_sq:.3e}")
print(f"  a_HVP^LO(μ) = {A_HVP_LO_MU*1e11:.0f}e-11")
print(f"  a_HVP^LO(e) = {A_HVP_LO_E:.3e} (= {A_HVP_LO_E*1e11:.4f}e-11)")

# Compute β_ghost from PT (sin²(theta_p, q_stat) × gamma_p for ghost primes p=11,13)
# Note: in ch21, sin² is the q_stat (vertex) version
beta_ghost = float(S2_STAT[11]) * float(GAMMA[11]) + float(S2_STAT[13]) * float(GAMMA[13])
gamma_7 = float(GAMMA[7])
factor_inactive = 1 - gamma_7
print(f"\n  β_ghost = sin²_11 · γ_11 + sin²_13 · γ_13")
print(f"          = {float(S2_STAT[11]):.5f} × {float(GAMMA[11]):.4f} + {float(S2_STAT[13]):.5f} × {float(GAMMA[13]):.4f}")
print(f"          = {beta_ghost:.5f}  (ch21 says 0.1039)")
print(f"  (1 - γ_7) = {factor_inactive:.4f}  (ch21 says 0.4045)")

# Apply ghost formula to electron:
A_GHOST_E = A_HVP_LO_E * beta_ghost * factor_inactive
A_GHOST_MU = A_HVP_LO_MU * beta_ghost * factor_inactive
print(f"\n  a_ghost(μ) = a_HVP^LO(μ) × β_ghost × (1 - γ_7) = {A_GHOST_MU*1e11:.2f}e-11")
print(f"               (ch21 says 286e-11)")
print(f"  a_ghost(e) = {A_GHOST_E:.3e} (= {A_GHOST_E*1e11:.6f}e-11)")

# Now compute a_e PT total
# QED contributions (Aoyama et al.):
ALPHA_PT = 1 / 137.0359988  # use 'PT-required' value derived earlier (matches a_e)
def a_e_qed(alpha):
    a_pi = alpha / math.pi
    return 0.5*a_pi - 0.328478965579193*a_pi**2 + 1.181241456587*a_pi**3 - 1.91206*a_pi**4 + 6.737*a_pi**5

A_E_QED = a_e_qed(ALPHA_PT)
A_E_EW = 0.030e-12
A_E_HAD_DD = 1.693e-12  # data-driven hadronic VP for electron (much smaller than muon)
A_E_LBL = 0.04e-12
# Add PT ghost contribution
A_E_PT_TOTAL = A_E_QED + A_E_EW + A_E_HAD_DD + A_E_LBL + A_GHOST_E

A_E_EXP = 1.15965218073e-3
print(f"\n  a_e (electron) total :")
print(f"    a_QED          = {A_E_QED:.13e}")
print(f"    a_EW           = {A_E_EW:.3e}")
print(f"    a_HVP (DD)     = {A_E_HAD_DD:.3e}")
print(f"    a_LBL          = {A_E_LBL:.3e}")
print(f"    a_ghost (PT)   = {A_GHOST_E:.3e}")
print(f"    -------------")
print(f"    a_e^PT_total   = {A_E_PT_TOTAL:.13e}")
print(f"    a_e^exp        = {A_E_EXP:.13e}")
print(f"    Δ              = {(A_E_PT_TOTAL - A_E_EXP):.3e}")
print(f"    Δ/a_e          = {(A_E_PT_TOTAL - A_E_EXP)/A_E_EXP:.3e}")

# Check muon too (sanity)
A_MU_EXP = 116592059e-11
A_MU_QED = 116584718.95e-11
A_MU_EW = 153.6e-11
A_MU_HVP_LO = A_HVP_LO_MU
A_MU_HVP_NLO = A_HVP_NLO_MU + A_HVP_NNLO_MU
A_MU_LBL = 92e-11
A_MU_PT_TOTAL = A_MU_QED + A_MU_EW + A_MU_HVP_LO + A_MU_HVP_NLO + A_MU_LBL + A_GHOST_MU
print(f"\n  Sanity check muon:")
print(f"    a_μ^PT_total = {A_MU_PT_TOTAL*1e11:.0f}e-11")
print(f"    ch21 says   = 116592058e-11")
print(f"    Fermilab    = {A_MU_EXP*1e11:.0f}e-11 ± 22e-11")
print(f"    Δ(PT-Ferm) / Ferm = {(A_MU_PT_TOTAL - A_MU_EXP)/A_MU_EXP:.3e}")

out = Path(__file__).parent.parent / "Results" / "q1q2q3_extended.txt"
out.parent.mkdir(exist_ok=True)
with open(out, 'w') as f:
    f.write(f"Q3+ a_e via ghost VP:\n")
    f.write(f"  β_ghost = {beta_ghost:.5f}\n")
    f.write(f"  a_ghost(e) = {A_GHOST_E:.3e}\n")
    f.write(f"  a_e^PT_total = {A_E_PT_TOTAL:.13e}\n")
    f.write(f"  a_e^exp = {A_E_EXP:.13e}\n")
    f.write(f"  Δ/exp = {(A_E_PT_TOTAL - A_E_EXP)/A_E_EXP:.3e}\n")
    f.write(f"\n  a_μ^PT_total = {A_MU_PT_TOTAL*1e11:.0f}e-11\n")
    f.write(f"  a_μ^exp Fermilab = {A_MU_EXP*1e11:.0f}e-11\n")
print(f"\nReport written to {out}")
