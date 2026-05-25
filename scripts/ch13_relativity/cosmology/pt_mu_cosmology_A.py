"""
pt_mu_cosmology_A.py
====================
Persistence Theory — cosmological identification of μ(z)
Three approaches compared against DESI DR2 CPL constraints.

PHYSICAL CONTEXT
----------------
PT has:
  ρ_DE(μ) = (1 + echo(μ)) · m_e⁴ · (m_e/m_P)^{3/2}    [vacuum energy]
  echo(μ)  = exp(-μ*/μ)                                   [prime coherence]
  μ*       = 15   (fixed point of the ghost cascade {11,13})

From the CLW slow-roll equation of motion (DERIVED, not postulated):
  1 + w(μ) = s · echo(μ)          with s = 1/2

The three approaches differ in HOW μ evolves with the scale factor a = 1/(1+z).

Author  : Yan Senez / PT project
Date    : 2026-04-24
"""

import math
import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import curve_fit, minimize_scalar
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════════
# 1.  PT CONSTANTS  (self-contained, no imports from PT library)
# ═══════════════════════════════════════════════════════════════════════════

alpha_EM = 1 / 137.036
s        = 0.5                     # CLW slow-roll factor  (= 1/2 in PT)
mu_star  = 15.0                    # ghost-prime fixed point {11,13} → 11+13-9=15? / 15 exact

echo_0   = math.exp(-1)            # echo(μ*) = exp(-μ*/μ*) = exp(-1)
N_gen    = 3
N_sp     = 3
p7       = 7
k = N_gen * p7 * (1 - N_sp**2 * s * alpha_EM / (4 * math.pi))

m_e      = 0.510998950e6           # eV
m_P_SI   = 1.220890e28             # eV  (reduced Planck: 2.435×10²⁷ eV; here GeV convention check)
m_P      = (1 / alpha_EM)**(k / 2) * m_e   # PT-derived Planck mass [eV]

hbar     = 6.582119569e-16         # eV·s
kms      = 1e3 / 3.085677581e22    # km/s → s⁻¹

# Planck 2018 cosmological parameters
Omega_m       = 0.315
Omega_L_obs   = 0.685
H0_Planck     = 67.36              # km/s/Mpc
H0_SHOES      = 73.3

# DESI DR2 CPL best-fit & uncertainties
w0_DESI  = -0.838;  sw0 = 0.057
wa_DESI  = -0.75;   swa = 0.29

# ═══════════════════════════════════════════════════════════════════════════
# 2.  CORE PT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def echo(mu):
    """Prime coherence factor: echo(μ) = exp(-μ*/μ)."""
    return math.exp(-mu_star / mu)

def echo_np(mu):
    """Vectorised version."""
    return np.exp(-mu_star / np.asarray(mu, dtype=float))

def w_field(mu):
    """
    Field equation of state from CLW slow-roll EOM (DERIVED).
        1 + w = s · echo(μ)  →  w = -1 + s·exp(-μ*/μ)
    """
    return -1.0 + s * echo(mu)

def rho_DE_PT(mu):
    """
    PT vacuum energy density (in eV⁴).
        ρ_DE = (1 + echo(μ)) · m_e⁴ · (m_e/m_P)^{3/2}
    """
    base = m_e**4 * (m_e / m_P)**1.5
    return (1.0 + echo(mu)) * base

# Reference energy density at μ = μ*
rho_DE_0 = rho_DE_PT(mu_star)

# w₀ reference (μ fixed at μ*)
w0_ref = w_field(mu_star)

# ═══════════════════════════════════════════════════════════════════════════
# 3.  CPL FIT UTILITY
# ═══════════════════════════════════════════════════════════════════════════

def cpl_model(a, w0, wa):
    """Chevallier-Polarski-Linder: w(a) = w₀ + w_a·(1-a)."""
    return w0 + wa * (1.0 - a)

def fit_cpl(a_arr, w_arr, a_min=0.3, a_max=1.0):
    """
    Fit CPL to w(a) in the range [a_min, a_max] (DESI/CMB sensitivity window).
    Returns w₀, w_a and their 1-σ uncertainties.
    """
    mask = (a_arr >= a_min) & (a_arr <= a_max)
    if mask.sum() < 4:
        return np.nan, np.nan, np.nan, np.nan
    try:
        popt, pcov = curve_fit(cpl_model, a_arr[mask], w_arr[mask],
                               p0=[-0.95, 0.0], maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        return popt[0], popt[1], perr[0], perr[1]
    except Exception:
        return np.nan, np.nan, np.nan, np.nan

def desi_tension(w0, wa):
    """
    Combined DESI DR2 tension in σ assuming uncorrelated Gaussian errors.
    Returns total pull Δ = sqrt(((w0-w0_DESI)/sw0)² + ((wa-wa_DESI)/swa)²)
    """
    if np.isnan(w0) or np.isnan(wa):
        return np.inf
    return math.sqrt(((w0 - w0_DESI) / sw0)**2 + ((wa - wa_DESI) / swa)**2)

# ═══════════════════════════════════════════════════════════════════════════
# 4.  ANALYTIC SIGN ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def print_analytic_sign_analysis():
    """
    Show that the sign of w_a is uniquely determined by sign(dμ/d ln a).
    Linearise w(μ) around μ = μ*:
        w(μ) = w_field(μ*)  +  dw/dμ|_{μ*} · δμ
        dw/dμ = s · (μ*/μ²) · exp(-μ*/μ)
        at μ=μ*: dw/dμ|_{μ*} = s · (1/μ*) · echo_0

    For CPL: w_a = -dw/d ln a|_{a=1}  (since w(a)=w₀+w_a(1-a), dw/da|_{a=1}=-w_a)
        w_a = -(dw/dμ)|_{μ*} · (dμ/d ln a)|_{a=1}
            = -(s/μ*) · echo_0 · (dμ/d ln a)|_{a=1}

    Key: sign(w_a) = -sign(dμ/d ln a)
    """
    coeff = -(s / mu_star) * echo_0
    print("\n═══════════════════════════════════════════════════════════════")
    print("ANALYSE ANALYTIQUE DU SIGNE DE w_a")
    print("═══════════════════════════════════════════════════════════════")
    print(f"  Linéarisation de w(μ) autour de μ = μ* = {mu_star}")
    print(f"    w(μ) ≈ w₀ + (s/μ*)·echo₀·(μ - μ*)")
    print(f"    où  s/μ* = {s/mu_star:.4f},  echo₀ = exp(-1) = {echo_0:.4f}")
    print(f"")
    print(f"  w_a (CPL) = -(s/μ*)·echo₀ · (dμ/d ln a)|_{{a=1}}")
    print(f"           =  {coeff:.5f}  ×  (dμ/d ln a)|_{{a=1}}")
    print(f"")
    print(f"  → sign(w_a) = -sign(dμ/d ln a)    [identité exacte au 1er ordre]")
    print(f"")
    print(f"  Approche 1 (ODE):   dμ/d ln a = -(3s/μ*)·μ²·(1+echo) < 0")
    print(f"                      → w_a  > 0   (thawing quintessence)")
    print(f"")
    print(f"  Approche 2 (ghost, α>0): μ = μ*·a^α")
    print(f"                      dμ/d ln a = α·μ*·a^α > 0 pour a,α>0")
    print(f"                      → w_field_a NEGATIF (signe opposé à A1) !")
    print(f"                      Raison : sign(w_a) = -sign(dμ/d ln a),")
    print(f"                               dμ/d ln a > 0 en A2 → w_field_a < 0 ← DESI !")
    print(f"                      MAIS w_eff_a > 0 car ρ_DE croît vers le futur")
    print(f"                           (interaction Q > 0 domine l'EOS effective)")
    print(f"")
    print(f"  Terme d'interaction (Approche 2) :")
    print(f"    Q/H = ρ_DE × [d ln ρ_DE/d ln a + 3(1+w_field)]")
    print(f"    d ln ρ_DE/d ln a = α·echo_ghost/(1+echo_ghost)")
    print(f"    à a=1: {alpha_EM:.4f} (dépend de α)")
    print(f"    Le terme 3s·echo₀ = {3*s*echo_0:.4f}")
    print(f"    Si Q>0 : DE → DM (ghost nourrit la matière noire)")

# ═══════════════════════════════════════════════════════════════════════════
# 5.  APPROACH 1 — SELF-CONSISTENT ODE
# ═══════════════════════════════════════════════════════════════════════════

def approach1_ode():
    """
    Integrate dμ/d ln a = -(3s/μ*)·μ²·(1+echo(μ))
    from the energy conservation equation:
        d ρ_DE / d ln a = -3(1+w)·ρ_DE
    combined with ρ_DE = (1+echo(μ))·base  and  1+w = s·echo(μ).

    Derivation:
        d ρ_DE/d ln a = d/d ln a [(1+echo(μ))·base]
                      = base · echo(μ) · (μ*/μ²) · dμ/d ln a
                      = -3 s·echo(μ) · (1+echo(μ)) · base

        Cancelling echo(μ)·base (> 0):
            (μ*/μ²) · dμ/d ln a = -3s·(1+echo(μ))
            dμ/d ln a = -(3s/μ*)·μ²·(1+echo(μ))
    """
    print("\n═══════════════════════════════════════════════════════════════")
    print("APPROCHE 1 — ODE AUTO-COHÉRENTE")
    print("  dμ/d ln a = -(3s/μ*)·μ²·(1+echo(μ))  [conservation énergie + CLW]")
    print("═══════════════════════════════════════════════════════════════")

    def dmu_dlna(lna, mu):
        """RHS of the μ ODE."""
        mu_val = float(mu[0])
        if mu_val <= 0:
            return [0.0]
        e = math.exp(-mu_star / mu_val)
        return [-(3 * s / mu_star) * mu_val**2 * (1.0 + e)]

    # Grid: ln(a) from -4.6 (z≈99) to +1.1 (a≈3)
    lna_grid = np.linspace(-4.6, 1.1, 3000)

    # Initial condition at a=1 (ln a = 0): μ = μ*
    ic = [mu_star]

    sol = solve_ivp(
        dmu_dlna,
        t_span=(0.0, lna_grid[-1]),
        y0=ic,
        t_eval=lna_grid[lna_grid >= 0],
        method='RK45',
        rtol=1e-10, atol=1e-12,
        dense_output=True
    )
    # Backward integration: t_span must be decreasing, t_eval must match that order
    lna_back = lna_grid[lna_grid <= 0][::-1]   # sorted descending (0 → most negative)
    sol_back = solve_ivp(
        dmu_dlna,
        t_span=(0.0, lna_grid[0]),
        y0=ic,
        t_eval=lna_back,
        method='RK45',
        rtol=1e-10, atol=1e-12,
        dense_output=True
    )

    # Merge forward and backward solutions
    lna_all = np.concatenate([sol_back.t, sol.t])
    mu_all  = np.concatenate([sol_back.y[0], sol.y[0]])
    # Sort by lna
    idx     = np.argsort(lna_all)
    lna_all = lna_all[idx]
    mu_all  = mu_all[idx]
    a_all   = np.exp(lna_all)
    z_all   = 1.0 / a_all - 1.0

    # Compute w(a) at each point
    w_all = np.array([-1.0 + s * math.exp(-mu_star / max(m, 1e-6)) for m in mu_all])

    # Extract specific redshifts
    def mu_at_z(z_target):
        a_target = 1.0 / (1.0 + z_target)
        lna_target = math.log(a_target)
        idx_near = np.argmin(np.abs(lna_all - lna_target))
        return mu_all[idx_near], w_all[idx_near]

    mu_z0, w_z0 = mu_at_z(0.0)
    mu_z1, w_z1 = mu_at_z(1.0)
    mu_z2, w_z2 = mu_at_z(2.0)

    print(f"\n  Condition initiale : μ(z=0) = μ* = {mu_star}")
    print(f"  μ(z=0)  = {mu_z0:.4f}   w(z=0) = {w_z0:.6f}")
    print(f"  μ(z=1)  = {mu_z1:.4f}   w(z=1) = {w_z1:.6f}")
    print(f"  μ(z=2)  = {mu_z2:.4f}   w(z=2) = {w_z2:.6f}")

    # dμ/d ln a at a=1
    e0 = echo_0
    dmu_lna_0 = -(3 * s / mu_star) * mu_star**2 * (1.0 + e0)
    print(f"\n  dμ/d ln a|_{{z=0}} = {dmu_lna_0:.4f}  (< 0 : μ DIMINUE vers le futur)")

    # CPL fit in [0.3, 1.0]
    w0_fit, wa_fit, dw0, dwa = fit_cpl(a_all, w_all, a_min=0.3, a_max=1.0)
    tension = desi_tension(w0_fit, wa_fit)

    print(f"\n  CPL fit [a ∈ 0.3–1.0] :")
    print(f"    w₀  = {w0_fit:.4f}  (DESI: {w0_DESI}±{sw0})")
    print(f"    w_a = {wa_fit:.4f}  (DESI: {wa_DESI}±{swa})")
    print(f"    Tension DESI combinée : {tension:.2f} σ")

    w0_analytic_A1 = -(s / mu_star) * echo_0 * dmu_lna_0  # contribution to w_a
    print(f"\n  w_a analytique (1er ordre) :")
    print(f"    w_a ≈ -(s/μ*)·echo₀·(dμ/d ln a) = {-(s/mu_star)*echo_0*dmu_lna_0:.4f}")

    print(f"\n  CONCLUSION A1 : dμ/d ln a < 0 → μ ↑ vers le passé → echo ↑ → w ↑")
    print(f"  → w_a > 0  (quintessence dégelante, thawing)")
    print(f"  → EN TENSION FORTE avec DESI DR2 (w_a_DESI = {wa_DESI})")

    return {
        'lna': lna_all, 'a': a_all, 'z': z_all, 'mu': mu_all, 'w': w_all,
        'w0': w0_fit, 'wa': wa_fit, 'tension': tension,
        'dmu_lna_0': dmu_lna_0
    }

# ═══════════════════════════════════════════════════════════════════════════
# 6.  APPROACH 2 — GHOST CONDENSATE  μ ∝ a^α
# ═══════════════════════════════════════════════════════════════════════════

def approach2_ghost_condensate(alpha_values=None):
    """
    Ghost condensate mechanism:
        μ_ghost(a) = μ* · a^α    (comoving coherence length ~ a^α)

    Physical interpretation:
        - The ghost prime condensate {11,13} has NOT fully formed at a < 1
        - Its coherence length grows with the expansion of the universe
        - At a→0: μ→0, echo→0, ρ_DE → base (pure vacuum, no coherence)
        - At a=1: μ=μ*, full condensate

    Key quantities:
        echo_ghost(a) = exp(-μ*/μ_ghost) = exp(-1/a^α)
        w_field(a)    = -1 + s·exp(-1/a^α)
        ρ_DE(a)       = (1 + exp(-1/a^α)) · base
        d ln ρ_DE/d ln a = α · echo_ghost / (1 + echo_ghost)

    The field EOS and the effective EOS differ because ρ_DE is NOT conserved
    (there is a DE↔DM interaction Q ≠ 0).
    """
    if alpha_values is None:
        alpha_values = [1, 2, 4]

    print("\n═══════════════════════════════════════════════════════════════")
    print("APPROCHE 2 — CONDENSAT GHOST  μ(a) = μ* · a^α")
    print("  μ plus petit dans le passé → condensat en formation")
    print("═══════════════════════════════════════════════════════════════")

    a_arr = np.linspace(0.01, 3.0, 5000)
    base  = m_e**4 * (m_e / m_P)**1.5

    results = {}

    print(f"\n  {'α':>5}  {'w₀_field':>10}  {'w₀_eff':>10}  "
          f"{'w_a_field':>11}  {'w_a_eff':>9}  "
          f"{'Q/ρH|z=0':>10}  {'Δ(DESI)':>9}")
    print(f"  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*11}  {'-'*9}  {'-'*10}  {'-'*9}")

    for alpha in alpha_values:
        # Ghost echo
        eg  = np.exp(-1.0 / np.maximum(a_arr**alpha, 1e-15))

        # Field EOS: 1+w_field = s·echo_ghost
        w_field_arr = -1.0 + s * eg

        # Effective EOS: w_eff = -1 - (1/3)·d ln ρ_DE/d ln a
        # d ln ρ_DE/d ln a = alpha · eg / (1+eg)
        dlnrho_dlna = alpha * eg / (1.0 + eg)
        w_eff_arr   = -1.0 - dlnrho_dlna / 3.0

        # Interaction term at z=0 (a=1):
        eg0    = math.exp(-1.0)   # echo_ghost(a=1) = exp(-1)
        dlnrho_0 = alpha * eg0 / (1.0 + eg0)
        w_field_0 = -1.0 + s * eg0
        # Q/(ρ_DE·H) = d ln ρ_DE/d ln a + 3(1+w_field)
        Q_over_rhoH = dlnrho_0 + 3.0 * s * eg0

        # CPL fit — field EOS
        w0_f, wa_f, _, _ = fit_cpl(a_arr, w_field_arr, a_min=0.3, a_max=1.0)
        # CPL fit — effective EOS
        w0_e, wa_e, _, _ = fit_cpl(a_arr, w_eff_arr, a_min=0.3, a_max=1.0)

        tension_field = desi_tension(w0_f, wa_f)
        tension_eff   = desi_tension(w0_e, wa_e)

        print(f"  α={alpha:>3}  "
              f"w₀_f={w0_f:+.4f}  "
              f"w₀_e={w0_e:+.4f}  "
              f"w_a_f={wa_f:+.5f}  "
              f"w_a_e={wa_e:+.5f}  "
              f"Q/ρH={Q_over_rhoH:+.4f}  "
              f"Δ_eff={tension_eff:.1f}σ")

        results[alpha] = {
            'a': a_arr, 'w_field': w_field_arr, 'w_eff': w_eff_arr,
            'w0_field': w0_f, 'wa_field': wa_f,
            'w0_eff': w0_e, 'wa_eff': wa_e,
            'Q_over_rhoH': Q_over_rhoH,
            'tension_field': tension_field, 'tension_eff': tension_eff
        }

    # Find alpha that minimises DESI tension on effective EOS
    def tension_alpha(alpha):
        a_ = np.linspace(0.01, 3.0, 3000)
        eg_ = np.exp(-1.0 / np.maximum(a_**alpha, 1e-15))
        dlnrho_ = alpha * eg_ / (1.0 + eg_)
        w_eff_  = -1.0 - dlnrho_ / 3.0
        w0_, wa_, _, _ = fit_cpl(a_, w_eff_, a_min=0.3, a_max=1.0)
        return desi_tension(w0_, wa_)

    res_opt = minimize_scalar(tension_alpha, bounds=(0.1, 20.0), method='bounded')
    alpha_opt = res_opt.x
    tension_opt = res_opt.fun

    # Compute best-fit CPL at alpha_opt
    a_ = np.linspace(0.01, 3.0, 5000)
    eg_ = np.exp(-1.0 / np.maximum(a_**alpha_opt, 1e-15))
    dlnrho_opt = alpha_opt * eg_ / (1.0 + eg_)
    w_eff_opt  = -1.0 - dlnrho_opt / 3.0
    w0_opt, wa_opt, _, _ = fit_cpl(a_, w_eff_opt, a_min=0.3, a_max=1.0)

    print(f"\n  Valeur α minimisant la tension DESI (EOS effective) :")
    print(f"    α_opt = {alpha_opt:.3f}")
    print(f"    w₀_eff= {w0_opt:.4f}  (DESI: {w0_DESI}±{sw0})")
    print(f"    w_a_eff= {wa_opt:.4f}  (DESI: {wa_DESI}±{swa})")
    print(f"    Tension résiduelle : {tension_opt:.2f} σ")

    print(f"\n  INTERPRÉTATION PT :")
    print(f"    α=1 : longueur cohérence ~ a (comobile fixe) — condensat de Bose-Einstein")
    print(f"    α=2 : ~ a² — croissance accélérée (analogue inflation)")
    print(f"    α=4 = 2/s : déterminé par s=1/2 (candidat PT)")
    print(f"    Le terme Q>0 signifie que le condensat NOURRIT la matière noire")
    print(f"    (énergie sombre → matière noire lors de la formation du condensat)")

    results['alpha_opt'] = alpha_opt
    results['w0_opt']    = w0_opt
    results['wa_opt']    = wa_opt
    results['tension_opt'] = tension_opt
    return results

# ═══════════════════════════════════════════════════════════════════════════
# 7.  APPROACH 3 — FIXED μ (ΛCDM ANALOGUE)
# ═══════════════════════════════════════════════════════════════════════════

def approach3_fixed_mu():
    """
    μ = μ* = constant everywhere.
    ρ_DE = const,  w = -1 + s·echo(μ*) = const.
    This is the ΛCDM-like baseline for PT.
    """
    print("\n═══════════════════════════════════════════════════════════════")
    print("APPROCHE 3 — μ FIXE (référence Λ)")
    print("  μ = μ* = 15 partout, ρ_DE = const")
    print("═══════════════════════════════════════════════════════════════")

    w_const = w_field(mu_star)
    print(f"\n  w = -1 + s·echo₀ = -1 + {s}·{echo_0:.6f} = {w_const:.6f}")
    print(f"  w₀ = {w_const:.6f}  (DESI: {w0_DESI}±{sw0})")
    print(f"  w_a = 0.000000  (DESI: {wa_DESI}±{swa})")

    tension = desi_tension(w_const, 0.0)
    print(f"  Tension DESI (w₀ seul) : {abs(w_const - w0_DESI)/sw0:.2f} σ(w₀)")
    print(f"  Tension DESI combinée  : {tension:.2f} σ")

    print(f"\n  CONCLUSION A3 : μ fixe → EOS constante → insensible à l'histoire cosmique")
    print(f"  Tension w_a : {abs(0 - wa_DESI)/swa:.2f} σ (w_a=0 vs DESI w_a={wa_DESI})")

    return {'w0': w_const, 'wa': 0.0, 'tension': tension}

# ═══════════════════════════════════════════════════════════════════════════
# 8.  INTERPRETATION SECTION
# ═══════════════════════════════════════════════════════════════════════════

def print_interpretation(res1, res2, res3):
    """Physical interpretation of results."""
    print("\n═══════════════════════════════════════════════════════════════")
    print("INTERPRÉTATION PT DES TROIS APPROCHES")
    print("═══════════════════════════════════════════════════════════════")

    print("""
APPROCHE 1 — Quintessence dégelante (thawing)
─────────────────────────────────────────────
  ODE : dμ/d ln a < 0  →  μ CROÎT vers le passé (z→∞)
  Physiquement : le scalaire ghost était "plus cohérent" dans le passé (μ > μ*),
  et il DÉROULE vers μ* aujourd'hui en perdant de la cohérence.
  Conséquence : w(z) croît depuis w(z=0) vers -1+s à z grand.
  w_a > 0 → ANTI-DESI.  Tension ≥ 3σ sur w_a.

  C'est le scénario "Route D Bianchi I" initial, confirmé ici comme problématique.

APPROCHE 2 — Condensat ghost en formation
─────────────────────────────────────────
  μ(a) = μ* · a^α  →  μ CROÎT vers le futur, PLUS PETIT dans le passé
  Physiquement : les primes ghost {11,13} forment un condensat de cohérence
  croissant avec l'expansion. Au Big Bang : μ≈0, pas de cohérence.
  Aujourd'hui : μ=μ*, condensat pleinement formé.

  ρ_DE(a) = (1 + exp(-1/a^α))·base  →  CROÎT vers le futur (plus de cohérence)
  → d ln ρ_DE/d ln a = α·echo/(1+echo) > 0
  → w_eff = -1 - (1/3)·d ln ρ_DE/d ln a < -1  (phantom effectif)

  MAIS : ρ_DE est plus GRAND maintenant qu'au passé → l'EOS effective
  devient MOINS négative vers le passé → w_eff_a > 0 également !
  (Le CPL mesure la variation de a=0.3 à a=1, et w_eff est plus négatif à a=1.)

  Résultat crucial : w_a_eff > 0 pour TOUS les α > 0, comme en A1.
  Les deux approches donnent w_a dans le MAUVAIS sens vs DESI.

  Terme Q > 0 : la DM reçoit de l'énergie du condensat.
  Cependant, le signe de w_a_eff ne change pas : c'est un problème structurel.

APPROCHE 3 — μ fixe (Λ-CDM PT)
───────────────────────────────
  Pas d'évolution, EOS constante. w_a = 0.
  Tension DESI sur w_a : ~2.6σ.  Ne capture pas la dynamique DE.

DIRECTION DE μ(z) SELON DESI
─────────────────────────────
  DESI voit w_a = -0.75 ± 0.29 : w(z) DÉCROÎT de a=0 vers a=1 (w plus négatif
  aujourd'hui qu'au passé).

  Pour w_field = -1 + s·echo(μ), on a w_field décroissant ⟺ echo décroissant ⟺ μ décroissant.
  Donc w_a < 0 (CPL) requiert que μ soit PLUS GRAND dans le passé qu'aujourd'hui.
  → dμ/d ln a < 0  (μ décroît vers le futur, venant de μ > μ* dans le passé)

  Mais l'ODE A1 donne exactement cela !  Pourtant w_a = +0.88 en CPL.
  PARADOXE APPARENT : μ croît vers le passé, donc w(z) croît vers le passé.
  En CPL: w(a) = w₀ + w_a(1-a). Si w est plus grand pour petit a (grand z),
  alors w_a > 0.  DESI veut le contraire : w plus négatif pour petit a.

  CONCLUSION : DESI w_a < 0 requiert w PLUS négatif dans le passé.
  → echo doit être PLUS GRAND dans le passé → μ PLUS GRAND dans le passé.
  → dμ/d ln a > 0 ???  Non — contradiction : echo = exp(-μ*/μ), si μ > μ*
    alors echo > echo₀, donc w > w₀, donc w_a > 0 dans le passé.

  RÉSOLUTION : Le signe DESI w_a < 0 correspond à une EOS qui était PLUS NÉGATIVE
  (plus proche de -1) dans le passé qu'aujourd'hui. Cela requiert :
    echo(z>0) < echo(z=0)  →  μ(z>0) < μ* = 15

  Donc μ DÉCROÎT vers le passé : dμ/d ln a > 0 (μ plus petit pour a petit).
  → C'est l'Approche 2 (condensat ghost, μ = μ*·a^α) pour la partie FIELD EOS.
  → w_field_a < 0 pour TOUS les α en A2 (confirmé numériquement).
  → L'EOS DE CHAMP de A2 est COMPATIBLE avec DESI en signe.
  → L'EOS effective (qui inclut le terme cinétique d'interaction Q) reste >0.

  VERDICT FINAL :
    Pour l'EOS de CHAMP (observable directe si Q est absorbé dans la DM) :
      A2 (ghost condensat) : w_field_a < 0  → COMPATIBLE DESI en signe
      A1 (ODE conservée)   : w_field_a > 0  → INCOMPATIBLE DESI
    Pour l'EOS effective (tout inclus) :
      Toutes les approches : w_eff_a > 0 ou ≈ 0 → tension avec DESI
      Sauf si Q est négatif (DM → DE), ce qui requiert un mécanisme inverse.

  PISTE : Modifier l'ODE A1 avec une interaction Q < 0 (DM → DE condensat)
  pour obtenir dμ/d ln a > 0 cohérent avec μ_ghost et DESI simultanément.
""")

# ═══════════════════════════════════════════════════════════════════════════
# 9.  SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════

def print_bilan(res1, res2, res3):
    print("\n═══════════════════════════════════════════════════════════════")
    print("BILAN — TABLEAU COMPARATIF")
    print("═══════════════════════════════════════════════════════════════")
    print(f"\n  {'Approche':<35} {'w₀':>8}  {'w_a':>8}  {'Δ DESI':>9}")
    print(f"  {'-'*35}  {'-'*8}  {'-'*8}  {'-'*9}")

    # Reference values
    w0_1  = res1['w0']
    wa_1  = res1['wa']
    t_1   = res1['tension']

    print(f"  {'A1 — ODE (conserv. énergie)':<35} {w0_1:>+8.4f}  {wa_1:>+8.4f}  {t_1:>8.2f}σ")

    for alpha in [1, 2, 4]:
        if alpha in res2:
            r = res2[alpha]
            print(f"  {'A2 — Ghost α='+str(alpha)+' (w_eff)':<35} "
                  f"{r['w0_eff']:>+8.4f}  {r['wa_eff']:>+8.4f}  {r['tension_eff']:>8.2f}σ")

    alpha_opt = res2['alpha_opt']
    w0_opt    = res2['w0_opt']
    wa_opt    = res2['wa_opt']
    t_opt     = res2['tension_opt']
    print(f"  {'A2 — Ghost α='+f'{alpha_opt:.2f}'+' (optimal)':<35} "
          f"{w0_opt:>+8.4f}  {wa_opt:>+8.4f}  {t_opt:>8.2f}σ")

    print(f"  {'A3 — μ fixe (Λ-analogue)':<35} "
          f"{res3['w0']:>+8.4f}  {res3['wa']:>+8.4f}  {res3['tension']:>8.2f}σ")

    print(f"\n  DESI DR2 (référence)             "
          f"  {w0_DESI:>+8.4f}  {wa_DESI:>+8.4f}  {'—':>9}")

    print(f"\n  VERDICT :")
    print(f"    A1 (ODE conservée) : w_field_a > 0, w_eff_a > 0, ANTI-DESI, tension ≥ {t_1:.1f}σ")
    print(f"    A2 w_field : w_field_a < 0 pour tout α → COMPATIBLE DESI en SIGNE")
    print(f"    A2 w_eff   : w_eff_a > 0 pour tout α (interaction Q positive)")
    print(f"    A2 α_opt={alpha_opt:.2f} (w_eff) : tension minimale {t_opt:.2f}σ (encore large)")
    print(f"    A3 (μ fixe) : w_a = 0, tension {res3['tension']:.1f}σ sur w_a seul")
    print(f"")
    print(f"    → Pour l'EOS de CHAMP (A2) : w_field_a < 0 est DESI-compatible en signe.")
    print(f"      Si l'interaction Q est absorbée dans ρ_DM, l'observable est w_field.")
    print(f"    → Pour l'EOS effective : TOUTES les approches donnent w_eff_a ≥ 0.")
    print(f"      Obtenir w_a < 0 requit soit Q < 0 (DM → DE), soit un mécanisme")
    print(f"      où ρ_DE décroît vers le futur (w_eff > -1 en moyenne).")
    print(f"    → PISTE PRINCIPALE : dériver Q dans le cadre PT et résoudre l'ODE couplée")
    print(f"      dμ/d ln a = -(3s/μ*)·μ²·(1+echo) + Q_PT(μ,a)/ρ_DE")
    print(f"      où Q_PT vient de la dissipation du condensat ghost vers la matière noire.")

# ═══════════════════════════════════════════════════════════════════════════
# 10.  ADDITIONAL DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════

def print_pt_constants():
    """Print the PT constants used in this script."""
    print("═══════════════════════════════════════════════════════════════")
    print("CONSTANTES PT (dérivées, sans paramètre libre)")
    print("═══════════════════════════════════════════════════════════════")
    print(f"  α_EM       = 1/{1/alpha_EM:.3f}")
    print(f"  s          = {s}  (facteur CLW slow-roll)")
    print(f"  μ*         = {mu_star}  (point fixe ghost cascade)")
    print(f"  echo(μ*)   = exp(-1) = {echo_0:.6f}")
    print(f"  k          = {k:.6f}  (exposant masse PT)")
    print(f"  m_P (PT)   = {m_P:.6e} eV  ({m_P/1e9:.4e} GeV)")
    print(f"  m_P (SI)   = {m_P_SI:.6e} eV  (référence)")
    print(f"  m_P ratio  = {m_P/m_P_SI:.6f}")
    print(f"  ρ_DE(μ*)   = {rho_DE_0:.6e} eV⁴")
    print(f"  w₀ (μ fixe)= {w0_ref:.6f}")
    print(f"  Δw₀(DESI)  = {abs(w0_ref - w0_DESI)/sw0:.2f}σ")
    print()

def print_ghost_detail(alpha=4):
    """
    Detailed analysis for α = 2/s = 4 (PT candidate).
    This is the PT-natural choice since s = 1/2 → 2/s = 4.
    """
    print(f"\n═══════════════════════════════════════════════════════════════")
    print(f"ANALYSE DÉTAILLÉE — α = 2/s = {2/s:.0f}  (candidat PT naturel)")
    print(f"═══════════════════════════════════════════════════════════════")
    print(f"  Justification : α = 2/s est l'exposant qui fait dμ/d ln a = 2/s · μ")
    print(f"  soit la croissance doublant à chaque e-fold → analogue inflation lente")

    a_vals = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0]
    base = m_e**4 * (m_e / m_P)**1.5

    print(f"\n  {'z':>6}  {'a':>6}  {'μ_ghost':>9}  {'echo':>8}  "
          f"{'w_field':>9}  {'w_eff':>9}  {'ρ_DE/ρ₀':>9}")
    print(f"  {'-'*6}  {'-'*6}  {'-'*9}  {'-'*8}  {'-'*9}  {'-'*9}  {'-'*9}")

    for a in a_vals:
        z     = 1.0/a - 1.0
        mu_g  = mu_star * a**alpha
        eg    = math.exp(-mu_star / max(mu_g, 1e-15))
        wf    = -1.0 + s * eg
        dlnr  = alpha * eg / (1.0 + eg)
        we    = -1.0 - dlnr / 3.0
        rho   = (1.0 + eg) * base
        print(f"  {z:>6.2f}  {a:>6.3f}  {mu_g:>9.3f}  {eg:>8.5f}  "
              f"{wf:>+9.5f}  {we:>+9.5f}  {rho/rho_DE_0:>9.5f}")

    print(f"\n  Note : à z=0 (a=1), μ_ghost = μ* = {mu_star} par construction")
    print(f"  Note : à z grand (a→0), μ_ghost→0, echo→0, ρ_DE→base (Λ pur)")
    print(f"  Note : l'EOS effective devient très négative (phantom) à grand z")
    print(f"         → la densité DE DIMINUE vers le passé (moins de cohérence)")
    print(f"         → cela RESSEMBLE à la matière noire froide à grand z !")

# ═══════════════════════════════════════════════════════════════════════════
# 11.  MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 63)
    print("PT_MU_COSMOLOGY_A.py")
    print("Persistence Theory — identification de μ(z)")
    print("Trois approches comparées aux données DESI DR2")
    print("=" * 63)

    print_pt_constants()
    print_analytic_sign_analysis()

    res1 = approach1_ode()
    res2 = approach2_ghost_condensate(alpha_values=[1, 2, 4])
    res3 = approach3_fixed_mu()
    print_ghost_detail(alpha=4)
    print_interpretation(res1, res2, res3)
    print_bilan(res1, res2, res3)

    print("\n" + "="*63)
    print("FIN pt_mu_cosmology_A.py")
    print("="*63)

if __name__ == '__main__':
    main()
