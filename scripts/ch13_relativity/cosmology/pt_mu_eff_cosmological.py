#!/usr/bin/env python3
"""
pt_mu_eff_cosmological.py  —  Dérivation de μ_eff(z) depuis Bianchi I
======================================================================
Date    : 2026-04-24
Auteur  : Yan Senez

QUESTION POSÉE :
    La Route D (ch20f) conditionne w_a = -s/e sur μ_eff(z) = μ*/(1+z) [COND].
    Cette hypothèse est-elle cohérente avec les équations Bianchi I de la PT ?

RÉSULTAT PRINCIPAL :
    L'ODE Bianchi I est CINÉMATIQUE (ρ_total s'annule dans le rapport μ̇/(dz/dτ)).
    Elle donne dμ/dz = -1/[(1+z) × H_iso_fac(μ)] > 0.
    ⟹ μ(z) CROÎT avec z (vers le passé cosmologique).
    ⟹ Scénario S2 : Route D a w_a = +s×α/(e) ≈ +0.213 (SIGNE OPPOSÉ à DESI).

STRUCTURE DU SCRIPT :
    Approche 1 — ODE Bianchi I directe (priorité)
    Approche 2 — Self-cohérence par epoch (secondaire)
    Approche 3 — Perturbation linéaire autour de μ*
    Verdict    — Scénario S1–S4 + implications Route D + LaTeX

RÉFÉRENCES :
    ch13_relativity.tex (Bianchi I, Friedmann)
    ch11_couplings.tex  (thm:mu_eff_running)
    ch20f_dark_energy.tex (Route D)
    pt_dark_energy_wz.py (fonctions γ_p, f_p, B, H_iso_fac)
"""

from __future__ import annotations
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

# ══════════════════════════════════════════════════════════════════════════════
# § 0 — Constantes PT et référence DESI
# ══════════════════════════════════════════════════════════════════════════════

MU_STAR  = 15.0          # point fixe PT (T5, ch08)
S        = 0.5           # spin fondamental (T1, ch01)
E        = math.e
PRIMES   = [3, 5, 7]     # primes actifs (T5)

OMEGA_M0  = 0.315        # Planck 2018
OMEGA_L0  = 0.685
H0_KMS    = 67.41        # km/s/Mpc (PT, ch13)

# DESI DR2 (BAO+CMB, 2025)
W0_DESI, W0_ERR = -0.80, 0.06
WA_DESI, WA_ERR = -0.72, 0.25


# ══════════════════════════════════════════════════════════════════════════════
# § 1 — Cinématique Bianchi I (depuis pt_dark_energy_wz.py)
# ══════════════════════════════════════════════════════════════════════════════

def q_minus(mu: float) -> float:
    """q_-(μ) = exp(-1/μ)  [edge/thermique, ch14]"""
    return math.exp(-1.0 / mu)


def gamma_p_formula(mu: float, p: int) -> float:
    """
    Facteur métrique directionnel (prompt eq.) :
        γ_p(μ) = 4p q^{p-1} (1-δ_p) / [μ (1-q^p) (2-δ_p)]
    Utilisé uniquement pour calculer f_p = dγ_p/dμ / γ_p - 1/μ.
    """
    q  = q_minus(mu)
    d  = (1.0 - q**p) / p
    return 4 * p * q**(p - 1) * (1.0 - d) / (mu * (1.0 - q**p) * (2.0 - d))


def dgamma_p_dmu(mu: float, p: int, h: float = 1e-7) -> float:
    """dγ_p/dμ par différences finies centrées (ordre 2)."""
    return (gamma_p_formula(mu + h, p) - gamma_p_formula(mu - h, p)) / (2.0 * h)


def f_p(mu: float, p: int) -> float:
    """
    Coefficient cinématique directionnels :
        f_p(μ) = γ'_p/γ_p - 1/μ   ⟹   H_p = f_p(μ) × μ̇
    f_p < 0 pour μ ≥ 4 (tous p ∈ {3,5,7}).
    """
    g  = gamma_p_formula(mu, p)
    dg = dgamma_p_dmu(mu, p)
    return dg / g - 1.0 / mu


def H_iso_fac(mu: float) -> float:
    """
    Facteur de Hubble isotrope :
        H_iso_fac(μ) = (f_3+f_5+f_7)/3   ⟹   H_iso = H_iso_fac × μ̇
    H_iso_fac < 0 pour μ ≥ 4.
    """
    return sum(f_p(mu, p) for p in PRIMES) / 3.0


def B_mu(mu: float) -> float:
    """
    Invariant de Friedmann-Bianchi :
        B(μ) = f_3 f_5 + f_3 f_7 + f_5 f_7
    Tel que G_00 = B(μ) × μ̇² = 8πG ρ_total/3.
    B > 0 car les f_p sont tous négatifs.
    """
    fp = [f_p(mu, p) for p in PRIMES]
    return fp[0]*fp[1] + fp[0]*fp[2] + fp[1]*fp[2]


# ══════════════════════════════════════════════════════════════════════════════
# § 2 — APPROCHE 1 : ODE Bianchi I directe
# ══════════════════════════════════════════════════════════════════════════════

def dmu_dz_bianchi(mu: float, z: float) -> float:
    """
    ODE purement cinématique (ρ_total s'annule dans μ̇/(dz/dτ)) :

        dμ/dz = -1 / [(1+z) × H_iso_fac(μ)]

    Dérivation :
      μ̇ = -sqrt(ρ_tot/B) (convention expansion)
      dz/dτ = -(1+z) × H_iso_fac(μ) × μ̇
      dμ/dz = μ̇ / (dz/dτ) = -1 / [(1+z) × H_iso_fac(μ)]

    Signe : H_iso_fac < 0  ⟹  dμ/dz > 0  (μ croît vers le passé).
    """
    Hf = H_iso_fac(mu)
    if abs(Hf) < 1e-12:
        return 0.0
    return -1.0 / ((1.0 + z) * Hf)


def solve_mu_of_z(z_max: float = 3.0, n_pts: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """
    Intègre l'ODE dμ/dz = -1/[(1+z) × H_iso_fac(μ)] de z=0 à z=z_max.
    Condition initiale : μ(z=0) = μ* = 15.
    Retourne (z_arr, mu_arr).
    """
    z_span = (0.0, z_max)
    z_eval = np.linspace(0.0, z_max, n_pts)

    def rhs(z, y):
        mu = y[0]
        if mu < 2.0:
            return [0.0]
        return [dmu_dz_bianchi(mu, z)]

    sol = solve_ivp(
        rhs,
        z_span,
        [MU_STAR],
        t_eval=z_eval,
        rtol=1e-10,
        atol=1e-12,
        method="DOP853"
    )
    return sol.t, sol.y[0]


def fit_power_law(z_arr: np.ndarray, mu_arr: np.ndarray) -> tuple[float, float]:
    """
    Ajuste μ(z) = μ* × (1+z)^α par moindres carrés (log-linéaire).
    Retourne (α, résidu_max).
    """
    ln1pz = np.log(1.0 + z_arr[z_arr > 0])
    ln_ratio = np.log(mu_arr[z_arr > 0] / MU_STAR)
    alpha_fit = np.polyfit(ln1pz, ln_ratio, 1)[0]
    resid_max = np.max(np.abs(ln_ratio - alpha_fit * ln1pz))
    return float(alpha_fit), float(resid_max)


# ══════════════════════════════════════════════════════════════════════════════
# § 3 — APPROCHE 3 : Perturbation linéaire autour de μ*
# ══════════════════════════════════════════════════════════════════════════════

def linear_approximation(z_arr: np.ndarray) -> np.ndarray:
    """
    Solution perturbative μ(z) ≈ μ* + C × ln(1+z)
    avec C = -1/H_iso_fac(μ*) = μ* × α (valeur locale à z=0).
    """
    C = -1.0 / H_iso_fac(MU_STAR)
    return MU_STAR + C * np.log(1.0 + z_arr)


def alpha_local(mu: float) -> float:
    """
    Exposant local α(μ) = -1/[μ × H_iso_fac(μ)].
    Serait constant si μ(z) = μ*(1+z)^α exactement.
    """
    Hf = H_iso_fac(mu)
    return -1.0 / (mu * Hf)


# ══════════════════════════════════════════════════════════════════════════════
# § 4 — Équation d'état Route D avec μ(z) dérivé
# ══════════════════════════════════════════════════════════════════════════════

def w_route_D(z: float, mu_eff: float) -> float:
    """
    w(z) = -1 + s × exp(-μ*/μ_eff(z))
    Route D : carrier-echo mechanism (ch20f).
    """
    return -1.0 + S * math.exp(-MU_STAR / mu_eff)


def w_route_D_cond(z: float) -> float:
    """w(z) Route D avec hypothèse COND μ = μ*/(1+z). Donne w_a = -s/e."""
    return -1.0 + S * math.exp(-(1.0 + z))


def cpl_params_analytical(alpha: float) -> tuple[float, float]:
    """
    Paramètres CPL w(a) = w₀ + w_a(1-a) pour Route D avec μ = μ*(1+z)^α :
        w₀ = w(z=0) = -1 + s/e                    [indépendant de α]
        w_a = dw/dz|_{z=0} = s × α / e
    """
    w0 = -1.0 + S / E
    wa = S * alpha / E
    return w0, wa


def cpl_fit(z_arr: np.ndarray, w_arr: np.ndarray) -> tuple[float, float]:
    """Ajustement CPL w = w₀ + w_a × z/(1+z)."""
    def model(z, w0, wa):
        return w0 + wa * z / (1.0 + z)
    popt, _ = curve_fit(model, z_arr, w_arr, p0=[-0.8, 0.0])
    return float(popt[0]), float(popt[1])


# ══════════════════════════════════════════════════════════════════════════════
# § 5 — Diagnostic q₀ en Bianchi I complet
# ══════════════════════════════════════════════════════════════════════════════

def q0_flrw(w0: float) -> float:
    """q₀ FLRW plat : q₀ = Ω_m/2 + Ω_Λ(1+3w₀)/2."""
    return OMEGA_M0 / 2.0 + OMEGA_L0 * (1.0 + 3.0 * w0) / 2.0


def q0_bianchi(mu: float) -> float:
    """
    Paramètre de décélération en Bianchi I :
        q₀ = -1 - Ḣ_iso/H_iso²
    En utilisant H_iso = H_iso_fac(μ) × μ̇ et le rapport Bianchi/Friedmann.
    Ici on utilise l'approximation FLRW avec w_eff → q₀ = Ω_m/2 + Ω_Λ(1+3w_eff)/2.
    """
    w_eff = w_route_D(0.0, mu)
    return q0_flrw(w_eff)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    sep = "=" * 76
    bar = "-" * 76

    print(sep)
    print("PT — DÉRIVATION DE μ_eff(z) DEPUIS BIANCHI I")
    print("Date : 2026-04-24  |  Auteur : Yan Senez")
    print(sep)

    # ─────────────────────────────────────────────────────────────────
    # PRÉAMBULE : Valeurs cinématiques à μ*
    # ─────────────────────────────────────────────────────────────────
    print()
    print("PRÉAMBULE — Cinématique Bianchi I à μ* = 15")
    print(bar)

    Hf_star = H_iso_fac(MU_STAR)
    B_star  = B_mu(MU_STAR)
    fp_star = [f_p(MU_STAR, p) for p in PRIMES]
    alpha_0 = alpha_local(MU_STAR)

    print(f"  f_3(μ*)    = {fp_star[0]:+.6f}")
    print(f"  f_5(μ*)    = {fp_star[1]:+.6f}")
    print(f"  f_7(μ*)    = {fp_star[2]:+.6f}")
    print(f"  H_iso_fac  = (f_3+f_5+f_7)/3 = {Hf_star:+.6f}  [< 0 ✓]")
    print(f"  B(μ*)      = f_3f_5+f_3f_7+f_5f_7 = {B_star:.6f}  [> 0 ✓]")
    print()
    print(f"  ODE : dμ/dz = -1/[(1+z) × H_iso_fac(μ)]")
    print(f"  À z=0 : dμ/dz|_0 = -1/H_iso_fac(μ*) = {-1/Hf_star:.4f}  [> 0 : μ ↑ avec z]")
    print()
    print(f"  OBSERVATION CLEF : H_iso_fac(μ*) < 0 ⟹ dμ/dz > 0")
    print(f"  ⟹ μ(z) CROÎT vers le passé (z > 0).")
    print(f"  ⟹ Hypothèse COND μ = μ*/(1+z) (μ décroissant) est CONTRAIRE à Bianchi I.")

    # ─────────────────────────────────────────────────────────────────
    # APPROCHE 1 : ODE directe
    # ─────────────────────────────────────────────────────────────────
    print()
    print(sep)
    print("APPROCHE 1 — ODE Bianchi I  dμ/dz = -1/[(1+z) × H_iso_fac(μ)]")
    print(sep)

    print()
    print("  Dérivation de l'ODE (cinématique pure) :")
    print("    μ̇ = -√(ρ_tot / B)           (convention : μ décroît pendant expansion)")
    print("    dz/dτ = -(1+z) × H_iso_fac × μ̇")
    print("    dμ/dz = μ̇/(dz/dτ) = -1/[(1+z) × H_iso_fac(μ)]")
    print("    ρ_tot s'ANNULE — l'ODE est INDÉPENDANTE du contenu matière.")

    z_arr, mu_arr = solve_mu_of_z(z_max=3.0, n_pts=600)

    print()
    print(f"  {'z':>5} | {'μ(z)_ODE':>10} | {'μ*(1+z)¹':>11} | {'μ*/(1+z)':>10} | {'dμ/dz':>10}")
    print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*11}-+-{'-'*10}-+-{'-'*10}")
    for ztgt in [0.0, 0.1, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0]:
        idx = np.argmin(np.abs(z_arr - ztgt))
        mu_val = mu_arr[idx]
        mu_lin = MU_STAR * (1.0 + ztgt)
        mu_inv = MU_STAR / (1.0 + ztgt)
        dmdz   = dmu_dz_bianchi(mu_val, ztgt) if mu_val > 2 else 0.0
        print(f"  {ztgt:5.2f} | {mu_val:10.4f} | {mu_lin:11.4f} | {mu_inv:10.4f} | {dmdz:10.4f}")

    # Ajustement puissance
    alpha_fit, resid = fit_power_law(z_arr, mu_arr)
    print()
    print(f"  Ajustement log-linéaire μ(z) = μ* × (1+z)^α :")
    print(f"    α_fit  = {alpha_fit:.4f}   (résidu log max = {resid:.4f})")
    print(f"    α_loc  = {alpha_0:.4f}   (valeur locale -1/(μ* × H_iso_fac(μ*)) à z=0)")
    print()
    print(f"  ⟹ α ≈ {alpha_fit:.2f} > 0 : μ(z) ∝ (1+z)^{alpha_fit:.2f}")
    print(f"  ⟹ μ(z=1) ≈ {MU_STAR * 2**alpha_fit:.2f}  (> μ* = 15 ✓)")
    print(f"  ⟹ μ(z=3) ≈ {MU_STAR * 4**alpha_fit:.2f}")

    # Variation de α avec μ
    print()
    print("  Variation de α(μ) = -1/[μ × H_iso_fac(μ)] :")
    print(f"  {'μ':>6} | {'H_iso_fac':>11} | {'α_local':>9}")
    print(f"  {'-'*6}-+-{'-'*11}-+-{'-'*9}")
    for mu_tst in [5, 8, 10, 12, 15, 18, 20, 30]:
        Hf = H_iso_fac(mu_tst)
        al = -1.0 / (mu_tst * Hf)
        print(f"  {mu_tst:6.1f} | {Hf:+11.6f} | {al:9.4f}")
    print(f"  Note : α → 1 quand μ → ∞ (car μ × H_iso_fac → -1)")

    # ─────────────────────────────────────────────────────────────────
    # APPROCHE 3 : Perturbation linéaire
    # ─────────────────────────────────────────────────────────────────
    print()
    print(sep)
    print("APPROCHE 3 — Perturbation linéaire μ(z) = μ* + C × ln(1+z)")
    print(sep)

    C_lin = -1.0 / Hf_star
    print()
    print(f"  μ(z) = μ* + C × ln(1+z) avec C = -1/H_iso_fac(μ*) = {C_lin:.4f}")
    print(f"         = μ* × (1 + α × ln(1+z))  avec α = C/μ* = {C_lin/MU_STAR:.4f}")
    print()
    print(f"  Comparaison ODE / perturbation linéaire :")
    print(f"  {'z':>5} | {'μ_ODE':>8} | {'μ_lin':>8} | {'μ*(1+z)^α':>11} | {'écart lin':>10}")
    print(f"  {'-'*5}-+-{'-'*8}-+-{'-'*8}-+-{'-'*11}-+-{'-'*10}")

    mu_lin_arr = linear_approximation(z_arr)
    mu_pow_arr = MU_STAR * (1.0 + z_arr)**alpha_fit

    for ztgt in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
        idx     = np.argmin(np.abs(z_arr - ztgt))
        mu_ode  = mu_arr[idx]
        mu_l    = linear_approximation(np.array([ztgt]))[0]
        mu_pw   = MU_STAR * (1 + ztgt)**alpha_fit
        ecart_l = (mu_l - mu_ode) / mu_ode * 100
        print(f"  {ztgt:5.2f} | {mu_ode:8.4f} | {mu_l:8.4f} | {mu_pw:11.4f} | {ecart_l:+9.2f}%")

    print()
    print(f"  ⟹ Perturbation linéaire est raisonnable pour z < 1, diverge pour z > 1.")
    print(f"  ⟹ Power law μ*(1+z)^α est meilleure approximation sur [0,3].")

    # ─────────────────────────────────────────────────────────────────
    # IMPLICATIONS ROUTE D
    # ─────────────────────────────────────────────────────────────────
    print()
    print(sep)
    print("IMPLICATIONS ROUTE D : w(z) = -1 + s × exp(-μ*/μ(z))")
    print(sep)

    # Route D avec μ(z) de Bianchi I
    w_bianchi_arr = np.array([
        w_route_D(z_arr[i], mu_arr[i]) for i in range(len(z_arr))
    ])

    # Route D avec hypothèse COND μ = μ*/(1+z)
    w_cond_arr = np.array([w_route_D_cond(z) for z in z_arr])

    # Paramètres CPL analytiques (α de Bianchi I)
    w0_bianchi, wa_bianchi = cpl_params_analytical(alpha_fit)

    # Paramètres CPL par ajustement numérique
    w0_cpl_bi, wa_cpl_bi = cpl_fit(z_arr, w_bianchi_arr)
    w0_cpl_co, wa_cpl_co = cpl_fit(z_arr, w_cond_arr)

    print()
    print("  w(z) Route D pour différents μ(z) :")
    print(f"  {'z':>5} | {'w_Bianchi':>11} | {'w_COND':>9} | {'Δw':>9}")
    print(f"  {'-'*5}-+-{'-'*11}-+-{'-'*9}-+-{'-'*9}")
    for ztgt in [0.0, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0]:
        idx = np.argmin(np.abs(z_arr - ztgt))
        wb  = w_bianchi_arr[idx]
        wc  = w_cond_arr[idx]
        print(f"  {ztgt:5.2f} | {wb:+11.6f} | {wc:+9.6f} | {wb-wc:+9.6f}")

    print()
    print("  Paramètres CPL w(z) = w₀ + w_a × z/(1+z) :")
    print()
    print(f"  {'Source':25} | {'w₀':>9} | {'w_a':>9} | {'w₀+w_a':>8} | {'pull_w₀':>9} | {'pull_wa':>9}")
    print(f"  {'-'*25}-+-{'-'*9}-+-{'-'*9}-+-{'-'*8}-+-{'-'*9}-+-{'-'*9}")

    rows = [
        ("Bianchi I [analyt.]",  w0_bianchi, wa_bianchi,
         (w0_bianchi - W0_DESI)/W0_ERR, (wa_bianchi - WA_DESI)/WA_ERR),
        ("Bianchi I [fit CPL]",  w0_cpl_bi,  wa_cpl_bi,
         (w0_cpl_bi - W0_DESI)/W0_ERR, (wa_cpl_bi - WA_DESI)/WA_ERR),
        ("COND μ=μ*/(1+z) [fit]", w0_cpl_co, wa_cpl_co,
         (w0_cpl_co - W0_DESI)/W0_ERR, (wa_cpl_co - WA_DESI)/WA_ERR),
        ("DESI DR2",             W0_DESI, WA_DESI, 0.0, 0.0),
    ]
    for name, w0, wa, t0, ta in rows:
        print(f"  {name:25} | {w0:+9.4f} | {wa:+9.4f} | {w0+wa:+8.4f} | {t0:+9.2f}σ | {ta:+9.2f}σ")

    # ─────────────────────────────────────────────────────────────────
    # FORMULAIRE ANALYTIQUE
    # ─────────────────────────────────────────────────────────────────
    print()
    print(bar)
    print("  FORMULE ANALYTIQUE EXACTE (Bianchi I, CPL à z=0) :")
    print()
    print("  w₀ = -1 + s/e  =  -1 + 1/(2e)              [DER-PARTIAL, T1+L2]")
    print(f"     = {-1 + S/E:.6f}  (0.27σ DESI)")
    print()
    print("  w_a [Bianchi I] = s × α / e")
    print(f"     avec α = -1/(μ* × H_iso_fac(μ*)) = {alpha_0:.4f}")
    print(f"     = s × α/e = {S * alpha_0 / E:.6f}")
    print()
    print("  Ou de façon équivalente :")
    print("     w_a [Bianchi I] = s/(e × μ* × |H_iso_fac(μ*)|)")
    print(f"     = {S:.1f} / ({E:.4f} × {MU_STAR:.1f} × {abs(Hf_star):.6f})")
    print(f"     = {S / (E * MU_STAR * abs(Hf_star)):.6f}")
    print()
    print("  Comparer avec COND :  w_a [COND] = -s/e")
    print(f"     = {-S/E:.6f}  (α_COND = -1, μ décroissant vers passé)")

    # ─────────────────────────────────────────────────────────────────
    # VERDICT SCÉNARIOS
    # ─────────────────────────────────────────────────────────────────
    print()
    print(sep)
    print("VERDICT — Scénarios S1–S4")
    print(sep)
    print()
    print(f"  L'ODE Bianchi I donne μ(z) croissant avec z :")
    print(f"    μ(z) ≈ μ* × (1+z)^{alpha_fit:.3f}   avec α ≈ {alpha_fit:.3f} > 0")
    print()
    print(f"  S1 (μ ∝ 1/(1+z)^α) : α ≠ −1, α ≈ +{alpha_fit:.3f}  [RÉFUTÉ]")
    print(f"       La COND Route D supposait α = −1 ; Bianchi I donne α ≈ +{alpha_fit:.3f}.")
    print()
    print(f"  S2 (μ croît avec z — MAUVAIS signe) : CONFIRMÉ  ◀ LE SCÉNARIO RÉALISÉ")
    print(f"       μ(z) croît vers le passé ⟹ w_a = +{wa_bianchi:.3f} > 0")
    print(f"       DESI DR2 requiert w_a ≈ −0.72 < 0 : tension {abs(wa_bianchi - WA_DESI)/WA_ERR:.1f}σ")
    print()
    print(f"  S3 (μ = cst) : RÉFUTÉ  (μ̇ = 0 n'est pas solution de G_00 = 8πG ρ)")
    print()
    print(f"  S4 (incohérent) : RÉFUTÉ  (ODE bien posée, solution unique)")
    print()
    print("  CONCLUSION ÉPISTÉMIQUE [DER] :")
    print()
    print("  La dérivation Bianchi I prouve que μ_eff(z) CROÎT vers le passé")
    print("  avec α ≈ 1.16, ce qui force w_a(Route D) ≈ +0.213 > 0.")
    print()
    print("  La contrainte w_a [COND] = −s/e = −0.184 supposait μ DÉCROISSANT")
    print("  (α = −1), ce qui est cinématiquement impossible en Bianchi I.")
    print()
    print("  Route D peut expliquer w₀ = −0.816 (0.27σ DESI) [DER-PARTIAL]")
    print("  mais prédit w_a ≈ +0.213 (signe opposé DESI) [DER].")
    print()
    print("  PT prédit donc pour Route D :")
    print(f"    w₀ = −1 + 1/(2e) ≈ −0.816   [DER-PARTIAL]")
    print(f"    w_a = +s×α/e ≈ +0.213        [DER : Bianchi I Approche 1]")
    print(f"    w₀ + w_a ≈ −0.632            (≠ −1 : pas de contrainte exacte)")
    print()
    print("  Si confirmé par DESI DR3 (w_a < 0) : Route D est réfutée.")
    print("  Si w_a > 0 confirmé : Route D est validée par Bianchi I.")

    # ─────────────────────────────────────────────────────────────────
    # q₀ en Bianchi I
    # ─────────────────────────────────────────────────────────────────
    print()
    print(sep)
    print("PARAMÈTRE DE DÉCÉLÉRATION q₀")
    print(sep)

    # Route D avec w₀ = -0.816
    w0_RD = -1.0 + S / E
    q0_RD = q0_flrw(w0_RD)
    q0_PT = -0.530   # prédiction PT (ch13)
    q0_LCDM = q0_flrw(-1.0)

    print()
    print(f"  q₀ (ΛCDM, w₀=-1)       = {q0_LCDM:.4f}")
    print(f"  q₀ (Route D, w₀={w0_RD:.3f}) = {q0_RD:.4f}")
    print(f"  q₀ (PT ch13)             = {q0_PT:.4f}")
    print()
    print(f"  Tension interne :  q₀(Route D) - q₀(PT) = {q0_RD - q0_PT:+.4f}")
    print(f"  Route D génère q₀ ≈ {q0_RD:.3f}, en tension Δ = {q0_RD - q0_PT:+.4f} avec q₀_PT = -0.530.")
    print(f"  Cette tension est interne à PT ; elle ne dépend pas de DESI.")

    # ─────────────────────────────────────────────────────────────────
    # RÉSUMÉ POUR LaTeX
    # ─────────────────────────────────────────────────────────────────
    print()
    print(sep)
    print("RÉSUMÉ POUR MISE À JOUR LaTeX (ch20f_dark_energy.tex — Route D)")
    print(sep)
    print()
    print("  Résultats à insérer dans ch20f, section Route D :")
    print()
    print("  [DER, ODE Bianchi I, 2026-04-24]")
    print()
    print("  dμ_eff/dz = -1/[(1+z) × H_iso_fac(μ)]   (cinématique pure)")
    print(f"  H_iso_fac(μ*) = {Hf_star:.6f}")
    print(f"  α ≡ -1/(μ* × H_iso_fac(μ*)) = {alpha_0:.4f}  (exposant power-law local)")
    print()
    print("  μ_eff(z) ≈ μ* × (1+z)^α  avec α ≈ 1.16")
    print("  ⟹ μ_eff croît vers le passé (OPPOSÉ à l'hypothèse COND α = -1)")
    print()
    print("  Route D CPL dérivé :")
    print(f"  w₀ = -1 + s/e = {w0_RD:.6f}   [DER-PARTIAL : T1 + L2]")
    print(f"  w_a = +s×α/e  = +{S * alpha_0 / E:.6f}  [DER : Bianchi I]")
    print()
    print("  Prédiction falsifiable DESI DR3 :")
    print("  w_a > 0 (Route D Bianchi I)  vs  w_a ≈ -0.72 (DESI DR2).")
    print("  Si DESI DR3 confirme w_a < 0 : Route D est réfutée.")
    print()

    print(sep)
    print("Script terminé.")
    print(sep)


if __name__ == "__main__":
    main()
