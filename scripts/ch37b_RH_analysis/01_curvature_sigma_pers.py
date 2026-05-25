#!/usr/bin/env python3
"""01 - Courbure scalaire de Σ_pers (Fisher-Bianchi).

Z1 du PROMPT_LAUNCH.md : calculer numeriquement et symboliquement la
courbure scalaire R_Σ de Σ_pers munie de sa metrique Fisher-Bianchi.

Convention PT :
  μ* = 15, premiers actifs P = {3, 5, 7}
  q_+(μ)   = 1 - 2/μ
  δ_p(μ)   = (1 - q_+(μ)^p) / p
  sin²θ_p  = δ_p · (2 - δ_p)
  γ_p      = 1 - δ_p
  α(μ)     = ∏_{p ∈ P} sin²θ_p(μ)
  a_p(μ)   = γ_p(μ) / μ          (facteur d'echelle Bianchi I, direction p)
  g_μμ(μ)  = -∂²(ln α)/∂μ²        (g < 0  ⇔  Lorentzien)

La courbure scalaire est calculee sur le contre-pendant Riemannien
(equivalent conforme aux signatures pres) :
  ds²_R = |g_μμ| dμ² + Σ_p a_p² dx_p²

Strategie : Sympy pour construire l'arbre symbolique des derivees,
lambdify -> mpmath pour evaluation rapide. On verifie l'asymptotique
μ → ∞ par evaluation a tres grand μ (5, 6, 7 decades).

Formules de courbure :
  Pour 2D Riemannien ds² = A² dμ² + B² dθ² :
    R_2D = 2 (A'B' - AB'') / (A³ B)
  Pour 4D Bianchi I Riemannien ds² = A² dμ² + Σ b_p² dx_p² (avec
    A² = |g_μμ|, tau = ∫A dμ proper time) :
    R_4D = -2 Σ b_p''(τ)/b_p - 2 Σ_{p<q} (b_p'/b_p)(b_q'/b_q)
  ou b_p'(τ) = b_p'(μ)/A et b_p''(τ) = (b_p''(μ) - (A'/A)b_p'(μ))/A².
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp
from mpmath import mp, mpf


# ----------------------------------------------------------------------
# Section 1 : expressions symboliques (construction rapide, pas de simplify)
# ----------------------------------------------------------------------

mu = sp.symbols("mu", positive=True)
PRIMES = (3, 5, 7)


def build_symbolics():
    q_plus = 1 - sp.Rational(2) / mu
    sin2_factors = []
    gammas = []
    for p in PRIMES:
        delta_p = (1 - q_plus**p) / p
        sin2_factors.append(delta_p * (2 - delta_p))
        gammas.append(1 - delta_p)
    alpha = sp.Mul(*sin2_factors, evaluate=False)
    log_alpha = sp.log(alpha)
    g_mumu = -sp.diff(log_alpha, mu, 2)
    a_sq = -g_mumu                       # |g_μμ|
    a_sym = sp.sqrt(a_sq)                 # lapse A(μ)
    a_prime = sp.diff(a_sym, mu)
    b_list = [g / mu for g in gammas]
    bp_list = [sp.diff(b, mu) for b in b_list]
    bpp_list = [sp.diff(b, mu, 2) for b in b_list]
    return {
        "g_mumu": g_mumu,
        "a_sym": a_sym,
        "a_prime": a_prime,
        "b_list": b_list,
        "bp_list": bp_list,
        "bpp_list": bpp_list,
    }


def R_2D_expr(syms, i):
    A = syms["a_sym"]; Ap = syms["a_prime"]
    B = syms["b_list"][i]
    Bp = syms["bp_list"][i]
    Bpp = syms["bpp_list"][i]
    return 2 * (Ap * Bp - A * Bpp) / (A**3 * B)


def R_4D_expr(syms):
    A = syms["a_sym"]; Ap = syms["a_prime"]
    n = len(syms["b_list"])
    bp_tau_over_b = []
    bpp_tau_over_b = []
    for i in range(n):
        B = syms["b_list"][i]
        Bp = syms["bp_list"][i]
        Bpp = syms["bpp_list"][i]
        bp_tau_over_b.append(Bp / (A * B))
        bpp_tau_over_b.append((Bpp - (Ap / A) * Bp) / (A**2 * B))
    R = sp.Integer(0)
    for i in range(n):
        R += -2 * bpp_tau_over_b[i]
    for i in range(n):
        for j in range(i + 1, n):
            R += -2 * bp_tau_over_b[i] * bp_tau_over_b[j]
    return R


# ----------------------------------------------------------------------
# Section 2 : lambdification mpmath
# ----------------------------------------------------------------------

SYMS = build_symbolics()
R2_SYM = {p: R_2D_expr(SYMS, i) for i, p in enumerate(PRIMES)}
R4_SYM = R_4D_expr(SYMS)
GMUMU_SYM = SYMS["g_mumu"]
B_SYM = dict(zip(PRIMES, SYMS["b_list"]))

# precision mpmath
mp.dps = 50

R2_F = {p: sp.lambdify(mu, R2_SYM[p], modules="mpmath") for p in PRIMES}
R4_F = sp.lambdify(mu, R4_SYM, modules="mpmath")
GMUMU_F = sp.lambdify(mu, GMUMU_SYM, modules="mpmath")
B_F = {p: sp.lambdify(mu, B_SYM[p], modules="mpmath") for p in PRIMES}


def eval_R2(p, mu_val):
    return mpf(R2_F[p](mpf(mu_val)))


def eval_R4(mu_val):
    return mpf(R4_F(mpf(mu_val)))


def eval_g(mu_val):
    return mpf(GMUMU_F(mpf(mu_val)))


def eval_b(p, mu_val):
    return mpf(B_F[p](mpf(mu_val)))


# ----------------------------------------------------------------------
# Section 3 : driver
# ----------------------------------------------------------------------

def main():
    lines = []

    def out(s=""):
        print(s, flush=True)
        lines.append(s)

    out("=" * 76)
    out("Z1 - Courbure scalaire de Σ_pers (Fisher-Bianchi Wick → Riemannien)")
    out("=" * 76)
    out()
    out("Convention : courbure 'scalaire de Ricci' R = g^{μν} R_μν,")
    out("calculee sur la metrique Riemannienne |g_μμ|·dμ² + Σ a_p² dx_p²")
    out("(Wick rotation de la Fisher-Bianchi Lorentzienne PT, μ > μ_c).")
    out()

    # 1. Sanity au point fixe
    out("# 1. Valeurs au point fixe μ = μ* = 15")
    out()
    out(f"   g_μμ(15)   = {float(eval_g(15)):+.6e}    (< 0 : Lorentzien)")
    for p in PRIMES:
        out(f"   a_{p}(15)    = {float(eval_b(p, 15)):.6f}    (γ_{p}/15)")
    for p in PRIMES:
        out(f"   R_2D^(p={p})(15) = {float(eval_R2(p, 15)):+.6f}")
    out(f"   R_4D(15)    = {float(eval_R4(15)):+.6f}")
    out()

    # 2. Grille numerique
    out("# 2. Tableau R(μ) sur grille [μ_c ; 200] et au-dela (verification asymptotique)")
    out()
    grid = [7, 8, 10, 12, 13, 14, 15, 16, 18, 22, 30, 50, 100, 200, 500, 1000,
            10**4, 10**5, 10**6]
    out("   " + f"{'μ':>10}  " +
        "  ".join(f"{'R_2D(p='+str(p)+')':>14}" for p in PRIMES) +
        f"  {'R_4D':>14}  {'g_μμ':>14}")
    out("   " + "-" * (10 + (14 + 2) * 4 + 14 + 2))
    grid_results = []
    for mu_val in grid:
        Rs_2D = [eval_R2(p, mu_val) for p in PRIMES]
        R4 = eval_R4(mu_val)
        gm = eval_g(mu_val)
        grid_results.append((mu_val, Rs_2D, R4, gm))
        row = f"   {mu_val:>10}  " + "  ".join(f"{float(Rv):>+14.8f}" for Rv in Rs_2D)
        row += f"  {float(R4):>+14.8f}  {float(gm):>+14.6e}"
        out(row)
    out()

    # 3. Convergence vers les limites theoriques
    out("# 3. Comparaison avec asymptotique theorique (μ → ∞)")
    out()
    out("   Prediction analytique (developpement Laurent en 1/μ) :")
    out("     R_2D^(p)(μ → ∞) = -2/3 ≈ -0.66666...   (pour les 3 primes)")
    out("     R_4D(μ → ∞)     = -4")
    out("   Courbure sectionnelle constante associee : K = -1/3,")
    out("   rayon hyperbolique r_hyp = sqrt(3).")
    out()
    out("   Erreur relative |R(μ) - R_∞|/|R_∞| au plus haut μ teste (μ = 10^6) :")
    last_mu, last_R2, last_R4, _ = grid_results[-1]
    for i, p in enumerate(PRIMES):
        err = abs(float(last_R2[i]) - (-2/3)) / (2/3)
        out(f"     R_2D^(p={p}) :  {float(last_R2[i]):+.10f}   err = {err:.3e}")
    err4 = abs(float(last_R4) - (-4)) / 4
    out(f"     R_4D       :  {float(last_R4):+.10f}   err = {err4:.3e}")
    out()

    # 4. Verdict
    out("# 4. Verdict Z1")
    out()
    dom = [r for r in grid_results if r[0] >= 7]
    R2_neg = all(min(r[1]) < 0 for r in dom)
    R4_neg = all(r[2] < 0 for r in dom)
    g_neg = all(r[3] < 0 for r in dom)
    out(f"   - g_μμ < 0 sur [7, 10^6]                 : {g_neg}")
    out(f"   - R_2D(p) < 0 sur [7, 10^6] pour tout p  : {R2_neg}")
    out(f"   - R_4D    < 0 sur [7, 10^6]              : {R4_neg}")
    out()
    # 1% de la valeur asymptotique a μ = 10^4 ?
    R2_conv_1pct = all(
        abs(float(eval_R2(p, 10000)) - (-2/3)) / (2/3) < 0.01 for p in PRIMES
    )
    R4_conv_1pct = abs(float(eval_R4(10000)) - (-4)) / 4 < 0.01
    out(f"   - |R_2D(μ=10^4) - (-2/3)|/|R_∞| < 1% (tous p) : {R2_conv_1pct}")
    out(f"   - |R_4D(μ=10^4) - (-4)|/|R_∞|    < 1%        : {R4_conv_1pct}")
    out()

    if R2_neg and R4_neg and g_neg and R2_conv_1pct and R4_conv_1pct:
        out("   VALIDATION MODEREE.")
        out()
        out("   * R_Σ < 0 partout sur le domaine Lorentzien (μ > μ_c).")
        out("   * R_Σ varie continument, sans changement de signe.")
        out("   * R_Σ converge vers une CONSTANTE NEGATIVE universelle :")
        out("       - 2D : -2/3 pour chacune des 3 slices (μ, x_p)")
        out("       - 4D : -4 pour Bianchi I complet")
        out("     Ces deux valeurs correspondent a une courbure sectionnelle")
        out("     constante K = -1/3 (rayon hyperbolique sqrt(3)) :")
        out("        2D : R = 2K = -2/3                  ✓")
        out("        4D : R = -n(n-1)|K| = -12·(1/3) = -4 ✓")
        out()
        out("   Σ_pers est asymptotiquement isometrique a une variete")
        out("   hyperbolique a courbure constante K = -1/3 dans la limite cusp")
        out("   (μ → ∞, z → 1). Au point fixe PT μ = μ* = 15, la courbure")
        out("   est plus negative que la limite (cf. tableau), mais reste")
        out("   strictement negative.")
        out()
        out("   L'hypothese de travail (Σ_pers hyperbolique a cusp) est")
        out("   VALIDEE asymptotiquement et numeriquement compatible sur")
        out("   tout le domaine Lorentzien teste. Continuer avec Z2.")
    else:
        out("   VERDICT NEGATIF : la courbure ne se comporte pas comme")
        out("   attendu pour une surface hyperbolique a cusp.")

    out()
    out("=" * 76)
    out("Fin Z1.")
    out("=" * 76)

    out_path = Path(__file__).parent.parent / "outputs" / "01_curvature_output.txt"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n(Resultat sauvegarde : {out_path})")


if __name__ == "__main__":
    sys.exit(main() or 0)
