"""
ÉTAPE E5a : test sur géométrie hyperbolique cuspidale de Σ_pers.

CONTEXTE (ch37b §17.1-17.2, théorèmes Z1, Z2) :
- Σ_pers : surface hyperbolique 4D, point isolé hyperelliptique de M_{0,30}^hyp
- Genus 14 (preprint A_PT)
- Régime asymptotique cusp : métrique
    ds²_R = (1/y²) [3 dy² + (4/225) Σ_p dx_p²]
  forme demi-espace hyperbolique H^4 de rayon r_hyp = √3
- Courbure scalaire : R_4D → -4 quand μ → ∞ (vs Bianchi I plat E4'')
- K_∞ = -1/3 = -1/|P|

Coordonnées :
  z ∈ (-1, 1) (PT_GeoFlow), μ(z) = μ*/(1-z²) = 15/(1-z²)
  y = 1/(1-z) coordonnée Poincaré, y ∈ (1/2, ∞) pour z ∈ (-1, 1)
  Cusp à y → ∞ (μ → ∞)
  Point fixe à z=0, μ=15, y=1

CALCUL : intégrer ⟨Δq^n⟩ sur le régime ASYMPTOTIQUE (limite cusp pur)
pour voir si R prend une valeur canonique (e.g. 17/16).
"""

import mpmath as mp

mp.mp.dps = 50

# =====================================================================
# §1 — Régime cusp pur (asymptotique μ → ∞)
# =====================================================================
print("="*70)
print("  §1. Régime cusp pur asymptotique (μ → ∞)")
print("="*70)
print("""
  En μ → ∞ :
    Δq(μ) ≈ 1/μ + 1/(2μ²) + O(1/μ³)
    Δq(μ)^n ≈ (1/μ)^n  (leading)

  Métrique cusp (ch37b Z2) en coordonnée y :
    ds²_R = (1/y²) [3 dy² + (4/225) Σ dx_p²]
    √g_R = (1/y²)² · √(3 · (4/225)³) = ... (dim 4)
    Volume sur cusp [y₀, ∞] :
      Vol_{cusp} = ∫√g_R dy d³x = const · ∫dy / y⁴ = const / (3 y₀³)
    Moments :
      ⟨Δq^n⟩ = ∫Δq^n √g_R dy d³x / Vol_{cusp}

  Avec μ(y) ≈ μ*·y/2 pour y grand (à corriger) :
    Δq(y) ≈ 2/(μ*·y)
    Δq^n(y) ≈ (2/μ*)^n · y^{-n}

  ⟨y^{-n}⟩ avec poids 1/y⁴ :
    ∫_{y₀}^∞ y^{-n} · y^{-4} dy / ∫_{y₀}^∞ y^{-4} dy
    = (1/(n+3)) y₀^{-(n+3)} / [(1/3) y₀^{-3}]
    = 3 / [(n+3) y₀^n]

  Donc :
    ⟨Δq²⟩_{cusp} = (2/μ*)² · 3/(5 y₀²)
    ⟨Δq⁴⟩_{cusp} = (2/μ*)⁴ · 3/(7 y₀⁴)

  Ratio R_{cusp} = ⟨Δq⁴⟩ / ⟨Δq²⟩²
    = [(2/μ*)⁴ · 3/(7 y₀⁴)] / [(2/μ*)⁴ · 9/(25 y₀⁴)]
    = (3/7) / (9/25) = 75/63 = 25/21
""")

R_cusp_asymp = mp.mpf(25)/21
print(f"  R_cusp asymptotique = 25/21 = {mp.nstr(R_cusp_asymp, 12)}")
print(f"  vs R_Bianchi-I (E4'') = 1.06208")
print(f"  vs cible 17/16 = {mp.nstr(mp.mpf(17)/16, 12)}")

lambda_H_cusp = 2 * R_cusp_asymp / 17
print(f"\n  λ_H prédit (formule note 16) = 2 · R_cusp / 17 = {mp.nstr(lambda_H_cusp, 12)}")
print(f"  vs PT 1/8 = 0.125")
print(f"  Écart = {mp.nstr(100*(lambda_H_cusp - mp.mpf(1)/8)/(mp.mpf(1)/8), 6)} %")

# =====================================================================
# §2 — Test alternatif : autres formules de fermeture possibles
# =====================================================================
print("\n" + "="*70)
print("  §2. Si R_cusp = 25/21 (hypothèse), quelle formule donne 1/8 ?")
print("="*70)

R = R_cusp_asymp
lambda_H_PT = mp.mpf(1)/8

# Combinaisons structurelles
print(f"\n  Tests sur R = {R}, cible λ_H = 1/8 = {lambda_H_PT}")
formulas = [
    ('R / (p_2+μ*) = R/17', R/17),
    ('R / (2·(p_2+μ*))', R/(2*17)),
    ('R · K_∞ / |P|', R * (mp.mpf(-1)/3) / 3),
    ('R · r_hyp² / (something)', R * 3 / 28.57),  # cherche 28.57 = ?
    ('R / (μ* + p_active² )', R/(mp.mpf(15) + 4)),  # 19?
    ('(R - 1) · |R_4D|', (R-1) * 4),  # exotic
]
for label, val in formulas:
    err = abs(val - lambda_H_PT) / lambda_H_PT * 100
    print(f"    {label:35s} = {mp.nstr(val, 12)}  écart à 1/8 = {mp.nstr(err, 6)} %")

# Si on cherche dénominateur D tel que R/D = 1/8 :
D = R * 8
print(f"\n  Dénominateur requis pour R/D = 1/8 : D = {mp.nstr(D, 12)}")
print(f"  Ce serait D = 8 · 25/21 = 200/21 ≈ 9.524")
print(f"  9.524 = ?  pas évident en constantes PT")

# =====================================================================
# §3 — Calcul plus réaliste : intégration sur intervalle [μ_7, μ_11]
# avec métrique hyperbolique exacte (pas seulement asymptotique)
# =====================================================================
print("\n" + "="*70)
print("  §3. Intégration sur intervalle [μ_7, μ_11] avec métrique cusp")
print("="*70)

# Mais : la coordonnée y = 1/(1-z) avec z² = 1 - μ*/μ requiert μ ≥ μ*
# Pour μ < μ* (notamment μ_7 = 11.63 < 15), z² < 0, y est imaginaire.
# Donc le paramétrage cusp NE COUVRE PAS l'intervalle complet.

print("""
  ATTENTION : le paramétrage cusp y = 1/(1-z) avec z² = 1 - μ*/μ
  ne couvre que μ ≥ μ* = 15. Pour μ < μ* (notamment μ_7 = 11.63),
  z est imaginaire et la métrique cusp n'est pas définie.

  ⇒ L'intervalle [μ_7, μ_11] = [11.63, 17.98] contient BOTH régime
    intérieur (μ < 15) ET régime cusp début (μ > 15).
  ⇒ La géométrie hyperbolique cuspidale s'applique uniquement à la
    partie μ ∈ [15, 17.98], donc demi-intervalle.
""")

# Intégration sur la partie cusp uniquement [μ*, μ_11]
mu_star = mp.mpf(15)
mu_11 = mp.mpf('17.9752246908774')

# Δq sur cette partie
def delta_q(mu):
    return mp.exp(-1/mu) - 1 + 2/mu

# Métrique cusp en variable y
def y_of_mu(mu):
    """μ = μ*/(1-z²), z = 1 - 1/y donc z² = 1 - 2/y + 1/y², μ*(1 - z²) = μ*(2/y - 1/y²)
       μ = μ* y²/(2y - 1)
       Inversement : 2y - 1 = μ* y²/μ ⟹ μ*y² - 2μy + μ = 0
       y = [2μ ± √(4μ² - 4μ*μ)] / (2μ*) = [μ ± √(μ² - μ*μ)] / μ*
       Branche positive : y = [μ + √(μ(μ-μ*))] / μ*"""
    return (mu + mp.sqrt(mu * (mu - mu_star))) / mu_star

print(f"  y(μ=15) = {y_of_mu(15)}")
print(f"  y(μ=17.98) = {y_of_mu(mu_11)}")

# Intégrales sur la métrique hyperbolique :
# Vol = ∫ (1/y²)·√[3 · (4/225)³] dy d³x avec d³x = (2π)³
# √[3 · (4/225)³] = √3 · (8/(225·15)) = ...
# Note : (4/225)^(3/2) = (4)^(3/2)/(225)^(3/2) = 8/3375

def sqrt_g_cusp(y):
    """√g_cusp = √[3 · (4/225)³] / y^4 = (8√3 / 3375) / y^4"""
    return 8 * mp.sqrt(3) / mp.mpf(3375) / y**4

def integrand_n(mu, n):
    y = y_of_mu(mu)
    # Jacobian : dy/dμ = ?
    # y = (μ + √(μ(μ-μ*))) / μ*
    # dy/dμ = (1 + (2μ-μ*)/(2√(μ(μ-μ*)))) / μ*
    if mu > mu_star:
        denom = 2 * mp.sqrt(mu * (mu - mu_star))
        dydmu = (1 + (2*mu - mu_star)/denom) / mu_star
    else:
        return 0  # singularity at mu=mu*
    return delta_q(mu)**n * sqrt_g_cusp(y) * dydmu * (2*mp.pi)**3

# Calcul pour μ ∈ (μ*, μ_11]
print("\n  Intégrales sur région cusp seulement (μ ∈ [μ*+ε, μ_11]) :")
eps = mp.mpf('0.01')
vol_cusp = mp.quad(lambda mu: integrand_n(mu, 0), [mu_star + eps, mu_11])
dq2_cusp = mp.quad(lambda mu: integrand_n(mu, 2), [mu_star + eps, mu_11])
dq4_cusp = mp.quad(lambda mu: integrand_n(mu, 4), [mu_star + eps, mu_11])

print(f"  Vol = {mp.nstr(vol_cusp, 12)}")
print(f"  ⟨Δq²⟩·Vol = {mp.nstr(dq2_cusp, 12)}")
print(f"  ⟨Δq⁴⟩·Vol = {mp.nstr(dq4_cusp, 12)}")
print(f"  ⟨Δq²⟩ = {mp.nstr(dq2_cusp/vol_cusp, 12)}")
print(f"  ⟨Δq⁴⟩ = {mp.nstr(dq4_cusp/vol_cusp, 12)}")
R_partial = (dq4_cusp/vol_cusp) / (dq2_cusp/vol_cusp)**2
print(f"\n  R_cusp partial = {mp.nstr(R_partial, 12)}")
print(f"  vs R_Bianchi-I = 1.06208")
print(f"  vs 17/16 = 1.0625")
print(f"  vs R_cusp asymp = 25/21 = {mp.nstr(R_cusp_asymp, 12)}")

# Test formule
lambda_H_partial = 2 * R_partial / 17
print(f"\n  λ_H = 2R/17 = {mp.nstr(lambda_H_partial, 12)}")
print(f"  vs 1/8 = 0.125, écart = {mp.nstr(100*(lambda_H_partial - mp.mpf(1)/8)/(mp.mpf(1)/8), 6)} %")

# =====================================================================
# §4 — Verdict E5a (session)
# =====================================================================
print("\n" + "="*70)
print("  §4. VERDICT E5a (cette session)")
print("="*70)
print("""
  GÉOMÉTRIE HYPERBOLIQUE CUSPIDALE :
    - Métrique exacte : ds² = (1/y²)·[3dy² + (4/225)Σdx²] (ch37b Z2)
    - Définie pour μ ≥ μ*, pas pour μ < μ*
    - Régime asymptotique pur : R_cusp → 25/21 = 1.190 (calcul ci-dessus)
    - Régime intermédiaire [μ*, μ_11] : R calculé numériquement

  COMPARAISON :
    R_Bianchi-I (E4'')         = 1.06208 → λ_H = 0.12495 (écart 0.04%)
    R_cusp asymp pur            = 25/21 ≈ 1.190 → λ_H = 0.1401 (écart 12%)
    R_cusp partial [μ*, μ_11]   = (voir ci-dessus)
    Cible 17/16 = 1.0625        → λ_H = 0.125 exact

  OBSERVATION CRUCIALE :
    Le régime cusp ne donne PAS 17/16. Il donne 25/21.
    Le régime Bianchi I tronqué (E4'') donne 1.062 ≈ 17/16 à 0.04 %.

    Donc :
    - Soit l'identification correcte est R_Bianchi-I = 17/16 (à correction
      structurelle près de 0.04%, possiblement géométrique)
    - Soit la formule λ_H = 2R/17 est CONTEXTE-DÉPENDANTE (différentes
      formules sur différents régimes)
    - Soit la coïncidence 17/8 est numérologique et le bon résultat est
      ailleurs.

  STATUS K4 après E5a :
    L'hypothèse H2 (géométrie hyperbolique complète donne R = 17/16 exact)
    est INFIRMÉE : le régime cusp asymp donne 25/21, pas 17/16.

    Cela ne ferme pas K4. Cela renforce l'idée que :
    - La formule λ_H = 2R/17 est spécifique à la géométrie Bianchi I × T³
    - L'écart 0.04 % résiduel n'est pas dû à la géométrie cuspidale
    - L'identification du préfacteur 1/17 reste structurellement intéressante
      mais sa dérivation propre nécessite plus de travail.

  K4 reste à [CONJ STRUCTURELLE FORTE] avec formule candidate sur
  Bianchi I × T³ uniquement.
""")
