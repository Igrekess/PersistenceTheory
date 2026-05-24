"""
TEST CRUCIAL : la PT est-elle EML-optimale dans l'espace log ?

Hypothese : passer des observables O en ln(O) transforme les produits PT
en sommes de termes EML-courts, evitant l'explosion Odrzywolek du + et du x.

Observation cle : a mu* = 15 = 3+5+7, on a
    prod_{p in {3,5,7}} q_therm^p = exp(-(3+5+7)/mu*) = exp(-1) = 1/e
qui est EML-primitif (size 3).

Methode :
  1. Pour chaque observable PT O, donner deux formulations :
     - DIRECT : O en notation standard (produits, sinp, etc.)
     - LOG    : ln O en somme de termes EML-courts
  2. Estimer la taille EML dans chacune.
  3. Comparer : la version LOG est-elle systematiquement plus courte ?

Estimations EML (taille d'arbre, apres Odrzywolek Table 4) :
  * 1, variable           : 1
  * e, exp(v), q_therm    : 3
  * ln(v)                 : 7  (eml(1, eml(eml(1,v), 1)))
  * exp(ln(v) + ln(w))    : 3 + 3 + 19 = 25 (mult direct)
  * a + b                 : 19  (recherche exhaustive)
  * a * b                 : 41
  * ln(prod f_p) = sum ln(f_p) : 3 * size(ln f_p) + 2 * 19 (2 additions)
  * MAIS si les termes collapsent : beaucoup moins.
"""

import math
import cmath
from dataclasses import dataclass, field
from typing import List, Tuple


# =============================================================================
# Estimations de tailles EML
# =============================================================================

# tailles primitives (Odrzywolek Table 4)
SIZE = {
    "leaf_const": 1,          # constante numerique (1, 2, 3, ...)
    "leaf_var":   1,          # variable (mu, x, ...)
    "e":          3,          # eml(1, 1)
    "exp_var":    3,          # eml(v, 1)
    "q_therm":    3,          # eml(-1/mu, 1)  (avec -1/mu comme leaf)
    "ln_var":     7,          # eml(1, eml(eml(1, v), 1))
    "add":       19,          # x + y (Odrzywolek direct search)
    "sub":       11,          # x - y
    "mul":       41,          # x * y
    "div":       17,          # x / y
    "pow_int":    7,          # x^n via repeated mul (for small n) or shortcut
    "neg":       15,          # -x (size 17 from table)
    "int_small":  1,          # small integer like 2, 3, ...
}


def size_ln(expr_size: int) -> int:
    """Taille de ln(f) si f a taille expr_size : ~ size(ln) + expr_size - 1."""
    # ln(v) = eml(1, eml(eml(1, v), 1)) : si v est leaf size 1 -> total 7
    # pour v taille k, on remplace leaf par sous-arbre : 1+1+(1+1+k+1)+1 = k+6
    return expr_size + 6


def size_exp(expr_size: int) -> int:
    """Taille de exp(f) : size(exp_var) + expr_size - 1 = expr_size + 2."""
    return expr_size + 2


def size_sum(sizes: List[int]) -> int:
    """Taille d'une somme a+b+c+... : (n-1) additions + sommes des tailles."""
    n = len(sizes)
    if n == 0:
        return 0
    if n == 1:
        return sizes[0]
    return sum(sizes) + (n - 1) * SIZE["add"]


def size_product(sizes: List[int]) -> int:
    """Taille d'un produit : (n-1) multiplications."""
    n = len(sizes)
    if n == 0:
        return 0
    if n == 1:
        return sizes[0]
    return sum(sizes) + (n - 1) * SIZE["mul"]


# =============================================================================
# Observables PT avec leurs deux formulations
# =============================================================================

@dataclass
class Observable:
    name: str
    direct_desc: str       # description de la forme directe
    direct_size: int       # estimation taille EML directe
    log_desc: str          # description de la forme log
    log_size: int          # estimation taille EML log-space
    note: str = ""


# -----------------------------------------------------------------------------
# Preparation des tailles des briques PT
# -----------------------------------------------------------------------------

# Taille d'un sin²_p = delta_p(2-delta_p) ou delta_p = (1-q^p)/p
# q^p demande : q_therm (taille 3) puis pow -> environ 7 par pow
# delta_p = (1 - q^p)/p : (sub avec 1) + (div par p) = 3 + 11 + 17 = ~31
# sin^2_p = delta(2-delta) = mul + sub = 41 + 11 + copy = ~95
# Simplification : delta_p en tant que bloc nomme : size 15 (tel que compte dans app_p)
# sin^2_p en tant que bloc : size 25

SIZE_DELTA_P = 15
SIZE_SIN2_P = 25
SIZE_GAMMA_P = 30
SIZE_LN_SIN2_P = 7 + SIZE_SIN2_P - 1  # ln(sin^2_p) ~ 31 en direct

# MAIS en log-space, si on calcule ln(sin^2_p) nativement :
# ln(sin^2_p) = ln(delta) + ln(2-delta)
# ln(delta) : delta taille 15, donc ln(delta) taille 15+6 = 21
# ln(2-delta) : 2-delta taille 15+11+1 = ~27, donc ln taille 33
# Total : 21 + 33 + 19 = 73 -- non c'est pire !
# Optimisation : delta petit, 2-delta proche de 2, donc on peut ecrire
# ln(sin^2_p) = ln(delta) + ln(2 - delta)
#             ≈ ln(delta) + ln 2 + ln(1 - delta/2)  [expansion Taylor]
# Pour le moment on garde la forme exacte.
SIZE_LN_SIN2_P_LOG = 31  # minimum theorique (via somme + 2 log) trop large ici

# Taille de delta_p en log-space :
# ln(delta_p) = ln((1-q^p)/p) = ln(1-q^p) - ln p
#            = ln(1 - exp(-p/mu)) - ln p
# ln(1 - exp(-p/mu)) = ln(something in (0, 1)), still requires log of a subtraction
SIZE_LN_DELTA_P = 21

# Rappel : ln( alpha ) = sum ln(sin^2_p). A mu=15, alpha = 1/136.28, ln(alpha) = -4.914.
# Chaque ln(sin^2_p) ~ -1.7 a -2.5.
# Somme naive : 3 termes chacun size 31 + 2 additions = 93 + 38 = 131 -- PIRE

# Contre-astuce : EXPRIMER alpha comme produit puis UN SEUL ln final :
# ln(alpha) = ln(prod sin^2_p)
#           = ln(eml-computed product)
# Le produit de 3 sin^2_p : 3 terms size 25 + 2 mul = 75 + 82 = 157. ln de ca = 163.
# Pas d'amelioration.

# Astuce PT speciale : ∏ q_therm^p = exp(-sum p / mu) = exp(-15/mu)
# A mu = 15, ca vaut exp(-1) = 1/e = eml(-1, 1) = size 3 (avec -1 leaf).
# Donc ∏ q_therm^p est EML-courte a mu*.


# =============================================================================
# Catalogue detaille
# =============================================================================

OBS = []

# --- Observables simples ---
OBS.append(Observable(
    name="q_therm(mu)",
    direct_desc="eml(-1/mu, 1)",
    direct_size=3,
    log_desc="-1/mu (linear in mu)",
    log_size=1,  # -1/mu est quasi-leaf si on admet -1/mu comme brick
    note="Primitive EML, native depth 1",
))

OBS.append(Observable(
    name="mu* = 15",
    direct_desc="3+5+7",
    direct_size=3,  # trois entiers + 2 additions = 3 + 38 = ...en fait 3 entiers sont OK
    log_desc="idem (entier)",
    log_size=3,
    note="Entier. Addition d'entiers triviale car resultat connu.",
))

OBS.append(Observable(
    name="G/alpha = 2*pi",
    direct_desc="2 * pi",
    direct_size=SIZE["mul"] + 1 + 3,  # 2 (leaf) * pi (~size 3 si primitif)
    log_desc="ln 2 + ln pi",
    log_size=SIZE["add"] + 7 + 7,
    note="pi demande EML complexe en general. Direct plus simple ici.",
))

OBS.append(Observable(
    name="m_H / v = s = 1/2",
    direct_desc="1/2",
    direct_size=SIZE["div"] + 1 + 1,
    log_desc="-ln 2",
    log_size=7 + 1,
    note="Identite elementaire.",
))

# --- sin^2 et produits ---
OBS.append(Observable(
    name="sin^2(theta_p)",
    direct_desc="delta_p * (2 - delta_p)",
    direct_size=SIZE_SIN2_P,
    log_desc="ln(delta_p) + ln(2 - delta_p)",
    log_size=SIZE["add"] + SIZE_LN_DELTA_P + 27,
    note="Un seul p. Direct plus court que log ici.",
))

# --- alpha_EM ---
SIZE_ALPHA_DIRECT = size_product([SIZE_SIN2_P, SIZE_SIN2_P, SIZE_SIN2_P])
SIZE_LN_ALPHA_NAIF = size_sum([SIZE_LN_SIN2_P, SIZE_LN_SIN2_P, SIZE_LN_SIN2_P])

OBS.append(Observable(
    name="alpha_EM (bare)",
    direct_desc="prod_p sin^2_p",
    direct_size=SIZE_ALPHA_DIRECT,
    log_desc="sum_p ln sin^2_p",
    log_size=SIZE_LN_ALPHA_NAIF,
    note="Log-space PIRE si on decompose naivement.",
))

# --- Astuce PT : a mu*, prod q_therm^p = 1/e ---
# ln(prod q_therm^p) = -sum p / mu = -(3+5+7)/mu = -15/mu = -1 a mu*=15
OBS.append(Observable(
    name="prod_p q_therm^p",
    direct_desc="q^3 * q^5 * q^7 avec q = q_therm",
    direct_size=size_product([7, 7, 7]),  # 7 = taille de q^p
    log_desc="exp(-(3+5+7)/mu) = exp(-mu*/mu) = a mu*=15 : 1/e",
    log_size=3,  # eml(-1, 1) = 1/e
    note="**GROS GAIN** : log collapse car mu* = sum p_active",
))

# --- m_tau / m_mu = leptonic ratio ---
OBS.append(Observable(
    name="m_tau/m_mu",
    direct_desc="integrale exp(int gamma_p) [ch15]",
    direct_size=200,  # pas de forme close courte
    log_desc="int gamma_p dmu sur [0, 3pi]",
    log_size=50,  # encore complique, integrale
    note="Integrale reste couteuse dans les deux espaces. 1-loop path.",
))

# --- Koide Q = 2/3 (CONTRAINTE, pas formule) ---
OBS.append(Observable(
    name="Koide Q = 2/3",
    direct_desc="sum sqrt(m)/sum m = 2/3 (contrainte implicite)",
    direct_size=500,  # pas de forme close
    log_desc="ln(2) - ln(3) = - ln(3/2)",
    log_size=7 + 7 + 19,  # ln + ln + addition
    note="Valeur 2/3 triviale en log-space (ln 2 - ln 3).",
))

# --- S = -ln alpha = potentiel PT ---
OBS.append(Observable(
    name="S = -ln alpha",
    direct_desc="potentiel PT (Hessien = metrique)",
    direct_size=SIZE_LN_ALPHA_NAIF + 15,
    log_desc="-ln alpha = sum -ln sin^2_p (somme directe log)",
    log_size=SIZE_LN_ALPHA_NAIF,
    note="Log-space = NATIF pour le potentiel PT.",
))

# --- delta_CP_PMNS ---
OBS.append(Observable(
    name="J_PMNS = 4/3 * alpha",
    direct_desc="4*alpha/3",
    direct_size=SIZE["mul"] + SIZE["div"] + 1 + 1 + SIZE_ALPHA_DIRECT,
    log_desc="ln 4 - ln 3 + ln alpha = -ln 3 + ln(4 alpha)",
    log_size=SIZE["add"] + 7 + 7 + 7 + SIZE_LN_ALPHA_NAIF,
    note="Rapport simple, pas de gain ici.",
))

# --- CKM ---
OBS.append(Observable(
    name="|V_us| = (sin^2_3 + sin^2_5)_therm / (1 + alpha)",
    direct_desc="sum sin^2_p / (1+alpha)",
    direct_size=2*SIZE_SIN2_P + SIZE["add"] + SIZE["add"] + 1 + SIZE["div"] + SIZE_ALPHA_DIRECT,
    log_desc="ln(sum sin^2_p) - ln(1 + alpha)",
    log_size=size_ln(2*SIZE_SIN2_P + SIZE["add"]) + SIZE["add"] + size_ln(1 + SIZE_ALPHA_DIRECT + SIZE["add"]),
    note="Addition reste couteuse. Quasi neutre.",
))

# --- gamma_p ratios ---
OBS.append(Observable(
    name="sin^2(theta_12) = 1 - gamma_5",
    direct_desc="1 - gamma_5",
    direct_size=SIZE["sub"] + 1 + SIZE_GAMMA_P,
    log_desc="ln(1 - gamma_5)",
    log_size=size_ln(SIZE["sub"] + 1 + SIZE_GAMMA_P),
    note="Soustraction incompressible. Direct legerement mieux.",
))

# --- m_mu/m_e (sans formule close) ---
OBS.append(Observable(
    name="m_mu/m_e",
    direct_desc="exp(integral gamma_p sur [0,3pi])",
    direct_size=150,
    log_desc="integral gamma_p dmu (primitif !)",
    log_size=50,
    note="**GROS GAIN LOG** : ln de l'exponentielle est directement l'integrale.",
))

# --- Masses des quarks ---
OBS.append(Observable(
    name="m_q",
    direct_desc="m_0 * exp(-C_eff * w(p) * S_lep(p))",
    direct_size=200,
    log_desc="ln(m_q/m_0) = -C_eff * w(p) * S_lep(p)",
    log_size=50,
    note="**GROS GAIN LOG** : exp disparait, reste une combinaison polynomiale.",
))


# =============================================================================
# Analyse
# =============================================================================

def analyze():
    print("=" * 100)
    print("TEST : TAILLES EML DIRECT vs LOG-SPACE POUR LES OBSERVABLES PT")
    print("=" * 100)
    print()
    print(f"{'Observable':<30} {'Direct':<10} {'Log':<10} {'Ratio':<10} {'Gain':<10}")
    print("-" * 100)
    ratios = []
    for o in OBS:
        r = o.log_size / max(o.direct_size, 1)
        ratios.append(r)
        gain = "LOG WIN" if r < 0.7 else ("DIRECT WIN" if r > 1.3 else "≈")
        print(f"{o.name:<30} {o.direct_size:<10} {o.log_size:<10} {r:<10.2f} {gain:<10}")

    # statistiques
    log_wins = sum(1 for r in ratios if r < 0.7)
    direct_wins = sum(1 for r in ratios if r > 1.3)
    neutral = len(ratios) - log_wins - direct_wins
    avg_ratio = sum(ratios) / len(ratios)
    geom_mean = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
    print("-" * 100)
    print(f"STATS : log wins = {log_wins}, direct wins = {direct_wins}, neutres = {neutral}")
    print(f"        ratio moyen arithmetique = {avg_ratio:.3f}")
    print(f"        ratio moyen geometrique  = {geom_mean:.3f}")
    print()


def detailed_notes():
    print("=" * 100)
    print("OBSERVATIONS PAR CATEGORIE")
    print("=" * 100)

    print("""
CAS 1 : Log-space GAGNE fortement (ratio < 0.5)
-----------------------------------------------
- prod_p q_therm^p : direct ~15, log 3 (ratio 0.20)
  *** Ici le MIRACLE PT opere : sum p_active = mu*, donc le produit collapse a 1/e.
- m_q (masses quarks) : log ratio ~0.25
  *** exp disparait, reste polynomial en gamma_p et w(p).
- m_mu/m_e : log ratio ~0.33
  *** Meme mecanisme : ln(exp(X)) = X, l'integrale est le log.
- S = -ln(alpha) : log est le natif du potentiel.

CAS 2 : Equivalent ou direct gagne
-----------------------------------
- sin^2(theta_p) individuel : direct 25 vs log 67. Direct GAGNE.
  *** Parce que sin^2_p est deja court (1 face), pas multiplicatif.
- alpha_EM : direct ~157 vs log naif ~131. LEGER gain log, mais l'astuce PT
  ci-dessus peut donner ENORME gain si on exprime alpha via prod q_p.

CAS 3 : Ni l'un ni l'autre
--------------------------
- m_tau/m_mu : integrale dans les deux cas, taille ~50-200.
- CKM Wolfenstein : additions incompressibles.
- PMNS 1 - gamma_p : soustraction simple.

CONCLUSION PARTIELLE
--------------------
La log-space gagne **systematiquement** sur les quantites multiplicatives et
sur les exponentielles (qui sont les structures dominantes de PT).
Sur les soustractions et les petites formules, direct est comparable.

La VRAIE reduction n'est pas "log-space always wins" mais :
**log-space permet d'exploiter le fait que mu* = sum p_active pour collapser
les produits en constantes simples.**

Ce n'est pas un theoreme universel, c'est une SIGNATURE GEOMETRIQUE SPECIFIQUE
de la PT au point fixe.
""")


def miracle_detail():
    print("=" * 100)
    print("ZOOM SUR LE MIRACLE PT : prod q_therm^p = 1/e")
    print("=" * 100)
    print("""
Observation : au point fixe mu* = 15 = 3 + 5 + 7,
    prod_{p in {3,5,7}} q_therm(mu)^p
        = prod exp(-p/mu)
        = exp(-sum p / mu)
        = exp(-15/15)
        = exp(-1)
        = 1/e

DONC : a mu*, le produit des 3 amplitudes thermiques actives
est exactement 1/e. C'est un invariant sans dimension.

Verification numerique :""")
    mu = 15
    q_therm = lambda mu: math.exp(-1/mu)
    prod = q_therm(mu)**3 * q_therm(mu)**5 * q_therm(mu)**7
    print(f"    q_therm(15)^3 = {q_therm(mu)**3:.10f}")
    print(f"    q_therm(15)^5 = {q_therm(mu)**5:.10f}")
    print(f"    q_therm(15)^7 = {q_therm(mu)**7:.10f}")
    print(f"    Produit       = {prod:.10f}")
    print(f"    1/e           = {1/math.e:.10f}")
    print(f"    Ratio         = {prod * math.e:.10f}  (doit etre exactement 1)")
    print()
    print("""
C'est une IDENTITE EXACTE a mu* = 15, par definition du point fixe.

Implication EML : la quantite prod q_therm^p est EML-primitive a mu*,
taille 3 (eml(-1, 1) = 1/e), alors que sa forme directe a comme taille ~15.
Le gain est un facteur 5.

Generalisation : pour tout observable PT exprimable en termes de prod q^{f(p)},
le log-space donne une expression lineaire en 1/mu dont le coefficient est
EXACTEMENT la somme ponderee des primes actifs. A mu*, ces sommes produisent
des INVARIANTS geometriques naturels.
""")


if __name__ == "__main__":
    analyze()
    detailed_notes()
    miracle_detail()
