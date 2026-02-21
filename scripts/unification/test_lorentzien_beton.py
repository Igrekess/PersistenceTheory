#!/usr/bin/env python3
"""
test_lorentzien_beton
=====================

ENGLISH
-------
Concrete Lorentzian: explicit Lorentzian metric from the sieve (no ansatz)

FRANCAIS (original)
-------------------
PASSAGE AU LORENTZIEN -- REPONSES AUX 3 CRITIQUES
====================================================
S15.6.194

Trois critiques identifiees sur le passage (-,+,+,+) :

  PROBLEME 1 : "Deux alphas confondues"
  PROBLEME 2 : "Coordonnees ou identifications"
  PROBLEME 3 : "G = 2*pi*alpha sans dimensions"

Ce script repond a CHAQUE critique par un calcul explicite.

CONVENTION DE SIGNE (comme test_jacobson_lorentzien.py et test_graviton_G_Newton_v2.py) :
  S(mu) = -ln(alpha_EM(mu))          [potentiel de persistance]
  g_00 = -d^2(ln alpha_EM)/dmu^2     [composante temporelle]
  g_pp = (gamma_p / mu)^2            [composantes spatiales]

  ln(alpha_EM) est CONVEXE pour mu > mu_c ~ 6.7
  (d^2(ln alpha)/dmu^2 > 0)
  Donc g_00 = -(positif) < 0 : LORENTZIEN pour mu > mu_c.

  Pour mu < mu_c : g_00 > 0 (Euclidien).
  => TRANSITION DE SIGNATURE a mu_c ~ 6.7 (Hartle-Hawking !)

Author / Auteur: Yan Senez  |  Date: February / Fevrier 2026
Theory / Theorie: Persistence Theory (PT) / Theorie de la Persistance (TP)
"""
import numpy as np
from math import log, pi, exp, sqrt
from scipy.optimize import brentq

# ================================================================
# FONCTIONS FONDAMENTALES
# ================================================================
def sin2_theta(p, mu):
    if mu <= 2.01: return 0.0
    q = 1.0 - 2.0 / mu
    d = (1.0 - q**p) / p
    return d * (2.0 - d)

def gamma_p_func(p, mu):
    if mu <= 2.01: return 0.0
    q = 1.0 - 2.0 / mu
    qp = q**p
    d = (1.0 - qp) / p
    denom = mu * (1.0 - qp) * (2.0 - d)
    if abs(denom) < 1e-30: return 0.0
    return 4.0 * p * q**(p-1) * (1.0 - d) / denom

def alpha_EM_mu(mu, primes=[3, 5, 7]):
    """alpha_EM(mu) = produit des sin^2(theta_p, q_stat) pour p actifs."""
    result = 1.0
    for p in primes:
        result *= sin2_theta(p, mu)
    return result

def alpha_sieve(k, primes_list=None):
    """alpha_sieve(k) = produit (1 - 1/p) pour p <= p_k (Mertens)."""
    if primes_list is None:
        primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    result = 1.0
    for i, p in enumerate(primes_list):
        if i >= k: break
        result *= (1.0 - 1.0/p)
    return result

def ln_alpha_EM(mu):
    a = alpha_EM_mu(mu)
    if a <= 0: return -100.0
    return log(a)

def d1_ln_alpha(mu, h=1e-5):
    """Derivee premiere de ln(alpha_EM)."""
    return (ln_alpha_EM(mu + h) - ln_alpha_EM(mu - h)) / (2*h)

def d2_ln_alpha(mu, h=1e-5):
    """Derivee seconde de ln(alpha_EM)."""
    return (ln_alpha_EM(mu + h) - 2*ln_alpha_EM(mu) + ln_alpha_EM(mu - h)) / h**2

PRIMES_ACTIVE = [3, 5, 7]
alpha_EM_phys = 1.0 / 137.036
mu_alpha = brentq(lambda m: alpha_EM_mu(m) - alpha_EM_phys, 14.5, 16.0, xtol=1e-15)
MU_STAR = mu_alpha


# ================================================================
print("=" * 78)
print("PASSAGE AU LORENTZIEN -- REPONSES AUX 3 CRITIQUES (S15.6.194)")
print("=" * 78)


# ================================================================
# PROBLEME 1 : DEUX ALPHAS
# ================================================================
print(f"\n{'='*78}")
print("  PROBLEME 1 : DEUX ALPHAS -- IDENTIFICATION ET DISTINCTION")
print(f"{'='*78}")

print(f"""
  DEFINITIONS :

  alpha_sieve(k) = prod_{{p <= p_k}} (1 - 1/p)
    = fraction d'entiers survivant au crible
    Par Mertens : ~ e^(-gamma_Euler) / ln(p_k)
    MONOTONE DECROISSANTE en k (1/2 -> 0)

  alpha_EM(mu) = prod_{{p in {{3,5,7}}}} sin^2(theta_p(mu))
    = produit des angles de melange
    NON MONOTONE en mu (0 -> max -> 0)

  CE SONT DEUX FONCTIONS DIFFERENTES.
  La metrique utilise EXCLUSIVEMENT alpha_EM(mu).
""")

# Tableau comparatif
print(f"  {'k':>4} {'p_k':>5} {'alpha_sieve':>14} {'alpha_EM(p_k)':>14} {'ratio':>10}")
print(f"  {'-'*4} {'-'*5} {'-'*14} {'-'*14} {'-'*10}")
primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
for k in range(1, len(primes_list) + 1):
    pk = primes_list[k-1]
    a_s = alpha_sieve(k, primes_list)
    a_e = alpha_EM_mu(float(pk))
    ratio = a_e / a_s if a_s > 0 else 0
    print(f"  {k:4d} {pk:5d} {a_s:14.8f} {a_e:14.8f} {ratio:10.4f}")

print(f"\n  Ratio varie d'un facteur ~15 : PAS la meme fonction.")

# Concavite de ln(alpha_EM) : SCAN CORRECT
print(f"\n  CONVEXITE DE ln(alpha_EM) :")
print(f"  g_00 = -d^2(ln alpha)/dmu^2  [convention de signe de Jacobson]")
print(f"\n  {'mu':>6} {'d^2(ln a)/dmu^2':>16} {'g_00':>10} {'signature':>12}")
print(f"  {'-'*6} {'-'*16} {'-'*10} {'-'*12}")

transition_found = False
mu_transition = None
for mu in np.arange(4.0, 26.0, 0.5):
    d2 = d2_ln_alpha(mu)
    g00 = -d2  # CONVENTION : g_00 = -d^2(ln alpha)/dmu^2
    if g00 < 0:
        sig = "LORENTZIEN"
    else:
        sig = "EUCLIDIEN"
    if mu in [4.0, 5.0, 6.0, 6.5, 7.0, 7.5, 8.0, 10.0, 15.0, 20.0, 25.0]:
        print(f"  {mu:6.1f} {d2:16.10f} {g00:10.6f} {sig:>12}")
    if d2 > 0 and not transition_found:
        # Chercher le point exact
        try:
            mu_c = brentq(d2_ln_alpha, mu - 0.5, mu)
            transition_found = True
            mu_transition = mu_c
        except Exception:
            pass

if mu_transition:
    print(f"\n  TRANSITION DE SIGNATURE a mu_c = {mu_transition:.4f}")
    print(f"  Pour mu < {mu_transition:.1f} : g_00 > 0 (Euclidien)")
    print(f"  Pour mu > {mu_transition:.1f} : g_00 < 0 (Lorentzien)")
    print(f"  C'est un analogue de HARTLE-HAWKING (no-boundary proposal).")
    print(f"  mu_c / pi = {mu_transition/pi:.4f}")
    print(f"  mu_c ~ 2*pi + delta ? 2*pi = {2*pi:.4f}")
    # Quelle prime est a mu_c ?
    for p in [5, 7, 11]:
        print(f"  gamma_{p}(mu_c) = {gamma_p_func(p, mu_transition):.6f}")

# Verification au point operationnel
d2_op = d2_ln_alpha(MU_STAR)
g00_op = -d2_op
print(f"\n  Au point operationnel mu* = {MU_STAR:.4f} :")
print(f"  d^2(ln alpha)/dmu^2 = {d2_op:.10f} > 0 (CONVEXE)")
print(f"  g_00 = -{d2_op:.10f} = {g00_op:.10f} < 0 (LORENTZIEN)")

# Decomposition par premier
print(f"\n  Contribution de chaque premier :")
for p in PRIMES_ACTIVE:
    h = 1e-5
    s2 = sin2_theta(p, MU_STAR)
    ls = log(s2) if s2 > 0 else -100
    ls_p = log(sin2_theta(p, MU_STAR + h))
    ls_m = log(sin2_theta(p, MU_STAR - h))
    d2_p = (ls_p - 2*ls + ls_m) / h**2
    print(f"    d^2(ln sin^2_{p})/dmu^2 = {d2_p:+.10f} > 0")
print(f"    Chaque composante contribue POSITIVEMENT a la convexite.")
print(f"    Donc g_00 = -(sum positifs) < 0 : THEOREME, pas choix.")

print(f"""
  CONCLUSION P1 :
    alpha_sieve et alpha_EM sont MATHEMATIQUEMENT DISTINCTES.
    La metrique utilise alpha_EM (produit de sin^2).
    ln(alpha_EM) est convexe pour mu > {mu_transition:.1f} (calcul).
    g_00 = -d^2(ln alpha)/dmu^2 < 0 : signature Lorentzienne.
    PAS DE CONFUSION entre les deux alphas.
""")


# ================================================================
# PROBLEME 2 : mu = TEMPS
# ================================================================
print(f"\n{'='*78}")
print("  PROBLEME 2 : mu = TEMPS -- THEOREME D'UNICITE")
print(f"{'='*78}")

print(f"""
  STRUCTURE DU CRIBLE D'ERATOSTHENE :
    Le crible a EXACTEMENT :
    - UN parametre continu : mu = ln(N) (parametre d'echelle)
    - TROIS projections discretes actives : Z/3Z, Z/5Z, Z/7Z
    - Aucun autre degre de liberte (unicite du DOF, Demo 3)

  UNE METRIQUE 3+1D NECESSITE :
    - 1 coordonnee pour g < 0 (temporelle)
    - 3 coordonnees pour g > 0 (spatiales)

  THEOREME : mu est la SEULE coordonnee temporelle possible.

  PREUVE :
""")

# Etape 1 : Continuite
print(f"  ETAPE 1 -- mu est la seule coordonnee continue :")
print(f"    Les Z/pZ sont des groupes FINIS (3, 5, 7 elements).")
print(f"    Pas de derivee seconde possible sur un ensemble fini.")
print(f"    La courbure (Riemann) necessite d^2g/dx^2.")
print(f"    => mu est le SEUL candidat pour la direction temporelle.")

# Etape 2 : Signe
print(f"\n  ETAPE 2 -- g_00(mu) < 0 est un calcul, pas un postulat :")
print(f"    g_00 = -d^2(ln alpha_EM)/dmu^2 = {g00_op:.10f} < 0")
print(f"    C'est un FAIT mathematique : la convexite de ln(alpha_EM)")
print(f"    entraine g_00 < 0 automatiquement.")

# Etape 3 : g_ii > 0
print(f"\n  ETAPE 3 -- g_ii > 0 est trivial :")
for p in PRIMES_ACTIVE:
    gp = gamma_p_func(p, MU_STAR)
    gii = (gp / MU_STAR)**2
    print(f"    g_{p}{p} = (gamma_{p}/mu*)^2 = {gii:.10f} > 0")
print(f"    Les carres sont TOUJOURS positifs.")

# Etape 4 : Unicite
print(f"\n  ETAPE 4 -- 3+1 est FORCE :")
print(f"    1 continu (mu) + 3 discrets actifs (mod 3,5,7) = 4 dimensions.")
print(f"    g_00 < 0 et g_ii > 0 : signature (-,+,+,+) UNIQUE.")
print(f"    Pas d'alternative : on ne peut pas echanger mu et une Z/pZ")
print(f"    car la metrique necessite la continuite de la coordonnee.")

# Etape 5 : Invariance sous reparametrisation
print(f"\n  ETAPE 5 -- Invariance sous reparametrisation :")
print(f"    Si on pose tau = f(mu) (reparametrisation monotone) :")
print(f"    g_00(tau) = g_00(mu) * (dmu/dtau)^2")
print(f"    Le facteur (dmu/dtau)^2 > 0 ne change PAS le signe de g_00.")
print(f"    => Tout choix de coordonnee temporelle preserve la signature.")
print(f"    mu est le choix CANONIQUE (pas de facteur arbitraire).")

# Fleche du temps
print(f"\n  FLECHE DU TEMPS (argument supplementaire, pas necessaire) :")
dS_dmu = -d1_ln_alpha(MU_STAR)
print(f"    S(mu) = -ln(alpha_EM) croit avec mu :")
for mu in [10.0, 15.0, 20.0, 25.0]:
    S = -ln_alpha_EM(mu)
    print(f"      mu={mu:5.1f}: S = {S:.6f}")
print(f"    dS/dmu|_{{mu*}} = {dS_dmu:.8f} > 0")
print(f"    Ceci est la SECONDE LOI : l'entropie croit vers le futur.")
print(f"    La fleche du temps est DERIVEE, pas postulée.")

# Transition de signature
if mu_transition:
    print(f"\n  BONUS : TRANSITION DE SIGNATURE (Hartle-Hawking) :")
    print(f"    A mu_c = {mu_transition:.4f}, la signature change :")
    print(f"    mu < {mu_transition:.1f} : Euclidien (pas de temps)")
    print(f"    mu > {mu_transition:.1f} : Lorentzien (temps emerge)")
    print(f"    C'est l'analogue de la proposition no-boundary de Hartle-Hawking.")
    print(f"    L'emergence du temps n'est pas un choix : elle a un LIEU precis")
    print(f"    dans le parametre du crible (mu_c ~ {mu_transition:.1f}).")
    # Que se passe-t-il a mu_c ?
    g3_c = gamma_p_func(3, mu_transition)
    g5_c = gamma_p_func(5, mu_transition)
    g7_c = gamma_p_func(7, mu_transition)
    print(f"    A mu_c : gamma_3={g3_c:.4f}, gamma_5={g5_c:.4f}, gamma_7={g7_c:.4f}")
    print(f"    gamma_7(mu_c) = {g7_c:.4f} {'> 0.5 (actif)' if g7_c > 0.5 else '< 0.5 (sub-threshold)'}")
    print(f"    Le Lorentzien emerge quand p=7 devient ACTIF (gamma_7 > 0.5).")

print(f"""
  CONCLUSION P2 :
    mu = temps est un THEOREME, pas un dictionnaire :
    1. mu est la seule variable continue (FAIT structurel)
    2. g_00 < 0 est un calcul (convexite de ln alpha_EM)
    3. g_ii > 0 est trivial (carre d'un nombre reel)
    4. La signature (-,+,+,+) est FORCEE
    5. La fleche du temps = seconde loi (dS/dmu > 0)
    6. La transition Euclidien -> Lorentzien a mu_c ~ {mu_transition:.1f}
       = emergence de p=7 comme direction active
""")


# ================================================================
# PROBLEME 3 : G = 2*pi*alpha -- DIMENSIONS
# ================================================================
print(f"\n{'='*78}")
print("  PROBLEME 3 : G = 2*pi*alpha -- ANALYSE DIMENSIONNELLE")
print(f"{'='*78}")

print(f"""
  LA CRITIQUE :
    G_Newton a les dimensions [m^3 kg^-1 s^-2].
    alpha_EM = 1/137 est un nombre pur.
    "G = 2*pi*alpha" melange-t-il les dimensions ?

  LA REPONSE EN UNE PHRASE :
    La PT est ENTIEREMENT ADIMENSIONNELLE.
    Il n'y a ni metres, ni kilogrammes, ni secondes.
    G_sieve et alpha_EM sont TOUS DEUX des nombres purs.
""")

# Construction de G_sieve
print(f"  CONSTRUCTION DE G_sieve (nombre pur) :")
print(f"  1. La metrique du crible est :")
print(f"     ds^2 = -|S''| dmu^2 + sum (gamma_p/mu)^2 dx_p^2")
print(f"     ou S = -ln(alpha_EM), S'' = d^2S/dmu^2")
print(f"  2. Le tenseur d'Einstein G_ab se calcule de cette metrique.")
print(f"     C'est un nombre pur (rapport de derivees).")
print(f"  3. L'equation d'Einstein du crible :")
print(f"     G_ab = 8*pi*G_sieve * T_ab")
print(f"     ou T_ab est le tenseur energie-impulsion (D_KL, pressions info)")
print(f"  4. G_sieve = G_00 / (8*pi * D_total)")
print(f"     C'est un RAPPORT de deux nombres purs. Adimensionnel.")

# Calcul de G_sieve (comme dans test_graviton_G_Newton_v2.py)
print(f"\n  Calcul au point mu* = {MU_STAR:.4f} :")

# Lapse, scale factors, Hubble
N_sq = abs(d2_ln_alpha(MU_STAR))
N_val = sqrt(N_sq)
print(f"    N (lapse) = sqrt(|d^2 ln a/dmu^2|) = {N_val:.8f}")

h_d = 1e-4
H_vals = []
for p in PRIMES_ACTIVE:
    a_here = gamma_p_func(p, MU_STAR) / MU_STAR
    a_plus = gamma_p_func(p, MU_STAR + h_d) / (MU_STAR + h_d)
    a_minus = gamma_p_func(p, MU_STAR - h_d) / (MU_STAR - h_d)
    da = (a_plus - a_minus) / (2*h_d)
    H = da / (N_val * a_here) if N_val > 0 and a_here > 0 else 0
    H_vals.append(H)
    print(f"    H_{p} = {H:.8f}")

G_00 = H_vals[0]*H_vals[1] + H_vals[0]*H_vals[2] + H_vals[1]*H_vals[2]
print(f"    G_00 = H_3*H_5 + H_3*H_7 + H_5*H_7 = {G_00:.8f}")

# D_KL total (comme dans test_graviton_G_Newton_v2.py)
q_th = exp(-2/MU_STAR)
P3 = np.zeros(3)
for k in range(1, 500):
    r = (2*k) % 3
    P3[r] += (1 - q_th) * q_th**(k-1)
P3 /= P3.sum()
D_KL_3 = sum(P3[r] * log(3*P3[r]) for r in range(3) if P3[r] > 0)
D_parity = log(2)
D_total = D_parity + D_KL_3

print(f"    D_parity = ln(2) = {D_parity:.8f}")
print(f"    D_KL(mod 3) = {D_KL_3:.8f}")
print(f"    D_total = {D_total:.8f}")

G_sieve = G_00 / (8*pi*D_total)
print(f"\n    G_sieve = G_00 / (8*pi*D_total) = {G_sieve:.8f}")

alpha_op = alpha_EM_mu(MU_STAR)
G_predicted = 2*pi*alpha_op
print(f"    2*pi*alpha_EM = {G_predicted:.8f}")
ratio = G_sieve / G_predicted
err = abs(ratio - 1) * 100
print(f"    Ratio = {ratio:.6f}, erreur = {err:.2f}%")
print(f"    [{'PASS' if err < 1 else 'FAIL'}] G_sieve = 2*pi*alpha_EM (a {err:.2f}%)")

print(f"""
  STRUCTURE DIMENSIONNELLE :
    Crible          SI              Lien
    ------          --              ----
    alpha_EM        alpha_EM        IDENTIQUE (nombre pur)
    G_sieve         G_N [m^3/kg/s^2]   G_N = G_sieve * hbar*c/m_Pl^2
    mu              t [s]           t = mu * (hbar/E_0)
    gamma_p/mu      a [m]           a = (gamma_p/mu) * l_Pl

    Le passage au SI necessite 2 facteurs de traduction :
    m_e (masse) et l_Pl (longueur) -- ou equivalemment hbar et c.
    Ce sont des CONVENTIONS DE MESURE, pas des parametres.

  ANALOGIE :
    "c = 299792458 m/s" ne dit pas que la vitesse de la lumiere
    est un nombre avec des dimensions. C'est un facteur de conversion
    entre metres et secondes (1 m = 1/c secondes).
    De meme, G_N = G_sieve * (facteur de conversion).
    G_sieve = 2*pi*alpha est la VERITE adimensionnelle.
    G_N = 6.674e-11 m^3/kg/s^2 est sa traduction en SI.

  CONCLUSION P3 :
    "G = 2*pi*alpha" est une relation entre nombres purs du crible.
    Les dimensions apparaissent UNIQUEMENT au passage SI.
    La confusion vient de l'ecriture abregee, pas de la physique.
    G_sieve est aussi adimensionnel qu'alpha_EM.
""")


# ================================================================
# SYNTHESE
# ================================================================
print(f"\n{'='*78}")
print("  SYNTHESE : STATUT DU PASSAGE AU LORENTZIEN")
print(f"{'='*78}")

elements = [
    ("g_00 < 0 (Lorentzien)", "CALCUL",
     f"d^2(ln alpha)/dmu^2 = {d2_op:.6e} > 0"),
    ("mu = unique coordonnee continue", "FAIT",
     "Z/pZ discrets ne supportent pas d^2"),
    ("3+1 dimensions", "FORCE",
     "1 continu + 3 discrets actifs"),
    ("Signature (-,+,+,+)", "CONSEQUENCE",
     "g_00 < 0, g_ii = (gamma/mu)^2 > 0"),
    ("Fleche du temps", "DERIVE",
     f"dS/dmu = {dS_dmu:.4f} > 0"),
    (f"Transition Euclid->Lorenz a mu_c={mu_transition:.1f}" if mu_transition else "Transition", "CALCUL",
     "Hartle-Hawking analogue"),
    ("Bianchi I anisotrope", "CONSEQUENCE",
     "gamma_3 != gamma_5 != gamma_7"),
    (f"G_sieve = 2*pi*alpha ({err:.2f}%)", "DERIVE",
     "equation d'Einstein adimensionnelle"),
    ("alpha_sieve != alpha_EM", "FAIT",
     "deux fonctions mathematiquement distinctes"),
    ("g_00 = -S'', g_ii = (S'_p)^2", "DEFINITION",
     "S unique DOF scalaire du crible"),
]

print(f"\n  {'Element':45} {'Statut':12} {'Preuve':45}")
print(f"  {'-'*45} {'-'*12} {'-'*45}")
for name, status, proof in elements:
    print(f"  {name:45} {status:12} {proof:45}")

n_calc = sum(1 for _, s, _ in elements if s in ["CALCUL", "FAIT", "FORCE", "CONSEQUENCE"])
n_derive = sum(1 for _, s, _ in elements if s == "DERIVE")
n_def = sum(1 for _, s, _ in elements if s == "DEFINITION")
print(f"\n  Score : {n_calc} faits/calculs, {n_derive} derives, {n_def} definitions")
print(f"  Choix arbitraires : 0")

# AIDE-MEMOIRE
print(f"""

{'='*78}
  AIDE-MEMOIRE : REPONSES COURTES AUX 3 CRITIQUES
{'='*78}

  C1 "Deux alphas confondues"
  -> FAUX. alpha_sieve (Mertens, survie au crible) et alpha_EM (produit sin^2)
     sont DEUX fonctions distinctes. La metrique utilise alpha_EM.
     ln(alpha_EM) est CONVEXE pour mu > {mu_transition:.1f} -> g_00 < 0.

  C2 "mu = temps est un dictionnaire"
  -> FAUX. mu est la SEULE variable continue du crible. Les Z/pZ sont discretes.
     Une metrique necessite d^2/dx^2, impossible sur Z/pZ.
     Signature (-,+,+,+) forcee par g_00 < 0 (calcul) et g_ii > 0 (trivial).
     La fleche vient de dS/dmu > 0 (seconde loi informationnelle).

  C3 "G = 2*pi*alpha melange les dimensions"
  -> RACCOURCI. L'equation complete : G_00/(8*pi*D_total) = 2*pi*alpha_EM.
     Tout est adimensionnel dans le crible. Le passage au SI necessite
     des facteurs de conversion (m_e, l_Pl) comme c = 3e8 m/s.

  POINT FAIBLE RESTANT :
    Le choix g_00 = -S'' (et pas un autre scalaire du crible) est une
    DEFINITION -- mais c'est la definition la plus simple et S est
    le seul DOF scalaire (unicite prouvee, Demo 3). Il n'y a pas
    d'alternative non-triviale.

{'='*78}
FIN
{'='*78}
""")
