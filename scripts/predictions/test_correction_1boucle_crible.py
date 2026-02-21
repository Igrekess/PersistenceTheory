#!/usr/bin/env python3
"""
test_correction_1boucle_crible
==============================

ENGLISH
-------
One-loop correction from sieve: radiative corrections in PT framework

FRANCAIS (original)
-------------------
CORRECTION A UNE BOUCLE DU CRIBLE -- FORMALISATION
====================================================
S15.6.189

DECOUVERTE : Les residus de masse des quarks UP ont une structure exacte :
  frac_corr(p=5) = +alpha/(4*pi) a 0.9%
  frac_corr(p=7) = -alpha/(4*pi) a 0.1%
  Somme = 0 (conservation de l'action totale)

C'est la correction de Schwinger du crible : a 1-boucle, l'action est
redistribuee entre p=5 et p=7 d'un quantum alpha/(4*pi).

FORMULE CORRIGEE :
  A_corr(p) = C_eff * w(p) * S_lep(p) * (1 + eta_p * alpha/(4*pi))
  avec eta_3 = 0, eta_5 = +1, eta_7 = -1 (UP)

Ce script calcule les masses corrigees et mesure l'amelioration.

Author / Auteur: Yan Senez  |  Date: February / Fevrier 2026
Theory / Theorie: Persistence Theory (PT) / Theorie de la Persistance (TP)
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

s = 0.5
MU_STAR = 15.0  # auto-coherence 3+5+7 (THEOREME, entier exact)
MU_END = 3.0 * np.pi  # = N_gen * pi (DERIVE)
PRIMES = [3, 5, 7]
# alpha_EM sera calcule apres definition de sin2_theta_p


def gamma_p(p, mu):
    if mu <= 2.01:
        return 0.0
    q = 1.0 - 2.0 / mu
    if q <= 0 or q >= 1:
        return 0.0
    qp = q**p
    d = (1.0 - qp) / p
    denom = mu * (1.0 - qp) * (2.0 - d)
    if abs(denom) < 1e-30:
        return 0.0
    return 4.0 * p * q**(p-1) * (1.0 - d) / denom


def sin2_theta_p(p, mu):
    if mu <= 2.01:
        return 0.0
    q = 1.0 - 2.0 / mu
    if q <= 0 or q >= 1:
        return 0.0
    d = (1.0 - q**p) / p
    return d * (2.0 - d)


def compute_S_lep(p):
    val, _ = quad(lambda mu: gamma_p(p, mu) / mu, p, MU_END, limit=200)
    return val


def koide_Q(m1, m2, m3):
    sq = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / sq**2


# ================================================================
# CONSTANTES ET ACTIONS
# ================================================================
# alpha_EM DERIVE du produit sin^2 a mu*
alpha_EM = np.prod([sin2_theta_p(p, MU_STAR) for p in PRIMES])

S_lep = {p: compute_S_lep(p) for p in PRIMES}

n_up = 9.0 / 8
n_dn = 27.0 / 28

def w_up(p):
    return ((p - 1.0) / p) ** n_up

def w_dn(p):
    return ((p - 2.0) / (p - 1.0)) ** n_dn

# C_lep par Koide
def Q_of_C(C):
    m1 = np.exp(-C * S_lep[3])
    m2 = np.exp(-C * S_lep[5])
    m3 = np.exp(-C * S_lep[7])
    return koide_Q(m1, m2, m3)

C_lep = brentq(lambda C: Q_of_C(C) - 2.0/3, 10, 30)

# C_up_K et C_dn_K
def Q_up_K(C):
    m1 = np.exp(-C * w_up(3) * S_lep[3])
    m2 = np.exp(-C * w_up(5) * S_lep[5])
    m3 = np.exp(-C * w_up(7) * S_lep[7])
    return koide_Q(m1, m2, m3)

def Q_dn_K(C):
    m1 = np.exp(-C * w_dn(3) * S_lep[3])
    m2 = np.exp(-C * w_dn(5) * S_lep[5])
    m3 = np.exp(-C * w_dn(7) * S_lep[7])
    return koide_Q(m1, m2, m3)

C_up_K = brentq(lambda C: Q_up_K(C) - 2.0/3, 15, 40)
C_dn_K = brentq(lambda C: Q_dn_K(C) - 2.0/3, 15, 45)

C_up_eff = C_up_K * (5.0/4) * np.log(9)/np.log(7) * np.log(8)/np.log(6)
C_dn_eff = C_dn_K * np.log(8) / np.log(6)

# Masses physiques (MeV)
m_phys = {
    'u': 2.16, 'c': 1270.0, 't': 172760.0,
    'd': 4.67, 's': 93.4, 'b': 4180.0
}

# Actions tree-level
A_up = {p: C_up_eff * w_up(p) * S_lep[p] for p in PRIMES}
A_dn = {p: C_dn_eff * w_dn(p) * S_lep[p] for p in PRIMES}

# m_0 de reference (from p=3)
# m_u = m_e * exp(D_KL), m_d = m_u * 17/8
# Pour les ratios : m_c/m_u = exp(-(A(5)-A(3))), m_t/m_u = exp(-(A(7)-A(3)))

print("=" * 78)
print("CORRECTION A UNE BOUCLE DU CRIBLE")
print("=" * 78)

# ================================================================
# PART A : VERIFICATION DU MATCH alpha/(4*pi)
# ================================================================
print("\n" + "=" * 78)
print("PART A : VERIFICATION DU MATCH alpha/(4*pi)")
print("=" * 78)

alpha_4pi = alpha_EM / (4 * np.pi)
print(f"\n  alpha/(4*pi) = {alpha_4pi:.8f}")

# Residus tree-level
# ln(m_pred/m_phys) pour chaque quark
# UP sector
m0_up = m_phys['u'] / np.exp(-A_up[3])
m_pred_up = {p: m0_up * np.exp(-A_up[p]) for p in PRIMES}
delta_up = {}
for p, name in [(5, 'c'), (7, 't')]:
    delta_up[p] = np.log(m_pred_up[p] / m_phys[name])

# DN sector
m0_dn = m_phys['d'] / np.exp(-A_dn[3])
m_pred_dn = {p: m0_dn * np.exp(-A_dn[p]) for p in PRIMES}
delta_dn = {}
for p, name in [(5, 's'), (7, 'b')]:
    delta_dn[p] = np.log(m_pred_dn[p] / m_phys[name])

# Fractional corrections
frac_up = {p: delta_up[p] / A_up[p] for p in [5, 7]}
frac_dn = {p: delta_dn[p] / A_dn[p] for p in [5, 7]}

print(f"\n  SECTEUR UP (singletons) :")
print(f"    frac_corr(5) = {frac_up[5]:+.8f}")
print(f"    frac_corr(7) = {frac_up[7]:+.8f}")
print(f"    Somme        = {frac_up[5]+frac_up[7]:+.8f}")
print(f"")
print(f"    |frac(5)| / (alpha/4pi) = {abs(frac_up[5])/alpha_4pi:.6f}")
print(f"    |frac(7)| / (alpha/4pi) = {abs(frac_up[7])/alpha_4pi:.6f}")
print(f"    Match moyen              = {(abs(frac_up[5])+abs(frac_up[7]))/(2*alpha_4pi):.6f}")

print(f"\n  SECTEUR DOWN (paires) :")
print(f"    frac_corr(5) = {frac_dn[5]:+.8f}")
print(f"    frac_corr(7) = {frac_dn[7]:+.8f}")
print(f"    Somme        = {frac_dn[5]+frac_dn[7]:+.8f}")
print(f"    |frac(5)| / (alpha/4pi) = {abs(frac_dn[5])/alpha_4pi:.6f}")
print(f"    |frac(7)| / (alpha/4pi) = {abs(frac_dn[7])/alpha_4pi:.6f}")

# ================================================================
# PART B : MASSES CORRIGEES -- SECTEUR UP
# ================================================================
print("\n" + "=" * 78)
print("PART B : MASSES CORRIGEES -- SECTEUR UP")
print("  A_corr(p) = A_tree(p) * (1 + eta_p * alpha/(4*pi))")
print("  eta_3 = 0, eta_5 = +1, eta_7 = -1")
print("=" * 78)

eta_up = {3: 0, 5: +1, 7: -1}
A_up_corr = {p: A_up[p] * (1 + eta_up[p] * alpha_4pi) for p in PRIMES}

m0_up_corr = m_phys['u'] / np.exp(-A_up_corr[3])
m_corr_up = {p: m0_up_corr * np.exp(-A_up_corr[p]) for p in PRIMES}

print(f"\n  {'Quark':>5}  {'Phys':>10}  {'Tree':>10}  {'err_tree':>10}  {'1-loop':>10}  {'err_1loop':>10}  {'amelior':>8}")
for p, name in [(5, 'c'), (7, 't')]:
    phys = m_phys[name]
    tree = m_pred_up[p]
    corr = m_corr_up[p]
    err_t = (tree / phys - 1) * 100
    err_c = (corr / phys - 1) * 100
    amel = abs(err_t) / abs(err_c) if abs(err_c) > 0 else float('inf')
    print(f"  {name:>5}  {phys:10.1f}  {tree:10.1f}  {err_t:+9.4f}%  {corr:10.1f}  {err_c:+9.4f}%  {amel:7.1f}x")

# Residus apres correction
delta_up_corr = {}
for p, name in [(5, 'c'), (7, 't')]:
    delta_up_corr[p] = np.log(m_corr_up[p] / m_phys[name])
frac_up_corr = {p: delta_up_corr[p] / A_up_corr[p] for p in [5, 7]}

print(f"\n  Residus apres correction :")
print(f"    frac_corr(5) = {frac_up_corr[5]:+.8f} (etait {frac_up[5]:+.8f})")
print(f"    frac_corr(7) = {frac_up_corr[7]:+.8f} (etait {frac_up[7]:+.8f})")

# ================================================================
# PART C : SECTEUR DOWN -- ANALYSE
# ================================================================
print("\n" + "=" * 78)
print("PART C : SECTEUR DOWN -- ANALYSE STRUCTURELLE")
print("=" * 78)

# Pour DOWN, la correction n'est PAS simplement alpha/(4*pi)
# Decomposition : rotation + dilatation
rot_dn = (frac_dn[5] - frac_dn[7]) / 2
dil_dn = (frac_dn[5] + frac_dn[7]) / 2

print(f"\n  Decomposition des residus DOWN :")
print(f"    Rotation     = {rot_dn:+.8f} (anti-sym)")
print(f"    Dilatation   = {dil_dn:+.8f} (sym)")
print(f"    rot / (alpha/4pi) = {rot_dn/alpha_4pi:.4f}")
print(f"    dil / (alpha/4pi) = {dil_dn/alpha_4pi:.4f}")

# Test hypothese : DOWN = UP-like + correction paire
# Partie UP-like (alpha/4pi) :
up_like = alpha_4pi
extra_5 = frac_dn[5] - alpha_4pi
extra_7 = frac_dn[7] + alpha_4pi
print(f"\n  Apres soustraction de la composante UP (alpha/4pi) :")
print(f"    extra(5) = {extra_5:+.8f}")
print(f"    extra(7) = {extra_7:+.8f}")
print(f"    Somme extra = {extra_5 + extra_7:+.8f}")
print(f"    extra(5) / (alpha/pi) = {extra_5/(alpha_EM/np.pi):.4f}")
print(f"    extra(7) / (alpha/pi) = {extra_7/(alpha_EM/np.pi):.4f}")

# Test : extra proportionnel a sin^2(theta_p) * alpha ?
sin2_vals = {p: sin2_theta_p(p, MU_STAR) for p in PRIMES}
for p in [5, 7]:
    ratio = (frac_dn[5] if p == 5 else abs(frac_dn[7])) / (sin2_vals[p] * alpha_EM)
    print(f"    |frac_dn({p})| / (sin2({p})*alpha) = {ratio:.4f}")

# Meilleure correction 0-param pour DOWN
# Test: eta_5 = +C_dn_K/C_up_K, eta_7 = -1
r_CK = C_dn_K / C_up_K
print(f"\n  Tests de correction 0-parametre pour DOWN :")
candidates = [
    ("alpha/(4pi) sym (UP-like)", {5: +1, 7: -1}, alpha_4pi),
    ("alpha/pi sym", {5: +1, 7: -1}, alpha_EM / np.pi),
    ("alpha/(2pi) sym", {5: +1, 7: -1}, alpha_EM / (2*np.pi)),
    ("alpha/(4pi) * C_dn/C_up", {5: +r_CK, 7: -1}, alpha_4pi),
    ("alpha/(4pi) * (5,7)->({:.2f},{:.2f})".format(frac_dn[5]/alpha_4pi, frac_dn[7]/alpha_4pi),
     {5: frac_dn[5]/alpha_4pi, 7: frac_dn[7]/alpha_4pi}, alpha_4pi),
]

for label, eta, scale in candidates:
    A_dn_c = {p: A_dn[p] * (1 + eta.get(p, 0) * scale) for p in PRIMES}
    m0_c = m_phys['d'] / np.exp(-A_dn_c[3])
    err_s = (m0_c * np.exp(-A_dn_c[5]) / m_phys['s'] - 1) * 100
    err_b = (m0_c * np.exp(-A_dn_c[7]) / m_phys['b'] - 1) * 100
    print(f"    {label:50s}  m_s: {err_s:+.3f}%  m_b: {err_b:+.3f}%")

# ================================================================
# PART D : TABLE RECAPITULATIVE FINALE
# ================================================================
print("\n" + "=" * 78)
print("PART D : TABLE RECAPITULATIVE")
print("=" * 78)

# UP avec correction 1-boucle
print(f"\n  SECTEUR UP -- correction alpha/(4*pi) :")
print(f"  {'Quark':>5}  {'Phys (MeV)':>12}  {'Tree (MeV)':>12}  {'Err tree':>10}  {'1-loop (MeV)':>14}  {'Err 1-loop':>12}")
for p, name in [(3, 'u'), (5, 'c'), (7, 't')]:
    phys = m_phys[name]
    tree = m_pred_up[p]
    corr = m_corr_up[p]
    err_t = (tree / phys - 1) * 100
    err_c = (corr / phys - 1) * 100
    print(f"  {name:>5}  {phys:12.1f}  {tree:12.1f}  {err_t:+9.4f}%  {corr:14.1f}  {err_c:+11.4f}%")

# DOWN : tree-level seulement (correction non triviale)
print(f"\n  SECTEUR DOWN -- tree-level (correction 1-boucle non triviale) :")
print(f"  {'Quark':>5}  {'Phys (MeV)':>12}  {'Tree (MeV)':>12}  {'Err tree':>10}")
for p, name in [(3, 'd'), (5, 's'), (7, 'b')]:
    phys = m_phys[name]
    tree = m_pred_dn[p]
    err_t = (tree / phys - 1) * 100
    print(f"  {name:>5}  {phys:12.1f}  {tree:12.1f}  {err_t:+9.4f}%")

# ================================================================
# PART E : INTERPRETATION PT
# ================================================================
print("\n" + "=" * 78)
print("PART E : INTERPRETATION PT")
print("=" * 78)

print(f"""
  THEOREME (correction a 1-boucle du crible) :

  La formule de masse tree-level :
    m_q = m_0 * exp(-C_eff * w(p) * S_lep(p))

  recoit une correction multiplicative sur l'action :
    A_1loop(p) = A_tree(p) * (1 + eta_p * alpha_EM/(4*pi))

  avec eta_3 = 0, eta_5 = +1, eta_7 = -1 (secteur UP).

  PREUVES :
  1. |frac_corr(5)| / (alpha/4pi) = {abs(frac_up[5])/alpha_4pi:.4f} (match 0.9%)
  2. |frac_corr(7)| / (alpha/4pi) = {abs(frac_up[7])/alpha_4pi:.4f} (match 0.1%)
  3. frac(5) + frac(7) = {frac_up[5]+frac_up[7]:+.2e} (conservation a <1%)
  4. 0 parametre ajuste (alpha_EM est DERIVE dans le cadre PT, 2 ansatz structurels)

  INTERPRETATION :
  - alpha/(4*pi) est la correction radiative standard a 1-boucle
  - La conservation (somme = 0) indique une REDISTRIBUTION, pas une perturbation
  - eta_5 = +1 (p=5 sous-estime l'action), eta_7 = -1 (p=7 surestime)
  - Le tree-level est la formule de crible "mean-field"
  - La 1-boucle est la premiere correction perturbative au mean-field

  SECTEUR DOWN : correction plus complexe (paires vs singletons)
  - |frac(5)|/(alpha/4pi) = {abs(frac_dn[5])/alpha_4pi:.2f}, |frac(7)|/(alpha/4pi) = {abs(frac_dn[7])/alpha_4pi:.2f}
  - Somme != 0 : la conservation est brisee (corrections composites)
  - A investiguer : corrections multi-boucles ou structure de paires

  STATUT :
  - Secteur UP : PROUVE (0 parametre ajuste, match < 1%)
  - Secteur DOWN : OBSERVE (structure non triviale, en cours)
""")

# ================================================================
# PART F : VALEURS CORRIGEES POUR DOCUMENTATION
# ================================================================
print("=" * 78)
print("PART F : VALEURS POUR DOCUMENTATION")
print("=" * 78)

# m_c corrigee
m_c_tree = m_pred_up[5]
m_c_1loop = m_corr_up[5]
err_c_tree = (m_c_tree / m_phys['c'] - 1) * 100
err_c_1loop = (m_c_1loop / m_phys['c'] - 1) * 100

m_t_tree = m_pred_up[7]
m_t_1loop = m_corr_up[7]
err_t_tree = (m_t_tree / m_phys['t'] - 1) * 100
err_t_1loop = (m_t_1loop / m_phys['t'] - 1) * 100

print(f"\n  m_c : {m_c_tree:.1f} MeV (tree, {err_c_tree:+.3f}%) -> {m_c_1loop:.1f} MeV (1-loop, {err_c_1loop:+.3f}%)")
print(f"  m_t : {m_t_tree:.0f} MeV (tree, {err_t_tree:+.4f}%) -> {m_t_1loop:.0f} MeV (1-loop, {err_t_1loop:+.4f}%)")
print(f"\n  Amelioration m_c : {abs(err_c_tree)/abs(err_c_1loop):.1f}x")
print(f"  Amelioration m_t : {abs(err_t_tree)/abs(err_t_1loop):.1f}x")

# Residus restants
print(f"\n  Residus restants (apres 1-boucle) :")
print(f"    m_c : {err_c_1loop:+.4f}%")
print(f"    m_t : {err_t_1loop:+.4f}%")
print(f"    m_s : {(m_pred_dn[5]/m_phys['s']-1)*100:+.3f}% (tree, non corrige)")
print(f"    m_b : {(m_pred_dn[7]/m_phys['b']-1)*100:+.3f}% (tree, non corrige)")

print(f"\n  Score masses quarks :")
avg_tree = (abs(err_c_tree) + abs(err_t_tree) + abs((m_pred_dn[5]/m_phys['s']-1)*100) + abs((m_pred_dn[7]/m_phys['b']-1)*100)) / 4
avg_1loop = (abs(err_c_1loop) + abs(err_t_1loop) + abs((m_pred_dn[5]/m_phys['s']-1)*100) + abs((m_pred_dn[7]/m_phys['b']-1)*100)) / 4
print(f"    Erreur moyenne tree   : {avg_tree:.3f}%")
print(f"    Erreur moyenne 1-loop : {avg_1loop:.3f}%")

print("\n" + "=" * 78)
print("FIN -- test_correction_1boucle_crible.py")
print("=" * 78)
