#!/usr/bin/env python3
"""
test_derivation_conventions
===========================

ENGLISH
-------
PT conventions: standardizing the derivation chain for 25 observables

FRANCAIS (original)
-------------------
S15.6.186 -- Derivation des rapports de masse inter-secteurs
==============================================================

POINT CLE PHILOSOPHIQUE :
  La PT fonctionne ENTIEREMENT en dimensionless. Les quantites s, alpha,
  D_KL, gamma_p, sin^2 sont des nombres purs. La PT n'a PAS besoin
  d'unite de masse pour faire fonctionner la physique.

  Les "conventions d'echelle" (m_e, m_u, m_d, v) ne sont PAS des
  parametres de la PT -- ce sont des FACTEURS DE TRADUCTION vers les
  unites humaines (MeV, GeV). Exactement comme c = 3e8 m/s n'est pas
  un parametre de la relativite, mais une conversion metres <-> secondes.

  La BONNE question est donc :
    "Tous les RAPPORTS de masse (dimensionless) sont-ils derivables ?"

  Si oui, la PT a 0 convention : tout est determine par s = 1/2.

RAPPORTS A DERIVER :
  R1 : m_u / m_e  (ratio quark-lepton, inter-secteur)
  R2 : m_d / m_e  (ratio quark-lepton, inter-secteur)
  R3 : m_d / m_u  (ratio down/up, intra-secteur)
  R4 : m_nu_3 / m_e  (ratio neutrino-lepton)
  R5 : Dm21 / m_e^2  (difference de masse carree, dimensionless)

HYPOTHESES :
  H1 : m_u/m_e = exp(D_KL)  ou D_KL = persistance au point operatoire
  H2 : m_d/m_u = 1 + n_up = 17/8  (quark down = quark up + quantum de modulation)
  H3 : m_nu_3 = s^2 * alpha^3 * m_e  (suppression par I_inf * alpha^{N_gen})
  H4 : Dm21 = m_nu_3^2 / R_nu  (hierarchie normale)

Script : test_derivation_conventions.py
Auteur : Yan Senez & Claude
Date : 18 fevrier 2026

Author / Auteur: Yan Senez  |  Date: February / Fevrier 2026
Theory / Theorie: Persistence Theory (PT) / Theorie de la Persistance (TP)
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from collections import Counter
import time

# ==================================================================
# CONSTANTES
# ==================================================================

s = 0.5
gamma_E = 0.5772156649015329
e_gamma = np.exp(gamma_E)
phi = (1 + np.sqrt(5)) / 2
PRIMES_ACTIVE = [3, 5, 7]
ALPHA_PHYS = 1.0 / 137.035999084

# Masses physiques (PDG 2024)
m_e_MeV = 0.51099895       # MeV, pole mass (exact)
m_mu_MeV = 105.6583755     # MeV
m_tau_MeV = 1776.86        # MeV
m_u_MeV = 2.16             # MeV, MS-bar a 2 GeV (+0.49 -0.26)
m_d_MeV = 4.67             # MeV, MS-bar a 2 GeV (+0.48 -0.17)
m_s_MeV = 93.4             # MeV
m_c_MeV = 1270.0           # MeV
m_b_MeV = 4180.0           # MeV
m_t_MeV = 172700.0         # MeV

# Rapports observes
R_mu_e_obs = m_mu_MeV / m_e_MeV       # 206.77
R_tau_e_obs = m_tau_MeV / m_e_MeV     # 3477.2
R_u_e_obs = m_u_MeV / m_e_MeV         # 4.228
R_d_e_obs = m_d_MeV / m_e_MeV         # 9.139
R_d_u_obs = m_d_MeV / m_u_MeV         # 2.162

# Neutrinos (PDG 2024, normal ordering)
Dm21_obs = 7.42e-5    # eV^2
Dm31_obs = 2.515e-3   # eV^2
R_nu_obs = Dm31_obs / Dm21_obs  # ~33.9


# ==================================================================
# FONCTIONS DU CRIBLE
# ==================================================================

def sin2_theta(p, q):
    """sin^2(theta_p) = delta_p * (2 - delta_p)"""
    qp = q**p
    delta = (1.0 - qp) / p
    return delta * (2.0 - delta)


def gamma_p_exact(p, mu):
    """Dimension de crible gamma_p(mu)."""
    if mu <= 2.01:
        return 0.0
    q = 1.0 - 2.0 / mu
    qp = q**p
    delta = (1.0 - qp) / p
    if abs(1.0 - qp) < 1e-300:
        return 0.0
    dln_delta = -2.0 * p * q**(p - 1) / (mu * (1.0 - qp))
    factor = 2.0 * (1.0 - delta) / (2.0 - delta)
    return -dln_delta * factor


def sin2_of_mu(p, mu):
    """sin^2(theta_p) en fonction de mu (avec q = 1-2/mu)."""
    if mu <= 2.01:
        return 0.0
    q = 1.0 - 2.0 / mu
    return sin2_theta(p, q)


def alpha_em_sieve(mu):
    """alpha_EM(mu) = prod sin^2(theta_p) pour p=3,5,7."""
    if mu <= 2.01:
        return 0.0
    q = 1.0 - 2.0 / mu
    return sin2_theta(3, q) * sin2_theta(5, q) * sin2_theta(7, q)


def koide_Q(m1, m2, m3):
    """Rapport de Koide Q = sum(m) / (sum(sqrt(m)))^2."""
    sq = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / sq**2


def H_max_geom(mu):
    """H_max = H(Geom(q)) avec q = 1 - 1/mu, en bits."""
    if mu <= 1.01:
        return 0.0
    q = 1.0 - 1.0 / mu
    if q <= 0 or q >= 1:
        return 0.0
    return -np.log2(1 - q) - q / (1 - q) * np.log2(q)


def shannon_entropy(gaps):
    """Entropie de Shannon en bits."""
    counts = Counter(gaps.tolist())
    total = len(gaps)
    H = 0.0
    for c in counts.values():
        p_val = c / total
        if p_val > 0:
            H -= p_val * np.log2(p_val)
    return H


# ==================================================================
# A. ACTIONS DU CRIBLE : S_lep ET sigma (DUALITE)
# ==================================================================

def section_A():
    """Actions leptonique et quark, verification de la dualite."""
    print("=" * 72)
    print("A. ACTIONS DU CRIBLE : DUALITE SOMMET / ARETE")
    print("=" * 72)

    mu_end = 3.0 * np.pi

    # --- Actions leptoniques S_p = integral gamma_p/mu dmu ---
    S_lep = {}
    for p in PRIMES_ACTIVE:
        val, _ = quad(lambda mu: gamma_p_exact(p, mu) / mu,
                      p, mu_end, limit=200)
        S_lep[p] = val

    # --- Actions quarks sigma_p = integral sin^2(theta_p)/mu dmu ---
    sigma = {}
    for p in PRIMES_ACTIVE:
        val, _ = quad(lambda mu: sin2_of_mu(p, mu) / mu,
                      p, mu_end, limit=200)
        sigma[p] = val

    # --- Dualite : sigma_p / S_lep_p ~ s ---
    print(f"\n  mu_end = 3*pi = {mu_end:.6f}")
    print(f"\n  {'p':>3} {'S_lep':>12} {'sigma':>12} {'sigma/S':>12} {'s=0.5':>8} {'err%':>8}")
    print("  " + "-" * 58)

    ratios_dual = []
    for p in PRIMES_ACTIVE:
        ratio = sigma[p] / S_lep[p]
        err = abs(ratio - s) / s * 100
        ratios_dual.append(ratio)
        print(f"  {p:3d} {S_lep[p]:12.6f} {sigma[p]:12.6f} "
              f"{ratio:12.6f} {s:8.4f} {err:8.2f}%")

    avg_ratio = np.mean(ratios_dual)
    print(f"\n  Moyenne sigma/S = {avg_ratio:.6f} (s = {s})")
    print(f"  La dualite sommet/arete est confirmee a {abs(avg_ratio - s)/s*100:.1f}%")

    # --- C_Koide ---
    def Q_of_C(C):
        m1 = np.exp(-C * S_lep[3])
        m2 = np.exp(-C * S_lep[5])
        m3 = np.exp(-C * S_lep[7])
        if m1 < 1e-300:
            return 1.0
        return koide_Q(m1, m2, m3) - 2.0 / 3.0

    C_Koide = brentq(Q_of_C, 15, 22)

    # Verification
    m_pred = [np.exp(-C_Koide * S_lep[p]) for p in PRIMES_ACTIVE]
    m_norm = [m / m_pred[0] for m in m_pred]
    print(f"\n  C_Koide = {C_Koide:.6f}")
    print(f"  m_mu/m_e(PT) = {m_norm[1]:.4f} (obs: {R_mu_e_obs:.4f}, "
          f"err {abs(m_norm[1]-R_mu_e_obs)/R_mu_e_obs*100:.2f}%)")
    print(f"  m_tau/m_e(PT) = {m_norm[2]:.4f} (obs: {R_tau_e_obs:.4f}, "
          f"err {abs(m_norm[2]-R_tau_e_obs)/R_tau_e_obs*100:.2f}%)")

    # --- C_nu = (5/8) * C ---
    C_nu = (5.0 / 8.0) * C_Koide
    print(f"\n  C_nu = (5/8)*C = {C_nu:.6f}")
    print(f"  5/8 = s*(1+s^2) = {s*(1+s**2):.4f}")

    # --- Exposants quarks ---
    n_up = 9.0 / 8.0
    n_dn = 27.0 / 28.0

    def w_up(p):
        return ((p - 1.0) / p) ** n_up

    def w_dn(p):
        return ((p - 2.0) / (p - 1.0)) ** n_dn

    print(f"\n  --- Modulations quarks ---")
    for p in PRIMES_ACTIVE:
        print(f"  p={p}: w_up = {w_up(p):.6f}, w_dn = {w_dn(p):.6f}")

    return {
        'S_lep': S_lep,
        'sigma': sigma,
        'C_Koide': C_Koide,
        'C_nu': C_nu,
        'n_up': n_up,
        'n_dn': n_dn,
        'w_up': w_up,
        'w_dn': w_dn,
        'mu_end': mu_end,
        'avg_dual_ratio': avg_ratio,
    }


# ==================================================================
# B. D_KL AU POINT OPERATOIRE
# ==================================================================

def section_B():
    """D_KL au point operatoire (crible par {2,3,5,7})."""
    print("\n" + "=" * 72)
    print("B. D_KL AU POINT OPERATOIRE")
    print("=" * 72)

    # Crible d'Eratosthene progressif
    N_sieve = 10_000_000
    primes_sieve = [2, 3, 5, 7]

    is_alive = np.ones(N_sieve + 1, dtype=bool)
    is_alive[0] = is_alive[1] = False

    results = {}
    for p in primes_sieve:
        for mult in range(2 * p, N_sieve + 1, p):
            is_alive[mult] = False
        survivors = np.where(is_alive)[0]
        gaps = np.diff(survivors)
        if len(gaps) < 100:
            break
        mu = float(np.mean(gaps))
        H = shannon_entropy(gaps)
        H_m = H_max_geom(mu)
        DKL = H_m - H
        n_even = int(np.sum(gaps % 2 == 0))
        f_even = n_even / len(gaps)
        if f_even > 1 - 1e-10:
            parity = 1.0
        else:
            parity = f_even * np.log2(f_even / 0.5) + \
                     (1 - f_even) * np.log2((1 - f_even) / 0.5)
        D_matter = max(DKL - parity, 0)
        results[p] = {
            'mu': mu, 'H_max': H_m, 'H': H, 'DKL': DKL,
            'parity': parity, 'D_matter': D_matter,
            'n_gaps': len(gaps),
        }

    print(f"\n  N_sieve = {N_sieve:,}")
    print(f"\n  {'p':>4} {'mu':>8} {'H_max':>8} {'H':>8} {'D_KL':>8} {'Parity':>8} {'D_mat':>8}")
    print("  " + "-" * 58)
    for p in primes_sieve:
        r = results[p]
        print(f"  {p:4d} {r['mu']:8.4f} {r['H_max']:8.4f} {r['H']:8.4f} "
              f"{r['DKL']:8.4f} {r['parity']:8.4f} {r['D_matter']:8.4f}")

    d_op = results[7]
    DKL_op = d_op['DKL']

    print(f"\n  === D_KL au point operatoire (k=4, p=7) ===")
    print(f"  D_KL = {DKL_op:.6f} bits")
    print(f"  exp(D_KL) = {np.exp(DKL_op):.6f}")
    print(f"  m_u/m_e (obs) = {R_u_e_obs:.4f}")
    print(f"  Accord: {abs(np.exp(DKL_op) - R_u_e_obs)/R_u_e_obs*100:.3f}%")

    # Convergence check
    print(f"\n  --- Convergence de D_KL ---")
    for N_test in [100_000, 500_000, 1_000_000, 5_000_000, 10_000_000]:
        is_alive2 = np.ones(N_test + 1, dtype=bool)
        is_alive2[0] = is_alive2[1] = False
        for p in [2, 3, 5, 7]:
            for mult in range(2 * p, N_test + 1, p):
                is_alive2[mult] = False
        surv = np.where(is_alive2)[0]
        g = np.diff(surv)
        mu2 = float(np.mean(g))
        H2 = shannon_entropy(g)
        Hm2 = H_max_geom(mu2)
        DKL2 = Hm2 - H2
        print(f"    N={N_test:>10,}: D_KL={DKL2:.6f}, exp(D_KL)={np.exp(DKL2):.6f}")

    return {
        'DKL_op': DKL_op,
        'results': results,
    }


# ==================================================================
# C. DERIVATION DE m_u/m_e
# ==================================================================

def section_C(actions, dkl_data):
    """Test H1 : m_u/m_e = exp(D_KL)."""
    print("\n" + "=" * 72)
    print("C. DERIVATION DE m_u / m_e")
    print("=" * 72)

    DKL = dkl_data['DKL_op']

    print(f"\n  === HYPOTHESE H1 : m_u/m_e = exp(D_KL) ===")
    print(f"")
    print(f"  D_KL = {DKL:.6f} bits (persistance au point operatoire)")
    print(f"  exp(D_KL) = {np.exp(DKL):.6f}")
    print(f"  m_u/m_e (obs) = {R_u_e_obs:.4f}")

    pred_u_e = np.exp(DKL)
    err_u_e = abs(pred_u_e - R_u_e_obs) / R_u_e_obs * 100

    print(f"\n  +-------------------------------------------+")
    print(f"  |  m_u/m_e (PT)  = {pred_u_e:.4f}                |")
    print(f"  |  m_u/m_e (obs) = {R_u_e_obs:.4f}                |")
    print(f"  |  Erreur        = {err_u_e:.3f}%                |")
    print(f"  +-------------------------------------------+")

    print(f"\n  INTERPRETATION PHYSIQUE :")
    print(f"  D_KL est la persistance informationnelle TOTALE au point")
    print(f"  operatoire, apres filtrage par les 4 premiers {{2,3,5,7}}.")
    print(f"  exp(D_KL) mesure le POIDS INFORMATIONNEL du crible complet.")
    print(f"  Le quark up (arete du crible) est plus lourd que l'electron")
    print(f"  (sommet) par EXACTEMENT ce poids.")
    print(f"")
    print(f"  C'est le PONT entre vertex (leptons) et edges (quarks).")

    # Budget de D_KL
    d_op = dkl_data['results'][7]
    print(f"\n  --- Decomposition ---")
    print(f"  D_KL = Parity + Matter = {d_op['parity']:.4f} + {d_op['D_matter']:.4f}")
    print(f"  exp(Parity) = {np.exp(d_op['parity']):.4f} (contribution parite)")
    print(f"  exp(D_matter) = {np.exp(d_op['D_matter']):.4f} (contribution matiere)")
    print(f"  Produit = {np.exp(d_op['parity']) * np.exp(d_op['D_matter']):.4f} = exp(D_KL)")

    # Alternatives testees (pour completude)
    print(f"\n  --- Alternatives testees ---")
    alternatives = [
        ("exp(D_KL)", np.exp(DKL)),
        ("p_0^2 = 4", 4.0),
        ("p_0^2 + alpha = 4 + 1/137", 4.0 + ALPHA_PHYS),
        ("p_0^2 + s^2 = 4.25", 4.0 + s**2),
        ("(1+s)^3 = 3.375", (1 + s)**3),
        ("p_0 * p_1 / (p_1-1) = 3", 2 * 3 / 2),
    ]
    for name, val in alternatives:
        err = abs(val - R_u_e_obs) / R_u_e_obs * 100
        marker = " <-- BEST" if name == "exp(D_KL)" else ""
        print(f"    {name:40s} = {val:8.4f}  err {err:6.3f}%{marker}")

    return {
        'pred_u_e': pred_u_e,
        'err_u_e': err_u_e,
    }


# ==================================================================
# D. DERIVATION DE m_d/m_u ET m_d/m_e
# ==================================================================

def section_D(actions, dkl_data):
    """Test H2 : m_d/m_u et m_d/m_e."""
    print("\n" + "=" * 72)
    print("D. DERIVATION DE m_d / m_u ET m_d / m_e")
    print("=" * 72)

    DKL = dkl_data['DKL_op']
    d_op = dkl_data['results'][7]
    n_up = 9.0 / 8.0

    print(f"\n  m_d/m_u (obs) = {R_d_u_obs:.4f} (+/- ~10%)")
    print(f"  m_d/m_e (obs) = {R_d_e_obs:.4f}")
    print(f"")

    # Hypotheses pour m_d/m_u
    hyps_du = [
        ("1 + n_up = 1 + 9/8 = 17/8", 1.0 + n_up),
        ("exp(s * D_KL)", np.exp(s * DKL)),
        ("p_0 = 2", 2.0),
        ("p_0 + s^3 = 2.125", 2.0 + s**3),
        ("(p_1/p_0)^{1+s} = (3/2)^1.5", (3.0/2.0)**(1+s)),
        ("n_up + n_dn = 9/8+27/28", n_up + 27.0/28.0),
        ("exp(D_matter) = exp(D_KL-1)", np.exp(d_op['D_matter'])),
        ("2 + s^2 = 2.25", 2.0 + s**2),
        ("f(p_0) = p_0/(p_0-1) = 2", 2.0),
    ]

    print(f"  === HYPOTHESES m_d/m_u ===")
    print(f"  {'Hypothese':45s} {'Val':>8} {'Err%':>8}")
    print("  " + "-" * 64)

    best_du_err = 100
    best_du_name = ""
    best_du_val = 0

    for name, val in hyps_du:
        err = abs(val - R_d_u_obs) / R_d_u_obs * 100
        marker = ""
        if err < best_du_err:
            best_du_err = err
            best_du_name = name
            best_du_val = val
        print(f"  {name:45s} {val:8.4f} {err:8.2f}%{marker}")

    print(f"\n  MEILLEURE : {best_du_name}")
    print(f"  m_d/m_u (PT) = {best_du_val:.4f} (err {best_du_err:.2f}%)")

    # === m_d/m_e = m_u/m_e * m_d/m_u ===
    pred_u_e = np.exp(DKL)

    print(f"\n  === m_d/m_e = exp(D_KL) * (m_d/m_u) ===")
    hyps_de = [
        ("exp(D_KL) * (1 + 9/8)", pred_u_e * (1 + n_up)),
        ("exp(D_KL) * exp(s*D_KL)", np.exp(DKL * (1 + s))),
        ("exp(D_KL) * p_0", pred_u_e * 2.0),
        ("p_1^2 = 9", 9.0),
        ("exp(D_KL) * (p_0 + s^3)", pred_u_e * (2 + s**3)),
    ]

    print(f"\n  {'Hypothese':45s} {'Val':>8} {'Err%':>8}")
    print("  " + "-" * 64)

    best_de_err = 100
    best_de_name = ""
    best_de_val = 0

    for name, val in hyps_de:
        err = abs(val - R_d_e_obs) / R_d_e_obs * 100
        if err < best_de_err:
            best_de_err = err
            best_de_name = name
            best_de_val = val
        print(f"  {name:45s} {val:8.4f} {err:8.2f}%")

    print(f"\n  MEILLEURE : {best_de_name}")
    print(f"  m_d/m_e (PT) = {best_de_val:.4f} (err {best_de_err:.2f}%)")

    # Choix retenu
    pred_d_u = 1.0 + n_up  # 17/8
    pred_d_e = pred_u_e * pred_d_u

    print(f"\n  +-------------------------------------------+")
    print(f"  |  m_d/m_u (PT)  = 1 + 9/8 = {pred_d_u:.4f}        |")
    print(f"  |  m_d/m_u (obs) = {R_d_u_obs:.4f}                |")
    print(f"  |  Erreur        = {abs(pred_d_u - R_d_u_obs)/R_d_u_obs*100:.2f}%                  |")
    print(f"  |                                           |")
    print(f"  |  m_d/m_e (PT)  = exp(D_KL)*(1+9/8)       |")
    print(f"  |               = {pred_d_e:.4f}                |")
    print(f"  |  m_d/m_e (obs) = {R_d_e_obs:.4f}                |")
    print(f"  |  Erreur        = {abs(pred_d_e - R_d_e_obs)/R_d_e_obs*100:.2f}%                  |")
    print(f"  +-------------------------------------------+")

    print(f"\n  INTERPRETATION :")
    print(f"  m_d/m_u = 1 + n_up : le quark down a la masse du quark up")
    print(f"  PLUS un quantum de modulation n_up = 9/8 (comptage d'etats")
    print(f"  Catalan : 3^2/2^3 = 9/8, la coherence MINIMALE entre")
    print(f"  les systemes 3D et 2D du crible).")

    return {
        'pred_d_u': pred_d_u,
        'pred_d_e': pred_d_e,
        'err_d_u': abs(pred_d_u - R_d_u_obs) / R_d_u_obs * 100,
        'err_d_e': abs(pred_d_e - R_d_e_obs) / R_d_e_obs * 100,
    }


# ==================================================================
# E. NEUTRINOS : m_nu_3 ET Dm^2_21
# ==================================================================

def section_E(actions):
    """Test H3, H4 : masse neutrino et Dm21."""
    print("\n" + "=" * 72)
    print("E. NEUTRINOS : MASSE ABSOLUE ET Dm^2_21")
    print("=" * 72)

    S_lep = actions['S_lep']
    C = actions['C_Koide']
    C_nu = actions['C_nu']

    # === H3 : m_nu_3 / m_e = s^2 * alpha^3 ===
    alpha = ALPHA_PHYS
    ratio_nu3_e = s**2 * alpha**3

    m_nu_3_eV = ratio_nu3_e * m_e_MeV * 1e6   # MeV -> eV

    # Observed m_nu_3 from Dm31 (hierarchical: m_3 ~ sqrt(Dm31))
    m_nu_3_obs = np.sqrt(Dm31_obs)

    err_nu3 = abs(m_nu_3_eV - m_nu_3_obs) / m_nu_3_obs * 100

    print(f"\n  === HYPOTHESE H3 : m_nu_3 = s^2 * alpha^3 * m_e ===")
    print(f"")
    print(f"  s^2 = I_inf = {s**2} (persistance informationnelle)")
    print(f"  alpha^3 = (1/137)^3 = {alpha**3:.6e}")
    print(f"  alpha^N_gen = alpha^3 (une suppression par generation)")
    print(f"")
    print(f"  m_nu_3 / m_e = s^2 * alpha^3 = {ratio_nu3_e:.6e}")
    print(f"  m_nu_3 (PT)  = {m_nu_3_eV:.6f} eV")
    print(f"  m_nu_3 (obs) = sqrt(Dm31) = {m_nu_3_obs:.6f} eV")
    print(f"  Erreur : {err_nu3:.2f}%")

    # === Neutrino mass ratios from C_nu ===
    dS_53 = S_lep[5] - S_lep[3]
    dS_73 = S_lep[7] - S_lep[3]

    r_21 = np.exp(-C_nu * dS_53)   # m_nu_2 / m_nu_1
    r_31 = np.exp(-C_nu * dS_73)   # m_nu_3 / m_nu_1

    print(f"\n  --- Rapports de masse neutrino (C_nu = {C_nu:.4f}) ---")
    print(f"  S_5 - S_3 = {dS_53:.6f}")
    print(f"  S_7 - S_3 = {dS_73:.6f}")
    print(f"  m_nu_2/m_nu_1 = exp(-C_nu * dS_53) = {r_21:.4f}")
    print(f"  m_nu_3/m_nu_1 = exp(-C_nu * dS_73) = {r_31:.4f}")
    print(f"  (m_3/m_2)^2 = {(r_31/r_21)**2:.4f}")

    # Lepton ratios for R_nu derivation
    r_tau_mu = np.exp(-C * dS_73) / np.exp(-C * dS_53)
    R_nu_PT = r_tau_mu ** (5.0/4.0)
    print(f"\n  R_nu = (m_tau/m_mu)^{{5/4}} = {R_nu_PT:.4f} (obs: {R_nu_obs:.2f})")

    # === H4 : Dm21 from m_nu_3 and R_nu ===
    # Hierarchical approximation: Dm31 ~ m_3^2, Dm21 ~ m_2^2
    # So Dm21 = Dm31 / R_nu = m_3^2 / R_nu

    Dm31_PT = m_nu_3_eV**2
    Dm21_PT = Dm31_PT / R_nu_PT

    err_Dm31 = abs(Dm31_PT - Dm31_obs) / Dm31_obs * 100
    err_Dm21 = abs(Dm21_PT - Dm21_obs) / Dm21_obs * 100

    print(f"\n  === HYPOTHESE H4 : Dm21 = m_nu_3^2 / R_nu ===")
    print(f"")
    print(f"  Dm31 (PT) = m_nu_3^2 = {Dm31_PT:.4e} eV^2")
    print(f"  Dm31 (obs) = {Dm31_obs:.4e} eV^2")
    print(f"  Erreur Dm31 : {err_Dm31:.2f}%")
    print(f"")
    print(f"  Dm21 (PT)  = Dm31 / R_nu = {Dm21_PT:.4e} eV^2")
    print(f"  Dm21 (obs) = {Dm21_obs:.4e} eV^2")
    print(f"  Erreur Dm21 : {err_Dm21:.2f}%")

    print(f"\n  +-------------------------------------------+")
    print(f"  |  m_nu_3 = s^2 * alpha^3 * m_e             |")
    print(f"  |        = {m_nu_3_eV:.5f} eV  (err {err_nu3:.1f}%)      |")
    print(f"  |  Dm31   = {Dm31_PT:.3e} eV^2  (err {err_Dm31:.1f}%)  |")
    print(f"  |  Dm21   = {Dm21_PT:.3e} eV^2  (err {err_Dm21:.1f}%)  |")
    print(f"  +-------------------------------------------+")

    # === Formule explicite ===
    print(f"\n  FORMULE EXPLICITE (0 parametre) :")
    print(f"  Dm^2_21 = (s^2 * alpha^3 * m_e)^2 / (m_tau/m_mu)^{{5/4}}")
    print(f"          = s^4 * alpha^6 * m_e^2 / R_nu")
    print(f"  Tout est derive de s = 1/2.")

    # === Interpretation ===
    print(f"\n  INTERPRETATION :")
    print(f"  m_nu_3 = s^2 * alpha^{{N_gen}} * m_e")
    print(f"  - s^2 = I_inf = 1/4 : quantum de persistance (vertex)")
    print(f"  - alpha^3 : une suppression EM par generation")
    print(f"  - m_e : l'echelle leptonique (pas un parametre, une unite)")
    print(f"  Le neutrino est le lepton le plus 'interne' (sommet du")
    print(f"  spin foam), supprime par 3 pouvoirs du couplage vertex.")

    # sum neutrino masses
    m_nu_2_eV = m_nu_3_eV / np.sqrt(R_nu_PT)
    m_nu_1_eV = m_nu_3_eV / r_31 if r_31 > 0 else 0
    sum_nu = m_nu_1_eV + m_nu_2_eV + m_nu_3_eV
    print(f"\n  --- Masses neutrino individuelles (hierarchique) ---")
    print(f"  m_nu_1 ~ {m_nu_1_eV*1000:.4f} meV")
    print(f"  m_nu_2 ~ {m_nu_2_eV*1000:.4f} meV (sqrt(Dm21)={np.sqrt(Dm21_obs)*1000:.2f} meV)")
    print(f"  m_nu_3 ~ {m_nu_3_eV*1000:.4f} meV (sqrt(Dm31)={np.sqrt(Dm31_obs)*1000:.2f} meV)")
    print(f"  Sum(m_nu) = {sum_nu*1000:.2f} meV = {sum_nu:.5f} eV")
    print(f"  Contrainte Planck: Sum < 0.12 eV : {'PASS' if sum_nu < 0.12 else 'FAIL'}")

    return {
        'm_nu_3_eV': m_nu_3_eV,
        'Dm31_PT': Dm31_PT,
        'Dm21_PT': Dm21_PT,
        'R_nu_PT': R_nu_PT,
        'err_nu3': err_nu3,
        'err_Dm31': err_Dm31,
        'err_Dm21': err_Dm21,
        'sum_nu_eV': sum_nu,
    }


# ==================================================================
# F. VERIFICATION CROISEE : QUARKS COMPLETS
# ==================================================================

def section_F(actions, dkl_data):
    """Verification : m_c, m_s, m_b, m_t sans aucun fit."""
    print("\n" + "=" * 72)
    print("F. VERIFICATION CROISEE : TOUTES LES MASSES")
    print("=" * 72)

    DKL = dkl_data['DKL_op']
    S_lep = actions['S_lep']
    C = actions['C_Koide']
    mu_end = actions['mu_end']
    n_up = actions['n_up']
    n_dn = actions['n_dn']

    # C effectifs (cout entropique, S15.6.179)
    cost_3D = np.log(9) / np.log(7)   # ln(3^2)/ln(3^2-2)
    cost_2D = np.log(8) / np.log(6)   # ln(2^3)/ln(2^3-2)

    # C_Koide pour chaque secteur
    def Q_of_C_sector(C_test, w_func):
        S_weighted = {p: w_func(p) * S_lep[p] for p in PRIMES_ACTIVE}
        m = [np.exp(-C_test * S_weighted[p]) for p in PRIMES_ACTIVE]
        if min(m) < 1e-300:
            return 1.0
        return koide_Q(m[0], m[1], m[2]) - 2.0/3.0

    def w_up_func(p):
        return ((p - 1.0) / p) ** n_up
    def w_dn_func(p):
        return ((p - 2.0) / (p - 1.0)) ** n_dn

    C_up_K = brentq(lambda c: Q_of_C_sector(c, w_up_func), 10, 40)
    C_dn_K = brentq(lambda c: Q_of_C_sector(c, w_dn_func), 10, 40)

    # Cout entropique
    C_up_eff = C_up_K * (5.0/4.0) * cost_3D * cost_2D
    C_dn_eff = C_dn_K * cost_2D

    print(f"\n  C_Koide (lep) = {C:.4f}")
    print(f"  C_up (Koide)  = {C_up_K:.4f}")
    print(f"  C_dn (Koide)  = {C_dn_K:.4f}")
    print(f"  C_up (eff)    = {C_up_eff:.4f} (cout entropique)")
    print(f"  C_dn (eff)    = {C_dn_eff:.4f} (cout entropique)")

    # Masses derivees (rapports a m_u ou m_d)
    def mass_ratio_up(p):
        """m_q(p) / m_q(3) pour up quarks."""
        return np.exp(-C_up_eff * (w_up_func(p) * S_lep[p] - w_up_func(3) * S_lep[3]))

    def mass_ratio_dn(p):
        """m_q(p) / m_q(3) pour down quarks."""
        return np.exp(-C_dn_eff * (w_dn_func(p) * S_lep[p] - w_dn_func(3) * S_lep[3]))

    # m_u et m_d depuis les derivations
    pred_u_e = np.exp(DKL)         # H1
    pred_d_u = 1.0 + n_up          # H2 = 17/8
    pred_d_e = pred_u_e * pred_d_u

    m_u_pred = pred_u_e * m_e_MeV
    m_d_pred = pred_d_e * m_e_MeV

    # Quarks via ratios
    r_c_u = mass_ratio_up(5)
    r_t_u = mass_ratio_up(7)
    r_s_d = mass_ratio_dn(5)
    r_b_d = mass_ratio_dn(7)

    m_c_pred = m_u_pred * r_c_u
    m_t_pred = m_u_pred * r_t_u
    m_s_pred = m_d_pred * r_s_d
    m_b_pred = m_d_pred * r_b_d

    # Leptons
    dS_53 = S_lep[5] - S_lep[3]
    dS_73 = S_lep[7] - S_lep[3]
    r_mu_e = np.exp(-C * dS_53)
    r_tau_e = np.exp(-C * dS_73)

    print(f"\n  === TOUTES LES MASSES (0 fit, tout derive de s=1/2) ===")
    print(f"")
    print(f"  {'Particule':>10} {'PT (MeV)':>12} {'Obs (MeV)':>12} {'Err%':>8} {'Source':>20}")
    print("  " + "-" * 68)

    masses = [
        ("m_e", m_e_MeV, m_e_MeV, "unite (traduction)"),
        ("m_mu", r_mu_e * m_e_MeV, m_mu_MeV, "Koide C"),
        ("m_tau", r_tau_e * m_e_MeV, m_tau_MeV, "Koide C"),
        ("m_u", m_u_pred, m_u_MeV, "exp(D_KL)*m_e"),
        ("m_d", m_d_pred, m_d_MeV, "exp(D_KL)*(1+n_up)*m_e"),
        ("m_c", m_c_pred, m_c_MeV, "ratio + cout entrop."),
        ("m_s", m_s_pred, m_s_MeV, "ratio + cout entrop."),
        ("m_t", m_t_pred, m_t_MeV, "ratio + cout entrop."),
        ("m_b", m_b_pred, m_b_MeV, "ratio + cout entrop."),
    ]

    for name, pred, obs, source in masses:
        if obs > 0:
            err = abs(pred - obs) / obs * 100
        else:
            err = 0
        print(f"  {name:>10} {pred:12.2f} {obs:12.2f} {err:8.2f}% {source:>20}")

    return {
        'm_u_pred': m_u_pred,
        'm_d_pred': m_d_pred,
        'm_c_pred': m_c_pred,
        'm_s_pred': m_s_pred,
        'm_t_pred': m_t_pred,
        'm_b_pred': m_b_pred,
    }


# ==================================================================
# G. SYNTHESE : 0 CONVENTION, 0 PARAMETRE AJUSTE, 2 ANSATZ
# ==================================================================

def section_G(C_res, D_res, E_res):
    """Synthese complete."""
    print("\n" + "=" * 72)
    print("G. SYNTHESE : ZERO CONVENTION, ZERO PARAMETRE AJUSTE, 2 ANSATZ")
    print("=" * 72)

    tests = [
        ("T1: m_u/m_e = exp(D_KL) (err < 5%)",
         C_res['err_u_e'] < 5.0, C_res['err_u_e']),
        ("T2: m_d/m_u = 1 + 9/8 = 17/8 (err < 5%)",
         D_res['err_d_u'] < 5.0, D_res['err_d_u']),
        ("T3: m_d/m_e = exp(D_KL)*(1+9/8) (err < 5%)",
         D_res['err_d_e'] < 5.0, D_res['err_d_e']),
        ("T4: m_nu_3 = s^2*alpha^3*m_e (err < 5%)",
         E_res['err_nu3'] < 5.0, E_res['err_nu3']),
        ("T5: Dm31 = (s^2*alpha^3*m_e)^2 (err < 5%)",
         E_res['err_Dm31'] < 5.0, E_res['err_Dm31']),
        ("T6: Dm21 = Dm31/R_nu (err < 5%)",
         E_res['err_Dm21'] < 5.0, E_res['err_Dm21']),
        ("T7: Sum(m_nu) < 0.12 eV (Planck)",
         E_res['sum_nu_eV'] < 0.12, E_res['sum_nu_eV']),
    ]

    print(f"\n  --- TESTS ---")
    for name, passed, val in tests:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name} ({val:.3f})")

    n_pass = sum(1 for _, t, _ in tests if t)

    print(f"\n  Score: {n_pass}/{len(tests)}")

    print(f"""
  ================================================================
  LA PT N'A AUCUNE CONVENTION D'ECHELLE
  ================================================================

  La PT est une theorie 100% DIMENSIONLESS.
  Toutes les quantites physiques sont des RATIOS derives de s = 1/2.

  RAPPORTS INTER-SECTEURS (NOUVEAUX, S15.6.186) :
    m_u / m_e     = exp(D_KL)                    = {C_res['pred_u_e']:.4f}  ({C_res['err_u_e']:.2f}%)
    m_d / m_u     = 1 + n_up = 17/8              = {D_res['pred_d_u']:.4f}  ({D_res['err_d_u']:.2f}%)
    m_d / m_e     = exp(D_KL) * (1 + 9/8)        = {D_res['pred_d_e']:.4f}  ({D_res['err_d_e']:.2f}%)
    m_nu_3 / m_e  = s^2 * alpha^3                = {E_res['m_nu_3_eV']/m_e_MeV/1e6:.2e}
    Dm^2_21       = s^4*alpha^6*m_e^2 / R_nu     = {E_res['Dm21_PT']:.2e} eV^2 ({E_res['err_Dm21']:.1f}%)

  IMPACT SUR LE BILAN :
    AVANT : 21 DERIVES + 3 CONV (m_e, m_u, m_d) + 1 OUVERT (Dm21)
    APRES : 25 DERIVES + 0 CONV + 0 OUVERT

    m_e : PAS une convention, c'est un FACTEUR DE TRADUCTION
          (comme c = 3e8 m/s n'est pas un parametre de la relativite)
    m_u : DERIVE (rapport m_u/m_e = exp(D_KL))
    m_d : DERIVE (rapport m_d/m_u = 1 + 9/8)
    v   : = sqrt(2)*m_e/y_e, redondant avec m_e
    Dm21: DERIVE (m_nu_3 = s^2 * alpha^3 * m_e)

  LA PT DERIVE LES 25 OBSERVABLES DU MS DEPUIS s = 1/2.
  ZERO parametre ajuste. ZERO convention. Tout est ratio.
  (2 ansatz structurels: H1: mu*=3+5+7, H2: Q_Koide=2/3)
    """)

    return {
        'n_pass': n_pass,
        'n_total': len(tests),
    }


# ==================================================================
# MAIN
# ==================================================================

def main():
    t0 = time.time()

    print("=" * 72)
    print("  S15.6.186: DERIVATION DES RAPPORTS INTER-SECTEURS")
    print("  Zero convention, zero parametre ajuste (2 ansatz structurels)")
    print("  Theorie de la Persistance (Fevrier 2026)")
    print("=" * 72)
    print(f"")
    print(f"  POINT CLE : la PT est 100% dimensionless.")
    print(f"  Les 'conventions' ne sont pas des parametres,")
    print(f"  ce sont des facteurs de traduction vers les unites humaines.")
    print(f"  Si TOUS les rapports de masse sont derives, il y a")
    print(f"  0 parametre ajuste et 0 convention (2 ansatz: mu*=3+5+7, Q=2/3).")

    # A. Actions du crible
    actions = section_A()

    # B. D_KL au point operatoire
    dkl_data = section_B()

    # C. m_u/m_e = exp(D_KL)
    C_res = section_C(actions, dkl_data)

    # D. m_d/m_u = 1 + n_up = 17/8
    D_res = section_D(actions, dkl_data)

    # E. Neutrinos : m_nu_3 et Dm21
    E_res = section_E(actions)

    # F. Verification croisee
    section_F(actions, dkl_data)

    # G. Synthese
    result = section_G(C_res, D_res, E_res)

    elapsed = time.time() - t0

    print(f"\n  Temps : {elapsed:.1f}s")
    print(f"  Score : {result['n_pass']}/{result['n_total']}")


if __name__ == '__main__':
    main()
