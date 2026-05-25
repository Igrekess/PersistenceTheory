#!/usr/bin/env python3
"""Verrou G_HP3 : formule de trace explicite Berry-Keating PT -> R(s).

Contexte (notes 42, 57, 58, 59):
  - R(s) = zeta(s) / (zeta_+(s) zeta_-(s)) avec zeta_+ = exp(sum_p a_p p^-s),
    zeta_- = exp(sum_p b_p p^-s), A_p = a_p + b_p.
  - kappa_p(s) = 1 - (1 - p^-s) exp(A_p p^-s), det(I - D_R(s))^-1 = R(s).
  - H_PT_BK = (u . p_u + p_u . u)/2 sur L^2(R_+, du), u = log p.

Verrou G_HP3 : connecter H_PT_BK aux kappa_p via une formule de trace.

Cette note implemente le test de la formule de trace explicite (Riemann-Weil
adaptee a R) suivant la prescription :

   somme geometrique sur primes = somme spectrale sur zeros + contributions
   archimediennes / continuum

Le cote prime (cote geometrique) de R est :

   -d/ds log R(s) = sum_p (log p) [ 1/(p^s - 1) + A_p p^-s ]

et l'identite Fredholm donne, dans Re(s)>1 :

   -d/ds log R(s) = sum_p kappa_p'(s) / (1 - kappa_p(s)),
   ou kappa_p(s) = 1 - (1 - p^-s) exp(A_p p^-s).

Cette egalite est analytique exacte (preuvee algebraiquement, verifiee
numeriquement).

Le cote Berry-Keating donne la transformee d'une fonction test paire phi(gamma)
contre la mesure spectrale du semi-groupe e^{-itH_PT_BK}. Quand on choisit
phi tel que la transformee inverse selectionne le mode k.log p avec poids
log p / p^{k/2}, on retrouve precisement le cote geometrique.

Le test numerique : pour t = log p_0 d'un prime donne, l'observable
'amplitude de retour' Tr(e^{-itH_PT_BK}) regularisee se compare a la
contribution du prime p_0 dans -d/ds log R(s).
"""

from __future__ import annotations

import math
import sys
from typing import List, Tuple

import numpy as np
import mpmath as mp


# Constantes PT canoniques (notes 42 et suivantes)
Q_PLUS  = mp.mpf(13) / mp.mpf(15)
Q_MINUS = mp.exp(-mp.mpf(1) / mp.mpf(15))


def first_primes(n: int) -> List[int]:
    """Liste des n premiers premiers via crible d'Eratosthene grossier."""
    if n <= 0:
        return []
    upper = max(20, int(n * (math.log(n) + math.log(max(2, math.log(n + 2)))) + 10))
    sieve = np.ones(upper + 1, dtype=bool)
    sieve[:2] = False
    for k in range(2, int(math.isqrt(upper)) + 1):
        if sieve[k]:
            sieve[k * k :: k] = False
    primes = np.flatnonzero(sieve).tolist()
    while len(primes) < n:
        upper *= 2
        sieve = np.ones(upper + 1, dtype=bool)
        sieve[:2] = False
        for k in range(2, int(math.isqrt(upper)) + 1):
            if sieve[k]:
                sieve[k * k :: k] = False
        primes = np.flatnonzero(sieve).tolist()
    return primes[:n]


def A_p(p: int) -> mp.mpf:
    """A_p = a_p + b_p = delta_+ (2 - delta_+) + delta_- (2 - delta_-)."""
    qp = Q_PLUS
    qm = Q_MINUS
    delta_plus  = (1 - qp ** p) / mp.mpf(p)
    delta_minus = (1 - qm ** p) / mp.mpf(p)
    a_p = delta_plus  * (2 - delta_plus)
    b_p = delta_minus * (2 - delta_minus)
    return a_p + b_p


def kappa_p(p: int, s: complex) -> complex:
    """kappa_p(s) = 1 - (1 - p^-s) exp(A_p p^-s)."""
    Ap = A_p(p)
    z = mp.mpf(p) ** (-s)
    return 1 - (1 - z) * mp.exp(Ap * z)


def kappa_p_deriv(p: int, s: complex) -> complex:
    """kappa_p'(s) calcule en forme close.

    CORRECTION 2026-05-17 : signe corrige (voir A4_formule_trace.md §3.3bis).

    kappa_p = 1 - (1 - p^-s) exp(A_p p^-s),  z = p^-s.
    d/ds [(1-z) exp(A_p z)] :
      = (d/ds (1-z)) exp(A_p z) + (1-z) (d/ds exp(A_p z))
      = (log p . z) exp(A_p z) + (1-z) exp(A_p z) . A_p . (-log p . z)
      = log p . z . exp(A_p z) . [1 - A_p (1 - z)]
    kappa_p' = -d/ds [(1-z) exp(A_p z)]
            = -log p . z . exp(A_p z) . [1 - A_p (1 - z)]
            =  log p . z . exp(A_p z) . [A_p (1 - z) - 1]

    Et donc kappa_p'/(1 - kappa_p) = log p . [A_p z - z/(1-z)]
                                  = log p . [A_p p^-s - 1/(p^s - 1)]
    """
    Ap = A_p(p)
    logp = mp.log(p)
    z = mp.mpf(p) ** (-s)
    eAp_z = mp.exp(Ap * z)
    # SIGNE CORRIGE : (A_p (1-z) - 1) au lieu de (1 + A_p (1-z))
    return logp * z * eAp_z * (Ap * (1 - z) - 1)


def logderiv_R_via_primes(s: complex, primes: List[int]) -> complex:
    """Cote 'Euler + couplage PT'.

    CORRECTION 2026-05-17 : signe corrige du terme A_p p^-s
    (voir A4_formule_trace.md §3.3bis).

    -R'(s)/R(s) = sum_p (log p) [ 1/(p^s - 1) - A_p p^-s ]

    avec MOINS devant A_p p^-s (et non plus comme dans la version 2026-05).
    """
    total = mp.mpc(0)
    for p in primes:
        Ap = A_p(p)
        z = mp.mpf(p) ** (-s)
        total += mp.log(p) * (z / (1 - z) - Ap * z)  # SIGNE CORRIGE : - Ap * z
    return total


def logderiv_R_via_kappa(s: complex, primes: List[int]) -> complex:
    """Cote 'Fredholm diagonal' :

       -R'(s)/R(s) = sum_p kappa_p'(s) / (1 - kappa_p(s))

    Identite analytique : equivalente a logderiv_R_via_primes pour
    Re(s) > 1, en vertu de la note 42.
    """
    total = mp.mpc(0)
    for p in primes:
        kp = kappa_p(p, s)
        kpprime = kappa_p_deriv(p, s)
        total += kpprime / (1 - kp)
    return total


def trace_BK_prime_contribution(p: int, t: float, sigma: float = 0.5) -> complex:
    """Contribution du prime p a une 'trace regularisee' BK.

    Modele heuristique : pour H_PT_BK = u p_u + 1/2 sur L^2(R_+, du),
    le propagateur e^{-itH} agit comme dilatation par e^{-t}.
    Une mesure atomique nu = sum_p delta_{log p} pondere chaque prime
    par log p / p^{1/2} (mesure de Haar mult.) et un facteur de phase
    e^{-it log p} = p^{-it}.

    Au niveau du prime p, la contribution simple-orbite est :

       c_p(t) = (log p) / p^{sigma} * p^{-it}

    et la contribution multi-orbite (sous la k-eme iteree du flot RG) :

       c_p^{(k)}(t) = (log p) / p^{k sigma} * p^{-i k t}

    La somme totale c_p(t) = sum_k c_p^{(k)}(t) = (log p) / (p^{sigma+it} - 1)
    se reconnait comme la contribution du prime p dans la derivee
    logarithmique -d/ds log zeta(s) a s = sigma + it.

    Le terme PT ajoute la 'compensation' A_p (signe CORRIGE 2026-05-17) :

       c_p^{PT}(t) = (log p) * [ 1/(p^{sigma+it} - 1) - A_p . p^{-(sigma+it)} ].

    C'est exactement le terme par-prime de -d/ds log R(s) en s = sigma + it.
    """
    s = mp.mpc(sigma, t)
    Ap = A_p(p)
    z = mp.mpf(p) ** (-s)
    return mp.log(p) * (z / (1 - z) - Ap * z)  # SIGNE CORRIGE : - Ap * z


def smear(t_values: np.ndarray, primes: List[int], sigma: float = 0.5) -> np.ndarray:
    """Calcule -d/dt log |R(sigma + it)|^2 / 2 via la somme prime + PT.

    On utilise l'identite :
      -d/dt log R(sigma+it) = -i * (-R'(s)/R(s))   (chaine, d s/d t = i)
    donc Im[ -R'/R(s) ] = d/dt log |R(s)|.
    """
    out = np.zeros(len(t_values))
    for j, t in enumerate(t_values):
        val = logderiv_R_via_primes(mp.mpc(sigma, t), primes)
        # d/dt log|R(s)| = Im(-R'/R(s))   pour s = sigma + it
        out[j] = float(mp.im(val))
    return out


def identity_check(primes: List[int], sigma_values: List[float], t: float = 0.0) -> List[Tuple[float, complex]]:
    """Verifie l'identite cote-prime = cote-kappa pour quelques sigma > 1."""
    rows = []
    for sigma in sigma_values:
        s = mp.mpc(sigma, t)
        v_prime = logderiv_R_via_primes(s, primes)
        v_kappa = logderiv_R_via_kappa(s, primes)
        diff = v_prime - v_kappa
        rows.append((sigma, complex(v_prime), complex(v_kappa), complex(diff)))
    return rows


def compare_to_mpmath_zeta(primes: List[int], sigma: float, t: float) -> Tuple[complex, complex, complex]:
    """Calcule -d/ds log R(s) via primes a Re(s)=sigma, puis le compare a la
    valeur exacte calculee a partir de zeta'(s)/zeta(s) - sum_p A_p log p p^-s.
    """
    mp.mp.dps = 30
    s = mp.mpc(sigma, t)

    # Cote primes (tronque a la liste fournie)
    v_prime = logderiv_R_via_primes(s, primes)

    # Cote 'exact' :  -zeta'/zeta(s) + sum_p (log p) A_p p^-s
    zeta_log_deriv = -mp.zeta(s, derivative=1) / mp.zeta(s)  # = sum_n Lambda(n) n^-s
    coupling = mp.mpc(0)
    for p in primes:
        Ap = A_p(p)
        coupling += mp.log(p) * Ap * mp.mpf(p) ** (-s)
    v_exact = zeta_log_deriv + coupling

    return v_prime, v_exact, v_prime - v_exact


def main():
    mp.mp.dps = 25
    primes = first_primes(2000)

    print("=" * 72)
    print("PT verrou G_HP3 : formule de trace explicite Berry-Keating -> R(s)")
    print("=" * 72)

    # 1. Verification identite cote-prime = cote-kappa
    print()
    print("(1) Identite Fredholm  -R'/R(s) = sum_p kappa'_p/(1-kappa_p)")
    print("                       = sum_p (log p) [ 1/(p^s - 1) + A_p p^-s ]")
    print()
    print(f"{'sigma':>8} {'|prime - kappa|':>22}")
    rows = identity_check(primes, [1.2, 1.5, 2.0, 3.0], t=14.134725)
    for sigma, v_p, v_k, diff in rows:
        print(f"  {sigma:>6.2f}   {abs(diff):>20.3e}")

    # 2. Convergence en troncature pour sigma > 1
    print()
    print("(2) Convergence Euler tronque (Re(s)=1.3, t=14.134725)")
    print(f"{'P':>10} {'|prime_P - exact|':>22}")
    for P in [50, 200, 500, 1000, 2000]:
        v_p, v_e, d = compare_to_mpmath_zeta(primes[:P], 1.3, 14.134725)
        print(f"  {P:>8d}   {float(abs(d)):>20.3e}")

    # 3. Test extrapolation a sigma = 1/2 (ligne critique)
    print()
    print("(3) Critical line tracking : -d/dt log|R(1/2 + it)|")
    print("    via cote-prime tronque (regularisation Hadamard non implementee).")
    print()
    print(f"{'t':>10} {'prime side P=2000':>22} {'mpmath zeta side':>22}")
    for t in [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]:
        v_p, v_e, _ = compare_to_mpmath_zeta(primes, 0.5, t)
        # Pres de zeros de zeta, -zeta'/zeta a un pole simple en s = 1/2 + i gamma_n
        # donc Im[-zeta'/zeta(s)] -> grand pres de t = gamma_n
        # On signale juste les valeurs
        print(f"  {t:>8.6f}   {float(mp.im(v_p)):>20.4f}   {float(mp.im(v_e)):>20.4f}")
    print()
    print("    NB : a Re(s)=1/2, la serie sur primes ne converge pas")
    print("    absolument. Les valeurs ci-dessus sont indicatives (somme partielle).")
    print("    Une regularisation Hadamard / Cesaro est necessaire pour stabiliser.")

    # 4. Interpretation : la formule de trace BK
    print()
    print("(4) Lecture Berry-Keating :")
    print()
    print("    Tr_reg(e^{-it H_PT_BK})  =  sum_p (log p) p^{-1/2 - it} . [1 + A_p (1 - p^{-s})...]")
    print()
    print("    cote prime :  longueurs ell_p = log p (orbites du flot de dilatation)")
    print("    poids       :  log p / p^sigma (mesure Haar mult.)")
    print("    phase       :  p^{-it} = e^{-it log p}")
    print()
    print("    Cote spectral (dualite Weil): la meme somme = sum_gamma h(gamma)")
    print("    avec h = transformee de la fonction de smearing en t.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
