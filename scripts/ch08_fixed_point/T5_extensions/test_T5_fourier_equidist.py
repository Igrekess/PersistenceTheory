"""
test_T5_fourier_equidist.py -- Decomposition de Fourier du boundary CRT
Contribution S15.6.302

OBJECTIF : Explorer la route trigonometrique vers le theoreme d'equidistribution.

IDEE CENTRALE :
  L'indicatrice de retrait chi(s) = 1_{s=0 mod p} se decompose en :
    chi(s) = (1/p) * sum_{t=0}^{p-1} omega^{t*s}   (omega = e^{2*pi*i/p})

  Le terme t=0 donne la prediction uniforme (1/p).
  Les termes t >= 1 sont les corrections oscillatoires.

  La sequence etendue (coprime a P_k dans [1, P_k*p]) est p COPIES EXACTES
  du pattern de niveau k. Cela implique des annulations specifiques dans
  les sommes exponentielles.

STRUCTURE :
  1. Calcul exact du boundary par position (contribution locale 4-gram)
  2. Decomposition de Fourier de l'indicatrice de retrait
  3. Verification : sum S(t) pour t>=1 = correction oscillatoire
  4. Recherche de patterns sin^2
  5. Analyse : le boundary est-il entierement determine par le niveau k ?
"""

import numpy as np
from math import prod, gcd
import unittest

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]


def sieve_survivors(prime_list):
    """Retourne les survivants dans [1, P] coprime a tous les p dans prime_list."""
    P = prod(prime_list)
    sieve = np.ones(P + 1, dtype=np.bool_)
    sieve[0] = False
    for p in prime_list:
        sieve[::p] = False
    return np.flatnonzero(sieve)


def gap_classes(survivors, P):
    """Calcule les classes de gaps (mod 3) pour une sequence cyclique."""
    n = len(survivors)
    gaps = np.empty(n, dtype=np.int64)
    gaps[:-1] = survivors[1:] - survivors[:-1]
    gaps[-1] = P + survivors[0] - survivors[-1]
    return gaps % 3


def compute_boundary_per_position(classes):
    """Pour chaque position j, calcule Delta_j = changement de n00 si j est retire.

    Quand on retire la position j :
    - gaps c_{j-1} et c_j fusionnent en c_new = (c_{j-1} + c_j) mod 3
    - 3 transitions sont affectees : a j-2, j-1, j (la transition a j disparait)

    Delta_j depend du 4-gramme (c_{j-2}, c_{j-1}, c_j, c_{j+1}).
    """
    N = len(classes)
    delta = np.zeros(N, dtype=np.int64)

    for j in range(N):
        cjm2 = classes[(j - 2) % N]
        cjm1 = classes[(j - 1) % N]
        cj   = classes[j]
        cjp1 = classes[(j + 1) % N]

        c_new = (cjm1 + cj) % 3

        # Transition a j-2 : (c_{j-2}, c_{j-1}) -> (c_{j-2}, c_new)
        old_00_jm2 = int(cjm2 == 0 and cjm1 == 0)
        new_00_jm2 = int(cjm2 == 0 and c_new == 0)

        # Transition a j-1 : (c_{j-1}, c_j) -> (c_new, c_{j+1})
        old_00_jm1 = int(cjm1 == 0 and cj == 0)
        new_00_jm1 = int(c_new == 0 and cjp1 == 0)

        # Transition a j : (c_j, c_{j+1}) -> DISPARAIT
        old_00_j = int(cj == 0 and cjp1 == 0)

        delta[j] = (new_00_jm2 - old_00_jm2) + (new_00_jm1 - old_00_jm1) - old_00_j

    return delta


# ================================================================
# Exploration pour k=3..7
# ================================================================

results = {}

for k in range(3, 8):
    pk = PRIMES[:k]  # premiers p_1, ..., p_k
    p_next = PRIMES[k]  # p_{k+1}
    P_k = prod(pk)
    P_ext = P_k * p_next

    if P_ext > 500_000_000:
        break

    # Niveau k : survivants et classes
    surv_k = sieve_survivors(pk)
    N_k = len(surv_k)
    cls_k = gap_classes(surv_k, P_k)
    n00_k = int(np.sum((cls_k == 0) & (np.roll(cls_k, -1) == 0)))

    # Sequence etendue : coprime a p_1...p_k dans [1, P_k*p]
    surv_ext = sieve_survivors(pk)  # On reconstruit dans [1, P_ext]
    # Plus efficace : p copies de surv_k
    surv_ext = np.concatenate([surv_k + m * P_k for m in range(p_next)])
    surv_ext.sort()
    M = len(surv_ext)  # = p * N_k
    cls_ext = gap_classes(surv_ext, P_ext)
    n00_ext = int(np.sum((cls_ext == 0) & (np.roll(cls_ext, -1) == 0)))

    # Positions retirees : celles divisibles par p_next
    removed_mask = (surv_ext % p_next == 0)
    removed_idx = np.flatnonzero(removed_mask)
    N_removed = len(removed_idx)

    # Niveau k+1 : survivants exacts
    surv_k1 = sieve_survivors(PRIMES[:k+1])
    N_k1 = len(surv_k1)
    cls_k1 = gap_classes(surv_k1, P_ext)
    n00_k1 = int(np.sum((cls_k1 == 0) & (np.roll(cls_k1, -1) == 0)))

    # Boundary exact (du code S15.6.301)
    bnd_exact = n00_k1 - (p_next - 3) * n00_k

    # Calcul du Delta par position sur la sequence etendue
    delta_ext = compute_boundary_per_position(cls_ext)

    # Somme des Delta aux positions retirees = changement REEL de n00
    sum_delta_removed = int(delta_ext[removed_idx].sum())

    # Verification : n00_{k+1} = n00_ext + sum_delta_removed ?
    # (si les retraits n'interagissent pas)
    n00_predicted = n00_ext + sum_delta_removed

    # Prediction uniforme VRAIE : (1/p) * sum de TOUS les delta
    sum_delta_all = int(delta_ext.sum())
    bnd_uniform_true = sum_delta_all / p_next

    # ===== DECOMPOSITION DE FOURIER =====
    omega = np.exp(2j * np.pi / p_next)

    # S(t) = sum_j omega^{t*s_j} * delta_j  pour j dans les positions retirees...
    # Non, S(t) = sum sur TOUTES les positions j de omega^{t*s_j} * delta_j
    # (car chi_j = (1/p) sum_t omega^{t*s_j})

    S = np.zeros(p_next, dtype=np.complex128)
    for t in range(p_next):
        if t == 0:
            S[0] = complex(sum_delta_all)
        else:
            phases = np.exp(2j * np.pi * t * surv_ext / p_next)
            S[t] = np.sum(phases * delta_ext)

    # Correction oscillatoire = (1/p) * sum_{t=1}^{p-1} S(t)
    correction = np.sum(S[1:]) / p_next

    # Verification : bnd = (1/p) * sum_t S(t) = S(0)/p + correction
    bnd_fourier = np.sum(S) / p_next

    # ===== ANALYSE DES INTERACTIONS =====
    # Distance minimale entre positions retirees consecutives
    if len(removed_idx) > 1:
        distances = np.diff(removed_idx)
        min_dist = int(distances.min())
        n_adjacent = int(np.sum(distances <= 2))
    else:
        min_dist = M
        n_adjacent = 0

    results[k] = {
        'p': p_next, 'N_k': N_k, 'M': M, 'N_removed': N_removed,
        'n00_k': n00_k, 'n00_ext': n00_ext, 'n00_k1': n00_k1,
        'bnd_exact': bnd_exact,
        'sum_delta_removed': sum_delta_removed,
        'n00_predicted': n00_predicted,
        'sum_delta_all': sum_delta_all,
        'bnd_uniform_true': bnd_uniform_true,
        'S': S, 'correction': correction,
        'bnd_fourier': bnd_fourier,
        'min_dist': min_dist, 'n_adjacent': n_adjacent,
        'delta_ext': delta_ext,
        'surv_ext': surv_ext,
        'removed_idx': removed_idx,
    }


# ================================================================
# Tests
# ================================================================

class TestFourierDecomposition(unittest.TestCase):
    """Decomposition de Fourier du boundary CRT."""

    def test_01_extended_structure(self):
        """La sequence etendue a bien p*n00_k transitions 00."""
        print("\n" + "=" * 120)
        print("STRUCTURE DE LA SEQUENCE ETENDUE")
        print("=" * 120)
        print(f"{'k':>3} {'p':>4} {'N_k':>8} {'M=p*N_k':>10} {'n00_k':>8} "
              f"{'n00_ext':>10} {'p*n00_k':>10} {'match':>6}")
        print("-" * 120)

        for k in sorted(results):
            r = results[k]
            p_n00 = r['p'] * r['n00_k']
            match = "YES" if r['n00_ext'] == p_n00 else "NO"
            print(f"{k:>3} {r['p']:>4} {r['N_k']:>8} {r['M']:>10} {r['n00_k']:>8} "
                  f"{r['n00_ext']:>10} {p_n00:>10} {match:>6}")
            self.assertEqual(r['n00_ext'], p_n00,
                msg=f"k={k}: n00_ext={r['n00_ext']} != p*n00_k={p_n00}")

        print("=" * 120)

    def test_02_boundary_vs_delta_sum(self):
        """Verifier si n00_{k+1} = n00_ext + sum(delta aux retraits).

        Si les retraits interagissent, la prediction additive peut diverger.
        """
        print("\n" + "=" * 120)
        print("BOUNDARY : PREDICTION ADDITIVE vs EXACT")
        print("=" * 120)
        print(f"{'k':>3} {'n00_k1':>10} {'n00_predicted':>14} {'match':>6} "
              f"{'sum_delta':>10} {'min_dist':>9} {'n_adj':>6}")
        print("-" * 120)

        for k in sorted(results):
            r = results[k]
            match = "YES" if r['n00_predicted'] == r['n00_k1'] else "NO"
            print(f"{k:>3} {r['n00_k1']:>10} {r['n00_predicted']:>14} {match:>6} "
                  f"{r['sum_delta_removed']:>10} {r['min_dist']:>9} {r['n_adjacent']:>6}")

        print()
        print("  EXPLICATION :")
        print("  - min_dist = distance min entre retraits consecutifs dans la sequence")
        print("  - n_adj = nombre de paires de retraits a distance <= 2 (interactions)")
        print("  - Si n_adj > 0, la prediction additive peut etre inexacte")
        print("=" * 120)

    def test_03_fourier_reconstruction(self):
        """Verifier que bnd_fourier = sum_delta_removed (reconstruction)."""
        print("\n" + "=" * 120)
        print("RECONSTRUCTION DE FOURIER")
        print("=" * 120)
        print(f"{'k':>3} {'bnd_exact':>10} {'sum_delta':>10} {'S(0)/p':>10} "
              f"{'correction':>12} {'bnd_fourier':>12} {'err':>10}")
        print("-" * 120)

        for k in sorted(results):
            r = results[k]
            S0_p = r['S'][0].real / r['p']
            corr = r['correction'].real
            bf = r['bnd_fourier'].real
            err = abs(bf - r['sum_delta_removed'])
            print(f"{k:>3} {r['bnd_exact']:>10} {r['sum_delta_removed']:>10} "
                  f"{S0_p:>10.2f} {corr:>12.4f} {bf:>12.4f} {err:>10.2e}")

        print()
        print("  S(0)/p = prediction uniforme VRAIE = (1/p) * sum de TOUS les delta")
        print("  correction = (1/p) * sum_{t>=1} S(t) = partie oscillatoire")
        print("  bnd_fourier = S(0)/p + correction = boundary reconstruit par Fourier")
        print("=" * 120)

    def test_04_fourier_modes(self):
        """Analyse des modes de Fourier S(t) pour t = 0, ..., p-1."""
        print("\n" + "=" * 120)
        print("MODES DE FOURIER S(t)")
        print("=" * 120)

        for k in sorted(results):
            r = results[k]
            p = r['p']
            S = r['S']
            print(f"\n  k={k}, p={p}:")
            print(f"  {'t':>4} {'|S(t)|':>12} {'Re(S(t))':>12} {'Im(S(t))':>12} "
                  f"{'phase/pi':>10}")
            print("  " + "-" * 60)

            for t in range(p):
                amp = abs(S[t])
                re = S[t].real
                im = S[t].imag
                if amp > 1e-10:
                    phase = np.angle(S[t]) / np.pi
                else:
                    phase = 0
                mark = " <-- DC" if t == 0 else ""
                print(f"  {t:>4} {amp:>12.4f} {re:>12.4f} {im:>12.4f} "
                      f"{phase:>10.4f}{mark}")

        print("=" * 120)

    def test_05_periodicity_cancellation(self):
        """Test de l'annulation par periodicite.

        Si delta_j est periodique avec periode N_k (p copies du pattern),
        alors S(t) = [sum_m omega^{t*m*P_k}] * [sum local], et le premier
        facteur = 0 pour t != 0.

        L'ecart a cette annulation vient des INTERACTIONS entre retraits
        adjacents et des effets de bord du pattern cyclique.
        """
        print("\n" + "=" * 120)
        print("TEST D'ANNULATION PAR PERIODICITE")
        print("=" * 120)
        print(f"{'k':>3} {'p':>4} {'S(0)':>12} {'|sum S(t>0)|':>14} "
              f"{'ratio correction/S(0)':>22} {'n_adj':>6}")
        print("-" * 120)

        for k in sorted(results):
            r = results[k]
            S0 = abs(r['S'][0])
            sum_St = abs(np.sum(r['S'][1:]))
            ratio = sum_St / S0 if S0 > 1e-10 else 0
            print(f"{k:>3} {r['p']:>4} {r['S'][0].real:>12.2f} {sum_St:>14.4f} "
                  f"{ratio:>22.6f} {r['n_adjacent']:>6}")

        print()
        print("  Si le pattern est EXACTEMENT periodique ET les retraits n'interagissent")
        print("  pas, la correction oscillatoire = 0 (annulation geometrique).")
        print("  Un ratio > 0 indique soit des interactions, soit des effets de bord.")
        print("=" * 120)

    def test_06_sin2_structure(self):
        """Chercher un pattern sin^2 dans les modes de Fourier.

        En PT, sin^2(pi*t/p) apparait naturellement dans les produits
        d'Euler et les sommes de Ramanujan.
        """
        print("\n" + "=" * 120)
        print("RECHERCHE DE PATTERN sin^2 DANS LES MODES")
        print("=" * 120)

        for k in sorted(results):
            r = results[k]
            p = r['p']
            S = r['S']

            print(f"\n  k={k}, p={p}:")
            print(f"  {'t':>4} {'|S(t)|':>12} {'sin^2(pi*t/p)':>14} "
                  f"{'|S|/sin^2':>12} {'Re(S)/sin^2':>12}")
            print("  " + "-" * 60)

            for t in range(1, p):
                amp = abs(S[t])
                sin2 = np.sin(np.pi * t / p) ** 2
                ratio_a = amp / sin2 if sin2 > 1e-10 else 0
                ratio_r = S[t].real / sin2 if sin2 > 1e-10 else 0
                print(f"  {t:>4} {amp:>12.4f} {sin2:>14.6f} "
                      f"{ratio_a:>12.4f} {ratio_r:>12.4f}")

        print("=" * 120)

    def test_07_boundary_from_level_k_only(self):
        """Verifier si le boundary est entierement determine par le niveau k.

        Si le boundary = sum_{j0=1}^{N_k} Delta_{j0} (somme sur le pattern
        de niveau k), alors il ne depend PAS de la distribution des retraits
        parmi les p copies -- seulement de la structure du pattern.

        Cela eliminerait le besoin d'un theoreme d'equidistribution !
        """
        print("\n" + "=" * 120)
        print("BOUNDARY DEPUIS LE NIVEAU k SEUL")
        print("=" * 120)

        for k in sorted(results):
            r = results[k]
            # Delta sur le pattern de niveau k
            surv_k = sieve_survivors(PRIMES[:k])
            P_k = prod(PRIMES[:k])
            cls_k = gap_classes(surv_k, P_k)
            delta_k = compute_boundary_per_position(cls_k)
            sum_delta_k = int(delta_k.sum())

            # Prediction : n00_{k+1} = p * n00_k + sum_delta_k
            # (si chaque position du pattern contribue exactement une fois)
            n00_predicted_k = r['p'] * r['n00_k'] + sum_delta_k

            # Compare avec le reel
            match = "YES" if n00_predicted_k == r['n00_k1'] else "NO"
            ecart = r['n00_k1'] - n00_predicted_k

            print(f"  k={k}: sum_delta_k = {sum_delta_k:>8}, "
                  f"p*n00_k + sum_delta_k = {n00_predicted_k:>10}, "
                  f"n00_k1 = {r['n00_k1']:>10}, match={match}, ecart={ecart:>6}")

        print()
        print("  Si match=YES, le boundary est determine par le niveau k SEUL")
        print("  et il n'y a PAS de phenomene d'equidistribution (pas d'amplification).")
        print("  L'amplification de S15.6.301 serait un artefact du modele 3-gramme.")
        print("=" * 120)

    def test_08_true_vs_3gram_uniform(self):
        """Comparer la prediction uniforme VRAIE avec le modele 3-gramme."""
        print("\n" + "=" * 120)
        print("PREDICTION UNIFORME VRAIE vs MODELE 3-GRAMME")
        print("=" * 120)

        for k in sorted(results):
            r = results[k]

            # Prediction vraie : sum de tous les delta / p
            true_unif = r['sum_delta_all'] / r['p']

            # Modele 3-gramme (de S15.6.301)
            surv_k = sieve_survivors(PRIMES[:k])
            P_k = prod(PRIMES[:k])
            N_k = len(surv_k)
            cls_k = gap_classes(surv_k, P_k)
            c_from = cls_k
            c_to = np.roll(cls_k, -1)
            c_to2 = np.roll(cls_k, -2)

            gram3 = np.zeros((3, 3, 3), dtype=np.int64)
            for a in range(3):
                ma = (c_from == a)
                for b in range(3):
                    mab = ma & (c_to == b)
                    for c in range(3):
                        gram3[a, b, c] = int((mab & (c_to2 == c)).sum())

            gain_L = int(gram3[0,0,0] + gram3[0,1,2] + gram3[0,2,1])
            gain_R = int(gram3[0,0,0] + gram3[1,2,0] + gram3[2,1,0])
            loss_ab = int(gram3[0,0,:].sum())
            loss_bc = int(gram3[:,0,0].sum())
            bnd_3gram = gain_L + gain_R - loss_ab - 2 * loss_bc

            # Boundary exact
            bnd_exact = r['bnd_exact']

            # Boundary = 3*n00_k + sum_delta_k (theorique)
            surv_k2 = sieve_survivors(PRIMES[:k])
            cls_k2 = gap_classes(surv_k2, P_k)
            delta_k2 = compute_boundary_per_position(cls_k2)
            sum_delta_k2 = int(delta_k2.sum())
            bnd_from_delta = 3 * r['n00_k'] + sum_delta_k2

            print(f"  k={k}: bnd_exact={bnd_exact:>10}, "
                  f"3*n00+sum_delta={bnd_from_delta:>10}, "
                  f"bnd_3gram={bnd_3gram:>10}, "
                  f"true_unif={true_unif:>10.1f}")

        print("=" * 120)

    def test_09_sin2_boundary_formula(self):
        """Tester si le boundary a une formule en sin^2.

        Hypothese PT : boundary / N_k ∝ sin^2(pi/p) ou 1 - cos(2*pi/p)
        ou une combinaison de termes sin^2.
        """
        print("\n" + "=" * 120)
        print("FORMULE sin^2 POUR LE BOUNDARY")
        print("=" * 120)
        print(f"{'k':>3} {'p':>4} {'bnd/N_k':>10} {'sin2(pi/p)':>12} "
              f"{'1/p':>8} {'2/(p-1)':>9} {'bnd/N_k / sin2':>14}")
        print("-" * 120)

        for k in sorted(results):
            r = results[k]
            p = r['p']
            bnd_per_N = r['bnd_exact'] / r['N_k']
            sin2 = np.sin(np.pi / p) ** 2
            inv_p = 1.0 / p
            two_pm1 = 2.0 / (p - 1)
            ratio = bnd_per_N / sin2 if sin2 > 1e-10 else 0
            print(f"{k:>3} {p:>4} {bnd_per_N:>10.4f} {sin2:>12.6f} "
                  f"{inv_p:>8.4f} {two_pm1:>9.4f} {ratio:>14.4f}")

        print("=" * 120)

    def test_10_synthesis(self):
        """Synthese de l'exploration trigonometrique."""
        print("\n" + "=" * 120)
        print("SYNTHESE S15.6.302 : EXPLORATION TRIGONOMETRIQUE")
        print("=" * 120)

        # Verifier la prediction "niveau k seul"
        level_k_works = True
        for k in sorted(results):
            r = results[k]
            surv_k = sieve_survivors(PRIMES[:k])
            P_k = prod(PRIMES[:k])
            cls_k = gap_classes(surv_k, P_k)
            delta_k = compute_boundary_per_position(cls_k)
            sum_delta_k = int(delta_k.sum())
            predicted = r['p'] * r['n00_k'] + sum_delta_k
            if predicted != r['n00_k1']:
                level_k_works = False

        print()
        if level_k_works:
            print("  DECOUVERTE MAJEURE :")
            print("  Le boundary est ENTIEREMENT determine par le pattern de niveau k.")
            print("  n00_{k+1} = p * n00_k + sum_{j=1}^{N_k} Delta_j(k)")
            print("  ou Delta_j(k) ne depend que du 4-gramme local au niveau k.")
            print()
            print("  CONSEQUENCE : il n'y a PAS besoin d'un theoreme d'equidistribution !")
            print("  L'amplification de S15.6.301 etait un artefact du modele 3-gramme.")
            print("  Le vrai probleme est : prouver que sum Delta_j > -3*n00_k + seuil.")
        else:
            print("  Le boundary n'est PAS determine par le niveau k seul.")
            print("  Il y a des interactions entre retraits ou des effets non-locaux.")
            print("  Le theoreme d'equidistribution est toujours necessaire.")

        # Chercher la meilleure formule sin^2
        print()
        print("  Analyse sin^2 :")
        for k in sorted(results):
            r = results[k]
            p = r['p']
            correction_ratio = abs(np.sum(r['S'][1:])) / abs(r['S'][0]) if abs(r['S'][0]) > 1e-10 else 0
            print(f"  k={k}: |correction|/|S(0)| = {correction_ratio:.6f}")

        print()
        print("=" * 120)


if __name__ == '__main__':
    unittest.main(verbosity=2)
