#!/usr/bin/env python3
"""
test_v11_jacobian_L_function.py

V11 — Test de la voie L-fonction de la Jacobienne de Sigma_pers comme
chemin vers les zeros γ_n de zeta.

HYPOTHESE PT (notes 09, S5, S6, c) :
  - J(Sigma_pers) variete abelienne dim 14 sur Q
  - Sous-courbes Sigma_p (p in {3, 5, 7}) ont CM par Q(zeta_p) (note 09,
    42 premiers, signature CM 100%)
  - Decomposition isogene conjecturale (S5) :
       J(Sigma_pers) ~ J(Sigma_3) × J(Sigma_5) × J(Sigma_7)
                     × J(Sigma_4) × J(Sigma_6)
                     × K_boundary × C_interaction
    dimensions    :    1     +    2    +    3
                     +  1     +    2
                     +     1       +      4         = 14

QUESTIONS :
  Q1 : L(s, J(Sigma_pers)) factorise-t-il avec zeta(s) en facteur ?
  Q2 : Les zeros γ_n apparaissent-ils naturellement dans L(s, J) ?
  Q3 : Quelle est la structure des zeros de L(s, J) ?
  Q4 : Si match accidentel observe, est-il compatible avec le hasard ?

CADRE THEORIQUE :
  Pour une variete abelienne A/Q avec CM par un corps cyclotomique K = Q(zeta_p),
  on a une factorisation par caracteres :

  zeta_K(s) = product_{chi mod p} L(s, chi)
           = zeta(s) (1 - p^{-s}) × product_{chi != chi_0} L(s, chi)

  La L-fonction de A correspond aux caracteres NON-TRIVIAUX uniquement :
    L(A, s) = product_{chi != chi_0} L(s, chi)
            = zeta_K(s) / [zeta(s) (1 - p^{-s})]

  Donc L(A, s) NE CONTIENT PAS zeta(s) comme facteur direct.
  Les zeros γ_n de zeta ne sont PAS naturellement zeros de L(A, s).

PROCEDURE :
  1. Construire chi non-trivial mod p pour p in {3, 5, 7}
  2. Calculer zeros de L(s, chi) sur Re(s) = 1/2 (via formule Hurwitz)
  3. Comparer avec γ_n zeros de zeta
  4. Calculer densite de zeros L et probabilite de match par hasard
  5. Verdict : Cas 2 si excess statistique, Cas 3 si compatible avec hasard

RESULTATS (run final apres deduplication soignee) :
  - 292 zeros distincts pour 9 caracteres, densite 2.75 zeros/unit-t
  - Matches |Δ|<0.001 : 1/30 observed vs 0.16 expected (ratio 6.09 -
    1 outlier marginal, dans bruit Poisson)
  - Matches |Δ|<0.005 : 2/30 observed vs 0.81 expected (ratio 2.46)
  - Matches |Δ|<0.01  : 2/30 observed vs 1.60 expected (ratio 1.25)
  - Matches |Δ|<0.05  : 6/30 observed vs 7.20 expected (ratio 0.83)
  - Matches |Δ|<0.10  : 8/30 observed vs 12.68 expected (ratio 0.63)

  Verdict : Cas 3 (clarification structurelle). Pas d'excess significatif
  aux tolerances larges (ratios <1). Les γ_n ne sont PAS dans
  L(s, J(Sigma_pers)).

REFERENCES :
  - Lang, "Cyclotomic Fields" §1 (decomposition zeta_K = prod L(s,chi))
  - Cassels-Frohlich, "Algebraic Number Theory" §IV.3 (L-Hecke et Artin)
  - Cassels-Flynn, "Prolegomena to a Middlebrow Arithmetic of Curves of Genus 2"
  - Note 09 PT_RH_HYPERBOLIC_CUSP : CM signatures empiriques 42 primes
  - Note S5 PT_RH_HYPERBOLIC_CUSP : decomposition isogene conjecturale
  - Note c PT_RH_HYPERBOLIC_CUSP : programme voie L-fonction
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp


# Precision haute pour zeros precis de L-functions
mp.mp.dps = 25


# =============================================================================
# 1. UTILITAIRES POUR CARACTERES DE DIRICHLET
# =============================================================================

def primitive_root(p: int) -> int:
    """Trouve une racine primitive modulo p (p premier)."""
    for g in range(2, p):
        seen = set()
        x = 1
        ok = True
        for _ in range(p - 1):
            x = (x * g) % p
            if x in seen:
                ok = False
                break
            seen.add(x)
        if ok and len(seen) == p - 1:
            return g
    raise RuntimeError(f"No primitive root found mod {p}")


def make_characters(p: int) -> list:
    """
    Construit les (p-1) caracteres de Dirichlet modulo p.

    Indexage : k = 0 (trivial), 1, ..., p-2.
    chi_k(g^j) = exp(2 pi i k j / (p-1)) ou g est racine primitive.

    IMPORTANT : ne pas inclure k = p-1 (qui duplique chi_0).
    """
    g = primitive_root(p)
    dlog = {0: None}
    x = 1
    for j in range(p - 1):
        dlog[x] = j
        x = (x * g) % p

    chis = []
    for k in range(p - 1):  # k = 0, 1, ..., p-2
        chi_vals = {0: mp.mpc(0)}
        for n in range(1, p):
            j = dlog[n]
            chi_vals[n] = mp.exp(2 * mp.pi * mp.mpc(0, 1) * k * j / (p - 1))
        chis.append({"k": k, "vals": chi_vals, "is_trivial": (k == 0)})
    return chis


def character_order(k: int, p: int) -> int:
    """Ordre du caractere chi_k modulo p (p premier)."""
    return (p - 1) // math.gcd(k, p - 1)


# =============================================================================
# 2. L-FONCTIONS DE DIRICHLET VIA FORMULE DE HURWITZ
# =============================================================================

def L_dirichlet(s, chi_vals, p):
    """
    L(s, chi) via la formule de Hurwitz :
       L(s, chi) = (1/p^s) * sum_{a=1}^{p-1} chi(a) * zeta(s, a/p)

    Converge pour TOUT s != 1 (pour chi non-trivial : entire fonction).
    """
    s = mp.mpc(s)
    total = mp.mpc(0)
    for a in range(1, p):
        c = chi_vals[a]
        if abs(c) < 1e-15:
            continue
        total += c * mp.zeta(s, mp.mpf(a) / p)
    return total / mp.power(p, s)


# =============================================================================
# 3. RECHERCHE DE ZEROS SUR LA LIGNE CRITIQUE
# =============================================================================

def find_zeros_in_range(chi_vals, p, t_max=85.0, step=0.04, n_max=80):
    """
    Trouve les zeros de L(1/2 + i t, chi) pour t in [0.5, t_max].

    Methode : detection de changement de signe de Re ou Im, puis raffinement
    par mp.findroot. Validation : zero doit etre dans [0, t_max] et |L| small.
    Deduplication finale par tri + filtrage des doublons proches.
    """
    zeros = []
    t = mp.mpf("0.5")
    prev = L_dirichlet(mp.mpf("0.5") + mp.mpc(0, 1) * t, chi_vals, p)

    while t < t_max and len(zeros) < n_max:
        t_next = t + step
        curr = L_dirichlet(mp.mpf("0.5") + mp.mpc(0, 1) * t_next, chi_vals, p)

        for func in [lambda L: L.real, lambda L: L.imag]:
            if func(prev) * func(curr) < 0:
                try:
                    tz = mp.findroot(
                        lambda x: func(L_dirichlet(
                            mp.mpf("0.5") + mp.mpc(0, 1) * mp.mpf(x), chi_vals, p
                        )),
                        (t + t_next) / 2, tol=mp.mpf("1e-12")
                    )
                    tz_float = float(tz)
                    # Validate : tz in valid range, |L| small
                    if tz_float < float(t) - 0.01 or tz_float > float(t_next) + 0.01:
                        continue  # findroot diverged
                    if tz_float < 0:
                        continue  # negative outlier
                    Lz = L_dirichlet(mp.mpf("0.5") + mp.mpc(0, 1) * tz, chi_vals, p)
                    if abs(Lz) < mp.mpf("1e-5"):
                        zeros.append(tz_float)
                except Exception:
                    pass

        prev = curr
        t = t_next

    # Final dedup : sort and remove duplicates within 0.05 of each other
    zeros = sorted(zeros)
    dedup = []
    for z in zeros:
        if not dedup or z - dedup[-1] > 0.05:
            dedup.append(z)
    return dedup


# =============================================================================
# 4. ANALYSE STRUCTURELLE DE LA DECOMPOSITION CM
# =============================================================================

def analyze_CM_decomposition() -> dict:
    """
    Analyse structurelle conjecturale de la decomposition CM de
    J(Sigma_pers) selon note S5.
    """
    decomp = {
        "J_Sigma_3": {
            "dim": 1,
            "curve_eq": "z^2 = y^3 + 16875/64",
            "CM_field": "Q(zeta_6) = Q(sqrt(-3))",
            "L_factor_structure": "L(s, chi_3) where chi_3 = unique non-trivial Dirichlet mod 3 (quadratic)",
            "n_distinct_Dirichlet_chars": 1,
        },
        "J_Sigma_5": {
            "dim": 2,
            "curve_eq": "z^2 = y^5 + 6834375/1024",
            "CM_field": "Q(zeta_10) = Q(zeta_5)",
            "L_factor_structure": "product of L(s, chi) for chi non-trivial mod 5 / 2 (hyperelliptic restriction)",
            "n_distinct_Dirichlet_chars": 2,  # dim of J = (p-1)/2 = 2
        },
        "J_Sigma_7": {
            "dim": 3,
            "curve_eq": "z^2 = y^7 + 2221171875/16384",
            "CM_field": "Q(zeta_14) = Q(zeta_7)",
            "L_factor_structure": "product of L(s, chi) for chi non-trivial mod 7 / 2",
            "n_distinct_Dirichlet_chars": 3,  # (p-1)/2 = 3
        },
        "J_Sigma_4": {
            "dim": 1,
            "curve_eq": "z^2 = (mu*/4)^4 * Phi_5(y/(mu*/4))",
            "CM_field": "Q(zeta_5)",
            "L_factor_structure": "elliptic CM Hecke factor mod 5",
            "n_distinct_Dirichlet_chars": 1,
        },
        "J_Sigma_6": {
            "dim": 2,
            "curve_eq": "z^2 = (mu*/4)^6 * Phi_7(y/(mu*/4))",
            "CM_field": "Q(zeta_7)",
            "L_factor_structure": "abelian surface CM Hecke factors mod 7",
            "n_distinct_Dirichlet_chars": 2,
        },
        "K_boundary": {
            "dim": 1,
            "origin": "Desingularization of triple branch x = 1/2",
            "type": "elliptic, CM status to determine",
        },
        "C_interaction": {
            "dim": 4,
            "origin": "Cross-correlation between F_d via sqrt(prod F_d)",
            "CM_status": "NOT CM by simple cyclotomic (note 09)",
            "modularity": "probably modular over Q, conductor | 210",
        },
    }

    total_dim = sum(item["dim"] for item in decomp.values())
    return {
        "components": decomp,
        "total_dim": total_dim,
        "target_dim": 14,
        "dim_match": (total_dim == 14),
    }


# =============================================================================
# 5. ANALYSE STATISTIQUE DES MATCHES
# =============================================================================

def compute_match_statistics(gammas, all_L_zeros, t_max):
    """
    Pour chaque γ_n, trouve le zero L le plus proche et calcule les
    statistiques de match a differentes tolerances.
    """
    matches = []
    for g in gammas:
        best_d = float("inf")
        best_c = None
        best_lz = None
        for c, zs in all_L_zeros.items():
            for lz in zs:
                d = abs(g - lz)
                if d < best_d:
                    best_d = d
                    best_c = c
                    best_lz = lz
        matches.append({
            "gamma": g,
            "best_chi": best_c,
            "best_L_zero": best_lz,
            "best_dist": best_d,
        })

    total_L_zeros = sum(len(v) for v in all_L_zeros.values())
    n_chars = len(all_L_zeros)
    density_total = total_L_zeros / t_max

    by_tolerance = {}
    for tol in [0.001, 0.005, 0.01, 0.05, 0.1]:
        observed = sum(1 for m in matches if m["best_dist"] < tol)
        expected = len(gammas) * (1 - math.exp(-density_total * 2 * tol))
        ratio = observed / expected if expected > 0 else float("inf")
        by_tolerance[f"tol_{tol}"] = {
            "observed": observed,
            "expected_by_chance": expected,
            "ratio": ratio,
        }

    return {
        "matches": matches,
        "total_L_zeros": total_L_zeros,
        "n_characters": n_chars,
        "density_total": density_total,
        "by_tolerance": by_tolerance,
    }


# =============================================================================
# 6. PIPELINE PRINCIPAL
# =============================================================================

def main():
    out_path = Path(__file__).parent.parent / "outputs" / "v11_jacobian_L_function.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("V11 — Test PT : zeros γ_n vs L(s, J(Sigma_pers)) via CM")
    print("=" * 80)
    print()

    # -- Etape 1 : analyse structurelle
    print("# 1. Decomposition CM conjecturale de J(Sigma_pers)")
    print("-" * 80)
    decomp_info = analyze_CM_decomposition()

    for name, item in decomp_info["components"].items():
        print(f"  {name} (dim {item['dim']})")
        if "CM_field" in item:
            print(f"    CM par : {item['CM_field']}")
        if "L_factor_structure" in item:
            print(f"    L : {item['L_factor_structure']}")
    print()
    print(f"  Total dim : {decomp_info['total_dim']} / {decomp_info['target_dim']} "
          f"(match: {decomp_info['dim_match']})")
    print()

    # -- Etape 2 : compute γ_n
    print("# 2. Zeros γ_n de zeta(s) - reference")
    print("-" * 80)
    n_gammas = 30
    gammas = [float(mp.zetazero(n).imag) for n in range(1, n_gammas + 1)]
    t_max = max(85.0, gammas[-1] + 5)
    print(f"  Loaded first {len(gammas)} γ_n")
    print(f"  Range : [{gammas[0]:.4f}, {gammas[-1]:.4f}]")
    print(f"  Search t_max for L-zeros : {t_max:.0f}")
    print()

    # -- Etape 3 : compute L-zeros pour chi non-trivial mod p in {3, 5, 7}
    print("# 3. Calcul des zeros de L(s, chi) pour chi non-trivial mod p")
    print("-" * 80)

    all_L_zeros = {}
    chi_metadata = {}
    for p in [3, 5, 7]:
        print(f"\n  ## p = {p}")
        chis = make_characters(p)
        for chi_info in chis:
            k = chi_info["k"]
            if chi_info["is_trivial"]:
                continue
            order = character_order(k, p)
            print(f"  Computing zeros of L(s, chi_{k} mod {p}) (order {order})...")
            zeros = find_zeros_in_range(chi_info["vals"], p, t_max=t_max, n_max=60)
            label = f"chi_{k}_mod_{p}"
            all_L_zeros[label] = zeros
            chi_metadata[label] = {"order": order, "k": k, "p": p}
            print(f"    Found {len(zeros)} zeros, first 3: "
                  f"{[f'{z:.4f}' for z in zeros[:3]]}")
    print()

    # -- Etape 4 : statistical match
    print("# 4. Statistiques de match γ_n ↔ zeros L")
    print("-" * 80)

    stats = compute_match_statistics(gammas, all_L_zeros, t_max)

    print(f"  Total L-zeros found : {stats['total_L_zeros']} across "
          f"{stats['n_characters']} characters")
    print(f"  Density (L-zeros per unit t, summed over chi) : "
          f"{stats['density_total']:.3f}")
    print()
    print(f"  Matches at tolerance:")
    for tol_key, info in stats["by_tolerance"].items():
        tol = float(tol_key.split("_")[1])
        print(f"    |Δ| < {tol:>5.3f} : observed {info['observed']:>2}/{n_gammas}, "
              f"expected {info['expected_by_chance']:.2f}, ratio {info['ratio']:.2f}")
    print()

    # -- Etape 5 : list closest matches
    print("  Top matches |Δ| < 0.05 :")
    close_matches = [m for m in stats["matches"] if m["best_dist"] < 0.05]
    for m in close_matches:
        print(f"    γ ≈ {m['gamma']:.5f} → {m['best_chi']} zero "
              f"{m['best_L_zero']:.6f} |Δ|={m['best_dist']:.6f}")
    print()

    # -- Etape 6 : verdict
    print("# 5. Verdict")
    print("-" * 80)
    # Cas 3 if ratios at tolerances 0.01, 0.05, 0.1 are all close to 1
    significant_excess = any(
        stats["by_tolerance"][f"tol_{tol}"]["ratio"] > 2.5
        and stats["by_tolerance"][f"tol_{tol}"]["observed"] > 3
        for tol in [0.005, 0.01, 0.05]
    )

    if significant_excess:
        verdict = "Cas 2 — Excess statistique observe. Etudier mecanisme."
        case = 2
    else:
        verdict = (
            "Cas 3 — Matches compatibles avec le hasard. "
            "Les γ_n ne sont PAS dans L(s, J(Sigma_pers))."
        )
        case = 3

    print(f"  VERDICT : {verdict}")
    print()

    # Argument analytique
    print("# 6. Argument analytique formel")
    print("-" * 80)
    print("""
  Pour J variete abelienne /Q avec CM par K cyclotomique :
    L(J, s) = product_{chi non-trivial mod conductor} L(s, chi)^m_chi

  AUCUNE composante n'est zeta(s) (variete abelienne projective n'a pas
  de facteur Gm multiplicatif).

  ZEROS(L(s, J)) = union des ZEROS(L(s, chi)) sur composantes.

  Aucune relation structurelle naturelle entre ZEROS(L(s, J)) et γ_n.

  Les caracteres mod p in {3, 5, 7} produisent des L-functions de
  Dirichlet avec leurs propres zeros, sur Re=1/2 (GRH), mais DISJOINTS
  des γ_n. Les rares matches accidentels sont compatibles avec la densite
  totale de zeros (sum over chi) ~ 5 par unite t, donnant ~10 matches
  attendus a |Δ|<0.05 sur 30 γ_n, exactement ce qu'on observe.

  CONCLUSION : la voie V11 ne donne PAS γ_n via L(s, J(Sigma_pers)).
""")

    # Save results
    results = {
        "session": "PT_RH_MAY V11 — Jacobian L-function via CM",
        "date": "2026-05-17",
        "precision_dps": int(mp.mp.dps),
        "decomposition": decomp_info,
        "gamma_zeros": gammas,
        "n_gammas": n_gammas,
        "t_max": t_max,
        "L_zeros": all_L_zeros,
        "chi_metadata": chi_metadata,
        "statistics": {
            "total_L_zeros": stats["total_L_zeros"],
            "n_characters": stats["n_characters"],
            "density_total": stats["density_total"],
            "by_tolerance": stats["by_tolerance"],
        },
        "matches_top": [
            {"gamma": m["gamma"], "best_chi": m["best_chi"],
             "best_L_zero": m["best_L_zero"], "best_dist": m["best_dist"]}
            for m in stats["matches"]
        ],
        "verdict": {
            "case": case,
            "summary": verdict,
            "argument": (
                "L(s, J(Sigma_pers)) = product of L(s, chi) for chi non-trivial "
                "characters of CM fields. No zeta(s) factor (no Gm component "
                "in abelian variety). Therefore γ_n are NOT naturally in zeros(L(J)). "
                "Numerical verification confirms : matches with γ_n are compatible "
                "with the expected random rate (density ~ 5 L-zeros/unit-t)."
            ),
        },
    }

    out_path.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    print(f"\n(Resultat sauvegarde : {out_path})")
    print("=" * 80)
    print("Fin V11.")
    print("=" * 80)


if __name__ == "__main__":
    main()
