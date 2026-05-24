"""
W7-bis — Validation étendue du Théorème W7-1 σ_crit^(k) = √(π(k+1)).

CONTEXTE
========
Le test W7 initial (p_max = 5×10⁵) a validé σ_crit^(1) = √(2π) à 0.73 %
mais n'a pas pu valider σ_crit^(2) = √(3π) à cause de la truncation de
turn 2 (0.2 % de couverture).

W7-bis sieve jusqu'à p_max = 1.535 × 10⁸ ≈ e^{6π} pour couvrir turn 2
INTÉGRALEMENT (~8.4 millions de primes).

OBJECTIF
========
Tester :
    σ_crit^(2) prédit = √(3π) ≈ 3.0700
    σ_crit^(2) observé = ?

Si validation à <1 %, le Théorème W7-1 est solidement établi.

OPTIMISATIONS
=============
- Sieve numpy optimisé pair/impair pour p_max = 1.5e8 (~30 sec, ~150MB)
- float64 plutôt que mpmath dps=40 : suffit pour σ ∈ [1.5, 5]
- Vectorisation numpy de toute la somme
- j=1 dominant pour primes ≥ 11 ; j=1..4 retenus pour cohérence avec W7

PROTOCOLE
=========
1. Sieve primes ≤ 1.535×10⁸ (couverture totale turn 2)
2. Binning par turn k = floor(ln p / 2π)
3. Pour σ ∈ {1.5, 2.0, √(2π), 3.0, √(3π), √(4π), 4.0, 5.0, √(5π)} :
   - Calculer T_G^turn_k(σ) pour k = 0, 1, 2 (et 3 partiel)
   - Comparer aux prédictions exactes √(π(k+1))
4. Reporter écart σ_crit observé vs prédit pour k=1, k=2, k=3
"""

import numpy as np
import json
import os
import time
from datetime import datetime

# -----------------------------------------------------------------------------
# Sieve numpy optimisé pair/impair
# -----------------------------------------------------------------------------

def sieve_numpy_fast(n):
    """Sieve d'Ératosthène vectorisé pair/impair. Retourne array de primes."""
    if n < 2:
        return np.array([], dtype=np.int64)
    # Sieve uniquement les indices impairs ≥ 3
    # sieve[i] représente le nombre 2i+1 (pour i ≥ 1)
    # Plus simple : sieve standard, on optimise par marquage 2-step
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    # Marquer tous les pairs ≥ 4 comme non-prime
    sieve[4::2] = False
    # Boucle impairs jusqu'à √n
    sqrt_n = int(np.sqrt(n)) + 1
    for i in range(3, sqrt_n + 1, 2):
        if sieve[i]:
            # Marquer multiples impairs de i à partir de i²
            sieve[i*i::2*i] = False
    return np.flatnonzero(sieve).astype(np.int64)

# -----------------------------------------------------------------------------
# Main parameters
# -----------------------------------------------------------------------------

# p_max = ceil(e^{6π}) = 153552935 ≈ 1.535e8
P_MAX = 153552935  # e^{6π} exact (couvre turn 2 entièrement)
TWO_PI = 2 * np.pi

print(f"Sieving primes ≤ {P_MAX:,} (≈ e^{{6π}})...")
t0 = time.time()
ALL_PRIMES = sieve_numpy_fast(P_MAX)
print(f"  Trouvé {len(ALL_PRIMES):,} primes en {time.time() - t0:.1f}s.")
print(f"  Premier : {ALL_PRIMES[0]}, dernier : {ALL_PRIMES[-1]:,}")

# Précompute log_p pour vectorisation
print("Précomputation log p...")
LOG_P = np.log(ALL_PRIMES.astype(np.float64))
SQRT_P = np.sqrt(ALL_PRIMES.astype(np.float64))

# Binning par turn
print("Binning par tour de spirale...")
TURN_INDEX = np.floor(LOG_P / TWO_PI).astype(np.int32)
TURNS_PRESENT = np.unique(TURN_INDEX)
print(f"  Tours présents : {TURNS_PRESENT.tolist()}")
for k in TURNS_PRESENT:
    mask = (TURN_INDEX == k)
    count = mask.sum()
    p_min_k = ALL_PRIMES[mask].min()
    p_max_k = ALL_PRIMES[mask].max()
    interval_low = np.exp(TWO_PI * k)
    interval_high = np.exp(TWO_PI * (k + 1))
    print(f"  Turn {k} : [e^{{{2*k}π}}, e^{{{2*(k+1)}π}}) = "
          f"[{interval_low:,.0f}, {interval_high:,.0f}) → "
          f"{count:,} primes (range {p_min_k:,} - {p_max_k:,})")

# Tour 2 totalement couverte ?
turn_2_count = (TURN_INDEX == 2).sum()
expected_turn_2 = (P_MAX / np.log(P_MAX)) - (np.exp(4 * np.pi) / (4 * np.pi))
print(f"\n  Tour 2 : {turn_2_count:,} primes observés "
      f"(PNT estime ~{expected_turn_2:,.0f}, OK ratio {turn_2_count/expected_turn_2:.2f})")

# -----------------------------------------------------------------------------
# T_G^turn_k vectorisé
# -----------------------------------------------------------------------------

def T_G_turn_k_vec(mask, sigma, j_max=4):
    """T_G^turn_k = 2 Σ_{p in turn k, j=1..j_max} (log p / p^{j/2}) g(j log p, σ).
    Tout en vectorisé numpy."""
    log_p = LOG_P[mask]
    p = ALL_PRIMES[mask].astype(np.float64)
    sqrt_p = SQRT_P[mask]

    norm = 1.0 / (2 * sigma * np.sqrt(np.pi))
    total = 0.0
    for j in range(1, j_max + 1):
        # g(j log p, σ) = norm * exp(-(j log p)² / (4σ²))
        jlogp = j * log_p
        g_j = norm * np.exp(-jlogp ** 2 / (4 * sigma ** 2))
        # term per prime : 2 * log_p / p^{j/2} * g_j
        if j == 1:
            p_pow = sqrt_p
        elif j == 2:
            p_pow = p
        elif j == 3:
            p_pow = p * sqrt_p
        else:
            p_pow = p ** (j / 2)
        contributions = 2 * log_p / p_pow * g_j
        total += contributions.sum()
    return float(total)

# -----------------------------------------------------------------------------
# Asymptotic predictions
# -----------------------------------------------------------------------------

def predicted_sigma_crit_k(k):
    """σ_crit^(k) = √(π(k+1))."""
    return np.sqrt(np.pi * (k + 1))

def predicted_ratio_turn_k(k, sigma):
    """Ratio asymptotique e^{πk} · exp(-π² k(k+1)/σ²) (approximation midpoint).
    Note : la prédiction EXACTE est l'intégrale ∫_{2πk}^{2π(k+1)} exp(-(u-σ²)²/(4σ²)) du
    rapportée à l'intégrale sur turn 0. Pour validation, on compare aux ratios mesurés."""
    return np.exp(np.pi * k) * np.exp(-np.pi**2 * k * (k+1) / sigma**2)

def predicted_ratio_exact_continuum(k, sigma):
    """Ratio EXACT continuum I_k(σ) / I_0(σ) avec intégrales gaussiennes.
    I_k(σ) = ∫_{2πk}^{2π(k+1)} exp(-(u-σ²)²/(4σ²)) du
    Utiliser fonction d'erreur via erf."""
    from scipy.special import erf
    sigma_u = sigma * np.sqrt(2)
    # I_k = (σ_u √π / 2) [erf((2π(k+1) - σ²)/σ_u) - erf((2πk - σ²)/σ_u)]
    def I_k(k_):
        a = (TWO_PI * k_ - sigma**2) / sigma_u
        b = (TWO_PI * (k_ + 1) - sigma**2) / sigma_u
        return erf(b) - erf(a)  # facteur (σ_u √π / 2) se simplifie dans le ratio
    return I_k(k) / I_k(0)

# -----------------------------------------------------------------------------
# Main scan
# -----------------------------------------------------------------------------

def main():
    SIGMAS = [
        1.5, 2.0, 2.4, 2.5, 2.5066,  # autour σ_crit^(1) = √(2π) = 2.5066
        2.8, 3.0, 3.07, 3.0700,       # autour σ_crit^(2) = √(3π) = 3.0700
        3.3, 3.5, 3.5449,             # autour σ_crit^(3) = √(4π) = 3.5449
        4.0, 5.0,
    ]
    SIGMAS = sorted(set(SIGMAS))

    TURNS = [int(k) for k in TURNS_PRESENT]
    masks = {k: (TURN_INDEX == k) for k in TURNS}

    # Tour 0 fully covered, 1 fully covered, 2 fully covered (with p_max = e^{6π})
    # Tour 3+ : empty (p_max = e^{6π})
    print()
    print("=" * 120)
    print("W7-bis — Validation étendue du Théorème W7-1 σ_crit^(k) = √(π(k+1))")
    print("=" * 120)
    print(f"p_max = {P_MAX:,}, {len(ALL_PRIMES):,} primes, "
          f"turns {TURNS} (turn 2 totalement couvert)")
    print()

    # Header
    print(f"{'σ':>9} " + " ".join(f"{'T_G^t'+str(k):>14}" for k in TURNS)
          + " " + " ".join(f"{'r_'+str(k)+'/0_obs':>12}" for k in TURNS if k > 0)
          + " " + " ".join(f"{'r_'+str(k)+'/0_cont':>12}" for k in TURNS if k > 0))
    print("-" * 120)

    results = {
        "metadata": {
            "date": datetime.now().isoformat(),
            "p_max": P_MAX,
            "n_primes": len(ALL_PRIMES),
            "turns_present": TURNS,
            "primes_per_turn": {str(k): int(masks[k].sum()) for k in TURNS},
            "j_max_powers": 4,
            "precision": "float64",
        },
        "predicted_sigma_crit": {
            str(k): float(predicted_sigma_crit_k(k)) for k in [1, 2, 3, 4]
        },
        "by_sigma": {}
    }

    for sigma in SIGMAS:
        T_turns = {k: T_G_turn_k_vec(masks[k], sigma) for k in TURNS}
        T_0 = T_turns[0]
        ratios_obs = {k: T_turns[k] / T_0 for k in TURNS if k > 0}
        ratios_cont = {k: float(predicted_ratio_exact_continuum(k, sigma))
                       for k in TURNS if k > 0}

        row = f"{sigma:>9.4f}"
        for k in TURNS:
            row += f" {T_turns[k]:>14.4e}"
        for k in TURNS:
            if k > 0:
                row += f" {ratios_obs[k]:>12.4f}"
        for k in TURNS:
            if k > 0:
                row += f" {ratios_cont[k]:>12.4f}"
        print(row)

        results["by_sigma"][str(sigma)] = {
            "sigma": float(sigma),
            "T_G_turns": {str(k): T_turns[k] for k in TURNS},
            "ratios_observed": {str(k): ratios_obs[k] for k in TURNS if k > 0},
            "ratios_continuum_predicted": {str(k): ratios_cont[k] for k in TURNS if k > 0},
        }

    # ------------------------------------------------------------------------
    # Find observed σ_crit for each turn (where ratio crosses 1) via bissection
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("SEUILS CRITIQUES — Validation σ_crit^(k) = √(π(k+1))")
    print("=" * 120)
    print(f"{'k':>5} {'σ_pred':>12} {'σ_obs':>12} {'écart':>10}")

    for k in TURNS:
        if k == 0:
            continue
        sigma_pred = float(predicted_sigma_crit_k(k))

        # Interpolation linéaire entre les σ pour trouver crossing
        sigmas_arr = sorted(SIGMAS)
        sigma_obs = None
        for i in range(len(sigmas_arr) - 1):
            s1, s2 = sigmas_arr[i], sigmas_arr[i+1]
            r1 = results["by_sigma"][str(s1)]["ratios_observed"][str(k)]
            r2 = results["by_sigma"][str(s2)]["ratios_observed"][str(k)]
            if (r1 - 1) * (r2 - 1) < 0:
                sigma_obs = s1 + (s2 - s1) * (1 - r1) / (r2 - r1)
                break

        ecart_str = "non traversé"
        if sigma_obs is not None:
            ecart = (sigma_obs - sigma_pred) / sigma_pred * 100
            ecart_str = f"{ecart:+.3f}%"
            results[f"sigma_crit_observed_k={k}"] = float(sigma_obs)
            results[f"sigma_crit_predicted_k={k}"] = sigma_pred
            results[f"sigma_crit_rel_error_k={k}"] = float((sigma_obs - sigma_pred) / sigma_pred)

        sigma_obs_str = f"{sigma_obs:.4f}" if sigma_obs is not None else "—"
        print(f"{k:>5} {sigma_pred:>12.4f} {sigma_obs_str:>12} {ecart_str:>10}")

    # ------------------------------------------------------------------------
    # Détail continuum exact vs observé pour chaque k, σ
    # ------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("RATIO observé / prédit continuum exact (Σ_pers PNT-density)")
    print("=" * 120)
    print(f"{'σ':>9} " + " ".join(f"{'k='+str(k):>14}" for k in TURNS if k > 0))
    print("-" * 120)
    for sigma in SIGMAS:
        row = f"{sigma:>9.4f}"
        for k in TURNS:
            if k > 0:
                r_obs = results["by_sigma"][str(sigma)]["ratios_observed"][str(k)]
                r_cont = results["by_sigma"][str(sigma)]["ratios_continuum_predicted"][str(k)]
                ratio = r_obs / r_cont if r_cont != 0 else float("inf")
                row += f" {ratio:>14.4f}"
        print(row)

    print()
    print("Si Théorème W7-1 strictement validé : ratio → 1 dans la limite PNT.")
    print("L'écart résiduel mesure la déviation discret-vs-continuum.")

    # Save
    out_dir = os.path.join(os.path.dirname(__file__), "..", "outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "w7bis_extended_sieve_results.json")

    # Convert numpy ints in metadata to native int for JSON
    def make_json_safe(obj):
        if isinstance(obj, dict):
            return {k: make_json_safe(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [make_json_safe(x) for x in obj]
        elif isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    with open(out_path, "w") as f:
        json.dump(make_json_safe(results), f, indent=2)
    print()
    print(f"Résultats sauvegardés : {os.path.relpath(out_path)}")

    return results


if __name__ == "__main__":
    main()
