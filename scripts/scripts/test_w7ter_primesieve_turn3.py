"""
W7-ter — Validation du Théorème W7-1 pour k=3 via primesieve.

OBJECTIF
========
Étendre la validation du Théorème W7-1 (σ_crit^(k) = √(π(k+1))) au cas k=3 :
    σ_crit^(3) prédit = √(4π) ≈ 3.5449

Turn 3 = [e^{6π}, e^{8π}) ≈ [1.535×10⁸, 7.896×10¹⁰).
Contient ~3.17 × 10⁹ premiers (PNT).

STRATÉGIE
=========
- W7-bis avec sieve numpy a couvert turns 0, 1, 2 (p_max = e^{6π})
- W7-ter : primesieve CLI subprocess streamé en chunks de range pour turn 3
- Sieve typique : 1 chunk = range de 1×10⁹ contiguës → ~5×10⁷ primes,
  ~10 sec/chunk, ~80 chunks total = ~13 min
- Accumulation per σ via streaming, mémoire bornée

PARSING
=======
- primesieve subprocess Popen avec stdout=PIPE binaire
- lecture par chunks ~100 MB
- np.frombuffer + parse via split + np.array (fast)

PROTOCOLE
=========
1. Charger résultats W7-bis pour T_G^turn_0/1/2 aux σ d'intérêt
2. Pour σ ∈ {3.0, 3.3, 3.5, √(4π), 3.6, 3.7, 4.0} :
   - Streamer turn 3 via primesieve
   - Accumuler T_G^turn_3(σ)
3. Comparer σ_crit^(3) observé à √(4π) ≈ 3.5449
"""

import numpy as np
import subprocess
import json
import os
import time
from datetime import datetime

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

# Turn 3 boundaries
TURN_3_LOW = int(np.ceil(np.exp(6 * np.pi)))    # ≈ 1.535×10⁸ (already known)
TURN_3_HIGH = int(np.floor(np.exp(8 * np.pi)))  # ≈ 7.896×10¹⁰

print(f"Turn 3 : [e^{{6π}}, e^{{8π}}) = [{TURN_3_LOW:,}, {TURN_3_HIGH:,})")
print(f"  Range : {TURN_3_HIGH - TURN_3_LOW:,} ≈ 7.9×10¹⁰")
print(f"  PNT estime ~{(TURN_3_HIGH/np.log(TURN_3_HIGH) - TURN_3_LOW/np.log(TURN_3_LOW)):,.0f} premiers")
print()

# Test σ : autour σ_crit^(3) = √(4π) ≈ 3.5449
SIGMAS = [3.0, 3.3, 3.5, 3.5449, 3.6, 3.7, 3.8, 4.0, 4.5]

CHUNK_RANGE = 1_000_000_000  # 10⁹ par chunk (≈ 5×10⁷ primes, ~400 MB int64)
TWO_PI = 2 * np.pi

# -----------------------------------------------------------------------------
# Streamed T_G accumulator for turn 3 via primesieve subprocess
# -----------------------------------------------------------------------------

def T_G_chunk(primes_arr, sigma, j_max=1):
    """Calcule T_G contribution pour un chunk de primes (vectorisé numpy).
    Pour primes ≥ 1.5×10⁸ dans turn 3 : j=1 dominant (j=2 → exp(-2·log_p²/...) << 1)."""
    log_p = np.log(primes_arr.astype(np.float64))
    sqrt_p = np.sqrt(primes_arr.astype(np.float64))
    norm = 1.0 / (2 * sigma * np.sqrt(np.pi))
    # j=1 only (suffit pour primes turn 3, contributions j≥2 négligeables)
    g_1 = norm * np.exp(-log_p ** 2 / (4 * sigma ** 2))
    contributions = 2 * log_p / sqrt_p * g_1
    return float(contributions.sum())


def stream_turn3_T_G(sigmas, chunk_range=CHUNK_RANGE):
    """Streame primes turn 3 via primesieve, accumule T_G^turn_3(σ) pour chaque σ."""
    accumulators = {sigma: 0.0 for sigma in sigmas}
    n_primes_total = 0
    n_chunks = 0

    current_low = TURN_3_LOW

    while current_low < TURN_3_HIGH:
        current_high = min(current_low + chunk_range, TURN_3_HIGH)
        chunk_low_str = str(current_low)
        chunk_high_str = str(current_high)

        t0 = time.time()
        # Capture stdout via subprocess (bytes mode)
        result = subprocess.run(
            ['primesieve', chunk_low_str, chunk_high_str, '--print', '--quiet'],
            capture_output=True,
            check=True,
        )
        t_sieve = time.time() - t0

        # Parse stdout bytes → numpy int64 array
        # np.fromstring with sep='\n' fonctionne sur bytes décodés
        t0 = time.time()
        primes_arr = np.fromstring(result.stdout.decode('ascii'), sep='\n', dtype=np.int64)
        t_parse = time.time() - t0

        n_primes_in_chunk = len(primes_arr)
        n_primes_total += n_primes_in_chunk

        # Compute T_G contributions pour chaque σ
        t0 = time.time()
        for sigma in sigmas:
            accumulators[sigma] += T_G_chunk(primes_arr, sigma)
        t_compute = time.time() - t0

        n_chunks += 1
        print(f"  Chunk {n_chunks:3d} [{current_low:,} → {current_high:,}] : "
              f"{n_primes_in_chunk:,} primes, "
              f"sieve {t_sieve:.1f}s + parse {t_parse:.1f}s + compute {t_compute:.1f}s = "
              f"{t_sieve+t_parse+t_compute:.1f}s")

        # Free memory before next chunk
        del primes_arr, result

        current_low = current_high

    print(f"\nTotal turn 3 : {n_primes_total:,} primes en {n_chunks} chunks.")
    return accumulators, n_primes_total


# -----------------------------------------------------------------------------
# Load W7-bis results pour T_G^turn_0
# -----------------------------------------------------------------------------

def load_T_G_turn_0(sigmas):
    """Charge T_G^turn_0(σ) depuis w7bis_extended_sieve_results.json ou recalcule."""
    out_dir = os.path.join(os.path.dirname(__file__), "..", "outputs")
    w7bis_path = os.path.join(out_dir, "w7bis_extended_sieve_results.json")
    if not os.path.exists(w7bis_path):
        raise FileNotFoundError(f"W7-bis output not found at {w7bis_path}. Run test_w7bis first.")
    with open(w7bis_path, "r") as f:
        w7bis = json.load(f)
    T_0 = {}
    for sigma in sigmas:
        # Match approximate σ (use closest available)
        sigma_str = str(sigma)
        if sigma_str in w7bis["by_sigma"]:
            T_0[sigma] = w7bis["by_sigma"][sigma_str]["T_G_turns"]["0"]
        else:
            # Need to compute or interpolate. For simplicity, recompute via primesieve sub-call.
            # Turn 0 = primes 2..523, fast.
            primes_turn0 = subprocess.run(
                ['primesieve', '2', '523', '--print', '--quiet'],
                capture_output=True, check=True
            )
            primes_arr = np.fromstring(primes_turn0.stdout.decode(), sep='\n', dtype=np.int64)
            T_0[sigma] = T_G_chunk(primes_arr, sigma, j_max=4)  # j_max=4 pour turn 0 (small primes)
    return T_0


# -----------------------------------------------------------------------------
# Predicted continuum ratio I_k / I_0 via erf
# -----------------------------------------------------------------------------

def predicted_ratio_exact_continuum(k, sigma):
    from scipy.special import erf
    sigma_u = sigma * np.sqrt(2)
    def I_k(k_):
        a = (TWO_PI * k_ - sigma**2) / sigma_u
        b = (TWO_PI * (k_ + 1) - sigma**2) / sigma_u
        return erf(b) - erf(a)
    return I_k(k) / I_k(0)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    print("=" * 100)
    print("W7-ter — Validation k=3 (σ_crit^(3) = √(4π) ≈ 3.5449) via primesieve")
    print("=" * 100)
    print(f"σ testés : {SIGMAS}")
    print(f"Chunk range : {CHUNK_RANGE:,} (≈ {CHUNK_RANGE / 2e7:.0f}M primes/chunk)")
    print()

    # Calcul T_G^turn_0 pour chaque σ (rapide, primes 2-523)
    print("Calcul T_G^turn_0 (recalcul rapide, primes 2-523)...")
    primes_turn0_result = subprocess.run(
        ['primesieve', '2', '523', '--print', '--quiet'],
        capture_output=True, check=True
    )
    primes_turn0 = np.fromstring(primes_turn0_result.stdout.decode(), sep='\n', dtype=np.int64)
    print(f"  Turn 0 : {len(primes_turn0)} primes (premier={primes_turn0[0]}, "
          f"dernier={primes_turn0[-1]})")

    T_G_turn0 = {}
    for sigma in SIGMAS:
        # Inclure j=1..4 pour turn 0 (small primes, j>1 contribue)
        log_p = np.log(primes_turn0.astype(np.float64))
        sqrt_p = np.sqrt(primes_turn0.astype(np.float64))
        p = primes_turn0.astype(np.float64)
        norm = 1.0 / (2 * sigma * np.sqrt(np.pi))
        total = 0.0
        for j in range(1, 5):
            jlogp = j * log_p
            g_j = norm * np.exp(-jlogp ** 2 / (4 * sigma ** 2))
            if j == 1:
                p_pow = sqrt_p
            elif j == 2:
                p_pow = p
            elif j == 3:
                p_pow = p * sqrt_p
            else:
                p_pow = p ** (j / 2)
            total += float((2 * log_p / p_pow * g_j).sum())
        T_G_turn0[sigma] = total
    print(f"  T_G^turn_0 calculé pour {len(SIGMAS)} valeurs σ.")
    print()

    # Streamage turn 3 via primesieve
    print(f"Streaming turn 3 via primesieve subprocess...")
    print(f"  Total chunks attendus : ~{(TURN_3_HIGH - TURN_3_LOW) // CHUNK_RANGE + 1}")
    print()

    t_start = time.time()
    T_G_turn3, n_primes_turn3 = stream_turn3_T_G(SIGMAS)
    t_total = time.time() - t_start

    print(f"\nTemps total turn 3 : {t_total:.1f}s ({t_total/60:.1f} min)")
    print(f"Primes turn 3 : {n_primes_turn3:,}")
    print()

    # Tableau résultats
    print("=" * 100)
    print(f"{'σ':>9} {'T_G^t0':>14} {'T_G^t3':>14} {'r_3/0_obs':>12} {'r_3/0_cont':>12} {'r_obs/cont':>12}")
    print("-" * 100)

    results = {
        "metadata": {
            "date": datetime.now().isoformat(),
            "turn_3_range": [TURN_3_LOW, TURN_3_HIGH],
            "n_primes_turn_3": n_primes_turn3,
            "chunk_range": CHUNK_RANGE,
            "time_seconds": t_total,
        },
        "predicted_sigma_crit_k3": float(np.sqrt(4 * np.pi)),
        "by_sigma": {}
    }

    for sigma in SIGMAS:
        T_0 = T_G_turn0[sigma]
        T_3 = T_G_turn3[sigma]
        r_obs = T_3 / T_0 if T_0 != 0 else float("inf")
        r_cont = float(predicted_ratio_exact_continuum(3, sigma))
        r_ratio = r_obs / r_cont if r_cont != 0 else float("inf")
        print(f"{sigma:>9.4f} {T_0:>14.4e} {T_3:>14.4e} "
              f"{r_obs:>12.4f} {r_cont:>12.4f} {r_ratio:>12.4f}")
        results["by_sigma"][str(sigma)] = {
            "sigma": sigma,
            "T_G_turn_0": T_0,
            "T_G_turn_3": T_3,
            "ratio_observed": r_obs,
            "ratio_continuum": r_cont,
            "ratio_obs_over_cont": r_ratio,
        }

    # Find σ_crit^(3) observé par interpolation
    print()
    print("=" * 100)
    print("SEUIL CRITIQUE σ_crit^(3)")
    print("=" * 100)

    sigma_crit_pred = float(np.sqrt(4 * np.pi))
    sigmas_sorted = sorted(SIGMAS)
    sigma_crit_obs = None
    for i in range(len(sigmas_sorted) - 1):
        s1, s2 = sigmas_sorted[i], sigmas_sorted[i+1]
        r1 = results["by_sigma"][str(s1)]["ratio_observed"]
        r2 = results["by_sigma"][str(s2)]["ratio_observed"]
        if (r1 - 1) * (r2 - 1) < 0:
            sigma_crit_obs = s1 + (s2 - s1) * (1 - r1) / (r2 - r1)
            break

    if sigma_crit_obs is not None:
        ecart = (sigma_crit_obs - sigma_crit_pred) / sigma_crit_pred * 100
        print(f"σ_crit^(3) prédit = √(4π) = {sigma_crit_pred:.4f}")
        print(f"σ_crit^(3) observé = {sigma_crit_obs:.4f}")
        print(f"Écart relatif : {ecart:+.3f}%")
        results["sigma_crit_observed"] = sigma_crit_obs
        results["sigma_crit_relative_error"] = (sigma_crit_obs - sigma_crit_pred) / sigma_crit_pred
    else:
        print("σ_crit^(3) non traversé dans la plage SIGMAS testée.")

    # Save
    out_dir = os.path.join(os.path.dirname(__file__), "..", "outputs")
    out_path = os.path.join(out_dir, "w7ter_primesieve_turn3_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nRésultats sauvegardés : {os.path.relpath(out_path)}")

    return results


if __name__ == "__main__":
    main()
