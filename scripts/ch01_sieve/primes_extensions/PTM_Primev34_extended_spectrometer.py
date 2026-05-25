"""
PTM_Primev34_extended_spectrometer.py

Extension du spectromètre 4D v32 sur trois axes simultanés:

A. Canaux quadratiques supplémentaires : mod 13, 17, 19, 23
   (en plus des existants 5, 7, 11)

B. Caractères d'ordre supérieur mod 7 :
   - Cubique (ordre 3)  : 3 directions cardinales sur S¹
   - Sextique (ordre 6) : 6 directions cardinales sur S¹
   (Z/7Z)* ≅ Z/6Z, donc ces caractères existent.

C. Moments d'ordre supérieur : variance, skewness per bin du biais Φ
   pour chaque canal quadratique. Objectif : détecter fluctuations
   invisibles à l'intégrale logarithmique.

Signature étendue (par sequence) :
    Σ_ext = ( λ_5, λ_7, λ_11, λ_13, λ_17, λ_19, λ_23,                 [7 scalaires R]
              arg Λ_4^(5), arg Λ_3^(7), arg Λ_6^(7),                  [3 phases S¹]
              var Φ^(5), var Φ^(7), var Φ^(11) )                      [3 moments R⁺]

Test : peut-on séparer des familles classiques qui collisionnent en 4D
(cousin vs safe) à partir des nouveaux canaux/caractères ?
"""
from __future__ import annotations
import math
import cmath
from typing import Callable, Dict, Iterable, List, Tuple

from PTM_Primev23_universal import sieve_with_spf, qr_set, find_generator
from PTM_Primev29_spectrometer import (
    gen_primes, gen_twin_primes, gen_sophie_germain
)
from PTM_Primev30_extended_library import (
    gen_cousin_primes, gen_sexy_primes, gen_safe_primes
)


S_CONST = 0.5            # normalisation (1/2)
LOG_START = math.log(100)
NUM_BINS = 30


def log_table_order_k(p: int, g: int, k: int) -> List[int]:
    """table[r] = log_g(r) mod k for r in 0..p-1 ; -1 if r=0.

    k must divide p-1 so that the character of order k is well-defined
    (non-trivial when k > 1 and k | p-1).
    """
    assert (p - 1) % k == 0, f"order {k} does not divide {p-1}"
    table = [-1] * p
    val = 1
    for j in range(p - 1):
        table[val] = j % k
        val = (val * g) % p
    return table


def _bin_edges(N_max: int, num_bins: int = NUM_BINS) -> List[float]:
    log_max = math.log(N_max)
    return [math.exp(LOG_START + i * (log_max - LOG_START) / num_bins)
            for i in range(num_bins + 1)]


def compute_lambda_quadratic_multi(integers: Iterable[int],
                                   channels: List[int],
                                   N_max: int,
                                   num_bins: int = NUM_BINS
                                   ) -> Dict[int, float]:
    """Compute lambda_p for each p in `channels` in a single pass.

    Also returns per-channel variance of Phi across bins (moment d'ordre 2
    empirique, Part C of v34).
    """
    edges = _bin_edges(N_max, num_bins)
    qrs = {p: qr_set(p) for p in channels}
    counts = {p: [[0, 0] for _ in range(num_bins)] for p in channels}
    n_total = 0
    for n in integers:
        if n < 7 or n > N_max:
            continue
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
            continue
        log_n = math.log(n)
        if log_n < LOG_START:
            continue
        idx = int((log_n - LOG_START) / (math.log(N_max) - LOG_START) * num_bins)
        if not 0 <= idx < num_bins:
            continue
        n_total += 1
        for p in channels:
            r = n % p
            if r == 0:
                continue
            counts[p][idx][0 if r in qrs[p] else 1] += 1
    results: Dict = {'n_total': n_total, 'channels': channels}
    for p in channels:
        phis, weights = [], []
        cum = 0.0
        for i in range(num_bins):
            nq, nn = counts[p][i]
            phi = math.log2(nn / nq) if (nq > 0 and nn > 0) else 0.0
            d_log = math.log(edges[i + 1]) - math.log(edges[i])
            cum += phi * d_log
            phis.append(phi)
            weights.append(d_log)
        results[p] = cum / S_CONST
        # variance pondérée de Phi across bins
        w_total = sum(weights) or 1.0
        mean = sum(p * w for p, w in zip(phis, weights)) / w_total
        var = sum(w * (p - mean) ** 2 for p, w in zip(phis, weights)) / w_total
        results[f'var_{p}'] = var
    return results


def compute_lambda_char(integers: Iterable[int],
                        p: int,
                        order: int,
                        N_max: int,
                        num_bins: int = NUM_BINS) -> Dict[str, float]:
    """Compute complex lambda using character of given `order` modulo `p`.

    Returns magnitude, phase (deg), and cumulative complex integral.
    """
    edges = _bin_edges(N_max, num_bins)
    g = find_generator(p)
    table = log_table_order_k(p, g, order)
    sums = [0j] * num_bins
    counts = [0] * num_bins
    for n in integers:
        if n < 7 or n > N_max:
            continue
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
            continue
        log_n = math.log(n)
        if log_n < LOG_START:
            continue
        idx = int((log_n - LOG_START) / (math.log(N_max) - LOG_START) * num_bins)
        if not 0 <= idx < num_bins:
            continue
        r = n % p
        if r == 0:
            continue
        k = table[r]
        psi = cmath.exp(2j * math.pi * k / order)
        sums[idx] += psi
        counts[idx] += 1
    cum = 0j
    for i in range(num_bins):
        if counts[i] > 0:
            mean_psi = sums[i] / counts[i]
            d_log = math.log(edges[i + 1]) - math.log(edges[i])
            cum += mean_psi * d_log
    mag = abs(cum)
    phase_deg = math.degrees(cmath.phase(cum)) if mag > 1e-12 else 0.0
    return {'magnitude': mag, 'phase_deg': phase_deg, 'cum': cum,
            'mod': p, 'order': order}


def signature_extended(factory: Callable[[], Iterable[int]],
                       N_max: int) -> Dict:
    """Extended signature : 7 quadratic channels + 3 higher-order phases + variances."""
    quad_channels = [5, 7, 11, 13, 17, 19, 23]
    quad = compute_lambda_quadratic_multi(factory(), quad_channels, N_max)
    # Higher-order characters (each needs one pass over the factory)
    q4_5 = compute_lambda_char(factory(), 5, 4, N_max)
    q3_7 = compute_lambda_char(factory(), 7, 3, N_max)
    q6_7 = compute_lambda_char(factory(), 7, 6, N_max)
    return {
        'n_total': quad['n_total'],
        'lambda': {p: quad[p] for p in quad_channels},
        'variance': {p: quad[f'var_{p}'] for p in quad_channels},
        'q4_5': q4_5, 'q3_7': q3_7, 'q6_7': q6_7,
    }


# ============================================================
# Theoretical predictions for cubic/sextic on primes
# ============================================================

def theoretical_phase_for_forbidden(p_mod: int, order: int,
                                    forbidden_residues: List[int]) -> float:
    """Predicted phase (deg) of complex lambda on P_A where A = (Z/pZ)* \\ forbidden.

    Formula derived in analogy with Theorem 5 (quartic):
      sum_{r in A} chi(r) = - sum_{r in forbidden} chi(r)
    so the phase = pi + arg(sum of chi over forbidden).
    """
    g = find_generator(p_mod)
    table = log_table_order_k(p_mod, g, order)
    total = 0j
    for r in forbidden_residues:
        if r % p_mod == 0:
            continue
        k = table[r % p_mod]
        total += cmath.exp(2j * math.pi * k / order)
    forbidden_sum = total  # sum over forbidden
    predicted = -forbidden_sum  # since sum over all non-zero is 0
    return math.degrees(cmath.phase(predicted)) if abs(predicted) > 1e-12 else 0.0


def main(N_max: int = 10_000_000):
    print("=" * 110)
    print(f"PTM_Primev34 — EXTENDED spectrometer")
    print(f"  A. Quadratic channels: 5, 7, 11 (old) + 13, 17, 19, 23 (new)")
    print(f"  B. Higher-order characters: quartic mod 5, cubic mod 7, sextic mod 7")
    print(f"  C. Variance of Phi per bin (moment 2)")
    print(f"  N_max = {N_max:,}")
    print("=" * 110)

    print("\nBuilding sieve...")
    is_prime, _ = sieve_with_spf(N_max)
    print("Done.\n")

    test_sequences = [
        ("primes (ref)",   lambda: gen_primes(is_prime, 7, N_max)),
        ("twin",           lambda: gen_twin_primes(is_prime, 7, N_max)),
        ("cousin",         lambda: gen_cousin_primes(is_prime, 7, N_max)),
        ("sexy",           lambda: gen_sexy_primes(is_prime, 7, N_max)),
        ("Sophie Germain", lambda: gen_sophie_germain(is_prime, 7, N_max)),
        ("safe",           lambda: gen_safe_primes(is_prime, 7, N_max)),
    ]

    sigs = {}
    for name, factory in test_sequences:
        print(f"Computing signature for '{name}'...")
        sigs[name] = signature_extended(factory, N_max)

    # ============ Report 1: 7 quadratic channels ============
    print("\n" + "=" * 110)
    print("QUADRATIC CHANNELS (7-dim)")
    print("=" * 110)
    header = f"{'sequence':>17} | {'count':>7} |" + \
             "".join(f"{f'λ_{p}':>9} |" for p in [5, 7, 11, 13, 17, 19, 23])
    print(header)
    print("-" * len(header))
    for name in sigs:
        s = sigs[name]
        line = f"{name:>17} | {s['n_total']:>7,} |"
        for p in [5, 7, 11, 13, 17, 19, 23]:
            line += f"{s['lambda'][p]:>+9.2f} |"
        print(line)

    # ============ Report 2: higher-order characters ============
    print("\n" + "=" * 110)
    print("HIGHER-ORDER CHARACTERS  (phase on S¹, magnitude = strength)")
    print("=" * 110)
    header = f"{'sequence':>17} |" + \
             f" {'|Λ_4⁽⁵⁾|':>9} {'arg Λ_4⁽⁵⁾':>11} |" + \
             f" {'|Λ_3⁽⁷⁾|':>9} {'arg Λ_3⁽⁷⁾':>11} |" + \
             f" {'|Λ_6⁽⁷⁾|':>9} {'arg Λ_6⁽⁷⁾':>11}"
    print(header)
    print("-" * len(header))
    for name in sigs:
        s = sigs[name]
        q45 = s['q4_5']; q37 = s['q3_7']; q67 = s['q6_7']
        line = f"{name:>17} |"
        line += f" {q45['magnitude']:>9.3f} {q45['phase_deg']:>+11.2f}° |"
        line += f" {q37['magnitude']:>9.3f} {q37['phase_deg']:>+11.2f}° |"
        line += f" {q67['magnitude']:>9.3f} {q67['phase_deg']:>+11.2f}°"
        print(line)

    # ============ Report 3: variance per channel ============
    print("\n" + "=" * 110)
    print("VARIANCE OF Φ PER CHANNEL  (higher moment = fluctuation strength)")
    print("=" * 110)
    header = f"{'sequence':>17} |" + \
             "".join(f"{f'var_{p}':>10} |" for p in [5, 7, 11, 13, 17, 19, 23])
    print(header)
    print("-" * len(header))
    for name in sigs:
        s = sigs[name]
        line = f"{name:>17} |"
        for p in [5, 7, 11, 13, 17, 19, 23]:
            line += f"{s['variance'][p]:>10.5f} |"
        print(line)

    # ============ Theoretical predictions ============
    print("\n" + "=" * 110)
    print("THEORETICAL CARDINAL PHASES  (forbidden-residue analogues of Thm. 5)")
    print("=" * 110)
    print("\n-- Quartic mod 5 (reference, should be 4 cardinals {0, +90, 180, -90}):")
    for r0 in (1, 2, 3, 4):
        phase = theoretical_phase_for_forbidden(5, 4, [r0])
        print(f"    forbid residue {r0} mod 5  → phase = {phase:+8.2f}°")

    print("\n-- Cubic mod 7 (3 cardinals: 0, +120, -120; each non-trivial residue class):")
    for r0 in range(1, 7):
        phase = theoretical_phase_for_forbidden(7, 3, [r0])
        print(f"    forbid residue {r0} mod 7  → phase = {phase:+8.2f}°")

    print("\n-- Sextic mod 7 (6 cardinals: 0, ±60, ±120, 180):")
    for r0 in range(1, 7):
        phase = theoretical_phase_for_forbidden(7, 6, [r0])
        print(f"    forbid residue {r0} mod 7  → phase = {phase:+8.2f}°")

    # ============ Disambiguation metric ============
    print("\n" + "=" * 110)
    print("PAIRWISE DISAMBIGUATION  (Euclidean distance on 10-dim signature)")
    print("=" * 110)

    def make_vec(s):
        vec = [s['lambda'][p] for p in [5, 7, 11, 13, 17, 19, 23]]
        # normalise phases to radians as 2D unit vectors
        for q in (s['q4_5'], s['q3_7'], s['q6_7']):
            phi = math.radians(q['phase_deg'])
            vec.extend([q['magnitude'] * math.cos(phi),
                        q['magnitude'] * math.sin(phi)])
        return vec

    names = list(sigs.keys())
    vecs = {name: make_vec(sigs[name]) for name in names}

    def dist(a, b):
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

    print(f"\n{'':>17} " + " ".join(f"{n:>15}" for n in names))
    for n1 in names:
        line = f"{n1:>17} "
        for n2 in names:
            line += f"{dist(vecs[n1], vecs[n2]):>15.3f} "
        print(line)

    print("\nDone.")


if __name__ == "__main__":
    main(N_max=10_000_000)
