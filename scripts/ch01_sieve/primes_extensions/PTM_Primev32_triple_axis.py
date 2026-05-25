"""
PTM_Primev32_triple_axis.py

TASK C : Spectromètre 3-canaux quadratique (mod 5, 7, 11) + caractère quartique mod 5.

Question : (lambda_5, lambda_7, lambda_11) suffit-il à disambiguer
twin / Sophie Germain / sexy / safe ? Sinon, le caractère quartique mod 5 le fait-il ?
"""
from __future__ import annotations
import math
import cmath
from typing import List, Set, Tuple
from PTM_Primev23_universal import sieve_with_spf, qr_set, find_generator
from PTM_Primev29_spectrometer import (
    gen_primes, gen_twin_primes, gen_sophie_germain
)
from PTM_Primev30_extended_library import (
    gen_cousin_primes, gen_sexy_primes, gen_safe_primes
)
from PTM_Primev23_universal import quartic_log_table as build_quartic_table


def compute_multi_lambda(integers, p_legendres: List[int],
                          N_max: int = 10_000_000, num_bins: int = 30,
                          log_start: float = math.log(100)) -> dict:
    """Compute lambda_p for each p in p_legendres in a single pass."""
    log_max = math.log(N_max)
    edges = [math.exp(log_start + i * (log_max - log_start) / num_bins)
             for i in range(num_bins + 1)]
    qrs = {p: qr_set(p) for p in p_legendres}
    counts = {p: [[0, 0] for _ in range(num_bins)] for p in p_legendres}
    n_total = 0

    for n in integers:
        if n < 7 or n > N_max:
            continue
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
            continue
        log_n = math.log(n)
        if log_n < log_start:
            continue
        bin_idx = int((log_n - log_start) / (log_max - log_start) * num_bins)
        if bin_idx < 0 or bin_idx >= num_bins:
            continue
        n_total += 1
        for p in p_legendres:
            r = n % p
            if r == 0:
                continue
            is_qr = r in qrs[p]
            counts[p][bin_idx][0 if is_qr else 1] += 1

    results = {}
    for p in p_legendres:
        cum = 0.0
        for i in range(num_bins):
            nq, nn = counts[p][i]
            phi = math.log2(nn / nq) if (nq > 0 and nn > 0) else 0.0
            d_log = math.log(edges[i + 1]) - math.log(edges[i])
            cum += phi * d_log
        results[p] = cum / 0.5
    results['n_total'] = n_total
    return results


def compute_quartic_lambda_mod5(integers, N_max: int = 10_000_000,
                                 num_bins: int = 30,
                                 log_start: float = math.log(100)) -> dict:
    """Compute complex lambda using quartic character mod 5."""
    log_max = math.log(N_max)
    edges = [math.exp(log_start + i * (log_max - log_start) / num_bins)
             for i in range(num_bins + 1)]
    g = find_generator(5)
    table = build_quartic_table(5, g)

    sums = [0j for _ in range(num_bins)]
    counts = [0 for _ in range(num_bins)]

    for n in integers:
        if n < 7 or n > N_max:
            continue
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
            continue
        log_n = math.log(n)
        if log_n < log_start:
            continue
        bin_idx = int((log_n - log_start) / (log_max - log_start) * num_bins)
        if bin_idx < 0 or bin_idx >= num_bins:
            continue
        r = n % 5
        if r == 0:
            continue
        k = table[r]
        psi = cmath.exp(1j * math.pi * k / 2)
        sums[bin_idx] += psi
        counts[bin_idx] += 1

    cum = 0j
    for i in range(num_bins):
        if counts[i] > 0:
            mean_psi = sums[i] / counts[i]
            d_log = math.log(edges[i + 1]) - math.log(edges[i])
            cum += mean_psi * d_log
    return {'cum_complex': cum, 'magnitude': abs(cum),
            'phase_deg': math.degrees(cmath.phase(cum)) if abs(cum) > 1e-12 else 0.0}


def main(N_max: int = 10_000_000):
    print("=" * 92)
    print(f"v32 TRIPLE-AXIS spectrometer (mod 5,7,11 + quartic mod 5), N_max = {N_max:,}")
    print("=" * 92)
    print("Building sieve...")
    is_prime, _ = sieve_with_spf(N_max)
    print("Done.")

    test_sequences = [
        ("primes (ref)", lambda: gen_primes(is_prime, 7, N_max)),
        ("twin (p, p+2)", lambda: gen_twin_primes(is_prime, 7, N_max)),
        ("cousin (p, p+4)", lambda: gen_cousin_primes(is_prime, 7, N_max)),
        ("sexy (p, p+6)", lambda: gen_sexy_primes(is_prime, 7, N_max)),
        ("Sophie Germain", lambda: gen_sophie_germain(is_prime, 7, N_max)),
        ("safe primes", lambda: gen_safe_primes(is_prime, 7, N_max)),
    ]

    p_channels = [5, 7, 11]
    print(f"\n{'sequence':>20} | {'count':>8} | {'lambda_5':>10} | {'lambda_7':>10} | "
          f"{'lambda_11':>10} | {'mag (q4_5)':>10} | {'phase q4_5':>10}")
    print("-" * 110)

    signatures = {}
    for name, factory in test_sequences:
        # Run once for quadratic, once for quartic (iterators consume)
        res_quad = compute_multi_lambda(factory(), p_channels, N_max=N_max)
        res_quartic = compute_quartic_lambda_mod5(factory(), N_max=N_max)
        signatures[name] = {
            **res_quad,
            'q4_mag': res_quartic['magnitude'],
            'q4_phase': res_quartic['phase_deg'],
        }
        print(f"{name:>20} | {res_quad['n_total']:>8,} | "
              f"{res_quad[5]:>+10.4f} | {res_quad[7]:>+10.4f} | {res_quad[11]:>+10.4f} | "
              f"{res_quartic['magnitude']:>10.4f} | {res_quartic['phase_deg']:>+10.2f}")

    # Disambiguation analysis
    print("\n" + "=" * 92)
    print("DISAMBIGUATION ANALYSIS")
    print("=" * 92)

    # Group by quadratic signs
    for name, sig in signatures.items():
        signs = (
            "+" if sig[5] > 0 else "-",
            "+" if sig[7] > 0 else "-",
            "+" if sig[11] > 0 else "-",
        )
        print(f"  {name:>20} : (sign5,sign7,sign11) = {signs}, "
              f"|q4_5|={sig['q4_mag']:.3f}, phase q4_5 = {sig['q4_phase']:+.1f}°")

    # Look at quartic phase as additional discriminator
    print("\nPHASE QUARTIQUE (mod 5) comme disambigateur final :")
    print("  La phase complexe ψ_4 distingue les résidus interdits MEME quand")
    print("  λ_5 quadratique est identique (ex: forbid 1 vs forbid 4 mod 5).")

    # Test theoretical predictions for forbidden residues mod 11
    print("\n" + "=" * 92)
    print("PREDICTIONS THEORIQUES MOD 11 (pour reference)")
    print("=" * 92)
    log_range = math.log(N_max) - math.log(100)
    qr11 = qr_set(11)
    for forbidden in range(1, 11):
        allowed = set(range(1, 11)) - {forbidden}
        n_qr = len(allowed & qr11)
        n_nqr = len(allowed - qr11)
        if n_qr > 0 and n_nqr > 0:
            phi = math.log2(n_nqr / n_qr)
            pred = phi * log_range / 0.5
            type_r = "QR" if forbidden in qr11 else "NQR"
            print(f"  forbid {forbidden:>2} mod 11 ({type_r}) : "
                  f"|allowed|=({n_qr} QR, {n_nqr} NQR), "
                  f"lambda_pred = {pred:+.2f}")


if __name__ == "__main__":
    main(N_max=10_000_000)
