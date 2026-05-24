"""
PTM_Primev41_inverse.py

Inverse spectroscopy: given a sequence S, compute Lambda_q^(p) across all
characters of order q | p-1 for several small p, then recover the subgroup
H_p such that {n mod p : n in S} ⊂ H_p (modulo noise).

Theorem (Inverse Spectroscopy):
  For S arithmetic sequence with residues mod p uniform on H_p ⊂ (Z/pZ)*,
  Λ_q^(p)(S) is non-zero iff χ_q is trivial on H_p, i.e. order(χ_q) | [G:H_p].

Algorithm:
  H_p = ⋂ {ker(χ_q) : |Λ_q^(p)(S)| > threshold}

Applications:
  1. k-th power sums: recover k from Im(x^k).
  2. Forbidden-residue families: recover the forbidden set.
  3. New sequences: classify their arithmetic structure.
"""
from __future__ import annotations
import math
import cmath
from typing import Callable, Dict, Iterable, List, Set, Tuple

from PTM_Primev23_universal import sieve_with_spf, find_generator
from PTM_Primev34_extended_spectrometer import log_table_order_k


LOG_START = math.log(100)
NUM_BINS = 30
THRESHOLD = 0.5  # |Λ_q^(p)| threshold to consider non-trivial signal


def _bin_edges(N_max: int, num_bins: int = NUM_BINS) -> List[float]:
    log_max = math.log(N_max)
    return [math.exp(LOG_START + i * (log_max - LOG_START) / num_bins)
            for i in range(num_bins + 1)]


def lambda_char(sequence: Iterable[int], p: int, order: int,
                N_max: int, num_bins: int = NUM_BINS) -> complex:
    """Compute Λ_order^(p)(S) = logarithmic integral of character χ_order mod p."""
    edges = _bin_edges(N_max, num_bins)
    log_max = math.log(N_max)
    g = find_generator(p)
    table = log_table_order_k(p, g, order)
    sums = [0j] * num_bins
    counts = [0] * num_bins
    for n in sequence:
        if n < 7 or n > N_max:
            continue
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
            continue
        log_n = math.log(n)
        if log_n < LOG_START:
            continue
        idx = int((log_n - LOG_START) / (log_max - LOG_START) * num_bins)
        if not 0 <= idx < num_bins:
            continue
        r = n % p
        if r == 0:
            continue
        sums[idx] += cmath.exp(2j * math.pi * table[r] / order)
        counts[idx] += 1
    cum = 0j
    for i in range(num_bins):
        if counts[i] > 0:
            cum += (sums[i] / counts[i]) * (math.log(edges[i+1]) - math.log(edges[i]))
    return cum


def divisors_of(n: int) -> List[int]:
    divs = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            divs.append(d)
            if d != n // d:
                divs.append(n // d)
        d += 1
    return sorted(divs)


def signature_full(sequence_factory: Callable[[], Iterable[int]],
                   primes: List[int], N_max: int) -> Dict[Tuple[int, int], complex]:
    """For each p in primes, compute Λ_q^(p) for ALL q > 1 dividing p-1."""
    sig = {}
    for p in primes:
        qs = [q for q in divisors_of(p - 1) if q > 1]
        for q in qs:
            sig[(p, q)] = lambda_char(sequence_factory(), p, q, N_max)
    return sig


def recover_coset(sequence_factory: Callable[[], Iterable[int]],
                  p: int, N_max: int,
                  threshold: float = THRESHOLD) -> Set[int]:
    """Recover support set H_p ⊂ (Z/pZ)* from phase-resolved character signals.

    Theory (refined): for sequence S concentrated on coset H ⊂ (Z/pZ)*, the
    character Λ_q^(p)(S) has:
      - magnitude |Λ| ≈ log_range when q | [G:H]  (χ_q constant on H)
      - magnitude |Λ| ≈ 0         when q ∤ [G:H]  (χ_q averages to 0 on H)
    The phase of Λ identifies WHICH coset of ker(χ_q) contains H.

    So for each q with non-trivial signal:
      coset_q = {r ∈ (Z/pZ)* : χ_q(r) = e^(i·phase)}
    Then H = ⋂_q coset_q.
    """
    g = find_generator(p)
    qs = [q for q in divisors_of(p - 1) if q > 1]
    H = set(range(1, p))  # start with full (Z/pZ)*
    for q in qs:
        val = lambda_char(sequence_factory(), p, q, N_max)
        if abs(val) > threshold:
            observed_phase = cmath.phase(val)  # in radians
            table = log_table_order_k(p, g, q)
            # Accept residues r whose character value is close to observed phase
            coset_q = set()
            for r in range(1, p):
                r_phase = 2 * math.pi * table[r] / q
                delta = (r_phase - observed_phase) % (2 * math.pi)
                delta = min(delta, 2 * math.pi - delta)
                if delta < math.pi / q + 0.1:  # tolerant
                    coset_q.add(r)
            H &= coset_q
    return H if H else set(range(1, p))


def classify_k_power_sums(signature: Dict[Tuple[int, int], complex]) -> Dict:
    """Given signature of a k-th power sum sequence, recover k.

    For n = a_1^k + ... + a_m^k restricted to (Z/pZ)*:
    The image of x→x^k has size (p-1)/gcd(k, p-1). Characters trivial on
    this image are those of order dividing gcd(k, p-1).

    Deduce gcd(k, p-1) per prime, then CRT to get k modulo lcm(p-1).
    """
    deductions = {}
    for (p, q), val in signature.items():
        if abs(val) > THRESHOLD:
            deductions.setdefault(p, set()).add(q)

    # For each p: signal orders q are those with q | gcd(k, p-1)
    # max such q tells us gcd(k, p-1) (if all divisors included).
    gcds = {}
    for p, qs in deductions.items():
        # gcd(k, p-1) must be a common multiple of all q in qs, and divisor of p-1
        candidates = [d for d in divisors_of(p - 1) if all(q | d for q in qs)]
        # pick the smallest consistent divisor
        gcds[p] = min(candidates) if candidates else 1
    return {'gcds': gcds, 'deductions': {p: sorted(qs) for p, qs in deductions.items()}}


# ---------- Generators for test sequences ----------

def gen_kth_powers(k: int, N_max: int):
    """Generate n = a^k with a positive integer, n ≤ N_max."""
    out = []
    a = 1
    while True:
        n = a ** k
        if n > N_max:
            break
        if n >= 7:
            out.append(n)
        a += 1
    return out


def gen_kth_power_sums_m3(is_prime: List[bool], k: int, N_max: int, C: int = 50):
    """Generate n = a^k + b^k + c^k with -C ≤ a ≤ b ≤ c ≤ C, 7 ≤ n ≤ N_max."""
    seen = set()
    for a in range(-C, C + 1):
        for b in range(a, C + 1):
            for c in range(b, C + 1):
                n = a**k + b**k + c**k
                if 7 <= n <= N_max and n not in seen:
                    seen.add(n)
    return sorted(seen)


def gen_primes_restricted(is_prime: List[bool], N_max: int,
                          residues_mod_5: Set[int]):
    out = []
    for n in range(7, N_max + 1):
        if is_prime[n] and n % 5 in residues_mod_5:
            out.append(n)
    return out


def gen_primes_restricted_modp(is_prime: List[bool], N_max: int,
                                p: int, allowed: Set[int]):
    out = []
    for n in range(7, N_max + 1):
        if is_prime[n] and n % p in allowed:
            out.append(n)
    return out


# ---------- Main demonstration ----------

def main(N_max: int = 1_000_000):
    print("=" * 100)
    print(f"PTM_Primev41 — INVERSE SPECTROSCOPY demonstration, N_max = {N_max:,}")
    print("=" * 100)
    is_prime, _ = sieve_with_spf(N_max)
    primes_list = [5, 7, 11, 13]  # for signature reading

    # ---------- Test 1: recover k from single k-th powers ----------
    # Theorem applies: n = a^k mod p uniform on Im(x^k), character trivial iff q | gcd(k, p-1)
    print("\n--- Test 1: recover k from single k-th powers (n = a^k) ---")
    print(f"{'k':>3} | {'nontrivial signals per p':>40} | {'expected gcd(k,p-1)':>30} | verdict")
    print("-" * 100)
    for k in (2, 3, 4, 5, 6):
        seq = gen_kth_powers(k, N_max)
        if len(seq) < 30:
            print(f"  k={k}: too few points ({len(seq)}), skipping")
            continue
        factory = lambda s=seq: s
        sig = signature_full(factory, primes_list, N_max)
        classif = classify_k_power_sums(sig)
        nontrivial = classif['deductions']
        gcds = classif['gcds']
        expected = {p: math.gcd(k, p - 1) for p in primes_list}
        consistent = gcds == expected
        print(f"{k:>3} | {str(nontrivial):>40} | "
              f"exp={expected} obs={gcds}  "
              f"{'OK' if consistent else 'MISMATCH'}")

    # ---------- Test 2: primes with restricted residue mod 5 ----------
    print("\n--- Test 2: recover forbidden residue set from restricted primes ---")
    test_restrictions = [
        ("primes == 1,4 mod 5 (QR)", {1, 4}),
        ("primes == 2,3 mod 5 (NQR)", {2, 3}),
        ("primes == 1 only mod 5", {1}),
        ("primes all residues mod 5", {1, 2, 3, 4}),
    ]
    for name, allowed in test_restrictions:
        seq = gen_primes_restricted(is_prime, N_max, allowed)
        if len(seq) < 200:
            print(f"  {name}: too few primes ({len(seq)}), skip")
            continue
        factory = lambda s=seq: s
        H_recovered = recover_coset(factory, 5, N_max)
        print(f"  {name}: {len(seq)} points, recovered H_5 = {sorted(H_recovered)}, "
              f"truth = {sorted(allowed)}")

    # ---------- Test 3: dense subgroup recovery mod 13 (cubic residues) ----------
    print("\n--- Test 3: recover coset on (Z/13Z)* via dense subsample ---")
    # Cubic residues mod 13 = {1, 5, 8, 12} (subgroup of order 4, index 3)
    # Quadratic residues mod 13 = {1,3,4,9,10,12} (subgroup of order 6, index 2)
    test_cases = [
        ("primes ≡ cubic residue mod 13", 13, {1, 5, 8, 12}),
        ("primes ≡ QR mod 13",            13, {1, 3, 4, 9, 10, 12}),
        ("primes ≡ 1 only mod 7",         7,  {1}),
        ("primes ≡ {1,6} mod 7 (squares)", 7, {1, 2, 4}),  # squares mod 7
    ]
    for name, p_mod, allowed in test_cases:
        seq = gen_primes_restricted_modp(is_prime, N_max, p_mod, allowed)
        if len(seq) < 500:
            print(f"  {name}: too few points ({len(seq)}), skip")
            continue
        factory = lambda s=seq: s
        H = recover_coset(factory, p_mod, N_max)
        ok = H == allowed
        print(f"  {name}: {len(seq)} pts,  recovered H_{p_mod} = {sorted(H)}, "
              f"truth = {sorted(allowed)}  {'OK' if ok else 'MISMATCH'}")

    # ---------- Test 4: sanity — all primes (no restriction) ----------
    print("\n--- Test 4: all primes sanity check ---")
    all_primes = [n for n in range(7, N_max + 1) if is_prime[n]]
    for p_mod in (5, 7, 11, 13):
        factory = lambda: all_primes
        H_p = recover_coset(factory, p_mod, N_max)
        full = set(range(1, p_mod))
        ok = H_p == full
        print(f"  all primes, p={p_mod}: {len(all_primes)} pts, H_{p_mod} = {sorted(H_p)}  "
              f"{'OK' if ok else 'MISMATCH (signal over threshold ⇒ spurious)'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
