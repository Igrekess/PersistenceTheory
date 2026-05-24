"""
PTM_Primev37_carmichael.py

Objectif
--------
Le spectromètre v34/v35 a révélé |Λ_4^(5)| = 3.35 à phase ≈ 0° sur les 34
Carmichael numbers dans la fenêtre log à N = 10^6. Par le Théorème 6
(forbidden-residue phase) avec k=4, p=5, générateur g=2, la phase 0°
correspond au résidu interdit r_0 = 4 (car chi_4(4) = chi_4(g^2) =
exp(2πi·2/4) = -1, donc -chi_4(4) = +1 ↔ phase 0°).

Prédiction : les Carmichael numbers évitent préférentiellement la classe
n ≡ 4 (mod 5).

Ce script :
  1. Génère tous les Carmichael ≤ 10^7 (Korselt).
  2. Tabule la distribution mod 5.
  3. Évalue la significativité statistique (p-value binomiale).
  4. Dérive l'argument théorique à partir de Korselt pour ω(n)=3.
  5. Teste la même signature sur les pseudo-premiers de Fermat base 2.

N'utilise PAS v33/v34/v35 comme dépendance pour ne pas les modifier.
"""
from __future__ import annotations
import math
from typing import Dict, Iterable, List, Tuple


# ============================================================
# Sieve utilities (local copy for autonomy)
# ============================================================

def sieve_with_spf(N: int) -> Tuple[List[bool], List[int]]:
    """Returns (is_prime, smallest_prime_factor) for [0..N]."""
    is_prime = [True] * (N + 1)
    if N >= 0:
        is_prime[0] = False
    if N >= 1:
        is_prime[1] = False
    spf = [0] * (N + 1)
    for i in range(2, N + 1):
        if spf[i] == 0:
            spf[i] = i
            for j in range(i * 2, N + 1, i):
                if spf[j] == 0:
                    spf[j] = i
                is_prime[j] = False
    return is_prime, spf


# ============================================================
# Carmichael generation (Korselt)
# ============================================================

def gen_carmichael(is_prime: List[bool], spf: List[int],
                   n_max: int) -> List[int]:
    """All Carmichael numbers in [2, n_max] via Korselt's criterion."""
    out: List[int] = []
    for n in range(9, n_max + 1):
        if n % 2 == 0:
            continue
        if is_prime[n]:
            continue
        # Factor via spf
        m = n
        primes_n: List[int] = []
        squarefree = True
        while m > 1:
            p = spf[m]
            primes_n.append(p)
            m //= p
            if m % p == 0:
                squarefree = False
                break
        if not squarefree:
            continue
        if len(primes_n) < 3:
            continue
        if any((n - 1) % (p - 1) != 0 for p in primes_n):
            continue
        out.append(n)
    return out


def factor_omega_and_primes(spf: List[int], n: int) -> List[int]:
    """Return distinct prime factors of n using spf."""
    primes_n: List[int] = []
    m = n
    while m > 1:
        p = spf[m]
        primes_n.append(p)
        while m % p == 0:
            m //= p
    return primes_n


# ============================================================
# Fermat pseudoprimes base 2 (A001567)
# ============================================================

def gen_fermat_pseudoprimes_base2(is_prime: List[bool], n_max: int) -> List[int]:
    """Composites n with 2^(n-1) ≡ 1 (mod n)."""
    out: List[int] = []
    for n in range(9, n_max + 1, 2):
        if is_prime[n]:
            continue
        if pow(2, n - 1, n) == 1:
            out.append(n)
    return out


# ============================================================
# Statistics
# ============================================================

def distribution_mod5(seq: List[int], exclude_div5: bool = True
                      ) -> Tuple[Dict[int, int], int]:
    """Return counts of n mod 5 and total (excluding class 0 if flag)."""
    counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    for n in seq:
        counts[n % 5] += 1
    total_nz = counts[1] + counts[2] + counts[3] + counts[4]
    return counts, total_nz


def _log_binom(n: int, k: int) -> float:
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def binomial_tail_pvalue(k_obs: int, n_total: int, p: float,
                         direction: str = "lower") -> float:
    """Two/one-sided binomial tail probability.

    direction='lower': P(X <= k_obs) under X ~ Binom(n_total, p).
    direction='upper': P(X >= k_obs).
    """
    assert direction in ("lower", "upper")
    if n_total == 0:
        return 1.0
    # Exact sum in log-space to keep precision for k up to n.
    total = 0.0
    rng = range(0, k_obs + 1) if direction == "lower" else range(k_obs, n_total + 1)
    for j in rng:
        lp = _log_binom(n_total, j) + j * math.log(p) + (n_total - j) * math.log(1 - p)
        total += math.exp(lp)
    return min(total, 1.0)


def chisq_uniform_4classes(counts: Dict[int, int]) -> Tuple[float, int]:
    """Chi-square for uniform distribution over classes {1,2,3,4}."""
    n = counts[1] + counts[2] + counts[3] + counts[4]
    if n == 0:
        return 0.0, 0
    expected = n / 4.0
    chi = sum((counts[r] - expected) ** 2 / expected for r in (1, 2, 3, 4))
    return chi, n


def chi2_pvalue_df3(chi2: float) -> float:
    """Upper-tail p-value for chi^2 with 3 degrees of freedom.

    Regularized upper incomplete gamma Q(3/2, x/2) ... but df=3 ⇒ Q(3/2, x/2).
    For df=3: P(chi^2 >= x) = erfc(sqrt(x/2)) + sqrt(2x/pi) * exp(-x/2)
    (derived from CDF of chi^2_3).
    """
    if chi2 <= 0:
        return 1.0
    x = chi2
    return math.erfc(math.sqrt(x / 2.0)) + math.sqrt(2.0 * x / math.pi) * math.exp(-x / 2.0)


# ============================================================
# Theoretical analysis : why Carmichael avoid class 4 mod 5
# ============================================================

def korselt_class_mod5_for_triple(p: int, q: int, r: int) -> int | None:
    """Given 3 distinct odd primes (p,q,r), test whether n = p*q*r
    satisfies Korselt ((p'-1) | (n-1) for all p' in {p,q,r}).
    If so, return n mod 5; else None.
    """
    if 5 in (p, q, r):
        return 0  # trivially class 0
    n = p * q * r
    if any((n - 1) % (x - 1) != 0 for x in (p, q, r)):
        return None
    return n % 5


def scan_small_carmichael_omega3(is_prime: List[bool], spf: List[int],
                                  n_max: int) -> Dict[int, int]:
    """For all Carmichael with ω=3 ≤ n_max, collect n mod 5 distribution
    restricted to gcd(n,5)=1."""
    counts = {1: 0, 2: 0, 3: 0, 4: 0, 0: 0}
    for n in gen_carmichael(is_prime, spf, n_max):
        pr = factor_omega_and_primes(spf, n)
        if len(pr) != 3:
            continue
        counts[n % 5] += 1
    return counts


# ============================================================
# Main analysis
# ============================================================

def analyse_sequence(name: str, seq: List[int]) -> None:
    counts, nz = distribution_mod5(seq)
    print(f"\n--- {name} ---")
    print(f"  Total n in sequence  : {len(seq):,}")
    print(f"  n ≡ 0 (mod 5) [div 5]: {counts[0]}")
    print(f"  Coprime-to-5 subset  : {nz}")
    print(f"  Distribution mod 5 over classes (1,2,3,4):")
    for r in (1, 2, 3, 4):
        pct = 100.0 * counts[r] / nz if nz > 0 else 0.0
        print(f"    class {r}: {counts[r]:>6,}   ({pct:5.2f} %)")
    if nz >= 20:
        chi2, n_used = chisq_uniform_4classes(counts)
        pval_uni = chi2_pvalue_df3(chi2)
        print(f"  χ²(df=3) vs uniform 1/4 : {chi2:.3f}  (p ≈ {pval_uni:.2e})")
        # One-sided tail: is class 4 depleted relative to uniform?
        p4_under = binomial_tail_pvalue(counts[4], nz, 0.25, "lower")
        p4_over  = binomial_tail_pvalue(counts[4], nz, 0.25, "upper")
        print(f"  Binom tail P(X<={counts[4]}|Binom({nz},0.25)) = {p4_under:.3e}")
        print(f"  Binom tail P(X>={counts[4]}|Binom({nz},0.25)) = {p4_over:.3e}")
    else:
        print(f"  Sample too small for chi²/binomial tests.")


def demonstrate_theory_omega3(is_prime: List[bool], spf: List[int],
                               n_max_scan: int) -> None:
    print("\n" + "=" * 80)
    print("THEORETICAL ANALYSIS  ω(n) = 3 Carmichael")
    print("=" * 80)

    # Collect all ω=3 Carmichael ≤ n_max_scan and tabulate mod 5
    counts = scan_small_carmichael_omega3(is_prime, spf, n_max_scan)
    print(f"\n  ω=3 Carmichael ≤ {n_max_scan:,}: "
          f"total {sum(counts.values())}, class 0 (div 5) = {counts[0]}")
    nz = sum(counts[r] for r in (1, 2, 3, 4))
    print(f"  Distribution over classes coprime-to-5 (n_total = {nz}):")
    for r in (1, 2, 3, 4):
        pct = 100.0 * counts[r] / nz if nz > 0 else 0.0
        print(f"    class {r}: {counts[r]:>5,} ({pct:5.2f} %)")

    # Symbolic argument
    print("\n  === Korselt → class-4 exclusion (informal argument) ===")
    print("  For n = p*q*r with p,q,r odd primes ≠ 5:")
    print("    - Korselt ⇒ (p-1)|(n-1), (q-1)|(n-1), (r-1)|(n-1).")
    print("    - In particular 2 | (n-1) (always true).")
    print("    - Key: many Carmichael have at least one factor ≡ 1 (mod 5)")
    print("      (so that p-1 ≡ 0 mod 5, hence 5 | (n-1), i.e. n ≡ 1 mod 5).")
    print("    - If no p_i ≡ 1 (mod 5), each p_i lies in {1,2,3,4}\\{1} = {2,3,4}")
    print("      mod 5, but since p_i ≠ 5, the p_i - 1 are in {1,2,3}.")
    print("    - lcm(p_i - 1) | (n-1) restricts (n-1) mod 5 in a way that rules")
    print("      out n-1 ≡ 3 (mod 5), i.e. n ≡ 4 (mod 5), with high prob.")

    # Numerical proof of concept: enumerate squarefree ω=3 triples with
    # max prime ≤ N_small and show the restriction.
    print("\n  === Brute small-triple census (p,q,r ≤ 200, distinct, odd, ≠5): ===")
    primes_small = [p for p in range(3, 201) if is_prime[p] and p != 5]
    counts_triple = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    total_valid = 0
    for i in range(len(primes_small)):
        for j in range(i + 1, len(primes_small)):
            for k in range(j + 1, len(primes_small)):
                p, q, r = primes_small[i], primes_small[j], primes_small[k]
                c = korselt_class_mod5_for_triple(p, q, r)
                if c is None:
                    continue
                total_valid += 1
                counts_triple[c] += 1
    print(f"  Found {total_valid} Korselt-valid ω=3 triples (p,q,r ≤ 200, ≠5).")
    if total_valid > 0:
        for r in range(5):
            pct = 100.0 * counts_triple[r] / total_valid
            print(f"    n ≡ {r} (mod 5) : {counts_triple[r]:>5,} ({pct:5.2f} %)")


def main(N_max: int = 10_000_000) -> None:
    print("=" * 90)
    print(f"PTM_Primev37 — Carmichael quartic-mod-5 alignment diagnosis")
    print(f"  Hypothesis: phase ≈ 0° in |Λ_4^(5)| ⇒ residue class 4 (mod 5) is avoided.")
    print(f"  N_max = {N_max:,}")
    print("=" * 90)

    print("\nBuilding sieve (this is the heavy step)...")
    is_prime, spf = sieve_with_spf(N_max)
    print(f"  Primes ≤ {N_max:,}: {sum(is_prime):,}")

    # Carmichael
    print("\nGenerating Carmichael numbers (Korselt)...")
    carm = gen_carmichael(is_prime, spf, N_max)
    print(f"  Found {len(carm):,} Carmichael in [0, {N_max:,}]")
    if carm:
        print(f"  First three: {carm[:3]}")
        print(f"  Last three : {carm[-3:]}")
    analyse_sequence("Carmichael (all ≤ N_max)", carm)

    # Carmichael restricted to ω = 3
    carm_w3 = [n for n in carm if len(factor_omega_and_primes(spf, n)) == 3]
    analyse_sequence("Carmichael (ω=3 only)", carm_w3)

    # Carmichael restricted to ω ≥ 4
    carm_w4plus = [n for n in carm if len(factor_omega_and_primes(spf, n)) >= 4]
    analyse_sequence("Carmichael (ω≥4)", carm_w4plus)

    # Fermat pseudoprimes base 2 (can be slow; use if N_max ≤ 10^7)
    if N_max <= 10_000_000:
        print("\nGenerating Fermat pseudoprimes base 2 ...")
        fpp = gen_fermat_pseudoprimes_base2(is_prime, N_max)
        print(f"  Found {len(fpp):,} Fermat pseudoprimes base 2.")
        analyse_sequence("Fermat pseudoprimes base 2", fpp)
    else:
        fpp = []

    # Theory
    demonstrate_theory_omega3(is_prime, spf, min(N_max, 1_000_000))

    # Conjecture
    print("\n" + "=" * 90)
    print("CONJECTURE (formalised)")
    print("=" * 90)
    print("""
  CONJECTURE (Carmichael mod 5 depletion).
  Let C(N) = #{Carmichael ≤ N, gcd(n,5) = 1}, and
      C_r(N) = #{n ∈ Carmichael, n ≤ N, n ≡ r (mod 5)}   for r ∈ {1,2,3,4}.
  Then
      lim sup_{N→∞}  C_4(N) / C(N)   <   1/4   strictly,
  and equivalently the spectrometer quartic-mod-5 amplitude
      |Λ_4^(5)| (Carmichael)  →  a positive limit  with phase ≈ 0°.

  A (stronger) precise form is
      C_4(N) / C(N)  ≍  (log log N)^{-α}   for some α ≥ 1,
  driven by the Korselt constraint 5 | (n-1) whenever some p_i ≡ 1 (mod 5).
""")
    print("=" * 90)
    print("Done.")


if __name__ == "__main__":
    import sys
    N = 10_000_000
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    main(N_max=N)
