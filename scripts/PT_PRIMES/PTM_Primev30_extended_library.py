"""
PTM_Primev30_extended_library.py

Bibliotheque etendue du lambda-spectrometer.

Prediction theorique pour subset arithmetique S coprime a 30 :

  Phi(S, bin) = log2 ( #{n in S : (n/5)=-1} / #{n in S : (n/5)=+1} )

Pour primes asymptotiquement uniformes sur residus allowed_set ⊂ {1,2,3,4} mod 5 :
  N_QR  = |{1,4} ∩ allowed_set|
  N_NQR = |{2,3} ∩ allowed_set|
  Phi   ≈ log2(N_NQR / N_QR)
  I     = Phi × log_range
  lambda = I / s

Donc pour un sous-ensemble de primes avec residu r mod 5 INTERDIT :
  forbid 1 (QR): allowed=(2,3,4), Phi=log2(2/1)=+1, lambda ≈ +23 (a N=10^7)
  forbid 4 (QR): meme calcul, lambda ≈ +23
  forbid 2 (NQR): allowed=(1,3,4), Phi=log2(1/2)=-1, lambda ≈ -23
  forbid 3 (NQR): meme calcul, lambda ≈ -23

Tests :
  twin primes (p, p+2) : forbid p ≡ 3 mod 5 -> lambda ≈ -23
  cousin primes (p, p+4) : forbid p ≡ 1 mod 5 -> lambda ≈ +23
  Sophie Germain (p, 2p+1) : forbid p ≡ 2 mod 5 -> lambda ≈ -23
"""
from __future__ import annotations
import math
from typing import List, Iterable, Tuple, Dict, Set
from PTM_Primev23_universal import sieve_with_spf, qr_set
from PTM_Primev29_spectrometer import (
    compute_lambda, legendre5, gen_primes, gen_twin_primes,
    gen_sophie_germain, gen_carmichael_small, gen_arithmetic_progression
)


# ============================
# Extended library
# ============================

class LambdaLibrary:
    """
    Bibliotheque de predictions lambda pour structures arithmetiques.
    """

    def __init__(self, lambda_primes: float, log_range: float, s: float = 0.5):
        self.lambda_primes = lambda_primes
        self.log_range = log_range
        self.s = s

    def multiplicative(self, factors: Tuple[int, ...]) -> float:
        """omega-strate fixed factors. lambda = prod((q/5)) * lambda(primes)."""
        prod = 1
        for q in factors:
            leg = legendre5(q)
            if leg == 0:
                return 0.0
            prod *= leg
        return prod * self.lambda_primes

    def forbidden_residue(self, allowed_set: Set[int]) -> float:
        """
        For primes (or omega=1 subset) restricted to allowed_set ⊂ {1,2,3,4} mod 5.
        Asymptotic uniform distribution on allowed_set.
        """
        qr = allowed_set & {1, 4}
        nqr = allowed_set & {2, 3}
        n_qr = len(qr)
        n_nqr = len(nqr)
        if n_qr == 0 or n_nqr == 0:
            # All in one class - Phi undefined (sequence highly biased)
            return float('inf') if n_qr == 0 else float('-inf')
        ratio = n_nqr / n_qr
        phi = math.log2(ratio)
        return phi * self.log_range / self.s

    def combined(self, factors: Tuple[int, ...], allowed_set: Set[int]) -> float:
        """
        Combination : multiplicative on factors + forbidden residue on free factor.
        For S = {n : n = q1*...*qk-1*m, m prime, m mod 5 in allowed_set}
        Phi = legendre_prod * Phi_forbidden
        lambda = legendre_prod * forbidden_residue_lambda
        """
        prod = 1
        for q in factors:
            leg = legendre5(q)
            if leg == 0:
                return 0.0
            prod *= leg
        return prod * self.forbidden_residue(allowed_set)

    def identify(self, observed_lambda: float, candidates: List[Tuple[str, float]]) -> List[Tuple[str, float, float]]:
        """
        Given observed lambda, return candidates sorted by closeness.
        Each tuple: (name, predicted_lambda, abs_distance).
        """
        scored = [(name, pred, abs(pred - observed_lambda)) for name, pred in candidates]
        return sorted(scored, key=lambda x: x[2])


# ============================
# Generators for new test cases
# ============================

def gen_cousin_primes(is_prime: List[bool], lo: int, hi: int):
    """p such that p+4 also prime."""
    for n in range(max(lo, 3), hi - 3):
        if is_prime[n] and is_prime[n + 4]:
            yield n


def gen_sexy_primes(is_prime: List[bool], lo: int, hi: int):
    """p such that p+6 also prime."""
    for n in range(max(lo, 3), hi - 5):
        if is_prime[n] and is_prime[n + 6]:
            yield n


def gen_safe_primes(is_prime: List[bool], lo: int, hi: int):
    """p such that (p-1)/2 is also prime (i.e., p = 2q+1 with q prime)."""
    for n in range(max(lo, 5), hi):
        if is_prime[n] and (n - 1) % 2 == 0 and is_prime[(n - 1) // 2]:
            yield n


def gen_primes_residue_mod(is_prime: List[bool], lo: int, hi: int,
                            modulus: int, residue: int):
    """Primes p ≡ residue mod modulus."""
    for n in range(max(lo, 2), hi + 1):
        if is_prime[n] and n % modulus == residue:
            yield n


# ============================
# Compute predicted forbidden residue from gap conditions
# ============================

def forbidden_residue_for_pair_condition(gap: int, mod_p: int = 5) -> Set[int]:
    """
    For sequence {p prime : p+gap also prime, both coprime to 5},
    we need (p+gap) mod 5 != 0, i.e., p mod 5 != -gap mod 5.
    Returns the SET of allowed residues (excludes the forbidden one).
    """
    forbidden = (-gap) % mod_p
    if forbidden == 0:
        return set(range(1, mod_p))  # gap is multiple of 5, no extra constraint
    return set(range(1, mod_p)) - {forbidden}


def forbidden_residue_for_safe_germain(scale: int = 2, shift: int = 1, mod_p: int = 5) -> Set[int]:
    """
    For sequence {p prime : scale*p + shift also prime}, we need
    scale*p + shift != 0 mod 5, i.e., p != -shift/scale mod p.
    Returns allowed residues.

    For Sophie Germain (scale=2, shift=1):
      2p+1 != 0 mod 5 -> p != -1/2 = -3 = 2 mod 5. Forbid 2.
    For safe primes (q=(p-1)/2 prime):
      We need q = (p-1)/2 != 0 mod 5, i.e., p-1 != 0 mod 10, i.e., p != 1 mod 5.
      (Use scale=1, shift=-1 with division by 2, but in mod 5 simpler:
       p != 1 mod 5 directly.)
    """
    inv_scale = pow(scale, -1, mod_p)
    forbidden = (-shift * inv_scale) % mod_p
    if forbidden == 0:
        return set(range(1, mod_p))
    return set(range(1, mod_p)) - {forbidden}


def forbidden_residue_for_safe_primes(mod_p: int = 5) -> Set[int]:
    """
    Safe primes : p prime AND (p-1)/2 prime.
    Constraint: (p-1)/2 != 0 mod 5, i.e., p != 1 mod 5.
    """
    return set(range(1, mod_p)) - {1}


# ============================
# Main
# ============================

def main(N_max: int = 10_000_000):
    print("=" * 92)
    print(f"v30 EXTENDED LIBRARY spectrometer, N_max = {N_max:,}")
    print("=" * 92)
    print("Building sieve...")
    is_prime, spf = sieve_with_spf(N_max)
    print("Done.")

    log_range = math.log(N_max) - math.log(100)
    print(f"\nLog range: {log_range:.3f}")

    # Compute lambda(primes) reference
    lam_primes, n_p = compute_lambda(gen_primes(is_prime, 7, N_max), N_max=N_max)
    print(f"lambda(primes) = {lam_primes:+.4f}  (count = {n_p:,})")

    lib = LambdaLibrary(lambda_primes=lam_primes, log_range=log_range)

    print("\n" + "=" * 92)
    print("PREDICTIONS THEORIQUES (forbidden residue)")
    print("=" * 92)
    for forbidden in [1, 2, 3, 4]:
        allowed = {1, 2, 3, 4} - {forbidden}
        pred = lib.forbidden_residue(allowed)
        print(f"  forbid r={forbidden} mod 5 (allowed={sorted(allowed)}): "
              f"lambda predicted = {pred:+.2f}")

    print("\n" + "=" * 92)
    print("MESURES + COMPARAISON BIBLIOTHEQUE")
    print("=" * 92)

    # Use lambda factories to recreate generators on demand (iterators consume!)
    test_cases = [
        ("primes (ref)", lambda: gen_primes(is_prime, 7, N_max),
         {'allowed_5': {1, 2, 3, 4}, 'desc': 'all residues, lambda_primes'}),
        ("twin primes (p, p+2)", lambda: gen_twin_primes(is_prime, 7, N_max),
         {'allowed_5': forbidden_residue_for_pair_condition(2), 'desc': 'forbid 3 mod 5'}),
        ("cousin primes (p, p+4)", lambda: gen_cousin_primes(is_prime, 7, N_max),
         {'allowed_5': forbidden_residue_for_pair_condition(4), 'desc': 'forbid 1 mod 5'}),
        ("sexy primes (p, p+6)", lambda: gen_sexy_primes(is_prime, 7, N_max),
         {'allowed_5': forbidden_residue_for_pair_condition(6), 'desc': 'gap=6 ≡ 1 mod 5, forbid 4'}),
        ("Sophie Germain (p, 2p+1)", lambda: gen_sophie_germain(is_prime, 7, N_max),
         {'allowed_5': forbidden_residue_for_safe_germain(2, 1), 'desc': 'forbid 2 mod 5'}),
        ("safe primes ((p-1)/2, p)", lambda: gen_safe_primes(is_prime, 7, N_max),
         {'allowed_5': forbidden_residue_for_safe_primes(), 'desc': 'forbid 1 mod 5'}),
    ]

    print(f"\n{'sequence':>30} | {'count':>8} | {'lambda obs':>11} | "
          f"{'lambda pred':>11} | {'diff':>8} | {'desc':>20}")
    print("-" * 100)
    for name, seq_factory, info in test_cases:
        lam, n = compute_lambda(seq_factory(), N_max=N_max)
        if 'allowed_5' in info:
            allowed = info['allowed_5']
            if allowed == {1, 2, 3, 4}:
                pred = lam_primes
            else:
                pred = lib.forbidden_residue(allowed)
        else:
            pred = 0.0
        diff = abs(pred - lam) if abs(pred) < 1e6 else float('inf')
        print(f"{name:>30} | {n:>8,} | {lam:>+11.4f} | {pred:>+11.4f} | "
              f"{diff:>8.4f} | {info['desc']:>20}")

    # Validation : Carmichael & sums of two squares
    print("\n" + "=" * 92)
    print("OTHER SEQUENCES (no a-priori prediction)")
    print("=" * 92)
    print(f"\n{'sequence':>30} | {'count':>8} | {'lambda obs':>11} | "
          f"{'comment':>40}")
    print("-" * 95)

    carms = gen_carmichael_small(N_max)
    lam, n = compute_lambda(carms, N_max=N_max)
    print(f"{'Carmichael numbers':>30} | {n:>8,} | {lam:>+11.4f} | "
          f"{'≥3 odd primes, complex structure':>40}")

    # Identification game
    print("\n" + "=" * 92)
    print("IDENTIFICATION GAME (extended library)")
    print("=" * 92)

    # Build extended library
    candidates = [("primes (no constraint)", lam_primes)]
    for forbidden in [1, 2, 3, 4]:
        allowed = {1, 2, 3, 4} - {forbidden}
        pred = lib.forbidden_residue(allowed)
        candidates.append((f"primes forbid {forbidden} mod 5", pred))
    for q in [7, 11, 13, 17, 19, 23, 29, 31]:
        candidates.append((f"biprimes spf={q}", lib.multiplicative((q,))))

    for name, seq_factory, info in test_cases[1:]:  # skip primes ref
        lam, n = compute_lambda(seq_factory(), N_max=N_max)
        ranked = lib.identify(lam, candidates)[:3]
        print(f"\n  {name} (n={n:,}, lam={lam:+.2f}) :")
        for i, (cname, cpred, cdist) in enumerate(ranked):
            marker = " <-- best" if i == 0 else ""
            print(f"    {i + 1}. {cname:<35} pred={cpred:+8.2f}  diff={cdist:6.2f}{marker}")


if __name__ == "__main__":
    main(N_max=10_000_000)
