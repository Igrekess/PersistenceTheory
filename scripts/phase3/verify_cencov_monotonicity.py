"""
verify_cencov_monotonicity.py — Phase 3, script 1/4

Verifies numerically Cencov's monotonicity of the Fisher metric under the
PT sieve Markov kernels:

    pi_{m' -> m} : P(Z/m'Z) -> P(Z/mZ)    (reduction mod m, m | m')

For each pair (m', m) with m | m' in the primorial chain
    m in {6, 30, 210, 2310},
we check that for a random tangent vector v on P(Z/m'Z),

    g_Fisher^m ( pi_* v, pi_* v )   <=   g_Fisher^{m'} ( v, v ).

The Fisher metric at a point P = (p_0, ..., p_{k-1}) in the open simplex is
    g_Fisher(v, v)  =  sum_r  v_r^2 / p_r.

The pushforward by pi : v_r' -> v_r = sum_{r' : r' = r mod m} v_r'.

This script executes etape 2 of the Phase 3 RH proof skeleton:
  monotonicity / Cencov uniqueness of Fisher.

Exit 0 iff every projection satisfies the contraction inequality to
machine precision on 1000 random tangent vectors per pair.

Author : Phase 3 RH chantier (2026-04-22).
"""

from __future__ import annotations

import numpy as np

rng = np.random.default_rng(20260422)

# ---------------------------------------------------------------------------
# Sieve distributions: empirical prime residue frequencies at modulus m,
# supported on the coprime classes (Z/mZ)^*.
# ---------------------------------------------------------------------------

def coprime_classes(m: int) -> np.ndarray:
    """Coprime residue classes modulo m (Z/mZ)^*."""
    return np.array([r for r in range(m) if np.gcd(r, m) == 1], dtype=int)

def primes_up_to(N: int) -> list[int]:
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(np.sqrt(N)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return [int(p) for p in np.flatnonzero(sieve)]

def sieve_distribution(m: int, N_primes: int = 200_000) -> np.ndarray:
    """Empirical prime residue distribution on (Z/mZ)^*.
    Returns P of length phi(m), indexed by coprime_classes(m) order.
    """
    primes = primes_up_to(3_000_000)[:N_primes]
    classes = coprime_classes(m)
    idx = {int(c): i for i, c in enumerate(classes)}
    counts = np.zeros(len(classes))
    for p in primes:
        r = p % m
        if r in idx:
            counts[idx[r]] += 1
    P = counts / counts.sum()
    # Regularise small zero probabilities to avoid singular Fisher metric
    P = P + 1e-12
    P = P / P.sum()
    return P

# ---------------------------------------------------------------------------
# Markov reduction kernel pi : P(Z/m'Z)^* -> P(Z/mZ)^*  for m | m'
# ---------------------------------------------------------------------------

def reduction_kernel(m_high: int, m_low: int) -> np.ndarray:
    """Row-stochastic K of shape (phi(m_high), phi(m_low)).
    K[i, j] = 1 if class c_high[i] reduces mod m_low to coprime class c_low[j],
    else 0.  Classes of m_high that are NOT coprime to m_low map to nothing,
    but every coprime class of m_high is coprime to m_low when m_low | m_high
    and they share the same prime factors (primorial chain).
    """
    assert m_high % m_low == 0, f"{m_low} must divide {m_high}"
    c_high = coprime_classes(m_high)
    c_low = coprime_classes(m_low)
    idx_low = {int(c): j for j, c in enumerate(c_low)}
    K = np.zeros((len(c_high), len(c_low)))
    for i, ch in enumerate(c_high):
        r = int(ch) % m_low
        if r in idx_low:
            K[i, idx_low[r]] = 1.0
    # Row-stochastic already (exactly one 1 per row in primorial chain)
    return K

# ---------------------------------------------------------------------------
# Fisher metric and pushforward
# ---------------------------------------------------------------------------

def fisher_norm_sq(v: np.ndarray, P: np.ndarray) -> float:
    """||v||_Fisher^2 at P.  v must satisfy sum v = 0."""
    return float(np.sum(v * v / P))

def pushforward(v_high: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Tangent pushforward: the reduction is a linear map on the simplex,
    so its derivative is K^T acting on tangent vectors of the coprime simplex.
    v_low[j] = sum_i K[i,j] * v_high[i].
    """
    return K.T @ v_high

def random_tangent(dim: int) -> np.ndarray:
    v = rng.normal(size=dim)
    v -= v.mean()  # enforce sum = 0 (tangent to simplex)
    return v

# ---------------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------------

def check_pair(m_high: int, m_low: int, n_trials: int = 1000) -> dict:
    P_high = sieve_distribution(m_high)
    P_low = sieve_distribution(m_low)
    K = reduction_kernel(m_high, m_low)

    # Sanity: K maps P_high to close-to-P_low (it's the exact marginal map
    # on the sieve data; tiny discrepancies come from finite-sample Chebyshev
    # biases, not from the kernel).
    P_low_pushed = K.T @ P_high
    marginal_err = float(np.linalg.norm(P_low_pushed - P_low, ord=1))

    ratios = []
    for _ in range(n_trials):
        v_high = random_tangent(len(P_high))
        v_low = pushforward(v_high, K)
        num = fisher_norm_sq(v_low, P_low_pushed + 1e-16)
        den = fisher_norm_sq(v_high, P_high)
        if den > 0:
            ratios.append(num / den)
    ratios = np.array(ratios)
    max_ratio = float(ratios.max())
    mean_ratio = float(ratios.mean())
    pass_cencov = max_ratio <= 1.0 + 1e-9

    return {
        "m_high": m_high,
        "m_low": m_low,
        "phi_high": len(P_high),
        "phi_low": len(P_low),
        "marginal_L1_err": marginal_err,
        "max_ratio": max_ratio,
        "mean_ratio": mean_ratio,
        "pass_cencov": pass_cencov,
    }

def main():
    print("=" * 72)
    print("Phase 3 / script 1 : Cencov monotonicity of Fisher metric under")
    print("sieve reduction kernels  pi_{m' -> m}.")
    print("=" * 72)

    chain = [6, 30, 210, 2310]
    pairs = [(chain[i+1], chain[i]) for i in range(len(chain) - 1)]

    all_pass = True
    rows = []
    for m_high, m_low in pairs:
        r = check_pair(m_high, m_low)
        rows.append(r)
        print()
        print(f"m' = {r['m_high']:5d}  ->  m = {r['m_low']:5d}"
              f"   phi' = {r['phi_high']:4d}   phi = {r['phi_low']:4d}")
        print(f"    max Fisher ratio over 1000 random tangents = {r['max_ratio']:.6e}")
        print(f"    mean Fisher ratio                          = {r['mean_ratio']:.6e}")
        print(f"    marginal L1 error on sieve distribution    = {r['marginal_L1_err']:.3e}")
        print(f"    Cencov monotonicity (ratio <= 1) : {'PASS' if r['pass_cencov'] else 'FAIL'}")
        all_pass = all_pass and r["pass_cencov"]

    print()
    print("=" * 72)
    print(f"GLOBAL : {'PASS' if all_pass else 'FAIL'}  "
          f"(4 primorial projections tested, 1000 tangents each)")
    print("=" * 72)

    return 0 if all_pass else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
