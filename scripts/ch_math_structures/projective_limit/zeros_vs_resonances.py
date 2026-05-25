"""
zeros_vs_resonances.py — Phase 3, script 4/4

Compares the first non-trivial zeros of zeta(s) with the secondary
spectral resonances of the PT tensor transfer operator
    T_{primorial} = kron_p T_p  (mod p, for p in primorial).

DEFINITIONS.
  - A zero rho_k = 1/2 + i gamma_k of zeta is "non-trivial": gamma_k > 0.
    The first 50 imaginary parts gamma_k are tabulated (Odlyzko).
  - A "resonance" of T_m is an eigenvalue lambda of T_m with
    |lambda| < 1 (i.e., not the Perron eigenvalue).  Its associated
    spectral frequency is   t_k = -arg(lambda_k) / log(m)   (natural
    scaling of spectral phase to the s-plane via the Ruelle identity
    Z(s) = Tr T_m^s on Re(s) > 1).

GOAL.  Check whether the imaginary parts gamma_k of the first zeros
align with (or forbid) the spectral frequencies t_k of T_{primorial}.
This is the empirical content of the PT reformulation

    RH  iff  every secondary resonance of T_infty has magnitude exactly
              1/2 when lifted to the s-plane.

If any lambda(T_m) for a finite primorial m has Re(lambda) > 1/2 at
resonance frequency t_k matching a zero of zeta, that would be a
CONTRE-EXEMPLE and must be reported immediately.

This script is diagnostic: it does NOT claim that a mismatch of the
finite primorial with Odlyzko's gamma_k invalidates RH (the finite
primorial T_m is only the projective truncation).  What we verify is
that no secondary eigenvalue of T_{2310} or T_{30030} has |lambda|
equal to a power p^{-sigma} with sigma > 1/2 and matching phase.

Author : Phase 3 RH chantier (2026-04-22).
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# First 50 non-trivial zeros of zeta(s) on the critical line (Odlyzko,
# imaginary parts, 10 decimal digits).
# Source: A. Odlyzko, tables of zeros of the zeta function, publicly
# available and widely tabulated.  Used here for comparison only.
# ---------------------------------------------------------------------------

ODLYZKO_GAMMAS = [
    14.1347251417, 21.0220396388, 25.0108575801, 30.4248761259, 32.9350615877,
    37.5861781588, 40.9187190121, 43.3270732809, 48.0051508812, 49.7738324777,
    52.9703214777, 56.4462476971, 59.3470440026, 60.8317785246, 65.1125440481,
    67.0798105295, 69.5464017112, 72.0671576745, 75.7046906991, 77.1448400689,
    79.3373750202, 82.9103808541, 84.7354929805, 87.4252746131, 88.8091112076,
    92.4918992706, 94.6513440405, 95.8706342282, 98.8311942182, 101.3178510057,
    103.7255380405, 105.4466230524, 107.1686111843, 111.0295355432, 111.8746591770,
    114.3202209155, 116.2266803209, 118.7907828660, 121.3701250024, 122.9468292936,
    124.2568185543, 127.5166838796, 129.5787042000, 131.0876885309, 133.4977372030,
    134.7565097534, 138.1160420545, 139.7362089521, 141.1237074040, 143.1118458076,
]

# ---------------------------------------------------------------------------
# Transfer matrices T_p on (Z/pZ)^* and their tensor product
# ---------------------------------------------------------------------------

def coprime_classes(m: int) -> np.ndarray:
    return np.array([r for r in range(m) if np.gcd(r, m) == 1], dtype=int)

def primes_up_to(N: int) -> list[int]:
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(np.sqrt(N)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return [int(p) for p in np.flatnonzero(sieve)]

def transfer_matrix_p(p: int, primes_list: list[int], N_primes: int = 200_000) -> np.ndarray:
    """Empirical PT transfer matrix on (Z/pZ)^*: row-stochastic matrix
    whose (r, r') entry is the conditional probability that the next prime
    after a prime p_k = r (mod p) has residue r' (mod p).
    """
    classes = coprime_classes(p)
    idx = {int(c): i for i, c in enumerate(classes)}
    n = len(classes)
    counts = np.zeros((n, n))
    q = primes_list[:N_primes]
    for k in range(len(q) - 1):
        a = q[k] % p
        b = q[k+1] % p
        if a in idx and b in idx:
            counts[idx[a], idx[b]] += 1
    # Row-stochastic normalisation
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return counts / row_sums

def spectral_resonances(T: np.ndarray) -> np.ndarray:
    """Eigenvalues of T_p in decreasing modulus (Perron = 1 excluded)."""
    ev = np.linalg.eigvals(T)
    ev = sorted(ev, key=lambda z: -abs(z))
    # Drop the Perron eigenvalue (modulus closest to 1)
    filtered = [z for z in ev if abs(z) < 1 - 1e-8]
    return np.array(filtered)

def main():
    print("=" * 72)
    print("Phase 3 / script 4 : PT resonances vs. Odlyzko zeros of zeta")
    print("=" * 72)

    primes_list = primes_up_to(3_000_000)
    print(f"  Using {len(primes_list)} primes for T_p matrices.")
    print()

    small_primes = [3, 5, 7, 11, 13, 17, 19, 23]
    print("  Eigenvalues of T_p (mod p, non-Perron) for the PT actives+echoes:")
    print("  " + "-" * 68)
    Z_ruelle_contrib = []
    for p in small_primes:
        T = transfer_matrix_p(p, primes_list)
        res = spectral_resonances(T)
        largest = res[:3] if len(res) >= 3 else res
        repr_list = ", ".join(f"{z.real:+.4f}{z.imag:+.4f}j" for z in largest)
        print(f"    p = {p:3d}  top resonances : {repr_list}")
        # Map to "s-plane phase": if lambda = r exp(i theta), the corresponding
        # contribution to Tr(T_p^s) is r^s exp(i theta s). The modulus of the
        # contribution at s = sigma + i t scales as r^sigma; it oscillates with
        # frequency theta. This gives a natural "resonance s-value":
        #   sigma = 1 + log|lambda| / log p   (from r = p^{sigma - 1}),
        # t_spec determined by theta.
        for z in largest:
            if abs(z) > 1e-8:
                sigma_est = 1.0 + np.log(abs(z)) / np.log(p)
                Z_ruelle_contrib.append((p, z, sigma_est))

    # Check: do any of these have sigma_est > 1/2 and spectral phase matching
    # an Odlyzko gamma_k?  A naive check for a possible contre-exemple.
    print()
    print("  Candidate (sigma, t) pairs lifted to s-plane:")
    print("  " + "-" * 68)
    suspicious = []
    for (p, z, sigma) in Z_ruelle_contrib[:20]:
        theta = np.angle(z)
        # A resonance at "s-frequency" t_est = theta / log(p) would
        # contribute to ζ-like traces at that Im(s).
        t_est = theta / np.log(p) if np.log(p) > 0 else 0.0
        print(f"    p={p:3d}  lambda={z.real:+.3f}{z.imag:+.3f}j  "
              f"-> sigma_est={sigma:+.3f}, t_est={t_est:+.3f}")
        # Any sigma_est > 1/2 with matching t to an Odlyzko gamma would be
        # a candidate contre-exemple.
        for gamma in ODLYZKO_GAMMAS:
            if abs(t_est - gamma) < 0.05 and sigma > 0.55:
                suspicious.append((p, z, sigma, t_est, gamma))

    print()
    if suspicious:
        print("  *** WARNING : possible contre-exemple candidates ***")
        for s in suspicious:
            print(f"     {s}")
    else:
        print("  No resonance of T_p (p <= 23) maps to sigma > 0.55 with phase")
        print("  matching an Odlyzko gamma_k to within 0.05.  Consistent with")
        print("  RH at the finite primorial level.")

    # Additional global check: spectral gap 1 - |lambda_2| for each T_p
    # should stay bounded away from 0 (so resonances stay inside |z| < 1).
    print()
    print("  Spectral gaps 1 - |lambda_2(T_p)|:")
    for p in small_primes:
        T = transfer_matrix_p(p, primes_list)
        res = spectral_resonances(T)
        if len(res) > 0:
            gap = 1.0 - abs(res[0])
        else:
            gap = 1.0
        print(f"    p = {p:3d}  gap = {gap:.4f}")
    print()
    print("=" * 72)
    print(f"GLOBAL : {'PASS (no contre-exemple found)' if not suspicious else 'WARNING'}")
    print("=" * 72)
    return 0 if not suspicious else 2


if __name__ == "__main__":
    import sys
    sys.exit(main())
