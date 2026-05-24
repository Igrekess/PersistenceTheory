"""
Vérification indépendante des claims spectraux T_5 / T_30 dans PT.

Construit T_5, T_30 empiriques à partir des résidus de premiers consécutifs.
Calcule les spectres et compare aux trois claims :
  - ch21 ch37b : "toutes non-Perron de T_30 sur |z| = 1/4" (DER-NUM)
  - ch03 T2     : "|λ_2(T_30)| = 1/4"
  - Script PT_v7 : T_5 modèle uniforme → eigenvalues {1, -1/3, -1/3, -1/3}
  - Subagent    : T_5 empirique → eigenvalues à |z| ≈ 0.13-0.17
"""

import numpy as np
from sympy import sieve

# Générer les premiers
print("Generating primes...")
sieve.extend_to_no(200_000)   # cible ~10^6 premiers (~1.3M en fait)
primes = np.array(list(sieve._list))
print(f"  Got {len(primes):,} primes, last = {primes[-1]:,}")

# Couper sur premiers > 5 pour avoir des résidus dans (Z/30Z)*
primes = primes[primes > 5]
print(f"  After filter p > 5: {len(primes):,}")

def build_T(p, primes, n_max=None):
    """Build empirical transition matrix on (Z/pZ)* from consecutive primes."""
    if n_max is None:
        n_max = len(primes) - 1
    primes = primes[:n_max + 1]
    residues = primes % p
    # restrict to coprime residues
    mask = np.gcd(residues, p) == 1
    valid_pairs_idx = np.where(mask[:-1] & mask[1:])[0]
    r_from = residues[:-1][mask[:-1] & mask[1:]]
    r_to = residues[1:][mask[:-1] & mask[1:]]
    # Liste des résidus coprime avec p
    reduced = sorted([r for r in range(p) if np.gcd(r, p) == 1])
    n = len(reduced)
    idx = {r: i for i, r in enumerate(reduced)}
    T = np.zeros((n, n))
    for rf, rt in zip(r_from, r_to):
        T[idx[rf], idx[rt]] += 1
    # Normaliser par ligne (stochasticité)
    row_sums = T.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    T = T / row_sums
    return T, reduced, len(r_from)

# ============================================================
# 1. T_5 empirique
# ============================================================
print("\n=== T_5 empirique (résidus (Z/5Z)* = {1,2,3,4}) ===")
T5, res5, n5_pairs = build_T(5, primes)
print(f"  Pairs used: {n5_pairs:,}")
print(f"  Residues: {res5}")
print(f"  Matrix T_5:")
print(f"  {T5}")
print(f"  diag(T_5) mean: {np.mean(np.diag(T5)):.4f}")
print(f"  diag(T_5) values: {np.diag(T5)}")
ev5 = np.linalg.eigvals(T5)
ev5_sorted = sorted(ev5, key=lambda e: -abs(e))
print(f"  Eigenvalues T_5: {[f'{e.real:+.4f}{e.imag:+.4f}i' for e in ev5_sorted]}")
print(f"  Moduli: {[f'{abs(e):.4f}' for e in ev5_sorted]}")

# ============================================================
# 2. T_30 empirique
# ============================================================
print("\n=== T_30 empirique (résidus (Z/30Z)* = {1,7,11,13,17,19,23,29}) ===")
T30, res30, n30_pairs = build_T(30, primes)
print(f"  Pairs used: {n30_pairs:,}")
print(f"  Residues: {res30}")
ev30 = np.linalg.eigvals(T30)
ev30_sorted = sorted(ev30, key=lambda e: -abs(e))
print(f"  Eigenvalues T_30 (8):")
for e in ev30_sorted:
    print(f"    {e.real:+.5f}{e.imag:+.5f}i   |.| = {abs(e):.5f}")
print(f"  Diag mean: {np.mean(np.diag(T30)):.5f}")

# ============================================================
# 3. T_3 ⊗ T_5 vs T_30
# ============================================================
print("\n=== Factorisation CRT : T_3 ⊗ T_5 vs T_30 ===")
T3, res3, n3 = build_T(3, primes)
print(f"  T_3 = {T3}")
print(f"  Diag T_3: {np.diag(T3)} (T1 dit diag=0 ; empirique?)")
ev3 = np.linalg.eigvals(T3)
print(f"  Eigenvalues T_3: {[f'{e.real:+.4f}{e.imag:+.4f}i' for e in sorted(ev3, key=lambda e:-abs(e))]}")

# T_3 ⊗ T_5
T3T5 = np.kron(T3, T5)
print(f"  T_3 ⊗ T_5 shape: {T3T5.shape}")
ev_kron = np.linalg.eigvals(T3T5)
ev_kron_sorted = sorted(ev_kron, key=lambda e: -abs(e))
print(f"  Eigenvalues T_3 ⊗ T_5:")
for e in ev_kron_sorted:
    print(f"    {e.real:+.5f}{e.imag:+.5f}i   |.| = {abs(e):.5f}")

# Comparer matrices : need to reorder rows/columns of T30 to match Kronecker basis
# Kronecker basis on residues: (r3, r5) ↔ via CRT to r30 ∈ (Z/30Z)*
# T_3 indexed by {1,2}, T_5 indexed by {1,2,3,4}
# Kronecker basis: (r3, r5) -> 4*(r3-1) + (r5-1) when restricting to non-zero
# CRT: r30 ≡ r3 mod 3, r30 ≡ r5 mod 5

# Build the CRT permutation
crt_order = []
for r3 in [1, 2]:
    for r5 in [1, 2, 3, 4]:
        # Find r30 in res30 such that r30 % 3 == r3 and r30 % 5 == r5
        for r30 in res30:
            if r30 % 3 == r3 and r30 % 5 == r5:
                crt_order.append(res30.index(r30))
                break
T30_reorder = T30[np.ix_(crt_order, crt_order)]

diff = T30_reorder - T3T5
print(f"\n  || T_30 (réordonné CRT) - T_3 ⊗ T_5 ||_F = {np.linalg.norm(diff):.5f}")
print(f"  || T_30 (réordonné CRT) - T_3 ⊗ T_5 ||_inf (max entry) = {np.max(np.abs(diff)):.5f}")
print(f"  Relative Frobenius : {np.linalg.norm(diff)/np.linalg.norm(T30_reorder):.5f}")

# ============================================================
# 4. J-intertwiner antidiagonal mod 5 (r ↔ 5-r)
# ============================================================
print("\n=== J-intertwiner antidiagonal mod 5 (r ↔ 5-r) ===")
# J sur résidus {1,2,3,4}: 1↔4, 2↔3
J5 = np.array([[0, 0, 0, 1],
               [0, 0, 1, 0],
               [0, 1, 0, 0],
               [1, 0, 0, 0]], dtype=float)
commutator = T5 @ J5 - J5 @ T5
print(f"  ||T_5 J - J T_5||_F = {np.linalg.norm(commutator):.5f}")
print(f"  ||T_5 J - J T_5||_inf = {np.max(np.abs(commutator)):.5f}")
print(f"  Relative: {np.linalg.norm(commutator)/np.linalg.norm(T5):.5f}")

# Vérifier si T_5 est presque-symétrique sous J: J T J^T
JTJ = J5 @ T5 @ J5  # J5 self-inverse
diff_JT = JTJ - T5
print(f"  ||J T_5 J^T - T_5||_F = {np.linalg.norm(diff_JT):.5f}")
print(f"  → Si ~0 alors J commute avec T_5 (T_5 J-symétrique)")

# ============================================================
# 5. Comparaison avec PT_v7 uniform model T_5
# ============================================================
print("\n=== PT_v7 modèle uniforme T_5 vs empirique ===")
T5_uniform = (np.ones((4, 4)) - np.eye(4)) / 3
print(f"  T_5 uniforme = (1-I)/3, eigenvalues = {sorted([abs(e) for e in np.linalg.eigvals(T5_uniform)], reverse=True)}")
print(f"  ||T_5 empirique - T_5 uniforme||_F = {np.linalg.norm(T5 - T5_uniform):.5f}")
print(f"  Relative: {np.linalg.norm(T5 - T5_uniform)/np.linalg.norm(T5):.5f}")
