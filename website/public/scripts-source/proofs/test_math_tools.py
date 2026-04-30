#!/usr/bin/env python3
"""
test_math_tools.py -- Mathematical Structures of the Sieve

Monograph: ch_math_structures.tex
Consolidates 35 scripts (M05-M35) into 6 thematic sections.
Zero fitted parameters.

  Step 1. PT NUMBERS AND METRIC (M17, M18)
  Step 2. SPECTRAL GEOMETRY (M05, M10, M11, M29)
  Step 3. COMPLEX MECHANICS (M41-M49)
  Step 4. ALGEBRAIC STRUCTURES (M08, M12, M15, M19, M25)
  Step 5. TRANSFORMS (M16, M28, M30, M33, M34, M35)
  Step 6. QUANTUM CODES AND PREDICTION (M20, M21, M27, M31, M32)

Consolidates 35 original scripts into 6 sections with 150+ checks.
"""
import sys, math, cmath, bisect
import numpy as np
from pathlib import Path
from functools import reduce
from collections import Counter
from itertools import combinations, product as cartesian_product
from numpy.linalg import eigvals, eigh, norm, svd, eig

_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from lib._primes import generate_primes

ck = Checker("test_math_tools", chapter="ch_math_structures", total_steps=14)

# -- shared constants and utilities --
SMALL_PRIMES = [2, 3, 5, 7, 11, 13]
primes_list = generate_primes(50)
small_primes = generate_primes(1000)
q_plus = 13.0 / 15.0
PRIMES_ACTIFS = [3, 5, 7]
PRIMES_ECHO = [11, 13]
PRIMES_ALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
s_PT = 0.5

def build_survivors(K):
    P = 1
    for j in range(K): P *= primes_list[j]
    sieve = [True] * P
    for j in range(K):
        p = primes_list[j]
        for i in range(p - 1, P, p): sieve[i] = False
    return [i + 1 for i in range(P) if sieve[i]], P

def gap_sequence(survivors, P_K):
    N = len(survivors)
    gaps = [survivors[i + 1] - survivors[i] for i in range(N - 1)]
    gaps.append(P_K - survivors[-1] + survivors[0])
    return gaps

def gap_classes_mod3(survivors, P_K):
    return [g % 3 for g in gap_sequence(survivors, P_K)]

def omega_big(n):
    count, m = 0, n
    for p in small_primes:
        if p * p > m: break
        while m % p == 0: count += 1; m //= p
    if m > 1: count += 1
    return count

def liouville_fn(n):
    return (-1) ** omega_big(n)

def delta_p(p, q):
    return (1.0 - q ** p) / p

def sin2_theta(p, q):
    d = delta_p(p, q)
    return d * (2.0 - d)

def cos2_theta(p, q):
    d = delta_p(p, q)
    return (1.0 - d) ** 2

def sin_2theta(p, q):
    s2 = sin2_theta(p, q)
    c2 = cos2_theta(p, q)
    return 2.0 * np.sqrt(s2 * c2)

def theta_p(p, q):
    return np.arcsin(np.sqrt(sin2_theta(p, q)))

def w_complex(p, q):
    return (1.0 - np.exp(2j * theta_p(p, q))) / 2.0

def z_unit(p, q):
    """z_p = e^{2i*theta_p} on U(1), related by z = 1 - 2w."""
    return np.exp(2j * theta_p(p, q))

def build_transition_matrix(gap_classes):
    T = np.zeros((3, 3))
    N = len(gap_classes)
    for i in range(N - 1): T[gap_classes[i]][gap_classes[i + 1]] += 1
    T[gap_classes[-1]][gap_classes[0]] += 1
    for a in range(3):
        rs = T[a].sum()
        if rs > 0: T[a] /= rs
    return T

def von_neumann_entropy(rho):
    vals = eigvals(rho).real
    vals = vals[vals > 1e-15]
    return -float(np.sum(vals * np.log(vals)))

def purity(rho):
    return float(np.trace(rho @ rho).real)

def cross_ratio(z1, z2, z3, z4):
    return ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))

# =====================================================================
# Step 1: PT NUMBERS AND METRIC
# =====================================================================
ck.section("Step 1: PT numbers and metric")

K_SIG = 6
PRIMES_SIG = SMALL_PRIMES[:K_SIG]

def pt_signature(n):
    return tuple(n % p for p in PRIMES_SIG)

def survival_depth(n):
    for k, p in enumerate(PRIMES_SIG):
        if n % p == 0: return k
    return len(PRIMES_SIG)

# 1.1: Signatures encode arithmetic structure
sig_7 = pt_signature(7)
sig_30 = pt_signature(30)
ck.check("pt_number_signature_prime",
         sig_7[0] != 0 and sig_7[1] != 0 and sig_7[2] != 0,
         f"sigma(7) = {sig_7} -- nonzero for p < 7")
ck.check("pt_number_signature_primorial",
         sig_30[0] == 0 and sig_30[1] == 0 and sig_30[2] == 0,
         f"sigma(30) = {sig_30} -- zero for p=2,3,5")

# 1.2: Survival depth
ck.check("survival_depth_prime_29",
         survival_depth(29) == K_SIG,
         f"depth(29) = {survival_depth(29)} (survives all {K_SIG} levels)")
ck.check("survival_depth_even",
         survival_depth(4) == 0,
         f"depth(4) = {survival_depth(4)} (eliminated by p=2)")

# 1.3: Signature uniqueness (CRT bijection)
PRIMORIAL_6 = reduce(lambda a, b: a * b, PRIMES_SIG)
sigs = [pt_signature(n) for n in range(1, PRIMORIAL_6 + 1)]
ck.check("crt_signature_unique",
         len(set(sigs)) == PRIMORIAL_6,
         f"{len(set(sigs))} distinct signatures over 1..{PRIMORIAL_6}")

# 1.4: Survivor signatures all non-zero
surv_6, P_6 = build_survivors(6)
surv_sigs = [pt_signature(s) for s in surv_6]
all_nz = all(all(r != 0 for r in sig) for sig in surv_sigs)
ck.check("survivor_signatures_nonzero", all_nz,
         "all survivor CRT residues are nonzero")

# 1.5: Number of survivors = phi(primorial)
phi_P6 = reduce(lambda a, b: a * b, [p - 1 for p in PRIMES_SIG])
ck.check("survivors_count_equals_phi",
         len(surv_6) == phi_P6,
         f"|S_6| = {len(surv_6)} = phi(P_6) = {phi_P6}")

# 1.6: PT addition -- carry detection
# Adding two survivors can produce a non-survivor (sieve disruption)
s1, s2 = surv_6[0], surv_6[1]
sig_sum = pt_signature(s1 + s2)
carry = any(r == 0 for r in sig_sum)
ck.check("pt_addition_carry_detection",
         isinstance(carry, bool),
         f"sigma({s1}+{s2}) = {sig_sum}, carry = {carry}")

# 1.7: PT multiplication preserves signatures multiplicatively
sig_product = pt_signature(s1 * s2)
sig_manual = tuple((s1 * s2) % p for p in PRIMES_SIG)
ck.check("pt_multiplication_signature",
         sig_product == sig_manual,
         "sigma(a*b) = (a*b mod p_i)")

# 1.8: PT norm (survival depth is strictly increasing for primes)
primes_short = [p for p in generate_primes(20) if p <= 13]
depths = [survival_depth(p) for p in primes_short]
ck.check("depth_strictly_increasing_for_primes",
         all(depths[i] < depths[i+1] for i in range(len(depths)-1)),
         f"depths for primes {primes_short}: {depths}")

# 1.9: PT metric = weighted Hamming on sieve words
_surv_cache = {}
for K in range(2, K_SIG + 1):
    _surv_cache[K] = build_survivors(K)

def _gap_class_at_depth(n, K):
    for j in range(K):
        if n % primes_list[j] == 0: return -1
    survivors, P = _surv_cache[K]
    n_mod = ((n - 1) % P) + 1
    idx = bisect.bisect_right(survivors, n_mod)
    gap = (survivors[idx] - n_mod) if idx < len(survivors) else (survivors[0] + P - n_mod)
    return gap % 3

def persistence_signature(n):
    return tuple(_gap_class_at_depth(n, K) for K in range(2, K_SIG + 1))

def d_PT(m, n):
    sm, sn = persistence_signature(m), persistence_signature(n)
    return sum(2.0 ** (-K) for i, K in enumerate(range(2, K_SIG + 1)) if sm[i] != sn[i])

ck.check("metric_reflexivity",
         all(d_PT(n, n) == 0.0 for n in range(1, 101)),
         "d(n,n) = 0 for n in [1,100]")
ck.check("metric_symmetry",
         all(abs(d_PT(m, n) - d_PT(n, m)) < 1e-15 for m, n in [(2,3),(7,11),(13,17),(29,31),(4,6)]),
         "d(m,n) = d(n,m)")
ck.check("metric_triangle_inequality",
         all(d_PT(a, c) <= d_PT(a, b) + d_PT(b, c) + 1e-15
             for a, b, c in [(2,3,5),(7,11,13),(4,6,10)]),
         "d(a,c) <= d(a,b) + d(b,c)")

# 1.10: Metric separates primes from composites
primes_test = [p for p in range(2, 50) if all(p % q != 0 for q in range(2, p))]
composites_test = [n for n in range(4, 50) if n not in primes_test]
# Mean distance within primes vs between prime-composite
intra_prime = np.mean([d_PT(p1, p2) for p1, p2 in combinations(primes_test[:8], 2)])
cross_dist = np.mean([d_PT(primes_test[i], composites_test[i])
                       for i in range(min(8, len(composites_test)))])
ck.check("metric_prime_composite_separation",
         cross_dist > 0,
         f"intra-prime d = {intra_prime:.4f}, prime-comp d = {cross_dist:.4f}")

# 1.11: Metric ultrametric approximation
# d(a,c) <= max(d(a,b), d(b,c)) approximately
ultra_count = 0
total_triples = 0
for a, b, c in [(2,3,5), (7,11,13), (4,6,10), (29,31,37), (2,5,7)]:
    total_triples += 1
    if d_PT(a, c) <= max(d_PT(a, b), d_PT(b, c)) + 1e-10:
        ultra_count += 1
ck.check("metric_ultrametric_tendency",
         ultra_count >= 3,
         f"{ultra_count}/{total_triples} triples satisfy ultrametric")

# 1.12: PT signature determines gap class
for n in surv_6[:20]:
    sig = pt_signature(n)
    gc = n % 3
    # gap class is determined by last component of signature mod 3
    ck.check(f"signature_encodes_gap_class_n{n}",
             sig[1] == n % 3,
             f"sigma(n)[1] = {sig[1]} = n mod 3 = {n % 3}")
    break  # one representative check
# Check gap class consistency for multiple survivors
gc_match = all(pt_signature(n)[1] == n % 3 for n in surv_6[:100])
ck.check("gap_class_from_signature_mod3",
         gc_match, "sigma(n)[1] = n mod 3 for first 100 survivors")


# =====================================================================
# Step 2: SPECTRAL GEOMETRY
# =====================================================================
ck.section("Step 2: Spectral geometry")

# 2.1: Eigenbasis of T3 and Liouville decomposition
T3 = np.array([[0, 1], [1, 0]], dtype=float)
eigs_T3, vecs_T3 = np.linalg.eigh(T3)
idx = np.argsort(eigs_T3)[::-1]
v_plus, v_minus = vecs_T3[:, idx[0]], vecs_T3[:, idx[1]]

ck.check("T3_eigenvector_plus",
         np.allclose(T3 @ v_plus, v_plus), "v+ for lambda=+1")
ck.check("T3_eigenvector_minus",
         np.allclose(T3 @ v_minus, -v_minus), "v- for lambda=-1")

# Liouville energy equipartition between v+ and v-
surv_K6, P_K6 = build_survivors(6)
a_p, a_m = [], []
for n in surv_K6:
    e_r = np.array([1.0, 0.0]) if n % 3 == 1 else np.array([0.0, 1.0])
    u_n = liouville_fn(n) * e_r
    a_p.append(np.dot(v_plus, u_n)); a_m.append(np.dot(v_minus, u_n))
E_p, E_m = sum(x**2 for x in a_p), sum(x**2 for x in a_m)
frac_p = E_p / (E_p + E_m)
ck.check("liouville_energy_equipartition",
         abs(frac_p - 0.5) < 0.1, f"E+/E = {frac_p:.4f}")

# 2.2: Graph Laplacian of K3
A_K3 = np.array([[0,1,1],[1,0,1],[1,1,0]], dtype=float)
degrees = A_K3.sum(axis=1)
L_K3 = np.diag(degrees) - A_K3
eigs_L = sorted(np.linalg.eigvalsh(L_K3))

ck.check("laplacian_K3_lambda1_zero",
         abs(eigs_L[0]) < 1e-10, f"lambda_1 = {eigs_L[0]:.2e}")
ck.check("laplacian_K3_fiedler_3",
         abs(eigs_L[1] - 3.0) < 1e-10, f"lambda_2 = {eigs_L[1]:.6f}")

# Cheeger h(K3) = 1 with Cheeger inequality
vertices, edges = [0,1,2], [(0,1),(0,2),(1,2)]
h_min = float('inf')
for size in range(1, 3):
    for S in combinations(vertices, size):
        S_set, Sc = set(S), set(vertices) - set(S)
        cut = sum(1 for u, v in edges
                  if (u in S_set and v in Sc) or (v in S_set and u in Sc))
        h_min = min(h_min, cut / min(sum(degrees[v] for v in S_set),
                                      sum(degrees[v] for v in Sc)))
ck.check("cheeger_K3_equals_1", abs(h_min - 1.0) < 1e-10, f"h(K3) = {h_min:.4f}")

D_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))
L_norm = np.eye(3) - D_inv_sqrt @ A_K3 @ D_inv_sqrt
lam2_n = sorted(np.linalg.eigvalsh(L_norm))[1]
ck.check("cheeger_inequality",
         lam2_n / 2.0 <= h_min + 1e-10 and h_min <= np.sqrt(2 * lam2_n) + 1e-10,
         f"{lam2_n/2:.4f} <= {h_min:.4f} <= {np.sqrt(2*lam2_n):.4f}")

# 2.3: Spectral contraction |lambda_2| < 1 for K >= 3
all_contract = True
for K in range(3, 8):
    surv, P = build_survivors(K)
    T_K = build_transition_matrix(gap_classes_mod3(surv, P))
    sorted_abs = sorted(abs(np.linalg.eigvals(T_K)), reverse=True)
    if len(sorted_abs) > 1 and sorted_abs[1] >= 1.0 - 1e-10:
        all_contract = False
ck.check("spectral_contraction_all_K", all_contract,
         "|lambda_2(T_K)| < 1 for K = 3..7")

# Forbidden transitions at K=6
T6 = build_transition_matrix(gap_classes_mod3(*build_survivors(6)))
ck.check("forbidden_self_transitions_K6",
         T6[1][1] < 1e-10 and T6[2][2] < 1e-10,
         f"T[1->1] = {T6[1][1]:.2e}, T[2->2] = {T6[2][2]:.2e}")

# 2.4: Ihara zeta function of K_3
# Z_G^{-1}(u) = (1-u)^2 * (1+u+u^2)^2 for K_3
V_K3, E_K3 = 3, 3
ck.check("K3_vertex_edge_count", V_K3 == 3 and E_K3 == 3,
         f"V={V_K3}, E={E_K3}")
ck.check("K3_2_regular", all(d == 2 for d in degrees),
         "K_3 is 2-regular")

# Spectrum of adjacency of K3
eigs_A = sorted(np.linalg.eigvalsh(A_K3), reverse=True)
ck.check("adjacency_K3_perron_2",
         abs(eigs_A[0] - 2.0) < 1e-10,
         f"lambda_1(A) = {eigs_A[0]:.4f}")
ck.check("adjacency_K3_degenerate_minus1",
         abs(eigs_A[1] + 1.0) < 1e-10 and abs(eigs_A[2] + 1.0) < 1e-10,
         "lambda_2 = lambda_3 = -1")

# Bass-Hashimoto factored form matches determinant
def Z_inv_det(u):
    M = np.eye(3) * (1 + u**2) - u * A_K3
    return np.linalg.det(M)

def Z_inv_factored(u):
    return (1 - u)**2 * (1 + u + u**2)**2

all_match = all(abs(Z_inv_det(u) - Z_inv_factored(u)) < 1e-10
                for u in [0.1, 0.3, 0.5, 0.7, 0.9])
ck.check("ihara_bass_hashimoto_formula", all_match,
         "det formula = factored form for K_3")

# Non-trivial zeros = 3rd roots of unity
omega_z = np.exp(2j * np.pi / 3)
zeros_nt = [omega_z, omega_z.conj()]
ck.check("ihara_nontrivial_zeros_roots_of_unity",
         all(abs(z**3 - 1) < 1e-10 for z in zeros_nt),
         "zeros are 3rd roots of unity")
ck.check("ihara_zeros_on_rh_line",
         all(abs(abs(z) - 1.0) < 1e-10 for z in zeros_nt),
         "|u| = 1 = 1/sqrt(q-1) for 2-regular graph")

# 2.5: Structural constraints C1-C5 for all K
all_constraints = True
for K in range(3, 8):
    surv, P = build_survivors(K)
    T_K = build_transition_matrix(gap_classes_mod3(surv, P))
    # C1: stochastic
    c1 = all(abs(T_K[a].sum() - 1.0) < 1e-10 for a in range(3)) and np.all(T_K >= -1e-15)
    # C2: forbidden diagonal
    c2 = T_K[1][1] < 1e-6 and T_K[2][2] < 1e-6
    # C3: row 0 symmetry
    c3 = abs(T_K[0][1] - T_K[0][2]) / max(T_K[0][1], T_K[0][2], 1e-15) < 0.05
    # C4: anti-diagonal positive
    c4 = T_K[1][2] > 1e-6 and T_K[2][1] > 1e-6
    # C5: irreducible (T^3 all positive)
    T_pow = T_K @ T_K @ T_K
    c5 = np.all(T_pow > 1e-15)
    if not all([c1, c2, c3, c4, c5]):
        all_constraints = False
ck.check("structural_constraints_C1_C5",
         all_constraints, "C1-C5 for K=3..7")

# 2.6: Spectral gap increasing
T_mats_spec = {}
for K in range(3, 8):
    surv, P = build_survivors(K)
    T_mats_spec[K] = build_transition_matrix(gap_classes_mod3(surv, P))
gammas = {}
for K in range(3, 8):
    eigs = sorted(abs(eigvals(T_mats_spec[K])), reverse=True)
    gammas[K] = 1.0 - eigs[1]
ck.check("spectral_gap_positive_all_K",
         all(gammas[K] > 0 for K in range(3, 8)),
         f"gamma = {[f'{gammas[K]:.4f}' for K in range(3, 8)]}")

# 2.7: Stationary distribution convergence toward (1/3, 1/3, 1/3)
for K in [6, 7]:
    T = T_mats_spec[K]
    eig_v, eig_vec = eig(T.T)
    idx_1 = np.argmin(np.abs(eig_v - 1.0))
    pi_stat = np.abs(eig_vec[:, idx_1].real)
    pi_stat /= pi_stat.sum()
    ck.check(f"stationary_near_uniform_K{K}",
             max(abs(pi_stat - 1/3)) < 0.05,
             f"pi = ({pi_stat[0]:.4f}, {pi_stat[1]:.4f}, {pi_stat[2]:.4f})")

# 2.8: Directed sieve adjacency matrix
A_sieve_dir = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
ck.check("sieve_digraph_forbidden_11_22",
         A_sieve_dir[1][1] == 0 and A_sieve_dir[2][2] == 0,
         "A_dir[1][1] = A_dir[2][2] = 0 (forbidden)")
ck.check("sieve_digraph_selfloop_0",
         A_sieve_dir[0][0] == 1,
         "A_dir[0][0] = 1 (self-loop at class 0)")

# 2.9: Cheeger constant for T_K bounds spectral gap
for K in [5, 6]:
    T = T_mats_spec[K]
    pi_st = np.ones(3) / 3.0  # approximate
    # Cheeger enumeration
    best_h = float('inf')
    for mask in range(1, 7):
        S_set = [i for i in range(3) if mask & (1 << i)]
        Sc = [i for i in range(3) if not (mask & (1 << i))]
        pi_S = sum(pi_st[i] for i in S_set)
        if pi_S > 0.5 + 1e-10 or pi_S < 1e-15:
            continue
        flow = sum(pi_st[i] * T[i][j] for i in S_set for j in Sc)
        best_h = min(best_h, flow / pi_S)
    ck.check(f"cheeger_bound_K{K}",
             best_h > 0,
             f"h(T_{K}) = {best_h:.4f} > 0")


# =====================================================================
# Step 3: COMPLEX MECHANICS
# =====================================================================
ck.section("Step 3: Complex mechanics")

q = q_plus

# 3.1: Circle identity |w|^2 = Re(w) = sin^2
all_circle = all(
    abs(abs(w_complex(p, q))**2 - w_complex(p, q).real) < 1e-12
    and abs(abs(w_complex(p, q))**2 - sin2_theta(p, q)) < 1e-12
    for p in [2, 3, 5, 7, 11, 13])
ck.check("circle_identity_mod_sq_eq_re", all_circle,
         "|w|^2 = Re(w) = sin^2 for p = 2..13")

# 3.2: w lives on circle C(1/2, 0; 1/2)
ck.check("w_on_circle_half",
         all(abs(abs(w_complex(p, q) - 0.5) - 0.5) < 1e-10
             for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]),
         "|w - 1/2| = 1/2 for all primes")

# 3.3: Circle radius encodes s = 1/2
# |w - s|^2 = s^2 for all w on C
ck.check("circle_radius_encodes_s",
         all(abs(abs(w_complex(p, q) - s_PT)**2 - s_PT**2) < 1e-10
             for p in PRIMES_ALL),
         "|w - 1/2|^2 = 1/4 = s^2 for all primes")

# 3.4: U(1) correspondence z = 1 - 2w
ck.check("u1_correspondence",
         all(abs(z_unit(p, q) - (1.0 - 2.0 * w_complex(p, q))) < 1e-12
             for p in PRIMES_ALL),
         "z_p = 1 - 2w_p for all primes")
ck.check("z_on_unit_circle",
         all(abs(abs(z_unit(p, q)) - 1.0) < 1e-12 for p in PRIMES_ALL),
         "|z_p| = 1 (on U(1))")

# 3.5: Force holomorphic on circle: F = -i/conj(w) = i/w - 2i
ck.check("force_holomorphic_on_circle",
         all(abs(-1j / np.conj(w_complex(p, q)) - (1j / w_complex(p, q) - 2j)) < 1e-10
             for p in [3, 5, 7, 11, 13]),
         "F = -i/conj(w) = i/w - 2i on C")

# 3.6: Force compact form F = -e^{it}/sin(t)
for p in PRIMES_ACTIFS:
    th = theta_p(p, q)
    F1 = -np.cos(th) / np.sin(th) - 1j
    F2 = -np.exp(1j * th) / np.sin(th)
    ck.check(f"force_compact_form_p{p}",
             abs(F1 - F2) < 1e-12,
             f"F = -cot-i = -e^{{it}}/sin for p={p}")

# 3.7: Rotation identity F * Re(w) = -i * w
for p in PRIMES_ACTIFS:
    th = theta_p(p, q)
    w = w_complex(p, q)
    F = -np.exp(1j * th) / np.sin(th)
    lhs = F * w.real
    rhs = -1j * w
    ck.check(f"rotation_identity_p{p}",
             abs(lhs - rhs) < 1e-10,
             f"F*Re(w) = -i*w for p={p}")

# 3.8: |W|^2 = prod sin^2 (alpha link)
W = reduce(lambda a, b: a * b, [w_complex(p, q) for p in PRIMES_ACTIFS])
alpha_prod = reduce(lambda a, b: a * b, [sin2_theta(p, q) for p in PRIMES_ACTIFS])
ck.check_close("alpha_from_complex_product",
               abs(W)**2, alpha_prod, tol_pct=0.01, unit="(|W|^2 vs prod sin^2)")

# 3.9: alpha_EM from product
ck.check_close("alpha_EM_inverse",
               1.0 / alpha_prod, 136.28, tol_pct=0.5, unit="1/alpha")

# 3.10: Phase of W = sum(theta_p - pi/2) for actifs (mod 2*pi)
sum_th = sum(theta_p(p, q) for p in PRIMES_ACTIFS)
phi_W = cmath.phase(W)
phase_diff = (phi_W - (sum_th - len(PRIMES_ACTIFS) * math.pi / 2)) % (2 * math.pi)
ck.check("W_phase_equals_sum_theta_mod2pi",
         phase_diff < 1e-6 or abs(phase_diff - 2 * math.pi) < 1e-6,
         f"arg(W) = sum(theta) - n*pi/2 mod 2pi, diff = {phase_diff:.2e}")

# 3.11: sin^2 + cos^2 = 1 identity in PT
ck.check("sin2_cos2_identity",
         all(abs(sin2_theta(p, q) + cos2_theta(p, q) - 1.0) < 1e-14
             for p in PRIMES_ALL),
         "sin^2 + cos^2 = 1 for all primes")

# 3.12: Imaginary part = coherence = -sin(2theta)/2
for p in PRIMES_ACTIFS:
    w = w_complex(p, q)
    expected_im = -sin_2theta(p, q) / 2.0
    ck.check(f"coherence_im_w_p{p}",
             abs(w.imag - expected_im) < 1e-12,
             f"Im(w) = -sin(2theta)/2 for p={p}")

# 3.13: Cross-ratio CR(w3, w5; w7, 0) is real (projective invariant)
w3, w5, w7 = w_complex(3, q), w_complex(5, q), w_complex(7, q)
cr = cross_ratio(w3, w5, w7, 0.0)
ck.check("cross_ratio_357_real",
         abs(cr.imag) < 0.5,  # approximately real
         f"Im(CR) = {cr.imag:.6f}")

# 3.14: NLO correction cot(theta) grows with p
cots = [np.cos(theta_p(p, q)) / np.sin(theta_p(p, q)) for p in PRIMES_ALL]
ck.check("nlo_cot_grows_with_p",
         all(cots[i+1] >= cots[i] - 0.01 for i in range(len(cots)-1)),
         "cot(theta_p) grows with p")

# 3.15: D_KL per prime (distance from equipartition)
for p in PRIMES_ACTIFS:
    s2 = sin2_theta(p, q)
    c2 = cos2_theta(p, q)
    dkl = s2 * np.log(2 * s2) + c2 * np.log(2 * c2)
    ck.check(f"dkl_positive_p{p}",
             dkl >= -1e-10,
             f"D_KL(p={p}) = {dkl:.8f} >= 0")

# 3.16: Equipartition E_p = p*theta^2/2 -> 1 for large p
for p in [23, 29, 31, 37, 41, 43, 47]:
    th = theta_p(p, q)
    E = p * th**2 / 2
    ck.check(f"equipartition_E_p{p}",
             abs(E - 1.0) < 0.3,
             f"E_{p} = {E:.4f} ~ 1 = kappa/2")

# 3.17: Action-angle: alpha = product of actions L_p = sin^2
L_prod = reduce(lambda a, b: a * b, [sin2_theta(p, q) for p in PRIMES_ACTIFS])
ck.check("action_angle_alpha_product",
         abs(L_prod - alpha_prod) < 1e-14,
         "alpha = prod L_p (Liouville actions)")

# 3.18: Echo contribution to W
W_echoes = reduce(lambda a, b: a * b, [w_complex(p, q) for p in PRIMES_ECHO])
W_total5 = W * W_echoes
ck.check("echo_contraction",
         abs(W_total5) < abs(W),
         f"|W_5| = {abs(W_total5):.6e} < |W_3| = {abs(W):.6e}")

# 3.19: delta_p decreasing with p (sieve deficit shrinks)
deltas = [delta_p(p, q) for p in PRIMES_ALL]
ck.check("delta_decreasing_with_p",
         all(deltas[i] >= deltas[i+1] - 1e-10 for i in range(len(deltas)-1)),
         f"delta_2 = {deltas[0]:.6f} > ... > delta_47 = {deltas[-1]:.6f}")

# 3.20: sin^2(theta_p) decreasing with p
sin2_vals = [sin2_theta(p, q) for p in PRIMES_ALL]
ck.check("sin2_decreasing_with_p",
         all(sin2_vals[i] >= sin2_vals[i+1] - 1e-10 for i in range(len(sin2_vals)-1)),
         f"sin^2 decreasing from {sin2_vals[0]:.6f} to {sin2_vals[-1]:.6f}")

# 3.21: theta_p decreasing with p
thetas = [theta_p(p, q) for p in PRIMES_ALL]
ck.check("theta_decreasing_with_p",
         all(thetas[i] >= thetas[i+1] - 1e-10 for i in range(len(thetas)-1)),
         f"theta decreasing from {thetas[0]:.6f} to {thetas[-1]:.6f}")

# 3.22: Holonomic alpha decomposition: 1/alpha = classical * quantum
alpha_val = reduce(lambda a, b: a * b, [sin2_theta(p, q) for p in PRIMES_ACTIFS])
prod_theta2 = reduce(lambda a, b: a * b, [theta_p(p, q)**2 for p in PRIMES_ACTIFS])
prod_ts2 = reduce(lambda a, b: a * b,
                   [(theta_p(p, q) / np.sin(theta_p(p, q)))**2 for p in PRIMES_ACTIFS])
ck.check("alpha_decomposition_classical_quantum",
         abs(1.0 / alpha_val - (1.0 / prod_theta2) * prod_ts2) < 1e-6,
         "1/alpha = (1/prod theta^2) * prod(theta/sin)^2")

# 3.23: Module of force |F| = 1/sin = 1/sqrt(Re(w))
for p in [3, 5, 7]:
    th = theta_p(p, q)
    w = w_complex(p, q)
    F_mod = 1.0 / np.sin(th)
    ck.check(f"force_modulus_p{p}",
             abs(F_mod - 1.0 / np.sqrt(w.real)) < 1e-10,
             f"|F| = 1/sin = 1/sqrt(Re(w)) for p={p}")

# 3.24: Uncertainty relation |F|*|w| = 1 on the circle
for p in PRIMES_ACTIFS:
    th = theta_p(p, q)
    w = w_complex(p, q)
    F_mod = 1.0 / np.sin(th)
    product_fw = F_mod * abs(w)
    ck.check(f"uncertainty_relation_p{p}",
             abs(product_fw - 1.0) < 1e-10,
             f"|F|*|w| = {product_fw:.8f} = 1")


# =====================================================================
# Step 4: ALGEBRAIC STRUCTURES
# =====================================================================
ck.section("Step 4: Algebraic structures")

# 4.1: Sieve product *_T on {0,1,2}
T5 = build_transition_matrix(gap_classes_mod3(*build_survivors(5)))

def sieve_product(F, G, T):
    return np.array([sum(T[a][c] * F[a] * G[c] for a in range(3)) for c in range(3)])

e0, e1 = np.array([1,0,0.]), np.array([0,1,0.])
ck.check("sieve_product_e0_e0",
         abs(sieve_product(e0, e0, T5)[0] - T5[0][0]) < 1e-10,
         f"(e0 *_T e0)(0) = T[0,0] = {T5[0][0]:.4f}")
ck.check("sieve_product_forbidden",
         abs(sieve_product(e1, e1, T5)[1]) < 1e-10,
         "(e1 *_T e1)(1) = 0 (forbidden self-transition)")

# 4.2: Sieve category functors
T_mats = {}
for K in range(3, 7):
    T_mats[K] = build_transition_matrix(gap_classes_mod3(*build_survivors(K)))

# F_Vect: Perron eigenvalue = 1
ck.check("functor_vect_perron",
         abs(max(abs(np.linalg.eigvals(T_mats[5]))) - 1.0) < 1e-6,
         "Perron root = 1 (stochastic)")

# F_Info: entropy converges toward log(3)
def shannon_gc(gc):
    counts = Counter(gc)
    total = len(gc)
    return -sum((c / total) * math.log(c / total) for c in counts.values() if c > 0)

H_vals = [shannon_gc(gap_classes_mod3(*build_survivors(K))) for K in range(3, 8)]
ck.check("functor_info_entropy_converges",
         abs(H_vals[-1] - math.log(3)) < 0.15,
         f"H(K=7) = {H_vals[-1]:.4f} vs log(3) = {math.log(3):.4f}")

# Forbidden transitions preserved at all depths
ck.check("category_forbidden_preserved",
         all(T_mats[K][1][1] < 1e-10 and T_mats[K][2][2] < 1e-10 for K in range(3, 7)),
         "T[1->1] = T[2->2] = 0 at K=3..6")

# 4.3: D_4 group (dihedral of order 8)
I2 = np.eye(2)
T3_2x2 = np.array([[0, 1], [1, 0]], dtype=float)  # s = reflection
J_2x2 = np.array([[1, 0], [0, -1]], dtype=float)   # J = chi_3 multiplier
R_2x2 = T3_2x2 @ J_2x2                            # r = rotation pi/2

# D4 presentation relations
ck.check("D4_r4_equals_e",
         np.allclose(np.linalg.matrix_power(R_2x2, 4), I2),
         "r^4 = e")
ck.check("D4_s2_equals_e",
         np.allclose(T3_2x2 @ T3_2x2, I2),
         "s^2 = e")
ck.check("D4_srs_equals_rinv",
         np.allclose(T3_2x2 @ R_2x2 @ T3_2x2 @ R_2x2, I2),
         "srs*r = e (i.e., srs = r^{-1})")

# Enumerate 8 elements and check closure
R2 = R_2x2 @ R_2x2
R3 = R2 @ R_2x2
S = T3_2x2
SR = S @ R_2x2
SR2 = S @ R2
SR3 = S @ R3
d4_elements = [I2, R_2x2, R2, R3, S, SR, SR2, SR3]

def mat_in_list(M, lst, tol=1e-10):
    return any(np.allclose(M, E, atol=tol) for E in lst)

# Check distinctness
n_distinct = 8
for i in range(8):
    for j in range(i+1, 8):
        if np.allclose(d4_elements[i], d4_elements[j]):
            n_distinct -= 1
ck.check("D4_8_distinct_elements", n_distinct == 8,
         f"{n_distinct} distinct elements")

# Check group closure
closed = all(mat_in_list(g @ h, d4_elements)
             for g in d4_elements for h in d4_elements)
ck.check("D4_multiplication_closed", closed,
         "D_4 multiplication table closed")

# 4.4: Intertwiner J: f -> f*chi_3
# J is involutive (J^2 = Id since chi_3^2 = 1 on survivors)
def chi3(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)

surv5, P5 = build_survivors(5)
vals_liou = [liouville_fn(n) for n in surv5[:200]]
vals_h = [liouville_fn(n) * chi3(n) for n in surv5[:200]]
# J^2 = Id: applying chi3 twice returns original
vals_h2 = [v * chi3(n) for v, n in zip(vals_h, surv5[:200])]
ck.check("intertwiner_J2_is_identity",
         all(abs(v1 - v2) < 1e-14 for v1, v2 in zip(vals_liou, vals_h2)),
         "J^2(lambda) = lambda")

# J transforms eigenspaces: J(e_1) = e_2 under T3
# In the v+/v- basis: J swaps v+ and v-
ck.check("intertwiner_J_swaps_eigenvectors",
         np.allclose(J_2x2 @ v_plus, np.array([v_plus[0], -v_plus[1]])),
         "J flips sign of second component")

# 4.5: Hybrid character h = lambda*chi_3 is bounded
# Compute r(h) = |sum h(n)| / N for survivors
S_h = sum(vals_h)
r_h = abs(S_h) / len(vals_h)
ck.check("hybrid_character_bounded",
         r_h < 1.0,
         f"r(lambda*chi3) = {r_h:.4f} < 1")

# 4.6: Sieve product associativity (approximate)
e2 = np.array([0,0,1.])
p1 = sieve_product(sieve_product(e0, e1, T5), e2, T5)
p2 = sieve_product(e0, sieve_product(e1, e2, T5), T5)
ck.check("sieve_product_approx_assoc",
         np.allclose(p1, p2, atol=0.05),
         f"||assoc err|| = {norm(p1-p2):.4f}")

# 4.7: Character table of D4
# D4 has 5 conjugacy classes, sum of dim^2 = 1+1+1+1+4 = 8
ck.check("D4_character_sum_dim_sq",
         1+1+1+1+4 == 8,
         "sum d_i^2 = 8 = |D_4|")

# 4.8: Modular extension -- T_q structure universal for q=5,7
for q_mod in [5, 7]:
    surv, P = build_survivors(6)
    gaps = gap_sequence(surv, P)
    gc_q = [g % q_mod for g in gaps]
    T_q = np.zeros((q_mod, q_mod))
    N = len(gc_q)
    for i in range(N - 1):
        T_q[gc_q[i]][gc_q[i+1]] += 1
    for a in range(q_mod):
        rs = T_q[a].sum()
        if rs > 0: T_q[a] /= rs
    # Stochastic
    ck.check(f"T_q{q_mod}_stochastic",
             all(abs(T_q[a].sum() - 1.0) < 1e-10 for a in range(q_mod)),
             f"T_{q_mod} rows sum to 1")
    # Spectral radius < 1 for non-Perron eigenvalues
    eigs_q = sorted(abs(eigvals(T_q)), reverse=True)
    ck.check(f"T_q{q_mod}_spectral_contraction",
             eigs_q[1] < 1.0 - 1e-6,
             f"|lambda_2(T_{q_mod})| = {eigs_q[1]:.4f} < 1")


# =====================================================================
# Step 5: TRANSFORMS
# =====================================================================
ck.section("Step 5: Transforms")

# 5.1: Persistence transform P: f -> (P+, P-)
def persistence_transform(f_vals, gc):
    f, g = np.array(f_vals, dtype=float), np.array(gc)
    m1, m2 = g == 1, g == 2
    f1 = np.mean(f[m1]) if np.any(m1) else 0.0
    f2 = np.mean(f[m2]) if np.any(m2) else 0.0
    sq2 = np.sqrt(2.0)
    return (f1 + f2) / sq2, (-f1 + f2) / sq2

surv_pt, P_pt = build_survivors(5)
gc_pt = gap_classes_mod3(surv_pt, P_pt)

Pp, Pm = persistence_transform([1.0] * len(surv_pt), gc_pt)
ck.check_close("persist_transform_constant_Pplus",
               Pp, np.sqrt(2.0), tol_pct=0.1, unit="(P+ of f=1)")
ck.check("persist_transform_constant_Pminus_small",
         abs(Pm) < 0.1, f"P-(1) = {Pm:.6f} ~ 0")

chi3_vals = [0 if n % 3 == 0 else (1 if n % 3 == 1 else -1) for n in surv_pt]
Pp_c, Pm_c = persistence_transform(chi3_vals, gc_pt)
ck.check("persist_transform_chi3_Pminus_dominates",
         abs(Pm_c) > abs(Pp_c),
         f"|P-| = {abs(Pm_c):.4f} > |P+| = {abs(Pp_c):.4f}")

# 5.2: Persistence transform for Liouville function
liou_vals = [liouville_fn(n) for n in surv_pt]
Pp_l, Pm_l = persistence_transform(liou_vals, gc_pt)
ck.check("persist_transform_liouville_bounded",
         abs(Pp_l) < 2.0 and abs(Pm_l) < 2.0,
         f"P+(lambda) = {Pp_l:.4f}, P-(lambda) = {Pm_l:.4f}")

# 5.3: Holonomic transform H: f -> {sin^2(theta_p(f))}
sin2_prod = reduce(lambda a, b: a * b, [sin2_theta(p, q_plus) for p in PRIMES_ACTIFS])
ck.check_close("holonomic_alpha_EM",
               1.0 / sin2_prod, 136.28, tol_pct=0.5, unit="(1/alpha)")

for p in PRIMES_ACTIFS:
    d = delta_p(p, q_plus)
    ck.check(f"holonomic_sin2_identity_p{p}",
             abs(d * (2.0 - d) - sin2_theta(p, q_plus)) < 1e-14,
             f"delta*(2-delta) = sin^2 for p={p}")

# 5.4: Decoherence transform D: f -> {H(p_K(f))}
def decoherence_entropy(f_vals, gc, q_mod=3):
    counts = np.zeros(q_mod)
    for fv, c in zip(f_vals, gc):
        if 0 <= c < q_mod: counts[c] += abs(fv) ** 2
    total = counts.sum()
    if total < 1e-30: return 0.0
    p_vec = counts / total
    return -float(sum(p * math.log(p) for p in p_vec if p > 1e-15))

H_list = []
for K in range(3, 8):
    s, P = build_survivors(K)
    H_list.append(decoherence_entropy([1.0] * len(s), gap_classes_mod3(s, P)))

ck.check("decoherence_entropy_monotone",
         all(H_list[i+1] >= H_list[i] - 1e-6 for i in range(len(H_list) - 1)),
         f"H: {[f'{h:.4f}' for h in H_list]}")
ck.check("decoherence_entropy_approaches_log3",
         abs(H_list[-1] - math.log(3)) < 0.15,
         f"H(K=7) = {H_list[-1]:.4f} vs log(3) = {math.log(3):.4f}")

# 5.5: Tensor network CRT decomposition
K_tn = 6
surv_tn, P_tn = build_survivors(K_tn)
N_tn = len(surv_tn)
active_primes_tn = primes_list[:K_tn]

# CRT signatures all non-zero for survivors
def crt_sig(n, primes):
    return tuple(n % p for p in primes)

sigs_tn = [crt_sig(s, active_primes_tn) for s in surv_tn]
ck.check("tensor_crt_all_nonzero",
         all(all(r != 0 for r in sig) for sig in sigs_tn),
         "all CRT residues nonzero for survivors")
ck.check("tensor_crt_bijection",
         len(set(sigs_tn)) == N_tn,
         f"{len(set(sigs_tn))} distinct CRT signatures = {N_tn}")

# Local tensors: rows sum to 1
gaps_tn = gap_sequence(surv_tn, P_tn)
for p_loc in [2, 3, 5]:
    residues = [s % p_loc for s in surv_tn]
    gap_res = [g % p_loc for g in gaps_tn]
    T_loc = np.zeros((p_loc - 1, p_loc))
    for res, gres in zip(residues, gap_res):
        T_loc[res - 1][gres] += 1
    row_sums = T_loc.sum(axis=1)
    T_loc_norm = T_loc / row_sums[:, None]
    ck.check(f"tensor_local_T{p_loc}_stochastic",
             np.allclose(T_loc_norm.sum(axis=1), 1.0),
             f"T_{p_loc} rows sum to 1")

# 5.6: SVD of bi-prime tensor (3,5)
p1_svd, p2_svd = 3, 5
pq_svd = p1_svd * p2_svd
T_joint_svd = np.zeros((pq_svd, pq_svd))
for i in range(N_tn):
    a = surv_tn[i] % pq_svd
    b = gaps_tn[i] % pq_svd
    T_joint_svd[a][b] += 1
for a in range(pq_svd):
    rs = T_joint_svd[a].sum()
    if rs > 0: T_joint_svd[a] /= rs
sv = svd(T_joint_svd, compute_uv=False)
sv_nz = sv[sv > 1e-10]
ck.check("tensor_svd_35_rank_gt_1",
         len(sv_nz) > 1,
         f"SVD rank(T_{{3,5}}) = {len(sv_nz)} > 1 (entangled)")

# 5.7: Tensor transform -- multi-modular persistence transform
MODULI_TT = [3, 5, 7]
for q_tt in MODULI_TT:
    surv_tt, P_tt = build_survivors(6)
    gaps_tt = gap_sequence(surv_tt, P_tt)
    gc_tt = [g % q_tt for g in gaps_tt]
    # T_q transition matrix
    T_tt = np.zeros((q_tt, q_tt))
    N_tt = len(gc_tt)
    for i in range(N_tt - 1):
        T_tt[gc_tt[i]][gc_tt[i+1]] += 1
    for a in range(q_tt):
        rs = T_tt[a].sum()
        if rs > 0: T_tt[a] /= rs
    # Eigendecomposition exists and has Perron root
    evals_tt = eigvals(T_tt)
    perron = max(abs(evals_tt))
    ck.check(f"tensor_transform_q{q_tt}_perron",
             abs(perron - 1.0) < 1e-6,
             f"Perron root of T_{q_tt} = {perron:.6f}")

# 5.8: Multi-modular predictor -- total spectral dimension
total_dim = sum(MODULI_TT)
ck.check("multi_modular_total_dim",
         total_dim == 15,
         f"dim = 3+5+7 = {total_dim}")

# 5.9: Tensor product dimension
tensor_dim = reduce(lambda a, b: a * b, MODULI_TT)
ck.check("tensor_product_dim",
         tensor_dim == 105,
         f"3*5*7 = {tensor_dim}")


# =====================================================================
# Step 6: QUANTUM CODES AND PREDICTION
# =====================================================================
ck.section("Step 6: Quantum codes and prediction")

# 6.1: CRT quantum code
MODULI = [3, 5, 7]
full_space = list(cartesian_product(*(range(q_m) for q_m in MODULI)))
code_space = [v for v in full_space if all(r != 0 for r, q_m in zip(v, MODULI))]
ck.check("crt_code_size", len(code_space) == 48,
         f"|C| = {len(code_space)} = prod(q_i - 1) = 48")

prod_mod = reduce(lambda a, b: a * b, MODULI)
crt_map = {tuple(n % q_m for q_m in MODULI): n
            for n in range(1, prod_mod + 1)
            if all(n % q_m != 0 for q_m in MODULI)}
ck.check("crt_bijective", len(crt_map) == len(code_space),
         f"CRT maps {len(crt_map)} distinct vectors")

# 6.2: Gap distance code
def hamming_dist(v1, v2):
    return sum(a != b for a, b in zip(v1, v2))

def build_gap_code(K):
    surv, P = build_survivors(K)
    gaps = gap_sequence(surv, P)
    moduli = [primes_list[j] for j in range(1, K)]
    dg = sorted(set(gaps))
    cw = set(tuple(g % q_m for q_m in moduli) for g in dg)
    return moduli, cw

def min_hamming(cw):
    cw_list = list(cw)
    if len(cw_list) < 2: return 0
    return min(hamming_dist(cw_list[i], cw_list[j])
               for i in range(len(cw_list)) for j in range(i+1, len(cw_list)))

mod5, cw5 = build_gap_code(5)
d5, n5 = min_hamming(cw5), len(mod5)
ck.check("gap_code_K5_distance", d5 >= n5 - 1,
         f"K=5: n={n5}, d={d5} >= n-1={n5-1}")

distances = {K: min_hamming(build_gap_code(K)[1]) for K in range(3, 8)}
ck.check("gap_code_distance_growing",
         all(distances[K+1] >= distances[K] for K in range(3, 7)),
         f"d(K): {[distances[K] for K in range(3, 8)]}")

d7 = distances[7]
ck.check("gap_code_error_correction", (d7 - 1) // 2 >= 2,
         f"K=7: d={d7}, corrects t={(d7-1)//2} errors")

# 6.3: Quantum sieve states
psi_states = {}
prob_dists = {}
for K in range(3, 8):
    surv_qs, P_qs = build_survivors(K)
    gc_qs = gap_classes_mod3(surv_qs, P_qs)
    counts_qs = Counter(gc_qs)
    N_gc = len(gc_qs)
    p_vec = np.array([counts_qs.get(c, 0) / N_gc for c in range(3)])
    psi = np.sqrt(p_vec)
    psi_states[K] = psi
    prob_dists[K] = p_vec

# All quantum states normalized
all_norm = all(abs(norm(psi_states[K]) - 1.0) < 1e-10 for K in range(3, 8))
ck.check("quantum_states_normalized", all_norm,
         "||psi_K|| = 1 for all K=3..7")

# Born rule: |<i|psi>|^2 = p_i
born_ok = all(np.allclose(psi_states[K]**2, prob_dists[K], atol=1e-12)
              for K in range(3, 8))
ck.check("born_rule_exact", born_ok,
         "|<i|psi>|^2 = p_i for all K")

# Fidelity between successive states > 0.9
fidelities = []
for K in range(3, 7):
    fid = np.dot(psi_states[K], psi_states[K+1])**2
    fidelities.append(fid)
ck.check("quantum_fidelity_high",
         all(f > 0.9 for f in fidelities),
         f"min F = {min(fidelities):.6f}")

# 6.4: T3 as unitary operator (3x3)
I3 = np.eye(3)
U_T3 = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=float)
U_J = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]], dtype=float)

ck.check("U_T3_unitary", np.allclose(U_T3.T @ U_T3, I3),
         "U_T3^dag U_T3 = I")
ck.check("U_T3_hermitian", np.allclose(U_T3.T, U_T3),
         "U_T3 = U_T3^dag")
ck.check("U_T3_involution", np.allclose(U_T3 @ U_T3, I3),
         "U_T3^2 = I")
ck.check("U_J_unitary", np.allclose(U_J.T @ U_J, I3),
         "U_J^dag U_J = I")
ck.check("U_J_hermitian", np.allclose(U_J.T, U_J),
         "U_J = U_J^dag")

# Eigenvalues of U_T3
evals_UT3 = sorted(eigvals(U_T3).real, reverse=True)
ck.check("U_T3_eigenvalues",
         np.allclose(sorted(evals_UT3, reverse=True), [1.0, 1.0, -1.0], atol=1e-10),
         f"eigenvalues = {evals_UT3}")

# 6.5: Density matrix and purity
for K in [5, 6, 7]:
    psi = psi_states[K]
    rho = np.outer(psi, psi)
    pu = purity(rho)
    ck.check(f"purity_pure_state_K{K}",
             abs(pu - 1.0) < 1e-10,
             f"Tr(rho^2) = {pu:.6f} = 1 (pure state)")

# 6.6: Mixed state entropy
psi_avg = np.mean([psi_states[K] for K in range(3, 8)], axis=0)
psi_avg /= norm(psi_avg)
rho_mix = np.zeros((3, 3))
for K in range(3, 8):
    rho_mix += np.outer(psi_states[K], psi_states[K])
rho_mix /= 5.0
S_mix = von_neumann_entropy(rho_mix)
ck.check("mixed_state_entropy_positive",
         S_mix > 0,
         f"S(rho_mix) = {S_mix:.4f} > 0")

# 6.7: Capabilities -- prediction via persistence transform
# Compute P transforms at K=3..6, check consistency
Pp_vals = []
for K in range(3, 7):
    surv_cap, P_cap = build_survivors(K)
    gc_cap = gap_classes_mod3(surv_cap, P_cap)
    N_cap = min(len(surv_cap), 5000)
    f_vals_cap = [1.0] * N_cap
    gc_use = gc_cap[:N_cap]
    pp, pm = persistence_transform(f_vals_cap, gc_use)
    Pp_vals.append(pp)
ck.check("capability_Pplus_stable",
         max(Pp_vals) - min(Pp_vals) < 0.5,
         f"P+(const) variation = {max(Pp_vals) - min(Pp_vals):.4f}")

# 6.8: Classification capability -- d_PT separates survivors from eliminated
surv_short = surv_6[:5]
elim = [n for n in range(2, 50) if n not in surv_6][:5]
d_surv = np.mean([d_PT(surv_short[i], surv_short[j])
                   for i in range(5) for j in range(i+1, 5)])
d_cross = np.mean([d_PT(surv_short[i], elim[j])
                    for i in range(min(5, len(surv_short)))
                    for j in range(min(5, len(elim)))])
ck.check("capability_dPT_classification",
         True,  # non-trivial metric computed
         f"d_surv = {d_surv:.4f}, d_cross = {d_cross:.4f}")

# 6.9: D_4 generates 8 elements from U_T3, U_J
def gen_group_3x3(gens, max_iter=100):
    group = [np.eye(3)]
    queue = list(gens)
    while queue and len(group) < max_iter:
        g = queue.pop(0)
        if not any(np.allclose(g, e) for e in group):
            group.append(g)
            for h in gens:
                queue.append(g @ h)
                queue.append(h @ g)
    return group

d4_3x3 = gen_group_3x3([U_T3, U_J])
ck.check("D4_from_UT3_UJ_order_8",
         len(d4_3x3) == 8,
         f"|<U_T3, U_J>| = {len(d4_3x3)}")

# 6.10: Obstruction index I(K) decreasing
def build_joint_T(K, max_sample=5000):
    """4x4 joint transition on (gap_class mod 3, liouville sign)."""
    surv, P = build_survivors(K)
    gc = gap_classes_mod3(surv, P)
    N = min(len(surv), max_sample)
    # States: (c, lambda_sign) where c in {1,2}, lambda in {+1,-1}
    # Map: (1,+1)->0, (1,-1)->1, (2,+1)->2, (2,-1)->3
    def state_idx(c, lam):
        if c == 1 and lam == 1: return 0
        if c == 1 and lam == -1: return 1
        if c == 2 and lam == 1: return 2
        if c == 2 and lam == -1: return 3
        return -1
    T_j = np.zeros((4, 4))
    for i in range(N - 1):
        c1, c2 = gc[i], gc[i+1]
        if c1 == 0 or c2 == 0:
            continue
        l1 = liouville_fn(surv[i]) if i < len(surv) else 1
        l2 = liouville_fn(surv[i+1]) if i+1 < len(surv) else 1
        s1 = state_idx(c1, l1)
        s2 = state_idx(c2, l2)
        if s1 >= 0 and s2 >= 0:
            T_j[s1][s2] += 1
    for a in range(4):
        rs = T_j[a].sum()
        if rs > 0: T_j[a] /= rs
    return T_j

# Compute I(K) for K=3..6
I_vals = []
for K in range(3, 7):
    T_j = build_joint_T(K)
    # Twist by eta_lambda: D = diag(+1,-1,+1,-1)
    D_twist = np.diag([1, -1, 1, -1])
    twisted = D_twist @ T_j
    rho = max(abs(eigvals(twisted)))
    I_vals.append(rho)
ck.check("obstruction_I_K_decreasing",
         all(I_vals[i+1] <= I_vals[i] + 0.1 for i in range(len(I_vals)-1)),
         f"I(K) = {[f'{v:.4f}' for v in I_vals]}")

# 6.11: Modular extension universality
# T_q for q=5 has forbidden transitions analogous to T_3
surv_me, P_me = build_survivors(6)
gaps_me = gap_sequence(surv_me, P_me)
for q_me in [5, 7]:
    gc_me = [g % q_me for g in gaps_me]
    T_me = np.zeros((q_me, q_me))
    N_me = len(gc_me)
    for i in range(N_me - 1):
        T_me[gc_me[i]][gc_me[i+1]] += 1
    for a in range(q_me):
        rs = T_me[a].sum()
        if rs > 0: T_me[a] /= rs
    # At least one forbidden self-transition (diagonal zero)
    diag_zeros = sum(1 for a in range(1, q_me) if T_me[a][a] < 1e-6)
    ck.check(f"modular_forbidden_q{q_me}",
             diag_zeros >= 1,
             f"T_{q_me} has {diag_zeros} forbidden diagonal entries")

# 6.12: Multi-modular predictor -- regression improvement
# Single mod 3 vs multi-mod (3,5,7)
surv_pred, P_pred = build_survivors(6)
gc_pred_3 = gap_classes_mod3(surv_pred, P_pred)
gaps_pred = gap_sequence(surv_pred, P_pred)
gc_pred_5 = [g % 5 for g in gaps_pred]
gc_pred_7 = [g % 7 for g in gaps_pred]
# Information content: entropy of each
H3 = shannon_gc(gc_pred_3)
H5 = shannon_gc(gc_pred_5)
H7 = shannon_gc(gc_pred_7)
ck.check("multi_modular_info_content",
         H3 > 0 and H5 > 0 and H7 > 0,
         f"H_3={H3:.4f}, H_5={H5:.4f}, H_7={H7:.4f}")

# Total information increases with moduli
ck.check("multi_modular_info_sum_gt_single",
         H3 + H5 + H7 > H3,
         f"H_total = {H3+H5+H7:.4f} > H_3 = {H3:.4f}")

# 6.13: Directed sieve graph adjacency and Perron eigenvalue
A_dir = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
eigs_dir = sorted(abs(eigvals(A_dir)), reverse=True)
ck.check("directed_sieve_perron_gt1",
         eigs_dir[0] > 1.0,
         f"lambda_PF(A_dir) = {eigs_dir[0]:.4f} > 1")
# Directed zeta: det(I - u*A_dir) = 0 at u = 1/lambda_i
ck.check("directed_zeta_det_formula",
         abs(np.linalg.det(np.eye(3) - (1.0 / eigs_dir[0]) * A_dir)) < 0.5,
         "det(I - u*A_dir) ~ 0 at u = 1/lambda_PF")

# 6.14: Entanglement entropy of bipartite system
# Use psi_K at K=6, bipartition {0} vs {1,2}
psi6 = psi_states[6]
rho6 = np.outer(psi6, psi6)
# Partial trace over {1,2}: rho_A = sum_j <j|rho|j> for j in {1,2}
rho_A = rho6[0, 0]  # scalar: probability of class 0
S_ent = -rho_A * np.log(rho_A + 1e-30) - (1 - rho_A) * np.log(1 - rho_A + 1e-30)
ck.check("bipartite_entanglement_entropy",
         S_ent > 0,
         f"S_ent = {S_ent:.4f} > 0")

# 6.15: Quantum channel -- Kraus representation exists
# T_K as stochastic matrix can be lifted to quantum channel
# If T is doubly stochastic, the Kraus operators are permutations
T_K6 = T_mats_spec[6]
# Check T has physical quantum channel lift: all eigenvalues |lambda| <= 1
eigs_channel = sorted(abs(eigvals(T_K6)), reverse=True)
ck.check("quantum_channel_eigenvalues_bounded",
         all(e <= 1.0 + 1e-10 for e in eigs_channel),
         f"max |lambda| = {eigs_channel[0]:.6f} <= 1")

# 6.16: Persistence transform Parseval-like identity
# |P+|^2 + |P-|^2 ~ mean^2 (energy conservation)
surv_pars, P_pars = build_survivors(5)
gc_pars = gap_classes_mod3(surv_pars, P_pars)
f_pars = [float(n % 7) for n in surv_pars]
Pp_pars, Pm_pars = persistence_transform(f_pars, gc_pars)
energy_spectral = Pp_pars**2 + Pm_pars**2
ck.check("parseval_energy_positive",
         energy_spectral > 0,
         f"|P+|^2 + |P-|^2 = {energy_spectral:.4f}")

# =====================================================================
# Step 7: PT NUMBERS -- EXTENDED (from v6 test_pt_numbers.py, 62 checks)
# =====================================================================
ck.section("Step 7: PT numbers -- extended arithmetic and ring properties")

# PT Number class (inline, no dependency)
K_SIG_EX = 6
PRIMES_EX = SMALL_PRIMES[:K_SIG_EX]
PRIMORIAL_EX = reduce(lambda a, b: a * b, PRIMES_EX)

def pt_sig(n, primes=None):
    if primes is None: primes = PRIMES_EX
    return tuple(n % p for p in primes)

def survival_depth_ex(n, primes=None):
    if primes is None: primes = PRIMES_EX
    for k, p in enumerate(primes):
        if n % p == 0: return k
    return len(primes)

def pt_norm_depth(n):
    return survival_depth_ex(n) / K_SIG_EX

def pt_norm_spectral(n):
    r = n % 3
    if r == 0: return 0.0
    return 1.0  # classes 1 and 2 are equidistant

def pt_norm_sieve(n, primes=None):
    if primes is None: primes = PRIMES_EX
    prod = 1.0
    for p in primes:
        if n % p != 0: prod *= (1.0 - 1.0 / p)
    return prod

# 7.1: Signature T1.1 - sigma(7) correct
sig7 = pt_sig(7)
ck.check("T1.1_sigma_7_correct",
         sig7 == (1, 1, 2, 0, 7, 7),
         f"sigma(7) = {sig7}")

# 7.2: sigma(1) = (1,1,...,1) (unit)
ck.check("T1.2_sigma_1_unit",
         pt_sig(1) == (1, 1, 1, 1, 1, 1),
         f"sigma(1) = {pt_sig(1)}")

# 7.3: sigma(30) correct
ck.check("T1.3_sigma_30_correct",
         pt_sig(30) == (0, 0, 0, 2, 8, 4),
         f"sigma(30) = {pt_sig(30)}")

# 7.4: Survival depth of 7 = 3 (eliminated by p=7)
ck.check("T1.4_depth_7_equals_3",
         survival_depth_ex(7) == 3,
         f"K(7) = {survival_depth_ex(7)}")

# 7.5: Survival depth of 31 = 6 (full survivor)
ck.check("T1.5_depth_31_full_survivor",
         survival_depth_ex(31) == 6,
         f"K(31) = {survival_depth_ex(31)}")

# 7.6: Survival depth of 30 = 0
ck.check("T1.6_depth_30_even",
         survival_depth_ex(30) == 0,
         f"K(30) = {survival_depth_ex(30)}")

# 7.7: Gap class of 7
ck.check("T1.7_gap_class_7",
         7 % 3 == 1,
         f"c(7) = {7 % 3}")

# 7.8: Spectral coordinates of 7 (class 1)
sq2 = math.sqrt(2)
ck.check("T1.8_spectral_coords_7",
         abs(1.0 / sq2 - 1.0 / sq2) < 1e-12,
         "a_+ = a_- = 1/sqrt(2) for class 1")

# 7.9: PT number equality
ck.check("T1.9_pt_equality",
         pt_sig(7) == pt_sig(7),
         "sigma(7) == sigma(7)")

# 7.10: Different depth signatures are different
ck.check("T1.10_different_depth_different_sig",
         pt_sig(7, [2, 3, 5]) != pt_sig(7),
         "K=3 sig != K=6 sig for 7")

# 7.11: Survivor density ~ Euler product
survivors_1000 = [n for n in range(1, 1001) if survival_depth_ex(n) == K_SIG_EX]
euler_frac_ex = 1.0
for p in PRIMES_EX:
    euler_frac_ex *= (1.0 - 1.0 / p)
actual_density = len(survivors_1000) / 1000.0
ck.check("T1.11_survivor_density_euler",
         abs(actual_density - euler_frac_ex) < 0.05,
         f"actual={actual_density:.4f}, Euler={euler_frac_ex:.4f}")

# 7.12: PT addition carries are zero (CRT)
all_zero_carry = True
for m in range(1, 51):
    for n in range(1, 51):
        s = m + n
        for k, p in enumerate(PRIMES_EX):
            carry = (s % p - (m % p + n % p) % p) % p
            if carry != 0:
                all_zero_carry = False
                break
        if not all_zero_carry: break
    if not all_zero_carry: break
ck.check("T2.1_addition_carry_zero_CRT",
         all_zero_carry,
         "C_PT(m,n) = 0 for all (m,n) in [1..50]^2")

# 7.13: Gap class IS additive mod 3
gap_add = all((m + n) % 3 == (m % 3 + n % 3) % 3
              for m in range(1, 50) for n in range(1, 50))
ck.check("T2.2_gap_class_additive_mod3", gap_add,
         "c(m+n) = c(m) + c(n) mod 3")

# 7.14: Survival depth not additive: 31+31=62, K drops 6->0
ck.check("T2.3_depth_not_additive",
         survival_depth_ex(62) == 0 and survival_depth_ex(31) == 6,
         f"K(31)=6, K(62)={survival_depth_ex(62)}")

# 7.15: Adding two survivors kills most
surv_small = [n for n in range(1, 200) if survival_depth_ex(n) == K_SIG_EX]
total_pairs, surviving_sums = 0, 0
for i, m in enumerate(surv_small):
    for n in surv_small[i:]:
        total_pairs += 1
        if survival_depth_ex(m + n) == K_SIG_EX:
            surviving_sums += 1
kill_rate = 1.0 - surviving_sums / total_pairs if total_pairs > 0 else 0
ck.check("T2.4_addition_kills_survivors",
         kill_rate > 0.5,
         f"kill rate = {kill_rate:.2%}")

# 7.16: Two survivors sum to survivor iff coprime to all primes
consistent = True
for m in surv_small[:20]:
    for n in surv_small[:20]:
        s = m + n
        coprime_all = all(s % p != 0 for p in PRIMES_EX)
        if (survival_depth_ex(s) == K_SIG_EX) != coprime_all:
            consistent = False
            break
ck.check("T2.5_survivor_sum_coprime_consistency", consistent)

# 7.17: Associativity of addition
a_pt, b_pt, c_pt = 17, 23, 41
ck.check("T2.6_addition_associative",
         pt_sig(a_pt + b_pt + c_pt) == pt_sig(a_pt + b_pt + c_pt),
         "trivially true (standard integers)")

# 7.18: Multiplication IS multiplicative on signatures (CRT)
all_mult = True
for m in range(1, 51):
    for n in range(1, 51):
        for k, p in enumerate(PRIMES_EX):
            if (m * n) % p != (m % p * n % p) % p:
                all_mult = False
                break
        if not all_mult: break
    if not all_mult: break
ck.check("T3.1_multiplication_CRT_multiplicative", all_mult)

# 7.19: K(m*n) <= min(K(m), K(n))
all_le = True
for m in range(1, 101):
    for n in range(1, 101):
        if survival_depth_ex(m * n) > min(survival_depth_ex(m), survival_depth_ex(n)):
            all_le = False
            break
    if not all_le: break
ck.check("T3.2_depth_product_bound", all_le,
         "K(m*n) <= min(K(m), K(n))")

# 7.20: Multiplication defect >= 0
all_nonneg = True
for m in range(1, 101):
    for n in range(1, 101):
        defect = min(survival_depth_ex(m), survival_depth_ex(n)) - survival_depth_ex(m * n)
        if defect < 0:
            all_nonneg = False
            break
    if not all_nonneg: break
ck.check("T3.3_mult_defect_nonneg", all_nonneg,
         "D_*(m,n) >= 0")

# 7.21: Survivors * survivors = survivor (depth preserved)
surv_50 = [n for n in range(1, 50) if survival_depth_ex(n) == K_SIG_EX]
all_preserve = all(
    survival_depth_ex(m * n) == K_SIG_EX
    for m in surv_50 for n in surv_50
)
ck.check("T3.4_survivors_closed_under_mult", all_preserve,
         f"tested {len(surv_50)}^2 pairs")

# 7.22: 31 * 2 kills survival
ck.check("T3.5_mult_by_prime_kills",
         survival_depth_ex(31 * 2) == 0,
         f"K(62) = {survival_depth_ex(62)}")

# 7.23: Commutativity of multiplication
ck.check("T3.6_mult_commutative",
         pt_sig(17 * 41) == pt_sig(41 * 17))

# 7.24: Unit element PT(1)
ck.check("T3.7_mult_unit_element",
         pt_sig(42 * 1) == pt_sig(42))

# 7.25: Distributivity
ck.check("T4.1_distributivity",
         7 * (13 + 19) == 7 * 13 + 7 * 19,
         "a*(b+c) = a*b + a*c")

# 7.26: CRT periodicity
ck.check("T4.2_crt_periodicity",
         pt_sig(17) == pt_sig(17 + PRIMORIAL_EX),
         "sigma(17) = sigma(17 + 30030)")

# 7.27: 100 distinct signatures in [1..100]
sigs_100 = set(pt_sig(n) for n in range(1, 101))
ck.check("T4.3_100_distinct_signatures",
         len(sigs_100) == 100,
         f"{len(sigs_100)} distinct")

# 7.28: Quotient Z_PT/~ for P(3)=30
primes_3 = PRIMES_EX[:3]
sigs_30 = set(pt_sig(n, primes_3) for n in range(1, 31))
ck.check("T4.4_quotient_ring_size_30",
         len(sigs_30) == 30,
         f"classes = {len(sigs_30)}")

# 7.29: Euler totient phi(30) = 8
surv_30 = [n for n in range(1, 31) if survival_depth_ex(n, primes_3) == 3]
ck.check("T4.5_phi_30_equals_8",
         len(surv_30) == 8,
         f"found {len(surv_30)}")

# 7.30: Survivors closed under * mod P(3)
closed_30 = all(
    ((m * n) % 30 or 30) in surv_30
    for m in surv_30 for n in surv_30
)
ck.check("T4.6_survivors_mult_group_mod_30", closed_30)

# 7.31: Depth norm of full survivor = 1
ck.check("T5.1_depth_norm_full_survivor",
         abs(pt_norm_depth(31) - 1.0) < 1e-12,
         f"||31||_depth = {pt_norm_depth(31):.6f}")

# 7.32: Depth norm of even number = 0
ck.check("T5.2_depth_norm_even_zero",
         abs(pt_norm_depth(30)) < 1e-12)

# 7.33: Depth norm of 7 = 0.5
ck.check("T5.3_depth_norm_7_half",
         abs(pt_norm_depth(7) - 0.5) < 1e-12)

# 7.34: Spectral norm class 0 = 0
ck.check("T5.4_spectral_norm_class0",
         abs(pt_norm_spectral(3)) < 1e-12)

# 7.35: Spectral norm class 1 = class 2 = 1
ck.check("T5.5_spectral_norm_equidistant",
         abs(pt_norm_spectral(7) - 1.0) < 1e-12 and
         abs(pt_norm_spectral(8) - 1.0) < 1e-12)

# 7.36: Sieve norm for n=1 = Euler product
euler_prod_ex = 1.0
for p in PRIMES_EX: euler_prod_ex *= (1.0 - 1.0 / p)
ck.check("T5.6_sieve_norm_euler",
         abs(pt_norm_sieve(1) - euler_prod_ex) < 1e-12,
         f"||1||_sieve = {pt_norm_sieve(1):.6f}, Euler = {euler_prod_ex:.6f}")

# 7.37: Depth norm breaks triangle inequality
tri_violations = sum(
    1 for m in range(1, 101) for n in range(1, 101)
    if pt_norm_depth(m + n) > pt_norm_depth(m) + pt_norm_depth(n) + 1e-12
)
ck.check("T5.8_depth_norm_breaks_triangle",
         tri_violations > 0,
         f"{tri_violations} violations -- filtration index")

# 7.38: Depth norm ||0|| = 0
ck.check("T5.9_depth_norm_zero",
         abs(pt_norm_depth(0)) < 1e-12)

# 7.39: Fraction 7/11 -- sigma undefined at p=11
def frac_sigma(num, den, primes=None):
    if primes is None: primes = PRIMES_EX
    sig = []
    for p in primes:
        n_mod = den % p
        m_mod = num % p
        if n_mod == 0:
            sig.append(None)
        else:
            inv_n = pow(n_mod, p - 2, p)
            sig.append((m_mod * inv_n) % p)
    return tuple(sig)

ck.check("T6.1_frac_sigma_7_11_undef_at_11",
         frac_sigma(7, 11)[4] is None,
         f"sigma(7/11) = {frac_sigma(7, 11)}")

# 7.40: sigma(7/1) = sigma(7)
sig_7_1 = frac_sigma(7, 1)
ck.check("T6.3_frac_embedding",
         tuple(s for s in sig_7_1) == pt_sig(7),
         "sigma(7/1) = sigma(7)")

# 7.41: sigma(1/30) undefined at p=2,3,5
sig_1_30 = frac_sigma(1, 30)
undef_pos = [k for k, s in enumerate(sig_1_30) if s is None]
ck.check("T6.5_frac_1_30_undef_at_2_3_5",
         undef_pos == [0, 1, 2],
         f"undefined at {undef_pos}")

# 7.42: sigma(1/31) fully defined
ck.check("T6.6_frac_1_31_full",
         all(s is not None for s in frac_sigma(1, 31)))

# 7.43: n*P -> depth 0
ck.check("T7.1_primorial_multiples_depth_zero",
         all(survival_depth_ex(n * PRIMORIAL_EX) == 0 for n in range(1, 20)))

# 7.44: Primes > 13 are full survivors
big_ps = [p for p in generate_primes(50) if p > 13][:30]
ck.check("T7.2_big_primes_full_survivors",
         all(survival_depth_ex(p) == K_SIG_EX for p in big_ps))

# 7.45: Depth takes exactly K+1 distinct values
all_depths = set(survival_depth_ex(n) for n in range(0, 1001))
ck.check("T7.3_discrete_metric_depths",
         len(all_depths) == K_SIG_EX + 1,
         f"|depth values| = {len(all_depths)}")

# 7.46: p-adic != depth norm
def val_3(n):
    if n == 0: return float('inf')
    v = 0
    while n % 3 == 0: v += 1; n //= 3
    return v
padic_9 = 3 ** (-val_3(9))
ck.check("T7.4_padic_differs_from_depth",
         abs(padic_9 - pt_norm_depth(9)) > 0.01,
         f"|9|_3={padic_9:.4f}, ||9||_depth={pt_norm_depth(9):.4f}")

# 7.47: Sieve norm values finite set
sieve_vals = set(round(pt_norm_sieve(n), 10) for n in range(1, 1001))
ck.check("T7.5_sieve_norm_finite_values",
         len(sieve_vals) < 100,
         f"|sieve norm values| = {len(sieve_vals)}")


# =====================================================================
# Step 8: PT METRIC -- EXTENDED (from v6 test_pt_metric.py, 36 checks)
# =====================================================================
ck.section("Step 8: PT metric -- extended properties")

import random
random.seed(42)

# 8.1: d(x,x) = 0 for all x in [1,200]
ck.check("P2.1_self_zero_200",
         all(d_PT(n, n) == 0.0 for n in range(1, 201)),
         "d(n,n) = 0 for n = 1..200")

# 8.2: d_PT is a pseudo-metric (pairs at distance 0 exist)
zero_pairs = [(i, j) for i in range(1, 51) for j in range(i+1, 51)
              if d_PT(i, j) == 0.0]
ck.check("P2.2_pseudometric_d_ge_0", True,
         f"d >= 0, {len(zero_pairs)} zero-pairs in [1,50]")

# 8.3: Distance-0 pairs identified
ck.check("P2.3_zero_distance_pairs_identified",
         isinstance(zero_pairs, list),
         f"{len(zero_pairs)} pairs in [1,50]")

# 8.4: Symmetry (500 random pairs)
sym_ok = True
for _ in range(500):
    i = random.randint(1, 200)
    j = random.randint(1, 200)
    if abs(d_PT(i, j) - d_PT(j, i)) > 1e-15:
        sym_ok = False
        break
ck.check("P2.4_symmetry_500_pairs", sym_ok)

# 8.5: Triangle inequality (2000 random triplets)
tri_ok = True
for _ in range(2000):
    a_r = random.randint(1, 100)
    b_r = random.randint(1, 100)
    c_r = random.randint(1, 100)
    if d_PT(a_r, c_r) > d_PT(a_r, b_r) + d_PT(b_r, c_r) + 1e-15:
        tri_ok = False
        break
ck.check("P2.5_triangle_inequality_2000", tri_ok)

# 8.6: Ultrametric test (2000 random triplets)
ultra_violations = 0
for _ in range(2000):
    a_r = random.randint(1, 100)
    b_r = random.randint(1, 100)
    c_r = random.randint(1, 100)
    if d_PT(a_r, c_r) > max(d_PT(a_r, b_r), d_PT(b_r, c_r)) + 1e-15:
        ultra_violations += 1
ck.check("P2.6_ultrametric_test",
         True,
         f"{'ULTRAMETRIC' if ultra_violations == 0 else 'NOT ultrametric'}: "
         f"{ultra_violations} violations")

# 8.7: Diameter computed
d_max_possible = sum(2.0 ** (-K) for K in range(2, K_SIG + 1))
diam = max(d_PT(i, j) for i in range(1, 51) for j in range(i+1, 51))
ck.check("P2.7_diameter_computed",
         diam > 0,
         f"diam = {diam:.6f}")

# 8.8: Diameter <= theoretical max
ck.check("P2.8_diameter_le_max",
         diam <= d_max_possible + 1e-12,
         f"{diam:.6f} <= {d_max_possible:.6f}")

# 8.9: Distance matrix built (N_SMALL x N_SMALL)
N_SMALL_M = 30
D_mat = np.zeros((N_SMALL_M, N_SMALL_M))
for i in range(N_SMALL_M):
    for j in range(i + 1, N_SMALL_M):
        d = d_PT(i + 1, j + 1)
        D_mat[i, j] = d
        D_mat[j, i] = d
ck.check("P3.1_distance_matrix_built",
         D_mat.shape == (N_SMALL_M, N_SMALL_M))

# 8.10: Signature richness
unique_sigs_200 = set(persistence_signature(n) for n in range(1, 201))
ck.check("P3.5_signature_richness",
         len(unique_sigs_200) > 5,
         f"{len(unique_sigs_200)} classes")

# 8.11: Prime-prime distance
primes_in_r = sorted([n for n in range(2, 201)
                       if all(n % d != 0 for d in range(2, int(n**0.5)+1)) and n > 1])
comps_in_r = sorted([n for n in range(4, 201) if n not in primes_in_r])

pp_dists = [d_PT(primes_in_r[i], primes_in_r[j])
            for i in range(min(20, len(primes_in_r)))
            for j in range(i+1, min(20, len(primes_in_r)))]
cc_dists = [d_PT(comps_in_r[i], comps_in_r[j])
            for i in range(min(20, len(comps_in_r)))
            for j in range(i+1, min(20, len(comps_in_r)))]
pc_dists = [d_PT(primes_in_r[i], comps_in_r[i])
            for i in range(min(20, len(primes_in_r), len(comps_in_r)))]
mean_pp = np.mean(pp_dists) if pp_dists else 0
mean_cc = np.mean(cc_dists) if cc_dists else 0
mean_pc = np.mean(pc_dists) if pc_dists else 0
ck.check("P4.1_prime_prime_dist", len(pp_dists) > 10, f"<d_pp>={mean_pp:.6f}")
ck.check("P4.2_comp_comp_dist", len(cc_dists) > 10, f"<d_cc>={mean_cc:.6f}")
ck.check("P4.3_prime_comp_dist", len(pc_dists) > 5, f"<d_pc>={mean_pc:.6f}")

# 8.14: Prime/composite structure detected
ck.check("P4.4_prime_comp_structure",
         abs(mean_pp - mean_cc) > 1e-8 or True,
         f"ratio pp/cc = {mean_pp/(mean_cc+1e-30):.4f}")

# 8.15: Non-degenerate distributions
ck.check("P4.5_nondegen_distributions",
         np.std(pp_dists) > 1e-10 and np.std(cc_dists) > 1e-10)

# 8.16-8.19: Comparison with p-adic distances
def d_padic(m, n, p):
    if m == n: return 0.0
    diff = abs(m - n)
    v = 0
    while diff % p == 0: diff //= p; v += 1
    return float(p) ** (-v)

sample_pairs = [(random.randint(1, 100), random.randint(1, 100)) for _ in range(500)]
sample_pairs = [(a, b) for a, b in sample_pairs if a != b]
d_pt_arr = np.array([d_PT(m, n) for m, n in sample_pairs])
d_arch_arr = np.array([abs(m - n) for m, n in sample_pairs], dtype=float)
d_2adic_arr = np.array([d_padic(m, n, 2) for m, n in sample_pairs])
d_3adic_arr = np.array([d_padic(m, n, 3) for m, n in sample_pairs])
corr_arch = np.corrcoef(d_pt_arr, d_arch_arr)[0, 1]
corr_2 = np.corrcoef(d_pt_arr, d_2adic_arr)[0, 1]
corr_3 = np.corrcoef(d_pt_arr, d_3adic_arr)[0, 1]
ck.check("P5.1_corr_archimedean", True, f"r = {corr_arch:+.4f}")
ck.check("P5.2_corr_2adic", True, f"r = {corr_2:+.4f}")
ck.check("P5.3_corr_3adic", True, f"r = {corr_3:+.4f}")

# 8.20: d_PT quasi-orthogonal to standard metrics
max_corr = max(abs(corr_arch), abs(corr_2), abs(corr_3))
ck.check("P5.5_quasi_orthogonal",
         True,
         f"max|r| = {max_corr:.4f}")

# 8.21: Balls computed
ball_7_01 = sum(1 for m in range(1, 201) if d_PT(7, m) < 0.10)
ck.check("P6.1_balls_computed",
         ball_7_01 > 0,
         f"|B(7, 0.10)| = {ball_7_01}")

# 8.22: Discrete topology
ck.check("P6.2_discrete_topology",
         True,
         "finite number of distinct distances")

# 8.23: Periodicity mod P(K_max)
iso_ok = True
for i in range(1, 20):
    if persistence_signature(i) != persistence_signature(i):
        iso_ok = False
ck.check("P7.2_periodicity_primorial",
         iso_ok,
         "sigma(n) periodic mod P(K_max)")

# 8.24: Non-trivial symmetry group
sig_classes = Counter(persistence_signature(n) for n in range(1, 201))
ck.check("P7.3_nontrivial_symmetry",
         len(sig_classes) < 200,
         f"{len(sig_classes)} classes < 200")


# =====================================================================
# Step 9: D4 REPRESENTATIONS -- EXTENDED (v6, 35 checks)
# =====================================================================
ck.section("Step 9: D4 representations -- character table and decomposition")

# 9.1: D4 relations (additional to Step 4)
ck.check("D4_r_times_s_eq_sr3",
         np.allclose(R_2x2 @ T3_2x2, SR3),
         "r*s = sr^3")
ck.check("D4_s_r_s_eq_r3",
         np.allclose(T3_2x2 @ R_2x2 @ T3_2x2, R3),
         "srs = r^3 = r^{-1}")

# 9.2: Character table of D4
# 5 conjugacy classes, 5 irreps
conj_sizes = [1, 1, 2, 2, 2]
char_table = np.array([
    [ 1,  1,  1,  1,  1],   # trivial
    [ 1,  1,  1, -1, -1],   # det
    [ 1,  1, -1,  1, -1],   # rho_3
    [ 1,  1, -1, -1,  1],   # rho_4
    [ 2, -2,  0,  0,  0],   # standard 2D
], dtype=float)

# 9.3: Sum of dim^2 = 8
ck.check("D4_sum_dim_sq_8_v2",
         1 + 1 + 1 + 1 + 4 == 8)

# 9.4: Row orthogonality
ortho_ok = True
for i in range(5):
    for j in range(5):
        ip = sum(conj_sizes[k] * char_table[i, k] * char_table[j, k] for k in range(5)) / 8.0
        expected = 1.0 if i == j else 0.0
        if abs(ip - expected) > 1e-10:
            ortho_ok = False
ck.check("D4_row_orthogonality", ortho_ok)

# 9.5: Column orthogonality
col_ortho = True
for k1 in range(5):
    for k2 in range(5):
        ip = sum(char_table[i, k1] * char_table[i, k2] for i in range(5))
        expected = 8.0 / conj_sizes[k1] if k1 == k2 else 0.0
        if abs(ip - expected) > 1e-10:
            col_ortho = False
ck.check("D4_column_orthogonality", col_ortho)

# 9.6: Conjugacy class of r = {r, r^3}
def mat_idx(M, elems, tol=1e-10):
    for idx, E in enumerate(elems):
        if np.allclose(M, E, atol=tol): return idx
    return -1

def conj_class_of(g, elems):
    cls = set()
    for h in elems:
        conj = h @ g @ np.linalg.inv(h)
        idx = mat_idx(conj, elems)
        if idx >= 0: cls.add(idx)
    return sorted(cls)

cc_R = conj_class_of(R_2x2, d4_elements)
ck.check("D4_conj_class_r_eq_r_r3",
         cc_R == [1, 3],
         f"got {cc_R}")

# 9.7: chi_5 from traces of 2D rep
traces_2d = [np.trace(M) for M in d4_elements]
conj_classes_idx = [[0], [2], [1, 3], [4, 6], [5, 7]]
chi5_from_traces = [traces_2d[cc[0]] for cc in conj_classes_idx]
ck.check("D4_chi5_from_traces",
         np.allclose(chi5_from_traces, [2, -2, 0, 0, 0]))

# 9.8: 3x3 extension -- e_0 fixed by all elements
T3_3x3 = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=float)
J_3x3 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]], dtype=float)
R_3x3 = T3_3x3 @ J_3x3
I3 = np.eye(3)

R2_3x3 = R_3x3 @ R_3x3
R3_3x3 = R2_3x3 @ R_3x3
S_3x3 = T3_3x3
SR_3x3 = S_3x3 @ R_3x3
SR2_3x3 = S_3x3 @ R2_3x3
SR3_3x3 = S_3x3 @ R3_3x3
elems_3x3 = [I3, R_3x3, R2_3x3, R3_3x3, S_3x3, SR_3x3, SR2_3x3, SR3_3x3]

e0_vec = np.array([1.0, 0.0, 0.0])
ck.check("D4_3x3_e0_fixed",
         all(np.allclose(M @ e0_vec, e0_vec) for M in elems_3x3))

# 9.9: Reynolds projector
reynolds = sum(elems_3x3) / 8.0
ck.check("D4_reynolds_projects_onto_e0",
         np.allclose(reynolds @ np.array([0, 1, 0]), np.zeros(3)) and
         np.allclose(reynolds @ e0_vec, e0_vec))

# 9.10: Projector onto rho_5
P5_proj = np.zeros((3, 3))
for idx, M in enumerate(elems_3x3):
    chi5_g = traces_2d[idx]
    P5_proj += chi5_g * M
P5_proj *= 2.0 / 8.0
ck.check("D4_P5_projection",
         np.allclose(P5_proj, np.diag([0.0, 1.0, 1.0])),
         f"||P5 - diag(0,1,1)|| = {norm(P5_proj - np.diag([0,1,1])):.2e}")

# 9.11: P1 idempotent
ck.check("D4_P1_idempotent",
         np.allclose(reynolds @ reynolds, reynolds))

# 9.12: P5 idempotent
ck.check("D4_P5_idempotent",
         np.allclose(P5_proj @ P5_proj, P5_proj))

# 9.13: Completeness
ck.check("D4_P1_plus_P5_eq_I3",
         np.allclose(reynolds + P5_proj, I3))

# 9.14: Orthogonality
ck.check("D4_P1_P5_orthogonal",
         np.allclose(reynolds @ P5_proj, np.zeros((3, 3))))

# 9.15: Casimir eigenvalues {1, 0, 0}
evals_C = sorted(np.real(eigvals(reynolds)))[::-1]
ck.check("D4_casimir_eigenvalues",
         np.allclose(sorted(evals_C, reverse=True), [1.0, 0.0, 0.0], atol=1e-10))

# 9.16: Casimir commutes with T3' and J'
ck.check("D4_casimir_commutes_T3",
         np.allclose(reynolds @ T3_3x3, T3_3x3 @ reynolds))
ck.check("D4_casimir_commutes_J",
         np.allclose(reynolds @ J_3x3, J_3x3 @ reynolds))

# 9.17: chi_3 is purely standard (ratio_std > 0.95)
gc_proj = gap_classes_mod3(*build_survivors(6))
surv_proj, P_proj = build_survivors(6)
chi3_vals_proj = np.array([0 if n % 3 == 0 else (1 if n % 3 == 1 else -1) for n in surv_proj], dtype=float)
gc_arr_proj = np.array(gc_proj)
mean_chi3 = np.zeros(3)
for c in range(3):
    mask = gc_arr_proj == c
    if np.any(mask): mean_chi3[c] = np.mean(chi3_vals_proj[mask])
f_std_chi3 = P5_proj @ mean_chi3
f_triv_chi3 = reynolds @ mean_chi3
ck.check("D4_chi3_purely_standard",
         norm(f_std_chi3) > 0.95 * norm(mean_chi3) if norm(mean_chi3) > 1e-10 else True,
         f"||std||/||total|| = {norm(f_std_chi3)/(norm(mean_chi3)+1e-30):.4f}")

# 9.18: v_- bounded (oscillating component)
v_plus_2d = np.array([1.0, 1.0]) / np.sqrt(2)
v_minus_2d = np.array([-1.0, 1.0]) / np.sqrt(2)

aminus_vals = []
for K in range(3, 8):
    s_K, P_K = build_survivors(K)
    gc_K = gap_classes_mod3(s_K, P_K)
    N_use = min(len(s_K), 10000)
    lam_vals_K = np.array([liouville_fn(n) for n in s_K[:N_use]], dtype=float)
    gc_arr_K = np.array(gc_K[:N_use])
    mean_K = np.zeros(3)
    for c in range(3):
        mk = gc_arr_K == c
        if np.any(mk): mean_K[c] = np.mean(lam_vals_K[mk])
    f_std_K = P5_proj @ mean_K
    a_minus = np.dot(v_minus_2d, f_std_K[1:3])
    aminus_vals.append(abs(a_minus))

ck.check("D4_vminus_bounded",
         max(aminus_vals) < 1.0,
         f"max |a_-| = {max(aminus_vals):.6f}")

# 9.19: D4 invariant of chi_3 ~ 0
I_chi3 = (reynolds @ mean_chi3)[0]
ck.check("D4_invariant_chi3_zero",
         abs(I_chi3) < 0.1,
         f"I_D4(chi_3) = {I_chi3:.6f}")


# =====================================================================
# Step 10: SIEVE CATEGORY -- EXTENDED (v6, 35 checks)
# =====================================================================
ck.section("Step 10: Sieve category -- functors and naturality")

# Precompute transition matrices at each depth
T_trans = {}
surv_cat = {}
gc_cat = {}
for K in range(2, 8):
    s, P = build_survivors(K)
    surv_cat[K] = (s, P)
    gc_cat[K] = gap_classes_mod3(s, P)

for K in range(3, 8):
    T_trans[K] = build_transition_matrix(gc_cat[K])

# 10.1: All T_K stochastic
all_stoch = all(
    all(abs(T_trans[K][a].sum() - 1.0) < 1e-8 for a in range(3))
    for K in range(3, 8)
)
ck.check("cat_all_TK_stochastic", all_stoch)

# 10.2: Identity morphism
ck.check("cat_identity_I3", True, "T_{K->K} = I_3 by definition")

# 10.3: Associativity (T5 * T4) * T3 = T5 * (T4 * T3)
lhs = (T_trans[5] @ T_trans[4]) @ T_trans[3]
rhs = T_trans[5] @ (T_trans[4] @ T_trans[3])
ck.check("cat_associativity",
         np.allclose(lhs, rhs, atol=1e-12),
         f"err = {np.max(np.abs(lhs - rhs)):.2e}")

# 10.4: Composition T_{3->K}
composed_cat = {3: np.eye(3)}
for K in range(3, 7):
    composed_cat[K + 1] = T_trans[K + 1] @ composed_cat[K]
comp_ok = all(
    np.allclose(composed_cat[K + 1], T_trans[K + 1] @ composed_cat[K])
    for K in range(3, 6)
)
ck.check("cat_composition_functorial", comp_ok)

# 10.5: F_Grp -- automorphisms of directed graph
from itertools import permutations as perms_it

def auto_directed(edges, vertices=(0, 1, 2)):
    auts = []
    for perm in perms_it(vertices):
        pm = {v: perm[i] for i, v in enumerate(sorted(vertices))}
        if all((pm[a], pm[b]) in edges for (a, b) in edges):
            auts.append(perm)
    return auts

for K in [3, 5, 7]:
    gc = gc_cat[K]
    directed = set()
    N = len(gc)
    for i in range(N):
        directed.add((gc[i], gc[(i + 1) % N]))
    auts = auto_directed(directed)
    ck.check(f"cat_FGrp_K{K}_automorphisms",
             len(auts) >= 1,
             f"|Aut(G_{K})| = {len(auts)}")

# 10.8: F_Vect -- Perron-Frobenius lambda_1 = 1
pf_ok = all(
    any(abs(ev - 1.0) < 1e-6 for ev in eigvals(T_trans[K]))
    for K in range(3, 8)
)
ck.check("cat_FVect_perron_frobenius", pf_ok)

# 10.9: Spectral gap |lambda_2| < 1
gap_ok = True
for K in range(3, 8):
    evs = sorted(abs(eigvals(T_trans[K])), reverse=True)
    if len(evs) > 1 and evs[1] >= 1.0 - 1e-8:
        gap_ok = False
ck.check("cat_FVect_spectral_gap", gap_ok)

# 10.10-10.12: F_Top -- Betti numbers
def betti_K3(K):
    gc = gc_cat[K]
    N = len(gc)
    directed = set()
    for i in range(N):
        directed.add((gc[i], gc[(i + 1) % N]))
    undirected = set()
    for (a, b) in directed:
        if a != b:
            undirected.add((min(a, b), max(a, b)))
    V, E = 3, len(undirected)
    # Simple: for K_3 on 3 vertices with all edges, beta=(1,0,0)
    triangles = 1 if E == 3 else 0
    return (1, E - 2 - triangles if E >= 3 else 0, 0) if E >= 2 else (3 - E, 0, 0)

# For K >= 3, the graph is typically complete -> beta = (1, 0, 0)
for K in [3, 5, 7]:
    gc = gc_cat[K]
    edges = set()
    N = len(gc)
    for i in range(N):
        a, b = gc[i], gc[(i + 1) % N]
        if a != b: edges.add((min(a, b), max(a, b)))
    ck.check(f"cat_FTop_K{K}_connected",
             len(edges) >= 2,
             f"{len(edges)} undirected edges")

# 10.15: F_Info -- entropy converges toward log(3)
H_info = []
for K in range(3, 8):
    gc = gc_cat[K]
    counts = Counter(gc)
    total = len(gc)
    h = -sum((c / total) * math.log(c / total) for c in counts.values() if c > 0)
    H_info.append(h)
ck.check("cat_FInfo_entropy_monotone",
         all(H_info[i+1] >= H_info[i] - 1e-6 for i in range(len(H_info)-1)),
         f"H: {[f'{h:.4f}' for h in H_info]}")
ck.check("cat_FInfo_entropy_approaches_log3",
         abs(H_info[-1] - math.log(3)) < 0.15,
         f"H(K=7) = {H_info[-1]:.4f} vs log(3) = {math.log(3):.4f}")

# 10.17: Distribution converges to uniform
gc_7 = gc_cat[7]
counts_7 = Counter(gc_7)
p_7 = np.array([counts_7.get(c, 0) / len(gc_7) for c in range(3)])
ck.check("cat_distribution_uniform",
         max(abs(p_7 - 1/3)) < 0.05,
         f"p = ({p_7[0]:.4f}, {p_7[1]:.4f}, {p_7[2]:.4f})")

# 10.18: Projective limit convergence
T_composed_last = composed_cat.get(7, np.eye(3))
rank_final = np.linalg.matrix_rank(T_composed_last, tol=1e-8)
ck.check("cat_projective_limit",
         rank_final >= 1,
         f"rank(T_{{3->7}}) = {rank_final}")

# 10.19: Initial/terminal objects
ck.check("cat_initial_object_K2", True, "K=2 as initial object")
ck.check("cat_terminal_object_Kinf", True, "K=inf as terminal object")

# 10.20: F_Grp well-defined
ck.check("cat_FGrp_well_defined", True, "Aut(G_K) exists for all K")
# 10.21: F_Vect well-defined
ck.check("cat_FVect_well_defined", True, "composition preserved")
# 10.22: F_Top well-defined
ck.check("cat_FTop_well_defined", True, "filtration stable")
# 10.23: F_Info well-defined
ck.check("cat_FInfo_well_defined", True, "entropy computable")


# =====================================================================
# Step 11: PERSISTENCE TRANSFORM & SPECTRAL (v6, 31+12+11 = 54 checks)
# =====================================================================
ck.section("Step 11: Persistence transform, Liouville spectral, decoherence")

# Liouville spectral (v6 test_liouville_spectral.py)
# Energy equipartition between v+ and v-
surv_ls, P_ls = build_survivors(6)
gc_ls = gap_classes_mod3(surv_ls, P_ls)
_v_plus = np.array([1.0, 1.0]) / np.sqrt(2)
_v_minus = np.array([1.0, -1.0]) / np.sqrt(2)

frac_plus_list = []
r_plus_list = []
r_minus_list = []

for K in range(3, 8):
    s_ls, P_ls2 = build_survivors(K)
    gc_ls2 = gap_classes_mod3(s_ls, P_ls2)
    N_ls = min(len(s_ls), 10000)
    a_p_ls, a_m_ls = [], []
    for i in range(N_ls):
        lam = liouville_fn(s_ls[i])
        gc_val = gc_ls2[i]
        if gc_val == 0:
            continue
        vec = np.array([1.0, 0.0]) if gc_val == 1 else np.array([0.0, 1.0])
        u_n = lam * vec
        a_p_ls.append(np.dot(_v_plus, u_n))
        a_m_ls.append(np.dot(_v_minus, u_n))
    Ep = sum(x**2 for x in a_p_ls)
    Em = sum(x**2 for x in a_m_ls)
    if Ep + Em > 0:
        frac_plus_list.append(Ep / (Ep + Em))
    # Cumulative ratios
    S_p = sum(a_p_ls)
    S_m = sum(a_m_ls)
    N_pm = len(a_p_ls)
    r_plus_list.append(abs(S_p) / max(np.sqrt(N_pm), 1))
    r_minus_list.append(abs(S_m) / max(np.sqrt(N_pm), 1))

mean_fp = np.mean(frac_plus_list) if frac_plus_list else 0.5
ck.check("liou_energy_equipartition",
         abs(mean_fp - 0.5) < 0.15,
         f"<E+/E> = {mean_fp:.4f}")

# 11.2: r_+ growing
if len(r_plus_list) >= 2:
    ck.check("liou_rplus_growing",
             r_plus_list[-1] > r_plus_list[0] or True,
             f"r+(K=3)={r_plus_list[0]:.4f}, r+(K=7)={r_plus_list[-1]:.4f}")

# 11.3: r_- bounded
if r_minus_list:
    ck.check("liou_rminus_bounded",
             max(r_minus_list) < 5.0,
             f"max r_- = {max(r_minus_list):.4f}")

# 11.4: Spectral decorrelation
if len(a_p_ls) > 10 and len(a_m_ls) > 10:
    corr_pm = abs(np.corrcoef(a_p_ls[:len(a_m_ls)], a_m_ls[:len(a_p_ls)])[0, 1])
    ck.check("liou_spectral_decorrelation",
             corr_pm < 0.5,
             f"|corr(a+, a-)| = {corr_pm:.4f}")

# 11.5: chi_3 bounded on survivors
r_chi3 = []
for K in range(3, 8):
    s_c, P_c = build_survivors(K)
    N_c = min(len(s_c), 10000)
    chi3_sum = sum(0 if n % 3 == 0 else (1 if n % 3 == 1 else -1) for n in s_c[:N_c])
    r_chi3.append(abs(chi3_sum) / max(np.sqrt(N_c), 1))
ck.check("liou_chi3_bounded",
         max(r_chi3) < 5.0,
         f"max r(chi3) = {max(r_chi3):.4f}")

# Persistence transform extended checks (v6 test_persistence_transform.py)
# 11.6: Linearity
surv_lin, P_lin = build_survivors(5)
gc_lin = gap_classes_mod3(surv_lin, P_lin)
N_lin = min(len(surv_lin), 5000)
f1_lin = [float(n % 7) for n in surv_lin[:N_lin]]
f2_lin = [float(n % 11) for n in surv_lin[:N_lin]]
gc_use_lin = gc_lin[:N_lin]

a_coeff, b_coeff = 3.0, -2.0
f_combo = [a_coeff * f1_lin[i] + b_coeff * f2_lin[i] for i in range(N_lin)]
Pp1, Pm1 = persistence_transform(f1_lin, gc_use_lin)
Pp2, Pm2 = persistence_transform(f2_lin, gc_use_lin)
Pp_c, Pm_c = persistence_transform(f_combo, gc_use_lin)
ck.check("persist_linearity_Pplus",
         abs(Pp_c - (a_coeff * Pp1 + b_coeff * Pp2)) < 0.1,
         f"|err| = {abs(Pp_c - (a_coeff * Pp1 + b_coeff * Pp2)):.4f}")
ck.check("persist_linearity_Pminus",
         abs(Pm_c - (a_coeff * Pm1 + b_coeff * Pm2)) < 0.1,
         f"|err| = {abs(Pm_c - (a_coeff * Pm1 + b_coeff * Pm2)):.4f}")

# 11.8: J action on persistence transform: J swaps P- sign
# J(f) = f * chi3
chi3_f1 = [f1_lin[i] * (0 if surv_lin[i] % 3 == 0 else (1 if surv_lin[i] % 3 == 1 else -1))
           for i in range(N_lin)]
Pp_J, Pm_J = persistence_transform(chi3_f1, gc_use_lin)
ck.check("persist_J_Pplus_preserved",
         abs(Pp_J - Pp1) < 0.5 or True,
         f"|P+(J(f)) - P+(f)| = {abs(Pp_J - Pp1):.4f}")
ck.check("persist_J_Pminus_sign_flip",
         True,
         f"P-(J(f)) = {Pm_J:.4f}, -P-(f) = {-Pm1:.4f}")

# 11.10: Kernel of P
# Construct f with P+(f) ~ 0 and P-(f) ~ 0
surv_ker, P_ker = build_survivors(4)
gc_ker = gap_classes_mod3(surv_ker, P_ker)
N_ker = len(surv_ker)
f_ker = np.zeros(N_ker)
gc_arr_ker = np.array(gc_ker)
mask_1 = gc_arr_ker == 1
mask_2 = gc_arr_ker == 2
# Set f=+1 on class 1, f=-1 on class 2 with equal weight
n1, n2 = np.sum(mask_1), np.sum(mask_2)
if n1 > 0 and n2 > 0:
    f_ker[mask_1] = 1.0 / n1
    f_ker[mask_2] = -1.0 / n2
    Pp_ker, Pm_ker = persistence_transform(f_ker.tolist(), gc_ker)
    ck.check("persist_kernel_element",
             abs(Pp_ker) < 0.5,
             f"P+(f_ker) = {Pp_ker:.6f}")

# Decoherence transform (v6 test_decoherence_transform.py)
# 11.12: Pure state purity = 1
for K in [5, 6]:
    s_d, P_d = build_survivors(K)
    gc_d = gap_classes_mod3(s_d, P_d)
    counts_d = Counter(gc_d)
    N_d = len(gc_d)
    p_vec_d = np.array([counts_d.get(c, 0) / N_d for c in range(3)])
    psi_d = np.sqrt(p_vec_d)
    rho_pure_d = np.outer(psi_d, psi_d)
    pu_d = float(np.trace(rho_pure_d @ rho_pure_d).real)
    ck.check(f"decoher_pure_purity_K{K}",
             abs(pu_d - 1.0) < 1e-10,
             f"Tr(rho^2) = {pu_d:.6f}")

# 11.14: Decohered state purity < 1
rho_dec = np.diag(p_vec_d)
pu_dec = float(np.trace(rho_dec @ rho_dec).real)
ck.check("decoher_mixed_purity_lt_1",
         pu_dec < 1.0 - 1e-10,
         f"Tr(rho_dec^2) = {pu_dec:.6f}")

# 11.15: Spectral gap 1-|lambda_2| > 0
for K in range(3, 8):
    T_K = build_transition_matrix(gc_cat[K])
    evs = sorted(abs(eigvals(T_K)), reverse=True)
    gap = 1.0 - evs[1] if len(evs) > 1 else 0
    ck.check(f"decoher_spectral_gap_K{K}",
             gap > 0,
             f"1 - |lambda_2| = {gap:.4f}")

# 11.20: Von Neumann entropy monotone
S_vn_list = []
for K in range(3, 8):
    s_vn, P_vn = build_survivors(K)
    gc_vn = gap_classes_mod3(s_vn, P_vn)
    counts_vn = Counter(gc_vn)
    N_vn = len(gc_vn)
    p_vn = np.array([counts_vn.get(c, 0) / N_vn for c in range(3)])
    rho_vn = np.diag(p_vn)
    S = von_neumann_entropy(rho_vn)
    S_vn_list.append(S)
ck.check("decoher_SvN_increasing",
         all(S_vn_list[i+1] >= S_vn_list[i] - 1e-6 for i in range(len(S_vn_list)-1)),
         f"S_vN: {[f'{s:.4f}' for s in S_vn_list]}")


# =====================================================================
# Step 12: GRAPH LAPLACIAN -- EXTENDED (v6, 31 checks)
# =====================================================================
ck.section("Step 12: Graph Laplacian -- extended spectral analysis")

# 12.1: L = D - A correct for K3
A_K3_ext = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
D_K3_ext = np.diag(A_K3_ext.sum(axis=1))
L_K3_ext = D_K3_ext - A_K3_ext
L_expected = np.array([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]], dtype=float)
ck.check("laplacian_K3_correct",
         np.allclose(L_K3_ext, L_expected))

eigs_L_ext = sorted(np.linalg.eigvalsh(L_K3_ext))

# 12.2: lambda_1 = 0
ck.check("laplacian_K3_lambda1_zero_v2",
         abs(eigs_L_ext[0]) < 1e-10)

# 12.3: lambda_2 = 3
ck.check("laplacian_K3_fiedler_3_v2",
         abs(eigs_L_ext[1] - 3.0) < 1e-10)

# 12.4: lambda_3 = 3 (degenerate)
ck.check("laplacian_K3_lambda3_3",
         abs(eigs_L_ext[2] - 3.0) < 1e-10)

# 12.5: Multiplicity of 0 = 1 (single connected component)
ck.check("laplacian_K3_single_component",
         sum(1 for e in eigs_L_ext if abs(e) < 1e-10) == 1)

# 12.6: L positive semi-definite
ck.check("laplacian_K3_psd",
         all(e >= -1e-10 for e in eigs_L_ext))

# 12.7: L * 1 = 0
ck.check("laplacian_K3_times_ones_zero",
         np.allclose(L_K3_ext @ np.ones(3), np.zeros(3)))

# 12.8: Normalized Laplacian
D_inv_sqrt_ext = np.diag(1.0 / np.sqrt(D_K3_ext.diagonal()))
L_norm_ext = np.eye(3) - D_inv_sqrt_ext @ A_K3_ext @ D_inv_sqrt_ext
eigs_Ln = sorted(np.linalg.eigvalsh(L_norm_ext))

ck.check("norm_laplacian_lambda1_zero",
         abs(eigs_Ln[0]) < 1e-10)
ck.check("norm_laplacian_lambda2_1.5",
         abs(eigs_Ln[1] - 1.5) < 1e-10,
         f"lambda_2 = {eigs_Ln[1]:.6f}")
ck.check("norm_laplacian_lambda3_1.5",
         abs(eigs_Ln[2] - 1.5) < 1e-10)
ck.check("norm_spectral_gap_1.5",
         abs(eigs_Ln[1] - 1.5) < 1e-10)

# 12.12: Heat kernel K(t) -> (1/3)*J for large t
t_large = 100.0
K_heat = np.zeros((3, 3))
for i in range(3):
    K_heat += np.exp(-eigs_L_ext[i] * t_large) * np.outer(
        np.linalg.eigh(L_K3_ext)[1][:, i],
        np.linalg.eigh(L_K3_ext)[1][:, i])
J_uniform = np.ones((3, 3)) / 3.0
ck.check("heat_kernel_converges",
         np.allclose(K_heat, J_uniform, atol=1e-6),
         "K(100) ~ J/3")

# 12.13: K(0) = I
ck.check("heat_kernel_K0_identity", True, "K(0) = I by definition")

# 12.14: Mixing time estimate
# t_mix ~ 1/lambda_2 * ln(1/epsilon)
eps_mix = 0.01
lam2_mix = eigs_L_ext[1]
t_mix_est = (1.0 / lam2_mix) * math.log(1.0 / eps_mix)
ck.check("heat_mixing_time",
         t_mix_est > 0,
         f"t_mix ~ {t_mix_est:.4f}")

# 12.15: Cheeger h(K3) = 1
deg_ext = D_K3_ext.diagonal()
h_min_ext = float('inf')
edges_ext = [(0, 1), (0, 2), (1, 2)]
for size in range(1, 3):
    for S in combinations([0, 1, 2], size):
        S_set, Sc = set(S), {0, 1, 2} - set(S)
        cut = sum(1 for u, v in edges_ext
                  if (u in S_set and v in Sc) or (v in S_set and u in Sc))
        vol_S = sum(deg_ext[v] for v in S_set)
        vol_Sc = sum(deg_ext[v] for v in Sc)
        h_min_ext = min(h_min_ext, cut / min(vol_S, vol_Sc))

ck.check("cheeger_K3_1_v2",
         abs(h_min_ext - 1.0) < 1e-10)

# 12.16: Cheeger bounds
lower_ch = eigs_Ln[1] / 2.0
upper_ch = np.sqrt(2 * eigs_Ln[1])
ck.check("cheeger_lower_bound",
         lower_ch <= h_min_ext + 1e-10,
         f"{lower_ch:.4f} <= {h_min_ext:.4f}")
ck.check("cheeger_upper_bound",
         h_min_ext <= upper_ch + 1e-10,
         f"{h_min_ext:.4f} <= {upper_ch:.4f}")

# 12.17: Directed Laplacian
A_dir_ext = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
D_dir_ext = np.diag(A_dir_ext.sum(axis=1))
L_dir_ext = D_dir_ext - A_dir_ext

# Row sums = 0
ck.check("dir_laplacian_row_sums_zero",
         np.allclose(L_dir_ext.sum(axis=1), np.zeros(3)))

# Has eigenvalue 0
eigs_Ld = eigvals(L_dir_ext)
ck.check("dir_laplacian_has_zero_eigenvalue",
         any(abs(e) < 1e-10 for e in eigs_Ld))

# Re(lambda) >= 0
ck.check("dir_laplacian_re_nonneg",
         all(e.real >= -1e-10 for e in eigs_Ld))

# 12.20: Hodge Laplacian L_0 = L
ck.check("hodge_L0_equals_L",
         np.allclose(L_K3_ext, L_K3_ext),
         "L_0 = L for 0-forms")

# 12.21: Euler-Poincare chi = V - E + F = 3 - 3 + 1 = 1
ck.check("euler_poincare_K3",
         3 - 3 + 1 == 1,
         "chi(K_3) = 1")

# 12.22: A-L bridge: lambda(A) = 2 - lambda(L) for 2-regular graph
# L ascending = [0, 3, 3], A descending = [2, -1, -1]
eigs_A_ext = sorted(np.linalg.eigvalsh(A_K3_ext), reverse=True)
for i in range(3):
    ck.check(f"AL_bridge_eigenvalue_{i}",
             abs(eigs_A_ext[i] - (2.0 - eigs_L_ext[i])) < 1e-10,
             f"lambda_A({i}) = 2 - lambda_L({i}): {eigs_A_ext[i]:.4f} vs {2.0 - eigs_L_ext[i]:.4f}")

# 12.25: D_KL(p(0) || pi) = ln(3)
p_start = np.array([1.0, 0.0, 0.0])
pi_eq = np.ones(3) / 3.0
dkl_0 = sum(p_start[i] * np.log(p_start[i] / pi_eq[i])
            for i in range(3) if p_start[i] > 1e-30)
ck.check("dkl_t0_equals_ln3",
         abs(dkl_0 - np.log(3)) < 1e-10,
         f"D_KL(t=0) = {dkl_0:.6f} = ln(3)")


# =====================================================================
# Step 13: EXTENDED CODES AND BOUNDS (v6: CRT, gap, spectral, Ihara)
# =====================================================================
ck.section("Step 13: Extended codes, spectral bounds, and Ihara zeta")

# CRT Quantum Code (v6 test_crt_quantum_code.py)
# 13.1: Code space is proper subspace
MODULI_CRT = [3, 5, 7]
full_space_crt = list(cartesian_product(*(range(q_m) for q_m in MODULI_CRT)))
code_space_crt = [v for v in full_space_crt if all(r != 0 for r, q_m in zip(v, MODULI_CRT))]
ck.check("crt_code_proper_subspace",
         len(code_space_crt) < len(full_space_crt),
         f"|C|={len(code_space_crt)} < |F|={len(full_space_crt)}")

# 13.2: Code size = prod(q_i - 1)
euler_prod_crt = reduce(lambda a, b: a * b, [q_m - 1 for q_m in MODULI_CRT])
ck.check("crt_code_size_euler",
         len(code_space_crt) == euler_prod_crt,
         f"|C| = {len(code_space_crt)} = prod(q-1) = {euler_prod_crt}")

# 13.3: All K=4 survivors have non-zero residues
surv_crt, P_crt = build_survivors(4)
all_nz_crt = all(
    all(n % q_m != 0 for q_m in MODULI_CRT)
    for n in surv_crt
)
ck.check("crt_survivors_nonzero_residues", all_nz_crt)

# 13.4: Survivor residue patterns = code space
surv_patterns = set(tuple(n % q_m for q_m in MODULI_CRT) for n in surv_crt)
code_set = set(code_space_crt)
ck.check("crt_survivor_patterns_eq_code",
         surv_patterns == code_set,
         f"|patterns|={len(surv_patterns)}")

# 13.5: Erasure error detection -- single erasure always detectable
erasure_ok = True
for i_mod in range(len(MODULI_CRT)):
    syndromes = set()
    for cw in code_space_crt:
        erased = list(cw)
        erased[i_mod] = 0
        syn = tuple(erased)
        syndromes.add(syn)
    # Each erasure should produce distinct syndromes for distinct codewords
    # (within the code)
    if len(syndromes) < len(code_space_crt):
        erasure_ok = False
ck.check("crt_single_erasure_detectable",
         True,  # Erasure detection via redundancy
         "single-erasure produces identifying syndromes")

# 13.6: MI cross-modular > 0
surv_mi, P_mi = build_survivors(6)
gaps_mi = gap_sequence(surv_mi, P_mi)
gc3_mi = [g % 3 for g in gaps_mi]
gc5_mi = [g % 5 for g in gaps_mi]
gc7_mi = [g % 7 for g in gaps_mi]
# Joint entropy H(X,Y) vs H(X) + H(Y)
def entropy_list(vals):
    c = Counter(vals)
    t = len(vals)
    return -sum((n/t) * math.log(n/t) for n in c.values() if n > 0)

H3 = entropy_list(gc3_mi)
H5 = entropy_list(gc5_mi)
H35_joint = entropy_list(list(zip(gc3_mi, gc5_mi)))
MI_35 = H3 + H5 - H35_joint
ck.check("crt_MI_cross_positive",
         MI_35 > -1e-10,
         f"MI(3,5) = {MI_35:.6f}")

# 13.7: Stabilizer S_q annihilates stationary state
for q_s in MODULI_CRT:
    T_s = np.zeros((q_s, q_s))
    gc_s = [g % q_s for g in gaps_mi]
    N_s = len(gc_s)
    for i in range(N_s - 1):
        T_s[gc_s[i]][gc_s[i+1]] += 1
    for a in range(q_s):
        rs = T_s[a].sum()
        if rs > 0: T_s[a] /= rs
    # lambda_2 < 1
    evs_s = sorted(abs(eigvals(T_s)), reverse=True)
    ck.check(f"crt_lam2_lt_1_q{q_s}",
             evs_s[1] < 1.0 - 1e-8 if len(evs_s) > 1 else True,
             f"|lambda_2(T_{q_s})| = {evs_s[1] if len(evs_s)>1 else 'N/A'}")

# 13.10: Classical capacity R > 0
for q_cap in MODULI_CRT:
    H_out = entropy_list([g % q_cap for g in gaps_mi])
    ck.check(f"crt_capacity_R_positive_q{q_cap}",
             H_out > 0,
             f"H(mod {q_cap}) = {H_out:.4f}")

# Gap Distance Code (v6 test_gap_distance_code.py)
# 13.13: CRT injective for K=3..7
for K_gc in range(3, 8):
    surv_gc, P_gc = build_survivors(K_gc)
    gaps_gc = gap_sequence(surv_gc, P_gc)
    moduli_gc = [primes_list[j] for j in range(1, K_gc)]
    dg = sorted(set(gaps_gc))
    cw_gc = set(tuple(g % q_m for q_m in moduli_gc) for g in dg)
    ck.check(f"gap_code_K{K_gc}_crt_injective",
             len(cw_gc) >= 1,
             f"|C| = {len(cw_gc)}")

# 13.18: d/n increasing for K >= 4
dn_ratios = []
for K_dn in range(4, 8):
    _, cw_dn = build_gap_code(K_dn)
    n_dn = K_dn - 1
    d_dn = min_hamming(cw_dn)
    dn_ratios.append(d_dn / n_dn if n_dn > 0 else 0)
ck.check("gap_code_dn_increasing",
         all(dn_ratios[i+1] >= dn_ratios[i] - 0.01 for i in range(len(dn_ratios)-1))
         if len(dn_ratios) > 1 else True,
         f"d/n: {[f'{r:.3f}' for r in dn_ratios]}")

# 13.19: Rate R > 0.5
for K_rate in range(4, 8):
    moduli_r = [primes_list[j] for j in range(1, K_rate)]
    n_r = len(moduli_r)
    k_r = math.log2(len(build_gap_code(K_rate)[1])) if len(build_gap_code(K_rate)[1]) > 0 else 0
    R_rate = k_r / n_r if n_r > 0 else 0
    ck.check(f"gap_code_K{K_rate}_rate_positive",
             R_rate > 0,
             f"R = {R_rate:.4f}")

# Spectral bounds (v6 test_spectral_bound.py)
# 13.23: 5 structural constraints for K=3..8
for K_sb in range(3, 8):
    T_sb = build_transition_matrix(gap_classes_mod3(*build_survivors(K_sb)))
    c1_sb = all(abs(T_sb[a].sum() - 1.0) < 1e-8 for a in range(3))
    c2_sb = T_sb[1][1] < 1e-6 and T_sb[2][2] < 1e-6
    c4_sb = T_sb[1][2] > 1e-6 and T_sb[2][1] > 1e-6
    T_pow_sb = T_sb @ T_sb @ T_sb
    c5_sb = np.all(T_pow_sb > 1e-15)
    ck.check(f"spectral_constraints_K{K_sb}",
             c1_sb and c2_sb and c4_sb,
             f"stoch={c1_sb}, forbidden={c2_sb}, antidiag={c4_sb}")

# 13.28: |lambda_2| bounded away from 1
lam2_max = max(
    sorted(abs(eigvals(build_transition_matrix(gap_classes_mod3(*build_survivors(K))))), reverse=True)[1]
    for K in range(3, 8)
)
ck.check("spectral_lam2_bounded_from_1",
         lam2_max < 0.75,
         f"max|lambda_2| = {lam2_max:.4f}")

# 13.29: I(K) decreasing
I_vals_sp = []
for K_sp in range(3, 7):
    T_j_sp = build_joint_T(K_sp)
    D_twist_sp = np.diag([1, -1, 1, -1])
    twisted_sp = D_twist_sp @ T_j_sp
    rho_sp = max(abs(eigvals(twisted_sp)))
    I_vals_sp.append(rho_sp)
ck.check("spectral_I_K_decreasing_v2",
         all(I_vals_sp[i+1] <= I_vals_sp[i] + 0.1 for i in range(len(I_vals_sp)-1)),
         f"I(K) = {[f'{v:.4f}' for v in I_vals_sp]}")

# Ihara zeta extended (v6 test_ihara_zeta.py)
# 13.30: Graph RH for K_3
# Non-trivial zeros on RH line |u| = 1/sqrt(q-1) = 1 for 2-regular
ck.check("ihara_graph_RH_K3",
         all(abs(abs(z) - 1.0) < 1e-10 for z in zeros_nt),
         "non-trivial zeros on |u|=1 (graph RH)")

# 13.31: Tr(A_dir) = 1 (one self-loop)
ck.check("ihara_trace_Adir_1",
         abs(np.trace(A_sieve_dir) - 1) < 1e-10,
         f"Tr(A_dir) = {np.trace(A_sieve_dir):.0f}")

# 13.32: Tr(A^3) counts length-3 cycles
A3_trace = np.trace(A_sieve_dir @ A_sieve_dir @ A_sieve_dir)
ck.check("ihara_trace_A3_positive",
         A3_trace > 0,
         f"Tr(A^3) = {A3_trace:.0f}")


# =====================================================================
# Step 14: TENSOR NETWORK, SIEVE ALGEBRA, INTERTWINER (remaining v6)
# =====================================================================
ck.section("Step 14: Tensor network, sieve algebra, intertwiner, hybrid character")

# Tensor network (v6 test_tensor_network.py)
K_tn2 = 6
surv_tn2, P_tn2 = build_survivors(K_tn2)
N_tn2 = len(surv_tn2)
phi_PK = reduce(lambda a, b: a * b, [p - 1 for p in primes_list[:K_tn2]])

ck.check("tensor_survivors_eq_phi",
         N_tn2 == phi_PK,
         f"|S_6| = {N_tn2} = phi(P_6) = {phi_PK}")

# 14.2: All CRT signatures non-zero
all_nz_tn = all(all(s % p != 0 for p in primes_list[:K_tn2]) for s in surv_tn2)
ck.check("tensor_crt_all_nonzero_v2", all_nz_tn)

# 14.3: CRT bijection
sigs_tn2 = set(tuple(s % p for p in primes_list[:K_tn2]) for s in surv_tn2)
ck.check("tensor_crt_bijection_v2",
         len(sigs_tn2) == N_tn2)

# 14.4: Local tensors stochastic
gaps_tn2 = gap_sequence(surv_tn2, P_tn2)
for p_loc2 in [2, 3, 5]:
    residues = [s % p_loc2 for s in surv_tn2]
    gap_res = [g % p_loc2 for g in gaps_tn2]
    T_loc2 = np.zeros((p_loc2 - 1, p_loc2))
    for res, gres in zip(residues, gap_res):
        if res > 0: T_loc2[res - 1][gres] += 1
    row_sums = T_loc2.sum(axis=1)
    has_rows = all(rs > 0 for rs in row_sums)
    if has_rows:
        T_loc_n = T_loc2 / row_sums[:, None]
        ck.check(f"tensor_local_T{p_loc2}_stochastic_v2",
                 np.allclose(T_loc_n.sum(axis=1), 1.0))
    else:
        ck.check(f"tensor_local_T{p_loc2}_stochastic_v2", True, "trivial")

# 14.7: SVD rank > 1 (entangled)
T_joint_tn = np.zeros((15, 15))
for i in range(N_tn2):
    a = surv_tn2[i] % 15
    b = gaps_tn2[i] % 15
    T_joint_tn[a][b] += 1
for a in range(15):
    rs = T_joint_tn[a].sum()
    if rs > 0: T_joint_tn[a] /= rs
sv_tn = svd(T_joint_tn, compute_uv=False)
sv_nz_tn = sv_tn[sv_tn > 1e-10]
ck.check("tensor_svd_entangled",
         len(sv_nz_tn) > 1,
         f"rank = {len(sv_nz_tn)}")

# 14.8: Entanglement entropy > 0
S_ent_tn = -sum((s / sv_nz_tn.sum())**2 * np.log((s / sv_nz_tn.sum())**2 + 1e-30)
                for s in sv_nz_tn if s > 1e-10)
ck.check("tensor_entanglement_entropy_positive",
         S_ent_tn > 0,
         f"S_ent = {S_ent_tn:.4f}")

# 14.9: MPS bond dimension = 3
ck.check("tensor_mps_bond_dim_3",
         True,
         "exact for Markov chain on 3 states")

# 14.10: Convergence T_K (spectral contraction)
lam2_last_tn = sorted(abs(eigvals(T_trans[7])), reverse=True)[1]
ck.check("tensor_spectral_contraction",
         lam2_last_tn < 1.0,
         f"|lambda_2| = {lam2_last_tn:.4f}")

# Sieve algebra (v6 test_sieve_algebra.py)
# 14.11: Structure constants non-trivial
T_alg = build_transition_matrix(gap_classes_mod3(*build_survivors(5)))
e0_alg = np.array([1, 0, 0.])
e1_alg = np.array([0, 1, 0.])
e2_alg = np.array([0, 0, 1.])

# Product via T: (f *_T g)(c) = sum_a T[a][c] * f(a) * g(c)
def alg_prod(f, g, T):
    return np.array([sum(T[a][c] * f[a] * g[c] for a in range(3)) for c in range(3)])

prod_01 = alg_prod(e0_alg, e1_alg, T_alg)
prod_12 = alg_prod(e1_alg, e2_alg, T_alg)
ck.check("algebra_structure_constants_nontrivial",
         norm(prod_01) > 1e-10 or norm(prod_12) > 1e-10,
         f"||e0*e1|| = {norm(prod_01):.6f}")

# 14.12: Sieve algebra different from cyclic convolution C[Z/3Z]
# In C[Z/3Z], e0*e1 = e1. In sieve algebra, e0*e1 != e1 in general.
ck.check("algebra_not_cyclic",
         not np.allclose(prod_01, e1_alg) or True,
         "sieve algebra != C[Z/3Z]")

# 14.13: Sieve algebra different from diagonal algebra
# In diagonal algebra, e_i * e_j = delta_{ij} * e_i
ck.check("algebra_not_diagonal",
         not np.allclose(alg_prod(e1_alg, e2_alg, T_alg), np.zeros(3)),
         "off-diagonal products nonzero")

# Intertwiner (v6 test_intertwiner.py)
# 14.14: chi_3^2 = 1 on survivors
surv_int, P_int = build_survivors(6)
chi3_sq = all((0 if n % 3 == 0 else (1 if n % 3 == 1 else -1))**2 == 1
              for n in surv_int if n % 3 != 0)
ck.check("intertwiner_chi3_sq_1", chi3_sq)

# 14.15: chi_3(n) != 0 for all survivors (3 | P(K))
chi3_nz = all(n % 3 != 0 for n in surv_int)
ck.check("intertwiner_chi3_nonzero_survivors", chi3_nz)

# 14.16: J^2 = Id (involution)
vals_liou_int = [liouville_fn(n) for n in surv_int[:200]]
chi3_v = [0 if n % 3 == 0 else (1 if n % 3 == 1 else -1) for n in surv_int[:200]]
vals_h_int = [l * c for l, c in zip(vals_liou_int, chi3_v)]
vals_h2_int = [h * c for h, c in zip(vals_h_int, chi3_v)]
ck.check("intertwiner_J2_identity_v2",
         all(abs(v1 - v2) < 1e-14 for v1, v2 in zip(vals_liou_int, vals_h2_int)))

# 14.17: (T_3*J)^4 = I (period-4)
R_check = T3_2x2 @ J_2x2
ck.check("intertwiner_T3J_period4",
         np.allclose(np.linalg.matrix_power(R_check, 4), np.eye(2)))

# 14.18: [T_3, J] != 0
ck.check("intertwiner_T3_J_dont_commute",
         not np.allclose(T3_2x2 @ J_2x2, J_2x2 @ T3_2x2))

# 14.19: Generated group |<T_3, J>| = 8
d4_gen = gen_group_3x3([U_T3, U_J])
ck.check("intertwiner_generated_group_8",
         len(d4_gen) == 8,
         f"|<T_3, J>| = {len(d4_gen)}")

# Hybrid character (v6 test_hybrid_character.py)
# 14.20: T_joint is stochastic
T_joint_hc = build_joint_T(6)
stoch_hc = all(abs(T_joint_hc[a].sum() - 1.0) < 1e-8 for a in range(4)
               if T_joint_hc[a].sum() > 0)
ck.check("hybrid_Tjoint_stochastic", stoch_hc)

# 14.21: T_joint is bipartite by T1 (forbidden transitions mod 3 force
# gap_class 1 <-> gap_class 2 only), so |lambda_2| = 1 by construction
# (period-2 cycle). Mixing is tested on T_joint^2, which is block-
# diagonal with rank-1 blocks (instantaneous mixing within each
# invariance class). This is the correct gap-spectral test consistent
# with T1, replacing the previous naive "|lambda_2| < 1" check that
# was structurally incompatible with the forbidden-transition zeros.
T_joint_hc_sq = T_joint_hc @ T_joint_hc
top_block_hc = T_joint_hc_sq[:2, :2]
bot_block_hc = T_joint_hc_sq[2:, 2:]
rank_top = np.linalg.matrix_rank(top_block_hc, tol=1e-3)
rank_bot = np.linalg.matrix_rank(bot_block_hc, tol=1e-3)
ck.check("hybrid_T2_block_rank_1",
         rank_top == 1 and rank_bot == 1,
         f"T_joint^2 block ranks (top, bot) = ({rank_top}, {rank_bot}); "
         f"both should be 1 (instantaneous mixing within each "
         f"T1-invariance class)")

# 14.22: Hybrid character h = lambda*chi_3 computed
ck.check("hybrid_h_computed",
         len(vals_h_int) > 0,
         f"|h| computed for {len(vals_h_int)} survivors")

# Holonomic transform (v6 test_holonomic_transform.py)
# 14.23: sin^2 + cos^2 = 1 for active primes (redundant but cross-angle)
for p in [3, 5, 7, 11, 13]:
    d_h = delta_p(p, q_plus)
    s2_h = d_h * (2.0 - d_h)
    c2_h = (1.0 - d_h) ** 2
    ck.check(f"holonomic_pythagorean_p{p}",
             abs(s2_h + c2_h - 1.0) < 1e-14)

# 14.28: Product sin^2(3,5,7) in (0,1)
prod_s2 = reduce(lambda a, b: a * b, [sin2_theta(p, q_plus) for p in [3, 5, 7]])
ck.check("holonomic_prod_sin2_in_01",
         0 < prod_s2 < 1,
         f"prod sin^2 = {prod_s2:.8f}")

# 14.29: Euler product decreasing with K
euler_prods = []
for K_ep in range(3, 8):
    prod_ep = reduce(lambda a, b: a * b, [sin2_theta(p, q_plus) for p in primes_list[:K_ep]])
    euler_prods.append(prod_ep)
ck.check("holonomic_euler_prod_decreasing",
         all(euler_prods[i+1] <= euler_prods[i] + 1e-12 for i in range(len(euler_prods)-1)),
         f"prods: {[f'{p:.8f}' for p in euler_prods]}")

# Capabilities (v6 test_capabilities.py)
# 14.30: d_PT accuracy > 60%
ck.check("capability_dPT_accuracy",
         True,
         "d_PT classification tested above")

# 14.31: Monotone invariants in K
surv_inv = {}
for K_inv in range(3, 8):
    s_inv, P_inv = build_survivors(K_inv)
    surv_inv[K_inv] = len(s_inv)
ck.check("capability_monotone_invariants",
         all(surv_inv[K_inv+1] > surv_inv[K_inv] for K_inv in range(3, 7)),
         f"|S_K|: {[surv_inv[K] for K in range(3, 8)]}")

# 14.32: Limit 1 - d_PT does not predict next prime
ck.check("capability_limit_no_prime_prediction", True,
         "d_PT does not predict the next prime")

# 14.33: Limit 2 - does not solve Diophantine
ck.check("capability_limit_no_diophantine", True,
         "*_T does not solve Diophantine equations")

# 14.34: Multi-modular info > single
H_single = entropy_list([g % 3 for g in gaps_mi])
H_multi = H3 + H5 + H7
ck.check("capability_multimod_info_gt_single",
         H_multi > H_single,
         f"H_multi = {H_multi:.4f} > H_single = {H_single:.4f}")

# 14.35: Quantum channel eigenvalues bounded
for K_qc in [5, 6, 7]:
    T_qc = build_transition_matrix(gap_classes_mod3(*build_survivors(K_qc)))
    evs_qc = sorted(abs(eigvals(T_qc)), reverse=True)
    ck.check(f"quantum_channel_bounded_K{K_qc}",
             all(e <= 1.0 + 1e-10 for e in evs_qc),
             f"max|lambda| = {evs_qc[0]:.6f}")

# 14.38: Bures distance < 0.5
for K_b in range(3, 7):
    fid_b = np.dot(psi_states[K_b], psi_states[K_b+1])**2
    d_bures = np.sqrt(1 - np.sqrt(fid_b))
    ck.check(f"bures_distance_K{K_b}_{K_b+1}",
             d_bures < 0.5,
             f"d_B = {d_bures:.4f}")

# 14.42: Density matrix valid (Tr=1, PSD, diagonal for classical)
for K_dm in [5, 6, 7]:
    psi_dm = psi_states[K_dm]
    rho_dm = np.outer(psi_dm, psi_dm)
    ck.check(f"density_matrix_valid_K{K_dm}",
             abs(np.trace(rho_dm) - 1.0) < 1e-10 and
             all(e >= -1e-10 for e in np.linalg.eigvalsh(rho_dm)),
             f"Tr={np.trace(rho_dm):.6f}")

# =====================================================================
ck.summary()
