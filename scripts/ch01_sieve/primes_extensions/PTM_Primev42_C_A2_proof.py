"""
PTM_Primev42_C_A2_proof.py

CONJECTURE C-A2 (post-v36) :
    Pour les biprimes n = q*m (spf = q), m prime, le ratio
        r_q(N) := λ_obs(biprimes spf=q, m twin-prime) / λ_pred(q, m twin-prime)
    converge vers 8/π² ≈ 0.810570  lorsque N → ∞.

Où λ_pred suit la formule v33 :
    λ_pred = phi * (log N_max - log 100) / s
    phi    = log2(n_NQR / n_QR)
    n_QR   = #{r ∈ {1,4} | ∃ m ∈ {1,2,4} mod 5 : (q*m) mod 5 = r}
    n_NQR  = #{r ∈ {2,3} | ∃ m ∈ {1,2,4} mod 5 : (q*m) mod 5 = r}
    s      = 0.5   (entropie PT)

Tâches :
  1. Mesurer r_q à N ∈ {10^6, 10^7, 10^8} sur q ∈ {7,11,...,97} (22 strates).
  2. Pearson ρ(log q, r_q) — indépendance q confirmée si |ρ| < 0.3.
  3. Étendre à ω=3 (triprimes n = q1*q2*m, m twin-prime), mesurer r_(q1,q2).
  4. Si r_q → 8/π² stable → C-A2 renforcée.
     Si r_q(ω=3) → (8/π²)² = 0.657 → T-O3 (Möbius universel).

Usage :
    python3 PTM_Primev42_C_A2_proof.py        # full: N up to 10^8
    python3 PTM_Primev42_C_A2_proof.py fast   # N up to 10^7 (stop there)
"""
from __future__ import annotations
import math
import statistics
import sys
import time
from typing import List, Tuple


# ================================================================
# Memory-efficient sieves
# ================================================================

def sieve_bytearray(N: int) -> bytearray:
    """Eratosthenes sieve as bytearray. is_p[i]=1 if prime."""
    is_p = bytearray(b'\x01') * (N + 1)
    is_p[0] = 0
    is_p[1] = 0
    i = 2
    while i * i <= N:
        if is_p[i]:
            is_p[i * i:N + 1:i] = b'\x00' * (((N - i * i) // i) + 1)
        i += 1
    return is_p


def sieve_spf_array(N: int) -> List[int]:
    """Smallest prime factor. O(N log log N) time, O(N) memory (int32-equivalent)."""
    spf = [0] * (N + 1)
    for i in range(2, N + 1):
        if spf[i] == 0:
            for j in range(i, N + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    return spf


# ================================================================
# Biprime generator n = q*m, m prime, m twin (m+2 prime)
# ================================================================

def gen_biprimes_spf_q_twin(is_p: bytearray, q: int, hi: int):
    """
    Yields n = q*m, m prime with m+2 prime, and spf(n)=q.
    Thus m > q (so m is strictly the 'second' factor).
    """
    m_max = hi // q
    # For spf(n)=q we need m > q (so q is smallest) AND m != q
    m_start = q + 1
    for m in range(max(m_start, 7), m_max + 1):
        if is_p[m] and m + 2 <= hi and is_p[m + 2]:
            n = q * m
            yield n


def gen_triprimes_q1q2_twin(is_p: bytearray, q1: int, q2: int, hi: int):
    """
    Yields n = q1*q2*m, m prime with m+2 prime, spf(n)=q1, 2nd factor=q2.
    Requires q1 < q2 < m.
    """
    assert q1 < q2
    bound = hi // (q1 * q2)
    m_start = q2 + 1
    for m in range(max(m_start, 7), bound + 1):
        if is_p[m] and m + 2 <= hi and is_p[m + 2]:
            yield q1 * q2 * m


# ================================================================
# λ computation (same scheme as v33)
# ================================================================

_QR5 = {1, 4}
_NQR5 = {2, 3}


def compute_lambda_mod5(integers, N_max: int, num_bins: int = 30,
                        log_start: float = math.log(100)) -> Tuple[float, int]:
    log_max = math.log(N_max)
    edges = [math.exp(log_start + i * (log_max - log_start) / num_bins)
             for i in range(num_bins + 1)]
    counts = [[0, 0] for _ in range(num_bins)]
    total = 0

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
        counts[bin_idx][0 if r in _QR5 else 1] += 1
        total += 1

    cum = 0.0
    for i in range(num_bins):
        nq, nn = counts[i]
        phi = math.log2(nn / nq) if (nq > 0 and nn > 0) else 0.0
        d_log = math.log(edges[i + 1]) - math.log(edges[i])
        cum += phi * d_log
    return cum / 0.5, total


def predict_lambda_biprime(q: int, N_max: int,
                           log_start: float = math.log(100)) -> float:
    """v33 formula, verbatim :
       phi = log2(n_NQR / n_QR) with n_residues from (q*m) mod 5, m in {1,2,4}
       λ_pred = phi * (log N_max - log 100) / 0.5
    """
    m_residues = {1, 2, 4}
    n_residues = set()
    for mr in m_residues:
        nr = (q * mr) % 5
        if nr != 0:
            n_residues.add(nr)
    n_qr = len(n_residues & _QR5)
    n_nqr = len(n_residues & _NQR5)
    if n_qr == 0 or n_nqr == 0:
        return float('inf') if n_qr == 0 else float('-inf')
    phi = math.log2(n_nqr / n_qr)
    return phi * (math.log(N_max) - log_start) / 0.5


def predict_lambda_triprime(q1: int, q2: int, N_max: int,
                            log_start: float = math.log(100)) -> float:
    """Like biprime but residue is (q1*q2*m) mod 5."""
    m_residues = {1, 2, 4}
    n_residues = set()
    q1q2_mod5 = (q1 * q2) % 5
    for mr in m_residues:
        nr = (q1q2_mod5 * mr) % 5
        if nr != 0:
            n_residues.add(nr)
    n_qr = len(n_residues & _QR5)
    n_nqr = len(n_residues & _NQR5)
    if n_qr == 0 or n_nqr == 0:
        return float('inf') if n_qr == 0 else float('-inf')
    phi = math.log2(n_nqr / n_qr)
    return phi * (math.log(N_max) - log_start) / 0.5


# ================================================================
# Main experiments
# ================================================================

PRIMES_STRATA = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
                 61, 67, 71, 73, 79, 83, 89, 97]  # 22 strata


def scan_r_q_at_N(is_p: bytearray, N: int, strata: List[int]) -> List[Tuple]:
    """Return [(q, lam_obs, lam_pred, r_q, count), ...]."""
    results = []
    for q in strata:
        lam_obs, cnt = compute_lambda_mod5(
            gen_biprimes_spf_q_twin(is_p, q, N), N_max=N)
        lam_pred = predict_lambda_biprime(q, N)
        if not math.isfinite(lam_pred) or abs(lam_pred) < 1e-9:
            r = float('nan')
        else:
            r = lam_obs / lam_pred
        results.append((q, lam_obs, lam_pred, r, cnt))
    return results


def pearson(xs, ys) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx * dy > 0 else 0.0


def print_scan_summary(title: str, scan: List[Tuple], target: float):
    print(f"\n{title}")
    print(f"  {'q':>3} | {'count':>9} | {'lam_obs':>10} | {'lam_pred':>10} | {'r_q':>8}")
    print("  " + "-" * 58)
    rs = []
    for q, lo, lp, r, c in scan:
        if math.isfinite(r):
            rs.append(r)
            print(f"  {q:>3} | {c:>9,} | {lo:>+10.4f} | {lp:>+10.4f} | {r:>8.4f}")
        else:
            print(f"  {q:>3} | {c:>9,} | {lo:>+10.4f} | {lp:>+10.4f} |   nan")
    if len(rs) >= 2:
        m = statistics.mean(rs)
        s = statistics.stdev(rs)
        print(f"  mean r_q = {m:.4f} ± {s:.4f}   (target 8/π² = {target:.6f})")
        print(f"  |mean − 8/π²| = {abs(m-target):.4f}   "
              f"({abs(m-target)/(s/math.sqrt(len(rs))):.2f} σ of mean)")
    return rs


def ratio_scaling(fit_points):
    """Fit r_q(N) = r_inf + A*(log N)^α * N^{-β} if possible.
    Simpler : for each N compute mean(r_q) and look at (mean - 8/π²)."""
    target = 8.0 / math.pi ** 2
    print("\nSCALING of mean r_q vs N :")
    print(f"  {'N':>11} | {'mean_r':>8} | {'std':>8} | {'gap_to_target':>13}")
    print("  " + "-" * 48)
    for N, rs in fit_points:
        if len(rs) < 2:
            continue
        m = statistics.mean(rs)
        s = statistics.stdev(rs)
        print(f"  {N:>11,} | {m:>8.4f} | {s:>8.4f} | {m-target:+13.5f}")


# ================================================================
# ω=3 extension
# ================================================================

def scan_triprime_ratios(is_p: bytearray, N: int,
                         q_pairs: List[Tuple[int, int]]) -> List[Tuple]:
    results = []
    for q1, q2 in q_pairs:
        lam_obs, cnt = compute_lambda_mod5(
            gen_triprimes_q1q2_twin(is_p, q1, q2, N), N_max=N)
        lam_pred = predict_lambda_triprime(q1, q2, N)
        if not math.isfinite(lam_pred) or abs(lam_pred) < 1e-9:
            r = float('nan')
        else:
            r = lam_obs / lam_pred
        results.append((q1, q2, lam_obs, lam_pred, r, cnt))
    return results


def print_triprime_scan(title: str, scan: List[Tuple], target: float):
    print(f"\n{title}  (target (8/π²)² = {target:.6f})")
    print(f"  {'q1':>3}*{'q2':<3} | {'count':>9} | {'lam_obs':>10} | {'lam_pred':>10} | {'r':>8}")
    print("  " + "-" * 62)
    rs = []
    for q1, q2, lo, lp, r, c in scan:
        if math.isfinite(r):
            rs.append(r)
            print(f"  {q1:>3}*{q2:<3} | {c:>9,} | {lo:>+10.4f} | {lp:>+10.4f} | {r:>8.4f}")
        else:
            print(f"  {q1:>3}*{q2:<3} | {c:>9,} | {lo:>+10.4f} | {lp:>+10.4f} |   nan")
    if len(rs) >= 2:
        m = statistics.mean(rs)
        s = statistics.stdev(rs)
        print(f"  mean r_ω3 = {m:.4f} ± {s:.4f}   (gap to (8/π²)² = {m-target:+.4f})")
    return rs


# ================================================================
# Main
# ================================================================

def main():
    target_8pi2 = 8.0 / math.pi ** 2
    target_sq = target_8pi2 ** 2

    fast = (len(sys.argv) > 1 and sys.argv[1] == "fast")
    Ns = [10 ** 6, 10 ** 7]
    if not fast:
        Ns.append(10 ** 8)

    print("=" * 92)
    print("v42 — C-A2 proof attempt : r_q(biprime, twin-m) → 8/π² ?")
    print(f"     target 8/π² = {target_8pi2:.6f}   (ω=3 target {target_sq:.6f})")
    print("=" * 92)

    print("\nλ_pred formula (v33 verbatim) :")
    print("    phi    = log2(n_NQR / n_QR)")
    print("    n_res  = { (q*m) mod 5 : m in {1,2,4} mod 5, residue != 0 }")
    print("    n_QR   = |n_res ∩ {1,4}|,  n_NQR = |n_res ∩ {2,3}|")
    print("    λ_pred = phi * (log N - log 100) / s,   s = 1/2")

    fit_points = []
    last_scan = None
    for N in Ns:
        t0 = time.time()
        print(f"\n--- Sieve N = {N:,} ---", flush=True)
        is_p = sieve_bytearray(N)
        t_sieve = time.time() - t0
        print(f"    sieve built in {t_sieve:.1f}s   primes <= N : {sum(is_p):,}",
              flush=True)

        t1 = time.time()
        scan = scan_r_q_at_N(is_p, N, PRIMES_STRATA)
        t_scan = time.time() - t1
        print(f"    strata scan in {t_scan:.1f}s")
        rs = print_scan_summary(f"Scan r_q at N = {N:,}", scan, target_8pi2)
        fit_points.append((N, rs))
        last_scan = (N, is_p, scan)

        # ω=3 test at the largest N to avoid redundant sieves
        if N == Ns[-1]:
            print("\n" + "=" * 92)
            print("ω=3 extension : triprimes n = q1*q2*m, m twin-prime")
            print("=" * 92)
            # Small q1*q2 so that N/(q1*q2) gives enough twin primes
            q_pairs = [(7, 11), (7, 13), (7, 17), (7, 19), (7, 23), (7, 29),
                       (11, 13), (11, 17), (11, 19), (13, 17), (13, 19),
                       (17, 19), (11, 23), (13, 23), (17, 23)]
            tri_scan = scan_triprime_ratios(is_p, N, q_pairs)
            tri_rs = print_triprime_scan(
                f"ω=3 triprime scan at N = {N:,}", tri_scan, target_sq)

        del is_p

    # Scaling analysis
    ratio_scaling(fit_points)

    # Pearson test on largest-N biprime scan
    if last_scan is not None:
        N, _, scan = last_scan
        pairs = [(math.log(q), r) for q, lo, lp, r, c in scan if math.isfinite(r)]
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        rho = pearson(xs, ys)
        print(f"\nPearson ρ(log q, r_q) at N = {N:,}  :  ρ = {rho:+.4f}")
        if abs(rho) < 0.3:
            print("  |ρ| < 0.3  →  q-independence compatible with constant r_q")
        else:
            print("  |ρ| ≥ 0.3  →  q-dependence present")

    # Verdict
    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    if fit_points:
        final_rs = fit_points[-1][1]
        if len(final_rs) >= 2:
            m = statistics.mean(final_rs)
            s = statistics.stdev(final_rs)
            gap = abs(m - target_8pi2)
            sigma_gap = gap / (s / math.sqrt(len(final_rs)))
            if sigma_gap < 2.0:
                print(f"  C-A2 RENFORCÉE : r_q = {m:.4f} ± {s:.4f} compatible "
                      f"avec 8/π² = {target_8pi2:.4f}  ({sigma_gap:.2f}σ)")
            elif sigma_gap < 3.5:
                print(f"  C-A2 à RAFFINER : r_q = {m:.4f} ± {s:.4f} dévie de "
                      f"{sigma_gap:.2f}σ de 8/π² — chercher terme correctif")
            else:
                print(f"  C-A2 INFIRMÉE : r_q = {m:.4f} s'écarte à "
                      f"{sigma_gap:.2f}σ de 8/π² — constante incorrecte")


if __name__ == "__main__":
    main()
