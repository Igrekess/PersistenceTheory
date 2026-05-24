"""
PTM_Primev40_A1_deep.py

Trancher PT-A1 : lambda_5(P) -> 2/pi est-ce
  (a) une vraie limite asymptotique,
  (b) une coincidence numerique a N=10^8 specifique,
  (c) un point fixe a identifier autrement (pondere, RMS, ...) ?

3 pistes :

  1. Sieve SEGMENTE pour N=10^9  (memoire ~O(sqrt(N)) + O(segment)) :
     calcul streame de lambda_5(P) et comparaison a 2/pi.
     Si budget wall-clock > ~30 min, on abandonne la piste 1.

  2. Reinterpretation 0.636 ~ 2/pi : scan lambda_5(P)(N) pour
     N in {5e7, 8e7, 1e8, 1.3e8, 2e8}. Si seul N=1e8 matche, c'est
     une resonance. Plus : tenter moyenne ponderee par log N ou sqrt N.

  3. Moyenne / RMS sur canaux p in {5,7,11,13,17,19,23} a N=1e8 :
     mean(lambda_p), mean(|lambda_p|), RMS(lambda_p) ; comparer a 2/pi.

Ne touche pas aux autres v.
"""
from __future__ import annotations
import math
import statistics
import sys
import time
from typing import Dict, Iterable, List, Tuple


# =============================================================
# SIEVE SEGMENTE : streaming primes sur [2, N]
# =============================================================

def _small_primes_upto(M: int) -> List[int]:
    """Eratosthenes simple, uniquement jusqu'a sqrt(N_max)."""
    if M < 2:
        return []
    sieve = bytearray(b"\x01") * (M + 1)
    sieve[0] = 0
    sieve[1] = 0
    i = 2
    while i * i <= M:
        if sieve[i]:
            sieve[i * i : M + 1 : i] = b"\x00" * (((M - i * i) // i) + 1)
        i += 1
    return [i for i in range(2, M + 1) if sieve[i]]


def segmented_sieve(N_max: int, segment_size: int = 10 ** 7) -> Iterable[int]:
    """Generator qui streame les premiers <= N_max.

    Memoire : sqrt(N_max) pour les petits premiers + segment_size octets
    pour le segment courant. A N=1e9 : sqrt=31623 premiers ~300 KB +
    segment 10 MB => < 15 MB total.
    """
    if N_max < 2:
        return
    limit = int(math.isqrt(N_max)) + 1
    small = _small_primes_upto(limit)

    # yield petits premiers directement
    for p in small:
        if p > N_max:
            return
        yield p

    lo = limit + 1
    while lo <= N_max:
        hi = min(lo + segment_size - 1, N_max)
        size = hi - lo + 1
        seg = bytearray(b"\x01") * size
        for p in small:
            if p * p > hi:
                break
            # premier multiple de p >= lo
            start = max(p * p, ((lo + p - 1) // p) * p)
            if start > hi:
                continue
            off = start - lo
            seg[off : size : p] = b"\x00" * (((size - off - 1) // p) + 1)
        for i in range(size):
            if seg[i]:
                yield lo + i
        lo = hi + 1


# =============================================================
# Lambda_5(P) en STREAMING (pas de sieve full en RAM)
# =============================================================

LOG_START = math.log(100.0)
NUM_BINS = 30
_QR5 = frozenset({1, 4})


def lambda5_streaming(primes_iter: Iterable[int], N_max: int,
                      num_bins: int = NUM_BINS) -> Tuple[float, int]:
    """Consomme un iterator de premiers, calcule lambda_5(P) a N_max.

    Retourne (lambda, n_primes_used).
    """
    log_max = math.log(N_max)
    edges = [math.exp(LOG_START + i * (log_max - LOG_START) / num_bins)
             for i in range(num_bins + 1)]
    counts = [[0, 0] for _ in range(num_bins)]
    n_used = 0
    scale = num_bins / (log_max - LOG_START)
    for n in primes_iter:
        if n < 7:
            continue
        if n > N_max:
            break
        r5 = n % 5
        if r5 == 0:
            continue
        # skip 2,3 automatiquement via n >= 7 et n premier
        log_n = math.log(n)
        if log_n < LOG_START:
            continue
        idx = int((log_n - LOG_START) * scale)
        if idx < 0 or idx >= num_bins:
            continue
        counts[idx][0 if r5 in _QR5 else 1] += 1
        n_used += 1
    cum = 0.0
    for i in range(num_bins):
        nq, nn = counts[i]
        phi = math.log2(nn / nq) if (nq > 0 and nn > 0) else 0.0
        d_log = math.log(edges[i + 1]) - math.log(edges[i])
        cum += phi * d_log
    return cum / 0.5, n_used


# =============================================================
# Lambda_p(P) multi-canal streaming
# =============================================================

def qr_set(p: int) -> frozenset:
    return frozenset(pow(k, 2, p) for k in range(1, p))


def lambdas_multi_streaming(primes_iter: Iterable[int], N_max: int,
                             p_channels: List[int],
                             num_bins: int = NUM_BINS) -> Dict[int, float]:
    log_max = math.log(N_max)
    edges = [math.exp(LOG_START + i * (log_max - LOG_START) / num_bins)
             for i in range(num_bins + 1)]
    qrs = {p: qr_set(p) for p in p_channels}
    counts = {p: [[0, 0] for _ in range(num_bins)] for p in p_channels}
    scale = num_bins / (log_max - LOG_START)
    for n in primes_iter:
        if n < 7:
            continue
        if n > N_max:
            break
        log_n = math.log(n)
        if log_n < LOG_START:
            continue
        idx = int((log_n - LOG_START) * scale)
        if idx < 0 or idx >= num_bins:
            continue
        for p in p_channels:
            r = n % p
            if r == 0:
                continue
            counts[p][idx][0 if r in qrs[p] else 1] += 1
    out: Dict[int, float] = {}
    for p in p_channels:
        cum = 0.0
        for i in range(num_bins):
            nq, nn = counts[p][i]
            phi = math.log2(nn / nq) if (nq > 0 and nn > 0) else 0.0
            d_log = math.log(edges[i + 1]) - math.log(edges[i])
            cum += phi * d_log
        out[p] = cum / 0.5
    return out


# =============================================================
# Scan lambda_5 a plusieurs N autour de 1e8 en re-utilisant un
# sieve bytearray unique (plus rapide que re-segmenter).
# =============================================================

def sieve_bytearray(N: int) -> bytearray:
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = 0
    is_p[1] = 0
    i = 2
    while i * i <= N:
        if is_p[i]:
            is_p[i * i : N + 1 : i] = b"\x00" * (((N - i * i) // i) + 1)
        i += 1
    return is_p


def scan_lambda5_multiN(is_p: bytearray, N_list: List[int],
                         num_bins: int = NUM_BINS) -> List[Tuple[int, float]]:
    out = []
    for N in N_list:
        log_max = math.log(N)
        edges = [math.exp(LOG_START + i * (log_max - LOG_START) / num_bins)
                 for i in range(num_bins + 1)]
        counts = [[0, 0] for _ in range(num_bins)]
        scale = num_bins / (log_max - LOG_START)
        n = 7
        while n <= N:
            if is_p[n]:
                r5 = n % 5
                if r5 != 0:
                    log_n = math.log(n)
                    if log_n >= LOG_START:
                        idx = int((log_n - LOG_START) * scale)
                        if 0 <= idx < num_bins:
                            counts[idx][0 if r5 in _QR5 else 1] += 1
            n += 1
        cum = 0.0
        for i in range(num_bins):
            nq, nn = counts[i]
            phi = math.log2(nn / nq) if (nq > 0 and nn > 0) else 0.0
            d_log = math.log(edges[i + 1]) - math.log(edges[i])
            cum += phi * d_log
        out.append((N, cum / 0.5))
    return out


# =============================================================
# Ponderation : moyenne ponderee par log N ou sqrt N
# =============================================================

def weighted_mean(pairs: List[Tuple[int, float]], weight_fn) -> float:
    ws = [weight_fn(N) for N, _ in pairs]
    vs = [lam for _, lam in pairs]
    s = sum(ws)
    return sum(w * v for w, v in zip(ws, vs)) / s if s > 0 else float("nan")


# =============================================================
# fmt
# =============================================================

def fmt_time(s: float) -> str:
    if s < 60:
        return f"{s:.1f}s"
    if s < 3600:
        return f"{s/60:.1f}min"
    return f"{s/3600:.2f}h"


# =============================================================
# Main
# =============================================================

def main():
    TWO_OVER_PI = 2.0 / math.pi
    t_global = time.time()

    # budget total : 30 min pour Piste 1
    BUDGET_SEC = 30 * 60

    print("=" * 96)
    print("PTM_Primev40_A1_deep — trancher PT-A1 (3 pistes)")
    print(f"  Target 2/pi = {TWO_OVER_PI:.8f}")
    print("=" * 96)

    # ----------------------------------------------------
    # PISTE 2 d'abord (rapide) : on construit un sieve 2e8
    # et on scanne N in {5e7, 8e7, 1e8, 1.3e8, 2e8}.
    # Coût : sieve 2e8 (~200 MB RAM, ~30s) + 5 scans.
    # ----------------------------------------------------
    print("\n" + "-" * 96)
    print("[PISTE 2] scan lambda_5(P)(N) pour N autour de 1e8")
    print("-" * 96)

    N_scan_set = sorted({
        5 * 10**7,
        8 * 10**7,
        10**8,
        13 * 10**7,   # 1.3e8
        2 * 10**8,
    })
    N_scan_max = N_scan_set[-1]

    t0 = time.time()
    print(f"  build bytearray sieve up to N={N_scan_max:,} ...", flush=True)
    is_p_scan = sieve_bytearray(N_scan_max)
    n_primes_scan = sum(is_p_scan)
    print(f"  done : {n_primes_scan:,} primes   ({fmt_time(time.time()-t0)})")

    t1 = time.time()
    scan_pairs = scan_lambda5_multiN(is_p_scan, N_scan_set)
    print(f"  scan done in {fmt_time(time.time()-t1)}")

    print(f"\n  {'N':>13} | {'lambda_5(P)':>12} | {'2/pi':>8} | {'dev':>10} | {'dev%':>7}")
    for N, lam in scan_pairs:
        dev = lam - TWO_OVER_PI
        dev_pct = 100.0 * dev / TWO_OVER_PI
        print(f"  {N:>13,} | {lam:>12.6f} | {TWO_OVER_PI:>8.4f} | "
              f"{dev:>+10.6f} | {dev_pct:>+6.2f}%")

    # resonance check : seul N=1e8 matche ?
    dev_abs = [(N, abs(lam - TWO_OVER_PI)) for N, lam in scan_pairs]
    N_best, best = min(dev_abs, key=lambda x: x[1])
    second_best = sorted(dev_abs, key=lambda x: x[1])[1][1]
    print(f"\n  meilleur match : N={N_best:,}  |dev|={best:.6f}")
    print(f"  2eme meilleur  : |dev|={second_best:.6f}")
    if N_best == 10**8 and second_best > 3 * best:
        print("  -> RESONANCE possible a N=1e8 (2eme meilleur >3x pire)")
    elif best < 0.01 and second_best < 0.03:
        print("  -> plateau ~2/pi sur TOUTE la fenetre, suggere limite asymptotique")
    else:
        print("  -> dispersion moderée ; ni resonance pure, ni plateau net")

    # moyennes ponderees
    w_log = weighted_mean(scan_pairs, lambda N: math.log(N))
    w_sqrt = weighted_mean(scan_pairs, lambda N: math.sqrt(N))
    w_uni = weighted_mean(scan_pairs, lambda N: 1.0)
    print(f"\n  moyenne uniforme    = {w_uni:.6f}  (dev 2/pi = {w_uni-TWO_OVER_PI:+.4f})")
    print(f"  moyenne log-weighted= {w_log:.6f}  (dev 2/pi = {w_log-TWO_OVER_PI:+.4f})")
    print(f"  moyenne sqrt-weight = {w_sqrt:.6f}  (dev 2/pi = {w_sqrt-TWO_OVER_PI:+.4f})")

    # ----------------------------------------------------
    # PISTE 3 : RMS et moyennes sur 7 canaux p a N=1e8
    # ----------------------------------------------------
    print("\n" + "-" * 96)
    print("[PISTE 3] statistiques sur 7 canaux p a N=1e8")
    print("-" * 96)
    t2 = time.time()
    # iter via is_p_scan restreint a <= 1e8
    def primes_le_1e8():
        for n in range(7, 10**8 + 1):
            if is_p_scan[n]:
                yield n

    p_channels = [5, 7, 11, 13, 17, 19, 23]
    lams_multi = lambdas_multi_streaming(primes_le_1e8(), 10**8, p_channels)
    print(f"  done in {fmt_time(time.time()-t2)}")

    print(f"\n  {'p':>4} | {'lambda_p(P)':>12} | {'|lambda|':>10} | {'lambda^2':>10}")
    vals = []
    for p in p_channels:
        v = lams_multi[p]
        vals.append(v)
        print(f"  {p:>4} | {v:>+12.6f} | {abs(v):>10.6f} | {v*v:>10.6f}")

    mean_raw = statistics.mean(vals)
    mean_abs = statistics.mean(abs(v) for v in vals)
    rms = math.sqrt(statistics.mean(v * v for v in vals))
    stdev_raw = statistics.stdev(vals)
    print(f"\n  mean(lambda_p)       = {mean_raw:+.6f}")
    print(f"  mean(|lambda_p|)     = {mean_abs:.6f}")
    print(f"  RMS(lambda_p)        = {rms:.6f}")
    print(f"  stdev(lambda_p)      = {stdev_raw:.6f}")
    print(f"  target 2/pi          = {TWO_OVER_PI:.6f}")
    print(f"\n  |RMS   - 2/pi| = {abs(rms - TWO_OVER_PI):.6f}  ({100*abs(rms-TWO_OVER_PI)/TWO_OVER_PI:.2f}%)")
    print(f"  |mean|λ| - 2/pi| = {abs(mean_abs - TWO_OVER_PI):.6f}  ({100*abs(mean_abs-TWO_OVER_PI)/TWO_OVER_PI:.2f}%)")

    if abs(rms - TWO_OVER_PI) < 0.05:
        print("  -> RMS ~ 2/pi : PT-A1 serait 'amplitude typique' (RMS), pas canal specifique")
    if abs(mean_abs - TWO_OVER_PI) < 0.05:
        print("  -> mean(|lambda|) ~ 2/pi : interpretation amplitude moyenne")

    # libere le sieve 2e8 avant Piste 1 (gros)
    del is_p_scan

    # ----------------------------------------------------
    # PISTE 1 : lambda_5(P) a N=1e9 via sieve segmente
    # ----------------------------------------------------
    print("\n" + "-" * 96)
    print("[PISTE 1] lambda_5(P) a N=1e9 via sieve segmente")
    print("-" * 96)
    elapsed = time.time() - t_global
    remaining = BUDGET_SEC - elapsed
    print(f"  elapsed so far : {fmt_time(elapsed)}  /  budget {fmt_time(BUDGET_SEC)}")
    if remaining < 60:
        print("  -> budget presque epuise, piste 1 SKIP")
        lam9 = None
    else:
        N9 = 10**9
        print(f"  computing lambda_5(P) streaming for N={N9:,}")
        print(f"  segment_size = 2e7 (expected time 10-30 min)")
        t3 = time.time()
        try:
            lam9, n_used = lambda5_streaming(
                segmented_sieve(N9, segment_size=2 * 10**7),
                N_max=N9,
            )
            dur = time.time() - t3
            print(f"  done in {fmt_time(dur)}   (primes used : {n_used:,})")
            dev9 = lam9 - TWO_OVER_PI
            print(f"\n  lambda_5(P) at N=1e9 = {lam9:.6f}")
            print(f"  2/pi                 = {TWO_OVER_PI:.6f}")
            print(f"  dev                  = {dev9:+.6f}  ({100*dev9/TWO_OVER_PI:+.3f}%)")
            if abs(dev9) < 0.01 * TWO_OVER_PI:
                print("  -> match < 1% : PT-A1 tres plausible comme VRAIE LIMITE")
            elif abs(dev9) > 0.10 * TWO_OVER_PI:
                print("  -> ecart > 10% : PT-A1 est COINCIDENCE 1e8 specifique")
            else:
                print("  -> ecart intermediaire : plateau faible, indecision")
        except KeyboardInterrupt:
            print(f"  INTERROMPU apres {fmt_time(time.time()-t3)}")
            lam9 = None
        except Exception as exc:
            print(f"  ERREUR : {exc}")
            lam9 = None

    # ----------------------------------------------------
    # VERDICT FINAL
    # ----------------------------------------------------
    print("\n" + "=" * 96)
    print("VERDICT FINAL PT-A1")
    print("=" * 96)

    # 3 critères
    lam8 = next(lam for N, lam in scan_pairs if N == 10**8)
    d8 = abs(lam8 - TWO_OVER_PI)
    # plateau check
    devs_scan = [abs(lam - TWO_OVER_PI) for _, lam in scan_pairs]
    plateau = max(devs_scan) < 0.05

    if lam9 is not None:
        d9 = abs(lam9 - TWO_OVER_PI)
        if d9 < 0.01 * TWO_OVER_PI and plateau:
            verdict = "VRAIE LIMITE ASYMPTOTIQUE"
        elif d9 > 0.10 * TWO_OVER_PI:
            verdict = "COINCIDENCE specifique a N=1e8"
        elif d8 < 0.01 and d9 > 0.03:
            verdict = "RESONANCE a N=1e8 (non asymptotique)"
        else:
            verdict = "PLATEAU faible, interpretation alternative necessaire"
    else:
        if plateau:
            verdict = "plateau ~2/pi sur [5e7,2e8], PT-A1 PLAUSIBLE sans N=1e9"
        else:
            verdict = "dispersion dans fenetre fine, PT-A1 probablement COINCIDENCE"

    if abs(rms - TWO_OVER_PI) < 0.03:
        verdict += "  ||  2/pi semble etre l'AMPLITUDE RMS 7-canaux (reinterpretation)"

    print(f"  {verdict}")
    print(f"\n  total time : {fmt_time(time.time()-t_global)}")
    print("=" * 96)


if __name__ == "__main__":
    main()
