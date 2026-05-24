"""
PTM_Primev20_proof.py
Demonstration numerique du theoreme PT :
  1. Bimodalite exacte delta_bar(r) = 9 - 2*(r/5) pour r coprime a 30
  2. Coherence residuelle des premiers : Phi(N) = log2(p_NQR/p_QR)
  3. Integrale cumulative int Phi d log N

Genere les donnees pour le papier LaTeX correspondant.
"""
from __future__ import annotations
import math
from typing import List, Dict, Tuple


def sieve_primes(N: int) -> List[bool]:
    """Sieve of Eratosthenes : retourne is_prime[0..N]."""
    sieve = [True] * (N + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, N + 1, i):
                sieve[j] = False
    return sieve


def legendre_5(n: int) -> int:
    """Symbole de Legendre (n/5). Renvoie +1 (QR), -1 (NQR), 0 (5 | n)."""
    r = n % 5
    if r in (1, 4):
        return +1
    if r in (2, 3):
        return -1
    return 0


# ============================
# Etape 1 : verification exacte du theoreme
# ============================

def verify_theorem_exact() -> Tuple[Dict, bool]:
    """
    Pour chaque r in {1, 7, 11, 13, 17, 19, 23, 29}, enumere toutes les
    classes (a%30, b%30) avec a, b coprimes a 15, a+b ≡ r mod 30, et
    verifie que mean_arc(r) = 9 - 2*(r/5) exactement.
    """
    coprime15 = [r for r in range(30) if math.gcd(r, 15) == 1]
    odd_coprime30 = [r for r in coprime15 if r % 2 == 1]

    print("\n" + "=" * 86)
    print("ETAPE 1 : verification exhaustive du theoreme de bimodalite")
    print("=" * 86)
    print(f"{'r':>3} | {'(r/5)':>6} | {'#classes':>9} | {'arcs':>20} | {'mean_arc':>9} | {'predit':>7}")
    print("-" * 86)

    results = {}
    for r in odd_coprime30:
        valid = []
        for a in coprime15:
            b = (r - a) % 30
            if math.gcd(b, 15) == 1:
                arc = min(abs(a - b), 30 - abs(a - b))
                valid.append((a, b, arc))
        arcs = [v[2] for v in valid]
        mean_arc = sum(arcs) / len(arcs) if arcs else float("nan")
        leg = legendre_5(r)
        predicted = 9 - 2 * leg
        results[r] = (mean_arc, predicted, leg, len(valid), arcs)
        arcs_str = "{" + ",".join(map(str, sorted(arcs))) + "}"
        print(f"{r:>3} | {leg:>+6} | {len(valid):>9} | {arcs_str:>20} | {mean_arc:>9.4f} | {predicted:>7}")

    all_match = all(abs(d[0] - d[1]) < 1e-9 for d in results.values())
    print("-" * 86)
    print(f"Theoreme verifie sur les 8 residus coprimes a 30 : {all_match}")
    return results, all_match


# ============================
# Etape 2 : integrale Phi(N) cumulative
# ============================

def cumulative_phi_integral(N_max: int, num_bins: int, log_start: float = math.log(100)
                            ) -> Tuple[List[Dict], float, float]:
    """
    Stratifie [exp(log_start), N_max] en num_bins bins logarithmiques.
    Pour chaque bin, calcule Phi_prime et Phi_compose, et integre.

    Phi(N) := log2 (#NQR / #QR) sur les n coprimes a 30 dans le bin.
    L'integrale cumulative mesure le depot informationnel coherent
    a travers la croissance de N.
    """
    is_prime_arr = sieve_primes(N_max)
    log_max = math.log(N_max)
    edges = [math.exp(log_start + i * (log_max - log_start) / num_bins) for i in range(num_bins + 1)]

    print("\n" + "=" * 130)
    print(f"ETAPE 2 : integrale Phi d log N stratifiee (N_max = {N_max:,}, {num_bins} bins log)")
    print("=" * 130)
    print(f"{'bin':>22} | {'n_p':>6} {'n_c':>6} | "
          f"{'p_NQR^p':>8} {'p_NQR^c':>8} | "
          f"{'Phi_p':>8} {'Phi_c':>8} {'Phi_diff':>9} | "
          f"{'d_logN':>7} | {'I_p_cum':>9} {'I_diff_cum':>11}")
    print("-" * 130)

    cum_p = 0.0
    cum_d = 0.0
    rows: List[Dict] = []
    for i in range(num_bins):
        lo = int(edges[i])
        hi = int(edges[i + 1])
        n_p_qr = n_p_nqr = n_c_qr = n_c_nqr = 0
        for n in range(max(lo, 7), hi + 1):
            if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
                continue
            leg = legendre_5(n)
            if is_prime_arr[n]:
                if leg == +1:
                    n_p_qr += 1
                else:
                    n_p_nqr += 1
            else:
                if leg == +1:
                    n_c_qr += 1
                else:
                    n_c_nqr += 1
        n_p = n_p_qr + n_p_nqr
        n_c = n_c_qr + n_c_nqr
        p_NQR_p = n_p_nqr / n_p if n_p > 0 else 0.0
        p_NQR_c = n_c_nqr / n_c if n_c > 0 else 0.0
        phi_p = math.log2(n_p_nqr / n_p_qr) if (n_p_qr > 0 and n_p_nqr > 0) else 0.0
        phi_c = math.log2(n_c_nqr / n_c_qr) if (n_c_qr > 0 and n_c_nqr > 0) else 0.0
        phi_diff = phi_p - phi_c
        d_log = math.log(hi) - math.log(lo)
        cum_p += phi_p * d_log
        cum_d += phi_diff * d_log
        rows.append({
            "lo": lo, "hi": hi, "n_p": n_p, "n_c": n_c,
            "p_NQR_p": p_NQR_p, "p_NQR_c": p_NQR_c,
            "phi_p": phi_p, "phi_c": phi_c, "phi_diff": phi_diff,
            "d_log": d_log, "cum_p": cum_p, "cum_d": cum_d,
        })
        label = f"[{lo:>7},{hi:>8}]"
        print(f"{label:>22} | {n_p:>6} {n_c:>6} | "
              f"{p_NQR_p:>8.4f} {p_NQR_c:>8.4f} | "
              f"{phi_p:>+8.4f} {phi_c:>+8.4f} {phi_diff:>+9.4f} | "
              f"{d_log:>7.3f} | {cum_p:>+9.4f} {cum_d:>+11.4f}")

    return rows, cum_p, cum_d


# ============================
# Etape 3 : decroissance asymptotique de Phi
# ============================

def fit_decay_law(rows: List[Dict]) -> Dict:
    """
    Ajuste Phi_diff(N) ~ C * (log N)^a / N^b par regression log-log.
    Predit-on b = 1/2 (Chebyshev) ?
    """
    import statistics
    # On garde les bins ou Phi_diff > 0 (positif systematique attendu)
    pts = [(math.sqrt(r["lo"] * r["hi"]), r["phi_diff"]) for r in rows if r["phi_diff"] > 0]
    if len(pts) < 3:
        return {"slope": float("nan"), "intercept": float("nan"), "n_points": len(pts)}
    log_x = [math.log(x) for x, y in pts]
    log_y = [math.log(y) for x, y in pts]
    n = len(log_x)
    mx = sum(log_x) / n
    my = sum(log_y) / n
    num = sum((x - mx) * (y - my) for x, y in zip(log_x, log_y))
    den = sum((x - mx) ** 2 for x in log_x)
    slope = num / den if den > 0 else float("nan")
    intercept = my - slope * mx
    return {"slope": slope, "intercept": intercept, "n_points": n}


# ============================
# Lancement et impression
# ============================

if __name__ == "__main__":
    print("=" * 86)
    print("DEMONSTRATION PT : bimodalite (theoreme exact) et coherence residuelle")
    print("=" * 86)

    # Etape 1 : verifier exhaustivement le theoreme
    results, all_match = verify_theorem_exact()

    # Etape 2 : integrale cumulative
    rows, cum_p, cum_d = cumulative_phi_integral(N_max=2_000_000, num_bins=20)

    # Etape 3 : pente decroissance
    fit = fit_decay_law(rows)

    print("\n" + "=" * 86)
    print("RESUME")
    print("=" * 86)
    print(f"Theoreme exact (bimodalite QR/NQR mod 5) : verifie = {all_match}")
    print(f"Integrale cumulative I_p     = int Phi_prime  d log N  = {cum_p:+.4f} bits")
    print(f"Integrale cumulative I_diff  = int (Phi_p - Phi_c) d log N = {cum_d:+.4f} bits")
    if not math.isnan(fit["slope"]):
        print(f"Decroissance Phi_diff ~ N^slope  : slope estime = {fit['slope']:+.4f}  "
              f"(prediction Chebyshev = -0.5)")
        print(f"  intercept (log)     = {fit['intercept']:+.4f}    n_points = {fit['n_points']}")

    # Sauvegarde des donnees pour LaTeX
    with open("v20_proof_data.txt", "w") as f:
        f.write(f"# Demonstration PT bimodalite + coherence residuelle\n")
        f.write(f"# Theoreme exact verifie = {all_match}\n")
        f.write(f"# I_p_cum = {cum_p:+.6f} bits\n")
        f.write(f"# I_diff_cum = {cum_d:+.6f} bits\n")
        f.write(f"# slope Phi_diff = {fit['slope']:+.6f} (Chebyshev predit -0.5)\n\n")
        f.write(f"# bin_lo  bin_hi  n_p  n_c  p_NQR_p  p_NQR_c  Phi_p  Phi_c  Phi_diff  cum_p  cum_d\n")
        for r in rows:
            f.write(
                f"{r['lo']:>9} {r['hi']:>9} {r['n_p']:>6} {r['n_c']:>6} "
                f"{r['p_NQR_p']:.6f} {r['p_NQR_c']:.6f} "
                f"{r['phi_p']:+.6f} {r['phi_c']:+.6f} {r['phi_diff']:+.6f} "
                f"{r['cum_p']:+.6f} {r['cum_d']:+.6f}\n"
            )
    print(f"\nDonnees sauvees dans v20_proof_data.txt")
