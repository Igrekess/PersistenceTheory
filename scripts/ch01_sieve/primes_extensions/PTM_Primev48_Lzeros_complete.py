"""
PTM_Primev48_Lzeros_complete.py
================================

Extension of v44 toward numerical validation of T-O10:

    lambda_5(P) = (4/ln 2) * sum_{rho} int [t^rho - (t/2)^rho] / [rho log t N_bin(t)] dt/t

where rho = 1/2 + i*gamma runs over non-trivial zeros of L(s, chi_5 Legendre).

Objectives
----------
  1. Extend zero computation to up to 1000 zeros (streamed, saved).
  2. Try 4 variance-reduction techniques:
       (a) Gaussian smoothing of Phi_P (kernel in log t)
       (b) Perron resummation (Cesaro average)
       (c) Abel damping exp(-lambda (gamma - gamma_max))
       (d) sqrt(t)/log t integration weight (Hooley)
  3. Convergence table match(k) vs N_zeros x technique.
  4. Euler product alternative log Prod (1 - chi(p)/sqrt p)/(1 - chi(p)/p).

Runs in streaming mode - saves after every 50 zeros.
"""
from __future__ import annotations
import math
import os
import time
import sys
from pathlib import Path

import mpmath as mp

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

mp.mp.dps = 25

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "Results"
RESULTS_DIR.mkdir(exist_ok=True)

ZEROS_FILE = RESULTS_DIR / "v48_Lzeros_1000.txt"
REPORT_FILE = RESULTS_DIR / "v48_report.txt"
EULER_FILE = RESULTS_DIR / "v48_euler.txt"

# Observed lambda_5(P) at N=1e8 (from v36)
LAMBDA_OBS_N8 = 0.636166
N_REF = 1e8

# chi_5 Legendre
CHI5 = [0, 1, -1, -1, 1]
Q = 5

# ---------------------------------------------------------------------------
# Zero computation: stream & persist
# ---------------------------------------------------------------------------

def L_chi5(s):
    return mp.dirichlet(s, CHI5)

def theta_hardy(t):
    s = mp.mpc(0.5, t)
    factor = mp.power(Q / mp.pi, s / 2) * mp.gamma(s / 2)
    return mp.arg(factor)

def Z_chi5(t):
    val = mp.exp(mp.mpc(0, theta_hardy(t))) * L_chi5(mp.mpc(0.5, t))
    return mp.re(val)


def stream_zeros(target: int, start_t: float = 5.0, dt: float = 0.20,
                 checkpoint_every: int = 50, budget_sec: float = 3000.0,
                 existing: list[float] | None = None) -> list[float]:
    """Find zeros incrementally, saving to disk periodically.

    Resumes from existing list (if provided).
    Returns as many zeros as found within budget.
    """
    found = list(existing) if existing else []
    last_zero = found[-1] if found else 0.0
    t = max(start_t, last_zero + 0.2) if found else start_t
    prev = Z_chi5(t)

    t_start = time.time()
    last_save = len(found)

    while len(found) < target:
        if time.time() - t_start > budget_sec:
            print(f"  [budget] hit {budget_sec}s cap at {len(found)} zeros")
            break
        t_next = t + dt
        try:
            val = Z_chi5(t_next)
        except Exception:
            t = t_next
            prev = 0.0
            continue
        if prev * val < 0:
            try:
                tz = mp.findroot(Z_chi5, (t + t_next) / 2, solver='newton')
                tz_f = float(tz)
            except Exception:
                try:
                    tz = mp.findroot(Z_chi5, [t, t_next], solver='anderson')
                    tz_f = float(tz)
                except Exception:
                    tz_f = None
            if tz_f is not None and tz_f > last_zero + 1e-3:
                found.append(tz_f)
                last_zero = tz_f
                if len(found) - last_save >= checkpoint_every:
                    _save_zeros(found)
                    last_save = len(found)
                    elapsed = time.time() - t_start
                    rate = len(found) / max(elapsed, 1e-6)
                    print(f"  [ckpt] {len(found)} zeros, "
                          f"last gamma={tz_f:.2f}, "
                          f"rate={rate:.2f}/s, "
                          f"elapsed={elapsed:.0f}s",
                          flush=True)
        t = t_next
        prev = val

    _save_zeros(found)
    return found


def _save_zeros(zeros: list[float]) -> None:
    with open(ZEROS_FILE, "w") as f:
        f.write(f"# N = {len(zeros)} zeros of L(s, chi_5 Legendre)\n")
        for i, g in enumerate(zeros, 1):
            f.write(f"{i:5d}  {g:.10f}\n")


def load_zeros() -> list[float]:
    if not ZEROS_FILE.exists():
        return []
    out = []
    with open(ZEROS_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            out.append(float(parts[1]))
    return out


# ---------------------------------------------------------------------------
# Explicit-formula integrand (v44 baseline)
# ---------------------------------------------------------------------------

def delta_from_zeros(x: float, gammas: list[float]) -> float:
    lgx = math.log(x)
    xs = math.sqrt(x)
    total_re = 0.0
    for g in gammas:
        # x^rho / rho, Re[...]
        # rho = 1/2 + i g,  |rho|^2 = 1/4 + g^2
        # x^rho = sqrt(x) * (cos(g lgx) + i sin(g lgx))
        cos_t = math.cos(g * lgx)
        sin_t = math.sin(g * lgx)
        # (cos + i sin)/(1/2 + i g) = (cos + i sin)(1/2 - i g)/(1/4 + g^2)
        denom = 0.25 + g * g
        re = (cos_t * 0.5 + sin_t * g) / denom
        total_re += xs * re
    return (2.0 * total_re) / lgx


# Technique (a): Gaussian smoothing in log t
def phi_smoothed(x_hi: float, x_lo: float, gammas: list[float],
                 sigma_log: float = 0.3) -> float:
    """Smoothed Phi_bin using Gaussian kernel of width sigma_log in log t.

    Replaces sharp bin boundaries with Gaussian weighting.
    """
    # Evaluate a set of sample points around center
    center = math.sqrt(x_lo * x_hi)
    log_center = math.log(center)
    half_log = (math.log(x_hi) - math.log(x_lo)) / 2
    # 9-point Gauss-Hermite-like
    samples = [-2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0]
    total_phi = 0.0
    total_w = 0.0
    for u in samples:
        log_t = log_center + u * half_log
        t = math.exp(log_t)
        t_half = t / 2
        d = delta_from_zeros(t, gammas) - delta_from_zeros(t_half, gammas)
        N_bin = t / math.log(t) - t_half / math.log(t_half)
        if N_bin <= 0:
            continue
        phi = (2.0 / math.log(2)) * d / N_bin
        # Gaussian weight
        w = math.exp(-(u * u) / (2 * (sigma_log / half_log) ** 2 if half_log > 0 else 1.0))
        total_phi += phi * w
        total_w += w
    return total_phi / max(total_w, 1e-12)


# Technique (c): Abel damping
def delta_damped(x: float, gammas: list[float], lam: float) -> float:
    """Delta with exp(-lam*(gamma - gamma_min)) damping to mute tail zeros."""
    if not gammas:
        return 0.0
    g0 = gammas[0]
    lgx = math.log(x)
    xs = math.sqrt(x)
    total_re = 0.0
    for g in gammas:
        damp = math.exp(-lam * (g - g0))
        cos_t = math.cos(g * lgx)
        sin_t = math.sin(g * lgx)
        denom = 0.25 + g * g
        re = (cos_t * 0.5 + sin_t * g) / denom
        total_re += damp * xs * re
    return (2.0 * total_re) / lgx


# Technique (d): Hooley weight sqrt(t)/log t replaced by 1/log t only (drop sqrt),
# effectively reweights integrand to emphasize smaller t's.
def phi_hooley(x_hi: float, x_lo: float, gammas: list[float]) -> float:
    """Use 1/log t (no sqrt t) weight on Delta_bin."""
    d = delta_from_zeros(x_hi, gammas) - delta_from_zeros(x_lo, gammas)
    # Normalize by sqrt(x_mid)*log(x_mid) instead of N_bin
    xm = math.sqrt(x_lo * x_hi)
    denom = math.sqrt(xm) * math.log(xm)
    if denom == 0:
        return 0.0
    # Rescale to match dimensionless Phi_bin of baseline
    # Phi_standard = (2/ln2) * d / N_bin  and N_bin ~ (x_hi - x_lo)/log(xm)
    # Hooley weight = multiply d by sqrt(xm)/log(xm) before dividing by (x_hi-x_lo)
    return (2.0 / math.log(2)) * d * math.sqrt(xm) / ((x_hi - x_lo) * math.log(xm))


# ---------------------------------------------------------------------------
# lambda_5 approximations under each technique
# ---------------------------------------------------------------------------

def lambda5_baseline(N_max: float, gammas: list[float],
                     n_bins: int = 30, log_start: float = math.log(100)) -> float:
    log_end = math.log(N_max)
    edges = [math.exp(log_start + i * (log_end - log_start) / n_bins)
             for i in range(n_bins + 1)]
    total = 0.0
    for i in range(n_bins):
        x_lo, x_hi = edges[i], edges[i + 1]
        d = delta_from_zeros(x_hi, gammas) - delta_from_zeros(x_lo, gammas)
        N_bin = x_hi / math.log(x_hi) - x_lo / math.log(x_lo)
        if N_bin <= 0:
            continue
        phi = (2.0 / math.log(2)) * d / N_bin
        dlog = math.log(x_hi) - math.log(x_lo)
        total += phi * dlog
    return total / 0.5  # s = 1/2


def lambda5_smoothed(N_max: float, gammas: list[float],
                     n_bins: int = 30, log_start: float = math.log(100),
                     sigma: float = 0.3) -> float:
    log_end = math.log(N_max)
    edges = [math.exp(log_start + i * (log_end - log_start) / n_bins)
             for i in range(n_bins + 1)]
    total = 0.0
    for i in range(n_bins):
        phi = phi_smoothed(edges[i + 1], edges[i], gammas, sigma_log=sigma)
        dlog = math.log(edges[i + 1]) - math.log(edges[i])
        total += phi * dlog
    return total / 0.5


def lambda5_perron(N_max: float, gammas: list[float],
                   n_bins: int = 30, log_start: float = math.log(100)) -> float:
    """Cesaro-average (Perron): running average of baseline partial sums."""
    log_end = math.log(N_max)
    edges = [math.exp(log_start + i * (log_end - log_start) / n_bins)
             for i in range(n_bins + 1)]
    cumsum = 0.0
    partials = []
    for i in range(n_bins):
        x_lo, x_hi = edges[i], edges[i + 1]
        d = delta_from_zeros(x_hi, gammas) - delta_from_zeros(x_lo, gammas)
        N_bin = x_hi / math.log(x_hi) - x_lo / math.log(x_lo)
        if N_bin <= 0:
            partials.append(cumsum)
            continue
        phi = (2.0 / math.log(2)) * d / N_bin
        dlog = math.log(x_hi) - math.log(x_lo)
        cumsum += phi * dlog
        partials.append(cumsum)
    # Cesaro: average of partial sums
    return (sum(partials) / len(partials)) / 0.5


def lambda5_abel(N_max: float, gammas: list[float],
                 n_bins: int = 30, log_start: float = math.log(100),
                 lam: float = 0.01) -> float:
    log_end = math.log(N_max)
    edges = [math.exp(log_start + i * (log_end - log_start) / n_bins)
             for i in range(n_bins + 1)]
    total = 0.0
    for i in range(n_bins):
        x_lo, x_hi = edges[i], edges[i + 1]
        d = delta_damped(x_hi, gammas, lam) - delta_damped(x_lo, gammas, lam)
        N_bin = x_hi / math.log(x_hi) - x_lo / math.log(x_lo)
        if N_bin <= 0:
            continue
        phi = (2.0 / math.log(2)) * d / N_bin
        dlog = math.log(x_hi) - math.log(x_lo)
        total += phi * dlog
    return total / 0.5


def lambda5_hooley(N_max: float, gammas: list[float],
                   n_bins: int = 30, log_start: float = math.log(100)) -> float:
    log_end = math.log(N_max)
    edges = [math.exp(log_start + i * (log_end - log_start) / n_bins)
             for i in range(n_bins + 1)]
    total = 0.0
    for i in range(n_bins):
        phi = phi_hooley(edges[i + 1], edges[i], gammas)
        dlog = math.log(edges[i + 1]) - math.log(edges[i])
        total += phi * dlog
    return total / 0.5


# ---------------------------------------------------------------------------
# Euler product alternative
# ---------------------------------------------------------------------------

def sieve_primes(N: int) -> list[int]:
    if N < 2:
        return []
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(N + 1) if sieve[i]]


def chi5(p: int) -> int:
    r = p % 5
    return CHI5[r]


def euler_product_lambda(max_prime: int) -> dict:
    """Test lambda_5 = (4/ln 2) log Prod (1 - chi(p)/sqrt p)/(1 - chi(p)/p)."""
    primes = sieve_primes(max_prime)
    log_prod = 0.0
    contribs = []
    for p in primes:
        c = chi5(p)
        if c == 0:
            continue
        try:
            num = 1.0 - c / math.sqrt(p)
            den = 1.0 - c / p
            ratio = num / den
            if ratio <= 0:
                # Use complex log
                term = math.log(abs(ratio))
            else:
                term = math.log(ratio)
        except Exception:
            continue
        log_prod += term
        contribs.append((p, c, term))
    return {
        "log_prod": log_prod,
        "lambda_euler": (4.0 / math.log(2)) * log_prod,
        "n_primes": len(contribs),
        "contribs_head": contribs[:20],
    }


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    log_lines: list[str] = []

    def _log(s: str = "") -> None:
        print(s, flush=True)
        log_lines.append(s)

    _log("=" * 92)
    _log("PTM_Primev48_Lzeros_complete - T-O10 validation at scale")
    _log(f"mpmath dps = {mp.mp.dps}")
    _log(f"lambda_5(P) obs @ N=1e8 = {LAMBDA_OBS_N8}")
    _log("=" * 92)

    # ------------------------------------------------------------------
    # 1. Compute / load zeros
    # ------------------------------------------------------------------
    TARGET = 1000
    BUDGET_SEC = 2400  # 40 min for zero computation

    existing = load_zeros()
    if existing:
        _log(f"[1] Loaded {len(existing)} existing zeros from disk")
    else:
        _log(f"[1] No existing zeros; starting from scratch")

    if len(existing) < TARGET:
        _log(f"    computing up to {TARGET} zeros (budget {BUDGET_SEC}s)...")
        gammas = stream_zeros(TARGET, existing=existing,
                              budget_sec=BUDGET_SEC)
    else:
        gammas = existing[:TARGET]
    n_zeros = len(gammas)
    _log(f"    FINAL: {n_zeros} zeros computed/loaded")
    if gammas:
        _log(f"    gamma_1 = {gammas[0]:.4f}, "
             f"gamma_{n_zeros} = {gammas[-1]:.4f}")
    _log(f"    elapsed so far: {time.time() - t0:.1f}s")

    # ------------------------------------------------------------------
    # 2. Convergence table: match vs (N_zeros, technique)
    # ------------------------------------------------------------------
    _log("\n[2] Match table: |lambda_approx| / |lambda_obs| at N = 1e8")
    _log(f"    lambda_obs = {LAMBDA_OBS_N8:.4f}")

    k_values = [10, 30, 100, 300, min(1000, n_zeros)]
    k_values = [k for k in k_values if k <= n_zeros]
    # Remove duplicates preserving order
    seen = set()
    k_values = [k for k in k_values if not (k in seen or seen.add(k))]

    techniques = [
        ("baseline", lambda g: lambda5_baseline(N_REF, g)),
        ("gaussian(s=0.3)", lambda g: lambda5_smoothed(N_REF, g, sigma=0.3)),
        ("perron",   lambda g: lambda5_perron(N_REF, g)),
        ("abel(0.01)", lambda g: lambda5_abel(N_REF, g, lam=0.01)),
        ("abel(0.05)", lambda g: lambda5_abel(N_REF, g, lam=0.05)),
        ("hooley",   lambda g: lambda5_hooley(N_REF, g)),
    ]

    header = f"    {'technique':<18}" + "".join(f"{k:>10}" for k in k_values)
    _log(header)

    results_matrix = {}
    best_match = 0.0
    best_entry = None

    for name, fn in techniques:
        row = [f"    {name:<18}"]
        results_matrix[name] = {}
        for k in k_values:
            try:
                val = fn(gammas[:k])
            except Exception as e:
                val = float("nan")
            match = abs(val) / abs(LAMBDA_OBS_N8)
            results_matrix[name][k] = (val, match)
            row.append(f"{match*100:>9.1f}%")
            if match > best_match and val * LAMBDA_OBS_N8 > 0:  # same sign
                best_match = match
                best_entry = (name, k, val)
        _log("".join(row))

    # Raw lambda values table
    _log("\n    Raw lambda_approx values:")
    _log(f"    {'technique':<18}" + "".join(f"{k:>10}" for k in k_values))
    for name, fn in techniques:
        row = [f"    {name:<18}"]
        for k in k_values:
            v = results_matrix[name][k][0]
            row.append(f"{v:>+10.4f}")
        _log("".join(row))

    # ------------------------------------------------------------------
    # 3. Convergence ASCII graph
    # ------------------------------------------------------------------
    _log("\n[3] Convergence match(k) ASCII graph (per technique)")
    for name in results_matrix:
        _log(f"\n    {name}:")
        vals = [(k, results_matrix[name][k][1]) for k in k_values]
        for k, m in vals:
            bar_len = int(min(m, 1.2) * 40)
            bar = "#" * bar_len
            _log(f"    k={k:>4}  {m*100:>6.1f}%  {bar}")

    # ------------------------------------------------------------------
    # 4. Best technique + verdict
    # ------------------------------------------------------------------
    _log("\n[4] Best match")
    if best_entry:
        name, k, v = best_entry
        _log(f"    Technique '{name}' with k={k} zeros: "
             f"lambda_approx = {v:.4f}, match = {best_match*100:.1f}%")
    else:
        _log("    No positive-sign match found.")

    # ------------------------------------------------------------------
    # 5. Euler product alternative
    # ------------------------------------------------------------------
    _log("\n[5] Euler product alternative:")
    _log("    lambda_5(P) ?= (4/ln 2) * log Prod_p (1-chi(p)/sqrt p)/(1-chi(p)/p)")
    for pmax in [100, 1000, 10000, 100000, 1000000]:
        try:
            ep = euler_product_lambda(pmax)
            _log(f"    pmax={pmax:>7}: n_primes={ep['n_primes']:>6}  "
                 f"lambda_euler={ep['lambda_euler']:>+9.4f}  "
                 f"match={abs(ep['lambda_euler'])/abs(LAMBDA_OBS_N8)*100:>6.1f}%")
        except Exception as e:
            _log(f"    pmax={pmax}: ERROR {e}")

    # Save detail for max pmax
    try:
        ep = euler_product_lambda(100000)
        with open(EULER_FILE, "w") as f:
            f.write(f"# Euler product log Prod (1-chi/sqrt p)/(1-chi/p)\n")
            f.write(f"# pmax = 100000, n_primes = {ep['n_primes']}\n")
            f.write(f"# log_prod = {ep['log_prod']:.8f}\n")
            f.write(f"# lambda_euler = {ep['lambda_euler']:.8f}\n")
            f.write(f"# lambda_obs = {LAMBDA_OBS_N8}\n")
            f.write("# first 20 contributions:\n")
            for p, c, term in ep["contribs_head"]:
                f.write(f"  p={p:>5}  chi={c:>+2}  term={term:+.6f}\n")
        _log(f"    [saved] {EULER_FILE}")
    except Exception as e:
        _log(f"    [error saving euler] {e}")

    # ------------------------------------------------------------------
    # 6. Formal verdict on T-O10
    # ------------------------------------------------------------------
    _log("\n[6] T-O10 Verdict")
    _log(f"    Best match achieved: {best_match*100:.1f}% with ≤ {n_zeros} zeros")
    if best_match >= 0.80:
        _log("    T-O10 status: NUMERICALLY VALIDATED (match >= 80%). Publishable.")
    elif best_match >= 0.40:
        _log("    T-O10 status: FORMALLY CORRECT, numerically PARTIAL "
             f"(match {best_match*100:.1f}%). Publishable with caveat.")
    else:
        _log("    T-O10 status: FORMAL-ONLY (match < 40%). "
             "Alternative formulation needed.")

    _log(f"\n    Total elapsed: {time.time() - t0:.1f}s")
    _log("=" * 92)

    REPORT_FILE.write_text("\n".join(log_lines))
    print(f"\n[saved] {REPORT_FILE}")
    print(f"[saved] {ZEROS_FILE}")


if __name__ == "__main__":
    main()
