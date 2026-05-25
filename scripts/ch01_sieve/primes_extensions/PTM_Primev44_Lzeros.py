"""
PTM_Primev44_Lzeros.py
======================

Bridge between PT invariant lambda_5(P) and non-trivial zeros of
Dirichlet L-functions L(s, chi) for chi mod 5, via the
Rubinstein-Sarnak explicit formula.

Theoretical setup (cf. report header of monograph App. O)
---------------------------------------------------------
For a primitive Dirichlet character chi mod q, the explicit formula
reads (under GRH):

    psi(x, chi) = sum_{n <= x} Lambda(n) chi(n)
                = - sum_{rho} x^{rho}/rho  +  small boundary terms

with Lambda the von Mangoldt function and rho = 1/2 + i gamma ranging
over the non-trivial zeros of L(s, chi).

Partial summation with weight chi(p) log p on primes (handling prime
powers is O(sqrt x)) gives

    theta(x, chi) = sum_{p <= x} chi(p) log p
                  = - sum_{rho} x^{rho}/rho + O(sqrt x).

For q = 5, Legendre character chi_5 (real, even, primitive), we have
chi_5(p) = +1 if p mod 5 in {1,4} (QR), -1 if in {2,3} (NQR).
Therefore

    #QR(x) - #NQR(x) ~ theta(x, chi_5) / log x
                     ~ -(1 / log x) * sum_{rho} x^{rho}/rho

Let Delta(x) := #NQR(x) - #QR(x) = +(1/log x) sum_{rho} x^{rho}/rho.

Our PT observable on a bin [x/2, x]:
    Phi_P^{(5)}(x) = log_2 ( #NQR(x,x/2) / #QR(x,x/2) )

Write n_QR = pi_QR(x,coprime,5) * (1 + e_+) and
      n_NQR = pi_NQR(x,coprime,5) * (1 + e_-).
pi_QR and pi_NQR are both ~ pi(x,coprime,5)/2 = N(x)/2.
The imbalance is Delta_bin(x) = #NQR[x/2,x] - #QR[x/2,x].

To first order:
    Phi_P^{(5)}(x) = log_2( (N/2 + Delta/2) / (N/2 - Delta/2) )
                   = log_2( (1 + Delta/N) / (1 - Delta/N) )
                   ~ (2 / ln 2) * (Delta_bin(x) / N_bin(x))

and the accumulated invariant

    lambda_5(P) = (1/s) * sum_bins Phi^{(5)} * d log x
                ~ (1/s) * (2/ln 2) * int (Delta_bin(x)/N_bin(x)) d log x
                ~ (1/s) * (2/ln 2) * sum_rho int (x^rho - (x/2)^rho) / (rho * log x * N_bin(x)) d log x

where we used
    Delta_bin(x) = theta(x,chi)/log x - theta(x/2,chi)/log(x/2)
                 ~ -(1/log x) [ sum_rho ( x^rho - (x/2)^rho ) / rho ]  (log smooth).

With s = 1/2 this yields the master formula used below (eq. (*)).

------------------------------------------------------------------------
Tasks implemented
-----------------
    (1) Compute first N_zeros zeros gamma_i of L(s, chi_5_Legendre) via
        mpmath Hardy-like function.
    (2) Evaluate the approximation lambda_5^{approx}(N) from a truncated
        explicit-formula sum, for N in {1e6, 1e7, 1e8}.
    (3) Study truncation effect (vary N_zeros).
    (4) Test whether sum 1/(rho * rho_bar) = sum 1/(1/4 + gamma^2) equals
        2/pi for chi_5 Legendre (candidate identity for the retired
        PT-A1 "lambda_5(P) = 2/pi" resonance).
    (5) Writes a full numerical report to Results/v44_lzeros.txt.

Dependencies: mpmath, numpy. Execution < 5 min for N_zeros = 50.
------------------------------------------------------------------------
"""
from __future__ import annotations
import math
import os
import time
from pathlib import Path

import mpmath as mp

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 25                              # working precision
RESULTS_PATH = Path(__file__).with_name("Results") / "v44_lzeros.txt"
RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)

# Observed values from v36 measurements (see project_pt_primes.md)
LAMBDA_OBS = {
    10**6: 0.671076,
    10**7: 0.615759,
    10**8: 0.636166,
}

# chi_5 Legendre character mod 5 (even, real, primitive)
CHI5 = [0, 1, -1, -1, 1]
Q = 5
A_PARITY = 0       # chi is even -> a = 0 in functional equation

# ---------------------------------------------------------------------------
# Part 1. L-function and zeros
# ---------------------------------------------------------------------------

def L_chi5(s):
    return mp.dirichlet(s, CHI5)

def theta_hardy(t):
    """Hardy theta for chi_5 even: arg( (q/pi)^{s/2} Gamma(s/2) ) at s=1/2+it."""
    s = mp.mpc(0.5, t)
    factor = mp.power(Q / mp.pi, s / 2) * mp.gamma(s / 2)
    return mp.arg(factor)

def Z_chi5(t):
    """Real-valued Z(t) on the critical line, sign changes = zeros."""
    val = mp.exp(mp.mpc(0, theta_hardy(t))) * L_chi5(mp.mpc(0.5, t))
    return mp.re(val)

def find_zeros_chi5(n_zeros: int = 50, t_min: float = 5.0,
                    dt: float = 0.25) -> list[float]:
    """Scan Z(t) for sign changes, refine with Newton/bisection."""
    found: list[float] = []
    t = t_min
    prev = Z_chi5(t)
    last_zero = 0.0
    while len(found) < n_zeros:
        t_next = t + dt
        val = Z_chi5(t_next)
        if prev * val < 0:
            try:
                tz = mp.findroot(Z_chi5, (t + t_next) / 2, solver='newton')
                tz_f = float(tz)
            except Exception:
                tz = mp.findroot(Z_chi5, [t, t_next], solver='anderson')
                tz_f = float(tz)
            if tz_f > last_zero + 1e-3:
                found.append(tz_f)
                last_zero = tz_f
        t = t_next
        prev = val
        # safety cap
        if t > 10_000:
            break
    return found


# ---------------------------------------------------------------------------
# Part 2. Explicit-formula approximation of Delta(x)
# ---------------------------------------------------------------------------

def delta_from_zeros(x: float, gammas: list[float]) -> float:
    """Compute -2 * Re( sum_{i} x^{1/2+i gamma_i}/(1/2+i gamma_i) ) / log x.

    Delta(x) = #NQR(x,coprime 5) - #QR(x,coprime 5)  (cumulative up to x)
    Under GRH with chi_5:
        Delta(x) ~ -(1/log x) * sum_{rho} x^rho/rho * 2  [real part only, complex conjugate pairs]
    """
    lgx = math.log(x)
    xs = math.sqrt(x)
    total = 0.0 + 0.0j
    for g in gammas:
        rho = complex(0.5, g)
        # x^rho = exp(rho * log x) = x^{1/2} * exp(i g log x)
        xrho = xs * complex(math.cos(g * lgx), math.sin(g * lgx))
        total += xrho / rho
    # Sum over rho_bar = complex conjugate -> take 2 * Re
    # Minus sign: Delta = #NQR - #QR but theta(x,chi) = #QR - #NQR weighted by log p
    # So #QR - #NQR ~ -(1/log x) sum x^rho/rho -> Delta(x) = +(1/log x) sum x^rho/rho * 2 Re
    return (2.0 * total.real) / lgx


def delta_bin(x_hi: float, gammas: list[float]) -> float:
    """Delta_bin = Delta(x_hi) - Delta(x_hi/2)  (in the bin [x_hi/2, x_hi])."""
    return delta_from_zeros(x_hi, gammas) - delta_from_zeros(x_hi / 2, gammas)


def phi_bin_theoretical(x_hi: float, gammas: list[float]) -> float:
    """Phi_P^{(5)}(bin) ~ (2/ln 2) * Delta_bin / N_bin.

    N_bin (# primes coprime to 30 in [x_hi/2, x_hi]) approx:
        pi(x_hi) - pi(x_hi/2), restricted to n mod 30 in {1,7,11,13,17,19,23,29}
        Density = phi(30)/30 = 8/30 = 4/15 of all naturals; but since we look
        at *primes*, and primes > 5 are all coprime to 30, we use
        N_bin ~ pi(x_hi) - pi(x_hi/2) ~ x_hi/(2 ln x_hi).

    Here we return the ratio Phi_bin multiplied by d log x per bin.
    """
    d = delta_bin(x_hi, gammas)
    # number of primes in (x_hi/2, x_hi] (PNT estimate)
    N_bin = x_hi / math.log(x_hi) - (x_hi / 2) / math.log(x_hi / 2)
    return (2.0 / math.log(2)) * d / N_bin


def lambda5_approx(N_max: float, gammas: list[float],
                   n_bins: int = 30, log_start: float = math.log(100)) -> float:
    """Approximate lambda_5(P) using explicit formula truncated to gammas.

    lambda_5(P) = (1/s) int_{log_start}^{log N_max} Phi_P^{(5)}(t) d log t
    discretized with n_bins equal log-intervals.
    """
    log_end = math.log(N_max)
    edges = [math.exp(log_start + i * (log_end - log_start) / n_bins)
             for i in range(n_bins + 1)]
    total = 0.0
    s_input = 0.5
    for i in range(n_bins):
        x_lo, x_hi = edges[i], edges[i + 1]
        # phi_bin evaluated with x_hi and x_hi/2 replaced by (x_lo, x_hi)
        # treat the bin directly: Delta_bin = Delta(x_hi) - Delta(x_lo)
        d = delta_from_zeros(x_hi, gammas) - delta_from_zeros(x_lo, gammas)
        # approximate N_bin as pi(x_hi) - pi(x_lo)
        N_bin = x_hi / math.log(x_hi) - x_lo / math.log(x_lo)
        if N_bin <= 0:
            continue
        phi = (2.0 / math.log(2)) * d / N_bin
        dlog = math.log(x_hi) - math.log(x_lo)
        total += phi * dlog
    return total / s_input


# ---------------------------------------------------------------------------
# Part 3. Truncation study
# ---------------------------------------------------------------------------

def truncation_study(all_gammas: list[float], N_max: float,
                     k_values: list[int]) -> list[tuple[int, float]]:
    return [(k, lambda5_approx(N_max, all_gammas[:k])) for k in k_values]


# ---------------------------------------------------------------------------
# Part 4. Zero-sum identity test (2/pi candidate)
# ---------------------------------------------------------------------------

def zero_sum_identities(gammas: list[float]) -> dict:
    """Compute several sums over the first zeros for comparison with 2/pi."""
    s1 = sum(1.0 / (0.25 + g * g) for g in gammas)          # sum 1/(rho rho_bar)
    s2 = sum(1.0 / g for g in gammas)                       # sum 1/gamma (divergent)
    s3 = sum(1.0 / (g * g) for g in gammas)                 # sum 1/gamma^2
    s4 = sum(g / (0.25 + g * g) ** 2 for g in gammas)
    # 2 times s1 because rho and rho_bar pair (zeros come in pairs ±gamma)
    # For chi_5 Legendre (real char): zeros symmetric under t -> -t,
    # so the "full" sum over non-trivial zeros is 2 * s1.
    return {
        "sum_1_over_rho_rhobar": s1,
        "sum_1_over_rho_full (x2)": 2 * s1,
        "sum_1_over_gamma_partial": s2,
        "sum_1_over_gamma2": s3,
        "sum_gamma_over_(1/4+g2)^2": s4,
        "target_2_over_pi": 2.0 / math.pi,
        "target_1_over_pi": 1.0 / math.pi,
        "gap_to_2pi (2*s1 - 2/pi)": 2 * s1 - 2.0 / math.pi,
    }


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    lines: list[str] = []

    def _log(msg: str = "") -> None:
        print(msg, flush=True)
        lines.append(msg)

    _log("=" * 92)
    _log("PTM_Primev44_Lzeros — explicit formula bridge for lambda_5(P)")
    _log(f"mpmath dps = {mp.mp.dps}")
    _log("=" * 92)

    # ------------------------------------------------------------------
    # 1. Find zeros
    # ------------------------------------------------------------------
    N_ZEROS = 60
    _log("\n[1] Computing first %d zeros of L(s, chi_5 Legendre) ..." % N_ZEROS)
    t1 = time.time()
    gammas = find_zeros_chi5(n_zeros=N_ZEROS)
    _log(f"    done in {time.time() - t1:.1f} s")
    _log(f"    first 15 gamma_i: " +
         ", ".join(f"{g:.4f}" for g in gammas[:15]))
    _log(f"    gamma_{N_ZEROS} = {gammas[-1]:.4f}")

    # ------------------------------------------------------------------
    # 2. Lambda_5 approx with all zeros, compare with observed
    # ------------------------------------------------------------------
    _log("\n[2] lambda_5^{approx}(N) from explicit formula vs observed")
    _log(f"    using first {N_ZEROS} zeros")
    _log(f"    {'N':>12}  {'lambda_obs':>12}  {'lambda_approx':>14}  {'diff':>10}")
    for N, obs in LAMBDA_OBS.items():
        lam = lambda5_approx(N, gammas)
        _log(f"    {N:>12,}  {obs:>12.4f}  {lam:>14.4f}  {lam - obs:>+10.4f}")

    # ------------------------------------------------------------------
    # 3. Truncation study at N = 1e7 and N = 1e8
    # ------------------------------------------------------------------
    _log("\n[3] Truncation study: lambda_5^{approx} vs number of zeros included")
    k_list = [1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60]
    for N in (1e7, 1e8):
        _log(f"\n    N = {N:.0e}  (observed = {LAMBDA_OBS[int(N)]:.4f})")
        _log(f"    {'k':>5}  {'lambda_approx':>14}")
        for k, v in truncation_study(gammas, N, k_list):
            _log(f"    {k:>5}  {v:>14.4f}")

    # ------------------------------------------------------------------
    # 4. Zero-sum identity candidates for 2/pi
    # ------------------------------------------------------------------
    _log("\n[4] Candidate identities linking zeros to 2/pi (retired C-A1)")
    ident = zero_sum_identities(gammas)
    for k, v in ident.items():
        _log(f"    {k:<34} = {v:.6f}")

    # Stability: sum_1_over_rho_rhobar as function of k
    _log("\n    convergence of sum_{i=1}^{k} 1/(1/4+gamma_i^2):")
    _log(f"    {'k':>5}  {'sum':>10}  {'2*sum':>10}")
    for k in (1, 5, 10, 20, 30, 40, 50, 60):
        s = sum(1.0 / (0.25 + g * g) for g in gammas[:k])
        _log(f"    {k:>5}  {s:>10.6f}  {2*s:>10.6f}")

    # ------------------------------------------------------------------
    # 5. Verdict on bridge
    # ------------------------------------------------------------------
    _log("\n[5] Verdict")
    lam_N7 = lambda5_approx(1e7, gammas)
    lam_N8 = lambda5_approx(1e8, gammas)
    obs_N7 = LAMBDA_OBS[10**7]
    obs_N8 = LAMBDA_OBS[10**8]
    err7 = abs(lam_N7 - obs_N7) / max(abs(obs_N7), 1e-9)
    err8 = abs(lam_N8 - obs_N8) / max(abs(obs_N8), 1e-9)
    _log(f"    relative error at N=1e7: {err7*100:.2f}%")
    _log(f"    relative error at N=1e8: {err8*100:.2f}%")
    if err7 < 0.05 and err8 < 0.05:
        _log("    VERDICT: explicit-formula truncation MATCHES observed lambda_5(P).")
    else:
        _log("    VERDICT: truncation to %d zeros INSUFFICIENT; main signal" % N_ZEROS)
        _log("             is drowned by bin-scale fluctuation or by the normalization.")
        _log("             (see discussion in v36; bin-edge resonance hypothesis)")

    _log(f"\n    Elapsed: {time.time() - t0:.1f} s")
    _log("=" * 92)

    RESULTS_PATH.write_text("\n".join(lines))
    print(f"\n[saved] {RESULTS_PATH}")


if __name__ == "__main__":
    main()
