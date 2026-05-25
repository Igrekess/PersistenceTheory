#!/usr/bin/env python3
"""
PT Cosmology Phase 2 - Chantier 5
Dark energy sensitivity test: can beta_super-echo scale-running
lift w_a from O(1e-3) to O(0.1) without breaking any IR observable?

Monograph baseline (ch14 eq:w_eff):
  w_0 = -1 - 2 / (e^gamma (ln N_0)^2 Omega_Lambda)
  w_a = same denominator * Omega_Lambda drift rate
  N_0 = 10^10, gamma = 0.5772 (Euler-Mascheroni)
  Current: w_0 ~= -1.003, w_a ~= -8.3e-4
  DESI DR2: w_0 = -0.75 +/- 0.067, w_a = -0.86 +/- 0.25 (prefers evolving)
  Pull vs PT: 3.78 sigma (w_0), 3.44 sigma (w_a)

This script:
  (i) reproduces the Mertens-drift baseline
  (ii) tests whether super-echo scale-dependent beta can source |w_a|
  (iii) tests whether keeping IR beta_echo=0.1039 (for a_mu) is compatible
       with a UV beta_super-echo boosting w_a
  (iv) produces a sensitivity matrix: w_0(N_0, beta_se) and w_a(N_0, beta_se)
  (v) verdict: can PT accommodate DESI DR2, or does it force falsification?
"""
from __future__ import annotations
import math


# ==== PT constants ====
GAMMA_EULER = 0.5772156649015329   # Euler-Mascheroni
OMEGA_LAMBDA_OBS = 0.685            # Planck 2018 CMB
OMEGA_M_OBS = 0.315
N0_DEFAULT = 1.0e10

# IR echo (UNIVERSAL, ch16 theorem ghost_exhaustion):
BETA_ECHO_IR = 0.1039               # IR, fixed by a_mu at 0.05 sigma

# UV super-echo (scale-dependent, super_ghost_integration_note.md):
# at m_b scale: beta_ghost(m_b; P_max=23) = 0.1598
# near cosmological horizon: TBD, this script explores it
BETA_SUPER_ECHO_UV = 0.1598        # canonical (m_b scale)


# ==== Mertens-drift baseline (ch14 w_eff) ====

def f_ghost(N: float) -> float:
    """F_ghost(N) = 1 - 2 / (e^gamma * ln N), Mertens-type."""
    return 1.0 - 2.0 / (math.exp(GAMMA_EULER) * math.log(N))


def ln_OmegaLambda_drift(N: float) -> float:
    """d ln Omega_Lambda / d ln a under PT baseline.
    ch14: 9.26e-3 at N=10^10."""
    f = f_ghost(N)
    # d F_ghost / d ln N = 2 / (e^gamma * (ln N)^2)
    dfghost = 2.0 / (math.exp(GAMMA_EULER) * math.log(N) ** 2)
    # Omega_Lambda ~ F_ghost * constant (bridge N propto a^3 -> d ln N = 3 d ln a)
    return 3.0 * dfghost / f


def w_0_baseline(N: float = N0_DEFAULT) -> float:
    """Mertens-drift prediction for w_0 at current epoch.
    From ch14: w_0 = -1 - 2 / (e^gamma (ln N)^2 Omega_Lambda)."""
    return -1.0 - 2.0 / (math.exp(GAMMA_EULER) * math.log(N) ** 2 * OMEGA_LAMBDA_OBS)


def w_a_baseline(N: float = N0_DEFAULT) -> float:
    """Mertens-drift prediction for w_a.
    ch14: w_a ~= -8.3e-4 at N=10^10."""
    # Approx: w_a = (d w / d a)|_0 ~= -3 * (w_0+1) * (2*ln(N)/ln^3(N) extra term)
    # Using the triple-log coefficient from ch14:
    # |w_a| / |w_0+1| = O(10^-3) structurally
    w0 = w_0_baseline(N)
    # empirical factor from ch14: w_a ~ 0.28 * (w0 + 1)
    return 0.28 * (w0 + 1.0)


# ==== Super-echo boost scenario ====

def w_with_super_echo(N: float, beta_se: float, mu_eff: float = 14.38) -> tuple[float, float]:
    """
    If beta_super-echo(mu) contributes to the vacuum energy density
    in addition to F_ghost, then:
      Omega_Lambda(a) = F_ghost(N(a)) + delta_Lambda(beta_se, mu(a))
    The super-echo term has a scale-dependent coefficient. Assume:
      delta_Lambda = beta_se * gamma_3 * ln(mu_eff / mu*)   (ansatz)
    where mu_eff(z=0) = 14.38 (ch21 running), mu* = 15.
    Then w_0 and w_a shift.

    Returns (w_0, w_a) under this ansatz.
    """
    # Base quantities
    f = f_ghost(N)
    # Super-echo contribution (drift-like, with 3 super-ghosts):
    # at z=0, mu_eff ~= 14.38 (< mu*=15), so ln(mu_eff/mu*) ~= -0.0424
    logterm = math.log(mu_eff / 15.0)
    # Coefficient: alpha_EM * gamma_3 (dimensional choice to stay O(0.01))
    ALPHA = 1.0 / 137.036
    GAMMA_3 = 0.80761
    delta_Lambda_se = beta_se * ALPHA * GAMMA_3 * logterm
    # Effective Omega_Lambda shift due to super-echo
    # w_0 modification: delta w_0 = -2 * (delta_Lambda_se / Omega_Lambda)
    delta_w0 = -2.0 * delta_Lambda_se / OMEGA_LAMBDA_OBS
    # w_a amplification if super-echo has residual d beta_se / d ln a:
    # Assume d beta_se / d ln a = beta_se * 3 (N propto a^3)
    d_beta_da = 3.0 * beta_se   # d ln beta / d ln a if linear
    delta_wa = -2.0 * beta_se * ALPHA * GAMMA_3 * d_beta_da / OMEGA_LAMBDA_OBS
    w0 = w_0_baseline(N) + delta_w0
    wa = w_a_baseline(N) + delta_wa
    return w0, wa


# ==== DESI DR2 tension function ====

DESI_W0 = -0.75
DESI_W0_ERR = 0.067
DESI_WA = -0.86
DESI_WA_ERR = 0.25


def desi_pull(w0: float, wa: float) -> tuple[float, float]:
    pull_w0 = (w0 - DESI_W0) / DESI_W0_ERR
    pull_wa = (wa - DESI_WA) / DESI_WA_ERR
    return pull_w0, pull_wa


# ==== Main ====

def main() -> None:
    print("=" * 72)
    print("PT COSMOLOGY PHASE 2 -- DARK ENERGY SENSITIVITY (Chantier 5)")
    print("=" * 72)

    N0 = N0_DEFAULT
    w0b = w_0_baseline(N0)
    wab = w_a_baseline(N0)
    print(f"\n Mertens-drift baseline (P20, exploratory):")
    print(f"  N_0            = {N0:.2e}")
    print(f"  F_ghost(N_0)   = {f_ghost(N0)*100:.3f} %")
    print(f"  w_0            = {w0b:.6f}")
    print(f"  w_a            = {wab:.6e}")
    pw0, pwa = desi_pull(w0b, wab)
    print(f"  Pull vs DESI DR2: w_0 {pw0:+.2f}σ, w_a {pwa:+.2f}σ")

    # Scan super-echo boost
    print()
    print("-" * 72)
    print("Sensitivity scan: w(z=0) under super-echo beta_se contribution")
    print("  (ansatz: delta_Lambda = beta_se * alpha_EM * gamma_3 * ln(mu_eff/mu*))")
    print()
    print("  beta_se    | w_0         | w_a         | pull_w0   | pull_wa")
    print("  ---------  | ----------  | ----------  | --------- | ---------")
    for beta_se in [0.0, 0.05, 0.1039, 0.16, 0.30, 0.60, 1.0, 2.0, 5.0, 50.0]:
        w0, wa = w_with_super_echo(N0, beta_se)
        pw0, pwa = desi_pull(w0, wa)
        print(f"  {beta_se:8.4f} | {w0:+.6f}  | {wa:+.6e} | {pw0:+.2f}σ   | {pwa:+.2f}σ")

    # Check: does any beta_se in the admissible range reach DESI DR2?
    print()
    print("-" * 72)
    print("Question: can beta_super-echo be boosted enough to match DESI DR2")
    print("  (target: w_0 ~ -0.75, w_a ~ -0.86) without breaking IR (a_mu)?")
    print()
    # Required delta_w0 to reach DESI
    needed_delta_w0 = DESI_W0 - w0b
    # From formula: delta_w0 = -2 * beta_se * alpha * gamma_3 * ln(14.38/15) / Omega_L
    # => beta_se = delta_w0 * Omega_L / (2 * alpha * gamma_3 * |ln(14.38/15)|)
    needed_beta = abs(needed_delta_w0) * OMEGA_LAMBDA_OBS / (
        2.0 * (1.0/137.036) * 0.80761 * abs(math.log(14.38/15.0))
    )
    print(f"  needed beta_se to reach DESI w_0 = -0.75:  ~{needed_beta:.1f}")
    print(f"  IR universal beta_echo (locked):           {BETA_ECHO_IR}")
    print(f"  UV super-echo @ m_b (from note):           {BETA_SUPER_ECHO_UV}")
    print()
    print(f"  Ratio needed/available_UV: {needed_beta/BETA_SUPER_ECHO_UV:.1f}x")
    print(f"  => Super-echo CANNOT source DESI DR2 magnitude by a factor {needed_beta/BETA_SUPER_ECHO_UV:.1f}x.")
    print(f"  => Consistent with ch21 rem:P20_DESI_DR2 conclusion: 'cannot be patched'.")

    # More refined: what if delta_Lambda is not logarithmic but linear in beta_se?
    print()
    print("-" * 72)
    print("Refined test: alternative ansatz (linear delta_Lambda = c * beta_se)")
    print()
    print("  coeff c required to reach DESI w_0:")
    # delta_w0 = -2 * c * beta_se / Omega_L
    for beta_se in [0.16, 0.60, 1.0]:
        c = abs(needed_delta_w0) * OMEGA_LAMBDA_OBS / (2.0 * beta_se)
        print(f"    beta_se = {beta_se:.3f}: c = {c:.3f}")
    print("  These c-values are order unity but NOT derivable from PT structure")
    print("  (no single ch13 quantity produces c ~ 0.1-1 by construction).")

    # Scale-dependent w(z) prediction
    print()
    print("-" * 72)
    print("PT w(z) prediction (Mertens drift only, ch14 baseline):")
    print()
    print("  z      | N(z)=N_0/(1+z)^3 | F_ghost | Omega_Lambda | w(z)")
    for z in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0]:
        Nz = N0 / (1.0 + z) ** 3
        fz = f_ghost(Nz)
        w_z = -1.0 - 2.0 / (math.exp(GAMMA_EULER) * math.log(Nz) ** 2 * OMEGA_LAMBDA_OBS)
        print(f"  {z:4.2f}   | {Nz:.3e}   | {fz*100:.2f}%   | (est.)       | {w_z:+.6f}")

    # Verdict
    print()
    print("=" * 72)
    print("VERDICT (Chantier 5):")
    print()
    print("  (1) PT Mertens-drift baseline (w_0 ~ -1.003, w_a ~ -8e-4)")
    print("      is 3.78 σ from DESI DR2 (w_0) and 3.44 σ (w_a).")
    print()
    print("  (2) Super-echo scale-dependent beta cannot close the gap:")
    print(f"      needs boost x{needed_beta/BETA_SUPER_ECHO_UV:.0f} above UV value, impossible without a new mechanism.")
    print()
    print("  (3) PT options for DESI DR3:")
    print("      (a) DR3 confirms DR2 at ≥5σ (w_0 ≈ -0.75, w_a ≈ -0.86)")
    print("          => P20 FALSIFIED (but only the exploratory P20, NOT the 28 Tier-1).")
    print("          PT retreats to: 'late-time cosmology outside sieve's IR scope.'")
    print("      (b) DR3 converges to Lambda-CDM (w_0 = -1, w_a = 0) at <1σ")
    print("          => P20 CONFIRMED at zero-parameter level. Rare across cosmologies.")
    print("      (c) DR3 intermediate (w_0 ~ -0.95, w_a ~ -0.2)")
    print("          => P20 marginal; super-echo mechanism research opened.")
    print()
    print("  (4) Conclusion: DO NOT promote P20 to Tier-1. The IR echo ")
    print("      (beta_echo = 0.104, a_mu, alpha_EM) is preserved regardless.")
    print("      The only structurally consistent PT path is option (b) or ")
    print("      acceptance that the sieve's IR derivation applies to the CMB")
    print("      epoch but not to late-time cosmology (scope boundary).")


if __name__ == "__main__":
    main()
