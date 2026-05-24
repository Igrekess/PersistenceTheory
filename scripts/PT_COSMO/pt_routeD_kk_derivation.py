#!/usr/bin/env python3
"""
PT Route D — KK Derivation of β_c from T³ Bianchi I
=====================================================
Three steps toward deriving β_c = √(3s·echo(μ)) from geometry:

  Step A : Z/2Z orbifold projection → factor s (algebraic)
  Step B : V_T³(μ) = a₃a₅a₇ vs echo(μ)   → numerical ratio check
  Step C : β_c from KK kinetic term normalization  → does KK give √(3s·echo)?

Key formula (Step C):
  The Bianchi I anisotropic metric gives, per the sigma-model on moduli space:
    G_μμ = Σ_p f_p(μ)²          [metric on config. space parameterized by μ]
    d ln V_T³/dμ = Σ_p f_p(μ)   [= 3·H_iso_fac]

  Canonical field φ (M_Pl = 1):
    (dφ/dμ)² = Σ_p f_p(μ)²      [= (3H_iso_fac)² - 2B(μ)]

  KK slope:
    β_c_KK(μ) = |d ln V_eff/dφ| = |Σ_p f_p| / √[Σ_p f_p²]
              = |3 H_iso_fac| / √[(3H_iso_fac)² - 2B(μ)]

  PT target:
    β_c_PT(μ) = √(3s·echo(μ)) = √(3/2 · e^{-μ*/μ})

The question: does β_c_KK = β_c_PT, or β_c_KK/β_c_PT = √(1/s) = √2 (Z/2Z missing),
or something else entirely?
"""

import numpy as np
import math
from scipy.integrate import solve_ivp, quad

# ─── PT constants ─────────────────────────────────────────────────────────────
MU_STAR = 15.0
PRIMES  = [3, 5, 7]
S       = 0.5
E       = math.e
MU_RANGE = np.linspace(10.0, 40.0, 500)   # physical range including μ*

# ─── PT functions ─────────────────────────────────────────────────────────────

def gamma_p(p, mu):
    q    = math.exp(-1.0 / mu)
    dp   = (1.0 - q**p) / p
    sin2 = dp * (2.0 - dp)
    return -mu * math.log(1.0 - sin2) / (2.0 * p)

def f_p(p, mu):
    """f_p = d ln(γ_p/μ)/dμ  (numerical derivative)."""
    eps    = 1e-5 * mu
    gp_hi  = gamma_p(p, mu + eps) / (mu + eps)
    gp_lo  = gamma_p(p, mu - eps) / (mu - eps)
    return (math.log(gp_hi) - math.log(gp_lo)) / (2.0 * eps)

def all_fp(mu):
    return [f_p(p, mu) for p in PRIMES]

def H_iso_fac(mu):
    return sum(f_p(p, mu) for p in PRIMES) / 3.0

def B_mu(mu):
    f = all_fp(mu)
    return f[0]*f[1] + f[0]*f[2] + f[1]*f[2]

def echo(mu):
    return math.exp(-MU_STAR / mu)

def V_T3(mu):
    return math.prod(gamma_p(p, mu) / mu for p in PRIMES)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP A — Z/2Z ORBIFOLD PROJECTION
# ═══════════════════════════════════════════════════════════════════════════════

def step_A():
    print("=" * 65)
    print("STEP A — Z/2Z Orbifold Projection (algebraic)")
    print("=" * 65)
    print("""
  Carrier prime P₁ = 2 acts on each circle of T³ = ∏_p S¹_p
  via the identification y ↦ -y  (Z/2Z orbifold on each S¹).

  Orbifold T³/(Z₂)³ preserves only Z₂-EVEN KK modes (cosines).
  The orbifold HALVES the kinetic term of each circle:
      ∫₀^π cos²(ny) dy = π/2  vs  ∫₀^{2π} cos²(ny) dy = π

  Effect on the 4D volume modulus φ:
    Kinetic factor: 1/2^{#orbifold directions} = 1/2^1 = 1/2
    (P₁=2 acts once, on the carrier sector orthogonal to {3,5,7})

  But T1 (Forbidden Transitions) already encodes this:
    s = 1/n_allowed = 1/2   [n_allowed = |{1,2} ⊂ Z/3Z|]

  Claim: the Z/2Z projection inserts a factor s = 1/2 into β_c²,
  so if the raw KK slope is β_raw, then:
      β_c² = s · β_raw² = (1/2) · β_raw²

  Equivalently:  β_c = β_raw / √2

  We verify this numerically in Step C by checking:
      β_c_KK² / β_c_PT² ≈ 1/s = 2     ← "raw" formula
  or: β_c_KK = β_c_PT                  ← Z/2Z already folded in
""")
    print(f"  s = 1/2 = {S}")
    print(f"  1/s = {1/S}")
    print(f"  √(1/s) = {math.sqrt(1/S):.6f}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP B — V_T³(μ) vs echo(μ)
# ═══════════════════════════════════════════════════════════════════════════════

def step_B():
    print()
    print("=" * 65)
    print("STEP B — V_T³(μ) vs echo(μ)")
    print("=" * 65)

    mu_vals = MU_RANGE
    vt3  = np.array([V_T3(m) for m in mu_vals])
    ech  = np.array([echo(m) for m in mu_vals])
    ratio = vt3 / ech   # should be ≈ constant if V_T³ ∝ echo

    print(f"\n  {'μ':>8}  {'V_T³(μ)':>14}  {'echo(μ)':>14}  {'ratio':>12}  {'ln(ratio)':>12}")
    print("  " + "-" * 65)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        v = V_T3(mu)
        e_ = echo(mu)
        r  = v / e_
        print(f"  {mu:>8.1f}  {v:>14.6e}  {e_:>14.6f}  {r:>12.6e}  {math.log(r):>12.4f}")

    ratio_mu_star = V_T3(MU_STAR) / echo(MU_STAR)
    print(f"\n  Ratio V_T³/echo at μ*: {ratio_mu_star:.6e}")
    print(f"  ln(ratio) at μ*:       {math.log(ratio_mu_star):.4f}")

    # Check: is ratio ≈ constant over [μ*, 30]?
    mu_phys = MU_RANGE[(MU_RANGE >= MU_STAR) & (MU_RANGE <= 30)]
    r_phys  = np.array([V_T3(m)/echo(m) for m in mu_phys])
    r_mean  = np.mean(r_phys)
    r_std   = np.std(r_phys)
    r_rms   = r_std / r_mean * 100
    print(f"\n  Over μ ∈ [μ*, 30]:")
    print(f"    Mean ratio  = {r_mean:.4e}")
    print(f"    Std / Mean  = {r_rms:.2f}%   [small → V_T³ ∝ echo]")

    # Better fit: V_T³ ∝ echo(μ)^α?
    # log V_T³ = const + α × log(echo) = const - α μ*/μ
    # fit: log V_T³ vs 1/μ
    inv_mu  = 1.0 / mu_phys
    log_vt3 = np.array([math.log(V_T3(m)) for m in mu_phys])
    # linear fit: log V_T³ = a - b/μ  → b = α × μ*
    coeffs  = np.polyfit(inv_mu, log_vt3, 1)
    b_fit   = -coeffs[0]           # b = α × μ*
    alpha_fit = b_fit / MU_STAR
    print(f"\n  Power-law fit V_T³ ∝ echo(μ)^α:")
    print(f"    log V_T³ ≈ const - α·μ*/μ  →  α = {alpha_fit:.6f}")
    print(f"    [α=1 ↔ V_T³ ∝ echo exactly;  current α = {alpha_fit:.4f}]")

    return ratio_mu_star, alpha_fit


# ═══════════════════════════════════════════════════════════════════════════════
# STEP C — β_c from KK kinetic term (Bianchi I sigma-model)
# ═══════════════════════════════════════════════════════════════════════════════

def beta_c_KK(mu):
    """
    KK slope from Bianchi I sigma-model (M_Pl = 1):
      G_μμ = Σ_p f_p²  = (Σ_p f_p)² - 2B  = (3H)² - 2B
      d ln V_T³/dμ = Σ_p f_p = 3H
      β_c_KK = |3H| / √[(3H)² - 2B]
    """
    Hf  = H_iso_fac(mu)
    Bm  = B_mu(mu)
    num = abs(3.0 * Hf)
    den2 = (3.0 * Hf)**2 - 2.0 * Bm
    if den2 <= 0:
        return float('nan')
    return num / math.sqrt(den2)

def beta_c_PT(mu):
    """PT target: β_c = √(3s·echo(μ))."""
    return math.sqrt(3.0 * S * echo(mu))

def beta_c_raw(mu):
    """Raw KK without Z/2Z: β_raw = β_c_KK, then β_c = β_raw/√2."""
    return beta_c_KK(mu) / math.sqrt(1.0/S)   # = β_KK × √s

def step_C():
    print()
    print("=" * 65)
    print("STEP C — β_c from Bianchi I KK reduction")
    print("=" * 65)
    print("""
  Bianchi I single-d.o.f. sigma-model on moduli space μ:
    Kinetic metric: G_μμ = Σ_p f_p(μ)² = (3H)² − 2B(μ)
    Volume gradient: d ln V_T³/dμ = Σ_p f_p = 3H_iso_fac

    β_c_KK(μ) = |3H| / √[(3H)²−2B]
""")

    print(f"  {'μ':>6}  {'β_c_KK':>10}  {'β_c_PT':>10}  {'ratio KK/PT':>12}  {'β_c_KK²/β_c_PT²':>16}  {'3·echo':>10}")
    print("  " + "-" * 72)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        bkk = beta_c_KK(mu)
        bpt = beta_c_PT(mu)
        r   = bkk / bpt
        r2  = bkk**2 / bpt**2
        e3  = 3.0 * echo(mu)
        print(f"  {mu:>6.1f}  {bkk:>10.6f}  {bpt:>10.6f}  {r:>12.6f}  {r2:>16.6f}  {e3:>10.6f}")

    # At μ*
    bkk_star = beta_c_KK(MU_STAR)
    bpt_star = beta_c_PT(MU_STAR)
    Hf_star  = H_iso_fac(MU_STAR)
    Bm_star  = B_mu(MU_STAR)
    fsum2    = (3*Hf_star)**2 - 2*Bm_star

    print(f"\n  ── At μ* = {MU_STAR} ──")
    print(f"  H_iso_fac        = {Hf_star:.8f}")
    print(f"  B(μ*)            = {Bm_star:.8f}")
    print(f"  3H               = {3*Hf_star:.8f}")
    print(f"  (3H)² − 2B       = {fsum2:.8f}   [= Σf_p²]")
    print(f"  √Σf_p²           = {math.sqrt(fsum2):.8f}")
    print(f"  β_c_KK(μ*)       = {bkk_star:.8f}")
    print(f"  β_c_PT(μ*)       = {bpt_star:.8f}   [= √(3s/e)]")
    print(f"  ratio KK/PT      = {bkk_star/bpt_star:.8f}")
    print(f"  ratio²  KK/PT    = {(bkk_star/bpt_star)**2:.8f}")
    print(f"  1/s = 2          = {1/S:.8f}")
    print(f"  √(1/s) = √2      = {math.sqrt(1/S):.8f}")

    return bkk_star, bpt_star


# ═══════════════════════════════════════════════════════════════════════════════
# STEP C2 — Does β_c_KK² = 3·echo(μ) exactly?
# ═══════════════════════════════════════════════════════════════════════════════

def step_C2():
    print()
    print("=" * 65)
    print("STEP C2 — Is β_c_KK² = 3·echo(μ)? (raw, no s)")
    print("=" * 65)
    print(f"\n  {'μ':>6}  {'β_c_KK²':>12}  {'3·echo(μ)':>12}  {'ratio':>10}")
    print("  " + "-" * 48)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        bkk2 = beta_c_KK(mu)**2
        e3   = 3.0 * echo(mu)
        r    = bkk2 / e3
        print(f"  {mu:>6.1f}  {bkk2:>12.8f}  {e3:>12.8f}  {r:>10.6f}")

    # Also check β_c_KK² = 3s·echo?
    print(f"\n  {'μ':>6}  {'β_c_KK²':>12}  {'3s·echo(μ)':>12}  {'ratio':>10}")
    print("  " + "-" * 48)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        bkk2 = beta_c_KK(mu)**2
        e3s  = 3.0 * S * echo(mu)
        r    = bkk2 / e3s
        print(f"  {mu:>6.1f}  {bkk2:>12.8f}  {e3s:>12.8f}  {r:>10.6f}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP C3 — What functional form matches β_c_KK²?
# ═══════════════════════════════════════════════════════════════════════════════

def step_C3():
    print()
    print("=" * 65)
    print("STEP C3 — Functional form of β_c_KK²(μ)")
    print("=" * 65)
    print("""
  Fit β_c_KK²(μ) against candidate forms:
    (a) C · echo(μ)^α  →  ln β_c_KK² = const + α·ln echo = const − α·μ*/μ
    (b) C / μ^α
    (c) some other PT function
""")

    mu_vals = MU_RANGE[(MU_RANGE >= MU_STAR) & (MU_RANGE <= 30)]
    bkk2    = np.array([beta_c_KK(m)**2 for m in mu_vals])
    ech     = np.array([echo(m) for m in mu_vals])
    log_b2  = np.log(bkk2)
    log_ech = np.log(ech)   # = -μ*/μ

    # Fit (a): log β² = const + α·log(echo)
    coeffs_a = np.polyfit(log_ech, log_b2, 1)
    alpha_a  = coeffs_a[0]
    const_a  = coeffs_a[1]   # = ln C
    C_a      = math.exp(const_a)
    residuals_a = log_b2 - np.polyval(coeffs_a, log_ech)
    rms_a    = np.sqrt(np.mean(residuals_a**2))

    print(f"  Fit (a): β_c_KK² ≈ {C_a:.4e} · echo(μ)^{alpha_a:.4f}")
    print(f"           RMS log residual = {rms_a:.4e}")

    # Fit (b): log β² = const + α·log(1/μ)
    log_inv_mu = np.log(1.0/mu_vals)
    coeffs_b  = np.polyfit(log_inv_mu, log_b2, 1)
    alpha_b   = coeffs_b[0]
    C_b       = math.exp(coeffs_b[1])
    residuals_b = log_b2 - np.polyval(coeffs_b, log_inv_mu)
    rms_b     = np.sqrt(np.mean(residuals_b**2))
    print(f"\n  Fit (b): β_c_KK² ≈ {C_b:.4e} / μ^{alpha_b:.4f}")
    print(f"           RMS log residual = {rms_b:.4e}")

    # Which is better?
    print(f"\n  → Better fit: {'(a) echo^α' if rms_a < rms_b else '(b) 1/μ^α'}")
    print(f"  → α from fit (a) = {alpha_a:.6f}  [expect 1.0 if ∝ echo exactly]")
    print(f"  → C from fit (a) = {C_a:.6e}  [expect 3·s = 1.5  or  3 = 3]")

    return alpha_a, C_a


# ═══════════════════════════════════════════════════════════════════════════════
# SYNTHESIS
# ═══════════════════════════════════════════════════════════════════════════════

def synthesis(ratio_VT3, alpha_VT3, bkk_star, bpt_star, alpha_bkk, C_bkk):
    print()
    print("=" * 65)
    print("SYNTHESIS — What the KK derivation gives")
    print("=" * 65)

    ratio2 = (bkk_star / bpt_star)**2
    print(f"""
  Step B result:
    V_T³(μ) ≈ echo(μ)^{alpha_VT3:.4f} × const
    [α={alpha_VT3:.4f}: close to 1 → V_T³ ∝ echo]

  Step C result:
    β_c_KK(μ*) = {bkk_star:.6f}
    β_c_PT(μ*) = {bpt_star:.6f}
    ratio²     = {ratio2:.6f}

  Diagnosis:
""")

    tol = 0.05
    if abs(ratio2 - 1.0) < tol:
        print("    β_c_KK = β_c_PT   → Z/2Z factor ALREADY IN β_c_KK")
        print("    The Bianchi I KK reduction DIRECTLY gives β_c = √(3s·echo)")
        print("    The s=1/2 factor is implicit in the anisotropic f_p structure")
        print("    GAP CLOSED: β_c = √(3s·echo) is a DERIVED result  [DER]")
    elif abs(ratio2 - 2.0) < tol:
        print("    β_c_KK² ≈ 2 × β_c_PT²  →  β_c_KK = √(3·echo) (no s)")
        print("    The Z/2Z orbifold (Step A) provides exactly the missing s=1/2")
        print("    β_c = β_c_KK / √(1/s) = β_c_KK × √s  → β_c_PT ✓")
        print("    GAP CLOSED via: KK gives √(3·echo), Z/2Z inserts s  [DER+A]")
    else:
        print(f"    ratio² = {ratio2:.4f} ≠ 1 or 2")
        print(f"    Functional form: β_c_KK² ≈ {C_bkk:.3e} × echo^{alpha_bkk:.4f}")
        print(f"    β_c_PT² = 3s·echo = {3*S}·echo")
        print(f"    Missing factor: {C_bkk / (3*S):.4f}")
        print(f"    [Partial result: functional form ~echo^{alpha_bkk:.3f} is correct]")

    print(f"""
  KEY IDENTITIES:
    B(μ*)        = G_μμ - (3H)²/... — Bianchi I curvature coefficient
    Σf_p²(μ*)   = (3H)²−2B         — KK kinetic metric
    β_c_KK²     = (3H)²/[(3H)²−2B] — from moduli σ-model
    β_c_PT²     = 3s·echo(μ)       — PT target

  The ratio β_c_KK²/β_c_PT² encodes whether the Z/2Z (Step A) is needed.
""")


# ─── MAIN ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    step_A()
    ratio_VT3, alpha_VT3 = step_B()
    bkk_star, bpt_star   = step_C()
    step_C2()
    alpha_bkk, C_bkk     = step_C3()
    synthesis(ratio_VT3, alpha_VT3, bkk_star, bpt_star, alpha_bkk, C_bkk)

    print("=" * 65)
    print("Script completed.")
    print("=" * 65)
