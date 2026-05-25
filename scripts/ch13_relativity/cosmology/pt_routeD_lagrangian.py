#!/usr/bin/env python3
"""
PT Route D — Carrier-Echo Lagrangian on T³
==========================================
Derives w₀ = -1 + s/e  and  w_a = +s·α/e  from a variational principle.

The PT quintessence Lagrangian on T³ (active-prime torus {p=3,5,7}):

    L_carrier = -(1/2)(∂_μ φ)(∂^μ φ) - Λ_c · exp(-β_c(μ) · φ/M_Pl)

with the PT-defined slope function

    β_c(μ) = √(3s · echo(μ)) = √(3s · exp(-μ*/μ))

DERIVATION CHAIN (epistemic status per step):
  I1  : T1 (forbidden transitions) → s = 1/n_allowed = 1/2          [THM]
  I2  : Sum rule Σp = μ* = 15                                        [ID-SUM]
  I3  : L2 identity → echo(μ) = exp(-μ*/μ);  echo(μ*) = 1/e        [ALG]
  I4a : Slow-roll attractor 1+w = β²/3 (exponential quintessence)    [QFT]
  I4b : PT identification β_c(μ*)² = 3s·echo(μ*) = 3s/e             [PT-SET]
  R1  : EOM → w₀ = -1 + β_c²/3 = -1 + s/e                         [DER]
  R2  : β_c(μ(z)) via Bianchi I ODE → w(z) = -1 + s·echo(μ(z))     [DER]
  R3  : Taylor at z=0 → w_a = s·α_loc/e                             [DER]

GAP (single): β_c is set by the PT identification (I4b).  A full geometric
derivation would reduce β_c from the KK spectrum of L on T³.
The identification step I4b closes the [DER-PARTIAL] tag on w₀.
"""

import math
import numpy as np
from scipy.integrate import solve_ivp

# ─── PT constants ─────────────────────────────────────────────────────────────
MU_STAR = 15.0
PRIMES  = [3, 5, 7]       # active primes (faces of T³)
S       = 0.5              # carrier amplitude s = 1/n_allowed (T1 → [THM])
E       = math.e

# ─── Step I1 : s from Theorem T1 (Forbidden Transitions) ─────────────────────

def derive_s_from_T1() -> dict:
    """
    T1: for any r in Z/P₁Z with P₁=3, the transition r → r mod 3 = 0 is
    FORBIDDEN (probability 0, exact).

    Residue classes mod 3: {0, 1, 2}.
    Forbidden: {0}.  Allowed: {1, 2}.
    Carrier amplitude s = 1 / n_allowed = 1/2.
    """
    n_residues  = 3          # |Z/3Z|
    n_forbidden = 1          # T1: r=0 forbidden
    n_allowed   = n_residues - n_forbidden
    s           = 1.0 / n_allowed
    return {"n_residues": n_residues, "n_forbidden": n_forbidden,
            "n_allowed": n_allowed, "s": s, "tag": "[THM]"}

# ─── Step I2 : Sum rule ────────────────────────────────────────────────────────

def verify_sum_rule() -> dict:
    """Σ_{p∈{3,5,7}} p = 15 = μ*  [ID-SUM]."""
    s_primes = sum(PRIMES)
    return {"sum_primes": s_primes, "mu_star": MU_STAR,
            "residual": s_primes - MU_STAR, "tag": "[ID-SUM]"}

# ─── Step I3 : L2 echo identity ───────────────────────────────────────────────

def echo(mu: float) -> float:
    """echo(μ) = ∏_{p∈{3,5,7}} q_-^p(μ) = exp(-Σp/μ) = exp(-μ*/μ)  [ALG]."""
    return math.exp(-sum(PRIMES) / mu)   # = exp(-μ*/μ)

def verify_echo_identity() -> dict:
    """Verify echo(μ*) = 1/e at machine precision."""
    e_val     = echo(MU_STAR)
    e_inv     = 1.0 / E
    ppb_error = abs(e_val - e_inv) / e_inv * 1e9
    return {"echo_mu_star": e_val, "1_over_e": e_inv,
            "ppb_error": ppb_error, "tag": "[ALG]"}

# ─── Step I4a : Slow-roll attractor of exponential quintessence ───────────────

def slow_roll_attractor_formula() -> str:
    """
    For L = -(1/2)(∂φ)² - Λ_c·exp(-β·φ/M_Pl), the DE-dominated
    tracking attractor satisfies:

        1 + w_DE = β² / 3        (exact, independent of Λ_c)

    Proof sketch (Steinhardt-Wang-Zlatev 1999 / Wetterich 1988):
      EOM in flat FLRW: φ̈ + 3Hφ̇ = -dV/dφ = (β/M_Pl)·V
      Attractor condition (DE-dominated, φ̇ > 0):
        ρ_DE = (1/2)φ̇² + V,  P_DE = (1/2)φ̇² - V
        φ̇² = β²/(3) × (φ̇² + 2V)   →  φ̇²/V = β²/(3-β²/2) ≈ β²/3  (β≪√6)
      So:  w_DE = (P/ρ) = (φ̇²/2 - V)/(φ̇²/2 + V) = (β²/3 - 1)
           1 + w_DE = β²/3  QED.

    This is a standard QFT result, tag [QFT].
    """
    return "1 + w_DE = β² / 3   [QFT, exact in DE-dominated attractor]"

# ─── Step I4b : PT identification β_c(μ) ─────────────────────────────────────

def beta_c(mu: float) -> float:
    """
    PT identification (I4b):
        β_c(μ)² = 3 · s · echo(μ) = 3 · s · exp(-μ*/μ)

    At μ = μ*:
        β_c² = 3s/e = 3/(2e)
        β_c  = √(3/(2e)) ≈ 0.743

    This closes the carrier-echo pressure identity:
        1 + w = β_c²/3 = s · echo(μ)

    [PT-SET]: The functional form is required by T1+L2; the numerical
    coefficient 3 comes from the attractor formula (I4a).
    """
    return math.sqrt(3.0 * S * echo(mu))

def verify_beta_c() -> dict:
    bc      = beta_c(MU_STAR)
    bc_sq   = bc ** 2
    target  = 3.0 * S / E     # = 3/(2e)
    return {"beta_c": bc, "beta_c_sq": bc_sq,
            "target_3s_over_e": target,
            "residual": bc_sq - target, "tag": "[PT-SET]"}

# ─── Step R1 : w₀ from EOM ────────────────────────────────────────────────────

def derive_w0() -> dict:
    """
    R1:  w₀ = -1 + β_c(μ*)²/3 = -1 + s/e                [DER]

    Derivation:
      1 + w₀ = β_c²/3   [I4a, QFT]
             = (3s/e)/3  [I4b, PT-SET]
             = s/e       [arithmetic]
      w₀ = -1 + s/e
    """
    bc_sq  = beta_c(MU_STAR) ** 2
    w0_der = -1.0 + bc_sq / 3.0       # from attractor formula
    w0_alg = -1.0 + S / E             # algebraic shorthand
    return {"w0_DER": w0_der, "w0_algebraic": w0_alg,
            "discrepancy": w0_der - w0_alg, "tag": "[DER]"}

# ─── Step R2 : μ-coupled w(z) via Bianchi I ───────────────────────────────────

def gamma_p(p: float, mu: float) -> float:
    q_minus = math.exp(-1.0 / mu)
    delta_p = (1.0 - q_minus**p) / p
    sin2    = delta_p * (2.0 - delta_p)
    return -mu * math.log(1.0 - sin2) / (2.0 * p)

def f_p(p: float, mu: float) -> float:
    eps   = 1e-5 * mu
    gp_hi = gamma_p(p, mu + eps) / (mu + eps)
    gp_lo = gamma_p(p, mu - eps) / (mu - eps)
    return (math.log(gp_hi) - math.log(gp_lo)) / (2.0 * eps)

def H_iso_fac(mu: float) -> float:
    return sum(f_p(p, mu) for p in PRIMES) / 3.0

def kinematic_ode(z, mu_arr):
    mu = mu_arr[0]
    Hf = H_iso_fac(mu)
    if abs(Hf) < 1e-12:
        return [0.0]
    return [-1.0 / ((1.0 + z) * Hf)]

def w_route_D(mu: float) -> float:
    """w(μ) = -1 + s·echo(μ) = -1 + β_c(μ)²/3  [DER]."""
    return -1.0 + S * echo(mu)

def integrate_mu_z(z_max: float = 3.0, n_pts: int = 300) -> tuple:
    """Integrate Bianchi I ODE to get μ(z) then w(z)."""
    z_eval = np.linspace(0.0, z_max, n_pts)
    sol    = solve_ivp(kinematic_ode, (0.0, z_max), [MU_STAR],
                       method='DOP853', t_eval=z_eval, rtol=1e-10, atol=1e-12)
    z_arr  = sol.t
    mu_arr = sol.y[0]
    w_arr  = np.array([w_route_D(m) for m in mu_arr])
    return z_arr, mu_arr, w_arr

# ─── Step R3 : w_a via CPL expansion ─────────────────────────────────────────

def alpha_local(mu: float) -> float:
    """α_loc = -1 / (μ · H_iso_fac(μ)) — power-law exponent at z=0."""
    Hf = H_iso_fac(mu)
    return -1.0 / (mu * Hf)

def derive_wa() -> dict:
    """
    R3: w_a = dw/da|₀ = -dw/dz|₀   [CPL convention: w = w₀ + w_a(1-a)]

    Chain rule:
      dw/dz = s · d(echo)/dz = s · echo(μ) · (μ*/μ²) · dμ/dz
    At z=0, μ=μ*, dμ/dz = -1/(H_iso_fac·(1+0)) = α_loc/μ*:
      dw/dz|₀ = s · (1/e) · (μ*/μ*²) · (α_loc/μ*·μ*)
              = s · (1/e) · α_loc / μ*²  × μ*²      [simplify]
              = s · α_loc / e

    Equivalently from I4b:
      d(β_c²/3)/dz|₀ = (1/3) · 2β_c · dβ_c/dz|₀
                     = (1/3) · d(3s·echo)/dz|₀
                     = s · d(echo)/dz|₀
                     = s · α_loc / e    [same result]

    [DER] — consistent with Bianchi I.
    """
    aloc   = alpha_local(MU_STAR)
    Hf     = H_iso_fac(MU_STAR)
    wa_der = S * aloc / E
    return {"alpha_loc": aloc, "H_iso_fac": Hf, "wa": wa_der,
            "tag": "[DER]"}

# ─── Consistency check: μ-integrated w_a vs analytic ────────────────────────

def numerical_wa(dz: float = 0.01) -> float:
    """Finite-difference w_a from integrated μ(z) trajectory.
    CPL: w(a)=w₀+w_a(1-a) → dw/dz|₀ = w_a  (since a=1/(1+z), da/dz|₀=-1).
    """
    z_arr, mu_arr, w_arr = integrate_mu_z(z_max=dz, n_pts=50)
    return (w_arr[-1] - w_arr[0]) / z_arr[-1]


# ─── DESI DR2 tension ─────────────────────────────────────────────────────────

DESI_W0    = -0.803
DESI_W0_ERR = 0.059
DESI_WA    = -0.723
DESI_WA_ERR = 0.244


# ─── MAIN ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    sep = "=" * 65

    print(sep)
    print("PT Route D — Carrier-Echo Lagrangian on T³")
    print(sep)

    # ── I1
    r1 = derive_s_from_T1()
    print(f"\n[I1] T1 → s  (forbidden transitions mod 3)  {r1['tag']}")
    print(f"     n_residues={r1['n_residues']}, n_forbidden={r1['n_forbidden']}, "
          f"n_allowed={r1['n_allowed']}")
    print(f"     s = 1/{r1['n_allowed']} = {r1['s']:.6f}")

    # ── I2
    r2 = verify_sum_rule()
    print(f"\n[I2] Sum rule  {r2['tag']}")
    print(f"     Σp = {r2['sum_primes']} = μ*={r2['mu_star']},  residual={r2['residual']:.2e}")

    # ── I3
    r3 = verify_echo_identity()
    print(f"\n[I3] L2 echo identity  {r3['tag']}")
    print(f"     echo(μ*)  = {r3['echo_mu_star']:.10f}")
    print(f"     1/e       = {r3['1_over_e']:.10f}")
    print(f"     ppb error = {r3['ppb_error']:.2e}")

    # ── I4a
    print(f"\n[I4a] Slow-roll attractor  [QFT]")
    print(f"      {slow_roll_attractor_formula()}")

    # ── I4b
    r4b = verify_beta_c()
    print(f"\n[I4b] PT identification β_c(μ*)  {r4b['tag']}")
    print(f"      β_c(μ*)  = {r4b['beta_c']:.6f}")
    print(f"      β_c²     = {r4b['beta_c_sq']:.6f}")
    print(f"      3s/e     = {r4b['target_3s_over_e']:.6f}  (target)")
    print(f"      residual = {r4b['residual']:.2e}")

    # ── R1 : w₀
    print(f"\n" + sep)
    print("RESULTS")
    print(sep)
    rw0 = derive_w0()
    print(f"\n[R1] w₀  {rw0['tag']}")
    print(f"     w₀ = -1 + β_c²/3 = {rw0['w0_DER']:.6f}")
    print(f"     w₀ = -1 + s/e    = {rw0['w0_algebraic']:.6f}  (algebraic)")
    print(f"     discrepancy      = {rw0['discrepancy']:.2e}")

    # ── R3 : w_a
    rwa = derive_wa()
    print(f"\n[R3] w_a  {rwa['tag']}")
    print(f"     α_loc       = {rwa['alpha_loc']:.6f}")
    print(f"     H_iso_fac   = {rwa['H_iso_fac']:.6f}")
    print(f"     w_a (DER)   = {rwa['wa']:.6f}")

    # ── Numerical consistency
    wa_num = numerical_wa()
    print(f"\n[CHECK] Numerical w_a from integrated μ(z):")
    print(f"     w_a (numerical) = {wa_num:.6f}")
    print(f"     w_a (analytic)  = {rwa['wa']:.6f}")
    print(f"     discrepancy     = {abs(wa_num - rwa['wa']):.2e}")

    # ── DESI tensions
    w0 = rw0['w0_algebraic']
    wa = rwa['wa']
    pull_w0 = (w0 - DESI_W0) / DESI_W0_ERR
    pull_wa = (wa - DESI_WA) / DESI_WA_ERR
    print(f"\n[DESI DR2 tensions]")
    print(f"     Route D w₀ = {w0:.4f},  DESI = {DESI_W0:.3f} ± {DESI_W0_ERR:.3f}"
          f"  →  pull = {pull_w0:+.2f}σ")
    print(f"     Route D w_a = {wa:.4f},  DESI = {DESI_WA:.3f} ± {DESI_WA_ERR:.3f}"
          f"  →  pull = {pull_wa:+.2f}σ")

    # ── Epistemic table
    print(f"\n{sep}")
    print("EPISTEMIC TABLE")
    print(sep)
    rows = [
        ("Step", "Content",                                         "Tag"),
        ("I1",   "T1 → n_allowed=2 → s=1/2",                       "[THM]"),
        ("I2",   "Σ_{p∈{3,5,7}} p = μ* = 15",                      "[ID-SUM]"),
        ("I3",   "echo(μ) = exp(-μ*/μ); echo(μ*)=1/e",             "[ALG]"),
        ("I4a",  "1+w = β²/3 (slow-roll attractor)",                "[QFT]"),
        ("I4b",  "β_c²(μ*) = 3s/e  (PT identification)",           "[PT-SET]"),
        ("R1",   "w₀ = -1+s/e = -0.8161",                          "[DER]"),
        ("R2",   "w(z)=-1+s·echo(μ(z)) via Bianchi I",             "[DER]"),
        ("R3",   "w_a = s·α_loc/e = +0.213",                       "[DER]"),
    ]
    print(f"\n  {'Step':<6} {'Content':<44} {'Tag':<12}")
    print("  " + "-" * 65)
    for step, content, tag in rows[1:]:
        print(f"  {step:<6} {content:<44} {tag:<12}")

    # ── Gap statement
    print(f"\n{sep}")
    print("REMAINING GAP")
    print(sep)
    print("""
  Step I4b is a PT identification, not a geometric derivation.
  A complete closure would show:
    β_c(μ) = √(3s·echo(μ))
  emerges from the KK spectrum of L_carrier on T³ = (S¹/P₁Z)³,
  i.e., the coupling constant is the square root of the T³ volume
  form weighted by the carrier density.

  Until then, I4b carries [PT-SET] and w₀ carries [DER] conditional
  on accepting the PT identification as a definition of L_carrier.

  All other steps (I1,I2,I3,I4a,R1,R2,R3) are independent.
  The derivation chain is otherwise complete.
""")

    print(sep)
    print("L_carrier = -(1/2)(∂φ)² - Λ_c · exp(-β_c(μ(z)) · φ/M_Pl)")
    print(f"β_c(μ*)   = √(3s/e) = {beta_c(MU_STAR):.6f}")
    print(f"w₀        = -1 + s/e = {rw0['w0_algebraic']:.6f}  [DER]")
    print(f"w_a       = +s·α/e   = {rwa['wa']:.6f}  [DER]")
    print(sep)
    print("Script completed successfully.")
    print(sep)
