#!/usr/bin/env python3
"""
PT Route D — Instanton Action Derivation
=========================================
Goal: derive S_p^BPS = p/μ for each active prime p ∈ {3,5,7}.
If S_p^BPS = p/μ, then ∏_p exp(−S_p) = echo(μ) = exp(−μ*/μ).

This script tests three routes:

  ROUTE I  — Worldline on Bianchi I torus:
             S_p = m_p × L_p,  L_p = 2π a_p(μ),  need m_p = p/(2πγ_p(μ))
             Check: is m_p a natural PT quantity? Is m_p ∝ p?

  ROUTE II — Thermal interpretation (T = μ):
             echo = ∏_p exp(−p/μ) = P_thermal({3,5,7} all excited at T=μ)
             Carrier field φ = −M_Pl × μ*/μ  (minus thermal entropy)
             Kinetic term: K_thermal = (μ*/μ²)²  vs  K_KK = Σf_p²
             Does K_thermal/K_KK → constant?

  ROUTE III — Worldline on prime-lattice circle of radius 1/μ:
             Action = (prime mass) × (circle circumference) = p × (1/μ)
             Requires: carrier couples to circles of radius 1/μ (dynamical)
             Check: is R=1/μ related to any Bianchi I quantity?

Key numerical checks:
  • m_p(μ) ≡ p/(2πγ_p(μ)): universality ratio m_p/p vs μ
  • K_thermal(μ)/K_KK(μ): constancy over [10,40]
  • If β_c from thermal kinetic term: β_c_thermal = |d ln echo/dφ_thermal|
"""

import math
import numpy as np

# ─── PT constants ─────────────────────────────────────────────────────────────
MU_STAR = 15.0
PRIMES  = [3, 5, 7]
S_      = 0.5
E_      = math.e

# ─── PT functions ─────────────────────────────────────────────────────────────

def gamma_p(p, mu):
    q    = math.exp(-1.0 / mu)
    dp   = (1.0 - q**p) / p
    sin2 = dp * (2.0 - dp)
    return -mu * math.log(1.0 - sin2) / (2.0 * p)

def a_p(p, mu):
    return gamma_p(p, mu) / mu

def f_p(p, mu):
    eps   = 1e-5 * mu
    hi    = gamma_p(p, mu+eps)/(mu+eps)
    lo    = gamma_p(p, mu-eps)/(mu-eps)
    return (math.log(hi) - math.log(lo))/(2*eps)

def H_iso_fac(mu):
    return sum(f_p(p, mu) for p in PRIMES) / 3.0

def K_KK(mu):
    """Bianchi I kinetic metric G_μμ = Σf_p²."""
    return sum(f_p(p, mu)**2 for p in PRIMES)

def echo(mu):
    return math.exp(-MU_STAR / mu)

def K_thermal(mu):
    """Thermal kinetic metric for φ = −M_Pl μ*/μ: (dφ/dμ)² = (μ*/μ²)²."""
    return (MU_STAR / mu**2)**2


# ══════════════════════════════════════════════════════════════════════════════
# ROUTE I — Worldline on Bianchi I torus
# ══════════════════════════════════════════════════════════════════════════════

def route_I():
    print("=" * 65)
    print("ROUTE I — Worldline on Bianchi I torus")
    print("=" * 65)
    print("""
  For each active prime p: instanton wraps S¹_p of circumference L_p = 2π a_p.
  BPS condition: S_p = m_p × L_p = p/μ
  → m_p(μ) = p / (2π γ_p(μ))

  If m_p ∝ p (universality), there is a carrier of UNIVERSAL mass g(μ) = m_p/p.
  Check: does m_p/p = 1/(2πγ_p) vary with p?
""")

    print(f"  {'μ':>5}  {'m₃/3':>10}  {'m₅/5':>10}  {'m₇/7':>10}  {'ratio m₃/m₅':>12}  {'ratio m₃/m₇':>12}")
    print("  " + "-" * 65)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        gp   = [gamma_p(p, mu) for p in PRIMES]
        mp   = [p/(2*math.pi*g) for p, g in zip(PRIMES, gp)]
        mpop = [m/p for m, p in zip(mp, PRIMES)]
        print(f"  {mu:>5.0f}  {mpop[0]:>10.6f}  {mpop[1]:>10.6f}  {mpop[2]:>10.6f}"
              f"  {mp[0]/mp[1]:>12.6f}  {mp[0]/mp[2]:>12.6f}")

    print(f"\n  Expected if m_p ∝ p:  m₃/m₅ = 3/5 = {3/5:.6f},  m₃/m₇ = 3/7 = {3/7:.6f}")

    # Check γ_p values
    print(f"\n  γ_p(μ) values at μ* = {MU_STAR}:")
    for p in PRIMES:
        gp = gamma_p(p, MU_STAR)
        ap = a_p(p, MU_STAR)
        mp = p / (2*math.pi*gp)
        Lp = 2*math.pi*ap
        Sp = mp * Lp
        print(f"    p={p}: γ_p={gp:.6f}, a_p={ap:.4e}, m_p={mp:.6f}, L_p={Lp:.4e}, S_p=m_p×L_p={Sp:.8f}, target={p/MU_STAR:.8f}")

    print(f"\n  VERDICT:")
    gp = [gamma_p(p, MU_STAR) for p in PRIMES]
    mp = [p/(2*math.pi*g) for p, g in zip(PRIMES, gp)]
    mpop = [m/p for m, p in zip(mp, PRIMES)]
    if max(mpop)/min(mpop) - 1 < 0.01:
        print("  m_p/p ≈ CONSTANT → universal carrier mass  [BPS universality holds]")
    else:
        spread = (max(mpop)/min(mpop)-1)*100
        print(f"  m_p/p varies by {spread:.1f}% → NOT universal  [BPS universality broken]")
        print(f"  Ratio γ_p/γ₃: {[gp[i]/gp[0] for i in range(3)]}")
        print(f"  Ratio p×γ₃/(p'×γ_p'): would need γ_p ∝ 1/p for universality")
        # Check if γ_p × p = const
        gpxp = [gp[i]*PRIMES[i] for i in range(3)]
        print(f"  γ_p × p: {gpxp}  [const ↔ γ_p ∝ 1/p]")


# ══════════════════════════════════════════════════════════════════════════════
# ROUTE II — Thermal interpretation: T = μ
# ══════════════════════════════════════════════════════════════════════════════

def route_II():
    print()
    print("=" * 65)
    print("ROUTE II — Thermal interpretation (T = μ, φ = −M_Pl μ*/μ)")
    print("=" * 65)
    print("""
  echo(μ) = ∏_p exp(−p/μ) = ∏_p P_p(T=μ)
  where P_p = Boltzmann weight for "prime p excited" at temperature T=μ.
  The product is the joint probability of {3,5,7} ALL excited (independent).

  Carrier field (entropic parameterization):
    φ(μ) = −M_Pl × μ*/μ   [φ<0, increasing toward μ→∞]
    → μ = −M_Pl μ*/φ  (inverse relation ✓)
    → echo = exp(−μ*/μ) = exp(+φ/M_Pl)   [growing with φ]

  Kinetic term:
    K_thermal(μ) = (dφ/dμ)² = (M_Pl μ*/μ²)²

  Compare with Bianchi I KK kinetic term:
    K_KK(μ) = Σ_p f_p(μ)²
""")

    print(f"  {'μ':>6}  {'K_thermal':>14}  {'K_KK':>12}  {'ratio T/KK':>12}  {'√ratio':>10}  {'√(s·echo)':>12}")
    print("  " + "-" * 68)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        Kt  = K_thermal(mu)
        Kk  = K_KK(mu)
        r   = Kt / Kk
        print(f"  {mu:>6.0f}  {Kt:>14.4e}  {Kk:>12.6f}  {r:>12.6f}  {math.sqrt(r):>10.6f}  {math.sqrt(S_*echo(mu)):>12.6f}")

    # At μ*
    mu = MU_STAR
    Kt = K_thermal(mu)
    Kk = K_KK(mu)
    print(f"\n  At μ* = {mu}:")
    print(f"    K_thermal = (μ*/μ²)² = (1/μ*)² = {Kt:.6e}")
    print(f"    K_KK      = Σf_p²    =            {Kk:.6e}")
    print(f"    ratio     = {Kt/Kk:.6f}")
    print(f"    √ratio    = {math.sqrt(Kt/Kk):.6f}")
    print(f"    3s·echo   = {3*S_*echo(mu):.6f}  [= β_c_PT²]")
    print(f"    ratio / (3s·echo) = {(Kt/Kk)/(3*S_*echo(mu)):.6f}")

    # β_c from thermal field
    print(f"""
  β_c from thermal field:
    V_thermal(φ) = Λ_c × echo(μ(φ)) = Λ_c × exp(+φ/M_Pl)
    β_c_thermal  = M_Pl |d ln V/dφ| = 1   [constant slope]
    → 1+w = β²/3 = 1/3 → w₀ = −2/3  ✗  (wrong)

  The thermal field gives the WRONG w₀ because:
    β_c from thermal = 1 (constant)
    β_c from PT      = √(3s·echo) (μ-dependent)

  Resolution: the thermal field IS the right carrier, but we need
  to RESCALE the kinetic term by K_KK/K_thermal to get the
  Bianchi I normalization. This rescaling multiplies β_c by:
    β_c_KK = β_c_thermal × √(K_thermal/K_KK)
            = 1 × √({Kt/Kk:.4f}) = {math.sqrt(Kt/Kk):.6f}

  But β_c_PT = {math.sqrt(3*S_*echo(mu)):.6f}  ← not equal!
""")


# ══════════════════════════════════════════════════════════════════════════════
# ROUTE III — Prime-lattice worldline (dynamical radius R = 1/μ)
# ══════════════════════════════════════════════════════════════════════════════

def route_III():
    print("=" * 65)
    print("ROUTE III — Prime-lattice worldline (R = 1/μ)")
    print("=" * 65)
    print("""
  Hypothesis: each active prime p creates a cycle of circumference L = 1/μ
  in the "prime-lattice" internal space (dual to the Bianchi I scale).

  Worldline of "mass" m_carrier = p on this cycle:
    S_p^BPS = m_p × L = p × (1/μ) = p/μ  ✓

  Product: ∏_p exp(−S_p) = ∏_p exp(−p/μ) = echo(μ)  ✓

  Question: where does R = 1/μ come from in PT / Bianchi I?
""")

    print(f"  {'μ':>6}  {'1/μ':>10}  {'√(Σa_p²)':>12}  {'1/(3H_iso)':>14}  {'1/√K_KK':>12}")
    print("  " + "-" * 58)
    for mu in [10, 12, 15, 18, 20, 25, 30, 40]:
        invmu  = 1.0/mu
        rms_ap = math.sqrt(sum(a_p(p, mu)**2 for p in PRIMES))
        H      = H_iso_fac(mu)
        inv3H  = abs(1/(3*H))
        sqrtK  = math.sqrt(K_KK(mu))
        print(f"  {mu:>6.0f}  {invmu:>10.6f}  {rms_ap:>12.4e}  {inv3H:>14.6f}  {sqrtK:>12.6f}")

    # Key check: is 1/μ related to Bianchi I?
    print(f"""
  Observation:
    1/μ* = {1/MU_STAR:.6f}
    1/|3H_iso(μ*)| = {1/abs(3*H_iso_fac(MU_STAR)):.6f}  [= 1/(3H) — the ODE denominator]
    1/√K_KK(μ*)   = {1/math.sqrt(K_KK(MU_STAR)):.6f}   [= 1/√Σf_p²]

  Are these close to 1/μ?
    1/(3H) / (1/μ) = {MU_STAR/abs(3*H_iso_fac(MU_STAR)):.4f}  ≈ μ/|3H| = α_loc = {-1/H_iso_fac(MU_STAR)/MU_STAR:.4f}⁻¹ × μ
    Recall α_loc = 1/(μ|H_iso|) ≈ 1.15 at μ*

  The Bianchi I ODE gives:  dμ/dz = −1/[(1+z)H_iso_fac]
  At z=0:  dμ/dz|₀ = 1/|H_iso| = α_loc × μ*/μ* × μ* = α_loc × μ*

  The "natural" length in the Bianchi I system is:
    L_Bianchi = |H_iso_fac|^{{−1}} = 1/(|H_iso|) ≈ μ/|3H·μ| = α_loc × μ

  For L = 1/μ: need α_loc × μ × (something) = 1/μ
  → (something) = 1/(α_loc × μ²)  [scale-dependent, not geometrically natural]
""")

    # Check the "natural" circumference = |H_iso|^{-1} gives what S_p?
    print(f"  If L = 1/|3H_iso(μ)| (natural Bianchi scale):")
    print(f"  {'μ':>6}  {'L=1/|3H|':>12}  {'m_p=p/μ×μ/L':>14}  {'S=m×L=p/μ ✓?':>14}")
    print("  " + "-" * 50)
    for mu in [10, 15, 20, 30]:
        L = 1/abs(3*H_iso_fac(mu))
        for p in PRIMES:
            m_p = (p/mu) / L
            S_p = m_p * L
            print(f"  {mu:>6.0f}  p={p}  L={L:.4f}  m_p={m_p:.6f}  S_p={S_p:.6f}  [target {p/mu:.6f}]")


# ══════════════════════════════════════════════════════════════════════════════
# KEY DIAGNOSTIC: γ_p × p — does it converge to a constant?
# ══════════════════════════════════════════════════════════════════════════════

def key_diagnostic():
    print()
    print("=" * 65)
    print("KEY DIAGNOSTIC — γ_p(μ) × p : towards BPS universality")
    print("=" * 65)
    print("""
  BPS universality requires m_p ∝ p, i.e., m_p/p = 1/(2πγ_p) = const.
  This holds iff γ_p(μ) = const (same for all p at fixed μ).
  But γ_p differ by definition. However:

  For LARGE μ: γ_p(μ) → 1/p (as computed: sin²(δ_p) ≈ 2/μ, γ_p ≈ 1/p).
  So: m_p/p = 1/(2πγ_p) → 1/(2π × 1/p) = p/(2π)  ← grows with p!

  For SMALL μ: need to compute numerically.

  The quantity γ_p × p:
""")
    print(f"  {'μ':>6}  {'γ₃×3':>10}  {'γ₅×5':>10}  {'γ₇×7':>10}  {'spread %':>10}")
    print("  " + "-" * 52)
    for mu in [5, 10, 15, 20, 30, 50, 100]:
        gp  = [gamma_p(p, mu) for p in PRIMES]
        gpxp = [g*p for g, p in zip(gp, PRIMES)]
        spread = (max(gpxp)/min(gpxp) - 1)*100
        print(f"  {mu:>6.0f}  {gpxp[0]:>10.6f}  {gpxp[1]:>10.6f}  {gpxp[2]:>10.6f}  {spread:>10.2f}%")

    # What about γ_p × √p?
    print(f"\n  The quantity γ_p × √p:")
    print(f"  {'μ':>6}  {'γ₃×√3':>10}  {'γ₅×√5':>10}  {'γ₇×√7':>10}  {'spread %':>10}")
    print("  " + "-" * 52)
    for mu in [5, 10, 15, 20, 30, 50, 100]:
        gp  = [gamma_p(p, mu) for p in PRIMES]
        gpxsp = [g*math.sqrt(p) for g, p in zip(gp, PRIMES)]
        spread = (max(gpxsp)/min(gpxsp) - 1)*100
        print(f"  {mu:>6.0f}  {gpxsp[0]:>10.6f}  {gpxsp[1]:>10.6f}  {gpxsp[2]:>10.6f}  {spread:>10.2f}%")

    # Find the exponent α such that γ_p × p^α = const
    print(f"\n  Finding α such that γ_p × p^α ≈ const at μ* = {MU_STAR}:")
    gp_star = [gamma_p(p, MU_STAR) for p in PRIMES]
    # γ_p ∝ p^{-α} → ln γ_p = const - α ln p
    log_gp = [math.log(g) for g in gp_star]
    log_pp = [math.log(p) for p in PRIMES]
    # Linear fit: log γ = A - α log p
    alpha = -(log_gp[2]-log_gp[0])/(log_pp[2]-log_pp[0])  # two-point
    print(f"    α (from p=3 and p=7) = {alpha:.4f}")
    alpha3 = -(log_gp[1]-log_gp[0])/(log_pp[1]-log_pp[0])
    alpha5 = -(log_gp[2]-log_gp[1])/(log_pp[2]-log_pp[1])
    print(f"    α (from p=3 and p=5) = {alpha3:.4f}")
    print(f"    α (from p=5 and p=7) = {alpha5:.4f}")
    alpha_mean = (alpha + alpha3 + alpha5) / 3
    print(f"    α mean = {alpha_mean:.4f}")
    print(f"\n    If α ≈ 1 → γ_p ∝ 1/p (BPS would need m_p ∝ p² — not p!)")
    print(f"    If α ≈ 0 → γ_p ≈ const (BPS universality m_p ∝ p holds)")

    gpxpalpha = [g * p**alpha_mean for g, p in zip(gp_star, PRIMES)]
    spread = (max(gpxpalpha)/min(gpxpalpha) - 1)*100
    print(f"    γ_p × p^α_mean: {[f'{v:.4f}' for v in gpxpalpha]}, spread = {spread:.2f}%")


# ══════════════════════════════════════════════════════════════════════════════
# SYNTHESIS
# ══════════════════════════════════════════════════════════════════════════════

def synthesis():
    print()
    print("=" * 65)
    print("SYNTHESIS — State of the instanton derivation")
    print("=" * 65)

    mu = MU_STAR
    gp   = [gamma_p(p, mu) for p in PRIMES]
    mp   = [p/(2*math.pi*g) for p, g in zip(PRIMES, gp)]
    Kt   = K_thermal(mu)
    Kk   = K_KK(mu)

    # α from γ_p ∝ p^{-α}
    log_gp = [math.log(g) for g in gp]
    log_pp = [math.log(p) for p in PRIMES]
    coeffs = np.polyfit(log_pp, log_gp, 1)
    alpha_gp = -coeffs[0]

    print(f"""
  WHAT ROUTE I SHOWS:
    m_p(μ*) = p/(2πγ_p) = {[f'{v:.4f}' for v in mp]}
    γ_p ∝ p^{{−{alpha_gp:.3f}}}  at μ*  [not exactly 1/p]
    → BPS universality (m_p ∝ p) requires γ_p = const, which fails.
    BUT: γ_p ∝ p^{{−{alpha_gp:.3f}}} means m_p ∝ p^{{1+{alpha_gp:.3f}}} = p^{{{1+alpha_gp:.3f}}}

  WHAT ROUTE II SHOWS:
    Thermal field φ = −M_Pl μ*/μ gives V = Λ exp(+φ/M_Pl) → β_c = 1 (wrong).
    K_thermal/K_KK = {Kt/Kk:.4f} (not constant over μ range).
    The thermal identification gives the right ECHO STRUCTURE
    but the wrong kinetic term normalization.

  WHAT ROUTE III SHOWS:
    If the "prime-lattice" circle has radius 1/μ (dynamical),
    then S_p = p × (1/μ) = p/μ  for any universal carrier mass m = p.
    But where does R=1/μ come from geometrically?
    Answer: 1/μ ≈ |H_iso_fac|/|3H| × ... — no clean Bianchi I identity.

  THE CORE DIFFICULTY:
    The Bianchi I torus has CURVED moduli space (γ_p not proportional to 1/p).
    BPS universality requires γ_p ∝ 1/p (or equivalently a_p ∝ p^{{−α}}).
    Numerically: γ_p ∝ p^{{−{alpha_gp:.3f}}} — close to 1/p^1 but not exact.

    IF γ_p were exactly ∝ 1/p, THEN:
      m_p = p/(2πγ_p) ∝ p/(2π × C/p) = p²/(2πC) ∝ p²  (not p)

    This means the BPS condition m_p = p/L_p requires L_p ∝ 1/p² (not 1/μ).

  THE MINIMAL GAP (final form):
    We need to show that the Bianchi I T³ has a DUAL description
    (via T-duality or modular transformation) in which:
      1. The "mass" of the carrier on circle p is m_p = p (= the prime itself)
      2. The "circumference" of circle p in the dual is L̃_p = 1/μ (universal)
      3. The dual action: S̃_p = p × (1/μ) = p/μ  ✓

    Condition 2 (L̃_p = 1/μ) is equivalent to:
      The dual radius R̃_p = 1/(2πμ)  [or R̃ = 1/μ if period = 2π]

    From Bianchi I:  a_p(μ) = γ_p(μ)/μ
    T-dual radius:   R̃_p = c_s²/a_p = c_s² × μ/γ_p(μ)

    For R̃_p = 1/μ:  c_s²/γ_p = 1  →  c_s = √γ_p  (PRIME-DEPENDENT!)
    The "string constant" c_s would need to depend on p → no universal T-duality.

  CONCLUSION:
    Standard T-duality does NOT give the required L̃_p = 1/μ (universal).
    The echo product structure requires a NON-STANDARD duality:
      either a MODULAR transformation specific to PT,
      or a DIFFERENT INTERNAL SPACE (not the Bianchi I T³).

  STATUS OF β_c² = 3 × s × echo:
    Factor 3:    from classical Bianchi I KK (d=3)     [DERIVED ✓]
    Factor s:    from Z/2Z / T1 projection             [ALGEBRAIC ✓]
    Factor echo: from BPS instanton product            [OPEN — requires
                  non-standard T-duality or prime-lattice modular structure]

  BEST CANDIDATE FOR NEXT STEP:
    Define the "PT modular transformation" T_PT on T³ that maps:
      a_p(μ) ↦ 1/(p × μ)   [modular dual of the Bianchi I metric]
    and show it generates the echo action S_p = p/μ.
    This would close the gap with a purely PT-algebraic argument.
""")
    print("=" * 65)
    print("Script completed.")
    print("=" * 65)


# ─── MAIN ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    route_I()
    route_II()
    route_III()
    key_diagnostic()
    synthesis()
