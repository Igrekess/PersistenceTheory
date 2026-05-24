#!/usr/bin/env python3
"""P2 - test the three candidate routes for deriving pi^2/240 from PT.

Route A: Pontryagin direct
    pi^2/240 ?= alpha_nu * 2*pi * (1 - beta_echo) ?
    (alpha_nu = prod sin^2(theta_p, q_+) over actives {3,5,7})

Route B: Fisher-Casimir
    Hess(ln Z_Ruelle) at the boundary
    (numerical, requires Z_Ruelle computation)

Route C: zeta(4) primorial split
    pi^2/240 = 3 zeta(4) / (8 pi^2)
    zeta(4) = prod_p (1 - p^{-4})^{-1}, split into actives {3,5,7} x echo
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# PT canonical values at mu* = 15
MU_STAR = 15
Q_PLUS = Fraction(13, 15)              # q_+ = 1 - 2/mu* = 13/15
Q_MINUS = math.exp(-1.0 / MU_STAR)     # q_- = exp(-1/mu*)

# Active primes
ACTIVES = [3, 5, 7]
# Ghost primes (cosmologically inactive but present)
ECHOES = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Target coefficient
TARGET_PI2_240 = math.pi ** 2 / 240.0


def delta_p(p: int, q_value: float) -> float:
    """Holonomy delta_p = (1 - q^p) / p."""
    return (1.0 - q_value ** p) / p


def sin2_theta(p: int, q_value: float) -> float:
    """sin^2(theta_p) = delta_p * (2 - delta_p)."""
    d = delta_p(p, q_value)
    return d * (2.0 - d)


def gamma_p_at_mu(p: int, mu: float, q_value: float) -> float:
    """Anomalous dimension gamma_p at scale mu."""
    q = q_value
    d = delta_p(p, q)
    if d == 0:
        return 0.0
    return (4.0 * p * q ** (p - 1) * (1.0 - d)) / (mu * (1.0 - q ** p) * (2.0 - d))


def euler_factor_inv(p: int, s: int) -> Fraction:
    """Inverse of Euler factor at p: (1 - p^{-s})."""
    return Fraction(1) - Fraction(1, p ** s)


def zeta_partial_product(primes: list[int], s: int) -> float:
    """Product of (1 - p^{-s})^{-1} over given primes."""
    total = 1.0
    for p in primes:
        total /= float(euler_factor_inv(p, s))
    return total


def route_A_pontryagin() -> dict:
    """Route A: pi^2/240 ?= alpha_nu * 2pi * (1 - beta_echo)."""
    q_plus = float(Q_PLUS)

    # alpha_nu = product sin^2(theta_p, q_+) over actives
    alpha_nu = 1.0
    for p in ACTIVES:
        alpha_nu *= sin2_theta(p, q_plus)

    # beta_echo = sum sin^2(theta_p) * gamma_p over {11, 13}
    beta_echo = 0.0
    for p in [11, 13]:
        beta_echo += sin2_theta(p, q_plus) * gamma_p_at_mu(p, MU_STAR, q_plus)

    candidate = alpha_nu * 2.0 * math.pi * (1.0 - beta_echo)
    return {
        "route": "A: Pontryagin alpha_nu * 2pi * (1 - beta_echo)",
        "alpha_nu": alpha_nu,
        "beta_echo": beta_echo,
        "1_minus_beta_echo": 1.0 - beta_echo,
        "candidate": candidate,
        "target": TARGET_PI2_240,
        "rel_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240,
        "ppm_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240 * 1e6,
    }


def route_A2_alternative() -> dict:
    """Route A variant: pi^2/240 ?= alpha_nu * 2pi (no echo correction)."""
    q_plus = float(Q_PLUS)
    alpha_nu = 1.0
    for p in ACTIVES:
        alpha_nu *= sin2_theta(p, q_plus)
    candidate = alpha_nu * 2.0 * math.pi
    return {
        "route": "A2: alpha_nu * 2pi (no echo)",
        "alpha_nu": alpha_nu,
        "candidate": candidate,
        "target": TARGET_PI2_240,
        "rel_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240,
        "ppm_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240 * 1e6,
    }


def route_C_zeta_split() -> dict:
    """Route C: pi^2/240 = 3 zeta(4) / (8 pi^2), zeta(4) decomposed."""
    # Exact: pi^2/240 = 3 * zeta(4) / (8 * pi^2)
    # zeta(4) = pi^4/90
    # Identity check: 3 * pi^4/90 / (8 pi^2) = pi^2 * 3 / 720 = pi^2 / 240. OK.

    zeta4_exact = math.pi ** 4 / 90.0

    # Split active vs echo
    zeta4_actives = zeta_partial_product(ACTIVES, 4)         # {3, 5, 7}
    zeta4_p2 = zeta_partial_product([2], 4)                   # parity channel
    zeta4_echo_p11plus = zeta_partial_product(ECHOES, 4)     # {11, 13, ...}

    zeta4_reconstructed = zeta4_actives * zeta4_p2 * zeta4_echo_p11plus

    # Tail correction (primes >= 47)
    tail = zeta4_exact / zeta4_reconstructed

    # Casimir coefficient via this decomposition
    coefficient = 3.0 * zeta4_exact / (8.0 * math.pi ** 2)

    # PT-meaning ratio: actives + p=2 / total
    pt_fraction = (zeta4_actives * zeta4_p2) / zeta4_exact

    return {
        "route": "C: pi^2/240 = 3 zeta(4) / (8 pi^2)",
        "zeta4_exact": zeta4_exact,
        "zeta4_actives_357": zeta4_actives,
        "zeta4_p2_parity": zeta4_p2,
        "zeta4_echo_p_ge_11_truncated": zeta4_echo_p11plus,
        "zeta4_reconstructed_truncated": zeta4_reconstructed,
        "tail_correction_factor": tail,
        "active_extended_fraction_of_zeta4": pt_fraction,
        "coefficient": coefficient,
        "target": TARGET_PI2_240,
        "rel_error": abs(coefficient - TARGET_PI2_240) / TARGET_PI2_240,
        "ppm_error": abs(coefficient - TARGET_PI2_240) / TARGET_PI2_240 * 1e6,
    }


def route_C2_pure_actives() -> dict:
    """Route C variant: ignore p=2 channel and echo. Use ONLY actives {3,5,7}."""
    zeta4_actives_only = zeta_partial_product(ACTIVES, 4)
    # Casimir candidate via actives-only:
    candidate = 3.0 * zeta4_actives_only / (8.0 * math.pi ** 2)
    return {
        "route": "C2: only actives {3,5,7} in zeta(4)",
        "zeta4_actives_only": zeta4_actives_only,
        "candidate": candidate,
        "target": TARGET_PI2_240,
        "rel_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240,
        "ppm_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240 * 1e6,
    }


def route_C3_with_p2() -> dict:
    """Route C variant: actives {3,5,7} + parity {2} (no echo)."""
    zeta4_with_p2 = zeta_partial_product([2] + ACTIVES, 4)
    candidate = 3.0 * zeta4_with_p2 / (8.0 * math.pi ** 2)
    return {
        "route": "C3: {2,3,5,7} (actives + parity, no echo)",
        "zeta4_2357": zeta4_with_p2,
        "candidate": candidate,
        "target": TARGET_PI2_240,
        "rel_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240,
        "ppm_error": abs(candidate - TARGET_PI2_240) / TARGET_PI2_240 * 1e6,
    }


def echo_tail_sum(s: int, p_min: int = 11, p_max: int = 100000) -> dict:
    """Compute sum_{p prime, p_min <= p <= p_max} 1/p^s exactly.

    Used to verify that the C3 residue equals the echo-tail contribution.
    """
    primes_in_range = []
    for n in range(p_min, p_max + 1):
        # quick primality
        if n < 2:
            continue
        if n == 2:
            primes_in_range.append(2)
            continue
        if n % 2 == 0:
            continue
        is_p = True
        d = 3
        while d * d <= n:
            if n % d == 0:
                is_p = False
                break
            d += 2
        if is_p:
            primes_in_range.append(n)

    tail_sum = sum(1.0 / (p ** s) for p in primes_in_range)
    log_tail = sum(-math.log1p(-1.0 / (p ** s)) for p in primes_in_range)
    multiplicative_correction = math.exp(log_tail)
    return {
        "s": s,
        "p_min": p_min,
        "p_max": p_max,
        "n_primes_in_tail": len(primes_in_range),
        "tail_sum_1_over_p_s": tail_sum,
        "tail_log_correction": log_tail,
        "tail_multiplicative_correction": multiplicative_correction,
        "tail_in_ppm": (multiplicative_correction - 1.0) * 1e6,
    }


def cross_check_1D_scalar() -> dict:
    """Cross-check Route C3 on the 1D scalar coefficient -pi/(24a).

    Standard derivation:
        -pi/(24a) = -zeta(2) / (4 pi a),
        zeta(2)  = pi^2 / 6 = prod_p (1 - p^{-2})^{-1}.

    PT-native truncation: keep only {2, 3, 5, 7} in the Euler product.
    """
    target = math.pi / 24.0
    zeta2_exact = math.pi ** 2 / 6.0

    # C3 truncation: actives + parity, no echo
    zeta2_C3 = zeta_partial_product([2] + ACTIVES, 2)
    candidate_C3 = zeta2_C3 / (4.0 * math.pi)

    # Variants
    zeta2_actives_only = zeta_partial_product(ACTIVES, 2)
    candidate_C2 = zeta2_actives_only / (4.0 * math.pi)

    # Exact reference
    candidate_exact = zeta2_exact / (4.0 * math.pi)

    return {
        "dimension": "1D scalar",
        "target_pi_over_24": target,
        "zeta2_exact": zeta2_exact,
        "zeta2_actives_only_357": zeta2_actives_only,
        "zeta2_C3_2357": zeta2_C3,
        "candidate_C2_actives_only": candidate_C2,
        "candidate_C3_2357": candidate_C3,
        "candidate_exact": candidate_exact,
        "C2_ppm_error": abs(candidate_C2 - target) / target * 1e6,
        "C3_ppm_error": abs(candidate_C3 - target) / target * 1e6,
        "exact_ppm_error": abs(candidate_exact - target) / target * 1e6,
    }


def prefactor_pt_identification() -> dict:
    """Test the PT-native identification of 3/(8 pi^2) prefactor in 3D Casimir.

    Hypothesis:
        3 / (8 pi^2) = N_spatial / (2^N_spatial * pi^N_spatial)
                     = 3 / (8 * pi^3) ... NOT MATCHING.

    Alternative:
        3 / (8 pi^2) = N_spatial / (2 * (2 pi)^2)
                     = 3 / (2 * 4 pi^2)
                     = 3 / (8 pi^2)  YES.

    So: prefactor = N_spatial / (2 * (2 pi)^N_spatial_in_loop)
                  with N_spatial_in_loop = 2 (transverse modes for parallel plates).

    Or: prefactor = (number of polarizations / 2) * (1 / (2 pi)^N_transverse).
    """
    N_spatial = 3                              # active primes {3, 5, 7}
    N_transverse = 2                           # transverse to plates
    N_polarizations = 2                        # EM has 2 polarizations
    candidate1 = N_spatial / (2.0 * (2.0 * math.pi) ** N_transverse)
    target = 3.0 / (8.0 * math.pi ** 2)
    return {
        "candidate_id": "N_spatial / (2 * (2 pi)^N_transverse)",
        "N_spatial": N_spatial,
        "N_transverse": N_transverse,
        "N_polarizations": N_polarizations,
        "candidate_value": candidate1,
        "target_3_over_8_pi2": target,
        "match": abs(candidate1 - target) < 1e-15,
        "interpretation": (
            "3 = N_spatial (number of active primes = number of spatial dimensions); "
            "(2 pi)^2 = (perimeter of S^1)^2 = limit of CRT torus area for transverse modes; "
            "factor 1/2 = standard zero-point energy normalization."
        ),
    }


def main() -> None:
    routes = [
        route_A_pontryagin(),
        route_A2_alternative(),
        route_C_zeta_split(),
        route_C2_pure_actives(),
        route_C3_with_p2(),
    ]

    # Question 1: echo tail sum verification
    echo_3D = echo_tail_sum(s=4, p_min=11)   # for 3D (zeta(4))
    echo_1D = echo_tail_sum(s=2, p_min=11)   # for 1D (zeta(2))

    # Question 2: prefactor identification
    prefactor = prefactor_pt_identification()

    # Question 3: 1D cross-check
    one_d = cross_check_1D_scalar()

    payload = {
        "status": "P2-coefficient-routes-test",
        "target_pi2_240": TARGET_PI2_240,
        "mu_star": MU_STAR,
        "q_plus_exact": "13/15",
        "q_minus": Q_MINUS,
        "actives": ACTIVES,
        "echo_pool": ECHOES,
        "routes": routes,
        "echo_tail_3D_zeta4": echo_3D,
        "echo_tail_1D_zeta2": echo_1D,
        "prefactor_pt_identification": prefactor,
        "cross_check_1D_scalar": one_d,
    }

    output = ROOT / "results" / "p2_coefficient_routes.json"
    output_md = ROOT / "results" / "p2_coefficient_routes.md"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    # Render markdown
    lines = [
        "# P2 - Coefficient Routes Test",
        "",
        f"Target: pi^2/240 = {TARGET_PI2_240:.16f}",
        f"mu* = {MU_STAR}, q_+ = 13/15, actives = {ACTIVES}",
        "",
        "## Summary",
        "",
        "| Route | Candidate | Target | Rel error | ppm error |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in routes:
        cand = r.get("candidate") or r.get("coefficient")
        lines.append(
            f"| {r['route']} | {cand:.10f} | {r['target']:.10f} | "
            f"{r['rel_error']:.3e} | {r['ppm_error']:.1f} |"
        )

    lines.extend(["", "## Detailed results", ""])
    for r in routes:
        lines.append(f"### {r['route']}")
        lines.append("")
        lines.append("```json")
        lines.append(json.dumps({k: v for k, v in r.items() if k != "route"}, indent=2))
        lines.append("```")
        lines.append("")

    # Append the 3 deeper-question answers
    lines.extend([
        "## Question 1: echo tail sum verification",
        "",
        "### 3D (zeta(4) tail)",
        "",
        "```json",
        json.dumps(echo_3D, indent=2),
        "```",
        "",
        "### 1D (zeta(2) tail)",
        "",
        "```json",
        json.dumps(echo_1D, indent=2),
        "```",
        "",
        "## Question 2: prefactor PT identification",
        "",
        "```json",
        json.dumps(prefactor, indent=2),
        "```",
        "",
        "## Question 3: cross-check on 1D scalar coefficient -pi/(24a)",
        "",
        "```json",
        json.dumps(one_d, indent=2),
        "```",
    ])
    output_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Wrote {output_md}")
    print()
    print("=" * 70)
    print(f"Target pi^2/240 = {TARGET_PI2_240:.10f}")
    print("=" * 70)
    for r in routes:
        cand = r.get("candidate") or r.get("coefficient")
        print(f"  {r['route']:50s} -> {cand:.10f}  err={r['ppm_error']:.1f} ppm")
    print()
    print("=" * 70)
    print("DEEPER QUESTIONS")
    print("=" * 70)
    print(f"Q1. 3D echo tail (zeta(4)):  multiplicative {echo_3D['tail_multiplicative_correction']:.6e}")
    print(f"    in ppm: {echo_3D['tail_in_ppm']:.1f} ppm  (vs C3 residue 131 ppm)")
    print(f"    1D echo tail (zeta(2)):  multiplicative {echo_1D['tail_multiplicative_correction']:.6e}")
    print(f"    in ppm: {echo_1D['tail_in_ppm']:.1f} ppm")
    print()
    print(f"Q2. prefactor 3/(8 pi^2) = {prefactor['target_3_over_8_pi2']:.10f}")
    print(f"    candidate {prefactor['candidate_id']} = {prefactor['candidate_value']:.10f}")
    print(f"    match: {prefactor['match']}")
    print()
    print(f"Q3. 1D scalar -pi/(24a), target |pi/24| = {one_d['target_pi_over_24']:.10f}")
    print(f"    C2 (actives only):   {one_d['candidate_C2_actives_only']:.10f}  err={one_d['C2_ppm_error']:.1f} ppm")
    print(f"    C3 (2357, no echo): {one_d['candidate_C3_2357']:.10f}  err={one_d['C3_ppm_error']:.1f} ppm")
    print(f"    exact (all primes):  {one_d['candidate_exact']:.10f}  err={one_d['exact_ppm_error']:.1f} ppm")


if __name__ == "__main__":
    main()
