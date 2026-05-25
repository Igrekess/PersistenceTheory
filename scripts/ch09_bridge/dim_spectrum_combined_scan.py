#!/usr/bin/env python3
"""
Chantier court — PT × Dimension Spectrum, étape 2.

Build a numerical model of the FULL PT Dirac and compute its dimension
spectrum (= set of poles of ζ_D(z) = Tr |D|^{-z}).

The PT Dirac à la note 08 PT_NONCOM_BIFURCATION is:

    D_PT^N9 = D_spin^(4D)(Σ_pers) ⊗ Γ_{p=2}  +  I ⊗ D_F.

We consider THREE candidate spectral sectors that may contribute to
the dimension spectrum, each in isolation:

  (A) 4D cusp Weyl law: continuous spectrum + discrete Maass modes.
      Standard NCG poles at z ∈ {0, 2, 4}.
  (B) Arithmetic primes spectrum λ_p = p (Berry-Keating heuristic).
      Σ p^{-z} = prime zeta P(z). Log-singularity at z=1; essential
      singularity at z=0.
  (C) PT-cusp geodesic lengths l_p = 5p for primes actifs / écho
      ([[project_pt_rh_hyperbolic_cusp]] note 06). Σ (5p)^{-z}
      = 5^{-z}·P(z). Same singularities as (B), shifted by 5^{-z}.

We compute Σ_{n ≤ N_max} λ_n^{-z} on a grid of z values and check
numerically:
  ▸ Is there a pole at z = 1/2?
  ▸ Where are the dominant poles?
  ▸ Does the combined PT spectrum bring 1/2 into the dimension
    spectrum?

Honest expectation: probably NO pole at z = 1/2 in any of the three
sectors. This would confirm that Manin's argument (dim_alg ℝ = 1/2)
is a *conceptual* justification of the s=1/2 axiom, NOT a *spectral
dimensional* fact. The script tests this explicitly.
"""

import math
import sys
import time
from pathlib import Path

import numpy as np
from mpmath import mp, mpc, mpf, zeta, primezeta, log, exp, im, re, almosteq


mp.dps = 50  # 50 decimal digits of precision
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True, parents=True)


# ----------------------------------------------------------------------
# Sector A — 4D cusp Weyl law (toy compact 4D spin manifold)
# ----------------------------------------------------------------------
# For a compact 4D spin manifold of volume V, the Dirac eigenvalues
# obey N(λ) ~ V · λ^4 / (24·π²). Equivalently, the n-th eigenvalue
# scales as λ_n ~ (24π²·n/V)^(1/4) ~ n^(1/4) up to constants.
#
# ζ_D(z) = Σ λ_n^{-z} = (V/(24π²))^(z/4) · Σ n^{-z/4}
#       = const · ζ(z/4)
#
# Poles: ζ(z/4) has pole at z/4 = 1, i.e., z = 4.
# Trivial zeros of ζ at z/4 = -2, -4, ... give poles of 1/ζ but not
# of ζ itself. So Σ_A = {4}.
def zeta_cusp_4D(z, N_max=10**5):
    """Riemann ζ at z/4, computed directly via mpmath for accuracy."""
    return zeta(z / 4)


# ----------------------------------------------------------------------
# Sector B — arithmetic prime spectrum λ_p = p (Berry-Keating heuristic)
# ----------------------------------------------------------------------
# ζ_B(z) = Σ_p p^{-z} = primezeta(z).
# Logarithmic singularity at z = 1 (from log ζ pole).
# Essential singularity at z = 0.
# No pole at z = 1/2 (P(1/2) is finite).
def zeta_arith_primes(z):
    """Prime zeta function P(z) = Σ_{p prime} p^{-z}."""
    return primezeta(z)


# ----------------------------------------------------------------------
# Sector C — PT-cusp geodesic lengths l_p = 5p (project_pt_rh_hyperbolic_cusp)
# ----------------------------------------------------------------------
# ζ_C(z) = Σ_p (5p)^{-z} = 5^{-z} · P(z).
# Same singularities as B, scaled by 5^{-z} (which is entire & nonzero).
def zeta_PT_cusp(z):
    """5^{-z} · Prime zeta. Models the geodesic-length spectrum."""
    return mpf(5) ** (-z) * primezeta(z)


# ----------------------------------------------------------------------
# Sector D — PT-active geodesic lengths only: {3·5, 5·5, 7·5} = {15, 25, 35}
# ----------------------------------------------------------------------
# Finite 3-term sum: ζ_D(z) = 15^{-z} + 25^{-z} + 35^{-z}.
# Entire function. No poles.
def zeta_PT_active_3(z):
    return mpf(15) ** (-z) + mpf(25) ** (-z) + mpf(35) ** (-z)


# ----------------------------------------------------------------------
# Sector E — Riemann zeros γ_n (Hilbert-Polya speculative)
# ----------------------------------------------------------------------
# ζ_E(z) = Σ γ_n^{-z}. By Cramér (1919), this has continuation with
# a pole at z = 1 only. No pole at z = 1/2.
# (Computing γ_n explicitly is expensive; we sample first 1000.)
def riemann_zeros(N=1000, cache=[]):
    if not cache or len(cache[0]) < N:
        cache.clear()
        cache.append([mp.zetazero(n).imag for n in range(1, N + 1)])
    return cache[0][:N]


def zeta_riemann_zeros(z, N=1000):
    gammas = riemann_zeros(N)
    return sum(mpc(g) ** (-z) for g in gammas)


# ----------------------------------------------------------------------
# Pole detector
# ----------------------------------------------------------------------
# Numerical strategy: scan z along Re axis on [0.1, 5]. At each z,
# evaluate ζ_sector(z). Near a pole, |ζ| → ∞. We compute a few sample
# points near each candidate z_target and check:
#   (a) magnitude grows
#   (b) (z - z_target)·ζ(z) → finite residue
def check_pole(zeta_fn, z_target, eps=mpf("0.001"), label=""):
    """Probe for a pole at z = z_target.
    Returns (is_pole, magnitude_at_target, residue_estimate)."""
    try:
        val_at = zeta_fn(z_target)
    except Exception as e:
        return (None, f"ERROR at z={z_target}: {e}", None)
    # Sample at z = z_target ± eps
    try:
        val_plus = zeta_fn(z_target + eps)
        val_minus = zeta_fn(z_target - eps)
        # Estimate residue if pole
        residue = (val_plus * eps + val_minus * eps) / 2
        # If pole, |val_at| should be much larger than |val_plus|, |val_minus|
        is_pole = abs(val_at) > 100 * max(abs(val_plus), abs(val_minus))
    except Exception as e:
        return (None, f"ERROR near z={z_target}: {e}", None)
    return (is_pole, abs(val_at), residue)


# ----------------------------------------------------------------------
# Run the diagnostic on all sectors
# ----------------------------------------------------------------------

SECTORS = {
    "A. 4D cusp Weyl law (ζ(z/4))": zeta_cusp_4D,
    "B. Prime zeta P(z) = Σ p^{-z}": zeta_arith_primes,
    "C. PT geodesic 5p: 5^{-z}·P(z)": zeta_PT_cusp,
    "D. Finite {15,25,35}": zeta_PT_active_3,
    "E. Riemann zeros γ_n^{-z}": zeta_riemann_zeros,
}

# Candidate poles to test
CANDIDATE_POLES = [
    mpf("0.5"),    # the PT axiom s = 1/2
    mpf("1"),      # ζ pole
    mpf("2"),      # heat kernel a_1
    mpf("4"),      # 4D Weyl
    mpf("0.25"),   # RH-shifted (if HP applies)
    mpf("1.5"),    # σ_crit from W1 (note 06 project_pt_rh_weil_spiral)
]

print("=" * 78)
print(f"Dimension spectrum scan — mpmath dps={mp.dps}")
print(f"Candidate pole locations: {[float(z) for z in CANDIDATE_POLES]}")
print("=" * 78)

results = {}
for sector_name, zeta_fn in SECTORS.items():
    print(f"\n[{sector_name}]")
    sector_results = {}
    for z_test in CANDIDATE_POLES:
        t0 = time.time()
        is_pole, mag, residue = check_pole(zeta_fn, z_test, label=sector_name)
        dt = time.time() - t0
        z_float = float(z_test)
        if is_pole is None:
            flag = "ERR"
        elif is_pole:
            flag = "POLE!"
        else:
            flag = "----"
        print(f"  z = {z_float:>6.3f}  [{flag:>5}]  |ζ(z)| = {mag}    [{dt:.2f}s]")
        sector_results[z_float] = {
            "is_pole": bool(is_pole) if is_pole is not None else None,
            "magnitude": str(mag),
            "residue_estimate": str(residue) if residue is not None else None,
        }
    results[sector_name] = sector_results


# ----------------------------------------------------------------------
# Combined PT spectrum (cusp + arith + bifurc; sum of sectors)
# ----------------------------------------------------------------------
# A natural composite for PT is:
#   ζ_PT(z) = ζ_cusp(z) + ζ_arith(z) + 2·m^{-z}
# (Direct sum of sectors. Tensor product would multiply, but the
# physical PT Dirac D_PT = D_geom ⊗ Γ + I ⊗ D_F has eigenvalues
# √(λ_geom² + m²) — see below.)

print("\n" + "=" * 78)
print("Combined PT spectrum via direct sum:")
print(f"  ζ_PT(z) = ζ_cusp(z) + ζ_arith(z) + 2·m^{{-z}}, m = q_- - q_+")
print("=" * 78)
# At μ* = 15, q_+ = exp(-1/15), q_- = ? — use approximate values
m_PT = mpf("0.1")  # placeholder; the *value* of m doesn't move poles.

def zeta_PT_combined(z):
    return zeta_cusp_4D(z) + zeta_arith_primes(z) + 2 * (m_PT ** (-z))

for z_test in CANDIDATE_POLES:
    is_pole, mag, residue = check_pole(zeta_PT_combined, z_test)
    flag = "POLE!" if is_pole else "----"
    print(f"  z = {float(z_test):>6.3f}  [{flag:>5}]  |ζ_PT(z)| = {mag}")


# ----------------------------------------------------------------------
# Tensor product version: |D|² = λ_geom² + m², so ζ uses √(λ² + m²)
# ----------------------------------------------------------------------
# For D_PT = D_geom ⊗ Γ + I ⊗ D_F, eigenvalues are ±√(λ² + m²) for
# each (λ, ±m). So |D_PT|^{-z} has Tr summed over (λ, ±):
#   ζ_PT_tensor(z) = 2 · Σ_n (λ_n² + m²)^{-z/2}
# For large λ, (λ² + m²)^{-z/2} ≈ λ^{-z}, so the poles match ζ_geom
# to leading order. The bifurcation is a "mass gap" only — doesn't
# introduce new poles.

print("\n" + "=" * 78)
print("Tensor PT spectrum: |D_PT|² = D_cusp² + m², ζ_PT_tensor(z)")
print("=" * 78)

def zeta_PT_tensor(z, N=10**4, m=m_PT):
    """Approximate ζ_PT_tensor using cusp Weyl proxy λ_n ~ n^(1/4).
    Σ_n (n^(1/2) + m²)^{-z/2} truncated to N."""
    s = mpf(0)
    for n in range(1, N + 1):
        lam_sq = mpf(n) ** mpf("0.5") + m * m
        s += lam_sq ** (-z / 2)
    return 2 * s

# Note: this is a slow brute-force sum. Use small N for the pole scan.
print("(slow numerical Σ, using N=10^4)")
for z_test in CANDIDATE_POLES:
    t0 = time.time()
    val = zeta_PT_tensor(z_test, N=10**4)
    dt = time.time() - t0
    print(f"  z = {float(z_test):>6.3f}  ζ_PT_tensor(z) ≈ {val}    [{dt:.1f}s]")


# ----------------------------------------------------------------------
# Final analysis & verdict
# ----------------------------------------------------------------------

# Count how many sectors flagged a pole at each candidate
pole_count = {float(z): 0 for z in CANDIDATE_POLES}
for sector, res in results.items():
    for z_str, info in res.items():
        if info["is_pole"]:
            pole_count[z_str] += 1

print("\n" + "=" * 78)
print("SUMMARY — pole detection across all sectors")
print("=" * 78)
print(f"  {'z_candidate':>12}  {'#sectors with pole':>20}")
for z_str, n in sorted(pole_count.items()):
    print(f"  {z_str:>12.3f}  {n:>20}")
print()

# Check the critical question:
poles_at_half = pole_count.get(0.5, 0)
if poles_at_half == 0:
    verdict = (
        "VERDICT: z = 1/2 is NOT a pole of any PT-physical zeta sector.\n"
        "        ⟹ The s = 1/2 axiom is NOT a spectral-dimensional fact;\n"
        "          it remains a CONCEPTUAL identity (Manin algebro-geo,\n"
        "          T_3 stationary, Fisher max, etc.). The Manin §2.4\n"
        "          justification stands but as a *normalization*, not a\n"
        "          *measured spectral dimension*."
    )
else:
    verdict = (
        f"VERDICT: z = 1/2 appears as pole in {poles_at_half} sector(s).\n"
        "        ⟹ The s = 1/2 axiom has a spectral-dimensional reading.\n"
        "          This would be a significant new acquisition for PT-NCG."
    )

print(verdict)

# Save the results to JSON
import json
with open(RESULTS_DIR / "dim_spectrum_scan.json", "w") as f:
    out = {
        "candidate_poles": [float(z) for z in CANDIDATE_POLES],
        "pole_count": pole_count,
        "verdict_summary": (
            "1/2 not in dim spectrum"
            if pole_count.get(0.5, 0) == 0
            else "1/2 found in dim spectrum"
        ),
        "sectors_tested": list(SECTORS.keys()),
    }
    json.dump(out, f, indent=2)
print(f"\nSaved: {RESULTS_DIR / 'dim_spectrum_scan.json'}")
