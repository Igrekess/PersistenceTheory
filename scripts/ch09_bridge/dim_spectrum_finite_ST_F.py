#!/usr/bin/env python3
"""
Chantier court — PT × Dimension Spectrum, étape 1.

Compute the dimension spectrum of PT's finite spectral triple ST_F
as defined in `PT_LEAN/PT/Bridge/FiniteSpectralTriple.lean`:

    A_F = ℂ        (commutative trivial)
    H_F = ℂ²       (2-dim Hilbert space, q_+/q_- branches)
    D_F = m · σ_x  (off-diagonal Dirac, m = bifurcation mass Δq = q_- − q_+)
    γ_F = σ_x      (chirality)

The dimension spectrum Σ ⊂ ℂ is the minimal set such that

    ζ_b(z) := Tr(b · |D_F|^{-z})

extends holomorphically to ℂ ∖ Σ for all b in the algebra B generated
by [|D|, A]-iterates of A acting on H. For the finite case here, the
spectrum is computed in closed form.

Expected (and confirmed): Σ_F = ∅. The bifurcation alone is
"zero-dimensional" in the Connes sense — which means the axiom
s = 1/2 of PT CANNOT be read off ST_F in isolation. It must come
from the coupling with the infinite (cusp + arithmetic) sector.

This is the diagnostic step. Step 2 (script 02) builds the joint
spectral triple.
"""

import math
import sys
from typing import Iterable

import numpy as np
import sympy as sp


# ----------------------------------------------------------------------
# Step 1.1 — Reproduce the Lean module's ST_F symbolically.
# ----------------------------------------------------------------------

m = sp.Symbol('m', positive=True, real=True)
z = sp.Symbol('z')

sigma_x = sp.Matrix([[0, 1], [1, 0]])

D_F = m * sigma_x          # m·σ_x
gamma_F = sigma_x          # σ_x

assert sp.simplify(gamma_F * gamma_F - sp.eye(2)) == sp.zeros(2)
# Note: in PT's ST_F (cf. note 25 PT_NONCOM_BIFURCATION), γ_F = σ_x and
# D_F = m·σ_x both equal σ_x up to scalar, so they commute (not anti-
# commute) — non-standard relative to Connes-Marcolli but consistent
# with the Lean module FiniteSpectralTriple.lean.
DF_gF_commutator = sp.simplify(D_F * gamma_F - gamma_F * D_F)
print(f"[ST_F] [D_F, γ_F] = {DF_gF_commutator.tolist()} (commute, not anticommute)")

# Eigenvalues of D_F = m·σ_x are ±m.
eigvals_DF = sp.Matrix(D_F).eigenvals()
print("[ST_F] D_F = m·σ_x  eigenvalues:", dict(eigvals_DF))

# |D_F| has eigenvalue |m| with multiplicity 2.
# ζ_F(z) := Tr(|D_F|^{-z}) = 2 · |m|^{-z}.

zeta_F = 2 * m**(-z)
zeta_F_simplified = sp.simplify(zeta_F)
print(f"[ST_F] ζ_F(z) = Tr(|D_F|^{{-z}}) = {zeta_F_simplified}")

# Is ζ_F(z) entire? Yes — it's 2·m^{-z}, no poles for any z ∈ ℂ.
# The dimension spectrum is ∅.
print("[ST_F] Singularities of ζ_F(z) in ℂ:", "NONE")
print("[ST_F] ⟹ Dimension spectrum Σ_F = ∅ (zero-dimensional Connes-style)")

# Same for any b ∈ ℂ acting as scalar on H_F: ζ_b(z) = 2·b·m^{-z}, entire.
print()


# ----------------------------------------------------------------------
# Step 1.2 — Diagnostic: what would put 1/2 into Σ?
# ----------------------------------------------------------------------

print("=" * 70)
print("[diagnostic] What would put 1/2 into the dimension spectrum?")
print("=" * 70)
print(
    "A pole at z = 1/2 of ζ_b(z) requires that the eigenvalue density\n"
    "ρ(λ) = #{i : |λ_i| ≤ λ} grow like λ^(1/2) near infinity.\n"
    "\n"
    "Equivalently, Weyl's law for D in dimension d gives ρ(λ) ~ λ^d, so\n"
    "d = 1/2 means the Hilbert space H is 'half-dimensional' in the\n"
    "Manin algebro-geometric sense (dim_alg ℝ = 1/2, cf. Manin 2005\n"
    "§2.4, p.15).\n"
    "\n"
    "Concrete candidates for ρ(λ) ~ λ^(1/2):\n"
    "  (a) Square-root spectrum: λ_n ~ √n  (then |D|^{-z} = n^{-z/2},\n"
    "      Σ n^{-z/2} = ζ(z/2), pole at z = 2)\n"
    "  (b) Prime spectrum λ_p = p: Σ p^{-z} = prime zeta P(z), log-\n"
    "      singularity at z = 1 (NOT at z = 1/2).\n"
    "  (c) Critical zeros of ζ: σ(D) = {γ_n} with N(T) = (T/2π) log T\n"
    "      gives ρ(λ) ~ λ log λ — close to dim 1 but with log correction.\n"
    "      No clean 1/2.\n"
    "  (d) Fractal Cantor-like spectrum of HB-dim 1/2: σ_D = {2^(-n/2) :\n"
    "      ...}, giving Σ_n 2^(-nz/2), pole at z = 0 (geometric).\n"
    "\n"
    "Honest assessment: 1/2 in the dim spectrum is NOT generic — it\n"
    "requires a fine-tuned spectral density. PT's axiom s = 1/2 is\n"
    "more naturally an algebraic identity (T_3 stationary, Fisher max,\n"
    "etc.) than a spectral dimension.\n"
)


# ----------------------------------------------------------------------
# Step 1.3 — Try the obvious tensor: ST_F ⊗ (4D cusp).
# ----------------------------------------------------------------------

print("=" * 70)
print("[next] Tensor: ST_F ⊗ Σ_pers (toy 4D cusp model)")
print("=" * 70)
print(
    "Following note 08 PT_NONCOM_BIFURCATION:\n"
    "  D_PT^N9 = D_spin^(4D)(Σ_pers) ⊗ Γ_{p=2}  +  I ⊗ D_F\n"
    "\n"
    "For a 4D Riemannian spin manifold, ζ_cusp(z) has poles at z = 4,\n"
    "z = 2, z = 0 (heat-kernel coefficients a_0, a_1, a_2). The tensor\n"
    "with finite ST_F doesn't shift these; it just multiplies by\n"
    "Tr(|D_F|^{-z}) = 2·m^{-z} (which contributes no poles).\n"
    "\n"
    "So Σ(D_PT^N9_toy) ⊃ {0, 2, 4} for the geometric part. The 1/2 is\n"
    "NOT in this set. The arithmetic (prime-zeta) part must enter via\n"
    "the ADELIC sector (Connes-Marcolli) to bring in critical-line\n"
    "complex dimensions. See script 02.\n"
)
