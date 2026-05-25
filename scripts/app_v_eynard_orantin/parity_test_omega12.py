r"""
parity_test_omega12.py — calcul efficace de omega_{1,2}^{S,cusp} via le kernel
EO standard, qui a la meme forme que pour notre cusp local.

Strategie : utiliser eo_recursion.omega_12 avec y = c*z^k (k = |S|, c symbolique
ou =1 pour aller plus vite). Le kernel K = 1/(2 y (z_0^2 - z^2)) est le meme
pour x = z^2/2 et x = (1+z^2)/2 (seule dx = z dz importe).

Une fois omega_{1,2}(z_0, z_1) obtenu comme Taylor en z_0 (rationnel en z_1),
on extrait l'ordre principal en z_0 -> 0 avec z_1 fixe = 1 par exemple.

Conjecture testee : alpha_z((1,2)) = (2k+5)/2, valuation z_0 = 1 - (2k+5) = -(2k+4).
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from eo_recursion import SpectralCurve, omega_12, omega_11, omega_03


def make_cusp_curve(k: int):
    """Curve y = z^k (on x = z^2/2 with simple branch at z=0)."""
    z = sp.Symbol('z')
    return SpectralCurve(
        x_of_z=z**2 / 2,
        y_of_z=z**k,
        z=z,
        branch_pts=[0],
        involution={0: -z},
    )


def extract_z0_valuation(expr: sp.Expr, z0: sp.Symbol) -> int:
    """Extract leading (most negative) power of z0 in expr at z0 -> 0.

    Pour eviter sp.series qui hangs, on factorise / met sous forme rationnelle
    et on prend la difference des deg_min de num et den.
    """
    e = sp.together(expr)
    n, d = sp.fraction(e)
    # Lowest power of z0 in numerator
    n_expanded = sp.expand(n)
    try:
        n_poly = sp.Poly(n_expanded, z0)
        n_min = min(m[0] for m in n_poly.monoms()) if n_poly.terms() else 0
    except sp.PolynomialError:
        n_min = 0
    try:
        d_poly = sp.Poly(sp.expand(d), z0)
        d_min = min(m[0] for m in d_poly.monoms()) if d_poly.terms() else 0
    except sp.PolynomialError:
        d_min = 0
    return n_min - d_min


def test_omega_12(k_values, z1_value=1):
    """Pour chaque k, calcule omega_{1,2}(z_0, z_1=z1_value) et son ordre en z_0."""
    print("=" * 76)
    print(f"  Test omega_{{1,2}} via EO recursion sur cusp y = z^k, z_1 = {z1_value}")
    print("=" * 76)
    print(f"  Conjecture : val_z0(omega_{{1,2}}) = -(2k+4)")
    print(f"  Predit : k=3 -> -10, k=4 -> -12")
    print()
    print(f"  {'k':<8}{'val_z0 calc':<15}{'val_z0 pred':<15}{'alpha_z':<10}{'match':<8}")
    print("-" * 76)

    z0 = sp.Symbol('z_0')
    z1 = sp.Symbol('z_1')

    for k in k_values:
        curve = make_cusp_curve(k)
        try:
            print(f"  k={k} : calcul omega_12 en cours...", flush=True)
            w12 = omega_12(curve, z0, z1)
            # Substituer z1 = z1_value
            w12_at_z1 = w12.subs(z1, z1_value)
            # Extraire valuation z0
            val = extract_z0_valuation(w12_at_z1, z0)
            val_pred = -(2 * k + 4)
            alpha_z = sp.Rational(1 - val, 2)
            match = "OUI" if val == val_pred else "NON"
            print(f"  {k:<8}{val:<15}{val_pred:<15}{str(alpha_z):<10}{match:<8}")
            if val != val_pred:
                print(f"    omega_12(z_0, {z1_value}) leading = "
                      f"{sp.simplify(w12_at_z1 * z0**(-val))}")
        except Exception as e:
            print(f"  k={k} : ERREUR {type(e).__name__}: {str(e)[:60]}")
            import traceback
            traceback.print_exc()


if __name__ == '__main__':
    # Test k=3 first (PT canonical case)
    test_omega_12([3])
    print()
    test_omega_12([4])
