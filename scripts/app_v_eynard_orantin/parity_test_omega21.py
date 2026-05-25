r"""
parity_test_omega21.py — test omega_{2,1}^{S,cusp} pour k=3 (cas PT canonique).

Conjecture : alpha_z((2,1)) = (3|S|+8)/2 = 17/2 pour k=3.
Valuation z_0 predite : 1 - 17 = -16.

omega_{2,1} = Res K(z_0, z) [omega_{1,2}(z, -z) + omega_{1,1}(z)*omega_{1,1}(-z)]

Nb : omega_12 et omega_11 sur la diagonale (z, -z) donnent des singularites
particulieres. Reduction symbolique probablement lourde.

Ne teste que k=3 (k pair degenere comme observe pour omega_{1,2}).
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from eo_recursion import SpectralCurve, omega_21


def make_cusp_curve(k: int):
    z = sp.Symbol('z')
    return SpectralCurve(
        x_of_z=z**2 / 2,
        y_of_z=z**k,
        z=z,
        branch_pts=[0],
        involution={0: -z},
    )


def extract_z0_valuation(expr: sp.Expr, z0: sp.Symbol) -> int:
    e = sp.together(expr)
    n, d = sp.fraction(e)
    try:
        n_poly = sp.Poly(sp.expand(n), z0)
        n_min = min(m[0] for m in n_poly.monoms()) if n_poly.terms() else 0
    except sp.PolynomialError:
        n_min = 0
    try:
        d_poly = sp.Poly(sp.expand(d), z0)
        d_min = min(m[0] for m in d_poly.monoms()) if d_poly.terms() else 0
    except sp.PolynomialError:
        d_min = 0
    return n_min - d_min


def main():
    print("=" * 76)
    print("  Test omega_{2,1} pour k=3 (PT canonique, y = z^3)")
    print("=" * 76)
    print(f"  Conjecture : val_z0 = -(3k + 7) = -16, alpha_z = 17/2")
    print()
    z0 = sp.Symbol('z_0')
    curve = make_cusp_curve(3)
    print("  Calcul omega_21(z_0) en cours...", flush=True)
    try:
        w21 = omega_21(curve, z0)
        val = extract_z0_valuation(w21, z0)
        val_pred = -16
        alpha_z = sp.Rational(1 - val, 2)
        match = "OUI" if val == val_pred else "NON"
        print(f"  val_z0 calc = {val}")
        print(f"  val_z0 pred = {val_pred}")
        print(f"  alpha_z     = {alpha_z}")
        print(f"  Match conjecture : {match}")
        if val != val_pred:
            print(f"  omega_{{2,1}}(z_0) leading = "
                  f"{sp.simplify(w21 * z0**(-val))} * z_0^{val}")
    except Exception as e:
        print(f"  ERREUR : {type(e).__name__}: {str(e)[:100]}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
