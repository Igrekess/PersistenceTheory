"""
validate_L_squared.py — Validation symbolique du Théorème L_S² universel.

Théorème (app_v Eynard-Orantin universelle, thm:app_eo_L_squared) :
    L_S² = 4^{|S|-1} / (mu_S^{|S|} * pi²)
pour les cinq solutions ACA listées dans le Tableau associé.

Vérifie que la formule générale et la forme simplifiée du tableau
coïncident pour chaque solution S, via arithmétique rationnelle exacte
(sympy.Rational), zéro précision flottante.

Référence: PT_MONOGRAPHY/appendices_fr/app_v_eynard_orantin_universelle.tex,
Tableau tab:app_eo_universality.

Statut: validation symbolique exacte, 5/5 PASS.
"""

from __future__ import annotations

from fractions import Fraction


def L_squared_formula(S: tuple[int, ...]) -> Fraction:
    """L_S² = 4^{|S|-1} / mu_S^{|S|} en rationnel exact (modulo pi²)."""
    abs_S = len(S)
    mu_S = sum(S)
    return Fraction(4 ** (abs_S - 1), mu_S ** abs_S)


def main() -> int:
    # Cinq solutions ACA (Active Connected Admissible) du Tableau
    # tab:app_eo_universality, app_v_eynard_orantin_universelle.tex
    # Format : (S, mu_S attendu, L_S² formule attendue numérateur, dénominateur)
    cases = [
        # (S, mu_S, numerateur, denominateur)  -- pour L_S² = num / (denom * pi²)
        ((1, 2, 3),    6,  16,    216),     # solution simple
        ((2, 3, 5),    10, 16,   1000),
        ((3, 5, 7),    15, 16,   3375),     # S^star = solution privilégiée PT
        ((2, 3, 5, 7), 17, 64,  83521),
        ((1, 3, 7, 9), 20, 64, 160000),
    ]

    passed = 0
    total = 0
    print("=" * 60)
    print("VALIDATION L_S² universel — app_v thm:app_eo_L_squared")
    print("=" * 60)

    for S, mu_expected, num_expected, denom_expected in cases:
        total += 1
        L2 = L_squared_formula(S)
        L2_expected = Fraction(num_expected, denom_expected)

        # Vérification mu_S
        mu_actual = sum(S)
        mu_ok = mu_actual == mu_expected

        # Vérification L_S² exact
        L2_ok = L2 == L2_expected

        status = "PASS" if (mu_ok and L2_ok) else "FAIL"
        if mu_ok and L2_ok:
            passed += 1

        S_str = "{" + ",".join(map(str, S)) + "}"
        print(
            f"  S={S_str:14s}  |S|={len(S)}  mu_S={mu_actual:3d}  "
            f"L_S²={L2.numerator}/{L2.denominator}*pi^{{-2}}  [{status}]"
        )

    print("=" * 60)
    print(f"TOTAL: {passed}/{total} PASS")
    print("=" * 60)
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
