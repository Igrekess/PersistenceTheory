"""
research_C2 — Step A: numerical scan of N(2) candidates and structural
verification of (p+1)^(p+1) - 1.

Goal of C2 (cf. ch23 sec:ch23_commitments, prop:Np in ch10):
    Promote the structural commitment "N(2) = 26 = (p+1)^(p+1) - 1"
    from [DER] to [THM] via a categorical universality argument
    that fixes the domain to [p+1] (rather than [p], [p-1], etc.).

This script establishes Step A:
  1. Confirm numerically that N=26 minimises the F(2) residual against
     the CODATA value of α_EM, in the window N ∈ [20, 31].
  2. Verify that 26 = 3^3 - 1 = (p+1)^(p+1) - 1 for p=2.
  3. Show that the alternatives in {20, 24, 25, 26, 27, 28, 30, 32}
     can be classified by their factorisation:
       - 26 = 3^3 - 1 (endofunctions on [3] minus identity) — PT
       - 27 = 3^3 (includes identity)
       - 24 = 4! (permutations on [4])
       - 25 = 5^2 (square)
       - 30 = 2·3·5 (primorial)
     Only 26 has the "endofunctions on [p+1] minus identity" reading.
  4. Report ppm errors for each candidate.

The remaining gap (Step B, treated in 02_categorical_universality.py)
is the categorical / universal-property argument for the domain choice
[p+1].
"""

from __future__ import annotations
import math
import sympy as sp
from fractions import Fraction


# ---------------------------------------------------------------------------
# 1. PT formula F(p) = sin²θ_p · cos²(θ_p/N(p)) · (μ* - p)/p²
# ---------------------------------------------------------------------------
MU_STAR = sp.Integer(15)  # T7 fixed point
P_PRIME = sp.Integer(2)   # binary channel
ALPHA_EM_CODATA = sp.Rational(1, 137035999084) * sp.Integer(10**9)  # exact ratio

# q = (μ* - 2)/μ* = 13/15 ; δ_p = (1 - q^p)/p
q_stat = sp.Rational(13, 15)


def delta_p(p: int) -> sp.Rational:
    """δ_p = (1 - q^p)/p with q = 13/15 = q_stat at μ* = 15."""
    return (1 - q_stat**p) / p


def cos_theta_p(p: int) -> sp.Rational:
    """cos θ_p = 1 - δ_p (exact rational since q_stat is rational)."""
    return 1 - delta_p(p)


def sin2_theta_p(p: int) -> sp.Rational:
    """sin²θ_p = δ_p (2 - δ_p) (exact rational)."""
    d = delta_p(p)
    return d * (2 - d)


def F_p_with_N(p: int, N: int) -> sp.Float:
    """F(p) = sin²θ_p · cos²(θ_p/N) · (μ* - p)/p²."""
    s2 = sin2_theta_p(p)
    cos_p = cos_theta_p(p)
    theta_p = sp.acos(cos_p)
    theta_over_N = theta_p / N
    cos2_factor = sp.cos(theta_over_N) ** 2
    depth = (MU_STAR - p) / sp.Integer(p) ** 2
    return s2 * cos2_factor * depth


# ---------------------------------------------------------------------------
# 2. Habillage: 1/α_EM = 1/α_bare + F(2) — solve for F(2) numerically
# ---------------------------------------------------------------------------
# α_bare from product of sin²θ_p over active primes {3, 5, 7}
def alpha_bare() -> sp.Float:
    """α_bare = Π sin²θ_p for p ∈ {3, 5, 7}."""
    out = sp.Integer(1)
    for p in (3, 5, 7):
        out = out * sin2_theta_p(p)
    return out


def predicted_inv_alpha(N_value: int) -> sp.Float:
    """1/α_EM = 1/α_bare + F(2; N)."""
    f2 = F_p_with_N(2, N_value).evalf(50)
    a_bare = alpha_bare().evalf(50)
    return (1 / a_bare + f2)


# ---------------------------------------------------------------------------
# 3. Scan N ∈ [20, 31] — confirm N=26 is the residual minimum
# ---------------------------------------------------------------------------
def scan_N_window(N_min: int = 20, N_max: int = 32) -> list[dict]:
    """Compute 1/α(N) and ppm error for each N candidate."""
    inv_alpha_codata = sp.Rational(137035999084, 10**9)  # CODATA value
    results = []
    for N in range(N_min, N_max + 1):
        inv_alpha_pred = predicted_inv_alpha(N)
        err_abs = sp.Abs(inv_alpha_pred - inv_alpha_codata)
        err_ppm = (err_abs / inv_alpha_codata * sp.Integer(10**6)).evalf(8)
        results.append({
            "N": N,
            "inv_alpha_pred": float(inv_alpha_pred),
            "ppm_error": float(err_ppm),
            "structural_reading": classify(N),
        })
    return results


def classify(N: int) -> str:
    """Classify each candidate N by its arithmetic reading."""
    classifications = {
        20: "= 4·5 (no pointed-set reading)",
        21: "= 3·7 (= p₁·p₃, but not endofunction count)",
        22: "= 2·11 (no sieve reading)",
        23: "prime (no combinatorial reading)",
        24: "= 4! (permutations, not endofunctions)",
        25: "= 5² (square, not End count)",
        26: "= 3³−1 = End([3])\\{id} — PT canonical",
        27: "= 3³ = End([3]) (includes identity, no leakage from id)",
        28: "= 4·7 (no sieve reading)",
        29: "prime",
        30: "= 2·3·5 (primorial, wrong counting)",
        31: "prime",
        32: "= 2⁵",
    }
    return classifications.get(N, "?")


# ---------------------------------------------------------------------------
# 4. Verify N(p) = (p+1)^(p+1) - 1 for small p
# ---------------------------------------------------------------------------
def N_PT_formula(p: int) -> int:
    """N(p) = (p+1)^(p+1) - 1 = |End([p+1])| - 1."""
    return (p + 1) ** (p + 1) - 1


def end_count(n: int) -> int:
    """Number of endofunctions on a set of size n."""
    return n ** n


# ---------------------------------------------------------------------------
# 5. Main verification — print everything
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 78)
    print("research_C2 / Step A — endofunctions scan and N(2) verification")
    print("=" * 78)

    # 4a. PT formula sanity check
    print("\n[1] PT formula N(p) = (p+1)^(p+1) - 1 for p ∈ {1, 2, 3, 5}:")
    for p in (1, 2, 3, 5):
        nval = N_PT_formula(p)
        ec = end_count(p + 1)
        print(f"    p={p}: N(p) = {p+1}^{p+1} - 1 = {ec} - 1 = {nval}"
              f"   (= |End([{p+1}])| − 1)")

    # 4b. Scan
    print("\n[2] Scan of N ∈ [20, 32] against CODATA α_EM:")
    print(f"    {'N':>3}  {'1/α_pred':>18}  {'|err|/CODATA (ppm)':>20}  reading")
    print(f"    {'-'*3}  {'-'*18}  {'-'*20}  {'-'*40}")
    rows = scan_N_window(20, 32)
    rows_sorted = sorted(rows, key=lambda r: r["ppm_error"])
    for r in rows:
        flag = "  ← PT" if r["N"] == 26 else "  ← min" if r["ppm_error"] == rows_sorted[0]["ppm_error"] else ""
        print(f"    {r['N']:>3}  {r['inv_alpha_pred']:>18.10f}  "
              f"{r['ppm_error']:>20.4f}  {r['structural_reading']}{flag}")

    n26 = next(r for r in rows if r["N"] == 26)
    nmin = rows_sorted[0]
    print(f"\n    N=26 ppm error: {n26['ppm_error']:.4f}")
    print(f"    Min ppm error : N={nmin['N']} → {nmin['ppm_error']:.4f}")
    if nmin["N"] == 26:
        print("    ✓ N=26 is BOTH the structurally motivated AND the residual min.")
    else:
        print(f"    ⚠ N=26 is NOT the residual min — closest is N={nmin['N']}")
        print("      (this is why C2 is [DER]: domain choice is structural,")
        print("       not residual-driven; the residual scan is a posteriori check)")

    # 4c. Catalogue of competing candidates with their domain interpretations
    print("\n[3] Domain interpretations for competing endofunction counts:")
    for n_dom in (1, 2, 3, 4, 5):
        ec = end_count(n_dom)
        print(f"    |End([{n_dom}])|       = {n_dom}^{n_dom} = {ec:>5}   "
              f"(remove id ⟹ {ec-1})")

    print("\n[4] Why [p+1] for p=2 (the structural commitment):")
    print("    The pointed-set (Z/2Z, 0, 2) is built as:")
    print("      - the 2 residue classes {[0], [1]} = Z/2Z   (2 elements)")
    print("      - the operator '2' (= prime performing elimination)  (1 element)")
    print("      ⟹ pointed-set has |Z/2Z| + 1 = 3 elements,")
    print("        i.e. domain = [p+1] = [3].")
    print("    Endofunctions on [3] count 3³ = 27; minus identity = 26.")

    print("\n[5] Categorical reading (ANSATZ — Step B target):")
    print("    [p+1] should be the UNIVERSAL pointed-set in the category")
    print("    Set_*^{Z/pZ-mod} (pointed sets equipped with Z/pZ-action +")
    print("    distinguished operator). This requires:")
    print("      (i)   defining the appropriate slice category")
    print("      (ii)  proving [p+1] is the initial / terminal object")
    print("      (iii) connecting to cyclic operads (Loday-Vallette 2012)")
    print("    Estimation: 6-12 months (open).")
    print()


if __name__ == "__main__":
    main()
