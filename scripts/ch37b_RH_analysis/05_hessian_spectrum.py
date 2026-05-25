"""
Script 05 — Spectre de la Hessienne D^2 F_m^PT(q^*), inversibilite (Z2a).

Le Hessien D^2 F_m^PT(q^*) est un operateur de multiplication sur C^0(K)
(cf. note 07 §2.4). On le represente comme la fonction H(mu) =
H_T5b(mu) + H_T1(mu) et on calcule :

  (i)   H_T5b(mu) = 2 (d_q Phi)^2 - 2 (mu - Phi) d_q^2 Phi
  (ii)  H_T1(mu) = 2 lambda mu^2
  (iii) min_K H, max_K H, conditionnement
  (iv)  dependance en epsilon et lambda

Validation : le Hessien diagonal donne directement les bornes BIFT du
theoreme Z2a (note 07 §2.6).

Usage :
    python 05_hessian_spectrum.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Callable

import numpy as np

ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]


def q_canonical(mu: float) -> float:
    return 1.0 - 2.0 / mu


def gamma_p_of_q(p: int, mu: float, q: float) -> float:
    if mu <= 2.01 or q <= 0.0 or q >= 1.0:
        return 0.0
    qp = q ** p
    d = (1.0 - qp) / p
    if d < 1e-15 or abs(2.0 - d) < 1e-15:
        return 0.0
    dln_delta = 2.0 * p * q ** (p - 1) / (mu * (1.0 - qp))
    factor = 2.0 * (1.0 - d) / (2.0 - d)
    return dln_delta * factor


def d_gamma_d_q(p: int, mu: float, q: float, h: float = 1e-6) -> float:
    """Numerical derivative d gamma_p / d q (pointwise in q)."""
    return (gamma_p_of_q(p, mu, q + h) - gamma_p_of_q(p, mu, q - h)) / (2.0 * h)


def d2_gamma_d_q2(p: int, mu: float, q: float, h: float = 1e-4) -> float:
    """Second derivative d^2 gamma_p / d q^2."""
    g_plus = gamma_p_of_q(p, mu, q + h)
    g_zero = gamma_p_of_q(p, mu, q)
    g_minus = gamma_p_of_q(p, mu, q - h)
    return (g_plus - 2.0 * g_zero + g_minus) / (h ** 2)


def S_eps(x: float, eps: float) -> float:
    return 0.5 * (1.0 + np.tanh(x / eps))


def S_eps_prime(x: float, eps: float) -> float:
    """S_eps'(x) = (1 / (2 eps)) sech^2(x/eps)."""
    return (1.0 / (2.0 * eps)) * (1.0 / np.cosh(x / eps)) ** 2


def S_eps_double_prime(x: float, eps: float) -> float:
    """S_eps''(x) = -(1/eps^2) sech^2(x/eps) tanh(x/eps)."""
    t = np.tanh(x / eps)
    sech2 = 1.0 / np.cosh(x / eps) ** 2
    return -(1.0 / eps ** 2) * sech2 * t


# Pointwise derivative d Phi / d q at q
def dPhi_dq(mu: float, q: float, m: int, eps: float) -> float:
    total = 0.0
    for p in ODD_PRIMES[:m]:
        g = gamma_p_of_q(p, mu, q)
        dg = d_gamma_d_q(p, mu, q)
        total += p * S_eps_prime(g - 0.5, eps) * dg
    return total


def d2Phi_dq2(mu: float, q: float, m: int, eps: float) -> float:
    """d^2 Phi / d q^2 (pointwise)."""
    total = 0.0
    for p in ODD_PRIMES[:m]:
        g = gamma_p_of_q(p, mu, q)
        dg = d_gamma_d_q(p, mu, q)
        d2g = d2_gamma_d_q2(p, mu, q)
        S1 = S_eps_prime(g - 0.5, eps)
        S2 = S_eps_double_prime(g - 0.5, eps)
        total += p * (S2 * dg ** 2 + S1 * d2g)
    return total


def Phi_value(mu: float, q: float, m: int, eps: float) -> float:
    total = 0.0
    for p in ODD_PRIMES[:m]:
        g = gamma_p_of_q(p, mu, q)
        total += p * S_eps(g - 0.5, eps)
    return total


# -------------------------------------------------------------------------
# Hessian density
# -------------------------------------------------------------------------

def hessian_density(mu: float, q: float, m: int, eps: float,
                    lam: float) -> dict:
    """Compute H(mu) = H_T5b(mu) + H_T1(mu) at q (which can be q^*)."""
    Phi = Phi_value(mu, q, m, eps)
    dPhi = dPhi_dq(mu, q, m, eps)
    d2Phi = d2Phi_dq2(mu, q, m, eps)
    H_T5b = 2.0 * dPhi ** 2 - 2.0 * (mu - Phi) * d2Phi
    H_T1 = 2.0 * lam * mu ** 2
    return {
        "mu": mu,
        "Phi": Phi,
        "dPhi_dq": dPhi,
        "d2Phi_dq2": d2Phi,
        "H_T5b": H_T5b,
        "H_T1": H_T1,
        "H_total": H_T5b + H_T1,
    }


def spectrum_on_K(m: int, eps: float, lam: float,
                  mu_grid: np.ndarray) -> dict:
    """Compute H(mu) over mu_grid at canonical q^*(mu)."""
    H_values = []
    for mu in mu_grid:
        q_star = q_canonical(float(mu))
        h = hessian_density(float(mu), q_star, m, eps, lam)
        H_values.append(h)

    H_T5b_arr = np.array([h["H_T5b"] for h in H_values])
    H_T1_arr = np.array([h["H_T1"] for h in H_values])
    H_arr = np.array([h["H_total"] for h in H_values])

    return {
        "m": m,
        "eps": eps,
        "lambda": lam,
        "H_T5b_min": float(np.min(H_T5b_arr)),
        "H_T5b_max": float(np.max(H_T5b_arr)),
        "H_T1_min": float(np.min(H_T1_arr)),
        "H_T1_max": float(np.max(H_T1_arr)),
        "H_min": float(np.min(H_arr)),
        "H_max": float(np.max(H_arr)),
        "kappa": float(np.max(H_arr) / np.max([np.min(H_arr), 1e-10])),
        "H_at_mu_15": next((h["H_total"] for h in H_values
                            if abs(h["mu"] - 15.0) < 0.1), None),
        "min_positive": bool(np.min(H_arr) > 0),
    }


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main() -> int:
    out_dir = Path(__file__).resolve().parent.parent / "outputs"
    out_dir.mkdir(exist_ok=True)

    mu_grid = np.linspace(3.0, 30.0, 100)  # coarser for speed
    m = 6  # m_K for K = [3, 30]

    print("=" * 72)
    print("Script 05 — Spectre de la Hessienne D^2 F_m^PT(q^*)")
    print("=" * 72)
    print(f"m = {m}, K = [3, 30] (100 points)")

    # Sweep over eps and lambda
    configurations = [
        (0.50, 1.0),
        (0.50, 1e4),
        (0.10, 1e4),
        (0.05, 1e4),
        (0.02, 1e4),
        (0.01, 1e4),
        (0.05, 1e2),
        (0.05, 1e6),
    ]

    print(f"\n{'eps':>6} {'lambda':>8} {'min H':>14} {'max H':>14} "
          f"{'kappa':>10} {'H(15)':>12} {'pos.?':>6}")
    print("-" * 78)
    results = []
    for eps, lam in configurations:
        r = spectrum_on_K(m, eps, lam, mu_grid)
        results.append(r)
        sign = "yes" if r["min_positive"] else "NO"
        print(f"{eps:>6.3f} {lam:>8.2e} {r['H_min']:>14.4e} "
              f"{r['H_max']:>14.4e} {r['kappa']:>10.2e} "
              f"{r['H_at_mu_15']:>12.4e} {sign:>6}")

    # Detailed dump at canonical config (eps=0.05, lambda=1e4)
    canonical = next(r for r in results
                     if abs(r["eps"] - 0.05) < 1e-6 and abs(r["lambda"] - 1e4) < 1)

    print("\n" + "=" * 72)
    print("Configuration canonique : eps = 0.05, lambda = 1e4")
    print("=" * 72)
    print(f"min H_T5b  : {canonical['H_T5b_min']:>14.4e}")
    print(f"max H_T5b  : {canonical['H_T5b_max']:>14.4e}")
    print(f"min H_T1   : {canonical['H_T1_min']:>14.4e}")
    print(f"max H_T1   : {canonical['H_T1_max']:>14.4e}")
    print(f"min H total: {canonical['H_min']:>14.4e}")
    print(f"max H total: {canonical['H_max']:>14.4e}")
    print(f"kappa      : {canonical['kappa']:>14.2e}")
    print(f"H(mu = 15) : {canonical['H_at_mu_15']:>14.4e}")

    # BIFT bounds
    r1 = 0.05  # heuristic from note 06
    r2 = r1 * canonical["H_min"]
    print("\nBIFT (Z2a) — bornes locales d'inversion :")
    print(f"  Rayon stabilite source (r_2) : ~{r2:.2e}")
    print(f"  Rayon stabilite cible  (r_1) : ~{r1:.4f}")

    # Asymptotic scaling check: H_T5b vs eps
    print("\n" + "=" * 72)
    print("Loi d'echelle : H_T5b ~ 1/eps^? (etude pour Z2b Nash-Moser)")
    print("=" * 72)
    eps_scan = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    H_T5b_at_15 = []
    for eps in eps_scan:
        h = hessian_density(15.0, q_canonical(15.0), m, eps, 0.0)
        H_T5b_at_15.append({"eps": eps, "H_T5b": h["H_T5b"],
                            "dPhi_dq": h["dPhi_dq"]})
        print(f"  eps = {eps:>6.3f}  H_T5b(15) = {h['H_T5b']:>14.4e}  "
              f"d Phi/d q (15) = {h['dPhi_dq']:>14.4e}")

    # Fit power law
    log_eps = np.log([s["eps"] for s in H_T5b_at_15])
    log_H = np.log([abs(s["H_T5b"]) + 1e-30 for s in H_T5b_at_15])
    if not np.any(np.isnan(log_H)):
        slope, intercept = np.polyfit(log_eps, log_H, 1)
        print(f"\n  Fit  log H_T5b = {slope:.3f} * log eps + {intercept:.3f}")
        print(f"  i.e. H_T5b ~ eps^({slope:.3f})  (Nash-Moser loss "
              f"order ~ {-slope:.2f})")

    # Save outputs
    payload = {
        "configurations": results,
        "canonical_config": canonical,
        "H_T5b_at_15_vs_eps": H_T5b_at_15,
    }
    out_file = out_dir / "05_hessian_spectrum.json"
    out_file.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\n[OK] Outputs written to {out_file.relative_to(out_dir.parent)}")

    # Verdict
    print("\n" + "=" * 72)
    print("VERDICT Z2a")
    print("=" * 72)
    all_positive = all(r["min_positive"] for r in results
                       if r["lambda"] >= 1e4)
    if all_positive:
        print("PASS  H(mu) > 0 sur tout K pour lambda >= 1e4 et eps testes.")
        print("PASS  Hessien inversible avec inverse borne (BIFT applicable).")
        print(f"PASS  Theoreme Z2a : inversion locale autour de q^* dans C^0.")
    else:
        print("PARTIAL  Certaines configurations donnent H <= 0.")
        for r in results:
            if not r["min_positive"]:
                print(f"  eps={r['eps']}, lambda={r['lambda']:.0e}, "
                      f"min H = {r['H_min']:.4e}")

    return 0 if all_positive else 1


if __name__ == "__main__":
    sys.exit(main())
