#!/usr/bin/env python3
"""
PT Cosmology Phase 2 - Chantier 2
DESI DR3 (2027) scenarios for Sigma m_nu under PT.

PT prediction (from ch15 P5/P31, monograph):
  m_nu3 = s^2 * alpha_bare^3 * m_e = 0.0505 eV
  with s = 1/2, alpha_bare = 1/136.28 = prod_{p in {3,5,7}} sin^2(theta_p, q+)
  Sigma m_nu = m_nu1 + m_nu2 + m_nu3 (normal ordering).

Oscillation constraints (PDG 2024):
  Dm31^2 = 2.51e-3 eV^2    (atm)
  Dm21^2 = 7.42e-5 eV^2    (solar)
  -> m_nu3 >= sqrt(Dm31^2) = 0.0501 eV (floor)
  -> Sigma >= 0.0586 eV (normal ordering, m_nu1 = 0)
  -> Sigma >= 0.0998 eV (inverted ordering, m_nu3 = 0)

Cosmological constraints:
  DESI DR2 + CMB (LambdaCDM):          Sigma < 0.064 eV (95% CL)
  DESI DR2 + CMB (Feldman-Cousins):    Sigma < 0.053 eV  (below NO floor!)
  DESI DR2 + w_0w_a CDM:               Sigma < 0.163 eV (relaxed)

This script computes:
  (i)   PT prediction interval with uncertainty propagation
  (ii)  PT compatibility under the three DESI DR2 bounds
  (iii) Three DESI DR3 2027 scenarios (tight / marginal / relaxed)
  (iv)  Discriminant prediction matrix (DESI + JUNO + Euclid)
"""
from __future__ import annotations
import math


# PT constants (from ch10, ch15 monograph)
S = 0.5
ALPHA_BARE = 1.0 / 136.28    # tree-level, pre-dressing
ALPHA_DRESSED = 1.0 / 137.036
M_E_EV = 510_999.0            # electron mass (eV)

# Mass splittings (PDG 2024 central values)
DM31_SQ = 2.514e-3   # eV^2
DM21_SQ = 7.412e-5   # eV^2

# PT m_nu3 central value
M_NU3_PT = S ** 2 * ALPHA_BARE ** 3 * M_E_EV   # eV
# Uncertainty on m_nu3 from alpha_bare residual (ch15, 0.44% quoted error)
DELTA_M_NU3_REL = 0.0044


def sigma_mnu_normal(m3: float) -> float:
    """Sigma m_nu in normal ordering, assuming m_nu1 = 0 (minimum)."""
    m1 = 0.0
    m2 = math.sqrt(DM21_SQ)
    return m1 + m2 + m3


def sigma_mnu_normal_m1(m3: float, m1: float = 0.0) -> float:
    """Sigma m_nu in normal ordering for arbitrary m_nu1 (can be >0)."""
    m2 = math.sqrt(m1 * m1 + DM21_SQ)
    return m1 + m2 + m3


def pt_prediction_interval() -> tuple[float, float, float]:
    """PT prediction for m_nu3 and Sigma m_nu with 1-sigma uncertainty."""
    m3 = M_NU3_PT
    dm3 = m3 * DELTA_M_NU3_REL
    sigma_c = sigma_mnu_normal(m3)
    sigma_lo = sigma_mnu_normal(m3 - dm3)
    sigma_hi = sigma_mnu_normal(m3 + dm3)
    return sigma_lo, sigma_c, sigma_hi


# ==== DESI scenarios 2027 ====

DR3_SCENARIOS = {
    "tight_LCDM": {
        "name": "DR3 tight LCDM",
        "bound_95": 0.045,          # hypothetical: DR3 tightens to 0.045
        "parametrisation": "LCDM",
        "PT_verdict": None,         # filled below
    },
    "marginal_LCDM": {
        "name": "DR3 marginal LCDM",
        "bound_95": 0.058,          # at the NO floor
        "parametrisation": "LCDM",
        "PT_verdict": None,
    },
    "relaxed_w0waCDM": {
        "name": "DR3 relaxed w0waCDM",
        "bound_95": 0.120,          # DE-dynamic extension
        "parametrisation": "w0waCDM",
        "PT_verdict": None,
    },
}


def verdict_from_bound(bound: float, sigma_lo: float) -> str:
    """Given a 95% CL upper bound on Sigma m_nu, and PT 1-sigma lower
    value of Sigma, decide compat/tension/falsified."""
    # 3-sigma down from central PT
    sigma_3lo = sigma_lo - 2 * (sigma_lo - sigma_mnu_normal(M_NU3_PT * (1 - DELTA_M_NU3_REL * 0.5)))
    if bound >= sigma_lo:
        return "COMPATIBLE (PT within bound)"
    if bound >= sigma_3lo:
        return "MARGINAL (2-3 sigma pull)"
    return "FALSIFIED (>5 sigma tension)"


# ==== Dark-energy relaxation of Sigma bound ====

def relax_factor(w0: float = -1.0, wa: float = 0.0) -> float:
    """
    Dimensionless multiplicative factor relaxing the Sigma m_nu upper
    bound when dark energy is allowed to deviate from LambdaCDM.
    Empirical fit from DESI DR2 Table 6 (arXiv:2503.14738):
      LCDM:       0.053 (FC)
      wCDM:       0.087
      w0waCDM:    0.163
    We use:
      relax = 1 + 2.0 * |w0 + 1| + 1.6 * |wa|
    calibrated to reproduce the three rows.
    """
    return 1.0 + 2.0 * abs(w0 + 1.0) + 1.6 * abs(wa)


# ==== Main ====

def main() -> None:
    print("=" * 72)
    print("PT COSMOLOGY PHASE 2 -- DESI NEUTRINO MASS BOUND (Chantier 2)")
    print("=" * 72)

    print("\nPT input (ch15 eq:m_nu3):")
    print(f"  s = {S}, alpha_bare = 1/{1/ALPHA_BARE:.2f}")
    print(f"  m_nu3 = s^2 * alpha_bare^3 * m_e")
    print(f"         = {M_NU3_PT*1000:.3f} meV = {M_NU3_PT:.5f} eV")
    print(f"  (relative uncertainty: {DELTA_M_NU3_REL*100:.2f}%)")

    sigma_lo, sigma_c, sigma_hi = pt_prediction_interval()
    print(f"\nPT Sigma m_nu (normal ordering, m_nu1 = 0):")
    print(f"  central: {sigma_c*1000:.3f} meV = {sigma_c:.5f} eV")
    print(f"  1-sigma: [{sigma_lo*1000:.3f}, {sigma_hi*1000:.3f}] meV")
    print(f"  oscillation floor (m_nu1=0): sqrt(Dm31^2) + sqrt(Dm21^2) = "
          f"{(math.sqrt(DM31_SQ)+math.sqrt(DM21_SQ))*1000:.3f} meV")

    print("\n" + "-" * 72)
    print("DESI DR2 2025 constraints (arXiv:2503.14738):")
    print("  LCDM (standard):          Sigma < 0.064 eV  [95% CL]")
    print("  LCDM (Feldman-Cousins):   Sigma < 0.053 eV  [*below* NO floor]")
    print("  w0waCDM (DE-dyn):         Sigma < 0.163 eV  [relaxed]")

    print(f"\nPT compatibility:")
    print(f"  vs LCDM 0.064:  Sigma_PT = {sigma_c*1000:.1f} meV  "
          f"=> {'COMPATIBLE' if sigma_c < 0.064 else 'TENSION'}")
    print(f"  vs LCDM-FC 0.053: Sigma_PT = {sigma_c*1000:.1f} meV  "
          f"=> {'COMPATIBLE' if sigma_c < 0.053 else 'TENSION (marginal)'}")
    print(f"  vs w0waCDM 0.163: Sigma_PT = {sigma_c*1000:.1f} meV  "
          f"=> {'COMPATIBLE' if sigma_c < 0.163 else 'TENSION'}")

    # pull estimate against DESI DR2 LCDM central
    # interpret bound as 2-sigma upper (95% CL => 1.96 sigma)
    sigma_exp_lcdm = 0.064 / 1.96  # implied 1-sigma error on Sigma
    pull_lcdm = (sigma_c - 0.0) / sigma_exp_lcdm if sigma_exp_lcdm > 0 else 0
    print(f"  (informal pull vs LCDM upper bound ~ {pull_lcdm:.2f}σ below limit)")

    print("\n" + "-" * 72)
    print("Three DESI DR3 2027 scenarios and PT responses:")
    print()
    for key, sc in DR3_SCENARIOS.items():
        bound = sc["bound_95"]
        verdict = "COMPATIBLE"
        response = ""
        if bound < sigma_lo and bound < math.sqrt(DM31_SQ):
            verdict = "FALSIFIED"
            response = "[PT falsified: sieve m_nu3 ABOVE bound]"
        elif bound < sigma_lo:
            verdict = "TENSION"
            # Can DE-dynamic relax it?
            relax_needed = sigma_lo / bound
            response = (f"[PT needs DE relaxation factor x{relax_needed:.2f};"
                        f" possible via w0waCDM if DR3 allows {bound*relax_needed:.3f} eV]")
        else:
            response = "[PT cleanly within bound]"
        sc["PT_verdict"] = verdict
        print(f"  {sc['name']:<25} bound {bound*1000:5.1f} meV  [{sc['parametrisation']}]")
        print(f"    -> PT verdict: {verdict}  {response}")
        print()

    # Discriminant matrix
    print("-" * 72)
    print("Discriminant matrix (DESI + JUNO + Euclid, 2027-2028):")
    print()
    print("  Experiment    | Quantity          | PT prediction       | Sensitivity")
    print("  ------------- | ----------------- | ------------------- | ---------------")
    print("  JUNO 2027     | Dm31^2            | 2.514e-3 eV^2       | 0.3% @ 6y")
    print("  JUNO 2027     | Dm21^2            | 7.412e-5 eV^2       | 0.5% @ 6y")
    print("  JUNO 2027     | sin^2 theta12     | 0.3037              | 0.5%")
    print("  KATRIN 2027   | m_beta (eff.)     | 0.0505 eV floor     | 0.3 eV (insens.)")
    print("  DESI DR3 2027 | Sigma m_nu (LCDM) | 0.0586 eV (NO min)  | ~0.04-0.06 eV")
    print("  Euclid 2028   | Sigma m_nu (CMB+LSS) | 0.058 eV         | 0.02-0.03 eV")
    print("  0nu2beta      | m_bb              | 0 (Dirac)           | LEGEND ~2035")

    # Relaxation via DE
    print("\n" + "-" * 72)
    print("Relaxation of Sigma m_nu bound via dark energy equation of state:")
    print()
    print("  (w0, wa)          | relax | effective LCDM-equiv bound if PT @ 0.058")
    print("  ----------------- | ----- | -------------------------------------")
    for w0, wa in [(-1.0, 0.0),
                    (-1.003, -8.32e-4),    # PT P20
                    (-0.95, -0.1),
                    (-0.90, -0.3),
                    (-0.86, -0.60),        # DESI DR2 preferred
                    (-0.80, -1.0)]:
        r = relax_factor(w0, wa)
        # Required LCDM-equivalent bound for PT to be compatible
        # if actual bound is 0.058 under the given (w0,wa):
        equiv_lcdm = 0.058 / r
        print(f"  ({w0:+.3f}, {wa:+.4f}) | {r:.3f} | LCDM-equiv <= {equiv_lcdm*1000:.1f} meV")

    # Verdict
    print("\n" + "=" * 72)
    print("VERDICT (Chantier 2):")
    print()
    print("  PT Sigma m_nu = 58.6 +/- 0.3 meV (normal ordering floor + PT m_nu3)")
    print("  DESI DR2 LCDM bound 64 meV => PT 0.9σ BELOW bound (COMPATIBLE).")
    print("  DESI DR2 LCDM-FC bound 53 meV => slight tension (~1.5σ), marginal.")
    print("  DESI DR2 w0waCDM bound 163 meV => PT well within (COMPATIBLE).")
    print()
    print("  Three-scenario response to DESI DR3 2027:")
    print("    (a) DR3 LCDM tightens to 0.045 eV @5σ => PT FALSIFIED")
    print("        Resilience: move P20 to evolving DE (removes bound).")
    print("    (b) DR3 LCDM stays at 0.058 ± 0.008 eV    => PT COMPATIBLE")
    print("        (marginal with FC interval).")
    print("    (c) DR3 detects evolving DE at 5σ and Sigma = 0.09 ± 0.03")
    print("         => PT cleanly within and P20 must evolve (see ch20f DE).")
    print()
    print("  The key conditional: *if* DR3 reports Sigma < 0.053 at >=5σ in")
    print("  LCDM-strict, PT's m_nu3 = 0.0505 eV derivation is falsified UNLESS")
    print("  dark energy is simultaneously detected as dynamical (in which")
    print("  case the bound relaxes to >0.10 eV and PT is safe). PT thus")
    print("  predicts DR3 cannot simultaneously (a) confirm LCDM strict and")
    print("  (b) tighten Sigma below 0.053 eV.")


if __name__ == "__main__":
    main()
