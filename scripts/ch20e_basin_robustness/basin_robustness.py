"""
basin_robustness.py — Chantier 5 du programme Class B Extensions

Teste la robustesse combinatoire du bassin d'attraction de mu* = 15 sous
ajout de contenu BSM. Itère la self-consistency map
    mu_{k+1} = sum( p : gamma_p(mu_k) > 1/2 )
avec le sous-ensemble des premiers actifs pondéré par le contenu BSM.

Méthode :
  - Modèle de référence SM : active primes = {3,5,7}, mu* = 15.
  - Pour chaque scénario BSM, on ajoute Delta_mu(BSM) à l'échelle initiale
    mu_0 (nombre de Weyl-fermions supplémentaires / Higgs extras, normalisé
    par la règle N_Weyl/gen = 4 N_c + 3 = 15, donc 1 Weyl = 1 unité mu).
  - On simule l'itération avec la map PT exacte (gamma_p via q = 1-2/mu,
    sin^2 holonomie, seuil s = 1/2) et on teste si mu* reste 15
    (Classe B compatible) ou sort du bassin [12,17] (Classe C).

Auteur : Class-B extension program (2026-04-22)
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Callable

# ---------------------------------------------------------------------------
# PT constants
# ---------------------------------------------------------------------------

MU_STAR_SM = 15
BASIN_LO, BASIN_HI = 12, 17
S_THRESHOLD = Fraction(1, 2)  # gamma_p > s = 1/2 threshold for activity

# Primes candidates pour activation (on prend p=3..251)
def sieve_primes(n_max: int) -> list[int]:
    """Simple Eratosthenes up to n_max."""
    sieve = [True] * (n_max + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n_max**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, n_max + 1, i):
                sieve[j] = False
    return [i for i, v in enumerate(sieve) if v]


PRIMES = sieve_primes(251)[1:]  # skip 2 (kinematic); start at 3


# ---------------------------------------------------------------------------
# PT holonomy : exact rationals
# ---------------------------------------------------------------------------

def q_stat(mu: int) -> Fraction:
    """q_+ = 1 - 2/mu (discrete)."""
    return Fraction(1) - Fraction(2, mu)


def delta_p(p: int, mu: int) -> Fraction:
    """delta_p = (1 - q^p)/p."""
    q = q_stat(mu)
    return (Fraction(1) - q ** p) / p


def sin2_theta(p: int, mu: int) -> Fraction:
    """sin^2 theta_p = delta_p (2 - delta_p)."""
    d = delta_p(p, mu)
    return d * (2 - d)


def gamma_p(p: int, mu: int) -> float:
    """
    gamma_p(mu) = -d ln(sin^2 theta_p) / d ln mu

    Closed form:
        gamma_p = [4 p q^{p-1} (1 - delta_p)] / [mu (1 - q^p)(2 - delta_p)]
    Evaluated in float (no issue since sin^2 in (0,1)).
    """
    q = float(q_stat(mu))
    if q <= 0 or q >= 1:
        return 0.0
    d = float(delta_p(p, mu))
    num = 4.0 * p * (q ** (p - 1)) * (1.0 - d)
    den = mu * (1.0 - q ** p) * (2.0 - d)
    if den <= 0:
        return 0.0
    return num / den


# ---------------------------------------------------------------------------
# Core iteration map
# ---------------------------------------------------------------------------

def active_primes(mu: int, p_max: int = 251) -> list[int]:
    """Return list of primes p in [3, p_max] with gamma_p(mu) > 1/2."""
    out = []
    for p in PRIMES:
        if p > p_max:
            break
        if gamma_p(p, mu) > 0.5:
            out.append(p)
    return out


def self_consistency_step(mu: int, p_max: int = 251) -> int:
    """One step of the unrestricted self-consistency map."""
    actives = active_primes(mu, p_max)
    return sum(actives)


def iterate_fixed_point(mu_0: int, p_max: int = 251, max_iter: int = 40,
                        verbose: bool = False) -> tuple[int, str, list[int]]:
    """
    Iterate the self-consistency map from mu_0.

    Returns
    -------
    (mu_final, regime, trajectory)
        regime in {"fixed", "collapse", "divergent"}
    """
    traj = [mu_0]
    mu = mu_0
    for k in range(max_iter):
        mu_new = self_consistency_step(mu, p_max)
        traj.append(mu_new)
        if verbose:
            print(f"  k={k}: mu={mu} -> mu_new={mu_new} "
                  f"actives={active_primes(mu, p_max)}")
        if mu_new == mu:
            return mu, "fixed", traj
        if mu_new == 0:
            return 0, "collapse", traj
        if mu_new > 6 * MU_STAR_SM:  # arbitrary cutoff for divergence
            return mu_new, "divergent", traj
        mu = mu_new
    return mu, "divergent", traj


# ---------------------------------------------------------------------------
# BSM scenario registry
# ---------------------------------------------------------------------------

@dataclass
class Scenario:
    name: str
    delta_weyl: int          # additional Weyl fermions in the counting
    delta_scalar: int        # additional real scalars (counted 1/2-unit each)
    pt_class: str            # "B", "C", "A-excluded" (for reference)
    rationale: str

    def perturbed_mu0(self, mu_sm: int = MU_STAR_SM) -> int:
        """Apply small perturbation to mu_0 before iteration."""
        # Convention : 1 Weyl fermion = +1 unit ; 1 real scalar = +1/2 unit
        delta = self.delta_weyl + self.delta_scalar // 2
        return mu_sm + delta


# Scénarios testés (cohérents avec ch20d)
SCENARIOS: list[Scenario] = [
    Scenario(
        name="SM only (reference)",
        delta_weyl=0, delta_scalar=0,
        pt_class="--",
        rationale="N_Weyl = 4 N_c + 3 = 15 per generation (PT baseline).",
    ),
    # Class B candidates
    Scenario(
        name="Axion / ALP (p=2 channel)",
        delta_weyl=0, delta_scalar=1,
        pt_class="B",
        rationale="1 real pseudo-scalar in p=2 info/anti-info; +1/2 unit.",
    ),
    Scenario(
        name="DM scalar singlet (Higgs portal)",
        delta_weyl=0, delta_scalar=1,
        pt_class="B",
        rationale="1 real scalar S (neutral, SM-singlet).",
    ),
    Scenario(
        name="Sterile neutrino eV (1 nu_R)",
        delta_weyl=2, delta_scalar=0,  # Majorana = 2 Weyl; Dirac = 2 Weyl
        pt_class="B",
        rationale="One sterile Majorana (eq 2 Weyl) per generation? Test 1 flavor.",
    ),
    Scenario(
        name="Sterile neutrinos 3x (full Majorana see-saw)",
        delta_weyl=6, delta_scalar=0,
        pt_class="B/C",
        rationale="3 heavy Majorana partners = 6 Weyl — potentially leaves basin.",
    ),
    Scenario(
        name="Light higgsino isolated",
        delta_weyl=4, delta_scalar=2,
        pt_class="B",
        rationale="Higgsino = 2 Dirac fermions = 4 Weyl + 2 complex scalars (=4 real).",
    ),
    Scenario(
        name="2HDM (second Higgs doublet)",
        delta_weyl=0, delta_scalar=4,
        pt_class="B",
        rationale="Second doublet = 4 real scalars (H+, H-, A0, H0).",
    ),
    Scenario(
        name="3HDM",
        delta_weyl=0, delta_scalar=8,
        pt_class="B/C",
        rationale="Two extra doublets = 8 real scalars.",
    ),
    # Class C / boundary
    Scenario(
        name="Natural SUSY partial (higgsino + stop + gluino)",
        delta_weyl=8, delta_scalar=12,
        pt_class="C",
        rationale="Higgsino (4W) + stop (12 real scalars LR+colour) + gluino (8W).",
    ),
    Scenario(
        name="Complete MSSM (1 generation)",
        delta_weyl=15, delta_scalar=30,
        pt_class="A",  # in the sense of excluded
        rationale="Full MSSM: N_Weyl doubles from 15 to 30.",
    ),
    # Arithmetic perturbations for robustness test
    Scenario(
        name="Tiny perturbation +1",
        delta_weyl=1, delta_scalar=0,
        pt_class="--",
        rationale="Robustness check: mu_0 = 16.",
    ),
    Scenario(
        name="Small perturbation +2",
        delta_weyl=2, delta_scalar=0,
        pt_class="--",
        rationale="mu_0 = 17 = upper basin boundary.",
    ),
    Scenario(
        name="Basin exit +3",
        delta_weyl=3, delta_scalar=0,
        pt_class="--",
        rationale="mu_0 = 18 = just outside basin [12,17].",
    ),
]


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def in_basin(mu: int) -> bool:
    return BASIN_LO <= mu <= BASIN_HI


def classify_result(mu_final: int, regime: str) -> str:
    if regime == "fixed" and mu_final == MU_STAR_SM:
        return "STABLE (mu*=15, PT intact)"
    if regime == "fixed" and in_basin(mu_final):
        return f"WITHIN BASIN (mu*={mu_final} in [12,17])"
    if regime == "collapse":
        return "COLLAPSE (mu->0)"
    if regime == "divergent":
        return "DIVERGENT (no stable phase)"
    return f"OTHER (mu={mu_final}, {regime})"


def run_scenario(s: Scenario, p_max: int = 251) -> dict:
    mu_0 = s.perturbed_mu0()
    mu_final, regime, traj = iterate_fixed_point(mu_0, p_max)
    return dict(
        name=s.name,
        pt_class_expected=s.pt_class,
        mu_0=mu_0,
        mu_final=mu_final,
        regime=regime,
        trajectory=traj,
        basin_respected=in_basin(mu_final) if regime == "fixed" else False,
        verdict=classify_result(mu_final, regime),
        rationale=s.rationale,
    )


def print_table(results: list[dict]) -> None:
    head = f"{'Scenario':<45} {'mu_0':>5} {'mu_final':>8} {'regime':<11} {'verdict'}"
    print(head)
    print("-" * len(head))
    for r in results:
        print(f"{r['name']:<45} {r['mu_0']:>5} {r['mu_final']:>8} "
              f"{r['regime']:<11} {r['verdict']}")


def main():
    print(f"=== Basin robustness under BSM perturbations (p_max=251) ===")
    print(f"SM baseline : mu* = {MU_STAR_SM}, basin = [{BASIN_LO}, {BASIN_HI}]")
    print(f"Threshold   : gamma_p > s = 1/2")
    print()

    results = [run_scenario(s) for s in SCENARIOS]
    print_table(results)

    # Breakdown
    stable = [r for r in results if r["mu_final"] == MU_STAR_SM and r["regime"] == "fixed"]
    within = [r for r in results if r["regime"] == "fixed" and in_basin(r["mu_final"])]
    divergent = [r for r in results if r["regime"] == "divergent"]

    print()
    print(f"Stable at mu*=15 : {len(stable)} / {len(results)}")
    print(f"Within basin     : {len(within)} / {len(results)}")
    print(f"Divergent        : {len(divergent)} / {len(results)}")

    # Breakage threshold search : number of extra Weyl fermions for basin exit
    print()
    print("=== Breakage threshold (extra Weyl fermions from SM baseline) ===")
    for d in range(0, 20):
        mu0 = MU_STAR_SM + d
        muf, reg, _ = iterate_fixed_point(mu0)
        tag = "OK" if reg == "fixed" and muf == MU_STAR_SM else (
            "WITHIN" if reg == "fixed" and in_basin(muf) else reg.upper()
        )
        marker = ">>>" if d == 3 else "   "
        print(f"  {marker} mu_0 = 15 + {d:2d} = {mu0:3d}  ->  "
              f"mu_final = {muf:4d}  [{reg:10s}]  {tag}")


if __name__ == "__main__":
    main()
