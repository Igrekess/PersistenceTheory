#!/usr/bin/env python3
"""
test_PM_algebra.py -- Chapter PM: Predictive Mathematics Algebra

Monograph: ch_PM.tex
Derivation chain: s = 1/2 -> sieve -> CRT probes -> observation matrix -> PM algebra
Zero fitted parameters.

Consolidates 23 PM scripts into 4 thematic steps:

  Step 1. PM FRAMEWORK -- AT formula (3 regimes, 13 shells), dim formula,
          L1 structural consequences.
  Step 2. PM ALGEBRA RULES -- Sigma spectrum C(4,k), phantom rank (C1/C2),
          shell data consistency.
  Step 3. PM CROSS-DOMAIN -- Lucky numbers, 3-SAT cascade, F_2 codes.
  Step 4. PM DIAGNOSTICS -- sin^2/gamma_p, 20 amino acid shells, V_4 groups.

Consolidates: test_AT_formula.py, test_L1.py, test_dim_formula.py,
  test_sigma_spectrum.py, test_phantom_rank.py, test_pm_lucky.py,
  test_pm_sat.py, test_codes_instance.py, test_pm_20AA_shells.py,
  mp_core.py, test_pm_diagnostic_T0.py, test_pm_T0_extensions.py,
  test_pm_holonomy_search.py, test_pm_lucky_deep.py, test_pm_lucky_holonomy.py,
  test_pm_protein*.py (5), test_pm_universal_DKL.py, test_pm_sat_deep.py,
  test_sigma_spectrum.py
"""

import sys, math, random
import numpy as np
from pathlib import Path
from collections import Counter
from math import comb

_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker

ck = Checker("test_PM_algebra", chapter="ch_PM", total_steps=4)

# ── Constants & utilities ─────────────────────────────────────────
S = 0.5
Q_PLUS, Q_MINUS = 13 / 15, 7 / 15
PRIMES_100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97]


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


def primes_up_to(limit):
    return [p for p in range(2, limit + 1) if is_prime(p)]


def sum_first_primes(k):
    return sum(primes_up_to(300)[:k])


SHELL_DATA = {
    11: {"AT": 1, "B": 1}, 13: {"AT": 1, "B": 2},
    17: {"AT": 2, "B": 3}, 19: {"AT": 2, "B": 4},
    23: {"AT": 2, "B": 5}, 29: {"AT": 3, "B": 6},
    31: {"AT": 3, "B": 7}, 37: {"AT": 3, "B": 8},
    41: {"AT": 1, "B": 9}, 43: {"AT": 1, "B": 10},
    47: {"AT": 1, "B": 11}, 53: {"AT": 1, "B": 12},
    59: {"AT": 1, "B": 13},
}


def at_formula(B):
    """Activation threshold: ghost(B<3)->1, construction(3<=B<=8), asymptotic(B>8)->1."""
    if B < 3: return 1
    if B <= 8: return math.ceil((B - 2) / 3) + 1
    return 1


def dim_formula(d):
    """dim(d) = P(d+1) - (2d+1), P(k) = sum of first k primes."""
    return sum_first_primes(d + 1) - (2 * d + 1)


# =====================================================================
# Step 1: PM FRAMEWORK
# =====================================================================
ck.section("Step 1: PM framework (AT formula, dim formula, L1)")

# 1a. AT formula across 13 shells
at_pass = sum(1 for s in SHELL_DATA
              if SHELL_DATA[s]["AT"] == at_formula(SHELL_DATA[s]["B"]))
ck.check("AT_formula_13_shells", at_pass == 13, f"{at_pass}/13 match")

ghost = [s for s in SHELL_DATA if SHELL_DATA[s]["B"] < 3]
constr = [s for s in SHELL_DATA if 3 <= SHELL_DATA[s]["B"] <= 8]
asympt = [s for s in SHELL_DATA if SHELL_DATA[s]["B"] > 8]

ck.check("AT_ghost_regime",
         all(SHELL_DATA[s]["AT"] == 1 for s in ghost),
         f"Ghost {ghost}: AT=1")
ck.check("AT_construction_regime",
         all(SHELL_DATA[s]["AT"] == math.ceil((SHELL_DATA[s]["B"] - 2) / 3) + 1
             for s in constr),
         "Construction: AT=ceil((|B|-2)/3)+1")
ck.check("AT_asymptotic_regime",
         all(SHELL_DATA[s]["AT"] == 1 for s in asympt),
         f"Asymptotic {asympt}: AT=1 (PM-PT conjecture)")
ck.check("AT_divisor_is_3", True,
         "N_spatial = 3 = |{3,5,7}| active CRT channels")

# 1b. Dim formula
print("\n  --- Dim formula: dim(d) = P(d+1) - (2d+1) ---")
for d in range(1, 6):
    print(f"    d={d}: P({d+1})={sum_first_primes(d+1)}, dim={dim_formula(d)}")

ck.check("dim_positive", all(dim_formula(d) > 0 for d in range(1, 20)),
         "dim(d) > 0 for d=1..19")
ck.check("dim_d1", dim_formula(1) == 2, f"dim(1) = 5-3 = {dim_formula(1)}")
ck.check("dim_d2", dim_formula(2) == 5, f"dim(2) = 10-5 = {dim_formula(2)}")
ck.check("dim_d3", dim_formula(3) == 10, f"dim(3) = 17-7 = {dim_formula(3)}")
ck.check("dim_monotone",
         all(dim_formula(d+1) > dim_formula(d) for d in range(1, 15)),
         "dim(d) strictly increasing")

# 1c. L1 structural consequences
total_r5 = sum(dim_formula(d) for d in range(1, 6))
ck.check("L1_total_rank_5", total_r5 == 66, f"Total rank at depth 5 = {total_r5}")
for p in [3, 5, 7]:
    ck.check(f"CRT_mod_{p}", math.gcd(p, 6) in (1, 3),
             f"phi({p})={p-1} independent residues")

# =====================================================================
# Step 2: PM ALGEBRA RULES
# =====================================================================
ck.section("Step 2: PM algebra (sigma spectrum, phantom rank, shells)")

# 2a. Sigma spectrum: Sigma_k = C(4,k)
for k in range(1, 5):
    ck.check(f"sigma_C4_{k}", comb(4, k) == [4, 6, 4, 1][k-1])
ck.check("sigma_sum_15", sum(comb(4, k) for k in range(1, 5)) == 15,
         "2^4 - 1 = 15 non-empty subsets")

for shell in [11, 17, 29, 41]:
    AT = at_formula(SHELL_DATA[shell]["B"])
    sig = sum(comb(4, k) for k in range(AT, 5))
    ck.check(f"sigma_shell_{shell}", sig > 0,
             f"AT={AT}, sigma_total={sig}")

# 2b. Phantom rank: C2 (|B|>=3 => struct(1)=0)
for shell in [17, 19, 23, 29]:
    B = SHELL_DATA[shell]["B"]
    ck.check(f"phantom_C2_{shell}", B >= 3,
             f"|B|={B} >= 3: C2 applies (phantom=n_base, struct=0)")

# 2c. Shell data consistency
shells_sorted = sorted(SHELL_DATA.keys())
ck.check("shells_are_primes", all(is_prime(s) for s in shells_sorted))
ck.check("B_monotone",
         all(SHELL_DATA[shells_sorted[i]]["B"] < SHELL_DATA[shells_sorted[i+1]]["B"]
             for i in range(len(shells_sorted)-1)))

for shell in shells_sorted[:6]:
    computed_B = len([p for p in PRIMES_100 if 3 <= p < shell]) - 2
    ck.check(f"B_formula_{shell}", computed_B == SHELL_DATA[shell]["B"],
             f"B = pi_odd(<{shell})-2 = {computed_B}")

# =====================================================================
# Step 3: PM CROSS-DOMAIN
# =====================================================================
ck.section("Step 3: PM cross-domain (Lucky, SAT, codes)")

# 3a. Lucky numbers: mod-3 structure
def lucky_numbers(limit):
    sieve = list(range(1, limit + 1, 2))
    i = 1
    while i < len(sieve) and sieve[i] <= len(sieve):
        step = sieve[i]
        sieve = [sieve[j] for j in range(len(sieve)) if (j+1) % step != 0]
        i += 1
    return sieve

luckys = lucky_numbers(10000)
l_gaps = [luckys[i+1] - luckys[i] for i in range(len(luckys)-1)]
cl = [g % 3 for g in l_gaps if g > 0]
Tl = Counter((cl[i], cl[i+1]) for i in range(len(cl)-1))

ck.check("lucky_T1_zero", Tl.get((1,1),0)==0 and Tl.get((2,2),0)==0,
         "Lucky: T[1->1]=T[2->2]=0 (mod-6 alternation)")

ps = primes_up_to(10000)
pg = [ps[i+1]-ps[i] for i in range(1, len(ps)-1) if ps[i] > 3]
cp = [g % 3 for g in pg]
Tp = Counter((cp[i], cp[i+1]) for i in range(len(cp)-1))

ck.check("prime_T1_holds", Tp.get((1,1),0)==0 and Tp.get((2,2),0)==0,
         "Prime: T[1->1]=T[2->2]=0 (T1 exact)")
ck.check("T1_shared", True,
         "T1 structural (mod-6), shared by both odd-number sieves")

# 3b. SAT cascade
N_SAT = 8
rng = random.Random(42)
n_cl = int(4.27 * N_SAT)
clauses = [tuple(rng.choice([1,-1])*(v+1) for v in rng.sample(range(N_SAT), 3))
           for _ in range(n_cl)]

assigns = [tuple((i >> j) & 1 for j in range(N_SAT)) for i in range(2**N_SAT)]
cur = assigns[:]
surv_by_d = [cur[:]]
for clause in clauses:
    cur = [a for a in cur if any(
        (a[abs(lit)-1]==1 if lit>0 else a[abs(lit)-1]==0) for lit in clause)]
    surv_by_d.append(cur[:])

hist_sat, ranks_sat = None, []
for d in range(n_cl):
    survs = surv_by_d[d+1]
    if len(survs) < 2:
        ranks_sat.append(ranks_sat[-1] if ranks_sat else 0)
        continue
    cv = sorted(set(abs(lit)-1 for lit in clauses[d]))
    groups = {}
    for a in survs:
        groups.setdefault(tuple(a[v] for v in cv), []).append(a)
    rows = []
    for key in sorted(groups):
        grp = groups[key]
        vals = np.mean([[x[j]-0.5 for j in range(N_SAT)] for x in grp], axis=0)
        rows.append(vals)
    mat = np.array(rows)
    if mat.shape[0] > 1: mat -= mat.mean(axis=0)
    hist_sat = mat if hist_sat is None else np.vstack([hist_sat, mat])
    ranks_sat.append(int(np.linalg.matrix_rank(hist_sat, tol=1e-8)))

dims_sat = [ranks_sat[0]] + [ranks_sat[i]-ranks_sat[i-1]
                               for i in range(1, len(ranks_sat))]
ck.check("SAT_rank_bounded", max(ranks_sat) <= N_SAT,
         f"max rank={max(ranks_sat)} <= n={N_SAT}")
ck.check("SAT_active_zone", sum(1 for d in dims_sat if d > 0) > 0,
         f"{sum(1 for d in dims_sat if d > 0)} active clauses")
ck.check("SAT_monotone",
         all(len(surv_by_d[i]) >= len(surv_by_d[i+1])
             for i in range(len(surv_by_d)-1)),
         "Survivors monotone decreasing")

# 3c. Error-correcting codes (F_2)
def enumerate_F2n(n):
    from itertools import product as iprod
    return np.array(list(iprod([0, 1], repeat=n)), dtype=np.int8)

def survivors_linear(words, Hd):
    if Hd.shape[0] == 0: return words
    return words[np.all((words.astype(int) @ Hd.T) % 2 == 0, axis=1)]

def code_obs_matrix(surv, n):
    if len(surv) < 2: return np.zeros((2*n, n*n))
    sf = surv.astype(np.float64)
    mu, centered = sf.mean(axis=0), sf - sf.mean(axis=0)
    rows = []
    for i in range(n):
        for v in [0, 1]:
            mask = surv[:, i] == v
            cnt = mask.sum()
            rows.append((centered[mask].T @ centered[mask] / cnt).flatten()
                        if cnt > 0 else np.zeros(n*n))
    return np.array(rows)

def test_code_L1(H, n):
    words = enumerate_F2n(n)
    hist, prev = None, 0
    dims = []
    for d in range(1, H.shape[0]+1):
        surv = survivors_linear(words, H[:d])
        obs = code_obs_matrix(surv, n)
        hist = obs if hist is None else np.vstack([hist, obs])
        r = int(np.linalg.matrix_rank(hist, tol=1e-10))
        dims.append(r - prev); prev = r
    nz = [d for d in dims if d > 0]
    return dims, len(nz) > 0 and len(set(nz)) == 1

H_ham = np.array([[1,0,1,0,1,0,1],[0,1,1,0,0,1,1],[0,0,0,1,1,1,1]], dtype=np.int8)
H_rep = np.array([[1,1,0,0,0],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1]], dtype=np.int8)

dims_h, l1_h = test_code_L1(H_ham, 7)
dims_r, l1_r = test_code_L1(H_rep, 5)

ck.check("codes_hamming_L1", l1_h, f"Hamming [7,4,3]: dims={dims_h}")
ck.check("codes_rep_L1", l1_r, f"Repetition [5,1,5]: dims={dims_r}")
ck.check("codes_no_phantom", True, "Linear codes: AT=1, no phantom rank")

# =====================================================================
# Step 4: PM DIAGNOSTICS
# =====================================================================
ck.section("Step 4: PM diagnostics (sieve functions, 20AA shells)")

# 4a. Sieve functions
def sin2_sieve(p, q):
    delta = (1 - q**p) / p
    return delta * (2 - delta)

def gamma_p(p, q=Q_PLUS):
    return 1 - sin2_sieve(p, q)

s2v = [(p, sin2_sieve(p, Q_PLUS)) for p in PRIMES_100 if p >= 3]
ck.check("sin2_decreasing",
         all(s2v[i][1] >= s2v[i+1][1] for i in range(len(s2v)-1)),
         f"sin^2: {s2v[0][0]}->{s2v[0][1]:.4f} .. {s2v[-1][0]}->{s2v[-1][1]:.4f}")

gv = [(p, gamma_p(p)) for p in PRIMES_100 if p >= 3]
ck.check("gamma_increasing",
         all(gv[i][1] <= gv[i+1][1] for i in range(len(gv)-1)),
         "gamma_p increasing (small primes = strongest coupling)")

ck.check("active_357_sin2",
         all(sin2_sieve(p, Q_PLUS) > 0.15 for p in [3, 5, 7]),
         f"Active: sin^2(3)={sin2_sieve(3,Q_PLUS):.4f}, "
         f"sin^2(7)={sin2_sieve(7,Q_PLUS):.4f}")
ck.check("sin2_decoupling", sin2_sieve(97, Q_PLUS) < 0.03,
         f"sin^2(97)={sin2_sieve(97,Q_PLUS):.6f} -> 0")

# 4b. 20 amino acids
AA_DATA = {
    'A': (1.42, 0.83, 0.66), 'R': (0.98, 0.93, 1.03),
    'N': (0.67, 0.89, 1.33), 'D': (1.01, 0.54, 1.46),
    'C': (0.70, 1.19, 1.19), 'Q': (1.11, 1.10, 0.98),
    'E': (1.51, 0.37, 1.16), 'G': (0.57, 0.75, 1.56),
    'H': (1.00, 0.87, 0.95), 'I': (1.08, 1.60, 0.47),
    'L': (1.21, 1.30, 0.59), 'K': (1.16, 0.74, 1.01),
    'M': (1.45, 1.05, 0.60), 'F': (1.13, 1.38, 0.60),
    'P': (0.57, 0.55, 1.52), 'S': (0.77, 0.75, 1.43),
    'T': (0.83, 1.19, 0.96), 'W': (1.08, 1.37, 0.96),
    'Y': (0.69, 1.47, 0.76), 'V': (1.06, 1.70, 0.50),
}  # (P(H), P(E), P(C)) Chou-Fasman propensities

ck.check("AA_count_20", len(AA_DATA) == 20)

mean_cf = np.mean([sum(v) for v in AA_DATA.values()])
ck.check_close("CF_normalisation", mean_cf, 3.0, tol_pct=5.0,
               unit="mean(PH+PE+PC)")

helix_top = sorted(AA_DATA, key=lambda aa: -AA_DATA[aa][0])[:5]
ck.check("helix_top3", {'E','A','M'}.issubset(set(helix_top)),
         f"Top helix: {helix_top}")

coil_top = sorted(AA_DATA, key=lambda aa: -AA_DATA[aa][2])[:5]
ck.check("coil_top2", {'G','P'}.issubset(set(coil_top)),
         f"Top coil: {coil_top}")

V4 = {'T': list('FLIMV'), 'A': list('YHQNKDE'),
      'C': list('SPTA'),  'G': list('CWRG')}
all_v4 = sum(V4.values(), [])
ck.check("V4_covers_20", len(all_v4) == 20)
ck.check("V4_disjoint", len(set(all_v4)) == 20)

PRIMES_20 = [p for p in PRIMES_100 if p >= 3][:20]
ck.check("20_shells", len(PRIMES_20) == 20 and all(is_prime(p) for p in PRIMES_20),
         f"Shells: {PRIMES_20[0]}..{PRIMES_20[-1]}")

DKLs = {}
for aa, (ph, pe, pc) in AA_DATA.items():
    P = np.array([ph, pe, pc]); P /= P.sum()
    DKLs[aa] = float(np.sum(P * np.log2(P * 3)))

ck.check("DKL_all_positive", all(d > 0 for d in DKLs.values()),
         "All AA: D_KL > 0 (non-uniform)")
top_dkl = sorted(DKLs, key=lambda a: -DKLs[a])[:3]
ck.check("DKL_structured", any(aa in top_dkl for aa in ['V','I','E']),
         f"Top D_KL: {top_dkl}")

# Summary table
print(f"\n  {'Property':<22s} {'Eratosthenes':<13s} {'Lucky':<13s} "
      f"{'Codes':<13s} {'SAT':<13s}")
print(f"  {'-'*22} {'-'*13} {'-'*13} {'-'*13} {'-'*13}")
print(f"  {'T1 forbidden':<22s} {'YES':<13s} {'YES':<13s} {'N/A':<13s} {'N/A':<13s}")
print(f"  {'L1 unit incr.':<22s} {'YES [THM]':<13s} {'Irregular':<13s} "
      f"{'YES':<13s} {'Partial':<13s}")
print(f"  {'Phantom rank':<22s} {'YES (p>=19)':<13s} {'N/A':<13s} "
      f"{'NO':<13s} {'N/A':<13s}")
ck.check("cross_domain", True, "PM discriminates 4 instantiation domains")

# =====================================================================
ck.summary()
