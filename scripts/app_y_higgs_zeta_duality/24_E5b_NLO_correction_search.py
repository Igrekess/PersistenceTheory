"""
E5b : Recherche STRUCTURELLE d'une correction PT NLO/1-loop/dressing
       qui expliquerait l'écart ε = 3.95 × 10⁻⁴ entre

       λ_H(K4) = 2R/17 = 0.124951    et    λ_H(PT) = 1/8 = 0.125

Le script #23 a déjà testé les constantes PT triviales (α_EM, 1/μ*², ...).
Ici on teste des STRUCTURES NLO PT, par analogie avec :
  - α_EM dressing (4 étapes : F(2), spiral, echo, 2-loop, ghost VP) ≈ 0.55 %
  - C_Koide NLO (tree 4/sin²₃ + NLO 1/21 + NNLO 5δ₃²/378) ≈ 0.26 % + 5e-3
  - β_ghost = Σ_{p∈{11,13}} sin²·γ ≈ 0.1039

Stratégie :
  1) bibliothèque PT exacte (sin², γ_p, β_ghost, Δq) en mpmath dps=50
  2) catalogue ~40 combinaisons NLO-plausibles (produits, ratios, Pythagore)
  3) ranking par |ε/candidate - 1| (match si < 1 %)
"""

import mpmath as mp
import json
from pathlib import Path

mp.mp.dps = 50

# ====================================================================
# 1. Constantes PT canoniques (dps=50)
# ====================================================================
mu_star  = mp.mpf(15)
s        = mp.mpf(1)/2
q_plus   = 1 - 2/mu_star                # q_+ = 13/15
q_minus  = mp.exp(-1/mu_star)           # q_-
Delta_q  = q_minus - q_plus             # bifurcation
N_br     = mp.mpf(2)
p2       = mp.mpf(2)

def delta_p(p, q):
    return (1 - q**p) / p

def sin2(p, q):
    d = delta_p(p, q)
    return d * (2 - d)

def gamma_p(p, q):
    d = delta_p(p, q)
    # γ_p = -d ln sin²(θ_p) / d ln μ
    # forme close au point fixe : 4 p q^{p-1} (1-d) / (μ (1-q^p)(2-d))
    return 4*p*q**(p-1)*(1-d) / (mu_star * (1-q**p) * (2-d))

# Calculs au point fixe (q_+)
s2_3 = sin2(3, q_plus)        # 0.21918...
s2_5 = sin2(5, q_plus)        # 0.19405...
s2_7 = sin2(7, q_plus)        # 0.17264...
s2_11 = sin2(11, q_plus)
s2_13 = sin2(13, q_plus)
g3   = gamma_p(3, q_plus)     # 0.80761...
g5   = gamma_p(5, q_plus)     # 0.69632...
g7   = gamma_p(7, q_plus)     # 0.59547...
g11  = gamma_p(11, q_plus)    # 0.42619...
g13  = gamma_p(13, q_plus)    # 0.35578...

alpha_EM_nu = s2_3 * s2_5 * s2_7              # 1/136.28
alpha_EM    = mp.mpf(1)/mp.mpf('137.0359990840')
beta_ghost  = s2_11*g11 + s2_13*g13           # ≈ 0.1039
delta_3     = (1 - q_plus**3) / 3             # = (1-(13/15)^3)/3

# ====================================================================
# 2. CIBLE
# ====================================================================
R       = mp.mpf('1.06208066520087022')
R_tgt   = mp.mpf(17)/16
eps     = 1 - R/R_tgt                          # 3.9467e-4

print("="*72)
print("  E5b — Recherche correction NLO pour fermer ε = 1 - R/(17/16)")
print("="*72)
print(f"\n  R = {mp.nstr(R, 15)}")
print(f"  17/16 = {mp.nstr(R_tgt, 8)}")
print(f"  ε := 1 - R/(17/16) = {mp.nstr(eps, 12)}")
print(f"  ε ≈ {mp.nstr(eps*100, 6)} %    (1/ε ≈ {mp.nstr(1/eps, 6)})")

print("\n  Constantes PT (dps=50) :")
print(f"    sin²(θ_3,q+) = {mp.nstr(s2_3,10)}    γ_3 = {mp.nstr(g3,10)}")
print(f"    sin²(θ_5,q+) = {mp.nstr(s2_5,10)}    γ_5 = {mp.nstr(g5,10)}")
print(f"    sin²(θ_7,q+) = {mp.nstr(s2_7,10)}    γ_7 = {mp.nstr(g7,10)}")
print(f"    sin²(θ_11)   = {mp.nstr(s2_11,10)}    γ_11= {mp.nstr(g11,10)}")
print(f"    sin²(θ_13)   = {mp.nstr(s2_13,10)}    γ_13= {mp.nstr(g13,10)}")
print(f"    Δq = q_-q_+ = {mp.nstr(Delta_q,10)}    Δq² = {mp.nstr(Delta_q**2,10)}")
print(f"    α_EM(nu) = {mp.nstr(alpha_EM_nu,10)}    α_EM = {mp.nstr(alpha_EM,10)}")
print(f"    β_ghost = {mp.nstr(beta_ghost,10)}    δ_3 = {mp.nstr(delta_3,10)}")

# ====================================================================
# 3. CATALOGUE de candidates NLO plausibles
# ====================================================================
# Stratégie : structures déjà attestées en PT pour les corrections 1-loop
# - β_ghost · (small) : analogue ghost-VP de l'habillage α_EM
# - sin²·γ : structure intrinsèque des contributions per-channel
# - cross-channel cube : analogue à Koide NLO 1/(p1 p3) mais carré ici
# - structures K4-natives : Δq^n / (p2+μ*)^k
# - 2-loop électrofaibles : (α_EM)^2 · k

cands = []
def add(label, val):
    if val == 0: return
    cands.append((label, val, eps/val, abs(eps/val - 1)))

# ---- (A) Structures ghost VP (analogues α_EM dressing) ----
add('β_ghost · α_EM',                beta_ghost * alpha_EM)
add('β_ghost · α_EM_nu',             beta_ghost * alpha_EM_nu)
add('β_ghost · s² · α_EM',           beta_ghost * s**2 * alpha_EM)
add('β_ghost / (μ*·17)',             beta_ghost / (mu_star*17))
add('β_ghost · Δq²',                 beta_ghost * Delta_q**2)
add('β_ghost · Δq² / 2',             beta_ghost * Delta_q**2 / 2)
add('β_ghost · Δq² / (p2+μ*)',       beta_ghost * Delta_q**2 / 17)
add('β_ghost / (4π · μ*)',           beta_ghost / (4*mp.pi*mu_star))
add('β_ghost · sin²₇',               beta_ghost * s2_7)
add('(β_ghost)² · 17',               beta_ghost**2 * 17)

# ---- (B) Boucle "top-like" PT : sin²_p · γ_p (analogue F2/2-loop) ----
add('sin²₃ · γ_3 · α_EM',            s2_3*g3*alpha_EM)
add('sin²₇ · γ_7 / 17²',             s2_7*g7 / 289)
add('Σ sin²·γ {3,5,7} · α_EM',      (s2_3*g3+s2_5*g5+s2_7*g7)*alpha_EM)
add('Σ sin²·γ {3,5,7} / (p2·μ*)²',  (s2_3*g3+s2_5*g5+s2_7*g7) / (p2*mu_star)**2)
add('Σ sin²·γ {3,5,7} / 17²',       (s2_3*g3+s2_5*g5+s2_7*g7) / 289)

# ---- (C) Cross-canal Koide-like : 1/(p_i · p_j · k) ----
add('1/(p1·p2·p3 · μ*)',             mp.mpf(1)/(3*5*7*15))
add('1/(p1²·p3²)',                   mp.mpf(1)/(9*49))
add('1/(p1·p3·17)',                  mp.mpf(1)/(3*7*17))
add('1/(p1·p2·p3·17)',               mp.mpf(1)/(3*5*7*17))
add('5·δ_3²/(18·21·17)',             5*delta_3**2/(18*21*17))   # NNLO Koide / 17

# ---- (D) Structures K4-natives : Δq^n / (p2+μ*)^k ----
add('Δq² / 17',                      Delta_q**2 / 17)
add('Δq² / (p2+μ*)²',                Delta_q**2 / 289)
add('Δq³ / 17',                      Delta_q**3 / 17)
add('Δq² · N_br / 17²',              Delta_q**2 * N_br / 289)
add('Δq² · μ* / 17³',                Delta_q**2 * mu_star / 17**3)
add('Δq⁴ / R',                       Delta_q**4 / R)
add('Δq² · α_EM',                    Delta_q**2 * alpha_EM)

# ---- (E) Boucle bosonique W/Z analogues : (sin²W)² · k ----
sin2W_tree = g7**2 / (g3**2 + g5**2 + g7**2)
add('(sin²W)² / 17',                 sin2W_tree**2 / 17)
add('sin²W · α_EM',                  sin2W_tree * alpha_EM)
add('sin²W · α_EM_nu',               sin2W_tree * alpha_EM_nu)
add('(sin²W)² · α_EM',               sin2W_tree**2 * alpha_EM)

# ---- (F) 2-loop EW : α_EM² · k ----
add('α_EM² · 4π',                    alpha_EM**2 * 4*mp.pi)
add('α_EM² · 2π · 17',               alpha_EM**2 * 2*mp.pi * 17)
add('α_EM² · μ*²',                   alpha_EM**2 * mu_star**2)
add('α_EM²_nu · μ*²',                alpha_EM_nu**2 * mu_star**2)
add('α_EM² · 17',                    alpha_EM**2 * 17)

# ---- (G) C_Koide-like exporté à K4 (NLO=1/21, NNLO=5δ₃²/378) divisés par échelle ----
add('1/21 · α_EM',                   alpha_EM/21)
add('1/21 / μ*²',                    1/(21*mu_star**2))
add('5δ_3²/(18·21)',                 5*delta_3**2/(18*21))
add('5δ_3²/(18·21) · μ*',            5*delta_3**2/(18*21)*mu_star)
add('1/(21·N_br)',                   1/(21*N_br))

# ---- (H) Structures "ε ≈ 1/2533" : décompositions arithmétiques ----
add('1/(p2 · μ*² · N_br · p1)',      1/(p2*mu_star**2*N_br*3))   # 1/1350 != cible
add('1/(p1·p2·p3·17 + 17)',          1/(3*5*7*17+17))
add('1/(17² · γ_3 · γ_5 · γ_7)',     1/(289 * g3*g5*g7))
add('1/(17² · μ*·s²)',               1/(289 * mu_star * s**2))

# ---- (I) BPS / cardinal Borsuk-Ulam · Δq^n ----
add('p2 · Δq^4 / (p2+μ*)²',         p2 * Delta_q**4 / 289)
add('p2 · Δq² / (p2+μ*)³',          p2 * Delta_q**2 / 17**3)

# ====================================================================
# 4. Ranking
# ====================================================================
cands_sorted = sorted(cands, key=lambda x: x[3])
print("\n" + "="*72)
print("  Catalogue NLO — classement par |ε/candidate − 1| (croissant)")
print("="*72)
print(f"  {'Candidate':40s}  {'Valeur':18s}  {'ε/cand':12s}  {'|écart-1|':10s}")
for label, val, ratio, dev in cands_sorted[:20]:
    print(f"  {label:40s}  {mp.nstr(val,10):18s}  {mp.nstr(ratio,8):12s}  {mp.nstr(dev*100,4):8s}%")

# ====================================================================
# 5. Verdict
# ====================================================================
print("\n" + "="*72)
print("  Verdict")
print("="*72)
best = cands_sorted[0]
print(f"\n  Meilleur match : '{best[0]}'")
print(f"    Valeur     : {mp.nstr(best[1],12)}")
print(f"    ε/cand     : {mp.nstr(best[2],8)}    (idéal = 1.0)")
print(f"    Déviation  : {mp.nstr(best[3]*100,4)} %")

# Critère honnête : match < 1 %
top_1pct = [c for c in cands_sorted if c[3] < 0.01]
top_5pct = [c for c in cands_sorted if c[3] < 0.05]
print(f"\n  # candidats < 1 % de déviation : {len(top_1pct)}")
print(f"  # candidats < 5 % de déviation : {len(top_5pct)}")

verdict = "AUCUN_MATCH"
if len(top_1pct) >= 1:
    verdict = "MATCH_TROUVE"
elif len(top_5pct) >= 1:
    verdict = "APPROCHANT"

print(f"\n  VERDICT : {verdict}")

# ====================================================================
# 6. Dump JSON
# ====================================================================
out = {
    'eps': str(eps),
    'eps_pct': float(eps*100),
    'one_over_eps': str(1/eps),
    'best_candidate': best[0],
    'best_value': str(best[1]),
    'best_ratio_eps_over_cand': float(best[2]),
    'best_dev_pct': float(best[3]*100),
    'verdict': verdict,
    'n_under_1pct': len(top_1pct),
    'n_under_5pct': len(top_5pct),
    'top_5': [
        {'label': l, 'value': str(v), 'ratio': float(r), 'dev_pct': float(d*100)}
        for (l,v,r,d) in cands_sorted[:5]
    ]
}
out_dir = Path("/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT_PROJECTS/PT_NONCOM_BIFURCATION/outputs")
out_path = out_dir / "24_E5b_output.json"
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n  → JSON sauvegardé : {out_path}")
