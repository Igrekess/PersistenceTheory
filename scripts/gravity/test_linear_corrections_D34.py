"""
test_linear_corrections_D34.py
==============================
Test numerique de la linearite des corrections et additivite des residus (Demo 34)

Structure verifiee :
  delta_total(mu_alpha) = delta_arith(mu*) + delta_EM(mu_alpha)
                        = 0.196%           + 0.730%
                        = 0.926%   [accord : 0.12%]

La pente d(G_00/RHS)/dmu est quasi-constante (variation < 3.7%) sur [mu_zero, mu_alpha].

Author: Yan Senez  |  Date: Fevrier 2026
"""

import numpy as np
from math import sqrt, log, pi, exp
from scipy.optimize import brentq

print("=" * 70)
print("D34 : Linearite des corrections PT -- TEST NUMERIQUE")
print("=" * 70)

# ==============================================================
# FONCTIONS (coherentes avec test_field_equations_D32.py v3)
# ==============================================================

ACTIVE = [3, 5, 7]

def q_stat(mu): return 1.0 - 2.0/mu

def sin2_p(p, mu):
    q = q_stat(mu); qp = q**p
    return (1 - qp)*(2*p - 1 + qp)/p**2

def gamma_p_analytic(p, mu):
    q = q_stat(mu); qp = q**p
    if abs(1 - qp) < 1e-15: return 0.0
    delta_p = (1 - qp)/p
    dln_delta = -2*p*q**(p-1)/(mu*(1 - qp))
    factor = 2*(1 - delta_p)/(2 - delta_p)
    return -dln_delta * factor

def alpha_mu(mu):
    a = 1.0
    for p in ACTIVE:
        a *= sin2_p(p, mu)
    return a

def S_total(mu): return -log(alpha_mu(mu))

def D_total_qtherm(mu):
    q = exp(-1.0/mu)
    P3 = np.zeros(3)
    for k in range(1, 1000):
        P3[(2*k) % 3] += (1 - q)*q**(k-1)
    P3 /= P3.sum()
    return log(2) + sum(P3[r]*log(3*P3[r]) for r in range(3) if P3[r] > 0)

def lapse_h(mu, h=1e-4):
    S_pp = (S_total(mu+h) - 2*S_total(mu) + S_total(mu-h))/h**2
    return sqrt(abs(S_pp))

def a_p(p, mu): return gamma_p_analytic(p, mu)/mu

def H_p_func(p, mu, hd=1e-4):
    N = lapse_h(mu, hd)
    ac = a_p(p, mu)
    da = (a_p(p, mu+hd) - a_p(p, mu-hd))/(2*hd)
    return da/(N*ac)

def G00(mu):
    Hs = {p: H_p_func(p, mu) for p in ACTIVE}
    return Hs[3]*Hs[5] + Hs[3]*Hs[7] + Hs[5]*Hs[7]

def RHS_F0(mu):
    return 16*pi**2 * alpha_mu(mu) * D_total_qtherm(mu)

def ratio(mu): return G00(mu) / RHS_F0(mu)

def residual_func(mu): return G00(mu) - RHS_F0(mu)

# ==============================================================
# ETAPE 1 : LES TROIS POINTS
# ==============================================================

print("\n--- Etape 1 : Les trois points ---")

mu_zero = brentq(residual_func, 14.9, 15.05, xtol=1e-10)
alpha_EM = 1/137.035999084
mu_alpha = brentq(lambda m: alpha_mu(m) - alpha_EM, 14.5, 16.0)

Delta1 = 15.0 - mu_zero
Delta2 = mu_alpha - 15.0

for mu_v, lbl in [(mu_zero, "mu_zero"), (15.0, "mu*"), (mu_alpha, "mu_alpha")]:
    r = ratio(mu_v)
    print(f"  {lbl:12s} ({mu_v:.5f}) : f={r:.7f}  residu={(r-1)*100:+.4f}%")

print(f"\n  Delta1 = mu* - mu_zero  = {Delta1:.5f}")
print(f"  Delta2 = mu_alpha - mu* = {Delta2:.5f}")
print(f"  Delta2/Delta1           = {Delta2/Delta1:.4f}")

# ==============================================================
# ETAPE 2 : PENTE QUASI-CONSTANTE
# ==============================================================

print("\n--- Etape 2 : Pente locale d(G_00/RHS)/dmu ---")
h_s = 0.03
pentes = []
for mu_v in [14.70, 14.85, 14.95, 15.00, 15.02, 15.04]:
    s = (ratio(mu_v+h_s) - ratio(mu_v-h_s)) / (2*h_s)
    pentes.append(s)
    print(f"  mu={mu_v:.4f}  pente={s:.5f} /mu")

variation = (max(pentes) - min(pentes)) / min(pentes) * 100
pente_moy = np.mean(pentes)
print(f"\n  Variation totale : {variation:.1f}%  (quasi-lineaire si < 5%)")
print(f"  Pente moyenne    : {pente_moy:.5f} /mu")
linearite_ok = variation < 5.0
print(f"  Linearite        : {'PASS' if linearite_ok else 'FAIL'} (<5% variation)")

# ==============================================================
# ETAPE 3 : ADDITIVITE
# ==============================================================

print("\n--- Etape 3 : Structure additive ---")

pente_ref = pentes[3]  # pente a mu*

r_zero = (ratio(mu_zero) - 1) * 100
r_star = (ratio(15.0) - 1) * 100
r_alpha = (ratio(mu_alpha) - 1) * 100

pred_arith = pente_ref * Delta1 * 100
pred_em    = pente_ref * Delta2 * 100
pred_total = pred_arith + pred_em

print(f"\n  Couche 1 : mu_zero -> mu*  (arithmetique D33)")
print(f"    Prediction : {pente_ref:.4f} * {Delta1:.5f} = {pred_arith:.4f}%")
print(f"    Mesure     : {r_star:.4f}%")
print(f"    Accord     : {abs(pred_arith - r_star)/abs(r_star)*100:.2f}%")

print(f"\n  Couche 2 : mu* -> mu_alpha  (electromagnetique D09)")
mesure_em = r_alpha - r_star
print(f"    Prediction : {pente_ref:.4f} * {Delta2:.5f} = {pred_em:.4f}%")
print(f"    Mesure     : {mesure_em:.4f}%")
print(f"    Accord     : {abs(pred_em - mesure_em)/abs(mesure_em)*100:.2f}%")

print(f"\n  Somme additive (mu_zero -> mu_alpha) :")
print(f"    Prediction : {pred_arith:.4f}% + {pred_em:.4f}% = {pred_total:.4f}%")
print(f"    Mesure     : {r_alpha:.4f}%")
accord_total = abs(pred_total - r_alpha)/abs(r_alpha)*100
print(f"    Accord     : {accord_total:.3f}%")
additif_ok = accord_total < 0.5
print(f"    Additivite : {'PASS' if additif_ok else 'FAIL'} (<0.5%)")

# ==============================================================
# ETAPE 4 : EQUIVALENCE D09
# ==============================================================

print("\n--- Etape 4 : Equivalence avec D09 ---")
delta_inv_alpha = 1/alpha_mu(mu_alpha) - 1/alpha_mu(15.0)
print(f"  alpha_nue = 1/{1/alpha_mu(15.0):.4f}  (mu*=15)")
print(f"  alpha_EM  = 1/{1/alpha_mu(mu_alpha):.6f}  (mu_alpha)")
print(f"  delta(1/alpha) = {delta_inv_alpha:.4f}  [D09 : 0.758]")
print(f"  Deplacement equivalent : Dmu_EM = {Delta2:.5f}")
print(f"  Dual D09/deplacement PT : {'PASS' if abs(delta_inv_alpha - 0.758) < 0.01 else 'CHECK'}")

# ==============================================================
# BILAN
# ==============================================================

print(f"\n{'='*70}")
print("BILAN")
print("="*70)
print(f"""
Trois points :
  mu_zero  = {mu_zero:.5f}   residu = {r_zero:.4f}%
  mu*      = 15.00000   residu = {r_star:.4f}%
  mu_alpha = {mu_alpha:.5f}   residu = {r_alpha:.4f}%

Linearite : variation de pente = {variation:.1f}%  {'[PASS]' if linearite_ok else '[FAIL]'}

Structure additive :
  delta_arith = {pred_arith:.4f}%  (couche arithmetique, D33)
  delta_EM    = {pred_em:.4f}%  (couche EM, D09 equivalent)
  Somme       = {pred_total:.4f}%  vs mesure {r_alpha:.4f}%  {'[PASS]' if additif_ok else '[FAIL]'}

Rapport Delta2/Delta1 = {Delta2/Delta1:.4f}  (la correction EM est 3.72x la correction arithmetique)

Interpretation :
  La fonction f(mu)=G_00/RHS est quasi-lineaire sur [mu_zero, mu_alpha].
  Les deux corrections (arithmetique D33, EM D09) s'accumulent de facon additive.
  La correction D09 correspond a un deplacement Dmu={Delta2:.5f} dans l'espace PT.
""")
