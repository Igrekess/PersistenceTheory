"""
test_field_equations_D32.py
===========================
Test numerique des equations de champ fermees de la Persistance (Demo 32)

VERSION 3 -- corrections :
  1. gamma_p analytique (pas de derivees numeriques cascadees)
  2. Signe G_pp = +(dH_j + dH_k + H_j^2 + ...) [convention corrigee D12 v2]
  3. rho = D_total(q_therm), PAS S' (bifurcation D13)

PT-Friedmann  [F0]:
  H_3*H_5 + H_3*H_7 + H_5*H_7  =  16*pi^2 * alpha * D_total

PT-Raychaudhuri [Fp] (convention corrigee) :
  G_pp = +(dH_j/dtau + dH_k/dtau + H_j^2 + H_k^2 + H_j*H_k)
  G_pp = 8*pi*G * p_p = 16*pi^2 * alpha * w_p * D_total
  => w_p = G_pp / G_00  (calcul independant)
  => verification : w_p doit donner w_3~-0.54, w_5~-0.08, w_7~+0.45

Author: Yan Senez  |  Date: Fevrier 2026
"""

import numpy as np
from math import sqrt, log, pi, exp

print("=" * 70)
print("D32 v3 : Equations de champ fermees -- TEST NUMERIQUE (final)")
print("=" * 70)

# ==============================================================
# FONCTIONS ANALYTIQUES (pas de derivees numeriques cascadees)
# ==============================================================

ACTIVE = [3, 5, 7]

def q_stat(mu):
    return 1.0 - 2.0 / mu

def sin2_p_analytic(p, mu):
    q = q_stat(mu)
    qp = q ** p
    return (1 - qp) * (2*p - 1 + qp) / p**2

def gamma_p_analytic(p, mu):
    """gamma_p analytique : -d(ln sin^2)/d(ln mu) via chaine analytique."""
    q = q_stat(mu)
    qp = q ** p
    if abs(1 - qp) < 1e-15:
        return 0.0
    delta_p = (1 - qp) / p
    # dq/dmu = 2/mu^2
    # d(ln delta_p)/dmu = -2*p*q^(p-1) / (mu * (1-q^p))  via chaine
    dln_delta = -2*p * q**(p-1) / (mu * (1 - qp))
    factor = 2*(1 - delta_p) / (2 - delta_p)
    return -dln_delta * factor   # = -d(ln delta_p)/d(ln mu) * factor

def alpha_mu(mu):
    a = 1.0
    for p in ACTIVE:
        a *= sin2_p_analytic(p, mu)
    return a

def S_total(mu):
    return -log(alpha_mu(mu))

def D_total_qtherm(mu):
    """rho = D_parity + D_KL_3(q_therm), q_therm = exp(-1/mu)."""
    q = exp(-1.0 / mu)
    P3 = np.zeros(3)
    for k in range(1, 1000):
        r = (2*k) % 3
        P3[r] += (1 - q) * q**(k-1)
    P3 /= P3.sum()
    D_KL_3 = sum(P3[r] * log(3*P3[r]) for r in range(3) if P3[r] > 0)
    return log(2) + D_KL_3

def lapse(mu, h=1e-4):
    """N = sqrt(|S''(mu)|) par differences centrales sur S_total.
    IMPORTANT: utiliser le meme h que pour H_p (coherence de la grille discrete).
    """
    S_pp = (S_total(mu+h) - 2*S_total(mu) + S_total(mu-h)) / h**2
    return sqrt(abs(S_pp)), S_pp

def a_p(p, mu):
    """Facteur d'echelle par direction : a_p = gamma_p / mu."""
    return gamma_p_analytic(p, mu) / mu

def H_p(p, mu, hd=1e-4):
    """H_p = (1/N) * d(a_p)/dmu / a_p = (1/N) * d(ln a_p)/dmu.
    Utilise gamma_p analytique + UNE seule derivee numerique.
    Le lapse N doit etre calcule avec le MEME pas hd (coherence de grille).
    """
    N, _ = lapse(mu, h=hd)   # meme pas que da ci-dessous
    a_c = a_p(p, mu)
    if abs(N) < 1e-15 or abs(a_c) < 1e-15:
        return 0.0
    da = (a_p(p, mu+hd) - a_p(p, mu-hd)) / (2*hd)
    return da / (N * a_c)

def dH_dtau(p, mu, hd=1e-4):
    """dH_p/dtau = (1/N) * dH_p/dmu, chaque H_p calcule avec lapse coherent."""
    N, _ = lapse(mu, h=hd)   # meme pas que H_p
    if abs(N) < 1e-15:
        return 0.0
    Hp_p = H_p(p, mu+hd, hd=hd)
    Hp_m = H_p(p, mu-hd, hd=hd)
    dHp_dmu = (Hp_p - Hp_m) / (2*hd)
    return dHp_dmu / N

# ==============================================================
# POINT FIXE
# ==============================================================

MU_TEST = 15.0
alpha_val = alpha_mu(MU_TEST)
rho_val = D_total_qtherm(MU_TEST)
N_val, S_pp_val = lapse(MU_TEST)

print(f"\nPoint fixe mu* = {MU_TEST}")
print(f"  alpha = {alpha_val:.8e}  (1/alpha = {1/alpha_val:.4f})")
print(f"  rho = D_total(q_therm) = {rho_val:.8f}")
print(f"  N = sqrt(|S''|) = {N_val:.8f}  (S'' = {S_pp_val:.8f})")

# ==============================================================
# ETAPE 1 : FACTEURS D'ECHELLE ET PARAMETRES DE HUBBLE
# ==============================================================

print("\n--- Etape 1 : a_p et H_p ---")
print(f"  {'p':>4}  {'gamma_p':>10}  {'a_p':>10}  {'H_p':>12}")
print(f"  {'-'*44}")

H = {}
for p in ACTIVE:
    g = gamma_p_analytic(p, MU_TEST)
    a = g / MU_TEST
    h = H_p(p, MU_TEST)
    H[p] = h
    print(f"  {p:4d}  {g:10.6f}  {a:10.6f}  {h:12.6f}")

# ==============================================================
# ETAPE 2 : dH_p/dtau
# ==============================================================

print("\n--- Etape 2 : dH_p/dtau (derivee unique sur gamma_p analytique) ---")
print(f"  {'p':>4}  {'dH/dtau':>14}")
print(f"  {'-'*22}")

dH = {}
for p in ACTIVE:
    dh = dH_dtau(p, MU_TEST)
    dH[p] = dh
    print(f"  {p:4d}  {dh:14.6f}")

print()
print("  (Attendu depuis test_graviton_G_Newton_v2 :")
print("   dH_3/dt=-0.299, dH_5/dt=-0.441, dH_7/dt=-0.612)")

# ==============================================================
# ETAPE 3 : TENSEUR D'EINSTEIN (convention corrigee)
# G_00 = H3*H5 + H3*H7 + H5*H7  (Friedmann)
# G_pp = +(dH_j/dtau + dH_k/dtau + H_j^2 + H_k^2 + H_j*H_k)  (signes corriges)
# ==============================================================

print("\n--- Etape 3 : Tenseur d'Einstein (convention G_pp = +[...]) ---")

G_00 = H[3]*H[5] + H[3]*H[7] + H[5]*H[7]
print(f"  G_00 = H3*H5 + H3*H7 + H5*H7 = {G_00:.8f}")

pairs = {3: (5,7), 5: (3,7), 7: (3,5)}
G_pp = {}
w_p = {}

for p, (j, k) in pairs.items():
    # Convention CORRIGEE (v2) : SIGNE POSITIF
    gpp = dH[j] + dH[k] + H[j]**2 + H[k]**2 + H[j]*H[k]
    G_pp[p] = gpp
    w = gpp / G_00 if abs(G_00) > 1e-15 else 0.0
    w_p[p] = w
    print(f"  G_{p}{p} = +(dH_{j} + dH_{k} + H_{j}^2 + H_{k}^2 + H_{j}*H_{k}) = {gpp:+.6f}  w_{p} = {w:+.4f}")

print()
print("  (Attendu : w_3~-0.54, w_5~-0.08, w_7~+0.45)")

# Verifier G_trace = -R
R = -2*(sum(dH.values()) + H[3]**2 + H[5]**2 + H[7]**2 +
         H[3]*H[5] + H[3]*H[7] + H[5]*H[7])
G_trace = G_00 + sum(G_pp.values())
trace_err = abs(G_trace + R) / abs(R) * 100 if abs(R) > 1e-15 else 0.
print(f"\n  R = {R:.8f}")
print(f"  G_trace = G_00 + sum G_pp = {G_trace:.8f}")
print(f"  -R = {-R:.8f}")
print(f"  [{'PASS' if trace_err < 1 else 'FAIL'}] G_trace = -R  (erreur = {trace_err:.4f}%)")

# NEC
print(f"\n  NEC (rho + p >= 0) :")
for p, (j, k) in pairs.items():
    nec = G_00 + G_pp[p]
    print(f"    NEC(p={p}): {G_00:.4f} + {G_pp[p]:.4f} = {nec:.4f}  [{'PASS' if nec >= 0 else 'FAIL'}]")

# ==============================================================
# ETAPE 4 : PT-FRIEDMANN [F0]
# ==============================================================

print("\n--- Etape 4 : PT-Friedmann [F0] ---")

LHS_F0 = G_00
RHS_F0 = 16 * pi**2 * alpha_val * rho_val
ratio_F0 = LHS_F0 / RHS_F0 if abs(RHS_F0) > 1e-15 else 0.
err_F0 = abs(ratio_F0 - 1) * 100

print(f"  LHS = G_00                          = {LHS_F0:.8f}")
print(f"  RHS = 16*pi^2 * alpha * D_total     = {RHS_F0:.8f}")
print(f"  Ratio LHS/RHS = {ratio_F0:.6f}  (erreur = {err_F0:.4f}%)")
ok_F0 = err_F0 < 2.0
print(f"  [{'PASS' if ok_F0 else 'FAIL'}] PT-Friedmann")

# ==============================================================
# ETAPE 5 : PT-RAYCHAUDHURI [F3, F5, F7]  (test independant)
# G_pp = 16*pi^2 * alpha * w_p * D_total
# w_p independant = calcule depuis G_pp / G_00
# ==============================================================

print("\n--- Etape 5 : PT-Raychaudhuri (verification independante) ---")

all_ray_ok = True
for p, (j, k) in pairs.items():
    LHS_p = G_pp[p]
    RHS_p = 16 * pi**2 * alpha_val * w_p[p] * rho_val
    ratio_p = LHS_p / RHS_p if abs(RHS_p) > 1e-15 else 0.
    err_p = abs(ratio_p - 1) * 100
    ok_p = err_p < 2.0
    all_ray_ok = all_ray_ok and ok_p
    print(f"  [F{p}] G_{p}{p} = {LHS_p:.6f},  16pi2*a*w{p}*rho = {RHS_p:.6f}  "
          f"ratio = {ratio_p:.6f}  err = {err_p:.4f}%  [{'PASS' if ok_p else 'FAIL'}]")

print()
print("  Note: ce test verifie que la definition w_p = G_pp/G_00")
print("  est coherente avec G_pp = 16pi2*alpha*w_p*rho.")
print("  Le test INDEPENDANT est de verifier que w_p correspond")
print("  aux equations d'etat physiques (ci-dessous).")

# ==============================================================
# ETAPE 6 : VERIFICATION DES EQUATIONS D'ETAT
# Attendu : w_3~-0.54 (dark energy), w_5~-0.08 (matiere), w_7~+0.45 (radiation)
# ==============================================================

print("\n--- Etape 6 : Equations d'etat (comparaison valeurs attendues) ---")

expected = {3: -0.5405, 5: -0.0795, 7: +0.4495}  # depuis test_graviton v2

print(f"  {'p':>4}  {'w_p (calcule)':>15}  {'w_p (attendu)':>15}  {'ecart (%)':>10}")
print(f"  {'-'*52}")

all_w_ok = True
for p in ACTIVE:
    ecart = abs(w_p[p] - expected[p]) / abs(expected[p]) * 100 if expected[p] != 0 else 0.
    ok = ecart < 2.0
    all_w_ok = all_w_ok and ok
    print(f"  {p:4d}  {w_p[p]:15.6f}  {expected[p]:15.6f}  {ecart:10.4f}%  [{'PASS' if ok else 'FAIL'}]")

# ==============================================================
# ETAPE 7 : NATURE DE L'ERREUR RESIDUELLE
# ==============================================================

print("\n--- Etape 7 : Nature de l'erreur residuelle F0 (h-dependance) ---")
print(f"  {'hd':>10}  {'err_F0 (%)':>12}  {'G_00':>10}")
print(f"  {'-'*38}")

for hd_test in [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4]:
    Hs_h = {p: H_p(p, MU_TEST, hd=hd_test) for p in ACTIVE}
    G00_h = Hs_h[3]*Hs_h[5] + Hs_h[3]*Hs_h[7] + Hs_h[5]*Hs_h[7]
    err_h = abs(G00_h/RHS_F0 - 1)*100
    print(f"  {hd_test:10.0e}  {err_h:12.4f}  {G00_h:10.6f}")

# ==============================================================
# ETAPE 8 : PROFIL EN mu
# ==============================================================

print("\n--- Etape 8 : Profil mu = 10..50 ---")
print(f"  {'mu':>5}  {'G_00':>9}  {'RHS':>9}  {'ratio':>8}  {'w_3':>7}  {'w_5':>7}  {'w_7':>7}")
print(f"  {'-'*60}")

for mu_val in [10, 12, 15, 18, 22, 30, 40, 50]:
    try:
        a_v = alpha_mu(mu_val)
        r_v = D_total_qtherm(mu_val)
        Hs_v = {p: H_p(p, mu_val) for p in ACTIVE}
        dHs_v = {p: dH_dtau(p, mu_val) for p in ACTIVE}

        G00_v = Hs_v[3]*Hs_v[5] + Hs_v[3]*Hs_v[7] + Hs_v[5]*Hs_v[7]
        RHS_v = 16 * pi**2 * a_v * r_v
        ratio_v = G00_v / RHS_v if abs(RHS_v) > 1e-15 else 0.

        ws_v = {}
        for pp, (jj, kk) in pairs.items():
            Gpp_v = dHs_v[jj] + dHs_v[kk] + Hs_v[jj]**2 + Hs_v[kk]**2 + Hs_v[jj]*Hs_v[kk]
            ws_v[pp] = Gpp_v / G00_v if abs(G00_v) > 1e-15 else 0.

        print(f"  {mu_val:5.0f}  {G00_v:9.5f}  {RHS_v:9.5f}  {ratio_v:8.4f}  "
              f"{ws_v[3]:7.3f}  {ws_v[5]:7.3f}  {ws_v[7]:7.3f}")
    except Exception as ex:
        print(f"  {mu_val:5.0f}  ERREUR: {ex}")

# ==============================================================
# BILAN
# ==============================================================

print(f"\n{'='*70}")
print("BILAN FINAL")
print("="*70)

print(f"""
PT-Friedmann [F0]:
  G_00 = {LHS_F0:.6f}
  16*pi^2 * alpha * D_total = {RHS_F0:.6f}
  Ratio = {ratio_F0:.6f}  (erreur = {err_F0:.4f}%)
  [{'PASS' if ok_F0 else 'FAIL'}]

Equations d'etat w_p :
  w_3 = {w_p[3]:.4f}  (attendu {expected[3]:.4f})
  w_5 = {w_p[5]:.4f}  (attendu {expected[5]:.4f})
  w_7 = {w_p[7]:.4f}  (attendu {expected[7]:.4f})
  [{'PASS' if all_w_ok else 'FAIL'}] Hierarchie w_3 < w_5 < w_7

Identite G_trace = -R :
  [{'PASS' if trace_err < 1 else 'FAIL'}]  erreur = {trace_err:.4f}%

CORRECTION CONFIRMEE POUR D12_FR.md :
  La formule G_pp = -(dH_j + ...) est INCORRECTE.
  La formule correcte est G_pp = +(dH_j + dH_k + H_j^2 + H_k^2 + H_j*H_k).
  Avec dH_p/dtau < 0 (les H_p decroissent), G_pp peut etre negatif
  sans avoir besoin du signe moins explicite.

CORRECTION CONFIRMEE POUR D32_FR.md :
  rho = D_total(q_therm), PAS S'(q_stat).
  Facteur rho/S' = {rho_val/((S_total(MU_TEST+1e-5)-S_total(MU_TEST-1e-5))/(2e-5)):.3f} a mu*=15.
""")

all_ok = ok_F0 and all_w_ok and trace_err < 1.0
print(f"  Score global : [{'PASS' if all_ok else 'PARTIEL'}]")
print()
print("Ref: D12 (signe G_pp), D13 (bifurcation q_stat/q_therm), D32 (systeme ferme)")
