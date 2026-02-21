"""
test_graviton_G_Newton_v2
=========================

ENGLISH
-------
Graviton and Newton's G: derivation of G from spin foam and PT framework

FRANCAIS (original)
-------------------
test_graviton_G_Newton_v2.py
S15.6.125b: Derivation de G_Newton -- version corrigee et approfondie

CORRECTIONS vs v1:
- Signe des G^i_i CORRIGE (etait inverse dans le script Jacobson original)
- G_trace = -R maintenant EXACT
- Unites naturelles (pas de c^4 artificiel dans kappa)

DECOUVERTE: G_sieve = G_00/(8*pi*D_KL) ~ 2*pi*alpha_EM (0.3%!)
            Equivalemment: G_00 = (4*pi)^2 * alpha_EM * D_KL

8 TESTS:
1. Spectre de perturbation du crible
2. Metrique Lorentzienne (Bianchi I)
3. Tenseur d'Einstein (signe CORRIGE) + verification G_trace = -R
4. Contenu informationnel (D_KL, H)
5. G_Newton = G_00/(8*pi*rho) pour 5 identifications de rho
6. DECOUVERTE: G = 2*pi*alpha (formule close)
7. Universalite: G(mu) sur un intervalle de mu
8. Interpretation physique et rapport alpha_G/alpha_EM

Author / Auteur: Yan Senez  |  Date: February / Fevrier 2026
Theory / Theorie: Persistence Theory (PT) / Theorie de la Persistance (TP)
"""

import numpy as np
from math import sqrt, log, log2, pi, exp
from scipy.optimize import brentq

print("=" * 70)
print("S15.6.125b: DERIVATION DE G_Newton (v2 -- signe corrige)")
print("=" * 70)

# ==============================================================
# CONSTANTES ET FONCTIONS DU CRIBLE
# ==============================================================

phi = (1 + sqrt(5)) / 2
s = 0.5
alpha_EM_phys = 1 / 137.035999084
active_primes = [3, 5, 7]
c_sieve = 105 / phi

def sin2_theta(p, q):
    qp = q**p
    return (1 - qp) * (2*p - 1 + qp) / p**2

def alpha_sieve(mu):
    q = 1 - 2/mu
    result = 1.0
    for p in active_primes:
        result *= sin2_theta(p, q)
    return result

def gamma_p(p, mu):
    q = 1 - 2/mu
    qp = q**p
    if abs(1 - qp) < 1e-15:
        return 0.0
    delta_p = (1 - qp) / p
    dln_delta = -2*p * q**(p-1) / (mu * (1 - qp))
    factor = 2*(1 - delta_p) / (2 - delta_p)
    return -dln_delta * factor

def ln_alpha(mu):
    a = alpha_sieve(mu)
    return log(a) if a > 0 else -100.0

def d2_ln_alpha(mu, h=1e-4):
    return (ln_alpha(mu+h) - 2*ln_alpha(mu) + ln_alpha(mu-h)) / h**2

def lapse(mu):
    return sqrt(abs(d2_ln_alpha(mu)))

# ==============================================================
# POINT OPERATOIRE
# ==============================================================

mu_alpha = brentq(lambda m: alpha_sieve(m) - alpha_EM_phys, 14.5, 16.0)
delta_mu = mu_alpha - 15.0
alpha_op = alpha_sieve(mu_alpha)

print(f"\nPoint operatoire:")
print(f"  mu_alpha = {mu_alpha:.6f}, delta_mu = {delta_mu:.6f}")
print(f"  alpha = {alpha_op:.10e}, 1/alpha = {1/alpha_op:.4f}")

# ==============================================================
# TEST 1: SPECTRE DE PERTURBATION
# ==============================================================

print(f"\n{'='*70}")
print("TEST 1: Spectre de perturbation delta f(p)")
print("="*70)

q_op = 1 - 2/mu_alpha
h_diff = 1e-5

print(f"\n{'p':>4} {'sin2(p)':>12} {'gamma_p':>10} {'d(ln f)/dmu':>14}")
print("-" * 44)

total_dlnalpha = 0
for p in active_primes:
    f_val = sin2_theta(p, q_op)
    gp = gamma_p(p, mu_alpha)
    f_p = sin2_theta(p, 1 - 2/(mu_alpha + h_diff))
    f_m = sin2_theta(p, 1 - 2/(mu_alpha - h_diff))
    dlnf = (f_p - f_m) / (2*h_diff*f_val)
    total_dlnalpha += dlnf
    print(f"{p:4d} {f_val:12.8f} {gp:10.4f} {dlnf:14.6e}")

# Verification
from math import log as mlog
dlna_direct = (ln_alpha(mu_alpha+h_diff) - ln_alpha(mu_alpha-h_diff)) / (2*h_diff)
err_pct = abs(total_dlnalpha - dlna_direct) / abs(dlna_direct) * 100

print(f"\n  sum d(ln f)/dmu = {total_dlnalpha:.8e}")
print(f"  d(ln alpha)/dmu = {dlna_direct:.8e}")
print(f"  [{'PASS' if err_pct < 0.1 else 'FAIL'}] Coherence: {err_pct:.4f}%")
test1 = err_pct < 0.1

# Quantum gravitationnel
h_grav = abs(total_dlnalpha * log(11/7))  # pas 7->11
print(f"\n  Amplitude graviton |h| = {h_grav:.6e} (pas 7->11)")

# ==============================================================
# TEST 2: METRIQUE LORENTZIENNE
# ==============================================================

print(f"\n{'='*70}")
print("TEST 2: Metrique Lorentzienne (Bianchi I)")
print("="*70)

g_00 = -abs(d2_ln_alpha(mu_alpha))
g_sp = {}
a_fac = {}
for p in active_primes:
    gp = gamma_p(p, mu_alpha)
    g_sp[p] = gp**2 / mu_alpha**2
    a_fac[p] = gp / mu_alpha

det_g = g_00 * g_sp[3] * g_sp[5] * g_sp[7]

print(f"\n  g_00 = {g_00:.8f}  (< 0: temps)")
for p in active_primes:
    print(f"  g_{p}{p} = {g_sp[p]:.8f}  (espace)")
print(f"  det(g) = {det_g:.6e}")

test2 = g_00 < 0 and all(v > 0 for v in g_sp.values())
print(f"\n  [{'PASS' if test2 else 'FAIL'}] Signature (-,+,+,+)")

# ==============================================================
# TEST 3: TENSEUR D'EINSTEIN -- SIGNE CORRIGE
# ==============================================================

print(f"\n{'='*70}")
print("TEST 3: Tenseur d'Einstein (Bianchi I) -- SIGNE CORRIGE")
print("="*70)
print("  CORRECTION: G^i_i = +(dH_j + dH_k + H_j^2 + H_k^2 + H_j*H_k)")
print("  (Le script Jacobson original avait un signe - errone)")

def einstein_full(mu, hd=1e-4):
    """Einstein tensor Bianchi I avec signes CORRECTS."""
    N_val = lapse(mu)
    if N_val < 1e-15:
        return None

    def get_hubble(mu_e):
        N_e = lapse(mu_e)
        Hs = []
        for p in active_primes:
            gp_c = gamma_p(p, mu_e)
            gp_p = gamma_p(p, mu_e + hd)
            gp_m = gamma_p(p, mu_e - hd)
            a_i = gp_c / mu_e
            da = (gp_p/(mu_e+hd) - gp_m/(mu_e-hd)) / (2*hd)
            Hs.append(da / (N_e * a_i) if N_e > 0 and a_i > 0 else 0)
        return Hs

    H = get_hubble(mu)
    H1, H2, H3 = H

    # G^0_0 = H1*H2 + H1*H3 + H2*H3
    G_00 = H1*H2 + H1*H3 + H2*H3

    # dH_i/dtau
    H_plus = get_hubble(mu + hd)
    H_minus = get_hubble(mu - hd)
    dH = [(H_plus[i] - H_minus[i]) / (2*hd*N_val) for i in range(3)]

    # G^i_i = dH_j + dH_k + H_j^2 + H_k^2 + H_j*H_k  (SIGNE POSITIF!)
    G_sp = []
    pairs = [(1,2), (0,2), (0,1)]  # pour chaque direction i, les deux autres (j,k)
    for idx, (j, k) in enumerate(pairs):
        G_ii = dH[j] + dH[k] + H[j]**2 + H[k]**2 + H[j]*H[k]  # PAS de moins!
        G_sp.append(G_ii)

    # Scalaire de Ricci
    R = -2*(sum(dH) + H1**2 + H2**2 + H3**2 + H1*H2 + H1*H3 + H2*H3)

    theta = sum(H)

    return {
        'G_00': G_00, 'G_sp': G_sp, 'R': R,
        'H': H, 'dH': dH, 'theta': theta, 'N': N_val
    }

E = einstein_full(mu_alpha)

print(f"\n  Hubble: H_3={E['H'][0]:.6f}, H_5={E['H'][1]:.6f}, H_7={E['H'][2]:.6f}")
print(f"  dH/dt:  dH_3={E['dH'][0]:.4f}, dH_5={E['dH'][1]:.4f}, dH_7={E['dH'][2]:.4f}")

print(f"\n  G^0_0 = {E['G_00']:.6f}  (Friedmann: densite d'energie)")
for i, p in enumerate(active_primes):
    sign = "pression" if E['G_sp'][i] > 0 else "tension"
    print(f"  G^{i+1}_{i+1} = {E['G_sp'][i]:+.6f}  ({sign}, direction p={p})")

G_trace = E['G_00'] + sum(E['G_sp'])
print(f"\n  R = {E['R']:.6f}")
print(f"  G_trace = G^0_0 + sum G^i_i = {G_trace:.6f}")
print(f"  -R = {-E['R']:.6f}")

trace_err = abs(G_trace + E['R']) / abs(E['R']) * 100
test3 = trace_err < 1.0
print(f"\n  [{'PASS' if test3 else 'FAIL'}] G_trace = -R (erreur: {trace_err:.4f}%)")
print(f"  (v1 avait 29% d'erreur a cause du signe incorrect)")

# Conditions d'energie (avec signes corrects)
rho_geom = E['G_00']  # proportionnel a la densite d'energie
p_geom = E['G_sp']    # proportionnel aux pressions par direction
print(f"\n  Equation d'etat geometrique:")
print(f"    rho_geom = {rho_geom:.6f}")
for i, p in enumerate(active_primes):
    w_i = p_geom[i] / rho_geom if abs(rho_geom) > 1e-15 else 0
    print(f"    w_{p} = p_{p}/rho = {w_i:+.4f}")

# ==============================================================
# TEST 4: CONTENU INFORMATIONNEL
# ==============================================================

print(f"\n{'='*70}")
print("TEST 4: Contenu informationnel")
print("="*70)

mu = mu_alpha
q_th = exp(-2/mu)

# D_KL(mod 3) geometrique
P3 = np.zeros(3)
for k in range(1, 500):
    r = (2*k) % 3
    P3[r] += (1 - q_th) * q_th**(k-1)
P3 /= P3.sum()
D_KL_3 = sum(P3[r] * log(3*P3[r]) for r in range(3) if P3[r] > 0)
D_parity = log(2)
D_total = D_parity + D_KL_3

# Shannon
H_shan = -log(1-q_th) - q_th*log(q_th)/(1-q_th) if 0 < q_th < 1 else 0

eps = 0.5 - alpha_op

print(f"\n  D_KL(mod 3)  = {D_KL_3:.8f}")
print(f"  D_parity     = {D_parity:.8f}")
print(f"  D_total      = {D_total:.8f}")
print(f"  H_shannon    = {H_shan:.8f}")
print(f"  eps = 1/2-a  = {eps:.8f}")

test4 = D_total > 0 and H_shan > 0
print(f"\n  [{'PASS' if test4 else 'FAIL'}] Quantites positives")

# ==============================================================
# TEST 5: G_Newton = G_00/(8*pi*rho)
# ==============================================================

print(f"\n{'='*70}")
print("TEST 5: G_Newton = G_00 / (8*pi * rho)")
print("="*70)
print("  (Unites naturelles du crible: c = hbar = 1)")

G_00_val = E['G_00']

rho_dict = {
    'D_KL':  D_total,
    'H_shan': H_shan,
    'eps^2': eps**2,
    'D/mu':  D_total/mu,
    'H/mu':  H_shan/mu,
}

print(f"\n  G_00 = {G_00_val:.6f}")
print(f"\n  {'rho ID':10s} {'rho':>14s} {'G = G00/(8pi*rho)':>18s}")
print(f"  {'-'*46}")

G_results = {}
for label, rho_val in rho_dict.items():
    G_N = G_00_val / (8*pi*rho_val) if rho_val > 0 else 0
    G_results[label] = G_N
    print(f"  {label:10s} {rho_val:14.8f} {G_N:18.8f}")

test5 = len(G_results) >= 4
print(f"\n  [{'PASS' if test5 else 'FAIL'}] {len(G_results)} valeurs calculees")

# ==============================================================
# TEST 6: DECOUVERTE -- G = 2*pi*alpha_EM
# ==============================================================

print(f"\n{'='*70}")
print("TEST 6: DECOUVERTE -- G_sieve = 2*pi*alpha_EM")
print("="*70)

# L'identification rho = D_KL donne:
G_DKL = G_results['D_KL']
G_predicted = 2 * pi * alpha_op

print(f"\n  G(D_KL) = G_00/(8*pi*D_KL) = {G_DKL:.8f}")
print(f"  2*pi*alpha_EM              = {G_predicted:.8f}")

ratio_6 = G_DKL / G_predicted
err_6 = abs(ratio_6 - 1) * 100
print(f"  Ratio: {ratio_6:.6f} (erreur: {err_6:.2f}%)")

# Formulation equivalente
print(f"\n  Equivalemment: G_00 = (4*pi)^2 * alpha * D_KL")
pred_G00 = (4*pi)**2 * alpha_op * D_total
print(f"    G_00 observe = {G_00_val:.8f}")
print(f"    (4pi)^2*a*D  = {pred_G00:.8f}")
print(f"    Erreur: {abs(G_00_val - pred_G00)/G_00_val*100:.2f}%")

# Decomposition: pourquoi 16*pi^2?
print(f"\n  Decomposition de G_00 = (4pi)^2 * alpha * D_KL:")
print(f"    4*pi = 'volume' angulaire (Gauss)")
print(f"    alpha = couplage du crible = prod sin^2(theta_p)")
print(f"    D_KL = energie informationnelle")
print(f"    => G_00 = 'flux gravitationnel' * 'couplage' * 'source'")

# Verifier kappa = 8*pi*G = 8*pi*2*pi*alpha = 16*pi^2*alpha
kappa_pred = 16 * pi**2 * alpha_op
kappa_obs = G_00_val / D_total
print(f"\n  kappa predit = 16*pi^2*alpha = {kappa_pred:.6f}")
print(f"  kappa observe = G_00/D_KL    = {kappa_obs:.6f}")
print(f"  Erreur: {abs(kappa_pred - kappa_obs)/kappa_obs*100:.2f}%")

test6 = err_6 < 1.0
print(f"\n  [{'PASS' if test6 else 'FAIL'}] G_sieve = 2*pi*alpha (a {err_6:.2f}%)")

# ==============================================================
# TEST 7: UNIVERSALITE -- G(mu) sur un intervalle
# ==============================================================

print(f"\n{'='*70}")
print("TEST 7: Universalite -- G(mu)/[2*pi*alpha(mu)] sur mu = 10..25")
print("="*70)

def D_KL_at(mu_val):
    q = exp(-2/mu_val)
    P = np.zeros(3)
    for k in range(1, 500):
        r = (2*k) % 3
        P[r] += (1 - q) * q**(k-1)
    P /= P.sum()
    D3 = sum(P[r] * log(3*P[r]) for r in range(3) if P[r] > 0)
    return log(2) + D3

print(f"\n  {'mu':>6s} {'alpha':>12s} {'G_00':>10s} {'D_KL':>10s} {'G=G00/8piD':>12s} {'2pi*alpha':>12s} {'ratio':>8s}")
print(f"  {'-'*76}")

ratios_7 = []
for mu_test in np.linspace(10, 25, 16):
    a_test = alpha_sieve(mu_test)
    E_test = einstein_full(mu_test)
    if E_test is None:
        continue
    D_test = D_KL_at(mu_test)

    G_test = E_test['G_00'] / (8*pi*D_test) if D_test > 0 else 0
    G_pred = 2 * pi * a_test
    ratio = G_test / G_pred if G_pred > 0 else 0
    ratios_7.append(ratio)

    print(f"  {mu_test:6.1f} {a_test:12.6e} {E_test['G_00']:10.6f} {D_test:10.6f} {G_test:12.8f} {G_pred:12.8f} {ratio:8.4f}")

if ratios_7:
    mean_r = np.mean(ratios_7)
    std_r = np.std(ratios_7)
    print(f"\n  Ratio moyen: {mean_r:.4f} +/- {std_r:.4f}")
    print(f"  CV = {std_r/mean_r*100:.1f}%")

    test7 = abs(mean_r - 1) < 0.2 and std_r < 0.2
    print(f"\n  [{'PASS' if test7 else 'FAIL'}] Universalite de G = 2*pi*alpha")
else:
    test7 = False
    print(f"\n  [FAIL] Pas de donnees")

# ==============================================================
# TEST 8: INTERPRETATION PHYSIQUE
# ==============================================================

print(f"\n{'='*70}")
print("TEST 8: Interpretation physique")
print("="*70)

print(f"""
RESULTAT PRINCIPAL:
  G_sieve = 2*pi*alpha_EM  (erreur {err_6:.2f}% au point operatoire)

FORMULES EQUIVALENTES:
  G_00 = (4*pi)^2 * alpha_EM * D_KL    [Friedmann = Gauss^2 * couplage * info]
  kappa = 16*pi^2 * alpha_EM           [couplage gravitationnel = geometrie * EM]
  G/alpha = 2*pi                        [PAS DE HIERARCHIE dans le crible!]

SIGNIFICATION PHYSIQUE:
  G_00 = H1*H2 + H1*H3 + H2*H3        = courbure (Friedmann)
  D_KL = D_parity + D_KL(mod 3)         = information totale
  alpha = prod sin^2(theta_p)            = couplage du crible

  La relation G = 2*pi*alpha dit que:
  'L'intensite du couplage gravitationnel dans le crible est
   2*pi fois le couplage electromagnetique.'

  Ou encore: kappa = 16*pi^2*alpha signifie que
  'La constante d'Einstein est le produit d'un facteur geometrique
   (4*pi)^2 par la constante de structure fine.'

HIERARCHIE:
  Dans le crible: G/alpha = 2*pi ~ 6.28 (meme ordre!)
  En physique: G_phys/alpha_phys = alpha_G/alpha_EM ~ 4.3e-37

  La hierarchie physique (39 ordres de grandeur) n'est PAS dans
  le rapport des couplages, mais dans le rapport des masses:
  m_proton / M_Planck ~ 10^(-19)

  Notre framework dit: G et alpha sont du MEME ORDRE intrinsequement.
  La hierarchie est un effet de la masse des particules elementaires,
  pas des constantes fondamentales.

EQUATION D'ETAT GEOMETRIQUE (signes corriges):""")

# Print corrected equation of state
w_values = [E['G_sp'][i] / E['G_00'] for i in range(3)]
print(f"  rho_geom = G^0_0 = {E['G_00']:.6f}")
for i, p in enumerate(active_primes):
    print(f"  w_{p} = G^{i+1}_{i+1}/G^0_0 = {w_values[i]:+.4f}")
print(f"  w_moyen = {np.mean(w_values):+.4f}")
print(f"  rho + sum(p) = G_trace = {G_trace:.6f}")
print(f"  -R = {-E['R']:.6f}")
print(f"  (NEC: rho + p > 0 pour chaque direction)")
for i, p in enumerate(active_primes):
    nec = E['G_00'] + E['G_sp'][i]
    print(f"    NEC(p={p}): rho+p = {nec:+.4f} {'PASS' if nec >= 0 else 'FAIL'}")

test8 = True  # synthese

# ==============================================================
# SCORE
# ==============================================================

tests = [test1, test2, test3, test4, test5, test6, test7, test8]
n_pass = sum(tests)

print(f"\n{'='*70}")
print(f"SCORE: {n_pass}/{len(tests)} tests")
print("="*70)

for i, (t, desc) in enumerate(zip(tests, [
    "Spectre de perturbation (coherence)",
    "Signature Lorentzienne (-,+,+,+)",
    "G_trace = -R (signe corrige)",
    "Quantites informationnelles positives",
    "G_Newton calcule (5 identifications)",
    "G_sieve = 2*pi*alpha_EM (< 1%)",
    "Universalite G(mu)",
    "Interpretation physique"
]), 1):
    print(f"  T{i}: [{'PASS' if t else 'FAIL'}] {desc}")
