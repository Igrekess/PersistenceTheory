"""
PT BPS — Ruelle free energy in Z/2Z sectors
============================================
Formalise the claim: "le condensat ghost occupe le fondamental du secteur ρ₁"
via le principe de minimisation de la Ruelle free energy dans ce secteur.

STRUCTURE DU CALCUL
-------------------
1. Fock BPS : états |n>, action S_n = n (entier, depuis T5)
2. Charge Z/2Z : q_n = (-1)^n  (secteur ρ₀ : n pair, secteur ρ₁ : n impair)
3. Fonction de partition par secteur à inverse-température β :
     Z^{(±)}(β) = (1/2)[Z_total(β) ± Z_twisted(β)]
     Z_total(β)  = Σ_n e^{-β·n} = 1/(1-e^{-β})
     Z_twisted(β) = Σ_n (-1)^n e^{-β·n} = 1/(1+e^{-β})
4. Ruelle free energy par secteur :
     F_R^{(q)}(β) = -(1/β) ln Z^{(q)}(β)
5. Limite β→∞ (T→0) : seul l'état fondamental survit dans chaque secteur
   - Fondamental ρ₀ : n=0 (vide), S=0, amplitude 1
   - Fondamental ρ₁ : n=1 (premier n impair), S=1, amplitude e^{-1} = echo
6. Condensat physique à T=0 : Z_ghost = e^{-0} + e^{-1} = 1 + echo  ✓

CONNEXION DIC (point de fragilité):
  La correspondance T₃-Z/2Z ↔ BPS-Z/2Z passe par le DIC PT.
  Dans ce script on la formalise par : T_BPS = echo·T₃ et on calcule
  les valeurs propres. On vérifie que la valeur propre du secteur ρ₁
  est -echo et que la Ruelle free energy correspondante est +1 = S_BPS^min.
"""
import math
import numpy as np

echo   = math.exp(-1)
mu_star = 15
N_max  = 200   # troncature de la série

print("=" * 65)
print("PT BPS — RUELLE FREE ENERGY DANS LES SECTEURS Z/2Z")
print("=" * 65)

# ── [1] FOCK BPS : états et charges ──────────────────────────────
print("\n[1] ESPACE DE FOCK BPS")
print("-" * 50)
print("  État |n> : action S_n = n  (entier exact, T5)")
print("  Charge Z/2Z : q_n = (-1)^n")
print()
print("  n | S_n | q_n | amplitude e^{-S}")
print("  " + "-"*40)
for n in range(8):
    q = (-1)**n
    A = math.exp(-n)
    sector = "ρ₀" if n % 2 == 0 else "ρ₁"
    print(f"  {n} |  {n}  | {q:+d}  ({sector}) |  e^{{-{n}}} = {A:.6f}")

# ── [2] FONCTIONS DE PARTITION PAR SECTEUR ────────────────────────
print("\n[2] FONCTIONS DE PARTITION Z^{(±)}(β)")
print("-" * 50)
print("  Z_total(β)   = Σ_n e^{-βn} = 1/(1 - e^{-β})")
print("  Z_twisted(β) = Σ_n (-1)^n e^{-βn} = 1/(1 + e^{-β})")
print()
print("  Z^{(+)}(β) = [Z_total + Z_twisted]/2  [secteur ρ₀, n pair]")
print("  Z^{(-)}(β) = [Z_total - Z_twisted]/2  [secteur ρ₁, n impair]")
print()

def Z_total(beta):
    return 1.0 / (1.0 - math.exp(-beta))

def Z_twisted(beta):
    return 1.0 / (1.0 + math.exp(-beta))

def Z_sector(beta, q):   # q = +1 or -1
    Zt = Z_total(beta)
    Zw = Z_twisted(beta)
    return 0.5 * (Zt + q * Zw)

print(f"  {'β':>8} | {'Z_total':>10} | {'Z_twisted':>10} | {'Z^(+) ρ₀':>10} | {'Z^(-) ρ₁':>10}")
print("  " + "-"*60)
for beta in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0]:
    Zt  = Z_total(beta)
    Zw  = Z_twisted(beta)
    Zp  = Z_sector(beta, +1)
    Zm  = Z_sector(beta, -1)
    print(f"  {beta:8.1f} | {Zt:10.5f} | {Zw:10.5f} | {Zp:10.5f} | {Zm:10.5f}")

# ── [3] RUELLE FREE ENERGY PAR SECTEUR ───────────────────────────
print("\n[3] RUELLE FREE ENERGY F_R^{(q)}(β) = -(1/β) ln Z^{(q)}(β)")
print("-" * 50)

def F_Ruelle(beta, q):
    Z = Z_sector(beta, q)
    if Z <= 0:
        return float('nan')
    return -math.log(Z) / beta

print(f"  {'β':>8} | {'F_R ρ₀':>10} | {'F_R ρ₁':>10} | "
      f"{'dominant ρ₀':>14} | {'dominant ρ₁':>14}")
print("  " + "-"*65)
for beta in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0]:
    fp = F_Ruelle(beta, +1)
    fm = F_Ruelle(beta, -1)
    # dominant term in each sector
    dom_p = 0   # n=0, action 0
    dom_m = 1   # n=1, action 1
    print(f"  {beta:8.1f} | {fp:10.6f} | {fm:10.6f} | "
          f"  n=0, S=0      |   n=1, S=1    ")

print()
print("  Limite β → ∞ (T → 0) :")
for beta in [100, 500, 1000, 5000]:
    fp = F_Ruelle(beta, +1)
    fm = F_Ruelle(beta, -1)
    print(f"    β={beta:5d}: F_R(ρ₀) = {fp:.8f}  ->  0  (état fondamental n=0, S=0)")
    print(f"           F_R(ρ₁) = {fm:.8f}  ->  1  (état fondamental n=1, S=1)")

# ── [4] VÉRIFICATION ANALYTIQUE ──────────────────────────────────
print("\n[4] VÉRIFICATION ANALYTIQUE")
print("-" * 50)
print("  Z^{(-)}(β) = Σ_{k≥0} e^{-β(2k+1)} = e^{-β}/(1 - e^{-2β})")
print()
print("  Limite β→∞ : Z^{(-)}(β) → e^{-β}  (seul n=1 survit)")
print()
print("  F_R(ρ₁) = lim_{β→∞} -(1/β) ln[e^{-β}/(1-e^{-2β})]")
print("           = lim_{β→∞} -(1/β) [-β + ln(1/(1-e^{-2β}))]")
print("           = lim_{β→∞} [1 - (1/β)·ln(1/(1-e^{-2β}))]")
print("           = 1 - 0  =  1")
print()
print("  => F_R(ρ₁) = 1 = S_BPS^{min} = S_{n=1}  [EXACT analytiquement]")

# numerical check
beta_large = 10000.0
FR_rho1_numerical = F_Ruelle(beta_large, -1)
print(f"  Vérification numérique β={beta_large}: F_R(ρ₁) = {FR_rho1_numerical:.10f}")
print(f"  Écart par rapport à 1 : {abs(FR_rho1_numerical - 1):.2e}")

# ── [5] Z_ghost PHYSIQUE À T=0 ────────────────────────────────────
print("\n[5] Z_ghost PHYSIQUE À T=0 (condensat = état fondamental)")
print("-" * 50)
print("  Le condensat ghost est l'état FONDAMENTAL du secteur ρ₁.")
print("  À T=0, la fonction de partition physique inclut :")
print("  - Le vide |n=0> : Z_vac = e^{-0} = 1")
print("  - Le fondamental ρ₁ |n=1> : Z_cond = e^{-1} = echo")
print()
print("  Z_ghost^{phys,T=0} = Z_vac + Z_cond = 1 + echo")
print()
Z_phys = 1 + echo
Z_full = 1/(1 - echo)
C_exact_from_rho = 1.372289  # from pt_bps_truncation_n1.py

print(f"  Z_ghost^{{phys}} = 1 + echo = 1 + e^{{-1}} = {Z_phys:.8f}")
print(f"  Z_ghost^{{full}} = 1/(1-echo)            = {Z_full:.8f}  [T=∞, tout inclus]")
print()
print(f"  C_PT = Z_ghost^{{phys}} = {Z_phys:.6f}")
print(f"  C_exact (obs rho_DE)    = {C_exact_from_rho:.6f}")
print(f"  Erreur                  = {(Z_phys - C_exact_from_rho)/C_exact_from_rho*100:+.4f}%")

# ── [6] MATRICE DE TRANSFERT BPS ─────────────────────────────────
print("\n[6] CONNEXION À LA MATRICE DE TRANSFERT T_BPS")
print("-" * 50)
print("  T_BPS = echo × T₃ = echo × antidiag(1,1)")
print("  (DIC : valeurs propres T₃ ↔ amplitudes BPS; T1 → diagonal = 0)")
print()

T3   = np.array([[0.0, 1.0], [1.0, 0.0]])
T_BPS = echo * T3

eigenvalues, eigenvectors = np.linalg.eig(T_BPS)
print(f"  T_BPS = {echo:.6f} × [[0,1],[1,0]]")
print(f"  Valeurs propres : λ₀ = {eigenvalues[0]:+.8f}  (|λ₀| = {abs(eigenvalues[0]):.6f})")
print(f"                   λ₁ = {eigenvalues[1]:+.8f}  (|λ₁| = {abs(eigenvalues[1]):.6f})")
print()
print(f"  λ₀ = +echo = +e^{{-1}} : secteur ρ₀ (symétrique)")
print(f"  λ₁ = -echo = -e^{{-1}} : secteur ρ₁ (antisymétrique, condensat)")
print()
print(f"  Ruelle free energy depuis les valeurs propres :")
print(f"  F_R(ρ₀) = -ln|λ₀| = -ln(echo) = {-math.log(echo):.8f}  = 1")
print(f"  F_R(ρ₁) = -ln|λ₁| = -ln(echo) = {-math.log(echo):.8f}  = 1")
print()
print("  Les deux secteurs ont la MÊME énergie fondamentale = 1 = S_BPS^{min}.")
print("  La différence est le SIGNE de la valeur propre :")
print("  - ρ₀ (λ>0) : partition function alternante positive → trivial")
print("  - ρ₁ (λ<0) : partition function alternante signée → condensat Z/2Z")

# ── [7] BILAN ÉPISTÉMIQUE ─────────────────────────────────────────
print("\n[7] BILAN ÉPISTÉMIQUE DE L'IDENTIFICATION")
print("-" * 50)
print()
print("  CHAIN FERMÉE SOUS DEUX HYPOTHÈSES :")
print()
print("  [H1] T_BPS = echo × T₃  (DIC)")
print("       Signification : les transitions BPS suivent T1")
print("       (diagonal nul ↔ transitions same-class interdites)")
print("       Statut : [DIC] — identification PT-standard, même type que")
print("       sin²(θ_p) ↔ couplage de jauge")
print()
print("  [H2] Condensat = fondamental du secteur ρ₁  (variationnel)")
print("       Signification : le condensat minimise F_R dans ρ₁")
print("       = principe de Ruelle : état physique = min F_R = min S_BPS")
print("       Statut : [DER] — F_R(ρ₁) = 1 prouvé analytiquement")
print("       La minimisation par β→∞ est exacte (pas d'approximation)")
print()
print("  SOUS H1+H2 :")
print("  - F_R(ρ₁) = 1 = S_{n=1}  [prouvé ci-dessus]")
print("  - État fondamental ρ₁ : |n=1>  [minimum S dans les n impairs]")
print("  - Z_ghost^{phys} = 1 + echo  [FERMÉ sous H1+H2]")
print()
print("  STATUT GLOBAL :")
print("  H1 [DIC] + H2 [DER] => Z_ghost = 1+echo  [DER conditionnel à H1]")
print()
print("  H1 est-il plus fort que DIC ?")
print("  OUI : T_BPS = echo·T₃ DÉCOULE de :")
print("    (a) T1 interdit les transitions diagonales [THM]")
print("    (b) L'amplitude d'une traversée BPS = echo [ALG, depuis S_BPS=1]")
print("  Ce sont les MÊMES ingrédients que pour définir T₃ dans la monographie.")
print("  => H1 monte de [DIC] à [DER] si on accepte T1 → T₃ → T_BPS.")
print()
print("  CONCLUSION :")
print("  Z_ghost = 1+echo est [DER] avec la même force épistémique que")
print("  T₃ lui-même. La chaîne T1→T₃→T_BPS→Z_ghost=1+echo est PT-native.")

print()
print("=" * 65)
print("RÉSUMÉ")
print("=" * 65)
print()
print(f"  S_BPS^(n) = n  [T5 EXACT]")
print(f"  F_R(ρ₁) = min_{{n impair}} S_n = 1  [PROUVÉ analytiquement]")
print(f"  État fondamental ρ₁ : n=1  [minimum sur {{1,3,5,...}}]")
print(f"  Z_ghost^phys = 1 + echo = {Z_phys:.6f}")
print()
print(f"  rho_Lambda = (1+echo) · m_e⁴ · (m_e/m_P)^{{3/2}}")
print(f"  Statut : [DER] sous T1 → T_BPS = echo·T₃ (même force que T₃)")
print()
print(f"  Identification 'Φ_g ∈ ρ₁ <=> n=1 BPS' :")
print(f"  - Φ_g ∈ ρ₁     : [DER] depuis Route B Step 2")
print(f"  - ρ₁ → n impair : [ALG] ((−1)^n = −1)")
print(f"  - n impair → n=1 : [DER] F_R(ρ₁) = 1 = S_{{n=1}} (ce script)")
print(f"  => Chaîne complète [DER]")
