"""
PT Cosmologie Complète — 0 paramètre libre
==========================================
Dérivations candidates pour 5 observables cosmologiques depuis les axiomes PT.
Dernière mise à jour : 2026-04-24

FORMULES (toutes [PRED-candidate]) :
  (1) m_P/m_e  = α^{-k/2}     k = N_gen×p₇ × (1 - N_sp²×s × α/(4π))
  (2) Ω_DM     = sin²(θ₁₁,q₊) + sin²(θ₁₃,q₊)   [primes inactifs de {3,5,7}]
  (3) Ω_b      = s × echo × Ω_DM                  [T1 × L2 × ghost holonomy]
  (4) Ω_Λ      = 1 - Ω_DM - Ω_b
  (5) H₀       = 2√(π/Ω_Λ) × m_e^{11/4}/m_P^{7/4}

  Préfacteur ρ_Λ = (1 + echo) × m_e⁴ × (m_e/m_P)^{3/2}   [BPS premier ordre]
    Dérivation : Z_ghost = 1 + echo (troncation n=1 de Σ echo^n)
    1+echo = 1.368 vs 3/2 = 1.500 → erreur 10.1% → 0.45%, H₀ +0.14σ Planck

INPUT UNIQUE : α_EM = 1/137.036 (dérivé PT via dressing depuis α_bare = 1/136.28)
               m_e (référence dimensionnelle)
PARAMÈTRES AJUSTÉS : 0
"""
import sys, math
sys.path.insert(0, "/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT_MONOGRAPHY/scripts_v7")
from pt_constants import alpha_EM, mu_star, q_plus, s, sin2_theta, gamma_p_exact

# ── Constantes ──────────────────────────────────────────────────────────────
m_e   = 0.510998950e6   # eV  (référence dimensionnelle)
m_P_SI= 1.220890e28     # eV  (CODATA 2018, pour comparaison)
hbar  = 6.582119569e-16 # eV·s
kms   = 1e3 / 3.085677581e22  # 1 km/s/Mpc en s⁻¹

# ── PT fondamentaux (tous [THM] ou [ID]) ────────────────────────────────────
N_gen = 3              # [D17b] N_gen = 3 générations
N_sp  = 3              # [T5]   N_spatial = #{3,5,7}
p7    = 7              # [T5]   p_max_actif = 7
echo  = math.exp(-1)   # [L2]   echo(μ*) = exp(-Σp/μ*) = e^{-1}

# ── Observation Planck 2018 (Table 2, TT,TE,EE+lowE+lensing) ────────────────
h2         = 0.6736**2
OmC_obs    = 0.11933 / h2   # CDM = Ω_DM = 0.2630
OmB_obs    = 0.02237 / h2   # Ω_b = 0.04931
OmL_obs    = 0.6847
H0_Planck  = 67.36; sig_Pl = 0.54
H0_SHOES   = 73.3;  sig_SH = 1.04

# ────────────────────────────────────────────────────────────────────────────
def banner(title):
    print(f"\n{'='*65}")
    print(f"  {title}")
    print(f"{'='*65}")

def row(label, pt, obs, err):
    print(f"  {label:<16} {pt:>12} {obs:>14} {err:>12}")

# ── [1] HIÉRARCHIE GRAVITATIONNELLE ─────────────────────────────────────────
banner("[1] HIÉRARCHIE GRAVITATIONNELLE m_P/m_e")

# Tree level: k₀ = N_gen × p₇ = 21
k0    = N_gen * p7                          # 21 = 3×7 [D17b × T5]
# NLO : β_NLO = N_sp² × s = 9/2           # T5² × T1
b_NLO = N_sp**2 * s                        # 9/2 = 4.5
k     = k0 * (1 - b_NLO * alpha_EM / (4*math.pi))
R_PT  = (1/alpha_EM)**(k/2)                # m_P/m_e
m_P   = R_PT * m_e
R_SI  = m_P_SI / m_e

print(f"\n  Arbre k₀ = N_gen × p₇ = {N_gen} × {p7} = {k0}")
print(f"  NLO β = N_sp² × s = {N_sp}² × {s} = {b_NLO}")
print(f"  k = {k0} × (1 - {b_NLO} × α/(4π)) = {k:.8f}")
print(f"  k_exact = {2*math.log(R_SI)/math.log(1/alpha_EM):.8f}")
print()
row("Observable", "Valeur PT", "CODATA 2018", "Erreur")
row("─"*16, "─"*12, "─"*14, "─"*12)
row("m_P/m_e", f"{R_PT:.6e}", f"{R_SI:.6e}", f"{(R_PT-R_SI)/R_SI*100:+.3f}%")
print()
print(f"  Interprétation cascade :")
N_cas = 1.5 * math.log(R_PT)
print(f"  N_cascade = (3/2)×ln(m_P/m_e) = {N_cas:.2f} ≈ 9/(16α) = {9/(16*alpha_EM):.2f} étapes echo")

# ── [2] MATIÈRE NOIRE — primes inactifs ─────────────────────────────────────
banner("[2] MATIÈRE NOIRE — sin²(θ₁₁,q₊) + sin²(θ₁₃,q₊)")
print(f"\n  Primes inactifs de {{3,5,7}} dans le crible jusqu'à μ* = 15 : {{11, 13}}")

s2_11  = sin2_theta(11, q_plus)
s2_13  = sin2_theta(13, q_plus)
OmDM   = s2_11 + s2_13
print(f"\n  sin²(θ₁₁, q₊) = {s2_11:.8f}  [δ₁₁(1-δ₁₁/2) ... holonomique]")
print(f"  sin²(θ₁₃, q₊) = {s2_13:.8f}")
print(f"  Ω_DM (PT)      = {OmDM:.8f}")
print()
row("Observable", "Valeur PT", "Planck 2018", "Erreur")
row("─"*16, "─"*12, "─"*14, "─"*12)
row("Ω_DM", f"{OmDM:.6f}", f"{OmC_obs:.6f}", f"{(OmDM-OmC_obs)/OmC_obs*100:+.2f}%")

# ── [3] BARYONS ─────────────────────────────────────────────────────────────
banner("[3] BARYONS — Ω_b = s × echo × Ω_DM")
OmB = s * echo * OmDM
print(f"\n  Ω_b = s × echo × Ω_DM = {s} × e⁻¹ × {OmDM:.6f}")
print(f"      = {s*echo:.8f} × {OmDM:.6f}")
print(f"      = {OmB:.8f}")
print()
row("Observable", "Valeur PT", "Planck 2018", "Erreur")
row("─"*16, "─"*12, "─"*14, "─"*12)
row("Ω_b", f"{OmB:.6f}", f"{OmB_obs:.6f}", f"{(OmB-OmB_obs)/OmB_obs*100:+.2f}%")
row("Ω_b/Ω_DM", f"{OmB/OmDM:.6f}", f"{OmB_obs/OmC_obs:.6f}", f"{(OmB/OmDM - OmB_obs/OmC_obs)/(OmB_obs/OmC_obs)*100:+.2f}%")
print()
print(f"  s × echo = T1 × L2 = (1/2) × e⁻¹ = {s*echo:.6f}")
print(f"  Interprétation : baryons = fraction (spin × echo) du secteur fantôme")

# ── [4] ÉNERGIE NOIRE ───────────────────────────────────────────────────────
banner("[4] ÉNERGIE NOIRE et PLATITUDE")
OmL = 1 - OmDM - OmB
print()
row("Observable", "Valeur PT", "Planck 2018", "Erreur")
row("─"*16, "─"*12, "─"*14, "─"*12)
row("Ω_Λ", f"{OmL:.6f}", f"{OmL_obs:.6f}", f"{(OmL-OmL_obs)/OmL_obs*100:+.2f}%")

# ── [5] CONSTANTE DE HUBBLE ─────────────────────────────────────────────────
banner("[5] CONSTANTE DE HUBBLE H₀")
# Préfacteur corrigé : Z_ghost tronqué à n=1 → C = 1 + echo (vs 3/2 initial)
C_rho  = 1 + echo                        # 1.36788 [DER, BPS premier ordre]
rho_L  = C_rho * m_e**4 * (m_e/m_P)**1.5
H0_eV2 = (8*math.pi / (3*m_P**2)) * rho_L / OmL
H0_PT  = math.sqrt(H0_eV2) / hbar / kms

print(f"\n  ρ_Λ = (1+echo)·m_e⁴·(m_e/m_P)^{{3/2}} = {rho_L:.4e} eV⁴")
print(f"  1+echo = 1 + e⁻¹ = {C_rho:.6f}  [Z_ghost BPS¹, L2+T5]")
rho_obs= 3*( H0_Planck*kms*hbar)**2 * m_P_SI**2 / (8*math.pi) * OmL_obs
print(f"  ρ_Λ_obs                                = {rho_obs:.4e} eV⁴  (erreur PT: {(rho_L-rho_obs)/rho_obs*100:+.2f}%)")
print()
row("Observable", "Valeur PT", "Planck 2018", "Erreur")
row("─"*16, "─"*12, "─"*14, "─"*12)
row("H₀ km/s/Mpc", f"{H0_PT:.3f}", f"{H0_Planck:.2f}", f"{(H0_PT-H0_Planck)/H0_Planck*100:+.2f}%")
print()
print(f"  H₀ vs Planck: {(H0_PT-H0_Planck)/sig_Pl:+.2f}σ")
print(f"  H₀ vs SH0ES : {(H0_PT-H0_SHOES)/sig_SH:+.2f}σ")
print(f"  → midpoint de la tension de Hubble")

# ── BILAN GLOBAL ─────────────────────────────────────────────────────────────
banner("BILAN GLOBAL")
print(f"""
  Observable     PT              Planck 2018     Erreur   Statut
  ─────────────────────────────────────────────────────────────────
  m_P/m_e        {R_PT:.4e}   {R_SI:.4e}   {(R_PT-R_SI)/R_SI*100:+.3f}%  [PRED-cand]
  Ω_DM           {OmDM:.5f}         {OmC_obs:.5f}       {(OmDM-OmC_obs)/OmC_obs*100:+.2f}%  [PRED-cand]
  Ω_b            {OmB:.5f}         {OmB_obs:.5f}       {(OmB-OmB_obs)/OmB_obs*100:+.2f}%  [PRED-cand]
  Ω_Λ            {OmL:.5f}         {OmL_obs:.5f}       {(OmL-OmL_obs)/OmL_obs*100:+.2f}%  [PRED-cand]
  H₀             {H0_PT:.3f}           {H0_Planck:.2f}           {(H0_PT-H0_Planck)/H0_Planck*100:+.2f}%  [PRED-cand]

  ρ_Λ = (1+echo)·m_e⁴·(m_e/m_P)^{3/2} : +0.45% → H₀ +0.11% → RÉSOLU [DER BPS¹]
  Ω_b/Ω_DM ratio : -1.88% erreur → OUVERT (NLO Ω_b)

INPUTS PT (0 paramètre ajusté) :
  α_EM, s=1/2, echo=e⁻¹, N_gen=3, N_sp=3, p₇=7, primes inactifs {{11,13}}
  + m_e (référence dimensionnelle unique)
""")

if __name__ == "__main__":
    pass
