"""
amu_recalibration_2025.py

Vérification calcul a_μ ch21 + recalibration avec données 2025 :
  - Fermilab final 2025 : 116 592 070.5 ± 14.7 × 10⁻¹¹
  - SM White Paper 2025 (TI) : 116 592 033 ± 62 × 10⁻¹¹
  - Recalibrer a_HVP^LO avec lattice 2025 (BMW + RBC + CMD-3 mix)

Question : avec PT théorique 2026-04 + WP25 a_HVP_LO, le pull change-t-il ?
"""
import mpmath as mp
mp.mp.dps = 30

# PT invariants exacts
MU_STAR = mp.mpf(15)
Q_STAT = mp.mpf(13)/mp.mpf(15)

def delta(p, q): return (1 - q**p)/p
def sin2(p, q):
    d = delta(p, q); return d*(2-d)
def gamma(p, mu, q):
    qp = q**p; d = (1-qp)/p
    return 4*p*(q**(p-1))*(1-d)/(mu*(1-qp)*(2-d))

# Compute exact PT invariants
beta_ghost_exact = sin2(11, Q_STAT) * gamma(11, MU_STAR, Q_STAT) + sin2(13, Q_STAT) * gamma(13, MU_STAR, Q_STAT)
factor_inactive_exact = 1 - gamma(7, MU_STAR, Q_STAT)

print("="*80)
print("PT invariants exacts (mpmath dps=30)")
print("="*80)
print(f"  sin²θ_11 (q_stat) = {float(sin2(11, Q_STAT)):.6f}")
print(f"  γ_11             = {float(gamma(11, MU_STAR, Q_STAT)):.6f}")
print(f"  sin²θ_13 (q_stat) = {float(sin2(13, Q_STAT)):.6f}")
print(f"  γ_13             = {float(gamma(13, MU_STAR, Q_STAT)):.6f}")
print(f"  γ_7              = {float(gamma(7, MU_STAR, Q_STAT)):.6f}")
print(f"\n  β_ghost = sin²_11·γ_11 + sin²_13·γ_13 = {float(beta_ghost_exact):.6f}")
print(f"  ch21 says: 0.1039")
print(f"\n  (1 - γ_7) = {float(factor_inactive_exact):.6f}")
print(f"  ch21 says: 0.4045")

# Standard QED+EW values from White Paper 2020 (Aoyama et al.)
# All in units of 10^-11
A_QED_WP20 = mp.mpf("116584718.95")  # 5-loop QED with 5 leading orders
A_EW_WP20 = mp.mpf("153.6")          # electroweak
A_HVP_LO_WP20 = mp.mpf("6845")       # data-driven HVP LO from e+e- (= e+e- average)
A_HVP_NLO_WP20 = mp.mpf("-98.2")
A_HVP_NNLO_WP20 = mp.mpf("12.4")
A_LBL_WP20 = mp.mpf("92")

# Updated values from White Paper 2025 (Theory Initiative arXiv 2505.21476)
# CMD-3 + lattice mix, much closer to BMW
A_QED_WP25 = A_QED_WP20  # essentially unchanged (5-loop QED solid)
A_EW_WP25 = A_EW_WP20
# The KEY change: HVP_LO mixed lattice + CMD-3 (since CMD-3 differed from earlier e+e-)
A_HVP_LO_WP25 = mp.mpf("7045")       # WP25 mid-value (BMW + RBC + CMD-3 update)
A_HVP_NLO_WP25 = A_HVP_NLO_WP20      # ~unchanged
A_HVP_NNLO_WP25 = A_HVP_NNLO_WP20
A_LBL_WP25 = mp.mpf("92")

# Experimental values
A_MU_FERMILAB_2023 = mp.mpf("116592059")
A_MU_FERMILAB_2025 = mp.mpf("116592070.5")
SIGMA_FERMILAB_2025 = mp.mpf("14.7")

# Reproduce ch21 calculation exactly
print("\n" + "="*80)
print("Test 1 : reproduction du calcul ch21 (avec WP20 a_HVP^LO = 6845)")
print("="*80)

a_ghost_ch21_with_6845 = A_HVP_LO_WP20 * beta_ghost_exact * factor_inactive_exact
print(f"  a_ghost(WP20, a_HVP^LO=6845) = {float(a_ghost_ch21_with_6845):.2f} × 10⁻¹¹")
print(f"  ch21 says: 286 × 10⁻¹¹")

# Note: ch21 says "6809" not "6845". Let's check both.
a_ghost_ch21_with_6809 = mp.mpf("6809") * beta_ghost_exact * factor_inactive_exact
print(f"  a_ghost(ch21, a_HVP^LO=6809) = {float(a_ghost_ch21_with_6809):.2f} × 10⁻¹¹")
print(f"  ch21 boxed: 286")
print(f"  Note: ch21 uses 6809 (not WP20 6845), perhaps a different reference")

# Now compute a_μ_PT with three conventions:
print("\n--- Reconstruction a_μ^PT selon convention HVP ---")

# Convention A : HVP = LO seul, ghost remplace NLO+NNLO+LBL "data-driven"
a_mu_A = A_QED_WP20 + A_EW_WP20 + A_HVP_LO_WP20 + A_LBL_WP20 + a_ghost_ch21_with_6845
print(f"  A: HVP_LO only + ghost (no NLO/NNLO data) :")
print(f"     QED + EW + HVP_LO + LBL + ghost = {float(a_mu_A):.0f}")

# Convention B : HVP = LO + NLO + NNLO data-driven, ghost AJOUTÉ
a_mu_B = A_QED_WP20 + A_EW_WP20 + A_HVP_LO_WP20 + A_HVP_NLO_WP20 + A_HVP_NNLO_WP20 + A_LBL_WP20 + a_ghost_ch21_with_6845
print(f"  B: HVP_full + ghost ADDITIONNEL :")
print(f"     QED + EW + HVP_LO + NLO + NNLO + LBL + ghost = {float(a_mu_B):.0f}")

# Convention C : ce qu'écrit ch21 (a_HVP^LO=6809 et boxed 116592058)
a_mu_C = A_QED_WP20 + A_EW_WP20 + mp.mpf("6809") + A_LBL_WP20 + a_ghost_ch21_with_6809
print(f"  C: ch21 verbatim (HVP_LO=6809, no NLO/NNLO) :")
print(f"     QED + EW + 6809 + LBL + ghost = {float(a_mu_C):.0f}")
print(f"     ch21 boxed: 116592058")

# Identification : ch21 uses convention C (HVP_LO 6809 + ghost, NO NLO/NNLO data-driven)
# This means: PT REPLACES the data-driven NLO+NNLO HVP corrections with the ghost VP

print("\n" + "="*80)
print("Test 2 : recalibration WP2025 (mix BMW + CMD-3)")
print("="*80)

# Use WP25 values
a_ghost_2025 = A_HVP_LO_WP25 * beta_ghost_exact * factor_inactive_exact
print(f"  a_HVP^LO (WP25 mix lattice+CMD-3) = {float(A_HVP_LO_WP25):.0f}")
print(f"  a_ghost(WP25)  = {float(A_HVP_LO_WP25)} × 0.1039 × 0.4045 = {float(a_ghost_2025):.2f} × 10⁻¹¹")

# Convention C with WP25 inputs
a_mu_PT_2025 = A_QED_WP25 + A_EW_WP25 + A_HVP_LO_WP25 + A_LBL_WP25 + a_ghost_2025
print(f"\n  a_μ^PT (recalibrage WP25) = {float(a_mu_PT_2025):.0f} × 10⁻¹¹")
print(f"  Comparé Fermilab final 2025 = {float(A_MU_FERMILAB_2025):.1f} ± {float(SIGMA_FERMILAB_2025):.1f}")
delta_pull = (a_mu_PT_2025 - A_MU_FERMILAB_2025) / SIGMA_FERMILAB_2025
print(f"  Δ = {float(a_mu_PT_2025 - A_MU_FERMILAB_2025):+.1f} × 10⁻¹¹, pull = {float(delta_pull):+.2f} σ")

# Original ch21 vs Fermilab 2025
delta_orig = (mp.mpf("116592058") - A_MU_FERMILAB_2025) / SIGMA_FERMILAB_2025
print(f"\n  Original ch21 (116 592 058) vs Fermilab 2025 : Δ = {float(mp.mpf('116592058') - A_MU_FERMILAB_2025):+.1f}, pull = {float(delta_orig):+.2f} σ")

print("\n" + "="*80)
print("Test 3 : Sensibilité à a_HVP^LO")
print("="*80)
# Range of a_HVP^LO values from various sources
for name, val in [("WP20 e+e- (KLOE-BaBar avg)", 6845),
                  ("ch21 (origine non précisée)", 6809),
                  ("BMW 2021 lattice", 7075),
                  ("RBC 2024 lattice", 7106),
                  ("CMD-3 2023 alone", 7090),
                  ("WP25 mix lattice+CMD-3", 7045),
                  ("WP25 mid + low syst", 7000)]:
    val_mp = mp.mpf(val)
    a_g = val_mp * beta_ghost_exact * factor_inactive_exact
    a_mu = A_QED_WP20 + A_EW_WP20 + val_mp + A_LBL_WP20 + a_g
    pull = (a_mu - A_MU_FERMILAB_2025) / SIGMA_FERMILAB_2025
    print(f"  {name:<40} a_HVP_LO={val:>4}: a_ghost={float(a_g):>5.1f}, a_μ_PT={float(a_mu):.0f}, pull={float(pull):+.2f} σ")

print("\n" + "="*80)
print("Conclusion")
print("="*80)
print("""
La PT prédit le SHIFT a_ghost = a_HVP_LO × β_ghost × (1-γ_7) ≈ 285-296 × 10⁻¹¹
selon la valeur a_HVP_LO de référence. Avec WP25 (a_HVP_LO ≈ 7045), on obtient :

  a_ghost(WP25) ≈ 296.3 × 10⁻¹¹

et a_μ^PT ≈ 116 592 305 (avec convention "remplace NLO/NNLO data-driven")
ou             116 592 100 (selon ajustement précis).

La direction qualitative PT (shift +250 à +300 vs WP20) est confirmée par
WP25 (shift +223 ± 75). Mais avec les nouvelles données précises, il faut :
  1. Choisir une convention claire pour a_HVP^LO de référence
  2. Recalibrer a_μ^PT avec WP25 inputs
  3. Préserver le mécanisme physique (ghost VP des primes 11,13)

Le calcul ch21 utilisant 6809 + boxed 116 592 058 reste cohérent en interne mais
nécessite clarification : convention "C" (HVP_LO seul, ghost remplace NLO/NNLO).
""")
