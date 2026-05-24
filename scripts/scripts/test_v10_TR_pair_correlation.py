"""
V10 — Test : Pair correlation des zéros γ_n de ζ vs Bergman kernel
       ω_{0,2}^pers de la courbe spectrale Σ_pers via Eynard-Orantin.

HYPOTHÈSE : la pair correlation des zéros (Montgomery 1973) suit le
sine-kernel GUE :
   K(x) = 1 - [sin(πx)/(πx)]²

Question PT : ω_{0,2}^pers évalué sur les paires (γ_n, γ_m) reproduit-il
ce kernel ? Si OUI = nouvelle correspondance PT ↔ Riemann (Cas 2 majeur).
Si NON = falsification disciplinée (Cas 4).

MÉTHODOLOGIE RIGOUREUSE :

(1) Unfolding (désinflation) des zéros :
    On définit γ̃_n = N(γ_n) ≈ (γ_n / 2π) [log(γ_n/2π) - 1] + 7/8
    Pour ces γ̃_n, l'espacement local est exactement 1 (par construction).

(2) R_2(x) empirique est définie comme :
    R_2(x) = lim_{N→∞} (1/N) Σ_{i≠j, i,j ≤ N} δ(x - (γ̃_i - γ̃_j))
    Plus précisément, pour un bin [a, b] :
    nombre de paires (i, j), i < j, avec x_ij = γ̃_j - γ̃_i ∈ [a, b]
    et puis on divise par (N × bin_width).

(3) Le pair correlation de Montgomery est :
    R_2(x) = 1 - [sin(πx)/(πx)]² + δ(x)
    Pour x > 0 (paires distinctes), on a R_2(x) = 1 - [sin(πx)/(πx)]²

(4) Comparaison avec :
    (a) Montgomery sine-kernel
    (b) Bergman kernel PT-modifié via substitution Mirzakhani

(5) Test ratio sur 2 décades : γ_n ∈ [10², 10³] et γ_n ∈ [10³, 10⁴]

DISCIPLINE :
- Test ratio sur 2 décades obligatoire (mémoire feedback)
- Distinguer "match" vs "coïncidence statistique"
- Si match suspect, RÉTRACTER (V6, V9, M1bis)
"""

import mpmath as mp
import numpy as np
from pathlib import Path
import json

mp.mp.dps = 30

# ============================================================================
# PARAMETRES
# ============================================================================

# Echelle 1: γ_n petits (γ ~ 10-500)
SCALE_LOW_MAX_N = 500

# Echelle 2: γ_n autour de γ ~ 1500-2500
SCALE_HIGH_START_N = 1000
SCALE_HIGH_END_N = 1700

# Plage de x pour pair correlation
X_MIN = 0.0
X_MAX = 3.0
N_BINS = 30  # 0.1 par bin pour avoir suffisamment de paires/bin

# PT canonical : primes actifs et leurs gamma_p
PRIMES = [3, 5, 7]
MU_STAR = mp.mpf(15)


# ============================================================================
# 1. CALCUL γ_p(μ*) — POUR LA PRÉDICTION PT
# ============================================================================

def q_plus(mu):
    """q_+(μ) = 1 - 2/μ"""
    return mp.mpf(1) - mp.mpf(2) / mp.mpf(mu)

def gamma_p_exact(p, mu):
    """γ_p(μ) = dimension anomale du canal p à l'échelle μ.

    Formule T6 : γ_p = 4p · q^(p-1) · (1-δ_p) / (μ · (1-q^p) · (2-δ_p))
    avec q = q_+(μ), δ_p = (1-q^p)/p.
    """
    p = mp.mpf(p)
    mu = mp.mpf(mu)
    q = q_plus(mu)
    qp = q**p
    one_minus_qp = mp.mpf(1) - qp
    delta = one_minus_qp / p
    numerator = 4 * p * (q**(p-1)) * (mp.mpf(1) - delta)
    denominator = mu * one_minus_qp * (mp.mpf(2) - delta)
    return numerator / denominator


GAMMA_P_AT_15 = {p: gamma_p_exact(p, MU_STAR) for p in PRIMES}
SUM_GAMMA = sum(GAMMA_P_AT_15.values())

print("=" * 78)
print(" V10 — TEST PAIR CORRELATION γ_n vs BERGMAN KERNEL ω_{0,2}^pers")
print("=" * 78)
print()
print(f"PT canonical (μ*=15, primes actifs {{3,5,7}}) :")
for p in PRIMES:
    print(f"  γ_{p} = {mp.nstr(GAMMA_P_AT_15[p], 12)}")
print(f"  Σ γ_p = {mp.nstr(SUM_GAMMA, 12)}")
print(f"  (1/2) Σ γ_p = {mp.nstr(SUM_GAMMA/2, 12)}")
print()
print(f"Comparaison transcendante :")
print(f"  π² = {np.pi**2:.6f}")
print(f"  2π²/3 = {2*np.pi**2/3:.6f}  [coefficient Mirzakhani]")
print(f"  Substitution Phase 2b : π²/3 → (1/4) Σγ_p, donc π² → (3/4) Σγ_p = {0.75*float(SUM_GAMMA):.6f}")
print(f"  Rapport (3/4 Σγ_p) / π² = {0.75*float(SUM_GAMMA)/(np.pi**2):.6f}")
print()


# ============================================================================
# 2. CHARGE LES ZÉROS DE RIEMANN
# ============================================================================

def load_zeros(n_start, n_end):
    """Charge γ_n pour n ∈ [n_start, n_end]."""
    n_count = n_end - n_start + 1
    print(f"  Chargement {n_count} γ_n (n ∈ [{n_start}, {n_end}])...")
    zeros = np.zeros(n_count)
    for i, n in enumerate(range(n_start, n_end + 1)):
        zn = mp.im(mp.zetazero(n))
        zeros[i] = float(zn)
        if (n - n_start) % 100 == 0 and n > n_start:
            print(f"    ... n={n}")
    return zeros


def unfold_zeros(gammas):
    """Unfolding par N(γ) = Riemann-von Mangoldt.

    N(γ) = (γ/2π) [log(γ/2π) - 1] + 7/8 + O(1/γ)
         = (γ/2π) log(γ/(2π e)) + 7/8

    Après unfolding, l'espacement moyen est 1 (par construction).
    """
    twopi = 2 * np.pi
    return (gammas / twopi) * np.log(gammas / (twopi * np.e)) + 7.0/8.0


# ============================================================================
# 3. KERNEL DE MONTGOMERY (GUE sine kernel)
# ============================================================================

def montgomery_kernel(x):
    """K_Mont(x) = 1 - [sin(πx) / (πx)]² (Montgomery 1973).

    Pair correlation des valeurs propres GUE (Gaussian Unitary Ensemble).
    Conjecturée pour les zéros γ̃_n unfolded de ζ.
    """
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    mask = x > 1e-12
    px = np.pi * x[mask]
    out[mask] = 1.0 - (np.sin(px) / px)**2
    # K(0) = 0 mathématiquement, mais en pratique on évite le bin x=0
    return out


# ============================================================================
# 4. KERNEL PT-MODIFIÉ (substitution Mirzakhani)
# ============================================================================

def pt_kernel(x, scheme='direct'):
    """Kernel PT via substitution Mirzakhani π²/3 → (1/4) Σγ_p.

    Schémas :
    - 'direct'   : K_PT = K_Mont (la substitution ne touche pas le sin-kernel)
    - 'pt_omega' : Substitution la fréquence π → ω_PT = √(3 Σγ_p / 4)
    - 'pt_factor': Substitution l'amplitude (3 Σγ_p / 4) / π²

    NOTE : Ces schémas sont des HYPOTHÈSES NON-DÉRIVÉES de TR rigoureuse.
    Le calcul exact de ω_{0,2}^pers pour Σ_pers genus 14 demanderait les
    Riemann theta functions, infaisable analytiquement. On teste donc des
    formes plausibles et on accepte la falsification si elles ne marchent pas.
    """
    x = np.asarray(x, dtype=float)
    sum_g = float(SUM_GAMMA)

    if scheme == 'direct':
        return montgomery_kernel(x)
    elif scheme == 'pt_omega':
        # Substitution de la fréquence : π → √((3/4) Σγ_p)
        omega_pt = np.sqrt(0.75 * sum_g)
        out = np.zeros_like(x)
        mask = x > 1e-12
        px = omega_pt * x[mask]
        out[mask] = 1.0 - (np.sin(px) / px)**2
        return out
    elif scheme == 'pt_factor':
        # Substitution de l'amplitude du terme oscillant
        # K_PT = 1 - [(3 Σγ_p / 4π²)] * [sin(πx)/(πx)]²
        ratio = 0.75 * sum_g / (np.pi**2)
        out = np.zeros_like(x)
        mask = x > 1e-12
        px = np.pi * x[mask]
        out[mask] = 1.0 - ratio * (np.sin(px) / px)**2
        return out
    else:
        raise ValueError(f"Unknown scheme: {scheme}")


# ============================================================================
# 5. PAIR CORRELATION EMPIRIQUE — VERSION RIGOUREUSE
# ============================================================================

def empirical_pair_correlation(gammas_unfolded, x_min=X_MIN, x_max=X_MAX,
                                n_bins=N_BINS, window_label="", verbose=True):
    """Calcule R_2(x) pour zéros unfolded γ̃_n.

    Pour des zéros unfolded (espacement local = 1), R_2(x) est défini par :
        R_2(x) = (1/N) · (nombre de paires (i,j), i ≠ j, |γ̃_i - γ̃_j| ∈ [x-Δ/2, x+Δ/2])
                  / Δ
    Pour x > 0, en prenant les paires ordonnées (i < j),
        R_2(x) = (2/N) · #{(i,j) : i < j, γ̃_j - γ̃_i ∈ bin centré sur x} / Δ

    L'espérance pour un processus de Poisson (zéros décorrélés) est R_2 = 1.

    Pour Riemann zeros, Montgomery conjecture
        R_2(x) → 1 - [sin(πx)/(πx)]² quand N → ∞.

    Args:
        gammas_unfolded : array des γ̃_n unfolded
        x_min, x_max : fenêtre de comparaison
        n_bins : nombre de bins

    Returns:
        x_centers : centres de bins
        R2_emp : pair correlation empirique
        counts : nombre brut de paires par bin
    """
    N = len(gammas_unfolded)
    bins = np.linspace(x_min, x_max, n_bins + 1)
    bin_width = bins[1] - bins[0]
    centers = (bins[:-1] + bins[1:]) / 2

    # Toutes les distances ordonnées x_ij = γ̃_j - γ̃_i pour i < j
    # ATTENTION : np.subtract.outer(a, b)[i,j] = a[i] - b[j]
    # Pour i < j, on a a[i] < a[j], donc subtract.outer[i,j] < 0
    # On veut |γ̃_i - γ̃_j| = γ̃_j - γ̃_i = -(a[i] - b[j])
    deltas = np.subtract.outer(gammas_unfolded, gammas_unfolded)
    iu = np.triu_indices(N, k=1)
    pos_deltas = -deltas[iu]  # tous positifs maintenant (γ̃_j - γ̃_i pour i < j)

    # Histogramme dans la fenêtre
    in_window = pos_deltas[(pos_deltas >= x_min) & (pos_deltas <= x_max)]
    counts, _ = np.histogram(in_window, bins=bins)

    # Normalisation rigoureuse :
    # R_2(x) dx = E[# paires (i,j) avec γ̃_j - γ̃_i ∈ [x, x+dx]] / N (paires ordonnées)
    # E[# paires ordonnées ∈ bin] = (N/2) · R_2(x_center) · bin_width
    # (à grande N pour x < N)
    # Donc R_2 ≈ 2 · counts / (N · bin_width)
    # MAIS : pour les zéros de Riemann, on est dans une fenêtre finie. La
    # densité moyenne est 1, donc dans [0, X] la moyenne de paires
    # par bin de largeur Δ est : N · Δ pour x petit comparé à N.
    # Donc R_2 = counts / (N · bin_width) si on compte les paires (i, j)
    # avec j > i. Il faut diviser par N pas N(N-1) car la moyenne d'un γ̃_j
    # à distance x d'un γ̃_i fixé est 1·dx pour Poisson, et nous sommons
    # sur tous les i (donne facteur N).

    R2_emp = counts / (N * bin_width)

    # Bruit Poisson sur R2
    noise_R2 = np.sqrt(counts) / (N * bin_width)

    if verbose:
        n_in_window = len(in_window)
        print(f"  [{window_label}] N={N} zeros, n_paires(i<j) ∈ [{x_min},{x_max}] = {n_in_window}")
        print(f"  N × bin_width = {N * bin_width:.1f}")

    return centers, R2_emp, counts, noise_R2


# ============================================================================
# 6. RUN — Echelle 1 (γ_n petits)
# ============================================================================

print("=" * 78)
print(" ECHELLE 1 : γ_n pour n ∈ [1, {}]".format(SCALE_LOW_MAX_N))
print("=" * 78)

zeros_scale_1 = load_zeros(1, SCALE_LOW_MAX_N)
unfolded_1 = unfold_zeros(zeros_scale_1)
print(f"  γ_1 = {zeros_scale_1[0]:.4f}, γ̃_1 = {unfolded_1[0]:.4f}")
print(f"  γ_{SCALE_LOW_MAX_N} = {zeros_scale_1[-1]:.4f}, γ̃_{SCALE_LOW_MAX_N} = {unfolded_1[-1]:.4f}")
print(f"  écart moyen unfolded = {(unfolded_1[-1] - unfolded_1[0]) / (SCALE_LOW_MAX_N - 1):.6f} (attendu : ≈ 1)")
print()

x1, R2_1, counts_1, noise_1 = empirical_pair_correlation(
    unfolded_1, window_label=f"n ∈ [1, {SCALE_LOW_MAX_N}]"
)

K_mont_1 = montgomery_kernel(x1)
K_pt_omega_1 = pt_kernel(x1, scheme='pt_omega')
K_pt_factor_1 = pt_kernel(x1, scheme='pt_factor')

print()
print(f"  {'x':>6} | {'R2_emp':>10} | {'K_Mont':>10} | {'K_PT_omega':>10} | {'K_PT_factor':>10} | {'noise':>8} | {'n':>6}")
print("  " + "-" * 86)
for i in range(0, len(x1)):
    print(f"  {x1[i]:>6.3f} | {R2_1[i]:>10.4f} | {K_mont_1[i]:>10.4f} | "
          f"{K_pt_omega_1[i]:>10.4f} | {K_pt_factor_1[i]:>10.4f} | {noise_1[i]:>8.4f} | {counts_1[i]:>6d}")

# Statistiques (poids inverse au bruit pour le χ²)
# NOTE: pour bins avec count très bas, le bruit Poisson n'est pas dans
# le régime gaussien, donc on impose un floor PLUS GRAND (= moyenne de noise)
# pour éviter que des bins quasi-vides dominent le chi-deux.
floor_noise_1 = max(0.05, np.mean(noise_1) * 0.5)
chi2_mont_1 = np.sum(((R2_1 - K_mont_1) / np.maximum(noise_1, floor_noise_1))**2)
chi2_pt_omega_1 = np.sum(((R2_1 - K_pt_omega_1) / np.maximum(noise_1, floor_noise_1))**2)
chi2_pt_factor_1 = np.sum(((R2_1 - K_pt_factor_1) / np.maximum(noise_1, floor_noise_1))**2)

print()
print(f"  χ²(Montgomery) = {chi2_mont_1:.2f}  (dof ~ {len(x1)})")
print(f"  χ²(PT omega)   = {chi2_pt_omega_1:.2f}")
print(f"  χ²(PT factor)  = {chi2_pt_factor_1:.2f}")
print(f"  χ²(PT)/χ²(Mont) — ratio omega = {chi2_pt_omega_1/chi2_mont_1:.3f}")
print(f"  χ²(PT)/χ²(Mont) — ratio factor = {chi2_pt_factor_1/chi2_mont_1:.3f}")


# ============================================================================
# 7. RUN — Echelle 2 (γ_n grands)
# ============================================================================

print()
print("=" * 78)
print(f" ECHELLE 2 : γ_n pour n ∈ [{SCALE_HIGH_START_N}, {SCALE_HIGH_END_N}]")
print("=" * 78)

zeros_scale_2 = load_zeros(SCALE_HIGH_START_N, SCALE_HIGH_END_N)
unfolded_2 = unfold_zeros(zeros_scale_2)
N2 = len(zeros_scale_2)
print(f"  γ_{SCALE_HIGH_START_N} = {zeros_scale_2[0]:.4f}, γ̃ = {unfolded_2[0]:.4f}")
print(f"  γ_{SCALE_HIGH_END_N} = {zeros_scale_2[-1]:.4f}, γ̃ = {unfolded_2[-1]:.4f}")
print(f"  écart moyen unfolded = {(unfolded_2[-1] - unfolded_2[0]) / (N2 - 1):.6f}")
print()

x2, R2_2, counts_2, noise_2 = empirical_pair_correlation(
    unfolded_2, window_label=f"n ∈ [{SCALE_HIGH_START_N}, {SCALE_HIGH_END_N}]"
)

K_mont_2 = montgomery_kernel(x2)
K_pt_omega_2 = pt_kernel(x2, scheme='pt_omega')
K_pt_factor_2 = pt_kernel(x2, scheme='pt_factor')

print()
print(f"  {'x':>6} | {'R2_emp':>10} | {'K_Mont':>10} | {'K_PT_omega':>10} | {'K_PT_factor':>10} | {'noise':>8} | {'n':>6}")
print("  " + "-" * 86)
for i in range(0, len(x2)):
    print(f"  {x2[i]:>6.3f} | {R2_2[i]:>10.4f} | {K_mont_2[i]:>10.4f} | "
          f"{K_pt_omega_2[i]:>10.4f} | {K_pt_factor_2[i]:>10.4f} | {noise_2[i]:>8.4f} | {counts_2[i]:>6d}")

floor_noise_2 = max(0.05, np.mean(noise_2) * 0.5)
chi2_mont_2 = np.sum(((R2_2 - K_mont_2) / np.maximum(noise_2, floor_noise_2))**2)
chi2_pt_omega_2 = np.sum(((R2_2 - K_pt_omega_2) / np.maximum(noise_2, floor_noise_2))**2)
chi2_pt_factor_2 = np.sum(((R2_2 - K_pt_factor_2) / np.maximum(noise_2, floor_noise_2))**2)

print()
print(f"  χ²(Montgomery) = {chi2_mont_2:.2f}  (dof ~ {len(x2)})")
print(f"  χ²(PT omega)   = {chi2_pt_omega_2:.2f}")
print(f"  χ²(PT factor)  = {chi2_pt_factor_2:.2f}")
print(f"  ratio omega = {chi2_pt_omega_2/chi2_mont_2:.3f}")
print(f"  ratio factor = {chi2_pt_factor_2/chi2_mont_2:.3f}")


# ============================================================================
# 8. TEST RATIO 2 DÉCADES + DIAGNOSTIC
# ============================================================================

print()
print("=" * 78)
print(" TEST RATIO 2 DÉCADES (discipline mémoire)")
print("=" * 78)

print(f"  Ratio χ²(Mont)/χ²(PT omega) sur 2 décades :")
print(f"    Echelle 1 : {chi2_mont_1/chi2_pt_omega_1:.4f}")
print(f"    Echelle 2 : {chi2_mont_2/chi2_pt_omega_2:.4f}")
print()
print(f"  Ratio χ²(Mont)/χ²(PT factor) sur 2 décades :")
print(f"    Echelle 1 : {chi2_mont_1/chi2_pt_factor_1:.4f}")
print(f"    Echelle 2 : {chi2_mont_2/chi2_pt_factor_2:.4f}")
print()

# RMSE empirique
rmse_mont_1 = np.sqrt(np.mean((R2_1 - K_mont_1)**2))
rmse_mont_2 = np.sqrt(np.mean((R2_2 - K_mont_2)**2))
rmse_pt_omega_1 = np.sqrt(np.mean((R2_1 - K_pt_omega_1)**2))
rmse_pt_omega_2 = np.sqrt(np.mean((R2_2 - K_pt_omega_2)**2))
rmse_pt_factor_1 = np.sqrt(np.mean((R2_1 - K_pt_factor_1)**2))
rmse_pt_factor_2 = np.sqrt(np.mean((R2_2 - K_pt_factor_2)**2))

print(f"  RMSE empirique :")
print(f"    Echelle 1 : RMSE vs Mont = {rmse_mont_1:.4f},  vs PT_omega = {rmse_pt_omega_1:.4f}")
print(f"    Echelle 2 : RMSE vs Mont = {rmse_mont_2:.4f},  vs PT_omega = {rmse_pt_omega_2:.4f}")
print(f"  Bruit Poisson moyen :")
print(f"    Echelle 1 : σ̄ = {np.mean(noise_1):.4f}")
print(f"    Echelle 2 : σ̄ = {np.mean(noise_2):.4f}")
print()

# Comparaison qualité Mont :
# Echelle 1 a moins de paires, plus de bruit → RMSE/σ peut être différent
# Si Montgomery est correct, on attend RMSE ~ σ (à un facteur constant)
print(f"  Ratio RMSE/σ̄ vs Mont :")
print(f"    Echelle 1 : {rmse_mont_1/np.mean(noise_1):.2f}")
print(f"    Echelle 2 : {rmse_mont_2/np.mean(noise_2):.2f}")
print(f"  Si Mont parfait, attendu ~ 1. Si > 2, déviation systématique présente.")


# ============================================================================
# 9. VERDICT
# ============================================================================

print()
print("=" * 78)
print(" VERDICT V10")
print("=" * 78)

# Critères basés sur RMSE et test ratio bruit Poisson (plus robuste que χ²)
# RMSE/σ proche de 1 = ajustement compatible avec bruit Poisson
# RMSE/σ > 3 = déviation systématique présente
# Ratio RMSE(Mont)/RMSE(PT) > 1.3 = PT préféré ; < 0.77 = Mont préféré

threshold = 1.3

ratio_omega_1 = rmse_mont_1 / rmse_pt_omega_1
ratio_omega_2 = rmse_mont_2 / rmse_pt_omega_2
ratio_factor_1 = rmse_mont_1 / rmse_pt_factor_1 if rmse_pt_factor_1 > 0 else 0
ratio_factor_2 = rmse_mont_2 / rmse_pt_factor_2 if rmse_pt_factor_2 > 0 else 0

print(f"  RMSE(Mont)/RMSE(PT omega) sur 2 décades :")
print(f"    Echelle 1 : {ratio_omega_1:.3f}")
print(f"    Echelle 2 : {ratio_omega_2:.3f}")
print(f"  RMSE(Mont)/RMSE(PT factor) sur 2 décades :")
print(f"    Echelle 1 : {ratio_factor_1:.3f}")
print(f"    Echelle 2 : {ratio_factor_2:.3f}")
print()

if (ratio_omega_1 > threshold and ratio_omega_2 > threshold):
    verdict = "PT_OMEGA_WINS"
    interpret = ("PT omega schéma BETTER que Mont sur 2 décades. "
                 "CAS 2 candidat : pair correlation suit substitution PT.")
elif (1/ratio_omega_1 > threshold and 1/ratio_omega_2 > threshold):
    verdict = "MONT_WINS_VS_PT_OMEGA"
    interpret = ("PT omega schéma FALSIFIÉ : Montgomery fitte beaucoup mieux. "
                 "Substitution π → √(3Σγ/4) ne décrit PAS la pair correlation.")
elif (ratio_factor_1 > threshold and ratio_factor_2 > threshold):
    verdict = "PT_FACTOR_WINS"
    interpret = ("PT factor schéma BETTER que Mont sur 2 décades. "
                 "CAS 2 candidat : amplitude PT meilleure.")
elif (1/ratio_factor_1 > threshold and 1/ratio_factor_2 > threshold):
    verdict = "MONT_WINS_VS_PT_FACTOR"
    interpret = ("PT factor schéma FALSIFIÉ : Montgomery fitte mieux.")
else:
    verdict = "INCONCLUSIVE_OR_DIRECT"
    interpret = ("Aucun schéma alternatif PT significativement préféré. "
                 "Montgomery est le bon kernel. PT direct (= Mont) est implicite. "
                 "Cas 3-4 : ne ferme pas HP.")

print(f"\n  Verdict : {verdict}")
print(f"  Interprétation : {interpret}")
print()

# Cas particulier : si Mont est BIEN ajusté, PT direct l'est aussi
# (puisqu'on a défini K_PT_direct = K_Mont). C'est l'interprétation
# minimale et la plus probable.
if rmse_mont_1 < 2 * np.mean(noise_1) and rmse_mont_2 < 2 * np.mean(noise_2):
    print(f"  Montgomery fitte BIEN les données (RMSE < 2σ Poisson).")
    print(f"  Conclusion supplémentaire : PT direct (= Mont) reproduit la pair correlation.")
    print(f"  Mais ne fournit pas de nouvelle prédiction PT-spécifique.")
elif rmse_mont_1 > 3 * np.mean(noise_1) or rmse_mont_2 > 3 * np.mean(noise_2):
    print(f"  Montgomery dévie significativement (RMSE > 3σ). Anomalie à investiguer.")


# ============================================================================
# 10. SAUVEGARDE
# ============================================================================

outputs_dir = Path("/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT_RH_MAY/outputs")
outputs_dir.mkdir(exist_ok=True)

results = {
    "verdict": verdict,
    "interpretation": interpret,
    "primes_actifs": PRIMES,
    "mu_star": float(MU_STAR),
    "sum_gamma_p": float(SUM_GAMMA),
    "ratio_3sumgamma_over_4pi2": 0.75 * float(SUM_GAMMA) / (np.pi**2),
    "scale_1": {
        "n_range": [1, SCALE_LOW_MAX_N],
        "gamma_range": [float(zeros_scale_1[0]), float(zeros_scale_1[-1])],
        "x_centers": x1.tolist(),
        "R2_emp": R2_1.tolist(),
        "K_montgomery": K_mont_1.tolist(),
        "K_pt_omega": K_pt_omega_1.tolist(),
        "K_pt_factor": K_pt_factor_1.tolist(),
        "noise_R2": noise_1.tolist(),
        "counts": counts_1.tolist(),
        "chi2_mont": float(chi2_mont_1),
        "chi2_pt_omega": float(chi2_pt_omega_1),
        "chi2_pt_factor": float(chi2_pt_factor_1),
        "rmse_vs_mont": float(rmse_mont_1),
        "noise_mean": float(np.mean(noise_1)),
    },
    "scale_2": {
        "n_range": [SCALE_HIGH_START_N, SCALE_HIGH_END_N],
        "gamma_range": [float(zeros_scale_2[0]), float(zeros_scale_2[-1])],
        "x_centers": x2.tolist(),
        "R2_emp": R2_2.tolist(),
        "K_montgomery": K_mont_2.tolist(),
        "K_pt_omega": K_pt_omega_2.tolist(),
        "K_pt_factor": K_pt_factor_2.tolist(),
        "noise_R2": noise_2.tolist(),
        "counts": counts_2.tolist(),
        "chi2_mont": float(chi2_mont_2),
        "chi2_pt_omega": float(chi2_pt_omega_2),
        "chi2_pt_factor": float(chi2_pt_factor_2),
        "rmse_vs_mont": float(rmse_mont_2),
        "noise_mean": float(np.mean(noise_2)),
    }
}

out_file = outputs_dir / "v10_pair_correlation_results.json"
with open(out_file, "w") as f:
    json.dump(results, f, indent=2)
print(f"\n  Résultats sauvegardés dans : {out_file}")
print()
print("=" * 78)
print(" FIN V10")
print("=" * 78)
