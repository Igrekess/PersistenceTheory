"""
V12 — Test : les zéros γ_n comme NŒUDS D'INTERFÉRENCE DESTRUCTIVE
       de la cascade de persistance sur la ligne critique,
       et symétrie T_3 → s ↔ 1-s.

HYPOTHÈSE MÉCANIQUE PT :
   La cascade de persistance projette chaque shell premier en oscillation
       ω_p(t) = sin(t · log p) · A_p(s)
   sur la ligne critique. L'onde spectrale résultante est
       R(1/2 + it) = ζ(1/2+it) / [ζ_+(1/2+it) · ζ_-(1/2+it)]
   et ses NŒUDS (zéros) sont les γ_n.

QUESTION FONDAMENTALE :
   La symétrie T_3 (involution antidiagonale mod 3) → s ↔ 1-s
   FORCE-T-ELLE les nœuds sur l'axe Re(s) = 1/2 ?

RÉPONSE STRUCTURELLE : NON, en général. Une fonction f symétrique sous
   s ↔ 1-s peut avoir des zéros par paires {s_0, 1-s_0} hors-axe.

OBJECTIF V12 :
   (A) Formaliser l'onde R(1/2+it) comme superposition explicite.
   (B) Identifier précisément le GAP entre symétrie T_3 et RH stricte.
   (C) Vérifier numériquement |R(σ+iτ)| sur grille 2D et identifier
       si les minima sont EXCLUSIVEMENT sur σ = 1/2.

DISCIPLINE :
   - Pas de "conjecture de Lehmer-like" inutilement chargée
   - Distinguer prouvé / heuristique
   - Test ratio sur 2 décades en τ
   - Si zéros hors-axe trouvés → contre-exemple à RH (impossible si RH vraie)
   - Si non → vérification empirique (sans preuve)
"""

import mpmath as mp
import numpy as np
from sympy import primerange
from pathlib import Path
import json
import time

mp.mp.dps = 25  # 25 chiffres significatifs (rapide mais précis)

# ============================================================================
# CONSTANTES PT
# ============================================================================

Q_PLUS = mp.mpf(13) / mp.mpf(15)
Q_MINUS = mp.exp(mp.mpf(-1) / mp.mpf(15))


def delta_p(p, q):
    """δ_p^σ = (1 - q^p)/p"""
    p = mp.mpf(p)
    return (mp.mpf(1) - q**p) / p


def a_p(p, q):
    """a_p^σ = δ_p^σ · (2 - δ_p^σ) = sin²(θ_p, q)"""
    d = delta_p(p, q)
    return d * (mp.mpf(2) - d)


def A_p(p):
    """A_p = a_p^+ + a_p^-"""
    return a_p(p, Q_PLUS) + a_p(p, Q_MINUS)


# ============================================================================
# A — FORMALISATION DE L'ONDE SPECTRALE R(s)
# ============================================================================

def log_zeta_pm(s, p_max):
    """log[ζ_+(s) ζ_-(s)] = Σ_p A_p · p^{-s}"""
    s = mp.mpc(s)
    total = mp.mpc(0)
    for p in primerange(2, p_max + 1):
        total += A_p(p) * mp.power(mp.mpc(p), -s)
    return total


def zeta_pm(s, p_max=2000):
    """ζ_+(s) · ζ_-(s) = exp(Σ A_p · p^{-s})"""
    return mp.exp(log_zeta_pm(s, p_max))


def R_spectral(s, p_max=2000):
    """R(s) = ζ(s) / [ζ_+(s) ζ_-(s)]"""
    s = mp.mpc(s)
    z = mp.zeta(s)
    zpm = zeta_pm(s, p_max)
    return z / zpm


def log_R_explicit(s, p_max=2000):
    """log R(s) = log ζ(s) - Σ A_p · p^{-s}

    Sur Re(s) > 1, on a aussi (preuve par produit d'Euler) :
        log ζ(s) = Σ_p Σ_{k≥1} p^{-ks}/k
    donc :
        log R(s) = Σ_p [ Σ_{k≥1} p^{-ks}/k - A_p · p^{-s} ]
                 = Σ_p p^{-s} · [(1/(1-p^{-s})) - 1 - A_p] + higher orders ...

    Au premier ordre en p^{-s} (terme dominant pour Re(s) près de 1/2,
    correction par les ordres supérieurs k≥2) :
        log R(s) ≈ Σ_p p^{-s} (1 - A_p)
    avec 1 - A_p ≈ 1 - 4/p pour p grand.

    Donc R(s) ~ exp[Σ_p (1 - A_p) p^{-s}]  pour le terme dominant.

    Sur Re(s) = 1/2, s = 1/2 + it, p^{-s} = p^{-1/2} · exp(-i t log p).
    """
    s = mp.mpc(s)
    z = mp.zeta(s)
    return mp.log(z) - log_zeta_pm(s, p_max)


# ============================================================================
# B — DÉCOMPOSITION AMPLITUDE / PHASE SUR LIGNE CRITIQUE
# ============================================================================

def harmonic_amplitudes_critical(t, p_max=2000):
    """Décomposition de log R(1/2 + it) en superposition harmonique.

    Returns dict :
        'log_R_zeta'  : log ζ(1/2+it) (depuis mpmath)
        'log_zeta_pm' : log[ζ_+ζ_-](1/2+it) (somme PT)
        'log_R'       : log R(1/2+it) = log ζ - log[ζ_+ζ_-]
        'harmonics'   : liste de [(p, amplitude, phase)] pour chaque prime
                         où amplitude = (1-A_p) · p^{-1/2}
                              phase = -t · log p

    Interprétation : chaque shell premier p contribue une oscillation
        c_p · cos(t log p - phi_p)
    avec c_p l'amplitude PT-modulée et phi_p déterminé par 1-A_p.
    """
    s = mp.mpc(mp.mpf("0.5"), t)
    harmonics = []
    for p in primerange(2, p_max + 1):
        amp = (mp.mpf(1) - A_p(p)) / mp.sqrt(mp.mpf(p))
        phase = -t * mp.log(p)
        harmonics.append((int(p), amp, phase))
    log_R_z = mp.log(mp.zeta(s))
    log_z_pm = log_zeta_pm(s, p_max)
    log_R_val = log_R_z - log_z_pm
    return {
        "log_R_zeta": log_R_z,
        "log_zeta_pm": log_z_pm,
        "log_R": log_R_val,
        "harmonics": harmonics,
    }


# ============================================================================
# C — TEST NUMÉRIQUE : SCAN |R(σ + iτ)| SUR GRILLE 2D
# ============================================================================

def scan_R_grid(sigma_range, tau_range, n_sigma=21, n_tau=200, p_max=1000):
    """Calcule |R(σ + iτ)| sur grille rectangulaire.

    Returns: (sigmas, taus, abs_R) avec abs_R[i, j] = |R(σ_i + iτ_j)|.
    """
    sigmas = np.linspace(sigma_range[0], sigma_range[1], n_sigma)
    taus = np.linspace(tau_range[0], tau_range[1], n_tau)
    abs_R = np.zeros((n_sigma, n_tau))
    t_start = time.time()
    for i, sigma in enumerate(sigmas):
        for j, tau in enumerate(taus):
            s = mp.mpc(mp.mpf(float(sigma)), mp.mpf(float(tau)))
            try:
                val = R_spectral(s, p_max=p_max)
                abs_R[i, j] = float(abs(val))
            except Exception as e:
                abs_R[i, j] = float("nan")
        if (i + 1) % 5 == 0 or i == n_sigma - 1:
            elapsed = time.time() - t_start
            print(f"  σ row {i+1}/{n_sigma} (σ={sigma:.3f}) done, "
                  f"elapsed {elapsed:.1f}s")
    return sigmas, taus, abs_R


def find_local_minima_2d(grid):
    """Identifie les minima locaux d'une grille 2D (par voisinage 8-connexe)."""
    minima = []
    n_rows, n_cols = grid.shape
    for i in range(1, n_rows - 1):
        for j in range(1, n_cols - 1):
            center = grid[i, j]
            is_min = True
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    if grid[i + di, j + dj] < center:
                        is_min = False
                        break
                if not is_min:
                    break
            if is_min:
                minima.append((i, j, center))
    return minima


# ============================================================================
# D — VÉRIFICATION RIGOUREUSE : Zéros connus de ζ vs zéros de R
# ============================================================================

def known_zeros(n_max=20):
    """Retourne les n_max premiers γ_n via mpmath."""
    return [float(mp.im(mp.zetazero(n))) for n in range(1, n_max + 1)]


def R_at_known_zeros(p_max=2000):
    """Vérifie que R(1/2 + i γ_n) ≈ 0 (puisque ζ(1/2+iγ_n) = 0 et ζ_+ζ_- ≠ 0)."""
    gammas = known_zeros(15)
    print(f"\nVérification : R(1/2 + i γ_n) = ζ(1/2+iγ_n) / [ζ_+ζ_-(1/2+iγ_n)]")
    print(f"     Théoriquement = 0/quelque chose ≠ 0 = 0")
    print(f"  n   |  γ_n        |  |R|         |  |ζ|         |  |ζ_+ζ_-|")
    print(f" {'-'*60}")
    for k, g in enumerate(gammas, 1):
        s = mp.mpc(mp.mpf("0.5"), mp.mpf(g))
        z = mp.zeta(s)
        zpm = zeta_pm(s, p_max)
        R_val = z / zpm
        print(f"  {k:2d}  |  {g:9.5f}  |  {float(abs(R_val)):.4e}  |  "
              f"{float(abs(z)):.4e}  |  {float(abs(zpm)):.6f}")


# ============================================================================
# E — TEST CRITIQUE : recherche de zéros HORS-AXE sur grille fine
# ============================================================================

def search_offaxis_zeros(sigma_range=(0.35, 0.65), tau_range=(10, 60),
                        n_sigma=25, n_tau=400, p_max=1500, abs_threshold=0.05):
    """Recherche systématique de minima |R(σ+iτ)| < seuil sur grille fine.

    Sur Re(s)=1/2 on s'attend à des zéros aux γ_n connus.
    On VÉRIFIE qu'aucun minimum significatif n'apparaît HORS de σ=0.5.

    Si on en trouve avec σ ≠ 0.5 et |R| << 1 → contre-exemple à RH.
    """
    print(f"\n[E] Recherche zéros hors-axe : grille {n_sigma}×{n_tau}")
    print(f"   σ ∈ [{sigma_range[0]}, {sigma_range[1]}], "
          f"τ ∈ [{tau_range[0]}, {tau_range[1]}]")
    print(f"   p_max = {p_max}, seuil |R| < {abs_threshold}")
    sigmas, taus, abs_R = scan_R_grid(sigma_range, tau_range,
                                       n_sigma, n_tau, p_max)
    minima = find_local_minima_2d(abs_R)
    print(f"\n   {len(minima)} minima locaux trouvés, seuil < {abs_threshold}:")
    significant = [(i, j, v) for (i, j, v) in minima if v < abs_threshold]
    onaxis = []
    offaxis = []
    # Indice central : σ = 0.5
    sigma_center_idx = np.argmin(np.abs(sigmas - 0.5))
    sigma_step = sigmas[1] - sigmas[0]
    for (i, j, v) in significant:
        sigma_min = sigmas[i]
        tau_min = taus[j]
        # Tolérance "sur axe" = 1 step de σ
        if abs(i - sigma_center_idx) <= 1:
            onaxis.append((sigma_min, tau_min, v))
        else:
            offaxis.append((sigma_min, tau_min, v))
    print(f"   Minima sur axe σ=0.5±{sigma_step:.3f} : {len(onaxis)}")
    for (s, t, v) in onaxis[:10]:
        print(f"     σ={s:.4f}, τ={t:.4f}, |R|={v:.5e}")
    print(f"   Minima HORS axe σ≠0.5 : {len(offaxis)}")
    for (s, t, v) in offaxis[:10]:
        print(f"     σ={s:.4f}, τ={t:.4f}, |R|={v:.5e}")
    return {
        "sigmas": sigmas.tolist(),
        "taus": taus.tolist(),
        "abs_R": abs_R.tolist(),
        "minima_onaxis": [(float(s), float(t), float(v)) for (s, t, v) in onaxis],
        "minima_offaxis": [(float(s), float(t), float(v)) for (s, t, v) in offaxis],
        "p_max": p_max,
        "abs_threshold": float(abs_threshold),
    }


# ============================================================================
# F — TEST DES PROFILS VERTICAUX |R(σ + iτ_fixed)|
# ============================================================================

def vertical_profile(tau_fixed, sigma_range=(0.1, 0.9), n_sigma=200, p_max=1500):
    """Profil de |R(σ + i·tau_fixed)| comme fonction de σ.

    Si γ_n est zéro de R sur Re=1/2, le profil à τ = γ_n doit montrer
    minimum à σ = 0.5.
    """
    sigmas = np.linspace(sigma_range[0], sigma_range[1], n_sigma)
    abs_R = np.zeros(n_sigma)
    for i, sigma in enumerate(sigmas):
        s = mp.mpc(mp.mpf(float(sigma)), mp.mpf(float(tau_fixed)))
        try:
            abs_R[i] = float(abs(R_spectral(s, p_max=p_max)))
        except Exception:
            abs_R[i] = float("nan")
    min_idx = np.argmin(abs_R)
    return sigmas, abs_R, sigmas[min_idx], abs_R[min_idx]


# ============================================================================
# G — TEST DE SYMÉTRIE EXPLICITE : R(s) vs R(1-s)
# ============================================================================

def test_functional_equation(t_values, p_max=2000):
    """Vérifie le lien entre R(s) et R(1-s) au regard de l'équation
    fonctionnelle de ζ (asymptotique).

    Pour ζ : ξ(s) = ξ(1-s) avec ξ(s) = (1/2)s(s-1)π^{-s/2}Γ(s/2)ζ(s).
    Donc ζ(1-s) = chi(s) · ζ(s) avec chi(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s).

    Pour R(s) = ζ(s)/(ζ_+(s)ζ_-(s)), il y a un facteur supplémentaire :
    R(1-s) = ζ(1-s)/(ζ_+(1-s)ζ_-(1-s))
           = chi(s) ζ(s) / (ζ_+(1-s) ζ_-(1-s))
           = chi(s) · [ζ_+(s)ζ_-(s) / ζ_+(1-s)ζ_-(1-s)] · R(s)
           = chi_PT(s) · R(s)
    avec chi_PT(s) = chi(s) · [ζ_+ζ_-(s) / ζ_+ζ_-(1-s)].

    Donc R(s) et R(1-s) sont liés par un facteur multiplicatif non nul
    sur la bande critique (par INC-1 : ζ_+ζ_- ≠ 0).
    Les zéros de R sont donc symétriques sous s ↔ 1-s (par paires).
    """
    print("\n[G] Test symétrie R(s) vs R(1-s)")
    print("   Vérifions que R(s) = 0 ⟺ R(1-s) = 0 (zéros symétriques)")
    print(f"  τ        | |R(1/2+iτ)| | |R(1/2-iτ)| | ratio  | chi_PT magnitude")
    print(f" {'-'*75}")
    for t in t_values:
        s = mp.mpc(mp.mpf("0.5"), mp.mpf(t))
        s_conj = mp.mpc(mp.mpf("0.5"), -mp.mpf(t))
        # 1 - s = 1 - 1/2 - it = 1/2 - it = s_conj sur l'axe
        R_s = R_spectral(s, p_max)
        R_1ms = R_spectral(s_conj, p_max)
        # Calcul chi(s)
        chi_s = mp.mpc(2) ** s * mp.pi ** (s - 1) * mp.sin(mp.pi * s / 2) * mp.gamma(1 - s)
        zpm_s = zeta_pm(s, p_max)
        zpm_1ms = zeta_pm(s_conj, p_max)
        chi_pt = chi_s * (zpm_s / zpm_1ms)
        print(f"  {float(t):6.3f}  |  {float(abs(R_s)):.4e}  |  "
              f"{float(abs(R_1ms)):.4e}  |  "
              f"{float(abs(R_s)/abs(R_1ms)):.4e}  |  {float(abs(chi_pt)):.4f}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print(" V12 — TEST : NŒUDS DE LA CASCADE PT vs SYMÉTRIE T_3")
    print("=" * 78)
    print()
    print(f"Précision mpmath dps = {mp.mp.dps}")
    print(f"q_+ = 13/15 = {float(Q_PLUS):.6f}")
    print(f"q_- = exp(-1/15) = {float(Q_MINUS):.6f}")
    print()

    # ========================================================================
    # ÉTAPE A — FORMALISATION ONDE SPECTRALE
    # ========================================================================
    print("─" * 78)
    print(" A — Décomposition harmonique de R(1/2 + it)")
    print("─" * 78)
    print()

    t_test = mp.mpf("14.135")  # proche du 1er zéro γ_1 = 14.1347...
    print(f"Test à t = {float(t_test)} (≈ γ_1 = 14.1347)")
    decomp = harmonic_amplitudes_critical(t_test, p_max=200)
    print(f"  log R          = {mp.nstr(decomp['log_R'], 12)}")
    print(f"  |R|            = {float(abs(mp.exp(decomp['log_R']))):.6e}")
    print(f"  log ζ          = {mp.nstr(decomp['log_R_zeta'], 12)}")
    print(f"  log[ζ_+ζ_-]    = {mp.nstr(decomp['log_zeta_pm'], 12)}")
    print()
    # Affiche les 10 plus grandes amplitudes (par valeur absolue 1-A_p/√p)
    harms_sorted = sorted(decomp['harmonics'],
                          key=lambda x: -float(abs(x[1])))
    print("  10 premières harmoniques (p, amp = (1-A_p)/√p, phase = -t log p):")
    for (p, amp, phase) in harms_sorted[:10]:
        print(f"     p={p:5d}, amp={float(amp):.5f}, phase={float(phase):.4f}")
    print()

    # ========================================================================
    # ÉTAPE D — VÉRIFICATION ZÉROS CONNUS
    # ========================================================================
    print("─" * 78)
    print(" D — Vérification : R(1/2 + i γ_n) ≈ 0 pour zéros connus de ζ")
    print("─" * 78)
    R_at_known_zeros(p_max=1500)
    print()

    # ========================================================================
    # ÉTAPE G — TEST SYMÉTRIE FONCTIONNELLE
    # ========================================================================
    print("─" * 78)
    print(" G — Symétrie R(s) vs R(1-s) sur axe (lien à équation fonctionnelle)")
    print("─" * 78)
    test_functional_equation([10.0, 14.135, 21.022, 25.011, 30.425],
                             p_max=1500)
    print()

    # ========================================================================
    # ÉTAPE F — PROFILS VERTICAUX à τ = γ_1, γ_2
    # ========================================================================
    print("─" * 78)
    print(" F — Profils verticaux |R(σ + iτ)| à τ = γ_1, γ_2, γ_3")
    print("─" * 78)
    print()
    gammas_known = known_zeros(5)
    profiles = {}
    for k, g in enumerate(gammas_known[:3], 1):
        print(f"\n  Profil à τ = γ_{k} = {g:.5f}:")
        sigmas_p, abs_R_p, sigma_min, val_min = vertical_profile(
            mp.mpf(g), sigma_range=(0.20, 0.80), n_sigma=121, p_max=1500)
        print(f"     σ argmin = {sigma_min:.4f}, |R| = {val_min:.4e}")
        # affiche valeurs autour de σ=0.5
        idx_center = np.argmin(np.abs(sigmas_p - 0.5))
        print(f"     |R(0.40+iτ)| = {abs_R_p[np.argmin(np.abs(sigmas_p-0.40))]:.4e}")
        print(f"     |R(0.45+iτ)| = {abs_R_p[np.argmin(np.abs(sigmas_p-0.45))]:.4e}")
        print(f"     |R(0.50+iτ)| = {abs_R_p[idx_center]:.4e}")
        print(f"     |R(0.55+iτ)| = {abs_R_p[np.argmin(np.abs(sigmas_p-0.55))]:.4e}")
        print(f"     |R(0.60+iτ)| = {abs_R_p[np.argmin(np.abs(sigmas_p-0.60))]:.4e}")
        profiles[f"gamma_{k}"] = {
            "gamma_value": g,
            "sigmas": sigmas_p.tolist(),
            "abs_R": abs_R_p.tolist(),
            "sigma_argmin": float(sigma_min),
            "val_min": float(val_min),
        }
    print()

    # ========================================================================
    # ÉTAPE E — RECHERCHE ZÉROS HORS-AXE
    # ========================================================================
    print("─" * 78)
    print(" E — Recherche zéros HORS-AXE : grille fine sur [0.35, 0.65] × [10, 50]")
    print("─" * 78)
    grid_results = search_offaxis_zeros(
        sigma_range=(0.35, 0.65),
        tau_range=(10, 50),
        n_sigma=15,
        n_tau=200,
        p_max=1500,
        abs_threshold=0.10,
    )
    print()

    # ========================================================================
    # SAUVEGARDE
    # ========================================================================
    out_path = Path("/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/"
                    "PT_RH_MAY/outputs/v12_nodes_symmetry_results.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "mpmath_dps": mp.mp.dps,
        "PT_constants": {
            "q_plus_str": "13/15",
            "q_minus_str": "exp(-1/15)",
        },
        "vertical_profiles": profiles,
        "grid_offaxis_search": grid_results,
    }
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    print(f"\nRésultats sauvegardés : {out_path}")

    # ========================================================================
    # CONCLUSION
    # ========================================================================
    print()
    print("=" * 78)
    print(" CONCLUSION V12")
    print("=" * 78)
    n_off = len(grid_results["minima_offaxis"])
    n_on = len(grid_results["minima_onaxis"])
    print(f"\n  • Minima locaux trouvés : {n_on} sur axe + {n_off} hors-axe")
    if n_off == 0:
        print(f"  → AUCUN minimum significatif hors-axe trouvé (à seuil 0.10)")
        print(f"    Vérification empirique cohérente avec RH SUR cette plage.")
    else:
        print(f"  → {n_off} minima hors-axe trouvés (= possibles candidats")
        print(f"    contre-exemple, OU artefacts de discrétisation à examiner)")
    print()
    print("  GAP IDENTIFIÉ (étape B de la note d'analyse):")
    print("    Symétrie T_3 → R(s)R(1-s) symétrique sous s↔1-s")
    print("    Mais : symétrie SEULE n'exclut PAS zéros par paires hors-axe")
    print("    Il faut un INGRÉDIENT POSITIVITÉ supplémentaire")
    print("    (Weil-positivity, ou rigidité spectrale dégénérée)")
    print()


if __name__ == "__main__":
    main()
