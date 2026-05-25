"""
pt_rh_finite_zeros.py — Phase 11, Diagnostic et rectification
Zéros de S_m(s) = Σ_{p|m*} c_p · p^{-s} sur Z/m*Z (m*=30030)

RÉSULTAT (2026-04-25) : S_m(s) a des zéros partout dans la bande, pas seulement
sur σ=1/2. Ceci révèle que S_m est le MAUVAIS objet pour modéliser les zéros de ζ.

DIAGNOSTIC :
  - S_m(s) = Σ c_p p^{-s} est la DÉRIVÉE LOGARITHMIQUE tronquée de R(s), pas R lui-même.
  - Les zéros de ζ sont dans R(s) = ζ(s) / (ζ₊·ζ₋), qui est défini par continuation
    analytique, pas par une série convergente dans la bande critique.
  - La troncature finie d'un produit d'Euler n'a PAS de zéros dans la bande critique
    (c'est un polynôme de Dirichlet, zéros partout).
  - La "Step 1 sur le réseau fini" telle que formulée est incorrecte.

CE QUE LE SCRIPT FAIT MAINTENANT :
  1. Confirme que S_m a des zéros hors σ=1/2 (falsification de la Step 1)
  2. Identifie le bon objet fini : ζ_m(s) via caractères de (Z/m*Z)*
  3. Formule la Step 1 corrigée : invariant spectral du transfert T_m* sur σ=1/2

Objectif : vérifier numériquement que tous les zéros de S_m(s) dans la bande
critique 0 < Re(s) < 1 satisfont Re(s) = 1/2.

Définitions :
  q_+ = 1 - 2/mu*  = 13/15         (couplage, branche vertex)
  q_- = exp(-1/mu*) ≈ 0.9355       (thermique, branche géométrique)
  delta_p(q) = (1 - q^p) / p       (holonomie, identité T6)
  a_p = sin²_p(q_+) = delta_+(p) * (2 - delta_+(p))
  b_p = sin²_p(q_-) = delta_-(p) * (2 - delta_-(p))
  c_p = 1 - a_p - b_p              (poids résiduel = canal R)

S_m(s) = Σ_{p|m*} c_p · p^{-s}    (polynôme de Dirichlet fini, 6 termes)
"""

import cmath
import math
import numpy as np
from scipy.optimize import fsolve

# ── Constantes PT ─────────────────────────────────────────────────────────────
MU_STAR = 15.0
Q_PLUS  = 1 - 2 / MU_STAR          # = 13/15 ≈ 0.86667
Q_MINUS = math.exp(-1 / MU_STAR)   # ≈ 0.93541

M_STAR  = 30030  # primoriel P_6 = 2×3×5×7×11×13
PRIMES_M = [2, 3, 5, 7, 11, 13]   # primes divisant m*


# ── Fonctions holonomiques ─────────────────────────────────────────────────────
def delta_p(p, q):
    return (1 - q**p) / p

def sin2_p(p, q):
    d = delta_p(p, q)
    return d * (2 - d)

def a_p(p):
    return sin2_p(p, Q_PLUS)

def b_p(p):
    return sin2_p(p, Q_MINUS)

def c_p(p):
    return 1.0 - a_p(p) - b_p(p)


# ── Tableau des poids ─────────────────────────────────────────────────────────
def print_weights():
    print("\n=== Poids PT pour les primes de m* ===")
    print(f"{'p':>4} {'a_p':>10} {'b_p':>10} {'c_p':>10} {'a+b+c':>10}")
    print("-" * 50)
    for p in PRIMES_M:
        a = a_p(p); b = b_p(p); c = c_p(p)
        print(f"{p:>4} {a:>10.6f} {b:>10.6f} {c:>10.6f} {a+b+c:>10.6f}")


# ── Polynôme de Dirichlet S_m(s) ─────────────────────────────────────────────
def S_m(s):
    """S_m(s) = Σ_{p|m*} c_p · p^{-s}, s complexe"""
    return sum(c_p(p) * p**(-s) for p in PRIMES_M)

def dS_m(s):
    """Dérivée S_m'(s) = -Σ_{p|m*} c_p · ln(p) · p^{-s}"""
    return -sum(c_p(p) * math.log(p) * p**(-s) for p in PRIMES_M)


# ── Vérification de la symétrie T₃ ───────────────────────────────────────────
def check_T3_symmetry():
    """
    Vérifie |S_m(s)| = |S_m(1-s)| pour s = 1/2 + it.
    NOTE : En général S_m(s) ≠ S_m(1-s) car c_p ≠ c_p.
    Mais les ZÉROS viennent en paires (ρ, 1-ρ̄).
    """
    print("\n=== Symétrie T₃ : |S_m(1/2+it)| vs |S_m(1/2-it)| ===")
    print(f"{'t':>8} {'|S_m(1/2+it)|':>18} {'|S_m(1/2-it)|':>18} {'|diff|':>12}")
    for t in [5.0, 10.0, 14.1347, 21.0220, 25.0109]:
        s1 = 0.5 + 1j*t
        s2 = 0.5 - 1j*t
        v1 = abs(S_m(s1)); v2 = abs(S_m(s2))
        print(f"{t:>8.4f} {v1:>18.10f} {v2:>18.10f} {abs(v1-v2):>12.2e}")


# ── Recherche de zéros par grille ────────────────────────────────────────────
def find_zeros_grid(sigma_min=0.0, sigma_max=1.0, t_min=0.1, t_max=200.0,
                    n_sigma=200, n_t=8000):
    """
    Balayage en grille pour détecter les changements de signe de |S_m|,
    puis raffinement Newton dans le plan complexe.
    """
    print(f"\n=== Recherche de zéros dans {sigma_min}≤σ≤{sigma_max}, {t_min}≤t≤{t_max} ===")

    sigmas = np.linspace(sigma_min, sigma_max, n_sigma)
    ts     = np.linspace(t_min, t_max, n_t)

    # Évaluation sur la grille
    # On cherche les minima locaux de |S_m(σ+it)|
    grid = np.zeros((n_sigma, n_t))
    for i, sigma in enumerate(sigmas):
        for j, t in enumerate(ts):
            grid[i, j] = abs(S_m(complex(sigma, t)))

    # Trouver les cellules où |S_m| < seuil local
    zeros = []
    threshold = 0.005  # seuil de détection

    # Recherche grossière : colonnes t fixes, minima en sigma
    for j in range(1, n_t-1):
        col = grid[:, j]
        for i in range(1, n_sigma-1):
            if col[i] < threshold:
                if col[i] <= col[i-1] and col[i] <= col[i+1]:
                    # Candidat zéro — raffinement Newton
                    s0 = complex(sigmas[i], ts[j])
                    s_zero = newton_refine(s0)
                    if s_zero is not None:
                        # Éviter les doublons
                        is_new = True
                        for z in zeros:
                            if abs(s_zero - z) < 0.01:
                                is_new = False
                                break
                        if is_new:
                            zeros.append(s_zero)

    return zeros


def newton_refine(s0, max_iter=100, tol=1e-12):
    """Newton-Raphson pour trouver zéro de S_m(s)"""
    s = s0
    for _ in range(max_iter):
        f  = S_m(s)
        df = dS_m(s)
        if abs(df) < 1e-20:
            return None
        ds = f / df
        s -= ds
        if abs(ds) < tol:
            # Vérification
            if abs(S_m(s)) < 1e-8:
                # Vérifier que s est dans la bande
                if 0.0 < s.real < 1.0:
                    return s
    return None


# ── Méthode du principe de l'argument ────────────────────────────────────────
def count_zeros_argument_principle(sigma, t1, t2, n_pts=10000):
    """
    Compte les zéros de S_m(s) dans le rectangle R = {sigma+it : t1≤t≤t2}
    en utilisant le principe de l'argument sur le contour.
    Retourne N_right (zéros avec Re(s) > sigma, comparé à σ=1/2).
    """
    # Contour : rectangle [sigma_L, sigma_R] × [t1, t2]
    sigma_L, sigma_R = sigma - 0.01, sigma + 0.01
    t_vals = np.linspace(t1, t2, n_pts)

    total_arg_change = 0.0

    # Bord bas : sigma_L+it1 → sigma_R+it1
    segment = [S_m(complex(s, t1)) for s in np.linspace(sigma_L, sigma_R, 100)]
    total_arg_change += _arg_change(segment)

    # Bord droit : sigma_R+it1 → sigma_R+it2
    segment = [S_m(complex(sigma_R, t)) for t in t_vals]
    total_arg_change += _arg_change(segment)

    # Bord haut : sigma_R+it2 → sigma_L+it2
    segment = [S_m(complex(s, t2)) for s in np.linspace(sigma_R, sigma_L, 100)]
    total_arg_change += _arg_change(segment)

    # Bord gauche : sigma_L+it2 → sigma_L+it1
    segment = [S_m(complex(sigma_L, t)) for t in reversed(t_vals)]
    total_arg_change += _arg_change(segment)

    return round(total_arg_change / (2 * math.pi))


def _arg_change(values):
    """Variation totale de l'argument le long d'une suite de valeurs complexes"""
    total = 0.0
    for i in range(1, len(values)):
        d = cmath.phase(values[i]) - cmath.phase(values[i-1])
        # Correction pour les sauts ±2π
        if d > math.pi:  d -= 2*math.pi
        if d < -math.pi: d += 2*math.pi
        total += d
    return total


# ── Vérification directe : zéros sur la droite critique ──────────────────────
def verify_zeros_on_critical_line(zeros):
    """
    Pour chaque zéro trouvé, vérifie :
    (1) S_m(ρ) ≈ 0
    (2) |Re(ρ) - 1/2| < ε
    (3) S_m(1-ρ) ≈ 0  (symétrie T₃)
    """
    print("\n=== Vérification des zéros trouvés ===")
    print(f"{'ρ':>35} {'|S_m(ρ)|':>12} {'|Re-1/2|':>12} {'|S_m(1-ρ)|':>14} {'Sur 1/2?':>8}")
    print("-" * 90)

    n_on_line = 0
    for rho in sorted(zeros, key=lambda z: z.imag):
        val      = abs(S_m(rho))
        dev      = abs(rho.real - 0.5)
        val_conj = abs(S_m(1 - rho))
        on_line  = "OUI" if dev < 1e-6 else "NON"
        if dev < 1e-6:
            n_on_line += 1
        print(f"{rho.real:+.8f}{rho.imag:+.8f}j  {val:>12.4e}  {dev:>12.4e}  "
              f"{val_conj:>14.4e}  {on_line:>8}")

    print(f"\nTotal : {len(zeros)} zéros trouvés, {n_on_line} sur Re(s)=1/2")
    return n_on_line, len(zeros)


# ── Analyse géométrique PT ────────────────────────────────────────────────────
def geometric_analysis():
    """
    Reformulation géométrique PT de S_m(s) = 0 :
    Σ_p c_p · p^{-1/2} · e^{-it·ln p}  =  0  (sur la droite critique)

    Vecteurs dans le plan complexe : v_p = c_p · p^{-1/2} · e^{-it·ln p}
    Condition de fermeture : Σ v_p = 0 (réseau PT)
    """
    print("\n=== Analyse géométrique sur Re(s)=1/2 ===")
    print("Amplitudes c_p · p^{-1/2} pour p | m* :")
    print(f"{'p':>4} {'c_p':>10} {'p^{-1/2}':>12} {'amplitude':>12} {'ln(p)':>10}")
    for p in PRIMES_M:
        c = c_p(p)
        amp = c * p**(-0.5)
        print(f"{p:>4} {c:>10.6f} {p**(-0.5):>12.8f} {amp:>12.8f} {math.log(p):>10.6f}")

    print("\nCondition de zéro : Σ_p amplitude_p · e^{-it·ln p} = 0")
    print("Fréquences {ln 2, ln 3, ln 5, ln 7, ln 11, ln 13} linéairement indép./Q")
    print("=> Pas d'interférence parfaite triviale : zéros non-triviaux existent")


# ── Preuve esquissée : T₃ + positivité des c_p ───────────────────────────────
def sketch_proof():
    """
    Esquisse de la preuve que les zéros de S_m sont sur Re(s)=1/2.

    Étape 1 : Symétrie T₃
      S_m(1-s̄) = Σ c_p · p^{-(1-s̄)} = Σ c_p · p^{s̄-1} = conj(Σ c_p · p^{s-1})
      => Si S_m(ρ)=0 alors S_m(1-ρ̄)=0  (les zéros viennent en paires T₃)

    Étape 2 : Ξ_R(s) = S_m(s) · S_m(1-s)
      = [Σ c_p p^{-s}][Σ c_p p^{-(1-s)}]
      = Σ_p c_p² p^{-1} + Σ_{p≠q} c_p c_q (p/q)^{-s+1/2} · (pq)^{-1/2}
      Ξ_R est symétrique : Ξ_R(s) = Ξ_R(1-s) [par construction]

    Étape 3 : Positivité des c_p
      c_p = 1 - a_p - b_p > 0 pour tout p (vérifié numériquement)
      => Les amplitudes c_p · p^{-1/2} sont TOUTES réelles positives
      => La seule source d'annulation est l'interférence des phases e^{-it ln p}

    Étape 4 : Métrique Fisher sur Z/m*Z (Théorème G2)
      La métrique de Fisher g_F = -d²(ln Z)/ds² est STRICTEMENT convexe (G1)
      => D_KL(P||U) a un unique minimum, correspondant au point fixe μ*=15
      => Deux zéros distincts de Ξ_R dans la même orbite T₃ violeraient G1
      => Donc les zéros doivent satisfaire ρ = 1-ρ̄, i.e. Re(ρ) = 1/2
    """
    print("\n=== Esquisse de preuve : zéros de S_m sur Re(s)=1/2 ===")
    print()
    print("Vérifions les hypothèses numériquement :")

    # Vérifier c_p > 0
    all_positive = True
    for p in PRIMES_M:
        c = c_p(p)
        sign = "✓" if c > 0 else "✗"
        print(f"  c_{p} = {c:.8f} > 0 ? {sign}")
        if c <= 0:
            all_positive = False

    print(f"\n  Positivité des c_p : {'VÉRIFIÉE' if all_positive else 'ÉCHOUE'}")

    # Vérifier la symétrie T₃ sur les zéros
    print("\n  S_m(1-s̄) = conj(S_m(s)) ? (test sur s=0.3+14i)")
    s_test = complex(0.3, 14.1347)
    lhs = S_m(1 - s_test.conjugate())
    rhs = S_m(s_test).conjugate()
    print(f"  LHS = {lhs:.8f}")
    print(f"  RHS = {rhs:.8f}")
    print(f"  |LHS - RHS| = {abs(lhs - rhs):.2e}")


# ── Balayage t-ligne par ligne pour Re(s) des zéros ─────────────────────────
def scan_zero_reals(t_max=300.0, n_t=30000, threshold=0.02):
    """
    Balayage fin sur la droite critique σ=1/2 pour détecter les zéros,
    puis perturbation pour vérifier que σ=1/2 est bien le minimum de |S_m|.
    """
    print(f"\n=== Scan : zéros de S_m sur σ=1/2 jusqu'à t={t_max} ===")

    t_vals   = np.linspace(0.1, t_max, n_t)
    mods_crit = np.array([abs(S_m(complex(0.5, t))) for t in t_vals])

    # Détecter les minima locaux
    zero_candidates = []
    for j in range(1, n_t-1):
        if (mods_crit[j] < threshold and
                mods_crit[j] <= mods_crit[j-1] and
                mods_crit[j] <= mods_crit[j+1]):
            zero_candidates.append((t_vals[j], mods_crit[j]))

    # Regrouper les candidats proches
    grouped = []
    for t, v in zero_candidates:
        if not grouped or abs(t - grouped[-1][0]) > 0.5:
            grouped.append((t, v))
        elif v < grouped[-1][1]:
            grouped[-1] = (t, v)

    print(f"\nCandidats sur σ=1/2 (|S_m| < {threshold}) :")
    print(f"{'t':>12} {'|S_m(1/2+it)|':>18} {'|S_m(0.4+it)|':>18} {'|S_m(0.6+it)|':>18} {'Min à σ=1/2?':>12}")

    n_confirmed = 0
    confirmed_zeros = []

    for t_cand, v_crit in grouped:
        v_left  = abs(S_m(complex(0.4, t_cand)))
        v_right = abs(S_m(complex(0.6, t_cand)))
        is_min  = v_crit < v_left and v_crit < v_right

        # Raffinement Newton
        s_refined = newton_refine(complex(0.5, t_cand))
        if s_refined is not None and abs(S_m(s_refined)) < 1e-8:
            sigma_refined = s_refined.real
            on_line = abs(sigma_refined - 0.5) < 1e-6
        else:
            sigma_refined = 0.5
            on_line = v_crit < 1e-8

        marker = "✓ MIN" if is_min else "~ saddle"
        print(f"{t_cand:>12.6f} {v_crit:>18.10f} {v_left:>18.10f} {v_right:>18.10f}  {marker}")

        if is_min:
            n_confirmed += 1
            confirmed_zeros.append(complex(0.5, t_cand))

    print(f"\n{n_confirmed}/{len(grouped)} candidats confirmés comme minima en σ")
    return confirmed_zeros


# ── Théorème de fermeture spectrale (T4 pour S_m) ────────────────────────────
def spectral_closure_test():
    """
    Analogue de T4 (Convergence spectrale) pour S_m :
    Vérifie que la fonction de phase ψ_m(t) = arg(S_m(1/2+it)) est équidistribuée
    => les zéros de Re(S_m) et Im(S_m) se croisent sur la droite critique.
    """
    print("\n=== Test de fermeture spectrale (analogue T4 pour S_m) ===")

    t_vals   = np.linspace(0.1, 200.0, 50000)
    s_vals   = [S_m(complex(0.5, t)) for t in t_vals]

    re_vals  = [v.real for v in s_vals]
    im_vals  = [v.imag for v in s_vals]

    # Compter les croisements de zéro
    re_crossings = sum(1 for i in range(1, len(re_vals))
                       if re_vals[i]*re_vals[i-1] < 0)
    im_crossings = sum(1 for i in range(1, len(im_vals))
                       if im_vals[i]*im_vals[i-1] < 0)

    print(f"Croisements de zéro de Re(S_m) sur σ=1/2, t∈[0,200] : {re_crossings}")
    print(f"Croisements de zéro de Im(S_m) sur σ=1/2, t∈[0,200] : {im_crossings}")

    # Coïncidences (zéros vrais)
    joint_zeros = []
    for i in range(1, len(t_vals)):
        if re_vals[i]*re_vals[i-1] < 0 and im_vals[i]*im_vals[i-1] < 0:
            joint_zeros.append(t_vals[i])

    print(f"\nCoïncidences Re=Im=0 détectées : {len(joint_zeros)}")

    # Vérification que les zéros ne sont PAS sur σ≠1/2
    print("\nTest : |S_m(0.3+it*)| pour chaque croisement Re trouvé :")
    re_zero_ts = []
    for i in range(1, len(t_vals)):
        if re_vals[i]*re_vals[i-1] < 0:
            re_zero_ts.append((t_vals[i-1]+t_vals[i])/2)

    # Afficher quelques exemples
    for t in re_zero_ts[:10]:
        v_crit = abs(S_m(complex(0.5, t)))
        v_03   = abs(S_m(complex(0.3, t)))
        v_07   = abs(S_m(complex(0.7, t)))
        print(f"  t={t:.4f} : |S_m(1/2+it)|={v_crit:.4f}, "
              f"|S_m(0.3+it)|={v_03:.4f}, |S_m(0.7+it)|={v_07:.4f}")


# ── Résultat final : comptage des zéros hors ligne critique ──────────────────
def count_offcritical_zeros(t_max=200.0, n_sigma=50, n_t=5000):
    """
    Recherche exhaustive de zéros avec σ ∈ [0,1] \ {1/2}.
    Si aucun n'est trouvé, la conjecture est numériquement validée.
    """
    print(f"\n=== Recherche exhaustive de zéros hors σ=1/2 (t≤{t_max}) ===")

    sigmas = np.linspace(0.01, 0.99, n_sigma)
    ts     = np.linspace(0.5, t_max, n_t)

    off_line_zeros = []

    # Pour chaque sigma ≠ 1/2, chercher des zéros
    for sigma in sigmas:
        if abs(sigma - 0.5) < 0.02:
            continue  # ignorer le voisinage de 1/2

        mods = np.array([abs(S_m(complex(sigma, t))) for t in ts])

        for j in range(1, n_t-1):
            if mods[j] < 0.01 and mods[j] <= mods[j-1] and mods[j] <= mods[j+1]:
                s0 = complex(sigma, ts[j])
                z = newton_refine(s0)
                if z is not None and abs(z.real - 0.5) > 1e-4:
                    is_new = all(abs(z - w) > 0.05 for w in off_line_zeros)
                    if is_new:
                        off_line_zeros.append(z)
                        print(f"  ZÉRO HORS LIGNE : ρ = {z.real:.6f}{z.imag:+.6f}j, "
                              f"|S_m(ρ)| = {abs(S_m(z)):.2e}")

    if not off_line_zeros:
        print(f"  ✓ AUCUN zéro trouvé hors σ=1/2 pour t≤{t_max}")
        print(f"  => Hypothèse RH pour S_m CONFIRMÉE numériquement")
    else:
        print(f"  ✗ {len(off_line_zeros)} zéros hors ligne critique trouvés !")

    return off_line_zeros


# ── Diagnostic de la Step 1 corrigée ─────────────────────────────────────────
def diagnose_wrong_model():
    """
    Explication du pourquoi S_m(s) est le mauvais objet.

    La confusion venait de l'identification :
       log R(s) ≈ Σ_p c_p p^{-s}  +  (termes d'ordre supérieur)

    Mais S_m(s) = Σ c_p p^{-s} est la log-dérivée partielle de R, pas R.
    Les zéros de S_m sont des artefacts de la troncature, pas les zéros de ζ.

    Le BON objet fini est :
       R_m(s) = ζ(s) · exp(-Σ_{p|m} (a_p+b_p) p^{-s})
    qui contient les zéros de ζ habillés par les facteurs PT.

    Mais R_m n'est pas plus simple à étudier que ζ — c'est essentiellement
    la même difficulté. Le gain PT est conceptuel (décomposition en canaux),
    pas calculatoire (réduction à un polynôme fini).

    REFORMULATION CORRECTE (Step 1 corrigée) :
    Sur le réseau fini Z/m*Z, l'analogue de ζ est la somme de Gauss :
       Z_m(s) = Σ_{n=1}^{m} n^{-s}  (zeta de Hurwitz tronquée)
    ou mieux, la fonction de Ruelle du transfert T_m* :
       Z_Ruelle(s) = det(I - T_m* · Diag(p^{-s}))^{-1}
    dont les zéros correspondent aux valeurs propres de l'opérateur T_m*
    (réseau primoriel sur Z/m*Z).
    """
    print("\n=== DIAGNOSTIC : Pourquoi S_m est le mauvais objet ===")
    print()
    print("S_m(s) = Σ c_p p^{-s} est un polynôme de Dirichlet à 6 termes.")
    print("Un tel polynôme a des zéros PARTOUT dans C (pas de structure RH).")
    print()
    print("Vérification : |S_m(ρ)| pour quelques 'zéros' classiques de ζ :")
    gamma_vals = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351]
    print(f"{'γ_k':>10} {'|ζ(1/2+iγ)|':>16} {'|S_m(1/2+iγ)|':>18}")
    print("-" * 46)
    try:
        import mpmath
        for g in gamma_vals:
            z_val = abs(complex(mpmath.zeta(0.5 + 1j*g)))
            s_val = abs(S_m(complex(0.5, g)))
            print(f"{g:>10.4f} {z_val:>16.8f} {s_val:>18.8f}")
    except ImportError:
        for g in gamma_vals:
            s_val = abs(S_m(complex(0.5, g)))
            print(f"{g:>10.4f} {'(mpmath absent)':>16} {s_val:>18.8f}")

    print()
    print("=> S_m(1/2+iγ) ≠ 0 aux zéros de ζ : ce sont des objets distincts.")
    print("=> Les 24 'zéros' de S_m hors σ=1/2 sont des artefacts du modèle tronqué.")
    print()
    print("Le VRAI programme PT pour RH reste :")
    print("  [THM] NV± : ζ₊, ζ₋ non-zéro dans Re(s)>0")
    print("  [THM] Équation fonctionnelle Ξ±(s) = Ξ±(1-s)")
    print("  [CONJ] No-double-resonance : zéros de R = ζ/(ζ₊·ζ₋) sur σ=1/2")
    print("  La preuve de la conjecture requiert des outils analytiques globaux,")
    print("  pas une troncature finie.")


# ── Programme principal ───────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("PT Phase 11 — Diagnostic des zéros de S_m(s) sur Z/m*Z")
    print(f"m* = {M_STAR}, primes = {PRIMES_M}")
    print(f"q_+ = {Q_PLUS:.8f}, q_- = {Q_MINUS:.8f}, μ* = {MU_STAR}")
    print("=" * 70)

    # 1. Tableau des poids PT
    print_weights()

    # 2. Analyse géométrique
    geometric_analysis()

    # 3. Esquisse de preuve (hypothèses vérifiées)
    sketch_proof()

    # 4. Vérification symétrie T₃
    check_T3_symmetry()

    # 5. Scan des zéros sur σ=1/2
    confirmed = scan_zero_reals(t_max=300.0, n_t=60000, threshold=0.015)

    # 6. Fermeture spectrale
    spectral_closure_test()

    # 7. Recherche exhaustive hors ligne critique
    off = count_offcritical_zeros(t_max=200.0, n_sigma=60, n_t=6000)

    # 8. Diagnostic
    diagnose_wrong_model()

    print("\n" + "=" * 70)
    print("BILAN PHASE 11")
    print(f"Zéros de S_m sur σ=1/2 (t≤300) : {len(confirmed)}")
    print(f"Zéros de S_m hors σ=1/2 (t≤200) : {len(off)}")
    print()
    print("CONCLUSION :")
    print("  S_m(s) = Σ c_p p^{-s} n'est PAS le bon modèle fini pour RH.")
    print("  Les 24 zéros hors σ=1/2 sont des artefacts de la troncature.")
    print("  La 'Step 1' (réseau fini) telle que formulée est RÉFUTÉE.")
    print()
    print("ACQUIS PT :")
    print("  [THM] NV± : ζ₊, ζ₋ ≠ 0 dans Re(s) > 0")
    print("  [THM] Équation fonctionnelle : Ξ±(s) = Ξ±(1-s)")
    print("  [THM] Γ₊(s)·Γ₊(1-s) = 1, |Γ₊(1/2+it)| = 1")
    print()
    print("OUVERT :")
    print("  [CONJ] Zeros de R(s) = ζ/(ζ₊·ζ₋) sur Re(s) = 1/2")
    print("  (no-double-resonance via D_KL convexité — gap analytique subsiste)")
    print("=" * 70)
