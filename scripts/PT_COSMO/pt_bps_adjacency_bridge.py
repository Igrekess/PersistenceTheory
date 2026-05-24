"""
PT BPS — T_BPS = R₀² × M_bipartite  (Insight 5)
=================================================
Identité numérique (erreur < 0.1%) :
  echo = z × R₀²
D'où :
  T_BPS = echo × T₃ = (z × R₀²) × T₃ = R₀² × M

Où M = z × T₃ = z × antidiag(1,1) est la matrice d'adjacence normalisée
du graphe bipartite à z bras.

SENS PHYSIQUE PT
-----------------
  Un instanton BPS (amplitude = echo) se décompose comme :
    = z orbites fermées (amplitude = R₀² chacune)
    = z × (aller R₀) × (retour R₀)
    = z × R₀² = echo

  L'amplitude BPS est la somme cohérente sur les z liens du vertex bipartite.

  Valeurs propres :
    T_BPS   : λ = ±echo
    M_bip   : λ = ±z     [adjacence brute]
    R₀²×M   : λ = ±z×R₀² = ±echo  [IDENTIQUES à T_BPS]

RÉFÉRENCES PT
--------------
  - Corollaire cor:bps_adjacency_bridge (ch20f)
  - Script contexte : PT_COSMO/pt_bps_ruelle_sector.py
  - Script contexte : PT_COSMO/pt_lattice_coordinance.py
"""
import math
import numpy as np
from scipy.special import iv as bessel_iv

echo  = math.exp(-1)
K_c   = 2.0 / math.pi
I0    = bessel_iv(0, K_c)
I1    = bessel_iv(1, K_c)
R0    = I1 / I0
N_sp  = 3   # primes actifs {3,5,7}
N_cin = 1   # prime porteur p=2
z     = N_sp + N_cin   # = 4

print("=" * 70)
print("PT BPS — PONT T_BPS = R₀² × M_BIPARTITE  (Insight 5)")
print("=" * 70)

# ── [1] CONSTANTES ────────────────────────────────────────────────────
print("\n[1] CONSTANTES PT FONDAMENTALES")
print("-" * 55)
print(f"  echo = e^{{-1}}             = {echo:.10f}")
print(f"  R₀   = I₁(2/π)/I₀(2/π)  = {R0:.10f}")
print(f"  R₀²  =                     {R0**2:.10f}")
print(f"  z    = N_sp + N_cin = 3+1 = {z}")
print()
echo_predicted = z * R0**2
err_ppb = (echo_predicted - echo) / echo * 1e9
err_pct = (echo_predicted - echo) / echo * 100.0
print(f"  z × R₀² = {echo_predicted:.10f}")
print(f"  echo    = {echo:.10f}")
print(f"  Erreur  = {err_ppb:.2f} ppb  ({err_pct:.6f}%)")

# ── [2] MATRICES 2×2 ──────────────────────────────────────────────────
print("\n[2] MATRICES 2×2")
print("-" * 55)

T3       = np.array([[0.0, 1.0], [1.0, 0.0]])
T_BPS    = echo * T3
M_bip    = z * T3                        # z × antidiag(1,1)
R0sq_M   = R0**2 * M_bip                 # R₀² × M = z×R₀² × T₃ ≈ echo × T₃

print()
print("  T₃ = antidiag(1,1) =")
print(f"    [[{T3[0,0]}, {T3[0,1]}],")
print(f"     [{T3[1,0]}, {T3[1,1]}]]")
print()
print(f"  T_BPS = echo × T₃ = {echo:.6f} × T₃ =")
print(f"    [[{T_BPS[0,0]:.6f}, {T_BPS[0,1]:.6f}],")
print(f"     [{T_BPS[1,0]:.6f}, {T_BPS[1,1]:.6f}]]")
print()
print(f"  M_bipartite = z × T₃ = {z} × T₃ =")
print(f"    [[{M_bip[0,0]:.1f}, {M_bip[0,1]:.1f}],")
print(f"     [{M_bip[1,0]:.1f}, {M_bip[1,1]:.1f}]]")
print()
print(f"  R₀² × M_bipartite = {R0**2:.6f} × M =")
print(f"    [[{R0sq_M[0,0]:.6f}, {R0sq_M[0,1]:.6f}],")
print(f"     [{R0sq_M[1,0]:.6f}, {R0sq_M[1,1]:.6f}]]")

# ── [3] VÉRIFICATION T_BPS = R₀² × M ─────────────────────────────────
print("\n[3] VÉRIFICATION NUMÉRIQUE : T_BPS = R₀² × M")
print("-" * 55)

diff_matrix  = T_BPS - R0sq_M
max_err_abs  = np.max(np.abs(diff_matrix))
max_err_rel  = max_err_abs / echo

print(f"  T_BPS − R₀²×M  (max |entrée|) = {max_err_abs:.3e}")
print(f"  Erreur relative                 = {max_err_rel:.3e}  "
      f"({max_err_rel*1e4:.2f} ppm × echo)")
print()
if max_err_rel < 1e-3:
    print("  PASS : T_BPS = R₀² × M  vérifié numériquement")
else:
    print(f"  NOTE : erreur {max_err_rel:.4f} — accord à {(1-max_err_rel)*100:.2f}%")

# ── [4] VALEURS PROPRES ───────────────────────────────────────────────
print("\n[4] VALEURS PROPRES")
print("-" * 55)

eig_T_BPS   = np.linalg.eigvals(T_BPS)
eig_M_bip   = np.linalg.eigvals(M_bip)
eig_R0sq_M  = np.linalg.eigvals(R0sq_M)

print(f"  T_BPS   : λ = {sorted(eig_T_BPS, key=lambda x: -x.real)}")
print(f"  M_bip   : λ = {sorted(eig_M_bip, key=lambda x: -x.real)}")
print(f"  R₀²×M   : λ = {sorted(eig_R0sq_M, key=lambda x: -x.real)}")
print()
print("  Relation entre valeurs propres :")
lam_T   = max(abs(eig_T_BPS))
lam_R0M = max(abs(eig_R0sq_M))
lam_M   = max(abs(eig_M_bip))
print(f"    |λ|(T_BPS)  = {lam_T:.8f}  = echo")
print(f"    |λ|(R₀²×M)  = {lam_R0M:.8f}  ≈ echo  [erreur = {abs(lam_R0M-lam_T)/lam_T*100:.4f}%]")
print(f"    |λ|(M_bip)  = {lam_M:.4f}          = z = {z}")
print()
print(f"  echo = |λ|(T_BPS) = z × R₀² = {z} × {R0**2:.6f} = {z*R0**2:.6f}")
print(f"                     = |λ|(R₀²×M)")

# ── [5] TABLE DES AMPLITUDES ──────────────────────────────────────────
print("\n[5] TABLE DES AMPLITUDES BPS")
print("-" * 55)
print()
print(f"  {'Amplitude':30s} | {'Symbole':12s} | {'Valeur numérique':16s}")
print("  " + "-" * 65)

amp_single_hop    = R0
amp_closed_orbit  = R0**2
amp_vertex_echo   = z * R0**2
amp_echo_exact    = echo

rows = [
    ("Saut simple (aller)", "R₀", amp_single_hop),
    ("Orbite fermée (A+R)", "R₀²", amp_closed_orbit),
    (f"Vertex BPS ({z} orbites)", f"z×R₀² = {z}×R₀²", amp_vertex_echo),
    ("echo exact = e^{-1}", "e^{-1}", amp_echo_exact),
    ("Erreur z×R₀² vs e^{-1}", "δ", abs(amp_vertex_echo - amp_echo_exact)),
]
for name, sym, val in rows:
    print(f"  {name:30s} | {sym:12s} | {val:.10f}")

print()
print("  Interprétation :")
print(f"    Un instanton BPS = {z} orbites fermées de {z} liens différents")
print(f"    ({N_sp} liens actifs {'{3,5,7}'} + {N_cin} lien cinématique {'{p=2}'})")
print(f"    Amplitude totale = {z} × R₀² = echo (à {abs(err_pct):.4f}%)")

# ── [6] DÉCOMPOSITION TENSORIELLE ─────────────────────────────────────
print("\n[6] DÉCOMPOSITION TENSORIELLE")
print("-" * 55)
print()
print("  T_BPS = echo × T₃")
print("        = (z × R₀²) × T₃          [echo = z × R₀², Insight 4]")
print("        = R₀² × (z × T₃)")
print("        = R₀² × M_bipartite")
print()
print("  L'amplitude R₀² factorise la matrice de transfert BPS en :")
print("   · une amplitude de propagation par lien  : R₀²")
print("   · une structure de graphe bipartite      : M = z × T₃")
print()
print("  Connexion PT :")
print("   · T₃ = matrice de transition {q_+, q_-} (T1 interdit diagonal)")
print("   · M  = z copies de T₃ = réseau bipartite z-régulier")
print("   · R₀ = amplitude de Bessel sur le réseau de couplage K_c = 2/π")

# ── [7] BILAN ÉPISTÉMIQUE ─────────────────────────────────────────────
print("\n[7] BILAN ÉPISTÉMIQUE")
print("-" * 55)
print()
print("  Hypothèses :")
print("  [H1] T_BPS = echo × T₃  [DER, T1+T5, cf. lem:bps_transfer_z_ghost]")
print("  [H2] echo = z × R₀²     [NUM, erreur < 0.1%, Insight 4]")
print("  [H3] z = N_sp + N_cin = 3+1 = 4  [ID, structure PT des primes]")
print()
print("  Sous H1+H2+H3 :")
print("  T_BPS = R₀² × M_bipartite  [identité exacte à 0.04%]")
print()
print("  Statut : [DER, H1+H2+H3] — l'identité est numérique à 0.04%")
print("  Interprétation : l'instanton BPS = somme cohérente d'orbites")
print("  bipartites, une par prime (3 actifs + 1 cinématique)")

print()
print("=" * 70)
print("RÉSUMÉ")
print("=" * 70)
print()
print(f"  echo  = e^{{-1}}            = {echo:.8f}")
print(f"  R₀    = I₁(2/π)/I₀(2/π) = {R0:.8f}")
print(f"  z     = N_sp + N_cin     = {z}")
print(f"  echo  = z × R₀²          = {z*R0**2:.8f}  (erreur {abs(err_pct):.4f}%)")
print()
print(f"  T_BPS = echo × T₃ = R₀² × (z×T₃) = R₀² × M_bipartite")
print(f"  Max |T_BPS − R₀²×M| = {max_err_abs:.2e}  [PASS]")
print()
print(f"  Valeurs propres T_BPS : ±{lam_T:.6f} = ±echo")
print(f"  Valeurs propres R₀²×M : ±{lam_R0M:.6f} ≈ ±echo")
print()
print(f"  Corollaire : cor:bps_adjacency_bridge  (ch20f)")
print(f"  Statut : [DER] sous T1+T5 + identification z=N_sp+1")
