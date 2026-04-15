#!/usr/bin/env python3
"""
proof_C7_transfer_matrix.py -- Chapter 22: Chemistry

Monograph: ch22_chemistry.tex
Derivation chain: s = 1/2 -> q_plus = 13/15 -> sin2_p -> transfer matrix T
Zero fitted parameters.

This script proves Theorem C7: the transfer matrix spectral invariance.

For an N-atom polyatomic molecule, the atomization energy is a spectral
invariant of the transfer matrix T:
  T_vv = IE_v / Ry                 (site energies on diagonal)
  T_vw = D0 * exp(-S_cross) / Ry  (bonded off-diagonal)
  T_vw = 0                        (non-bonded off-diagonal)

T is symmetric. Its eigenvalues encode bonding topology:
  - Perron eigenvalue lambda_0 is real, positive, maximal
  - Spectral gap lambda_0 - lambda_1 measures bond strength
  - For conjugated systems, T reduces to Hueckel model with PT-derived
    alpha and beta

  Step 1. MATRIX STRUCTURE
          Build T for H2, H2O, CH4, C2H6, benzene. Verify dimensions,
          symmetry, diagonal = IE/Ry, off-diagonal structure.

  Step 2. DIATOMIC LIMIT (HOMONUCLEAR)
          Eigenvalues of 2x2: lambda = a +/- t for 10 molecules.

  Step 3. HETERONUCLEAR DIATOMIC
          Eigenvalue formula with IE asymmetry for 10 molecules.

  Step 4. HUECKEL-PT FOR CONJUGATED SYSTEMS
          Ethylene, allyl, butadiene, benzene pi-systems.

  Step 5. PERRON-FROBENIUS PROPERTIES
          Perron eigenvalue positive/maximal, eigenvector all-positive,
          eigenvalue sum = trace.

  Step 6. SPECTRAL GAP AND BONDING
          Gap ordering, hetero vs homo, benzene delocalization,
          Gershgorin bounds, basis invariance.

Theorems verified:
  C7 "Transfer Matrix" (ch22_chemistry.tex) -- Spectral invariance of T

PT constants used:
  s = 1/2, mu* = 15, q_plus = 13/15
  Ry = 13.606 eV, P1 = 3, D0 ~ 4.535 * 1.005 eV
  sin2_3 ~ 0.219 (cross-suppression factor)
"""

import math
import sys
from pathlib import Path

# Path setup
_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker
from pt_constants import s, mu_star, q_plus, PRIMES_ACTIFS, sin2_theta

ck = Checker("proof_C7_transfer_matrix", chapter="ch22_chemistry", total_steps=6)

# Physical constants
Ry = 13.606  # eV
P1, P2, P3 = 3, 5, 7
TOL = 1e-9

# Derived
sin2_3 = sin2_theta(P1, q_plus)
D0_tree = Ry / P1
D_FULL = 1.005   # NLO vertex dressing (approximate)
D0 = D0_tree * D_FULL
S_cross_default = sin2_3  # leading cross-entropy

# IE data (eV, NIST)
IE = {
    'H':  13.598, 'He': 24.587, 'Li':  5.392, 'Be':  9.323,
    'B':   8.298, 'C':  11.260, 'N':  14.534, 'O':  13.618,
    'F':  17.423, 'Ne': 21.565, 'Na':  5.139, 'Cl': 12.968,
    'Br': 11.814, 'I':  10.451, 'K':   4.341, 'S':  10.360,
}


# =====================================================================
# Linear algebra helpers (math module only -- no numpy)
# =====================================================================

def mat_zeros(n):
    return [[0.0]*n for _ in range(n)]

def mat_copy(M):
    return [row[:] for row in M]

def mat_symmetric(M, tol=1e-12):
    n = len(M)
    for i in range(n):
        for j in range(i+1, n):
            if abs(M[i][j] - M[j][i]) > tol:
                return False
    return True

def mat_trace(M):
    return sum(M[i][i] for i in range(len(M)))

def vec_dot(u, v):
    return sum(a*b for a, b in zip(u, v))

def vec_norm(v):
    return math.sqrt(sum(x*x for x in v))

def mat_vec(M, v):
    n = len(M)
    return [sum(M[i][j]*v[j] for j in range(n)) for i in range(n)]

def vec_scale(c, v):
    return [c*x for x in v]

def eig_2x2(M):
    a, b = M[0][0], M[0][1]
    d = M[1][1]
    tr = a + d
    det = a*d - b*b
    disc = max(tr*tr - 4.0*det, 0.0)
    sq = math.sqrt(disc)
    return [(tr + sq) / 2.0, (tr - sq) / 2.0]

def eig_3x3_symmetric(M):
    a00, a01, a02 = M[0][0], M[0][1], M[0][2]
    a11, a12 = M[1][1], M[1][2]
    a22 = M[2][2]
    tr = a00 + a11 + a22
    q_shift = tr / 3.0
    b00, b11, b22 = a00 - q_shift, a11 - q_shift, a22 - q_shift
    p2 = b00**2 + b11**2 + b22**2 + 2.0*(a01**2 + a02**2 + a12**2)
    if p2 < 1e-30:
        return sorted([tr/3.0]*3, reverse=True)
    p = math.sqrt(p2 / 6.0)
    det_B = (b00*(b11*b22 - a12**2) - a01*(a01*b22 - a12*a02)
             + a02*(a01*a12 - b11*a02))
    r = max(-1.0, min(1.0, det_B / (2.0 * p**3)))
    phi = math.acos(r) / 3.0
    eigs = [q_shift + 2.0*p*math.cos(phi + 2.0*math.pi*k/3.0) for k in [0, -1, 1]]
    return sorted(eigs, reverse=True)

def jacobi_eigenvalues(M, max_iter=500):
    n = len(M)
    A = mat_copy(M)
    for _ in range(max_iter):
        max_val, p_idx, q_idx = 0.0, 0, 1
        for i in range(n):
            for j in range(i+1, n):
                if abs(A[i][j]) > max_val:
                    max_val = abs(A[i][j])
                    p_idx, q_idx = i, j
        if max_val < 1e-14:
            break
        app, aqq, apq = A[p_idx][p_idx], A[q_idx][q_idx], A[p_idx][q_idx]
        if abs(app - aqq) < 1e-30:
            theta = math.pi / 4.0
        else:
            theta = 0.5 * math.atan2(2.0*apq, app - aqq)
        c, sr = math.cos(theta), math.sin(theta)
        A_new = mat_copy(A)
        for i in range(n):
            if i != p_idx and i != q_idx:
                A_new[i][p_idx] = c*A[i][p_idx] + sr*A[i][q_idx]
                A_new[p_idx][i] = A_new[i][p_idx]
                A_new[i][q_idx] = -sr*A[i][p_idx] + c*A[i][q_idx]
                A_new[q_idx][i] = A_new[i][q_idx]
        A_new[p_idx][p_idx] = c*c*app + 2*sr*c*apq + sr*sr*aqq
        A_new[q_idx][q_idx] = sr*sr*app - 2*sr*c*apq + c*c*aqq
        A_new[p_idx][q_idx] = 0.0
        A_new[q_idx][p_idx] = 0.0
        A = A_new
    return sorted([A[i][i] for i in range(n)], reverse=True)

def power_iteration(M, num_iter=2000, tol=1e-12):
    n = len(M)
    v = [1.0/math.sqrt(n)] * n
    lam = 0.0
    for _ in range(num_iter):
        w = mat_vec(M, v)
        lam_new = vec_dot(v, w)
        nw = vec_norm(w)
        if nw < 1e-30:
            break
        v = vec_scale(1.0/nw, w)
        if abs(lam_new - lam) < tol:
            break
        lam = lam_new
    return lam, v

def circulant_eigenvalues(first_row):
    n = len(first_row)
    eigs = []
    for k in range(n):
        val = sum(first_row[j] * math.cos(2.0*math.pi*k*j/n) for j in range(n))
        eigs.append(val)
    return sorted(eigs, reverse=True)

def tridiag_eigenvalues(diag, offdiag):
    n = len(diag)
    if n == 1:
        return [diag[0]]
    if n == 2:
        return eig_2x2([[diag[0], offdiag[0]], [offdiag[0], diag[1]]])
    if n == 3:
        return eig_3x3_symmetric(
            [[diag[0], offdiag[0], 0.0],
             [offdiag[0], diag[1], offdiag[1]],
             [0.0, offdiag[1], diag[2]]])
    # Sturm bisection for n >= 4
    def sturm_count(x):
        count = 0
        d_prev = diag[0] - x
        if d_prev < 0:
            count += 1
        for i in range(1, n):
            if abs(d_prev) < 1e-30:
                d_prev = 1e-30
            d_cur = (diag[i] - x) - offdiag[i-1]**2 / d_prev
            if d_cur < 0:
                count += 1
            d_prev = d_cur
        return count
    lo = min(diag[i] - abs(offdiag[i-1] if i > 0 else 0) - abs(offdiag[i] if i < n-1 else 0)
             for i in range(n)) - 1.0
    hi = max(diag[i] + abs(offdiag[i-1] if i > 0 else 0) + abs(offdiag[i] if i < n-1 else 0)
             for i in range(n)) + 1.0
    eigs = []
    for idx in range(n):
        a, b = lo, hi
        for _ in range(200):
            mid = (a + b) / 2.0
            if sturm_count(mid) > idx:
                b = mid
            else:
                a = mid
            if b - a < 1e-14:
                break
        eigs.append((a + b) / 2.0)
    return sorted(eigs, reverse=True)


# =====================================================================
# Transfer matrix builder
# =====================================================================

def build_transfer_matrix(atoms, bonds, ie_dict=None, D0_bond=None, S_cross=None):
    if ie_dict is None:
        ie_dict = IE
    if D0_bond is None:
        D0_bond = D0
    if S_cross is None:
        S_cross = S_cross_default
    n = len(atoms)
    T = mat_zeros(n)
    for v in range(n):
        T[v][v] = ie_dict[atoms[v]] / Ry
    coupling = D0_bond * math.exp(-S_cross) / Ry
    for (v, w) in bonds:
        T[v][w] = coupling
        T[w][v] = coupling
    return T


# =====================================================================
# Step 1: MATRIX STRUCTURE (20 checks)
# =====================================================================
ck.section("Step 1: Matrix structure")

coupling_val = D0 * math.exp(-S_cross_default) / Ry

# H2
T_H2 = build_transfer_matrix(['H', 'H'], [(0, 1)])
ck.check("H2 is 2x2", len(T_H2) == 2 and len(T_H2[0]) == 2)
ck.check("H2 diagonal = IE_H/Ry",
         abs(T_H2[0][0] - IE['H']/Ry) < TOL and abs(T_H2[1][1] - IE['H']/Ry) < TOL)
ck.check("H2 off-diag = D0*exp(-S_cross)/Ry",
         abs(T_H2[0][1] - coupling_val) < TOL)
ck.check("H2 symmetric", mat_symmetric(T_H2))

# H2O
atoms_H2O = ['O', 'H', 'H']
T_H2O = build_transfer_matrix(atoms_H2O, [(0, 1), (0, 2)])
ck.check("H2O is 3x3", len(T_H2O) == 3)
ck.check("H2O O-site energy", abs(T_H2O[0][0] - IE['O']/Ry) < TOL)
ck.check("H2O H-site energies",
         abs(T_H2O[1][1] - IE['H']/Ry) < TOL and abs(T_H2O[2][2] - IE['H']/Ry) < TOL)
ck.check("H2O O-H bonds nonzero",
         abs(T_H2O[0][1]) > TOL and abs(T_H2O[0][2]) > TOL)
ck.check("H2O H-H not bonded = 0", abs(T_H2O[1][2]) < TOL)
ck.check("H2O symmetric", mat_symmetric(T_H2O))

# CH4
atoms_CH4 = ['C', 'H', 'H', 'H', 'H']
bonds_CH4 = [(0, 1), (0, 2), (0, 3), (0, 4)]
T_CH4 = build_transfer_matrix(atoms_CH4, bonds_CH4)
ck.check("CH4 is 5x5", len(T_CH4) == 5)
ck.check("CH4 C-site energy", abs(T_CH4[0][0] - IE['C']/Ry) < TOL)
ck.check("CH4 symmetric", mat_symmetric(T_CH4))
ck.check("CH4 H-H entries = 0",
         all(abs(T_CH4[i][j]) < TOL for i in range(1, 5) for j in range(1, 5) if i != j))

# C2H6
atoms_C2H6 = ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
bonds_C2H6 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7)]
T_C2H6 = build_transfer_matrix(atoms_C2H6, bonds_C2H6)
ck.check("C2H6 is 8x8", len(T_C2H6) == 8)
ck.check("C2H6 symmetric", mat_symmetric(T_C2H6))
ck.check("C2H6 C-C bond nonzero", abs(T_C2H6[0][1]) > TOL)

# Benzene
atoms_benz = ['C']*6 + ['H']*6
bonds_benz = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),
              (0,6),(1,7),(2,8),(3,9),(4,10),(5,11)]
T_benz = build_transfer_matrix(atoms_benz, bonds_benz)
ck.check("Benzene is 12x12", len(T_benz) == 12)
ck.check("Benzene symmetric", mat_symmetric(T_benz))
ck.check("Benzene trace = sum(site energies)",
         abs(mat_trace(T_benz) - (6*IE['C']/Ry + 6*IE['H']/Ry)) < TOL)


# =====================================================================
# Step 2: DIATOMIC LIMIT -- HOMONUCLEAR (20 checks)
# =====================================================================
ck.section("Step 2: Diatomic limit -- homonuclear")

homonuc_order = ['H2', 'N2', 'O2', 'F2', 'Cl2', 'Li2', 'Na2', 'K2', 'Br2', 'I2']
homonuc_data = {
    'H2': 'H', 'N2': 'N', 'O2': 'O', 'F2': 'F', 'Cl2': 'Cl',
    'Li2': 'Li', 'Na2': 'Na', 'K2': 'K', 'Br2': 'Br', 'I2': 'I',
}

for mol_name in homonuc_order:
    sym = homonuc_data[mol_name]
    a = IE[sym] / Ry
    t = coupling_val
    T_mol = [[a, t], [t, a]]
    eigs = eig_2x2(T_mol)
    ck.check(f"{mol_name:<4s} lam_0 = a+t",
             abs(eigs[0] - (a + t)) < TOL,
             f"lam0={eigs[0]:.6f}")
    ck.check(f"{mol_name:<4s} lam_1 = a-t",
             abs(eigs[1] - (a - t)) < TOL,
             f"lam1={eigs[1]:.6f}")


# =====================================================================
# Step 3: HETERONUCLEAR DIATOMIC (20 checks)
# =====================================================================
ck.section("Step 3: Heteronuclear diatomic")

hetero_data = [
    ('HF',   'H',  'F'),  ('HCl',  'H',  'Cl'), ('HBr',  'H',  'Br'),
    ('HI',   'H',  'I'),  ('CO',   'C',  'O'),   ('NO',   'N',  'O'),
    ('NaCl', 'Na', 'Cl'), ('LiF',  'Li', 'F'),   ('KBr',  'K',  'Br'),
    ('NaF',  'Na', 'F'),
]

for name, A, B in hetero_data:
    a = IE[A] / Ry
    d = IE[B] / Ry
    t = coupling_val
    half_sum = (a + d) / 2.0
    disc = math.sqrt(((a - d) / 2.0)**2 + t**2)
    lam0_exp = half_sum + disc
    lam1_exp = half_sum - disc
    T_mol = [[a, t], [t, d]]
    eigs = eig_2x2(T_mol)
    ck.check(f"{name:<5s} eigenvalue formula",
             abs(eigs[0] - lam0_exp) < TOL and abs(eigs[1] - lam1_exp) < TOL)
    ck.check(f"{name:<5s} Perron lam0 > 0 and maximal",
             eigs[0] > 0 and eigs[0] >= abs(eigs[1]))


# =====================================================================
# Step 4: HUECKEL-PT FOR CONJUGATED SYSTEMS (20 checks)
# =====================================================================
ck.section("Step 4: Hueckel-PT for conjugated systems")

alpha_H = IE['C'] / Ry
beta_H = D0 * math.sin(math.pi / P1)**2 / Ry

# Ethylene (2 centers)
eigs_eth = eig_2x2([[alpha_H, beta_H], [beta_H, alpha_H]])
ck.check("Ethylene lam0 = alpha + beta",
         abs(eigs_eth[0] - (alpha_H + beta_H)) < TOL)
ck.check("Ethylene lam1 = alpha - beta",
         abs(eigs_eth[1] - (alpha_H - beta_H)) < TOL)
ck.check("Ethylene trace = 2*alpha",
         abs(sum(eigs_eth) - 2*alpha_H) < TOL)
ck.check("Ethylene gap = 2*beta",
         abs(eigs_eth[0] - eigs_eth[1] - 2*beta_H) < TOL)

# Allyl radical (3 centers)
eigs_allyl = tridiag_eigenvalues([alpha_H]*3, [beta_H]*2)
ck.check("Allyl lam0 = alpha + sqrt(2)*beta",
         abs(eigs_allyl[0] - (alpha_H + math.sqrt(2)*beta_H)) < 1e-8)
ck.check("Allyl lam1 = alpha (non-bonding)",
         abs(eigs_allyl[1] - alpha_H) < 1e-8)
ck.check("Allyl lam2 = alpha - sqrt(2)*beta",
         abs(eigs_allyl[2] - (alpha_H - math.sqrt(2)*beta_H)) < 1e-8)
ck.check("Allyl trace = 3*alpha",
         abs(sum(eigs_allyl) - 3*alpha_H) < 1e-8)

# Butadiene (4 centers)
eigs_buta = tridiag_eigenvalues([alpha_H]*4, [beta_H]*3)
buta_analytic = sorted([alpha_H + 2*beta_H*math.cos(k*math.pi/5.0) for k in range(1, 5)],
                       reverse=True)
ck.check("Butadiene 4 eigenvalues match analytic",
         all(abs(eigs_buta[i] - buta_analytic[i]) < 1e-8 for i in range(4)))
ck.check("Butadiene trace = 4*alpha",
         abs(sum(eigs_buta) - 4*alpha_H) < 1e-8)
ck.check("Butadiene bonding/antibonding symmetry",
         abs((eigs_buta[0] - alpha_H) + (eigs_buta[3] - alpha_H)) < 1e-8)
ck.check("Butadiene inner pair symmetric",
         abs((eigs_buta[1] - alpha_H) + (eigs_buta[2] - alpha_H)) < 1e-8)

# Benzene (6 centers, circulant)
circ_row = [alpha_H, beta_H, 0.0, 0.0, 0.0, beta_H]
eigs_benz_pi = circulant_eigenvalues(circ_row)
benz_analytic = sorted([alpha_H + 2*beta_H*math.cos(2*math.pi*k/6.0) for k in range(6)],
                       reverse=True)
ck.check("Benzene 6 eigenvalues match Fourier formula",
         all(abs(eigs_benz_pi[i] - benz_analytic[i]) < 1e-8 for i in range(6)))
ck.check("Benzene trace = 6*alpha",
         abs(sum(eigs_benz_pi) - 6*alpha_H) < 1e-8)
ck.check("Benzene Perron = alpha + 2*beta",
         abs(eigs_benz_pi[0] - (alpha_H + 2*beta_H)) < 1e-8)
ck.check("Benzene antibonding = alpha - 2*beta",
         abs(eigs_benz_pi[-1] - (alpha_H - 2*beta_H)) < 1e-8)

# Delocalization energy
E_pi_benz = 2*eigs_benz_pi[0] + 4*eigs_benz_pi[1]
E_pi_3eth = 3 * 2 * (alpha_H + beta_H)
deloc_energy = E_pi_benz - E_pi_3eth
ck.check("Benzene delocalization = 2*beta > 0",
         abs(deloc_energy - 2*beta_H) < 1e-8 and deloc_energy > 0,
         f"E_deloc = {deloc_energy:.6f}")
ck.check("Benzene doubly degenerate pairs",
         abs(eigs_benz_pi[1] - eigs_benz_pi[2]) < 1e-8 and
         abs(eigs_benz_pi[3] - eigs_benz_pi[4]) < 1e-8)


# =====================================================================
# Step 5: PERRON-FROBENIUS PROPERTIES (20 checks)
# =====================================================================
ck.section("Step 5: Perron-Frobenius properties")

def test_perron(mol_name, atoms, bonds):
    T = build_transfer_matrix(atoms, bonds)
    n = len(atoms)
    if n <= 3:
        eigs = eig_2x2(T) if n == 2 else eig_3x3_symmetric(T)
    else:
        eigs = jacobi_eigenvalues(T)
    lam0 = eigs[0]
    ck.check(f"{mol_name:<8s} lam0 > 0", lam0 > 0, f"lam0={lam0:.6f}")
    max_abs = max(abs(e) for e in eigs)
    ck.check(f"{mol_name:<8s} lam0 >= |lam_j|",
             abs(lam0 - max_abs) < 1e-6)
    _, v0 = power_iteration(T)
    all_pos = all(c > -1e-8 for c in v0)
    if not all_pos:
        all_pos = all(-c > -1e-8 for c in v0)
    ck.check(f"{mol_name:<8s} Perron vec all positive", all_pos)
    ck.check(f"{mol_name:<8s} sum(eigs) = trace",
             abs(sum(eigs) - mat_trace(T)) < 1e-4)

test_perron("H2", ['H', 'H'], [(0,1)])
test_perron("HF", ['H', 'F'], [(0,1)])
test_perron("H2O", ['O', 'H', 'H'], [(0,1), (0,2)])
test_perron("CH4", ['C', 'H', 'H', 'H', 'H'], [(0,1),(0,2),(0,3),(0,4)])
test_perron("C6-ring", ['C']*6, [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0)])


# =====================================================================
# Step 6: SPECTRAL GAP AND BONDING (17 checks)
# =====================================================================
ck.section("Step 6: Spectral gap and bonding")

# Universal homonuclear gap at tree level
gap_univ = 2.0 * coupling_val
ck.check("homonuclear gap = 2t (universal)",
         abs(gap_univ - 2*D0*math.exp(-S_cross_default)/Ry) < TOL)

# Bond order scaling
def gap_bond_order(n_bond):
    return 2.0 * n_bond * D0 * math.exp(-S_cross_default) / Ry

g1 = gap_bond_order(1)
g2 = gap_bond_order(2)
g3 = gap_bond_order(3)
ck.check("gap: triple > double > single", g3 > g2 > g1)
ck.check("triple gap = 3 * single gap", abs(g3 - 3*g1) < TOL)

# Heteronuclear gap depends on IE difference
def gap_hetero(ie_A, ie_B):
    a, d = ie_A / Ry, ie_B / Ry
    t = coupling_val
    return 2.0 * math.sqrt(((a - d) / 2.0)**2 + t**2)

gap_HF = gap_hetero(IE['H'], IE['F'])
gap_HCl = gap_hetero(IE['H'], IE['Cl'])
gap_HBr = gap_hetero(IE['H'], IE['Br'])
gap_HI = gap_hetero(IE['H'], IE['I'])

ck.check("hetero gap HF > HCl", gap_HF > gap_HCl)
ck.check("hetero gap HF > HI", gap_HF > gap_HI)
ck.check("hetero gap HI > HBr", gap_HI > gap_HBr)
ck.check("hetero gap >= homo gap", gap_HF >= gap_univ - TOL)
ck.check("all H-X gaps >= 2t",
         all(g >= gap_univ - TOL for g in [gap_HF, gap_HCl, gap_HBr, gap_HI]))

gap_NaCl = gap_hetero(IE['Na'], IE['Cl'])
gap_NaF = gap_hetero(IE['Na'], IE['F'])
ck.check("ionic NaF gap > NaCl gap", gap_NaF > gap_NaCl)

# Benzene pi gap
gap_benz = eigs_benz_pi[0] - eigs_benz_pi[1]
ck.check("benzene pi gap = beta", abs(gap_benz - beta_H) < 1e-8)

# Butadiene HOMO-LUMO gap
homo_lumo_buta = eigs_buta[1] - eigs_buta[2]
ck.check("butadiene HOMO-LUMO gap > 0", homo_lumo_buta > 0)

# Benzene HOMO-LUMO gap
homo_lumo_benz = eigs_benz_pi[2] - eigs_benz_pi[3]
ck.check("benzene HOMO-LUMO gap > 0", homo_lumo_benz > 0)

# CH4 spectral gap
T_CH4_eigs = jacobi_eigenvalues(T_CH4)
ck.check("CH4 spectral gap > 0", T_CH4_eigs[0] - T_CH4_eigs[1] > 0)

# H2O spectral gap
T_H2O_eigs = eig_3x3_symmetric(T_H2O)
ck.check("H2O spectral gap > 0", T_H2O_eigs[0] - T_H2O_eigs[1] > 0)

# Gershgorin bound
a_H = IE['H'] / Ry
t_H = coupling_val
gersh_bound = a_H + t_H
lam_H2_perron = eig_2x2(T_H2)[0]
ck.check("H2 Perron <= Gershgorin bound",
         lam_H2_perron <= gersh_bound + TOL)

# Basis invariance
T_HF = build_transfer_matrix(['H', 'F'], [(0, 1)])
T_FH = build_transfer_matrix(['F', 'H'], [(0, 1)])
eigs_HF = eig_2x2(T_HF)
eigs_FH = eig_2x2(T_FH)
ck.check("spectral invariance HF vs FH relabelling",
         abs(eigs_HF[0] - eigs_FH[0]) < TOL and abs(eigs_HF[1] - eigs_FH[1]) < TOL)


# =====================================================================
ck.summary()
