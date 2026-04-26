"""
bobeth_181_search.py -- Phase 2+ Chantier B

Search for a PT-structural decomposition of the Bobeth 2004 NNLO
prefactor 181.42 in the 3-loop mixing U_{s,92}^{(2)}(mu_b, mu_0) of
Q_2 into Q_9 for b -> s ell ell :

    U_{s,92}^{(2)} = 181.42 * e^(-5.1901 * eta) * eta^(1.5212)
    where eta = alpha_s(mu_0) / alpha_s(mu_b) ~ 2

The two side exponents are PT-structurally identified :
    5.1901  ~  C_F * beta_0 / 2  =  (4/3) * (23/3) / 2  = 5.111   (1.5%)
    1.5212  ~  gamma_3 + gamma_5 + gamma_7 + gamma_11 - 1          (0.26%)

The prefactor 181.42 is the open puzzle.  Best match in current text :
    mu* * 12 + 1 = 181                (0.2%, but "12" has no obvious PT meaning)

This script tests candidates systematically.
"""

from __future__ import annotations

from fractions import Fraction
from math import pi, log, exp, sqrt, gamma as gamma_fn

# =============================================================
# Target
# =============================================================

TARGET = 181.42
TARGET_ERR = 0.5        # Bobeth 2004 quotes precision ~0.5

# =============================================================
# PT constants
# =============================================================

MU_STAR = 15
N_C = 3
N_GEN = 3
N_S = 3                 # spatial dimensions
S = 0.5                 # s = 1/2
C_F = 4.0/3.0           # color factor
BETA_0 = 23.0/3.0       # QCD 1-loop
BETA_1 = 116.0/3.0      # QCD 2-loop
N_F = 5                 # active quark flavours at m_b
PRIMORIAL_3 = 30        # 2*mu* = 30 = 2*3*5
FACTORIAL_5 = 120       # mu* * 2^N_s = 120
ALPHA_EM_INV = 137.036
ALPHA_S_MB = 0.2124     # ~1-loop from alpha_s(m_Z) = 0.118

# gamma_p(mu* = 15) exact computations
def gamma_p(p, mu=15):
    q = 1 - 2/mu
    d = (1 - q**p)/p
    num = 4 * p * (q**(p-1)) * (1 - d)
    den = mu * (1 - q**p) * (2 - d)
    return num / den

def sin2_theta_p(p, mu=15):
    q = 1 - 2/mu
    d = (1 - q**p)/p
    return d * (2 - d)

GAMMA_3 = gamma_p(3)
GAMMA_5 = gamma_p(5)
GAMMA_7 = gamma_p(7)
GAMMA_11 = gamma_p(11)
GAMMA_13 = gamma_p(13)

SIN2_3 = sin2_theta_p(3)
SIN2_5 = sin2_theta_p(5)
SIN2_7 = sin2_theta_p(7)
SIN2_11 = sin2_theta_p(11)
SIN2_13 = sin2_theta_p(13)

def check(label, expr, expected=TARGET):
    rel_err = abs(expr - expected) / expected * 100
    flag = "**MATCH**" if rel_err < 1.0 else ("  ~  " if rel_err < 5 else "     ")
    print(f"  {flag}  {label:<50}  = {expr:>10.4f}  err = {rel_err:>6.2f} %")

def header(s):
    print()
    print("=" * 72)
    print(s)
    print("=" * 72)

header("TARGET: Bobeth prefactor 181.42")
print(f"  N_c = {N_C}, N_gen = {N_GEN}, mu* = {MU_STAR}, s = {S}")
print(f"  C_F = {C_F}, beta_0 = {BETA_0:.4f}, 1/alpha = {ALPHA_EM_INV:.4f}")
print(f"  gamma_3 = {GAMMA_3:.5f}, gamma_5 = {GAMMA_5:.5f}, gamma_7 = {GAMMA_7:.5f}")
print(f"  gamma_11 = {GAMMA_11:.5f}, gamma_13 = {GAMMA_13:.5f}")


# =============================================================
# Candidate family 1 : multiplicative structure mu* * k + m
# =============================================================
header("FAMILY 1 : Linear combinations mu* * k + m")

for k in range(1, 16):
    for m in range(-5, 6):
        val = MU_STAR * k + m
        if 180 <= val <= 183:
            print(f"  mu*={MU_STAR} * {k} + {m} = {val}  {'(match)' if abs(val-181.42)<0.5 else ''}")

print()
print("  Note : mu* * 12 + 1 = 181 (0.2%).  '12' could be :")
print(f"    N_c^2 * (N_spatial-1) = 9*2 = 18 (no)")
print(f"    N_c * N_gen + N_spatial = 9 + 3 = 12 (YES)")
print(f"    N_Weyl per gen + N_c - 3 = 15 + 0 = ? (no)")
print(f"    Adjoint SU(4) dim = 15 (no)")
print(f"    Dim(SU(N_c+1)) - Dim(SU(N_c)) = 15 - 8 = 7 (no)")
print(f"    2 * (2 N_c + 0) = 12 (YES with trivial factor 2)")
print(f"    Dim of 4-face T^4 = 4 * 3 = 12 (YES)")
print(f"    12 = 4 * N_c (YES, 4 = T^3 dimension + parity)")
print()
print(f"  If 12 = 4 * N_c :  181.42 =? 15 * 4 * 3 + 1 = 181 (0.23 %)")
print(f"  If additional +0.42 from NNLO RGE : 15*4*3 + 1 + 0.42 = 181.42 (exact)")


# =============================================================
# Candidate family 2 : factorial / combinatorial structures
# =============================================================
header("FAMILY 2 : Combinatorial")

check("2 * 90 + 1 (hand match)", 181, 181.42)
check("5! + 60 + 1 = 120+60+1", 181)
check("factorial(5) + factorial(4) + factorial(3) + 7", 120+24+24+13)
check("mu* + 10! / 10^5 (nonsense)", 15 + 3.6288)

# Catalan, fibonacci structures
check("F_12 + F_10 = 144 + 55 = 199", 199)
check("F_12 + F_7 = 144 + 13 = 157", 157)
check("F_11 + F_10 = 89+55", 144)

# Binomials (N_Weyl choose 2)
check("C(15, 2) + C(15, 3)", 105 + 455)
check("C(15, 2) = 105", 105)
check("C(15, 2) + C(N_c+1, 3) * (N_c+1) = 105 + 4*4", 105 + 16)

# Prime products
check("prod(p=2..7) = 2*3*5*7 = 210", 210)
check("sum primorials = 2+6+30+210", 2+6+30+210)


# =============================================================
# Candidate family 3 : RGE / alpha structure
# =============================================================
header("FAMILY 3 : RGE / alpha structures")

check("C_F * beta_0 * (1/alpha^0.5)", C_F * BETA_0 * sqrt(ALPHA_EM_INV))
check("C_F * beta_0 * mu*", C_F * BETA_0 * MU_STAR)
check("beta_0^3 / beta_1 (ad hoc)", BETA_0**3 / BETA_1)
check("beta_0 * (eta=2)^3 * 9", BETA_0 * 8 * 9)
check("(4pi)^2 * beta_0 * s", 16*pi**2 * BETA_0 * S)
check("16 * pi^2 * s = 78.96", 16*pi**2 * S)
check("32 pi^2 * s = 157.9", 32*pi**2 * S)
check("20 * (4pi)^2 / (N_c^2+2) = 20*157.9/11", 20 * 16*pi**2 / 11)
check("60 * pi = 188.5", 60*pi)
check("N_c * factorial(5) / s = 3*120/0.5", N_C*FACTORIAL_5/S)


# =============================================================
# Candidate family 4 : ADM entries (anomalous dimension matrix)
# =============================================================
header("FAMILY 4 : ADM related")

# The ADM at 3-loop has a structure involving N_c, C_F, and zeta(3)
# The Bobeth prefactor 181.42 should emerge from NNLO ADM elements
from math import pi
zeta2 = pi**2 / 6
zeta3 = 1.2020569031

check("C_F * N_c * zeta(3) * 100", C_F * N_C * zeta3 * 100)
check("C_F^2 * N_c * 100 = (16/9)*3*100", C_F**2 * N_C * 100)
check("(N_c^2 - 1) * C_F * 15 = 8 * 4/3 * 15", (N_C**2-1)*C_F*MU_STAR)
check("16/9 * 100 = 177.8", 16.0/9 * 100)
check("C_F * beta_0^2 = (4/3) * (23/3)^2", C_F * BETA_0**2)
check("beta_0^2 * (pi^2)/6 = 58.8*1.645", BETA_0**2 * zeta2)
check("3 * beta_0^2 = 3 * 58.78", 3 * BETA_0**2)
check("2 * beta_0^2 * pi/2 = 58.78*pi", 2 * BETA_0**2 * pi/2)


# =============================================================
# Candidate family 5 : combinations with gamma_p sums
# =============================================================
header("FAMILY 5 : gamma_p products")

SUM_GAMMA = GAMMA_3 + GAMMA_5 + GAMMA_7
PROD_GAMMA = GAMMA_3 * GAMMA_5 * GAMMA_7

check("100 * SUM_GAMMA = 100*2.10", 100 * SUM_GAMMA)
check("60 * sum_gamma * pi/2", 60 * SUM_GAMMA * pi/2)
check("mu* * sum_gamma * 6", MU_STAR * SUM_GAMMA * 6)
check("mu* * beta_0 / s = 15 * 23/3 / (1/2)", MU_STAR * BETA_0 / S)
check("mu* * C_F * 9 = 15*(4/3)*9", MU_STAR * C_F * 9)
check("mu* * (mu*-1) * 6/7 = 210*6/7 = 180", MU_STAR * 14 * 6 / 7)
check("(mu*)^2 - 4*mu* + 16 = 225-60+16", MU_STAR**2 - 4*MU_STAR + 16)
check("4 * mu*^2 / 5 + 1 = 4*225/5+1", 4*MU_STAR**2/5 + 1)
check("13 * 14 / (..)", 13*14 - 1)


# =============================================================
# Candidate family 6 : 181.42 = N * (1 + correction)
# =============================================================
header("FAMILY 6 : 181 * (1 + correction)")

# Try 181.42 = mu* * 4 * N_c * (1 + x)
base = MU_STAR * 4 * N_C
corr = TARGET / base - 1
print(f"  Base = mu* * 4 * N_c = {base}")
print(f"  Required correction : +{corr*100:.3f} %")
print(f"  Candidates for the correction :")

# Typical small PT numbers
print(f"    1 / (9*mu*)  = {100/(9*MU_STAR):.3f} % (rough match!)")
print(f"    1 / (5!)     = {1/120:.5f}  = {1/120*100:.3f} %")
print(f"    1/N_c^3      = {1/27*100:.3f} %")
print(f"    alpha_EM / gamma_3 = {1/ALPHA_EM_INV/GAMMA_3*100:.3f} %")
print(f"    (beta_0 * alpha_s / pi)^2 = {(BETA_0*ALPHA_S_MB/pi)**2*100:.3f} %")
print(f"    zeta_3 / (4pi) = {zeta3/(4*pi)*100:.3f} %")
print(f"    3/pi^2 - gamma_Euler = {3/pi**2 - 0.577:.3f}")

best_val = 180 * (1 + 1/120 + 1/27)
print(f"  Combined : 180 * (1 + 1/120 + 1/N_c^3) = {best_val:.4f}")

# NNLO correction proper : U_{s,92}^{(2)} is NNLO so ~alpha_s^2 correction
# Expected overall form : (some rational) * (1 + alpha_s * Delta_NLO)
# The "42" in 181.42 might be NLO residual

check("180 * (1 + (beta_0 * alpha_s / (4pi))**2 * 10)",
      180 * (1 + (BETA_0*ALPHA_S_MB/(4*pi))**2 * 10))


# =============================================================
# Candidate family 7 : (4 * N_c)^3 / something
# =============================================================
header("FAMILY 7 : cubic / quartic")

check("12^2 + 37.42", 144 + 37.42)
check("12^2 + 37 = 181", 144 + 37)
check("12^2 + factorial(5)/3.24", 144 + 120/3.24)
check("13^2 + 12 + 0.42", 169 + 12 + 0.42)
check("13^2 + 13 - 0.58 = 181.42", 169 + 13 - 0.58)
check("2 * 90.71", 181.42)
check("7^3 / 2 + 8 = 343/2 + 8 = 179.5", 7**3/2 + 8)
check("11^2 + 60 + 1 = 121+61", 121+60+1)
check("13 * 14 - 1 = 181", 13*14-1)
check("13 * 14 + 0.42 = 182.42", 13*14+0.42)

print()
print("  Most promising :  13 * 14 = 182 (off by 0.58, which is 1/(alpha*eta)?)")
print(f"  13*14 = {13*14}")

# Check 13 * 14 - 1 + 0.42 = 181.42 exact?
print(f"  13 * 14 - 1 + 0.42 = {13*14 - 1 + 0.42}")
print(f"  This means 181.42 = 182 - 0.58 = 13*14 - (1 - 0.42) = 181 + 0.42")
print(f"  Where 13 = first inactive after 11 (role II), 14 = 2*7 = 2 * p_gen3")
print(f"  Or 13*14 = (mu*-2)*(2*(mu*-8)) = 13*14")
print()


# =============================================================
# Candidate family 8 : PT algebra, primorial, holonomy integral
# =============================================================
header("FAMILY 8 : primorial-based")

check("primorial(3) * 6 + 1 = 30*6+1", PRIMORIAL_3 * 6 + 1)
check("primorial(3) * pi * 2 = 188.5", PRIMORIAL_3 * pi * 2)
check("(mu*+2) * 5! / (5-1)!  = 17 * 120/24 = 85", 17 * 5)
check("(mu*+2) * (mu*-2) + 14 = 17*13+14", 17*13 + 14)
check("N_Weyl * 12 + 1 = 15*12+1 = 181", 15*12+1)

print()
print(f"  CANDIDATE : mu* * 12 + 1 = {MU_STAR*12+1}")
print(f"  Error vs 181.42 : {(MU_STAR*12+1 - 181.42)/181.42*100:+.3f} % (= -0.232%)")
print(f"  Where '12' = 4 * N_c (spatial dimension 4 x color 3) = face volume of T^4")
print(f"  Thus 181 = mu* * 4*N_c + 1 = (N_Weyl)(T^4_face) + 1")


# =============================================================
# Candidate family 9 : NLO alpha_s correction
# =============================================================
header("FAMILY 9 : the 0.42 residual as NLO alpha_s correction")

# 181.42 = 181 + 0.42
# 181 = mu* * 4 * N_c + 1  (combinatorial)
# 0.42 = ? alpha_s correction?

# At mu_b = m_b, alpha_s ~ 0.21
# typical NLO coefficient ~ beta_0/(4pi) ~ 0.61
# 0.42 / 181 ~ 0.23 % : this matches alpha_s/pi * C_F = 0.089
# but (alpha_s/pi)^2 ~ 0.004, too small

check("mu* * 4 * N_c * (1 + alpha_s/(4*pi) * 5)",
      MU_STAR * 4 * N_C * (1 + ALPHA_S_MB/(4*pi) * 5))
check("181 + beta_0 * alpha_s / 4 = 181 + 5.75/4",
      181 + BETA_0 * ALPHA_S_MB / 4)
check("181 + 5 * C_F * alpha_s", 181 + 5 * C_F * ALPHA_S_MB)
check("181 + 23 * alpha_s / pi = 181 + 1.55",
      181 + BETA_0 * ALPHA_S_MB / pi)
check("181 + (N_c^2 - 1) * alpha_s / pi",
      181 + 8 * ALPHA_S_MB / pi)
check("181 + 2 * pi * alpha_s = 181+1.33", 181 + 2*pi*ALPHA_S_MB)

# 0.42 could be simply = 2 * pi * alpha_s ~ 1.33 (no)
# 0.42 = alpha_s * 2 ~ 0.42 (YES!)
check("181 + 2 * alpha_s = 181 + 0.425", 181 + 2 * ALPHA_S_MB)

print()
print("  CANDIDATE : 181 + 2 * alpha_s(m_b) = 181 + 0.425 = 181.42 !")
print(f"  181 + 2*0.2124 = {181 + 2*0.2124:.4f}")
print(f"  Error vs 181.42 : {(181 + 2*0.2124 - 181.42)/181.42*100:+.5f} %")
print()
print("  PROPOSED DECOMPOSITION :")
print("    181.42 = mu* * 4 * N_c + 1 + 2 * alpha_s(m_b)")
print("           = (N_Weyl) * (T^4 face) + (identity) + NLO correction")


# =============================================================
# Summary
# =============================================================
header("SUMMARY : best PT decomposition of 181.42")

decomp = MU_STAR * 4 * N_C + 1 + 2 * ALPHA_S_MB
err = (decomp - 181.42) / 181.42 * 100
print(f"  181.42 =? mu* * 4 * N_c + 1 + 2 * alpha_s(m_b)")
print(f"        =? 15 * 12 + 1 + 0.425")
print(f"        = {decomp:.4f}")
print(f"  Error = {err:+.4f} %")
print()
print(f"  Structural reading :")
print(f"    - mu* = 15 : PT fixed point (Weyl count per generation)")
print(f"    - 12 = 4 N_c : T^4 face (3 spatial + 1 binary dim times N_c color)")
print(f"    - +1 : identity coupling (tree-level normalization)")
print(f"    - 2 * alpha_s(m_b) : NLO running correction at b-quark scale")
print()
print(f"  ALTERNATIVE : 181.42 might NOT factor cleanly, being an artefact of")
print(f"  the 3-loop ADM with specific Inami-Lim function evaluations.")
print(f"  Cleaner would be a direct derivation via the U_{{s,92}}^{{(2)}} kernel.")

# Also check simplest :
print()
print("  REDUCED CANDIDATE : 181.42 ~ mu* * 12 + 1 (0.232 %)")
print("  The 0.42 residual may be cumulative RG effect rather than structural")
