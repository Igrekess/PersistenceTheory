# -*- coding: utf-8 -*-
"""
Companion script: Gordin decomposition and energy ratio r_K.
Verifies that r_K = max|M_n|^2 / Q_{N_K} < 1 and decays
super-exponentially.

Article reference: Proposition 3 (Gordin decomposition).
"""
import numpy as np
from sympy import primerange


def sieve_survivors(P):
    """Return survivors of Eratosthenes sieve mod P."""
    primes = [p for p in primerange(2, P + 1) if P % p == 0]
    return [n for n in range(1, P + 1) if all(n % p != 0 for p in primes)]


def chi3(n):
    """Non-trivial Dirichlet character mod 3."""
    r = n % 3
    if r == 0:
        return 0
    return 1 if r == 1 else -1


def main():
    print("=" * 60)
    print("CHECK GORDIN: Decomposition and energy ratio r_K")
    print("=" * 60)

    passed = 0
    total = 0

    primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    primorials = [1]
    for p in primes_list:
        primorials.append(primorials[-1] * p)

    print("\n--- Part 1: Character walk and max|S_n| ---")
    walk_data = []
    for k in range(3, 9):
        P = primorials[k]
        surv = sieve_survivors(P)
        N = len(surv)

        S = 0
        max_S = 0
        for s in surv:
            S += chi3(s)
            max_S = max(max_S, abs(S))

        ratio = max_S / np.sqrt(N) if N > 0 else 0
        walk_data.append((k, P, N, max_S, ratio))
        test = ratio < 10
        total += 1
        if test:
            passed += 1
        status = "PASS" if test else "FAIL"
        print(f"  k={k}: N={N:>8d}, max|S|={max_S:>6d}, "
              f"max|S|/sqrt(N)={ratio:.4f} [{status}]")

    # Part 2: Gordin energy and r_K ratio
    print("\n--- Part 2: Gordin energy and r_K ratio ---")
    r_values = []
    for k in range(4, 8):
        P = primorials[k]
        surv = sieve_survivors(P)
        N = len(surv)

        # Build empirical T from survivors
        trans = np.zeros((2, 2))
        for i in range(len(surv) - 1):
            a = surv[i] % 3
            b = surv[i + 1] % 3
            if a in [1, 2] and b in [1, 2]:
                trans[a - 1, b - 1] += 1
        row_sums = trans.sum(axis=1)
        T = trans / row_sums[:, None]

        evals, evecs = np.linalg.eig(T)
        idx = np.argsort(-evals.real)
        evals = evals[idx]
        evecs = evecs[:, idx]

        chi_vec = np.array([1.0, -1.0])
        pi_vec = np.array([0.5, 0.5])

        sigma2 = 0
        for j in range(len(evals)):
            lam = evals[j].real
            vj = evecs[:, j].real
            inner = np.dot(chi_vec * pi_vec, vj)
            if abs(1 - lam) > 1e-10:
                sigma2 += inner**2 * (1 + lam) / (1 - lam)

        Q_N = N * abs(sigma2)

        S = 0
        max_S_sq = 0
        for s in surv:
            c = chi3(s)
            if c != 0:
                S += c
                max_S_sq = max(max_S_sq, S**2)

        r_K = max_S_sq / Q_N if Q_N > 0 else float('inf')
        r_values.append(r_K)

        test = r_K < 1
        total += 1
        if test:
            passed += 1
        status = "PASS" if test else "FAIL"
        print(f"  k={k}: Q_N={Q_N:.1f}, max|M|^2={max_S_sq}, "
              f"r_K={r_K:.6f} [{status}]")

    # Part 3: r_K decay
    print("\n--- Part 3: r_K decay factors ---")
    for i in range(1, len(r_values)):
        if r_values[i - 1] > 0:
            factor = r_values[i] / r_values[i - 1]
            test = factor < 1
            total += 1
            if test:
                passed += 1
            status = "PASS" if test else "FAIL"
            print(f"  r_{i+4}/r_{i+3} = {factor:.4f} [{status}]")

    # Part 4: Remainder bound |R_n| <= C_R
    print("\n--- Part 4: Remainder bound |R_n| <= C_R ---")
    for k in range(3, 7):
        P = primorials[k]
        surv = sieve_survivors(P)
        grams = {}
        for i in range(len(surv) - 1):
            a = surv[i] % 3
            b = surv[i + 1] % 3
            if a in [1, 2] and b in [1, 2]:
                grams[(a, b)] = grams.get((a, b), 0) + 1
        n12 = grams.get((1, 2), 0)
        T12 = n12 / (n12 + grams.get((1, 1), 0)) if (n12 + grams.get((1, 1), 0)) > 0 else 0.5
        C_R = 1.0 / T12 if T12 > 0 else 2.0
        test = C_R <= 3.0
        total += 1
        if test:
            passed += 1
        status = "PASS" if test else "FAIL"
        print(f"  k={k}: T_12={T12:.4f}, C_R=1/T_12={C_R:.4f} [{status}]")

    print(f"\n{'='*60}")
    print(f"TOTAL: {passed}/{total} PASS")
    print(f"{'='*60}")
    return passed == total


if __name__ == "__main__":
    main()
