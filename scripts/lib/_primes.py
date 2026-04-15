#!/usr/bin/env python3
"""
_primes.py — Prime number generation with fast C fallback.

Uses primesieve (C library) when available, falls back to pure-Python
sieve of Eratosthenes. Cross-platform, no external dependencies required.
"""

import math


def generate_primes(n):
    """Return a list of the first n primes.

    Uses primesieve (C library) if installed for speed on large n,
    otherwise falls back to a pure-Python Eratosthenes sieve.

    Parameters
    ----------
    n : int
        Number of primes to generate.

    Returns
    -------
    list[int]
        The first n prime numbers [2, 3, 5, 7, 11, ...].
    """
    # Fast path: use primesieve C library if available
    try:
        import primesieve
        it = primesieve.Iterator()
        return [it.next_prime() for _ in range(n)]
    except ImportError:
        pass

    # Pure Python fallback: sieve of Eratosthenes
    if n <= 0:
        return []
    if n < 6:
        limit = 15
    else:
        # Upper bound: p_n ~ n * (ln(n) + ln(ln(n))) for n >= 6
        limit = int(n * (math.log(n) + math.log(math.log(n))) * 1.3) + 200

    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    primes = [i for i in range(2, limit + 1) if sieve[i]]

    # If estimate was too low, extend
    if len(primes) < n:
        limit2 = int(limit * 1.5) + 1000
        sieve2 = [True] * (limit2 + 1)
        sieve2[0] = sieve2[1] = False
        for i in range(2, int(limit2 ** 0.5) + 1):
            if sieve2[i]:
                for j in range(i * i, limit2 + 1, i):
                    sieve2[j] = False
        primes = [i for i in range(2, limit2 + 1) if sieve2[i]]

    return primes[:n]
