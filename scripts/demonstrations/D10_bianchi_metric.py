#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Theorem D10: Bianchi I metric with scale factors a_p = gamma_p / mu.
GENUINE TEST: Compute gamma_p from first principles via numerical differentiation.
"""
import numpy as np

def sin2(p, mu):
    q = 1.0 - 2.0/mu
    delta = (1.0 - q**p) / p
    return delta * (2.0 - delta)

def gamma_p(p, mu, dmu=0.001):
    """Anomalous dimension: gamma_p = -d(ln sin^2)/d(ln mu)"""
    s_plus = sin2(p, mu + dmu)
    s_minus = sin2(p, mu - dmu)
    dlns = (np.log(s_plus) - np.log(s_minus)) / (2*dmu)
    dlnmu = 1.0 / mu
    return -dlns / (1.0/mu)  # = -mu * d(ln sin^2)/d(mu)

mu = 15.0
active = [3, 5, 7]
print('Bianchi I scale factors at mu* = {}:'.format(mu))
print('{:>5s} {:>10s} {:>10s}'.format('p', 'gamma_p', 'a_p'))

gammas = {}
for p in active:
    g = gamma_p(p, mu)
    a = g / mu
    gammas[p] = g
    print('{:5d} {:10.4f} {:10.6f}'.format(p, g, a))

# Verify hierarchy: gamma_3 > gamma_5 > gamma_7 > 1/2
assert gammas[3] > gammas[5] > gammas[7], 'FAIL: hierarchy broken'
assert all(gammas[p] > 0.5 for p in active), 'FAIL: some gamma < 1/2'
# Verify gamma_11 < 1/2
g11 = gamma_p(11, mu)
print('\ngamma_11 = {:.4f} (should be < 1/2)'.format(g11))
assert g11 < 0.5, 'FAIL: gamma_11 = {:.4f} >= 1/2'.format(g11)

print('\nD10 VERIFIED: Bianchi I metric computed, hierarchy gamma_3>gamma_5>gamma_7>1/2>gamma_11.')
