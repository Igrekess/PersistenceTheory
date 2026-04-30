#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_qg_topology_changing_foams_PT.py -- PT topology-changing foam sector
=========================================================================

This script tests the topology-changing sector forced by the physical Perron
projector Pi_0 = |1><pi|.

On the physical quotient, each connected boundary component carries the
rank-one Frobenius algebra

    eta(1) = |1>
    eps(v) = <pi|v>
    m(v,w) = |1> <pi|v <pi|w
    Delta(v) = |1> tensor |1> <pi|v

These four maps generate cups, caps, pair-of-pants mergers/splits, and
handles.  The tests verify the Frobenius/TQFT relations:

  - unit/counit on the physical quotient;
  - associativity, commutativity, coassociativity, cocommutativity;
  - Frobenius relation;
  - specialness m Delta = Pi_0, hence handles are identity physically;
  - genus and pair-of-pants decompositions give the same amplitude;
  - disjoint union factorizes.

This closes topology-changing foams in the projected PT physical sector.
Higher Fourier/RG decorations and black-hole/ringdown precision skeletons
are tested in the neighbouring QG scripts.
"""

import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def sieve_primes(n_max):
    is_prime = [True] * (n_max + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n_max**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n_max + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n_max + 1) if is_prime[i]]


def build_T_gaps(gaps, modulus):
    residues = [g % modulus for g in gaps]
    classes = sorted(set(residues))
    idx = {c: i for i, c in enumerate(classes)}
    counts = np.zeros((len(classes), len(classes)), dtype=float)
    for a, b in zip(residues, residues[1:]):
        counts[idx[a], idx[b]] += 1.0

    T = np.zeros_like(counts)
    for i in range(counts.shape[0]):
        row_sum = counts[i].sum()
        if row_sum > 0:
            T[i] = counts[i] / row_sum
    return T, classes


def stationary_left(T):
    eigvals, eigvecs = np.linalg.eig(T.T)
    idx = np.argmin(np.abs(eigvals - 1.0))
    pi = np.real(eigvecs[:, idx])
    if pi.sum() < 0:
        pi = -pi
    pi = np.maximum(pi, 0.0)
    return pi / pi.sum()


def mat_norm(A):
    return float(np.linalg.norm(A))


def check(name, condition, detail=""):
    global n_pass, n_total
    n_total += 1
    status = "PASS" if condition else "FAIL"
    if condition:
        n_pass += 1
    print(f"  T{n_total:02d} [{status}] {name}")
    if detail:
        print(f"       {detail}")


def frobenius_maps(T):
    n = T.shape[0]
    eta = np.ones((n, 1), dtype=complex)
    eps = stationary_left(T).reshape(1, n).astype(complex)
    Pi = eta @ eps
    # m: H tensor H -> H, Delta: H -> H tensor H
    m = np.kron(eps, eps)
    m = eta @ m
    Delta = np.kron(eta, eta) @ eps
    return eta, eps, Pi, m, Delta


def random_vec(n, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, 1)) + 1j * rng.normal(size=(n, 1))


def surface_amplitude(T, incoming, outgoing, genus=0):
    """Connected physical surface amplitude from n incoming to m outgoing."""
    eta, eps, Pi, m, Delta = frobenius_maps(T)
    in_weight = 1.0 + 0j
    for v in incoming:
        in_weight *= (eps @ v)[0, 0]
    out_state = np.array([[1.0 + 0j]])
    for _ in outgoing:
        out_state = np.kron(out_state, eta)
    # Handles act as m Delta = Pi, hence identity on the physical line.
    handle_weight = (eps @ Pi @ eta)[0, 0] ** genus
    return handle_weight * in_weight * out_state


print("=" * 82)
print("  PT QUANTUM GRAVITY: TOPOLOGY-CHANGING FOAMS")
print("=" * 82)

primes = [p for p in sieve_primes(100000) if p >= 5]
gaps = [b - a for a, b in zip(primes, primes[1:])]
T6, classes6 = build_T_gaps(gaps, 6)
T30, classes30 = build_T_gaps(gaps, 30)

print(f"\n  T_6  shape = {T6.shape}, classes = {classes6}")
print(f"  T_30 shape = {T30.shape}, classes = {classes30}")
print()

n_pass = 0
n_total = 0

for label, T in (("T6", T6), ("T30", T30)):
    n = T.shape[0]
    I = np.eye(n, dtype=complex)
    eta, eps, Pi, m, Delta = frobenius_maps(T)
    v = random_vec(n, 100 + n)
    w = random_vec(n, 200 + n)
    u = random_vec(n, 300 + n)

    print(f"--- {label}: boundary dimension {n} ---")

    norm_eta = (eps @ eta)[0, 0]
    pi_err = mat_norm(Pi @ Pi - Pi)
    check(f"{label}: normalized Perron Frobenius line",
          abs(norm_eta - 1.0) < 1e-12 and pi_err < 1e-12,
          f"<pi|1>={norm_eta:.12f}, ||Pi^2-Pi||={pi_err:.2e}")

    left_unit = m @ np.kron(eta, I)
    right_unit = m @ np.kron(I, eta)
    unit_err = mat_norm(left_unit - Pi) + mat_norm(right_unit - Pi)
    check(f"{label}: pair-of-pants unit laws hold physically",
          unit_err < 1e-12,
          f"||m(eta,id)-Pi||+||m(id,eta)-Pi||={unit_err:.2e}")

    left_counit = np.kron(eps, I) @ Delta
    right_counit = np.kron(I, eps) @ Delta
    counit_err = mat_norm(left_counit - Pi) + mat_norm(right_counit - Pi)
    check(f"{label}: split counit laws hold physically",
          counit_err < 1e-12,
          f"||(.eps,id)Delta-Pi||+||(id,eps)Delta-Pi||={counit_err:.2e}")

    assoc_left = m @ np.kron(m, I)
    assoc_right = m @ np.kron(I, m)
    assoc_err = mat_norm(assoc_left - assoc_right)
    check(f"{label}: merger associativity is topology independent",
          assoc_err < 1e-12,
          f"||m(m,id)-m(id,m)||={assoc_err:.2e}")

    coassoc_left = np.kron(Delta, I) @ Delta
    coassoc_right = np.kron(I, Delta) @ Delta
    coassoc_err = mat_norm(coassoc_left - coassoc_right)
    check(f"{label}: splitter coassociativity is topology independent",
          coassoc_err < 1e-12,
          f"||(Delta,id)Delta-(id,Delta)Delta||={coassoc_err:.2e}")

    swap = np.zeros((n * n, n * n), dtype=complex)
    for i in range(n):
        for j in range(n):
            swap[j * n + i, i * n + j] = 1.0
    comm_err = mat_norm(m @ swap - m) + mat_norm(swap @ Delta - Delta)
    check(f"{label}: boundary exchange symmetry",
          comm_err < 1e-12,
          f"||m tau-m||+||tau Delta-Delta||={comm_err:.2e}")

    frob_left = Delta @ m
    frob_mid = np.kron(m, I) @ np.kron(I, Delta)
    frob_right = np.kron(I, m) @ np.kron(Delta, I)
    frob_err = mat_norm(frob_left - frob_mid) + mat_norm(frob_left - frob_right)
    check(f"{label}: Frobenius relation holds",
          frob_err < 1e-12,
          f"Frobenius relation error={frob_err:.2e}")

    handle = m @ Delta
    handle_err = mat_norm(handle - Pi)
    check(f"{label}: handle is physically trivial",
          handle_err < 1e-12,
          f"||m Delta - Pi||={handle_err:.2e}")

    sphere = (eps @ eta)[0, 0]
    torus = (eps @ handle @ eta)[0, 0]
    genus3 = (eps @ np.linalg.matrix_power(handle, 3) @ eta)[0, 0]
    genus_err = abs(sphere - 1.0) + abs(torus - 1.0) + abs(genus3 - 1.0)
    check(f"{label}: closed connected topology is genus-invariant physically",
          genus_err < 1e-12,
          f"Z(S2)={sphere:.3f}, Z(T2)={torus:.3f}, Z(g=3)={genus3:.3f}")

    merge_12_then_3 = m @ np.kron(m @ np.kron(v, w), u)
    merge_1_then_23 = m @ np.kron(v, m @ np.kron(w, u))
    pair_err = mat_norm(merge_12_then_3 - merge_1_then_23)
    check(f"{label}: two pair-of-pants decompositions agree on states",
          pair_err < 1e-10,
          f"state-level pair-of-pants error={pair_err:.2e}")

    split_12_then_3 = np.kron(Delta, I) @ Delta @ v
    split_1_then_23 = np.kron(I, Delta) @ Delta @ v
    split_err = mat_norm(split_12_then_3 - split_1_then_23)
    check(f"{label}: two splitting decompositions agree on states",
          split_err < 1e-10,
          f"state-level split error={split_err:.2e}")

    amp_connected = surface_amplitude(T, [v, w], [0, 1, 2], genus=2)
    amp_decomp = surface_amplitude(T, [m @ np.kron(v, w)], [0, 1, 2], genus=2)
    decomp_err = mat_norm(amp_connected - amp_decomp)
    check(f"{label}: general connected foam amplitude is decomposition independent",
          decomp_err < 1e-10,
          f"connected amplitude decomposition error={decomp_err:.2e}")

    amp_a = surface_amplitude(T, [v], [0], genus=0)
    amp_b = surface_amplitude(T, [w], [0, 1], genus=1)
    disjoint = np.kron(amp_a, amp_b)
    combined = surface_amplitude(T, [v, w], [0, 1, 2], genus=1)
    # The rank-one physical theory factorizes, so disconnected union equals
    # the connected physical amplitude with the same projected boundary weights.
    disjoint_err = mat_norm(disjoint - combined)
    check(f"{label}: disjoint union factorizes in the physical sector",
          disjoint_err < 1e-10,
          f"disjoint-union factorization error={disjoint_err:.2e}")

    # Negative control: a random multiplication does not satisfy Frobenius.
    rng = np.random.default_rng(900 + n)
    m_raw = rng.normal(size=(n, n * n)) + 1j * rng.normal(size=(n, n * n))
    Delta_raw = rng.normal(size=(n * n, n)) + 1j * rng.normal(size=(n * n, n))
    raw_err = mat_norm(Delta_raw @ m_raw - np.kron(m_raw, I) @ np.kron(I, Delta_raw))
    check(f"{label}: topology invariance is not automatic",
          raw_err > 1e-3,
          f"negative-control Frobenius error={raw_err:.2e}")

print("\n" + "=" * 82)
print(f"  SCORE: {n_pass}/{n_total} PASS")
print("  CLOSED HERE: topology-changing foams in the projected physical sector,")
print("               via the rank-one Frobenius/TQFT structure forced by Pi0.")
print("  RELATED: Fourier/RG decorations and black-hole/ringdown precision skeletons")
print("           are closed in the neighbouring QG scripts.")
print("=" * 82)

sys.exit(0 if n_pass == n_total else 1)
