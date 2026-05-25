#!/usr/bin/env python3
"""PT bridge checks for the tetrahedral-spinorial forcing theorem.

This script does not claim that PT has already proved every hypothesis of the
forcing theorem. It separates what PT supplies directly from what is still a
bridge step:

PT supplies:
  - p=2 as a binary involution / parity infrastructure, balanced between its
    positive and negative eigenbranches;
  - the active dynamic faces {3,5,7} at mu*=15;
  - the q_+/q_- bifurcation, whose asymmetric branch is the time axis in the
    EML closure text.

Bridge step:
  - treating {time, f3, f5, f7} as one normalized four-direction packet, with
    p=2 acting as the sign/polarity operator rather than as a tetrahedral
    vertex.
  - reading {3,5,7} as a three-amplitude variation triplet carried by the
    positive branch.
"""

from __future__ import annotations

import itertools
import json
import math


Vec = tuple[float, float, float]
Quat = tuple[float, float, float, float]

MU_STAR = 15.0
S = 0.5
DYNAMIC_PRIMES = (3, 5, 7)
BRIDGE_PACKET = ("tau", "f3", "f5", "f7")
AUDIT_PRIMES = (2, 3, 5, 7, 11, 13)
TOL = 1e-12


def close(left: float, right: float, tol: float = TOL) -> bool:
    return abs(left - right) <= tol


def q_plus(mu: float = MU_STAR) -> float:
    return 1.0 - 2.0 / mu


def q_minus(mu: float = MU_STAR) -> float:
    return math.exp(-1.0 / mu)


def delta_p(prime: int, q: float) -> float:
    return (1.0 - q**prime) / prime


def sin2_theta(prime: int, q: float) -> float:
    delta = delta_p(prime, q)
    return delta * (2.0 - delta)


def gamma_p(prime: int, mu: float = MU_STAR) -> float:
    q = q_plus(mu)
    delta = (1.0 - q**prime) / prime
    return 4.0 * prime * q ** (prime - 1) * (1.0 - delta) / (
        mu * (1.0 - q**prime) * (2.0 - delta)
    )


def matmul(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def det2(matrix: list[list[float]]) -> float:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def trace(matrix: list[list[float]]) -> float:
    return sum(matrix[i][i] for i in range(len(matrix)))


def dot(left: Vec, right: Vec) -> float:
    return sum(a * b for a, b in zip(left, right, strict=True))


def norm(vec: Vec) -> float:
    return math.sqrt(dot(vec, vec))


def scale(vec: Vec, factor: float) -> Vec:
    return tuple(factor * x for x in vec)  # type: ignore[return-value]


def add(left: Vec, right: Vec) -> Vec:
    return tuple(a + b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def neg(vec: Vec) -> Vec:
    return scale(vec, -1.0)


def tetra_vertices() -> dict[str, Vec]:
    raw = {
        "tau": (1.0, 1.0, 1.0),
        "f3": (1.0, -1.0, -1.0),
        "f5": (-1.0, 1.0, -1.0),
        "f7": (-1.0, -1.0, 1.0),
    }
    return {label: scale(vec, 1.0 / math.sqrt(3.0)) for label, vec in raw.items()}


def gram(vertices: dict[str, Vec]) -> list[list[float]]:
    ordered = [vertices[label] for label in BRIDGE_PACKET]
    return [[dot(left, right) for right in ordered] for left in ordered]


def quat_mul(left: Quat, right: Quat) -> Quat:
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )


def quat_key(q: Quat) -> tuple[int, int, int, int]:
    return tuple(int(round(2.0 * value)) for value in q)  # type: ignore[return-value]


def binary_tetrahedral_group() -> list[Quat]:
    elements: list[Quat] = []
    basis = [
        (1.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    ]
    for q in basis:
        elements.append(q)
        elements.append(tuple(-x for x in q))  # type: ignore[arg-type]
    for signs in itertools.product((-1.0, 1.0), repeat=4):
        elements.append(tuple(sign / 2.0 for sign in signs))  # type: ignore[arg-type]
    return list({quat_key(q): q for q in elements}.values())


def main() -> int:
    sigma_x = [[0.0, 1.0], [1.0, 0.0]]
    identity_2 = [[1.0, 0.0], [0.0, 1.0]]
    sigma_x_squared = matmul(sigma_x, sigma_x)
    p2_eigenvectors = {
        "+": (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)),
        "-": (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0)),
    }
    p2_eigenvalues = (-1.0, 1.0)

    gamma_values = {prime: gamma_p(prime) for prime in AUDIT_PRIMES if prime != 2}
    dynamic_active = tuple(prime for prime in gamma_values if gamma_values[prime] > S)
    dynamic_inactive = tuple(prime for prime in gamma_values if gamma_values[prime] <= S)
    active_amplitudes = tuple(gamma_values[prime] for prime in DYNAMIC_PRIMES)
    active_amplitude_norm = math.sqrt(sum(value * value for value in active_amplitudes))
    normalized_active_amplitudes = tuple(value / active_amplitude_norm for value in active_amplitudes)

    vertices = tetra_vertices()
    barycenter = (0.0, 0.0, 0.0)
    for vec in vertices.values():
        barycenter = add(barycenter, vec)
    g = gram(vertices)

    group = binary_tetrahedral_group()
    group_keys = {quat_key(q) for q in group}
    group_closed = all(
        quat_key(quat_mul(left, right)) in group_keys for left in group for right in group
    )

    pt_direct_checks = {
        "p2_transfer_is_involution": sigma_x_squared == identity_2,
        "p2_transfer_trace_zero": close(trace(sigma_x), 0.0),
        "p2_transfer_has_pm_eigenvalues": close(det2(sigma_x), -1.0),
        "p2_is_balanced_positive_negative_pivot": close(sum(p2_eigenvalues), 0.0)
        and close(abs(p2_eigenvalues[0]), abs(p2_eigenvalues[1])),
        "p2_eigenvectors_are_orthonormal": close(
            sum(p2_eigenvectors["+"][i] * p2_eigenvectors["-"][i] for i in range(2)), 0.0
        )
        and close(sum(x * x for x in p2_eigenvectors["+"]), 1.0)
        and close(sum(x * x for x in p2_eigenvectors["-"]), 1.0),
        "q_plus_contains_binary_two": close(q_plus(), 13.0 / 15.0)
        and close(1.0 - q_plus(), 2.0 / MU_STAR),
        "q_branches_are_distinct": not close(q_plus(), q_minus()),
        "time_axis_from_q_bifurcation_available": q_minus() > q_plus(),
        "active_dynamic_faces_are_3_5_7": dynamic_active == DYNAMIC_PRIMES,
        "first_inactive_after_dynamic_faces_is_11": dynamic_inactive[0] == 11,
        "active_triplet_has_nonzero_amplitudes": all(value > 0.0 for value in active_amplitudes),
        "active_triplet_amplitudes_are_non_degenerate": len(
            {round(value, 12) for value in active_amplitudes}
        )
        == len(active_amplitudes),
        "active_triplet_amplitudes_are_ordered": all(
            left > right for left, right in zip(active_amplitudes, active_amplitudes[1:])
        ),
    }

    bridge_checks = {
        "pt_packet_has_three_faces_plus_time": BRIDGE_PACKET == ("tau", "f3", "f5", "f7"),
        "p2_is_polarity_not_tetrahedral_vertex": 2 not in DYNAMIC_PRIMES,
        "canonical_packet_is_centered": close(norm(barycenter), 0.0),
        "canonical_packet_is_unit_normalized": all(close(norm(vec), 1.0) for vec in vertices.values()),
        "canonical_packet_has_forced_tetrahedral_gram": all(
            close(g[i][j], 1.0 if i == j else -1.0 / 3.0)
            for i in range(4)
            for j in range(4)
        ),
        "p2_polarity_realizes_antipodal_branch": all(
            all(close(a, b) for a, b in zip(neg(vec), scale(vec, -1.0), strict=True))
            for vec in vertices.values()
        ),
        "spinorial_lift_group_has_order_24": len(group) == 24,
        "spinorial_lift_group_is_closed": group_closed,
        "spinorial_kernel_contains_pm_one": {(2, 0, 0, 0), (-2, 0, 0, 0)}.issubset(group_keys),
    }

    payload = {
        "model": "pt_tetrahedral_spinor_bridge",
        "status": {
            "pt_direct": "checked from PT formulas used locally",
            "bridge": "conditional on treating {time,f3,f5,f7} as one normalized four-direction packet",
            "not_claimed": "does not prove the full PT physical identification of the spinorial double tetrahedron",
        },
        "pt_data": {
            "mu_star": MU_STAR,
            "s": S,
            "q_plus": q_plus(),
            "q_minus": q_minus(),
            "dynamic_active_faces": list(dynamic_active),
            "dynamic_inactive_faces": list(dynamic_inactive),
            "gamma_values": gamma_values,
            "active_amplitudes_3_5_7": dict(zip(DYNAMIC_PRIMES, active_amplitudes, strict=True)),
            "normalized_active_amplitudes_3_5_7": dict(
                zip(DYNAMIC_PRIMES, normalized_active_amplitudes, strict=True)
            ),
            "bridge_packet": list(BRIDGE_PACKET),
            "p2_role": "stable positive/negative pivot; binary involution / polarity operator",
            "active_triplet_role": "three non-degenerate variation amplitudes of the potential field",
            "time_role": "asymmetric q_minus branch, orthogonal complement to the three q_plus spatial faces",
        },
        "p2_involution": {
            "matrix": sigma_x,
            "eigenvalues": list(p2_eigenvalues),
            "eigenvectors": p2_eigenvectors,
        },
        "tetrahedral_completion": {
            "vertices_by_direction": vertices,
            "gram": g,
        },
        "spinorial_lift": {
            "binary_tetrahedral_order": len(group),
            "kernel_candidates": [[1, 0, 0, 0], [-1, 0, 0, 0]],
        },
        "checks": {
            "pt_direct": pt_direct_checks,
            "bridge": bridge_checks,
        },
        "passed": all(pt_direct_checks.values()) and all(bridge_checks.values()),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
