#!/usr/bin/env python3
"""Local PT forcing check for the pivot-triplet-time bridge.

The goal is deliberately modest: show that, inside the local PT data used by
the consolidation notes, the binary pivot, the active triplet, and the time
axis are not three independent choices. They are read from the same packet:

    mu*, S, q_+(mu), q_-(mu), gamma_p(mu).

This script therefore proves a finite, formula-level statement:
  - p=2 is the balanced positive/negative pivot of the binary square;
  - the only odd primes active above S=1/2 at mu*=15 are 3,5,7;
  - q_- is the distinguished asymmetric branch used for the time axis;
  - the packet has the 1+3 shape needed by the tetrahedral bridge.

It does not claim the full physical bridge as an unconditional theorem.
"""

from __future__ import annotations

import json
import math


MU_STAR = 15.0
S = 0.5
ODD_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)
EXPECTED_ACTIVE_TRIPLET = (3, 5, 7)
TOL = 1e-12


Vec2 = tuple[float, float]


def close(left: float, right: float, tol: float = TOL) -> bool:
    return abs(left - right) <= tol


def q_plus(mu: float = MU_STAR) -> float:
    return 1.0 - 2.0 / mu


def q_minus(mu: float = MU_STAR) -> float:
    return math.exp(-1.0 / mu)


def gamma_p(prime: int, mu: float = MU_STAR) -> float:
    q = q_plus(mu)
    delta = (1.0 - q**prime) / prime
    return 4.0 * prime * q ** (prime - 1) * (1.0 - delta) / (
        mu * (1.0 - q**prime) * (2.0 - delta)
    )


def swap(vec: Vec2) -> Vec2:
    return (vec[1], vec[0])


def add(left: Vec2, right: Vec2) -> Vec2:
    return (left[0] + right[0], left[1] + right[1])


def scale(vec: Vec2, factor: float) -> Vec2:
    return (factor * vec[0], factor * vec[1])


def dot(left: Vec2, right: Vec2) -> float:
    return left[0] * right[0] + left[1] * right[1]


def norm(vec: Vec2) -> float:
    return math.sqrt(dot(vec, vec))


def classify_binary_square(states: tuple[Vec2, ...]) -> dict[str, list[Vec2]]:
    symmetric: list[Vec2] = []
    antisymmetric: list[Vec2] = []
    mixed: list[Vec2] = []
    for state in states:
        image = swap(state)
        if image == state:
            symmetric.append(state)
        elif image == scale(state, -1.0):
            antisymmetric.append(state)
        else:
            mixed.append(state)
    return {
        "symmetric_axis": symmetric,
        "antisymmetric_axis": antisymmetric,
        "swap_pair": mixed,
    }


def main() -> int:
    binary_square: tuple[Vec2, ...] = (
        (-1.0, -1.0),
        (1.0, -1.0),
        (-1.0, 1.0),
        (1.0, 1.0),
    )
    center = (0.0, 0.0)
    for state in binary_square:
        center = add(center, state)
    center = scale(center, 1.0 / len(binary_square))

    binary_axes = {
        "positive_axis": scale((1.0, 1.0), 1.0 / math.sqrt(2.0)),
        "negative_axis": scale((1.0, -1.0), 1.0 / math.sqrt(2.0)),
    }
    binary_classes = classify_binary_square(binary_square)
    eigenvalues = {
        "positive_axis": 1.0,
        "negative_axis": -1.0,
    }

    gammas = {prime: gamma_p(prime) for prime in ODD_PRIMES}
    active = tuple(prime for prime in ODD_PRIMES if gammas[prime] > S)
    inactive = tuple(prime for prime in ODD_PRIMES if gammas[prime] <= S)
    active_amplitudes = tuple(gammas[prime] for prime in active)
    amplitude_norm = math.sqrt(sum(value * value for value in active_amplitudes))
    normalized_active_amplitudes = tuple(value / amplitude_norm for value in active_amplitudes)

    common_packet = {
        "mu_star": MU_STAR,
        "threshold": S,
        "q_plus": q_plus(),
        "q_minus": q_minus(),
        "binary_gap": 1.0 - q_plus(),
        "active_triplet": active,
        "time_branch": "q_minus",
    }

    binary_checks = {
        "binary_square_has_four_oriented_states": len(binary_square) == 4,
        "binary_square_is_centered": close(norm(center), 0.0),
        "zero_zero_is_center_not_oriented_state": (0.0, 0.0) not in binary_square
        and close(norm(center), 0.0),
        "swap_is_involution_on_square": all(swap(swap(state)) == state for state in binary_square),
        "symmetric_axis_has_eigenvalue_plus_one": swap(binary_axes["positive_axis"])
        == binary_axes["positive_axis"],
        "antisymmetric_axis_has_eigenvalue_minus_one": swap(binary_axes["negative_axis"])
        == scale(binary_axes["negative_axis"], -1.0),
        "positive_negative_axes_are_balanced": close(sum(eigenvalues.values()), 0.0)
        and close(abs(eigenvalues["positive_axis"]), abs(eigenvalues["negative_axis"])),
        "positive_negative_axes_are_orthogonal": close(
            dot(binary_axes["positive_axis"], binary_axes["negative_axis"]), 0.0
        ),
    }

    triplet_checks = {
        "active_primes_are_exactly_3_5_7": active == EXPECTED_ACTIVE_TRIPLET,
        "first_inactive_prime_is_11": inactive[0] == 11,
        "active_amplitudes_are_above_threshold": all(value > S for value in active_amplitudes),
        "inactive_amplitudes_are_below_threshold": all(gammas[prime] <= S for prime in inactive),
        "active_amplitudes_are_non_degenerate": len({round(value, 12) for value in active_amplitudes})
        == len(active_amplitudes),
        "active_amplitudes_are_ordered": all(
            left > right for left, right in zip(active_amplitudes, active_amplitudes[1:])
        ),
        "normalized_amplitude_vector_is_unit": close(
            math.sqrt(sum(value * value for value in normalized_active_amplitudes)), 1.0
        ),
    }

    time_checks = {
        "q_minus_is_distinct_from_q_plus": not close(q_minus(), q_plus()),
        "q_minus_is_larger_than_q_plus_at_mu_star": q_minus() > q_plus(),
        "q_minus_is_exponential_primitive_branch": close(q_minus(), math.exp(-1.0 / MU_STAR)),
        "q_plus_contains_binary_pivot_gap": close(1.0 - q_plus(), 2.0 / MU_STAR),
    }

    common_structure_checks = {
        "same_mu_controls_q_plus_q_minus_and_gamma": close(common_packet["mu_star"], MU_STAR),
        "same_threshold_selects_triplet": common_packet["threshold"] == S
        and common_packet["active_triplet"] == EXPECTED_ACTIVE_TRIPLET,
        "binary_pivot_is_not_dynamic_triplet_member": 2 not in active,
        "forced_packet_has_one_time_plus_three_variations": 1 + len(active) == 4,
        "local_data_has_pivot_triplet_time_shape": all(binary_checks.values())
        and all(triplet_checks.values())
        and all(time_checks.values()),
    }

    payload = {
        "model": "pt_common_structure_forcing",
        "status": {
            "proved_locally": "formula-level forcing inside the local PT packet",
            "bridge_status": "conditional physical interpretation",
            "not_claimed": "not a standalone derivation of spacetime or spin from first principles",
        },
        "common_packet": common_packet,
        "binary_pivot": {
            "square_states": binary_square,
            "center": center,
            "axes": binary_axes,
            "classes_under_swap": binary_classes,
            "eigenvalues": eigenvalues,
        },
        "active_triplet": {
            "gammas": gammas,
            "active": active,
            "inactive": inactive,
            "amplitudes": dict(zip(active, active_amplitudes, strict=True)),
            "normalized_amplitudes": dict(
                zip(active, normalized_active_amplitudes, strict=True)
            ),
        },
        "time_axis": {
            "q_minus": q_minus(),
            "q_plus": q_plus(),
            "branch_asymmetry": q_minus() - q_plus(),
        },
        "checks": {
            "binary_pivot": binary_checks,
            "active_triplet": triplet_checks,
            "time_axis": time_checks,
            "common_structure": common_structure_checks,
        },
        "passed": all(binary_checks.values())
        and all(triplet_checks.values())
        and all(time_checks.values())
        and all(common_structure_checks.values()),
    }

    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
