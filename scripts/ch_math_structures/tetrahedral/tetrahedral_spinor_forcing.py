#!/usr/bin/env python3
"""Finite checks for the tetrahedral and spinorial forcing theorem.

The theorem has three layers:

1. A normalized four-state packet with zero barycenter and equal pairwise
   correlations is forced to be a regular tetrahedron.
2. A free antipodal branch involution forces the signed double tetrahedron.
3. The orientation-preserving tetrahedral rotations lift to the binary
   tetrahedral group in unit quaternions, giving the spinorial double cover.
"""

from __future__ import annotations

import itertools
import json
import math


Vec = tuple[float, float, float]
Quat = tuple[float, float, float, float]
Perm = tuple[int, int, int, int]

TOL = 1e-12


def close(left: float, right: float, tol: float = TOL) -> bool:
    return abs(left - right) <= tol


def dot(left: Vec, right: Vec) -> float:
    return sum(a * b for a, b in zip(left, right, strict=True))


def norm(vec: Vec) -> float:
    return math.sqrt(dot(vec, vec))


def scale(vec: Vec, factor: float) -> Vec:
    return tuple(factor * x for x in vec)  # type: ignore[return-value]


def add(left: Vec, right: Vec) -> Vec:
    return tuple(a + b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def sub(left: Vec, right: Vec) -> Vec:
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def neg(vec: Vec) -> Vec:
    return scale(vec, -1.0)


def tetra_vertices() -> list[Vec]:
    raw = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    return [scale(vec, 1.0 / math.sqrt(3.0)) for vec in raw]


def gram(values: list[Vec]) -> list[list[float]]:
    return [[dot(left, right) for right in values] for left in values]


def rank_symmetric(matrix: list[list[float]], tol: float = 1e-10) -> int:
    """Small Gaussian rank helper for exact-shaped Gram matrices."""
    rows = [row[:] for row in matrix]
    n_rows = len(rows)
    n_cols = len(rows[0])
    rank = 0
    for col in range(n_cols):
        pivot = None
        for row in range(rank, n_rows):
            if abs(rows[row][col]) > tol:
                pivot = row
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][col]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for row in range(n_rows):
            if row != rank and abs(rows[row][col]) > tol:
                factor = rows[row][col]
                rows[row] = [a - factor * b for a, b in zip(rows[row], rows[rank], strict=True)]
        rank += 1
    return rank


def forced_gram(size: int = 4) -> list[list[float]]:
    off_diagonal = -1.0 / (size - 1)
    return [
        [1.0 if i == j else off_diagonal for j in range(size)]
        for i in range(size)
    ]


def quat_mul(left: Quat, right: Quat) -> Quat:
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )


def quat_conj(q: Quat) -> Quat:
    a, b, c, d = q
    return (a, -b, -c, -d)


def quat_norm(q: Quat) -> float:
    return math.sqrt(sum(value * value for value in q))


def quat_to_vec(q: Quat) -> Vec:
    return (q[1], q[2], q[3])


def vec_to_quat(vec: Vec) -> Quat:
    return (0.0, vec[0], vec[1], vec[2])


def rotate(q: Quat, vec: Vec) -> Vec:
    return quat_to_vec(quat_mul(quat_mul(q, vec_to_quat(vec)), quat_conj(q)))


def rounded_quat(q: Quat) -> tuple[int, int, int, int]:
    return tuple(int(round(2.0 * value)) for value in q)  # type: ignore[return-value]


def rounded_vec(vec: Vec) -> tuple[int, int, int]:
    return tuple(int(round(math.sqrt(3.0) * value)) for value in vec)  # type: ignore[return-value]


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
    unique = {}
    for q in elements:
        unique[rounded_quat(q)] = q
    return list(unique.values())


def permutation_parity(perm: Perm) -> int:
    inversions = 0
    for i, j in itertools.combinations(range(len(perm)), 2):
        if perm[i] > perm[j]:
            inversions += 1
    return inversions % 2


def induced_permutation(q: Quat, vertices: list[Vec]) -> Perm:
    lookup = {rounded_vec(vec): idx for idx, vec in enumerate(vertices)}
    image = []
    for vec in vertices:
        key = rounded_vec(rotate(q, vec))
        image.append(lookup[key])
    return tuple(image)  # type: ignore[return-value]


def quat_key(q: Quat) -> tuple[int, int, int, int]:
    return rounded_quat(q)


def canonical_rotation_key(q: Quat) -> tuple[int, int, int, int]:
    key = quat_key(q)
    neg_key = tuple(-x for x in key)
    return min(key, neg_key)


def main() -> int:
    vertices = tetra_vertices()
    g_actual = gram(vertices)
    g_forced = forced_gram(4)
    barycenter = (0.0, 0.0, 0.0)
    for vec in vertices:
        barycenter = add(barycenter, vec)

    group = binary_tetrahedral_group()
    group_keys = {quat_key(q) for q in group}
    multiplication_closed = all(
        quat_key(quat_mul(left, right)) in group_keys
        for left in group
        for right in group
    )

    permutations = [induced_permutation(q, vertices) for q in group]
    unique_permutations = sorted(set(permutations))
    even_permutations = {
        tuple(perm)
        for perm in itertools.permutations(range(4))
        if permutation_parity(perm) == 0
    }
    rotation_keys = {canonical_rotation_key(q) for q in group}

    kernel = [q for q in group if induced_permutation(q, vertices) == (0, 1, 2, 3)]

    checks = {
        "forced_gram_matches_tetrahedron": all(
            close(g_actual[i][j], g_forced[i][j])
            for i in range(4)
            for j in range(4)
        ),
        "gram_has_rank_three": rank_symmetric(g_forced) == 3,
        "barycenter_zero": close(norm(barycenter), 0.0),
        "all_vertices_unit": all(close(norm(vec), 1.0) for vec in vertices),
        "all_edges_equal": len(
            {
                round(norm(sub(vertices[i], vertices[j])), 12)
                for i, j in itertools.combinations(range(4), 2)
            }
        )
        == 1,
        "antipodal_double_has_eight_vertices": len({rounded_vec(vec) for vec in vertices + [neg(v) for v in vertices]}) == 8,
        "binary_tetrahedral_group_has_order_24": len(group) == 24,
        "binary_tetrahedral_group_is_unit": all(close(quat_norm(q), 1.0) for q in group),
        "binary_tetrahedral_group_is_closed": multiplication_closed,
        "central_minus_one_present": (-2, 0, 0, 0) in group_keys,
        "central_minus_one_acts_trivially_on_visible_tetrahedron": induced_permutation((-1.0, 0.0, 0.0, 0.0), vertices)
        == (0, 1, 2, 3),
        "kernel_of_visible_action_is_plus_minus_one": {quat_key(q) for q in kernel}
        == {(2, 0, 0, 0), (-2, 0, 0, 0)},
        "visible_rotation_group_has_order_12": len(unique_permutations) == 12,
        "visible_rotation_group_is_A4": set(unique_permutations) == even_permutations,
        "su2_to_so3_is_two_to_one_on_this_group": len(rotation_keys) == 12,
    }

    payload = {
        "model": "tetrahedral_spinor_forcing",
        "theorem_layers": [
            "regular simplex forcing by Gram matrix",
            "antipodal double by free metric involution",
            "spinorial lift by binary tetrahedral group 2T -> A4",
        ],
        "tetrahedron": {
            "vertices": vertices,
            "barycenter": barycenter,
            "forced_gram": g_forced,
            "actual_gram": g_actual,
            "gram_rank": rank_symmetric(g_forced),
            "edge_length": norm(sub(vertices[0], vertices[1])),
        },
        "spinorial_lift": {
            "binary_tetrahedral_order": len(group),
            "visible_rotation_order": len(unique_permutations),
            "kernel_keys": sorted(quat_key(q) for q in kernel),
            "unique_visible_permutations": [list(perm) for perm in unique_permutations],
        },
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
