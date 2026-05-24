#!/usr/bin/env python3
"""Finite checks for the bipolar tetrahedral-potential hypothesis.

The hypothesis tested here is deliberately modest:

    binary branch + free involution + symmetric normalization

does not merely duplicate a four-probe potential tetrahedron. It realizes it as
two antipodal tetrahedra on S^2. The quotient by the involution is the original
four-probe tetrahedral packet.
"""

from __future__ import annotations

import itertools
import json
import math
from collections.abc import Iterable


Vec = tuple[float, float, float]
Face = tuple[int, ...]
SignedFace = tuple[tuple[int, int], ...]

TOL = 1e-12


def dot(left: Vec, right: Vec) -> float:
    return sum(a * b for a, b in zip(left, right, strict=True))


def norm(vec: Vec) -> float:
    return math.sqrt(dot(vec, vec))


def scale(vec: Vec, factor: float) -> Vec:
    return tuple(factor * x for x in vec)  # type: ignore[return-value]


def neg(vec: Vec) -> Vec:
    return scale(vec, -1.0)


def sub(left: Vec, right: Vec) -> Vec:
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def close(left: float, right: float, tol: float = TOL) -> bool:
    return abs(left - right) <= tol


def tetra_positive() -> list[Vec]:
    raw = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    return [scale(vec, 1.0 / math.sqrt(3.0)) for vec in raw]


def pairwise(values: list[Vec]) -> list[tuple[int, int, float, float]]:
    return [
        (i, j, dot(values[i], values[j]), norm(sub(values[i], values[j])))
        for i, j in itertools.combinations(range(len(values)), 2)
    ]


def powerset(values: Iterable[int]) -> list[Face]:
    items = tuple(values)
    return [
        tuple(combo)
        for cardinality in range(len(items) + 1)
        for combo in itertools.combinations(items, cardinality)
    ]


def signed_powerset(values: Iterable[tuple[int, int]]) -> list[SignedFace]:
    items = tuple(values)
    return [
        tuple(combo)
        for cardinality in range(len(items) + 1)
        for combo in itertools.combinations(items, cardinality)
    ]


def faces(face: Face) -> list[Face]:
    return powerset(face)


def signed_faces(face: SignedFace) -> list[SignedFace]:
    return [
        tuple(combo)
        for cardinality in range(len(face) + 1)
        for combo in itertools.combinations(face, cardinality)
    ]


def is_simplicial_complex(family: set[Face]) -> bool:
    return all(subface in family for face in family for subface in faces(face))


def is_signed_simplicial_complex(family: set[SignedFace]) -> bool:
    return all(subface in family for face in family for subface in signed_faces(face))


def base_potential(required_size: int) -> set[Face]:
    return {face for face in powerset(range(4)) if len(face) < required_size}


def project_signed(face: SignedFace) -> Face:
    return tuple(sorted({idx for idx, _sign in face}))


def signed_lift(base_family: set[Face]) -> set[SignedFace]:
    signed_vertices = tuple((idx, sign) for idx in range(4) for sign in (-1, 1))
    return {face for face in signed_powerset(signed_vertices) if project_signed(face) in base_family}


def involution(face: SignedFace) -> SignedFace:
    return tuple(sorted((idx, -sign) for idx, sign in face))


def minimal_non_faces_signed(family: set[SignedFace]) -> list[SignedFace]:
    signed_vertices = tuple((idx, sign) for idx in range(4) for sign in (-1, 1))
    non_faces = [face for face in signed_powerset(signed_vertices) if face not in family]
    return sorted(
        face
        for face in non_faces
        if not any(subface not in family for subface in signed_faces(face) if subface != face)
    )


def signed_permutation_matrices() -> list[tuple[tuple[int, ...], tuple[int, int, int]]]:
    return [
        (perm, signs)
        for perm in itertools.permutations(range(3))
        for signs in itertools.product((-1, 1), repeat=3)
    ]


def apply_signed_permutation(vec: Vec, perm: tuple[int, ...], signs: tuple[int, int, int]) -> Vec:
    raw = (vec[perm[0]] * signs[0], vec[perm[1]] * signs[1], vec[perm[2]] * signs[2])
    return raw


def rounded(vec: Vec) -> tuple[int, int, int]:
    return tuple(int(round(x * math.sqrt(3.0))) for x in vec)  # type: ignore[return-value]


def main() -> int:
    positive = tetra_positive()
    negative = [neg(vec) for vec in positive]
    union = positive + negative

    union_set = {rounded(vec) for vec in union}
    cube_set = set(itertools.product((-1, 1), repeat=3))

    signed_mats = signed_permutation_matrices()
    preserves_union = []
    preserves_positive = []
    swaps_positive_negative = []
    positive_set = {rounded(vec) for vec in positive}
    negative_set = {rounded(vec) for vec in negative}

    for perm, signs in signed_mats:
        image_union = {rounded(apply_signed_permutation(vec, perm, signs)) for vec in union}
        image_positive = {rounded(apply_signed_permutation(vec, perm, signs)) for vec in positive}
        if image_union == union_set:
            preserves_union.append((perm, signs))
        if image_positive == positive_set:
            preserves_positive.append((perm, signs))
        if image_positive == negative_set:
            swaps_positive_negative.append((perm, signs))

    potential_payloads = []
    for required_size in (1, 2, 3, 4):
        base = base_potential(required_size)
        lift = signed_lift(base)
        minimal_signed = minimal_non_faces_signed(lift)
        projected_minimal = {project_signed(face) for face in minimal_signed}
        potential_payloads.append(
            {
                "required_size": required_size,
                "base_f_vector": [
                    sum(1 for face in base if len(face) == cardinality)
                    for cardinality in range(5)
                ],
                "signed_lift_size": len(lift),
                "minimal_signed_non_faces": [
                    [[idx, sign] for idx, sign in face] for face in minimal_signed[:12]
                ],
                "minimal_signed_non_face_count": len(minimal_signed),
                "projected_minimal_non_faces": [list(face) for face in sorted(projected_minimal)],
                "checks": {
                    "base_is_simplicial": is_simplicial_complex(base),
                    "lift_is_simplicial": is_signed_simplicial_complex(lift),
                    "lift_is_involution_invariant": all(involution(face) in lift for face in lift),
                    "minimal_non_faces_project_to_required_size": all(
                        len(face) == required_size for face in projected_minimal
                    ),
                },
            }
        )

    checks = {
        "positive_vertices_on_unit_sphere": all(close(norm(vec), 1.0) for vec in positive),
        "negative_is_antipodal_image": all(
            all(close(a, b) for a, b in zip(negative[idx], neg(positive[idx]), strict=True))
            for idx in range(4)
        ),
        "positive_is_regular_tetrahedron": all(
            close(product, -1.0 / 3.0) and close(distance, math.sqrt(8.0 / 3.0))
            for _i, _j, product, distance in pairwise(positive)
        ),
        "negative_is_regular_tetrahedron": all(
            close(product, -1.0 / 3.0) and close(distance, math.sqrt(8.0 / 3.0))
            for _i, _j, product, distance in pairwise(negative)
        ),
        "union_is_cube_vertices_on_sphere": union_set == cube_set,
        "central_involution_has_no_fixed_vertex": all(positive[idx] != negative[idx] for idx in range(4)),
        "signed_symmetry_group_has_order_48": len(preserves_union) == 48,
        "branch_preserving_subgroup_has_order_24": len(preserves_positive) == 24,
        "branch_swapping_coset_has_order_24": len(swaps_positive_negative) == 24,
        "all_potential_lifts_pass": all(all(payload["checks"].values()) for payload in potential_payloads),
    }

    payload = {
        "model": "bipolar_tetrahedral_potential",
        "interpretation": {
            "claim": "binary free involution plus symmetric normalization realizes a tetrahedral potential as two antipodal tetrahedra on S^2",
            "caution": "the number 2 alone does not force this; the involutive antipodal normalization is an explicit hypothesis",
        },
        "geometry": {
            "positive_vertices": positive,
            "negative_vertices": negative,
            "positive_pairwise_dot_distance": pairwise(positive),
            "cross_dot_values": sorted(
                {round(dot(pos, negv), 12) for pos in positive for negv in negative}
            ),
            "union_rounded_vertices": sorted(union_set),
        },
        "symmetry": {
            "signed_permutation_group_order": len(preserves_union),
            "branch_preserving_order": len(preserves_positive),
            "branch_swapping_order": len(swaps_positive_negative),
        },
        "potential_lifts": potential_payloads,
        "checks": checks,
        "passed": all(checks.values()),
    }

    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
