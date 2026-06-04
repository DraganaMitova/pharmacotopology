from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from typing import Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import (
    CouplingConstraint,
    CouplingDataset,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


COUPLING_NEGATIVE_CONTROL_KIND = "external_coupling_negative_controls_v0"

EXTERNAL_COUPLING_CONTROL_NAMES = (
    "external_shuffled_same_row_same_separation",
    "external_confidence_permuted",
    "external_cross_row_swapped",
    "external_random_long_range_same_count",
    "external_low_confidence_tail",
)


@dataclass(frozen=True)
class CouplingControlDataset:
    control_name: str
    control_kind: str
    dataset: CouplingDataset
    constraint_count: int
    coordinate_truth_used_to_build_constraints: bool = False
    native_truth_used_before_coupling_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "control_name": self.control_name,
            "control_kind": self.control_kind,
            "constraint_count": self.constraint_count,
            "coordinate_truth_used_to_build_constraints": (
                self.coordinate_truth_used_to_build_constraints
            ),
            "native_truth_used_before_coupling_selection": (
                self.native_truth_used_before_coupling_selection
            ),
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _constraint_id(control_name: str, row_id: str, i: int, j: int, rank: int) -> str:
    digest = hashlib.sha256(
        f"{control_name}:{row_id}:{i}:{j}:{rank}".encode("utf-8")
    ).hexdigest()
    return f"{control_name}_{digest[:12]}"


def _row_maps(
    rows: Sequence[RealCoordinateVisualRow],
) -> tuple[dict[str, RealCoordinateVisualRow], dict[str, tuple[CouplingConstraint, ...]]]:
    row_by_id = {row.row_id: row for row in rows}
    return row_by_id, {row.row_id: () for row in rows}


def _with_constraints(
    source: CouplingDataset,
    constraints: Sequence[CouplingConstraint],
) -> CouplingDataset:
    return replace(
        source,
        coordinate_truth_used_to_build_constraints=False,
        native_truth_used_before_coupling_selection=False,
        raw_sequence_exposed=False,
        constraints=tuple(constraints),
    )


def _bounded_pair(length: int, separation: int, seed: int) -> tuple[int, int]:
    separation = max(1, min(separation, length - 1))
    start_count = max(1, length - separation)
    i = seed % start_count + 1
    return (i, i + separation)


def _same_row_same_separation(
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
) -> tuple[CouplingConstraint, ...]:
    row_by_id, _ = _row_maps(rows)
    output: list[CouplingConstraint] = []
    for rank, constraint in enumerate(dataset.constraints, start=1):
        row = row_by_id[constraint.row_id]
        seed = constraint.i * 37 + constraint.j * 101 + rank * 17
        i, j = _bounded_pair(row.sequence_length, constraint.sequence_separation, seed)
        if (i, j) == constraint.pair():
            i, j = _bounded_pair(
                row.sequence_length,
                constraint.sequence_separation,
                seed + constraint.sequence_separation + 1,
            )
        output.append(
            replace(
                constraint,
                constraint_id=_constraint_id(
                    "external_shuffled_same_row_same_separation",
                    constraint.row_id,
                    i,
                    j,
                    rank,
                ),
                i=i,
                j=j,
                normalized_separation=round((j - i) / row.sequence_length, 6),
                coordinate_truth_used_to_build_constraint=False,
                native_truth_used_before_coupling_selection=False,
                structure_model_used=False,
            )
        )
    return tuple(output)


def _confidence_permuted(dataset: CouplingDataset) -> tuple[CouplingConstraint, ...]:
    constraints = tuple(dataset.constraints)
    if not constraints:
        return ()
    confidences = [constraint.confidence for constraint in constraints]
    shift = max(1, len(confidences) // 3)
    rotated = confidences[shift:] + confidences[:shift]
    return tuple(
        replace(
            constraint,
            constraint_id=_constraint_id(
                "external_confidence_permuted",
                constraint.row_id,
                constraint.i,
                constraint.j,
                rank,
            ),
            confidence=round(rotated[rank - 1], 6),
            coordinate_truth_used_to_build_constraint=False,
            native_truth_used_before_coupling_selection=False,
            structure_model_used=False,
        )
        for rank, constraint in enumerate(constraints, start=1)
    )


def _cross_row_swapped(
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
) -> tuple[CouplingConstraint, ...]:
    row_by_id = {row.row_id: row for row in rows}
    grouped = dataset.constraints_by_row_id()
    row_ids = [row.row_id for row in rows]
    output: list[CouplingConstraint] = []
    rank = 1
    for index, target_row_id in enumerate(row_ids):
        source_row_id = row_ids[(index + 1) % len(row_ids)]
        source_row = row_by_id[source_row_id]
        target_row = row_by_id[target_row_id]
        for constraint in grouped.get(source_row_id, ()):
            scaled_i = round(constraint.i / source_row.sequence_length * target_row.sequence_length)
            scaled_j = round(constraint.j / source_row.sequence_length * target_row.sequence_length)
            i = max(1, min(target_row.sequence_length - 1, scaled_i))
            j = max(i + 1, min(target_row.sequence_length, scaled_j))
            output.append(
                replace(
                    constraint,
                    row_id=target_row.row_id,
                    source_accession=target_row.source_accession,
                    constraint_id=_constraint_id(
                        "external_cross_row_swapped",
                        target_row.row_id,
                        i,
                        j,
                        rank,
                    ),
                    i=i,
                    j=j,
                    sequence_separation=j - i,
                    normalized_separation=round((j - i) / target_row.sequence_length, 6),
                    coordinate_truth_used_to_build_constraint=False,
                    native_truth_used_before_coupling_selection=False,
                    structure_model_used=False,
                )
            )
            rank += 1
    return tuple(output)


def _random_long_range_same_count(
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
) -> tuple[CouplingConstraint, ...]:
    row_by_id = {row.row_id: row for row in rows}
    output: list[CouplingConstraint] = []
    for rank, constraint in enumerate(dataset.constraints, start=1):
        row = row_by_id[constraint.row_id]
        min_separation = min(max(24, constraint.sequence_separation), row.sequence_length - 1)
        seed = int(
            hashlib.sha256(f"{constraint.constraint_id}:{rank}".encode("utf-8")).hexdigest()[:8],
            16,
        )
        i, j = _bounded_pair(row.sequence_length, min_separation, seed)
        output.append(
            replace(
                constraint,
                constraint_id=_constraint_id(
                    "external_random_long_range_same_count",
                    constraint.row_id,
                    i,
                    j,
                    rank,
                ),
                i=i,
                j=j,
                sequence_separation=j - i,
                normalized_separation=round((j - i) / row.sequence_length, 6),
                coordinate_truth_used_to_build_constraint=False,
                native_truth_used_before_coupling_selection=False,
                structure_model_used=False,
            )
        )
    return tuple(output)


def _low_confidence_tail(dataset: CouplingDataset) -> tuple[CouplingConstraint, ...]:
    return tuple(
        replace(
            constraint,
            constraint_id=_constraint_id(
                "external_low_confidence_tail",
                constraint.row_id,
                constraint.i,
                constraint.j,
                rank,
            ),
            confidence=round(max(0.001, constraint.confidence * 0.20), 6),
            coordinate_truth_used_to_build_constraint=False,
            native_truth_used_before_coupling_selection=False,
            structure_model_used=False,
        )
        for rank, constraint in enumerate(dataset.constraints, start=1)
    )


def generate_external_coupling_negative_controls(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
) -> Mapping[str, CouplingControlDataset]:
    builders = {
        "external_shuffled_same_row_same_separation": _same_row_same_separation(
            rows,
            dataset,
        ),
        "external_confidence_permuted": _confidence_permuted(dataset),
        "external_cross_row_swapped": _cross_row_swapped(rows, dataset),
        "external_random_long_range_same_count": _random_long_range_same_count(
            rows,
            dataset,
        ),
        "external_low_confidence_tail": _low_confidence_tail(dataset),
    }
    return {
        name: CouplingControlDataset(
            control_name=name,
            control_kind=COUPLING_NEGATIVE_CONTROL_KIND,
            dataset=_with_constraints(dataset, constraints),
            constraint_count=len(constraints),
        )
        for name, constraints in builders.items()
    }
