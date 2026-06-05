from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Mapping, Sequence

from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


EVOLUTIONARY_COUPLING_LAYER_KIND = (
    "locked_safe_coupling_constraint_layer_v1"
)
COUPLING_CONSTRAINT_KIND = "safe_residue_pair_coupling_constraint_v1"
COUPLING_ASSESSMENT_KIND = "coupling_preserved_closure_state_v1"


@dataclass(frozen=True)
class CouplingConstraint:
    row_id: str
    source_accession: str
    constraint_id: str
    i: int
    j: int
    sequence_separation: int
    normalized_separation: float
    confidence: float
    constraint_class: str
    source_kind: str
    coordinate_truth_used_to_build_constraint: bool
    native_truth_used_before_coupling_selection: bool = False
    structure_model_used: bool = False
    raw_sequence_exposed: bool = False
    raw_score: float = 0.0
    apc_corrected_score: float = 0.0
    rank: int = 0
    rank_fraction: float = 0.0

    def pair(self) -> tuple[int, int]:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CouplingDataset:
    layer_kind: str
    constraint_kind: str
    source_benchmark_file: str
    source_constraint_kind: str
    coupling_source_kind: str
    coordinate_truth_used_to_build_constraints: bool
    native_truth_used_before_coupling_selection: bool
    external_evolutionary_couplings_used: bool
    raw_sequence_exposed: bool
    constraints: tuple[CouplingConstraint, ...]
    structure_model_used_before_coupling_selection: bool = False

    @property
    def per_constraint_coordinate_truth_used(self) -> bool:
        return any(
            constraint.coordinate_truth_used_to_build_constraint
            for constraint in self.constraints
        )

    @property
    def coordinate_truth_tainted(self) -> bool:
        return (
            self.coordinate_truth_used_to_build_constraints
            or self.per_constraint_coordinate_truth_used
        )

    @property
    def per_constraint_native_truth_used(self) -> bool:
        return any(
            constraint.native_truth_used_before_coupling_selection
            for constraint in self.constraints
        )

    @property
    def native_truth_tainted(self) -> bool:
        return (
            self.native_truth_used_before_coupling_selection
            or self.per_constraint_native_truth_used
        )

    @property
    def per_constraint_structure_model_used(self) -> bool:
        return any(constraint.structure_model_used for constraint in self.constraints)

    @property
    def structure_model_tainted(self) -> bool:
        return (
            self.structure_model_used_before_coupling_selection
            or self.per_constraint_structure_model_used
        )

    @property
    def oracle_constraint_control(self) -> bool:
        return (
            self.coordinate_truth_tainted
            or self.native_truth_tainted
            or self.structure_model_tainted
        )

    @cached_property
    def _constraints_by_row_id(self) -> dict[str, tuple[CouplingConstraint, ...]]:
        grouped: dict[str, list[CouplingConstraint]] = {}
        for constraint in self.constraints:
            grouped.setdefault(constraint.row_id, []).append(constraint)
        return {
            row_id: tuple(
                sorted(
                    row_constraints,
                    key=lambda item: (-item.confidence, item.i, item.j),
                )
            )
            for row_id, row_constraints in grouped.items()
        }

    def constraints_by_row_id(self) -> dict[str, tuple[CouplingConstraint, ...]]:
        return self._constraints_by_row_id


@dataclass(frozen=True)
class CouplingClosureAssessment:
    row_id: str
    source_accession: str
    event_id: str
    direct_coupling_count: int
    direct_coupling_confidence: float
    direct_support_score: float
    future_coupling_count: int
    future_preserved_count: int
    future_preservation_score: float
    blocked_future_count: int
    blocked_future_confidence: float
    blocked_future_pressure: float
    coupling_selectivity_score: float
    constraint_pairs_total: int
    coordinate_truth_used_to_build_constraints: bool
    native_truth_used_before_coupling_selection: bool
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _score(value: float) -> float:
    return round(value, 6)


def _constraint_id(row_id: str, i: int, j: int) -> str:
    encoded = f"{row_id}:{i}:{j}".encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _constraint_from_raw(
    raw: Mapping[str, Any],
    *,
    source_kind: str,
) -> CouplingConstraint:
    i = int(raw["i"])
    j = int(raw["j"])
    if j <= i:
        raise ValueError("coupling constraints must use i < j")
    confidence = _rounded(float(raw["confidence"]))
    row_id = str(raw["row_id"])
    sequence_separation = j - i
    return CouplingConstraint(
        row_id=row_id,
        source_accession=str(raw["source_accession"]),
        constraint_id=str(raw.get("constraint_id") or _constraint_id(row_id, i, j)),
        i=i,
        j=j,
        sequence_separation=sequence_separation,
        normalized_separation=_rounded(float(raw["normalized_separation"])),
        confidence=confidence,
        constraint_class=str(raw.get("constraint_class", "coupling_anchor")),
        source_kind=source_kind,
        coordinate_truth_used_to_build_constraint=bool(
            raw.get("coordinate_truth_used_to_build_constraint", False)
        ),
        native_truth_used_before_coupling_selection=bool(
            raw.get("native_truth_used_before_coupling_selection", False)
        ),
        structure_model_used=bool(raw.get("structure_model_used", False)),
        raw_sequence_exposed=bool(raw.get("raw_sequence_exposed", False)),
        raw_score=_score(float(raw.get("raw_score", 0.0))),
        apc_corrected_score=_score(float(raw.get("apc_corrected_score", 0.0))),
        rank=int(raw.get("rank", 0)),
        rank_fraction=_rounded(float(raw.get("rank_fraction", 0.0))),
    )


def load_coupling_dataset(path: Path) -> CouplingDataset:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{path} must contain a JSON object")
    constraints_raw = parsed.get("constraints", [])
    if not isinstance(constraints_raw, list):
        raise ValueError("constraints must be a list")
    source_kind = str(parsed.get("coupling_source_kind", "unknown"))
    constraints = tuple(
        _constraint_from_raw(raw, source_kind=source_kind)
        for raw in constraints_raw
        if isinstance(raw, Mapping)
    )
    return CouplingDataset(
        layer_kind=str(parsed.get("layer_kind", EVOLUTIONARY_COUPLING_LAYER_KIND)),
        constraint_kind=str(parsed.get("constraint_kind", COUPLING_CONSTRAINT_KIND)),
        source_benchmark_file=str(parsed.get("source_benchmark_file", "")),
        source_constraint_kind=str(parsed.get("source_constraint_kind", "")),
        coupling_source_kind=source_kind,
        coordinate_truth_used_to_build_constraints=bool(
            parsed.get("coordinate_truth_used_to_build_constraints", False)
        ),
        native_truth_used_before_coupling_selection=bool(
            parsed.get("native_truth_used_before_coupling_selection", False)
        ),
        external_evolutionary_couplings_used=bool(
            parsed.get("external_evolutionary_couplings_used", False)
        ),
        raw_sequence_exposed=bool(parsed.get("raw_sequence_exposed", False)),
        constraints=constraints,
        structure_model_used_before_coupling_selection=bool(
            parsed.get("structure_model_used_before_coupling_selection", False)
            or parsed.get("structure_model_used", False)
        ),
    )


def validate_coupling_dataset(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
) -> None:
    row_by_id = {row.row_id: row for row in rows}
    if dataset.layer_kind != EVOLUTIONARY_COUPLING_LAYER_KIND:
        raise ValueError(f"unsupported coupling layer kind: {dataset.layer_kind}")
    if dataset.constraint_kind != COUPLING_CONSTRAINT_KIND:
        raise ValueError(f"unsupported constraint kind: {dataset.constraint_kind}")
    if dataset.raw_sequence_exposed:
        raise ValueError("coupling dataset must not expose raw sequence text")
    for constraint in dataset.constraints:
        if constraint.raw_sequence_exposed:
            raise ValueError(
                f"coupling constraint must not expose raw sequence text: "
                f"{constraint.constraint_id}"
            )
        if constraint.structure_model_used:
            raise ValueError(
                f"coupling constraint must not use structure models before "
                f"selection: {constraint.constraint_id}"
            )
        row = row_by_id.get(constraint.row_id)
        if row is None:
            raise ValueError(f"constraint row not in benchmark: {constraint.row_id}")
        if constraint.source_accession != row.source_accession:
            raise ValueError(f"source accession mismatch for {constraint.row_id}")
        if constraint.i < 1 or constraint.j > row.sequence_length:
            raise ValueError(f"constraint outside sequence bounds: {constraint}")
        if constraint.sequence_separation != constraint.j - constraint.i:
            raise ValueError(f"invalid constraint separation: {constraint}")
        expected = _rounded((constraint.j - constraint.i) / row.sequence_length)
        if abs(constraint.normalized_separation - expected) > 0.000001:
            raise ValueError(f"invalid normalized separation: {constraint}")


def _event_residues(event: NucleusClosureEvent) -> set[int]:
    return set(range(event.segment_a_start, event.segment_a_end + 1)) | set(
        range(event.segment_b_start, event.segment_b_end + 1)
    )


def _overlap_ratio(left: NucleusClosureEvent, right: NucleusClosureEvent) -> float:
    left_residues = _event_residues(left)
    right_residues = _event_residues(right)
    return len(left_residues & right_residues) / max(
        1,
        min(len(left_residues), len(right_residues)),
    )


def _arcs_cross(left: NucleusClosureEvent, right: NucleusClosureEvent) -> bool:
    return (
        left.segment_a_start < right.segment_a_start < left.segment_b_start < right.segment_b_start
    ) or (
        right.segment_a_start < left.segment_a_start < right.segment_b_start < left.segment_b_start
    )


def compatible_future_event(
    locked: NucleusClosureEvent,
    candidate: NucleusClosureEvent,
) -> bool:
    if locked.event_id == candidate.event_id:
        return False
    if _overlap_ratio(locked, candidate) >= 0.75:
        return False
    if _arcs_cross(locked, candidate) and _overlap_ratio(locked, candidate) > 0.0:
        return False
    if candidate.geometry_violation_cost > 0.55:
        return False
    if candidate.frustration_cost > 0.72:
        return False
    return True


def _confidence_sum(constraints: Sequence[CouplingConstraint]) -> float:
    return sum(constraint.confidence for constraint in constraints)


def assess_coupling_closure(
    *,
    event: NucleusClosureEvent,
    row_constraints: Sequence[CouplingConstraint],
    row_events: Sequence[NucleusClosureEvent],
    row_events_by_pair: Mapping[tuple[int, int], Sequence[NucleusClosureEvent]] | None = None,
    coordinate_truth_used_to_build_constraints: bool,
    native_truth_used_before_coupling_selection: bool,
) -> CouplingClosureAssessment:
    region_pairs = set(event.candidate_region_pairs())
    event_residues = _event_residues(event)
    direct = tuple(
        constraint for constraint in row_constraints if constraint.pair() in region_pairs
    )
    remaining = tuple(
        constraint for constraint in row_constraints if constraint.pair() not in region_pairs
    )
    future_preserved: list[CouplingConstraint] = []
    blocked: list[CouplingConstraint] = []
    for constraint in remaining:
        if row_events_by_pair is None:
            covering_candidates = tuple(
                candidate
                for candidate in row_events
                if constraint.pair() in candidate.candidate_region_pairs()
            )
        else:
            covering_candidates = tuple(row_events_by_pair.get(constraint.pair(), ()))
        covering_future = [
            candidate
            for candidate in covering_candidates
            if compatible_future_event(event, candidate)
        ]
        if covering_future:
            future_preserved.append(constraint)
        elif constraint.i in event_residues or constraint.j in event_residues:
            blocked.append(constraint)

    direct_confidence = _confidence_sum(direct)
    future_confidence = _confidence_sum(future_preserved)
    blocked_confidence = _confidence_sum(blocked)
    all_confidence = max(1.0, _confidence_sum(row_constraints))
    future_pool = max(1, len(remaining))
    direct_support_score = _rounded(direct_confidence / 3.0)
    future_preservation_score = _rounded(
        (direct_confidence + future_confidence) / all_confidence
    )
    blocked_future_pressure = _rounded(blocked_confidence / all_confidence)
    coupling_selectivity_score = _score(
        0.54 * direct_support_score
        + 0.34 * future_preservation_score
        + 0.16 * event.registry_support
        + 0.12 * event.contact_cluster_gain
        + 0.08 * event.secondary_structure_compatibility
        - 0.30 * blocked_future_pressure
        - 0.10 * event.loop_entropy_cost
        - 0.08 * event.geometry_violation_cost
        - 0.06 * event.frustration_cost
    )
    return CouplingClosureAssessment(
        row_id=event.row_id,
        source_accession=event.source_accession,
        event_id=event.event_id,
        direct_coupling_count=len(direct),
        direct_coupling_confidence=_score(direct_confidence),
        direct_support_score=direct_support_score,
        future_coupling_count=future_pool,
        future_preserved_count=len(future_preserved),
        future_preservation_score=future_preservation_score,
        blocked_future_count=len(blocked),
        blocked_future_confidence=_score(blocked_confidence),
        blocked_future_pressure=blocked_future_pressure,
        coupling_selectivity_score=coupling_selectivity_score,
        constraint_pairs_total=len(row_constraints),
        coordinate_truth_used_to_build_constraints=(
            coordinate_truth_used_to_build_constraints
        ),
        native_truth_used_before_coupling_selection=(
            native_truth_used_before_coupling_selection
        ),
    )


def assess_coupling_closures(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    events: Sequence[NucleusClosureEvent],
    dataset: CouplingDataset,
) -> tuple[CouplingClosureAssessment, ...]:
    validate_coupling_dataset(rows=rows, dataset=dataset)
    constraints_by_row = dataset.constraints_by_row_id()
    events_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in events:
        events_by_row.setdefault(event.row_id, []).append(event)
    pairs_by_row = {
        row_id: {constraint.pair() for constraint in row_constraints}
        for row_id, row_constraints in constraints_by_row.items()
    }
    events_by_row_pair: dict[str, dict[tuple[int, int], list[NucleusClosureEvent]]] = {}
    for row_id, row_events in events_by_row.items():
        wanted_pairs = pairs_by_row.get(row_id, set())
        pair_map: dict[tuple[int, int], list[NucleusClosureEvent]] = {}
        for event in row_events:
            for pair in event.candidate_region_pairs():
                if pair in wanted_pairs:
                    pair_map.setdefault(pair, []).append(event)
        events_by_row_pair[row_id] = pair_map
    assessments: list[CouplingClosureAssessment] = []
    for event in events:
        assessments.append(
            assess_coupling_closure(
                event=event,
                row_constraints=constraints_by_row.get(event.row_id, ()),
                row_events=tuple(events_by_row.get(event.row_id, ())),
                row_events_by_pair=events_by_row_pair.get(event.row_id, {}),
                coordinate_truth_used_to_build_constraints=(
                    dataset.coordinate_truth_tainted
                ),
                native_truth_used_before_coupling_selection=(
                    dataset.native_truth_tainted
                ),
            )
        )
    return tuple(assessments)
