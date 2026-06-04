from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent


CLOSURE_GEOMETRY_KIND = "sequence_only_closure_geometry_v1"


@dataclass(frozen=True)
class ClosureCompatibility:
    left_event_id: str
    right_event_id: str
    row_id: str
    compatibility_label: str
    overlap_residue_count: int
    overlap_ratio: float
    crossing_arcs: bool
    steric_risk_proxy: float

    def to_dict(self) -> dict[str, object]:
        return {
            "row_id": self.row_id,
            "left_event_id": self.left_event_id,
            "right_event_id": self.right_event_id,
            "compatibility_label": self.compatibility_label,
            "overlap_residue_count": self.overlap_residue_count,
            "overlap_ratio": self.overlap_ratio,
            "crossing_arcs": self.crossing_arcs,
            "steric_risk_proxy": self.steric_risk_proxy,
        }


def rounded_unit(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def event_residue_indices(event: NucleusClosureEvent) -> frozenset[int]:
    return frozenset(
        tuple(range(event.segment_a_start, event.segment_a_end + 1))
        + tuple(range(event.segment_b_start, event.segment_b_end + 1))
    )


def event_arc(event: NucleusClosureEvent) -> tuple[float, float]:
    left = (event.segment_a_start + event.segment_a_end) / 2
    right = (event.segment_b_start + event.segment_b_end) / 2
    return (left, right)


def arcs_cross(left: NucleusClosureEvent, right: NucleusClosureEvent) -> bool:
    left_a, left_b = event_arc(left)
    right_a, right_b = event_arc(right)
    return (left_a < right_a < left_b < right_b) or (
        right_a < left_a < right_b < left_b
    )


def event_overlap_ratio(
    left: NucleusClosureEvent,
    right: NucleusClosureEvent,
) -> tuple[int, float]:
    left_residues = event_residue_indices(left)
    right_residues = event_residue_indices(right)
    shared = len(left_residues & right_residues)
    denominator = max(1, min(len(left_residues), len(right_residues)))
    return shared, rounded_unit(shared / denominator)


def steric_risk_proxy(
    left: NucleusClosureEvent,
    right: NucleusClosureEvent,
) -> float:
    span_pressure = max(0.0, left.normalized_span - 0.55) + max(
        0.0, right.normalized_span - 0.55
    )
    geometry_pressure = left.geometry_violation_cost + right.geometry_violation_cost
    entropy_pressure = (left.loop_entropy_cost + right.loop_entropy_cost) / 2
    return rounded_unit(
        0.44 * span_pressure
        + 0.28 * geometry_pressure
        + 0.18 * entropy_pressure
        + (0.10 if arcs_cross(left, right) else 0.0)
    )


def closure_compatibility(
    left: NucleusClosureEvent,
    right: NucleusClosureEvent,
) -> ClosureCompatibility:
    shared, overlap_ratio = event_overlap_ratio(left, right)
    crossing = arcs_cross(left, right)
    risk = steric_risk_proxy(left, right)
    if left.row_id != right.row_id:
        label = "different_row_not_compared"
    elif overlap_ratio >= 0.50:
        label = "overlapping"
    elif crossing and risk >= 0.22:
        label = "topologically_conflicting"
    elif risk >= 0.34:
        label = "sterically_risky"
    elif shared > 0:
        label = "competing"
    elif (
        abs(left.normalized_span - right.normalized_span) >= 0.55
        and max(left.normalized_span, right.normalized_span) >= 0.70
    ):
        label = "domain_inconsistent"
    else:
        label = "compatible"
    return ClosureCompatibility(
        left_event_id=left.event_id,
        right_event_id=right.event_id,
        row_id=left.row_id if left.row_id == right.row_id else "",
        compatibility_label=label,
        overlap_residue_count=shared,
        overlap_ratio=overlap_ratio,
        crossing_arcs=crossing,
        steric_risk_proxy=risk,
    )


def compatibility_rows(
    selected_events: Sequence[NucleusClosureEvent],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in selected_events:
        by_row.setdefault(event.row_id, []).append(event)
    for row_id in sorted(by_row):
        row_events = sorted(
            by_row[row_id],
            key=lambda event: (
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for left_index, left in enumerate(row_events):
            for right in row_events[left_index + 1 :]:
                rows.append(closure_compatibility(left, right).to_dict())
    return rows

