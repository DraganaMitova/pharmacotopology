from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent


NUCLEUS_DECOY_FALSIFICATION_KIND = "matched_nucleus_decoy_falsification_v1"


@dataclass(frozen=True)
class NucleusDecoyMatch:
    row_id: str
    source_accession: str
    real_event_id: str
    decoy_event_id: str
    real_normalized_span: float
    decoy_normalized_span: float
    real_hydrophobic_burial_gain: float
    decoy_hydrophobic_burial_gain: float
    real_contact_cluster_gain: float
    decoy_contact_cluster_gain: float
    real_native_positive_after_scoring: bool
    decoy_native_positive_after_scoring: bool
    native_truth_used_before_decoy_matching: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def decoy_distance(
    real: NucleusClosureEvent,
    candidate: NucleusClosureEvent,
) -> tuple[float, str]:
    distance = (
        1.00 * abs(real.normalized_span - candidate.normalized_span)
        + 0.45
        * abs(real.hydrophobic_burial_gain - candidate.hydrophobic_burial_gain)
        + 0.35 * abs(real.contact_cluster_gain - candidate.contact_cluster_gain)
        + 0.20 * abs(real.registry_support - candidate.registry_support)
        + 0.10 * abs(real.sequence_span - candidate.sequence_span) / max(
            real.sequence_length,
            1,
        )
    )
    return (round(distance, 6), candidate.event_id)


def matched_decoys_for_selected_events(
    *,
    selected_events: Sequence[NucleusClosureEvent],
    candidate_events: Sequence[NucleusClosureEvent],
) -> tuple[NucleusDecoyMatch, ...]:
    selected_ids = {event.event_id for event in selected_events}
    candidates_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in candidate_events:
        if event.event_id not in selected_ids:
            candidates_by_row.setdefault(event.row_id, []).append(event)

    matches: list[NucleusDecoyMatch] = []
    for real in selected_events:
        row_candidates = candidates_by_row.get(real.row_id, [])
        if row_candidates:
            decoy = min(row_candidates, key=lambda candidate: decoy_distance(real, candidate))
        else:
            decoy = real
        matches.append(
            NucleusDecoyMatch(
                row_id=real.row_id,
                source_accession=real.source_accession,
                real_event_id=real.event_id,
                decoy_event_id=decoy.event_id,
                real_normalized_span=real.normalized_span,
                decoy_normalized_span=decoy.normalized_span,
                real_hydrophobic_burial_gain=real.hydrophobic_burial_gain,
                decoy_hydrophobic_burial_gain=decoy.hydrophobic_burial_gain,
                real_contact_cluster_gain=real.contact_cluster_gain,
                decoy_contact_cluster_gain=decoy.contact_cluster_gain,
                real_native_positive_after_scoring=(
                    real.native_contact_count_after_scoring > 0
                ),
                decoy_native_positive_after_scoring=(
                    decoy.native_contact_count_after_scoring > 0
                ),
            )
        )
    return tuple(matches)


def decoy_native_overlap_rate(matches: Sequence[NucleusDecoyMatch]) -> float:
    if not matches:
        return 0.0
    return round(
        sum(1 for match in matches if match.decoy_native_positive_after_scoring)
        / len(matches),
        6,
    )


def real_native_positive_rate(matches: Sequence[NucleusDecoyMatch]) -> float:
    if not matches:
        return 0.0
    return round(
        sum(1 for match in matches if match.real_native_positive_after_scoring)
        / len(matches),
        6,
    )


def real_vs_decoy_enrichment_ratio(matches: Sequence[NucleusDecoyMatch]) -> float:
    decoy_rate = decoy_native_overlap_rate(matches)
    real_rate = real_native_positive_rate(matches)
    if decoy_rate == 0.0:
        return 0.0
    return round(real_rate / decoy_rate, 6)

