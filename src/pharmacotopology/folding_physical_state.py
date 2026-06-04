from __future__ import annotations

from dataclasses import asdict, dataclass

from pharmacotopology.folding_burial_frustration import BurialFrustrationPacket
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent


PHYSICAL_CLOSURE_STATE_KIND = "coarse_sequence_only_physical_closure_state_v1"


@dataclass(frozen=True)
class PhysicalClosureState:
    state_id: str
    row_id: str
    source_accession: str
    sequence_hash: str
    event_id: str
    segment_a_start: int
    segment_a_end: int
    segment_b_start: int
    segment_b_end: int
    normalized_span: float
    loop_strain: float
    steric_clash_score: float
    burial_gain: float
    unsatisfied_polar_penalty: float
    future_frustration_score: float
    physical_state_score: float
    state_build_success: bool
    state_build_failure_reason: str
    native_contact_count_after_scoring: int
    native_long_range_contact_count_after_scoring: int
    native_truth_used_before_physical_scoring: bool = False
    native_label_attached_after_physical_scoring: bool = True
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(value, 6)


def physical_state_score(
    *,
    event: NucleusClosureEvent,
    burial: BurialFrustrationPacket,
) -> float:
    return _rounded(
        burial.burial_gain
        - 0.75 * burial.loop_strain
        - 0.85 * burial.steric_clash_score
        - 0.75 * burial.unsatisfied_polar_penalty
        - 0.65 * burial.future_frustration_score
        + 0.18 * event.registry_support
        + 0.12 * event.contact_cluster_gain
    )


def physical_state_from_event(
    *,
    event: NucleusClosureEvent,
    burial: BurialFrustrationPacket,
) -> PhysicalClosureState:
    return PhysicalClosureState(
        state_id=f"{event.row_id}:{event.event_id}:physical",
        row_id=event.row_id,
        source_accession=event.source_accession,
        sequence_hash=event.sequence_hash,
        event_id=event.event_id,
        segment_a_start=event.segment_a_start,
        segment_a_end=event.segment_a_end,
        segment_b_start=event.segment_b_start,
        segment_b_end=event.segment_b_end,
        normalized_span=event.normalized_span,
        loop_strain=burial.loop_strain,
        steric_clash_score=burial.steric_clash_score,
        burial_gain=burial.burial_gain,
        unsatisfied_polar_penalty=burial.unsatisfied_polar_penalty,
        future_frustration_score=burial.future_frustration_score,
        physical_state_score=physical_state_score(event=event, burial=burial),
        state_build_success=True,
        state_build_failure_reason="",
        native_contact_count_after_scoring=event.native_contact_count_after_scoring,
        native_long_range_contact_count_after_scoring=(
            event.native_long_range_contact_count_after_scoring
        ),
    )

