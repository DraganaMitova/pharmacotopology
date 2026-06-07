from __future__ import annotations

from dataclasses import asdict, dataclass

from pharmacotopology.folding_physical_state import PhysicalClosureState


FUTURE_FRUSTRATION_GATE_KIND = "coarse_future_frustration_gate_v1"
FUTURE_FRUSTRATION_LIMIT = 0.26


@dataclass(frozen=True)
class FutureFrustrationAssessment:
    state_id: str
    event_id: str
    row_id: str
    future_frustration_score: float
    future_closure_path_preserved: bool
    rejection_reason: str
    native_truth_used_before_future_gate: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def assess_future_frustration(
    state: PhysicalClosureState,
    *,
    limit: float = FUTURE_FRUSTRATION_LIMIT,
) -> FutureFrustrationAssessment:
    preserved = state.future_frustration_score <= limit
    return FutureFrustrationAssessment(
        state_id=state.state_id,
        event_id=state.event_id,
        row_id=state.row_id,
        future_frustration_score=state.future_frustration_score,
        future_closure_path_preserved=preserved,
        rejection_reason="passed" if preserved else "future_frustration_rejection",
    )
