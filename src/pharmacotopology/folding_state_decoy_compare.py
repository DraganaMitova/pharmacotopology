from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

from pharmacotopology.folding_nucleus_decoy_falsification import NucleusDecoyMatch
from pharmacotopology.folding_physical_state import PhysicalClosureState


PHYSICAL_STATE_DECOY_COMPARE_KIND = "physical_state_matched_decoy_compare_v1"


@dataclass(frozen=True)
class PhysicalStateDecoyComparison:
    row_id: str
    source_accession: str
    real_state_id: str
    decoy_state_id: str
    real_event_id: str
    decoy_event_id: str
    real_physical_state_score: float
    decoy_physical_state_score: float
    real_beats_decoy_by_score: bool
    real_native_positive_after_scoring: bool
    decoy_native_positive_after_scoring: bool
    native_truth_used_before_physical_scoring: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def physical_state_decoy_comparisons(
    *,
    matches: Sequence[NucleusDecoyMatch],
    real_states: Sequence[PhysicalClosureState],
    decoy_states: Sequence[PhysicalClosureState],
) -> tuple[PhysicalStateDecoyComparison, ...]:
    state_by_event = {
        state.event_id: state for state in tuple(real_states) + tuple(decoy_states)
    }
    output: list[PhysicalStateDecoyComparison] = []
    for match in matches:
        real = state_by_event[match.real_event_id]
        decoy = state_by_event[match.decoy_event_id]
        output.append(
            PhysicalStateDecoyComparison(
                row_id=match.row_id,
                source_accession=match.source_accession,
                real_state_id=real.state_id,
                decoy_state_id=decoy.state_id,
                real_event_id=real.event_id,
                decoy_event_id=decoy.event_id,
                real_physical_state_score=real.physical_state_score,
                decoy_physical_state_score=decoy.physical_state_score,
                real_beats_decoy_by_score=(
                    real.physical_state_score > decoy.physical_state_score
                ),
                real_native_positive_after_scoring=(
                    real.native_contact_count_after_scoring > 0
                ),
                decoy_native_positive_after_scoring=(
                    decoy.native_contact_count_after_scoring > 0
                ),
            )
        )
    return tuple(output)


def mean_physical_score(states: Sequence[PhysicalClosureState]) -> float:
    if not states:
        return 0.0
    return round(sum(state.physical_state_score for state in states) / len(states), 6)


def real_vs_decoy_physical_enrichment_ratio(
    comparisons: Sequence[PhysicalStateDecoyComparison],
) -> float:
    if not comparisons:
        return 0.0
    real_mean = sum(item.real_physical_state_score for item in comparisons) / len(
        comparisons
    )
    decoy_mean = sum(item.decoy_physical_state_score for item in comparisons) / len(
        comparisons
    )
    if decoy_mean == 0.0:
        return 0.0
    return round(real_mean / decoy_mean, 6)


def real_beats_decoy_score_rate(
    comparisons: Sequence[PhysicalStateDecoyComparison],
) -> float:
    if not comparisons:
        return 0.0
    return round(
        sum(1 for item in comparisons if item.real_beats_decoy_by_score)
        / len(comparisons),
        6,
    )

