from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

from pharmacotopology.folding_evolutionary_constraints import (
    CouplingClosureAssessment,
)
from pharmacotopology.folding_nucleus_decoy_falsification import NucleusDecoyMatch


COUPLING_DECOY_FALSIFICATION_KIND = "coupling_supported_decoy_falsification_v1"


@dataclass(frozen=True)
class CouplingDecoyComparison:
    row_id: str
    source_accession: str
    real_event_id: str
    decoy_event_id: str
    real_coupling_selectivity_score: float
    decoy_coupling_selectivity_score: float
    real_direct_support_score: float
    decoy_direct_support_score: float
    real_future_preservation_score: float
    decoy_future_preservation_score: float
    real_blocked_future_pressure: float
    decoy_blocked_future_pressure: float
    real_beats_decoy_by_coupling_score: bool
    real_native_positive_after_scoring: bool
    decoy_native_positive_after_scoring: bool
    coordinate_truth_used_to_build_constraints: bool
    native_truth_used_before_coupling_selection: bool
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def coupling_decoy_comparisons(
    *,
    matches: Sequence[NucleusDecoyMatch],
    assessments: Sequence[CouplingClosureAssessment],
) -> tuple[CouplingDecoyComparison, ...]:
    assessment_by_event = {assessment.event_id: assessment for assessment in assessments}
    output: list[CouplingDecoyComparison] = []
    for match in matches:
        real = assessment_by_event[match.real_event_id]
        decoy = assessment_by_event[match.decoy_event_id]
        output.append(
            CouplingDecoyComparison(
                row_id=match.row_id,
                source_accession=match.source_accession,
                real_event_id=match.real_event_id,
                decoy_event_id=match.decoy_event_id,
                real_coupling_selectivity_score=real.coupling_selectivity_score,
                decoy_coupling_selectivity_score=decoy.coupling_selectivity_score,
                real_direct_support_score=real.direct_support_score,
                decoy_direct_support_score=decoy.direct_support_score,
                real_future_preservation_score=real.future_preservation_score,
                decoy_future_preservation_score=decoy.future_preservation_score,
                real_blocked_future_pressure=real.blocked_future_pressure,
                decoy_blocked_future_pressure=decoy.blocked_future_pressure,
                real_beats_decoy_by_coupling_score=(
                    real.coupling_selectivity_score
                    > decoy.coupling_selectivity_score
                ),
                real_native_positive_after_scoring=(
                    match.real_native_positive_after_scoring
                ),
                decoy_native_positive_after_scoring=(
                    match.decoy_native_positive_after_scoring
                ),
                coordinate_truth_used_to_build_constraints=(
                    real.coordinate_truth_used_to_build_constraints
                    or decoy.coordinate_truth_used_to_build_constraints
                ),
                native_truth_used_before_coupling_selection=(
                    real.native_truth_used_before_coupling_selection
                    or decoy.native_truth_used_before_coupling_selection
                ),
            )
        )
    return tuple(output)


def real_vs_decoy_coupling_enrichment_ratio(
    comparisons: Sequence[CouplingDecoyComparison],
) -> float:
    if not comparisons:
        return 0.0
    real_mean = sum(
        item.real_coupling_selectivity_score for item in comparisons
    ) / len(comparisons)
    decoy_mean = sum(
        item.decoy_coupling_selectivity_score for item in comparisons
    ) / len(comparisons)
    if decoy_mean == 0.0:
        return 0.0
    return round(real_mean / decoy_mean, 6)


def real_beats_decoy_coupling_score_rate(
    comparisons: Sequence[CouplingDecoyComparison],
) -> float:
    if not comparisons:
        return 0.0
    return round(
        sum(1 for item in comparisons if item.real_beats_decoy_by_coupling_score)
        / len(comparisons),
        6,
    )
