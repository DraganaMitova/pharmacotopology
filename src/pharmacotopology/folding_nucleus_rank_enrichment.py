from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Sequence

from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent


NUCLEUS_RANK_ENRICHMENT_KIND = "nucleus_graph_rank_enrichment_v1"
RANK_ENRICHMENT_CUTOFFS = (10, 25, 50)


@dataclass(frozen=True)
class RankEnrichmentRow:
    row_id: str
    source_accession: str
    cutoff: int
    top_rank_event_count: int
    top_rank_native_positive_rate: float
    baseline_native_positive_rate: float
    rank_enrichment: float
    native_truth_used_before_ranking: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(value, 6)


def native_positive_rate(events: Sequence[NucleusClosureEvent]) -> float:
    if not events:
        return 0.0
    return _rounded(
        sum(1 for event in events if event.native_contact_count_after_scoring > 0)
        / len(events)
    )


def rank_enrichment_rows_for_row(
    *,
    row_id: str,
    source_accession: str,
    candidate_events: Sequence[NucleusClosureEvent],
    score_function: Callable[[NucleusClosureEvent], float],
    cutoffs: Sequence[int] = RANK_ENRICHMENT_CUTOFFS,
) -> tuple[RankEnrichmentRow, ...]:
    ranked = sorted(
        candidate_events,
        key=lambda event: (
            -score_function(event),
            event.segment_a_start,
            event.segment_b_start,
            event.event_id,
        ),
    )
    baseline_rate = native_positive_rate(ranked)
    rows: list[RankEnrichmentRow] = []
    for cutoff in cutoffs:
        top = tuple(ranked[: min(cutoff, len(ranked))])
        top_rate = native_positive_rate(top)
        enrichment = top_rate / baseline_rate if baseline_rate else 0.0
        rows.append(
            RankEnrichmentRow(
                row_id=row_id,
                source_accession=source_accession,
                cutoff=cutoff,
                top_rank_event_count=len(top),
                top_rank_native_positive_rate=top_rate,
                baseline_native_positive_rate=baseline_rate,
                rank_enrichment=_rounded(enrichment),
            )
        )
    return tuple(rows)


def mean_rank_enrichment_at(
    rows: Sequence[RankEnrichmentRow],
    cutoff: int,
) -> float:
    values = [row.rank_enrichment for row in rows if row.cutoff == cutoff]
    if not values:
        return 0.0
    return _rounded(sum(values) / len(values))


def mean_native_positive_top_rank_rate(
    rows: Sequence[RankEnrichmentRow],
    cutoff: int,
) -> float:
    values = [
        row.top_rank_native_positive_rate
        for row in rows
        if row.cutoff == cutoff
    ]
    if not values:
        return 0.0
    return _rounded(sum(values) / len(values))
