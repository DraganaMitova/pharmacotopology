from __future__ import annotations

from dataclasses import asdict, dataclass
from statistics import median
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_contact_topology import (
    BETA_RESIDUES,
    ContactTopologyPrediction,
    VisualMechanismRow,
)
from pharmacotopology.folding_energy_landscape import EnergyLandscapePacket
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    normalized_contact_pairs,
)


CONTACT_GAP_ANALYSIS_SIGNATURE_KIND = "native_gap_analysis_after_prediction_v1"
LONG_RANGE_THRESHOLD = 24
SHORT_RANGE_THRESHOLD = 12


@dataclass(frozen=True)
class ContactCluster:
    cluster_id: str
    cluster_kind: str
    contact_count: int
    separation_class: str
    center_i: float
    center_j: float
    min_i: int
    max_i: int
    min_j: int
    max_j: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ContactGapPacket:
    row_id: str
    gap_analysis_signature_kind: str
    missed_native_contact_count: int
    false_predicted_contact_count: int
    native_cluster_miss_count: int
    false_contact_cluster_count: int
    missed_native_contact_clusters: tuple[ContactCluster, ...]
    false_predicted_contact_clusters: tuple[ContactCluster, ...]
    beta_pairing_contact_recall: float | None
    beta_pairing_native_contact_count: int
    beta_pairing_true_positive_count: int
    first_native_contact_rank: int | None
    median_native_contact_rank: float | None
    first_native_contact_step: float | None
    closure_timing_label: str
    premature_compaction: bool
    failure_mechanism: str

    def to_safe_dict(self) -> dict[str, object]:
        return {
            "row_id": self.row_id,
            "gap_analysis_signature_kind": self.gap_analysis_signature_kind,
            "missed_native_contact_count": self.missed_native_contact_count,
            "false_predicted_contact_count": self.false_predicted_contact_count,
            "native_cluster_miss_count": self.native_cluster_miss_count,
            "false_contact_cluster_count": self.false_contact_cluster_count,
            "missed_native_contact_clusters": summarize_clusters(
                self.missed_native_contact_clusters
            ),
            "false_predicted_contact_clusters": summarize_clusters(
                self.false_predicted_contact_clusters
            ),
            "beta_pairing_contact_recall": (
                ""
                if self.beta_pairing_contact_recall is None
                else self.beta_pairing_contact_recall
            ),
            "beta_pairing_native_contact_count": (
                self.beta_pairing_native_contact_count
            ),
            "beta_pairing_true_positive_count": (
                self.beta_pairing_true_positive_count
            ),
            "first_native_contact_rank": (
                "" if self.first_native_contact_rank is None else self.first_native_contact_rank
            ),
            "median_native_contact_rank": (
                "" if self.median_native_contact_rank is None else self.median_native_contact_rank
            ),
            "first_native_contact_step": (
                "" if self.first_native_contact_step is None else self.first_native_contact_step
            ),
            "closure_timing_label": self.closure_timing_label,
            "premature_compaction": self.premature_compaction,
            "failure_mechanism": self.failure_mechanism,
        }


def analyze_contact_gap(
    *,
    row: VisualMechanismRow,
    prediction: ContactTopologyPrediction,
    metrics: ContactMetricPacket,
    energy: EnergyLandscapePacket,
    visible_partial_success: bool,
) -> ContactGapPacket:
    native = set(normalized_contact_pairs(row.native_contact_pairs))
    predicted = set(normalized_contact_pairs(prediction.predicted_contact_pairs))
    missed_native = sorted(native - predicted)
    false_predicted = sorted(predicted - native)
    missed_clusters = cluster_contact_pairs(missed_native, cluster_kind="missed_native")
    false_clusters = cluster_contact_pairs(false_predicted, cluster_kind="false_predicted")
    beta_recall, beta_native_count, beta_tp_count = beta_pairing_contact_recall(
        row=row,
        predicted_pairs=predicted,
    )
    timing = closure_timing_analysis(prediction=prediction, native_pairs=native)
    premature = _premature_compaction(
        metrics=metrics,
        energy=energy,
        false_cluster_count=len(false_clusters),
    )
    failure_mechanism = classify_failure_mechanism(
        row=row,
        metrics=metrics,
        beta_pairing_recall=beta_recall,
        premature_compaction=premature,
        visible_partial_success=visible_partial_success,
    )
    return ContactGapPacket(
        row_id=row.row_id,
        gap_analysis_signature_kind=CONTACT_GAP_ANALYSIS_SIGNATURE_KIND,
        missed_native_contact_count=len(missed_native),
        false_predicted_contact_count=len(false_predicted),
        native_cluster_miss_count=len(missed_clusters),
        false_contact_cluster_count=len(false_clusters),
        missed_native_contact_clusters=missed_clusters,
        false_predicted_contact_clusters=false_clusters,
        beta_pairing_contact_recall=beta_recall,
        beta_pairing_native_contact_count=beta_native_count,
        beta_pairing_true_positive_count=beta_tp_count,
        first_native_contact_rank=timing["first_native_contact_rank"],
        median_native_contact_rank=timing["median_native_contact_rank"],
        first_native_contact_step=timing["first_native_contact_step"],
        closure_timing_label=str(timing["closure_timing_label"]),
        premature_compaction=premature,
        failure_mechanism=failure_mechanism,
    )


def cluster_contact_pairs(
    pairs: Iterable[ContactPair],
    *,
    cluster_kind: str,
    proximity: int = 4,
) -> tuple[ContactCluster, ...]:
    remaining = list(normalized_contact_pairs(pairs))
    clusters: list[list[ContactPair]] = []
    while remaining:
        seed = remaining.pop(0)
        cluster = [seed]
        changed = True
        while changed:
            changed = False
            for pair in tuple(remaining):
                if _near_cluster(pair, cluster, proximity=proximity):
                    cluster.append(pair)
                    remaining.remove(pair)
                    changed = True
        clusters.append(sorted(cluster))
    return tuple(
        _cluster_packet(
            cluster,
            cluster_id=f"{cluster_kind}_{index:02d}",
            cluster_kind=cluster_kind,
        )
        for index, cluster in enumerate(clusters, start=1)
    )


def summarize_clusters(clusters: Sequence[ContactCluster]) -> str:
    return ";".join(
        (
            f"{cluster.cluster_id}:n={cluster.contact_count}:"
            f"{cluster.separation_class}:"
            f"{cluster.center_i:.1f}-{cluster.center_j:.1f}"
        )
        for cluster in clusters
    )


def beta_pairing_contact_recall(
    *,
    row: VisualMechanismRow,
    predicted_pairs: Iterable[ContactPair],
) -> tuple[float | None, int, int]:
    axis = row.truth_axes.get("secondary_structure_axis", "weak_or_unknown")
    if axis not in {"beta_rich", "alpha_beta_mixed"}:
        return None, 0, 0
    native_targets = tuple(
        pair
        for pair in normalized_contact_pairs(row.native_contact_pairs)
        if _is_beta_pairing_target(row.sequence, pair)
    )
    if not native_targets:
        return None, 0, 0
    predicted = set(normalized_contact_pairs(predicted_pairs))
    true_positive = sum(1 for pair in native_targets if pair in predicted)
    return round(true_positive / len(native_targets), 6), len(native_targets), true_positive


def closure_timing_analysis(
    *,
    prediction: ContactTopologyPrediction,
    native_pairs: Iterable[ContactPair],
) -> dict[str, object]:
    native = set(normalized_contact_pairs(native_pairs))
    ranks = [
        rank
        for rank, candidate in enumerate(prediction.candidates, start=1)
        if candidate.pair() in native
    ]
    if not ranks:
        return {
            "first_native_contact_rank": None,
            "median_native_contact_rank": None,
            "first_native_contact_step": None,
            "closure_timing_label": "no_native_contact_overlap",
        }
    candidate_count = max(len(prediction.candidates), 1)
    first_rank = min(ranks)
    first_step = round(first_rank / candidate_count, 6)
    median_rank = round(float(median(ranks)), 6)
    if first_step <= 0.25:
        label = "early_native_closure"
    elif first_step <= 0.60:
        label = "mid_native_closure"
    else:
        label = "late_native_closure"
    return {
        "first_native_contact_rank": first_rank,
        "median_native_contact_rank": median_rank,
        "first_native_contact_step": first_step,
        "closure_timing_label": label,
    }


def classify_failure_mechanism(
    *,
    row: VisualMechanismRow,
    metrics: ContactMetricPacket,
    beta_pairing_recall: float | None,
    premature_compaction: bool,
    visible_partial_success: bool,
) -> str:
    axes = row.truth_axes
    if visible_partial_success:
        return "visible_partial_success"
    if axes.get("order_axis") == "disordered_flexible":
        return "disorder_over_collapse"
    if axes.get("environment_axis") == "membrane_like":
        return "membrane_mis_topology"
    if axes.get("architecture_axis") == "multidomain_or_segmented":
        return "architecture_fragmentation"
    if beta_pairing_recall is not None and beta_pairing_recall < 0.35:
        return "bad_beta_pairing"
    if premature_compaction:
        return "premature_compaction"
    if metrics.short_range_contact_recall > metrics.long_range_contact_recall:
        return "local_only_collapse"
    return "native_gap_unresolved"


def _near_cluster(
    pair: ContactPair,
    cluster: Sequence[ContactPair],
    *,
    proximity: int,
) -> bool:
    return any(
        abs(pair[0] - other[0]) <= proximity
        and abs(pair[1] - other[1]) <= proximity
        for other in cluster
    )


def _cluster_packet(
    cluster: Sequence[ContactPair],
    *,
    cluster_id: str,
    cluster_kind: str,
) -> ContactCluster:
    lefts = [pair[0] for pair in cluster]
    rights = [pair[1] for pair in cluster]
    separations = [pair[1] - pair[0] for pair in cluster]
    return ContactCluster(
        cluster_id=cluster_id,
        cluster_kind=cluster_kind,
        contact_count=len(cluster),
        separation_class=_dominant_separation_class(separations),
        center_i=round(sum(lefts) / len(lefts), 3),
        center_j=round(sum(rights) / len(rights), 3),
        min_i=min(lefts),
        max_i=max(lefts),
        min_j=min(rights),
        max_j=max(rights),
    )


def _dominant_separation_class(separations: Sequence[int]) -> str:
    labels = [_separation_class(separation) for separation in separations]
    return max(sorted(set(labels)), key=labels.count)


def _separation_class(separation: int) -> str:
    if separation <= SHORT_RANGE_THRESHOLD:
        return "short_range"
    if separation < LONG_RANGE_THRESHOLD:
        return "medium_range"
    return "long_range"


def _is_beta_pairing_target(sequence: str, pair: ContactPair) -> bool:
    separation = pair[1] - pair[0]
    if separation < SHORT_RANGE_THRESHOLD:
        return False
    left = sequence[pair[0] - 1]
    right = sequence[pair[1] - 1]
    if left in BETA_RESIDUES or right in BETA_RESIDUES:
        return True
    return separation >= LONG_RANGE_THRESHOLD


def _premature_compaction(
    *,
    metrics: ContactMetricPacket,
    energy: EnergyLandscapePacket,
    false_cluster_count: int,
) -> bool:
    return (
        metrics.false_contact_rate >= 0.85
        and metrics.native_contact_recall < 0.25
        and false_cluster_count >= 3
        and energy.trajectory_contact_gain_monotonicity >= 0.80
    )


def gap_rows_from_packets(
    packets: Sequence[ContactGapPacket],
) -> list[dict[str, object]]:
    return [packet.to_safe_dict() for packet in packets]


def gap_packet_by_row_id(
    packets: Sequence[ContactGapPacket],
) -> Mapping[str, ContactGapPacket]:
    return {packet.row_id: packet for packet in packets}
