from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_structure_benchmark import (
    LABEL_BENCHMARK_KIND_V0,
    ORDER_CONTROL_KINDS,
    STRUCTURE_BENCHMARK_KIND,
    StructureEvidenceRow,
    _control_sequence,
    _signature_delta_mean,
    load_real10_with_structure_evidence,
)
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    CHARGED_AMINO_ACIDS,
    CONTACT_PROXY_DIMENSIONS,
    DISORDER_PROMOTING_AMINO_ACIDS,
    FOLDING_TOPOLOGY_DIMENSIONS,
    HYDROPHOBIC_AMINO_ACIDS,
    FoldingReferenceExample,
    FoldingTopologySignature,
    classify_broad_fold_topology,
    contact_map_proxy_similarity,
    normalize_fold_class,
    normalize_sequence,
    sequence_features,
    signature_to_dict,
)


ORDER_AWARE_BENCHMARK_KIND = "order_aware_folding_topology_benchmark"
POSITIVE_AMINO_ACIDS = frozenset("KRH")
NEGATIVE_AMINO_ACIDS = frozenset("DE")
BREAKERS = frozenset("PG")
HELIX_PRESSURE_AMINO_ACIDS = frozenset("ALEMQKRH")
BETA_PRESSURE_AMINO_ACIDS = frozenset("VIFYWTC")
ORDER_FEATURE_KEYS: tuple[str, ...] = (
    "hydrophobic_cluster_density",
    "hydrophobic_cluster_max_fraction",
    "hydrophobic_cluster_periodicity",
    "hydrophobic_interruption_rate",
    "charge_blockiness",
    "opposite_charge_pairing_potential",
    "same_charge_repulsion_pressure",
    "breaker_boundary_alignment",
    "breaker_hydrophobic_interruption",
    "cysteine_bridge_potential",
    "local_helix_pressure_mean",
    "local_beta_pressure_mean",
    "local_disorder_pressure_mean",
    "segment_boundary_contrast",
    "long_range_closure_potential",
    "contact_prior_density",
    "contact_prior_mean_weight",
    "contact_prior_max_weight",
    "contact_prior_average_order",
    "contact_prior_crossing_fraction",
)


@dataclass(frozen=True)
class Cluster:
    start: int
    end: int
    length: int
    kind: str


@dataclass(frozen=True)
class SegmentProfile:
    index: int
    start: int
    end: int
    length: int
    hydrophobic_fraction: float
    aromatic_fraction: float
    positive_fraction: float
    negative_fraction: float
    charged_fraction: float
    breaker_fraction: float
    cysteine_fraction: float
    disorder_fraction: float
    helix_pressure: float
    beta_pressure: float
    collapse_pressure: float
    charge_frustration: float
    boundary_pressure: float


@dataclass(frozen=True)
class ContactPriorEdge:
    protein_id: str
    segment_i: int
    segment_j: int
    segment_i_start: int
    segment_i_end: int
    segment_j_start: int
    segment_j_end: int
    sequence_distance: int
    sequence_distance_norm: float
    edge_weight: float
    hydrophobic_compatibility: float
    opposite_charge_compatibility: float
    same_charge_repulsion: float
    cysteine_bridge_potential: float
    aromatic_anchor_potential: float
    breaker_penalty: float
    disorder_penalty: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class OrderAwarePrediction:
    protein_id: str
    sequence_length: int
    topology_signature: FoldingTopologySignature
    fold_class: str
    order_features: dict[str, float]
    contact_prior_edges: tuple[ContactPriorEdge, ...]
    order_topology_score: float
    contact_prior_signal: float


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _mean(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _fraction(sequence: str, alphabet: frozenset[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in alphabet) / len(sequence)


def _clusters(sequence: str, alphabet: frozenset[str], kind: str) -> tuple[Cluster, ...]:
    clusters: list[Cluster] = []
    start: int | None = None
    for index, residue in enumerate(sequence):
        if residue in alphabet and start is None:
            start = index
        elif residue not in alphabet and start is not None:
            clusters.append(Cluster(start=start, end=index, length=index - start, kind=kind))
            start = None
    if start is not None:
        clusters.append(
            Cluster(start=start, end=len(sequence), length=len(sequence) - start, kind=kind)
        )
    return tuple(clusters)


def _profile_for_segment(index: int, start: int, segment: str) -> SegmentProfile:
    length = max(len(segment), 1)
    hydrophobic = _fraction(segment, HYDROPHOBIC_AMINO_ACIDS)
    aromatic = _fraction(segment, AROMATIC_AMINO_ACIDS)
    positive = _fraction(segment, POSITIVE_AMINO_ACIDS)
    negative = _fraction(segment, NEGATIVE_AMINO_ACIDS)
    charged = _fraction(segment, CHARGED_AMINO_ACIDS)
    breaker = _fraction(segment, BREAKERS)
    cysteine = _fraction(segment, frozenset("C"))
    disorder = _fraction(segment, DISORDER_PROMOTING_AMINO_ACIDS)
    helix = _clamp(
        _fraction(segment, HELIX_PRESSURE_AMINO_ACIDS) * 0.85
        + hydrophobic * 0.15
        - breaker * 0.55
    )
    beta = _clamp(
        _fraction(segment, BETA_PRESSURE_AMINO_ACIDS) * 0.90
        + aromatic * 0.20
        - breaker * 0.28
    )
    collapse = _clamp(hydrophobic * 0.70 + aromatic * 0.20 + cysteine * 0.10)
    charge_frustration = _clamp(
        abs(positive - negative) * 0.85 + charged * 0.20
    )
    boundary = _clamp(breaker * 0.48 + disorder * 0.32 + charge_frustration * 0.20)
    return SegmentProfile(
        index=index,
        start=start,
        end=start + len(segment),
        length=length,
        hydrophobic_fraction=round(hydrophobic, 6),
        aromatic_fraction=round(aromatic, 6),
        positive_fraction=round(positive, 6),
        negative_fraction=round(negative, 6),
        charged_fraction=round(charged, 6),
        breaker_fraction=round(breaker, 6),
        cysteine_fraction=round(cysteine, 6),
        disorder_fraction=round(disorder, 6),
        helix_pressure=round(helix, 6),
        beta_pressure=round(beta, 6),
        collapse_pressure=round(collapse, 6),
        charge_frustration=round(charge_frustration, 6),
        boundary_pressure=round(boundary, 6),
    )


def _segments(sequence: str, segment_size: int = 12) -> tuple[SegmentProfile, ...]:
    return tuple(
        _profile_for_segment(
            index=index,
            start=start,
            segment=sequence[start : start + segment_size],
        )
        for index, start in enumerate(range(0, len(sequence), segment_size))
    )


def _adjacent_contrast(segments: Sequence[SegmentProfile]) -> float:
    if len(segments) < 2:
        return 0.0
    contrasts = []
    for left, right in zip(segments, segments[1:]):
        contrasts.append(
            (
                abs(left.hydrophobic_fraction - right.hydrophobic_fraction)
                + abs(left.positive_fraction - right.positive_fraction)
                + abs(left.negative_fraction - right.negative_fraction)
                + abs(left.breaker_fraction - right.breaker_fraction)
                + abs(left.helix_pressure - right.helix_pressure)
                + abs(left.beta_pressure - right.beta_pressure)
                + abs(left.disorder_fraction - right.disorder_fraction)
            )
            / 7
        )
    return round(_mean(contrasts), 6)


def _cluster_periodicity(clusters: Sequence[Cluster], sequence_length: int) -> float:
    if len(clusters) < 3 or sequence_length <= 1:
        return 0.0
    centers = [(cluster.start + cluster.end) / 2 for cluster in clusters]
    gaps = [right - left for left, right in zip(centers, centers[1:])]
    mean_gap = _mean(gaps)
    if mean_gap == 0:
        return 0.0
    mean_abs_deviation = _mean(abs(gap - mean_gap) for gap in gaps)
    return _rounded(1.0 - mean_abs_deviation / max(mean_gap, 1.0))


def _breaker_interruption_rate(sequence: str, hydrophobic_clusters: Sequence[Cluster]) -> float:
    if not hydrophobic_clusters:
        return 0.0
    interrupted = 0
    for left, right in zip(hydrophobic_clusters, hydrophobic_clusters[1:]):
        gap = sequence[left.end : right.start]
        if any(residue in BREAKERS for residue in gap):
            interrupted += 1
    return round(interrupted / max(len(hydrophobic_clusters) - 1, 1), 6)


def _charge_blockiness(sequence: str) -> float:
    charge_symbols = [
        "+"
        if residue in POSITIVE_AMINO_ACIDS
        else "-"
        if residue in NEGATIVE_AMINO_ACIDS
        else "0"
        for residue in sequence
    ]
    charged = [symbol for symbol in charge_symbols if symbol != "0"]
    if len(charged) < 2:
        return 0.0
    same_neighbors = sum(
        1 for left, right in zip(charged, charged[1:]) if left == right
    )
    return round(same_neighbors / (len(charged) - 1), 6)


def _cysteine_bridge_potential(sequence: str) -> float:
    positions = [index for index, residue in enumerate(sequence) if residue == "C"]
    if len(positions) < 2:
        return 0.0
    distances = [
        abs(right - left) / max(len(sequence), 1)
        for left_index, left in enumerate(positions)
        for right in positions[left_index + 1 :]
    ]
    return _rounded(_mean(min(distance * 1.8, 1.0) for distance in distances))


def _contact_prior_edges(
    protein_id: str,
    segments: Sequence[SegmentProfile],
    sequence_length: int,
) -> tuple[ContactPriorEdge, ...]:
    edges: list[ContactPriorEdge] = []
    for left_index, left in enumerate(segments):
        for right in segments[left_index + 2 :]:
            center_left = (left.start + left.end) / 2
            center_right = (right.start + right.end) / 2
            distance = int(abs(center_right - center_left))
            distance_norm = distance / max(sequence_length, 1)
            distance_factor = _clamp(distance_norm * 1.8)
            hydrophobic = min(left.collapse_pressure, right.collapse_pressure)
            opposite_charge = (
                left.positive_fraction * right.negative_fraction
                + left.negative_fraction * right.positive_fraction
            )
            same_charge = (
                left.positive_fraction * right.positive_fraction
                + left.negative_fraction * right.negative_fraction
            )
            cysteine = min(left.cysteine_fraction, right.cysteine_fraction) * distance_factor
            aromatic = min(left.aromatic_fraction, right.aromatic_fraction)
            breaker_penalty = (left.breaker_fraction + right.breaker_fraction) / 2
            disorder_penalty = (left.disorder_fraction + right.disorder_fraction) / 2
            raw_weight = (
                hydrophobic * 0.42
                + opposite_charge * 1.15
                + cysteine * 0.34
                + aromatic * 0.28
                + distance_factor * 0.16
                - same_charge * 0.40
                - breaker_penalty * 0.24
                - disorder_penalty * 0.18
            )
            weight = _rounded(raw_weight)
            if weight <= 0.05:
                continue
            edges.append(
                ContactPriorEdge(
                    protein_id=protein_id,
                    segment_i=left.index,
                    segment_j=right.index,
                    segment_i_start=left.start + 1,
                    segment_i_end=left.end,
                    segment_j_start=right.start + 1,
                    segment_j_end=right.end,
                    sequence_distance=distance,
                    sequence_distance_norm=round(distance_norm, 6),
                    edge_weight=weight,
                    hydrophobic_compatibility=round(hydrophobic, 6),
                    opposite_charge_compatibility=round(opposite_charge, 6),
                    same_charge_repulsion=round(same_charge, 6),
                    cysteine_bridge_potential=round(cysteine, 6),
                    aromatic_anchor_potential=round(aromatic, 6),
                    breaker_penalty=round(breaker_penalty, 6),
                    disorder_penalty=round(disorder_penalty, 6),
                )
            )
    return tuple(sorted(edges, key=lambda edge: edge.edge_weight, reverse=True))


def _edge_crossing_fraction(edges: Sequence[ContactPriorEdge]) -> float:
    if len(edges) < 2:
        return 0.0
    crossings = 0
    reviewed = 0
    for index, left in enumerate(edges):
        for right in edges[index + 1 :]:
            reviewed += 1
            if left.segment_i < right.segment_i < left.segment_j < right.segment_j:
                crossings += 1
    return round(crossings / reviewed, 6) if reviewed else 0.0


def extract_order_aware_features(
    sequence: str,
    *,
    protein_id: str = "sequence",
) -> tuple[dict[str, float], tuple[ContactPriorEdge, ...]]:
    normalized = normalize_sequence(sequence)
    segments = _segments(normalized)
    hydrophobic_clusters = _clusters(
        normalized,
        HYDROPHOBIC_AMINO_ACIDS,
        "hydrophobic",
    )
    breaker_clusters = _clusters(normalized, BREAKERS, "breaker")
    edges = _contact_prior_edges(protein_id, segments, len(normalized))
    possible_edges = max(len(segments) * (len(segments) - 1) / 2, 1.0)
    top_edges = edges[: max(1, min(12, len(edges)))]
    mean_weight = _mean(edge.edge_weight for edge in edges)
    mean_top_weight = _mean(edge.edge_weight for edge in top_edges)
    long_range_weight = _mean(
        edge.edge_weight
        for edge in edges
        if edge.sequence_distance_norm >= 0.30
    )
    opposite_pairing = _mean(edge.opposite_charge_compatibility for edge in edges)
    same_repulsion = _mean(edge.same_charge_repulsion for edge in edges)
    hydrophobic_edge = _mean(edge.hydrophobic_compatibility for edge in edges)
    contact_order = _mean(
        edge.sequence_distance_norm * edge.edge_weight for edge in edges
    )
    breaker_boundary_alignment = _mean(
        min(cluster.length / 4, 1.0)
        for cluster in breaker_clusters
    )
    features = {
        "hydrophobic_cluster_count": float(len(hydrophobic_clusters)),
        "hydrophobic_cluster_density": round(
            len(hydrophobic_clusters) / max(len(normalized) / 18, 1.0),
            6,
        ),
        "hydrophobic_cluster_max_fraction": round(
            max((cluster.length for cluster in hydrophobic_clusters), default=0)
            / max(len(normalized), 1),
            6,
        ),
        "hydrophobic_cluster_periodicity": _cluster_periodicity(
            hydrophobic_clusters,
            len(normalized),
        ),
        "hydrophobic_interruption_rate": _breaker_interruption_rate(
            normalized,
            hydrophobic_clusters,
        ),
        "charge_blockiness": _charge_blockiness(normalized),
        "opposite_charge_pairing_potential": round(opposite_pairing, 6),
        "same_charge_repulsion_pressure": round(same_repulsion, 6),
        "breaker_density": _fraction(normalized, BREAKERS),
        "breaker_boundary_alignment": round(breaker_boundary_alignment, 6),
        "breaker_hydrophobic_interruption": _breaker_interruption_rate(
            normalized,
            hydrophobic_clusters,
        ),
        "cysteine_bridge_potential": _cysteine_bridge_potential(normalized),
        "local_helix_pressure_mean": round(
            _mean(segment.helix_pressure for segment in segments),
            6,
        ),
        "local_beta_pressure_mean": round(
            _mean(segment.beta_pressure for segment in segments),
            6,
        ),
        "local_disorder_pressure_mean": round(
            _mean(segment.disorder_fraction for segment in segments),
            6,
        ),
        "local_collapse_pressure_mean": round(
            _mean(segment.collapse_pressure for segment in segments),
            6,
        ),
        "segment_boundary_contrast": _adjacent_contrast(segments),
        "long_range_closure_potential": round(long_range_weight, 6),
        "contact_prior_edge_count": float(len(edges)),
        "contact_prior_density": round(len(edges) / possible_edges, 6),
        "contact_prior_mean_weight": round(mean_weight, 6),
        "contact_prior_top_weight": round(mean_top_weight, 6),
        "contact_prior_max_weight": round(
            max((edge.edge_weight for edge in edges), default=0.0),
            6,
        ),
        "contact_prior_average_order": round(contact_order, 6),
        "contact_prior_hydrophobic_weight": round(hydrophobic_edge, 6),
        "contact_prior_crossing_fraction": _edge_crossing_fraction(top_edges),
    }
    return features, edges


def predict_order_aware_topology_signature(
    sequence: str,
    *,
    protein_id: str = "sequence",
) -> OrderAwarePrediction:
    normalized = normalize_sequence(sequence)
    composition = sequence_features(normalized)
    order_features, edges = extract_order_aware_features(
        normalized,
        protein_id=protein_id,
    )
    helix = order_features["local_helix_pressure_mean"]
    beta = order_features["local_beta_pressure_mean"]
    disorder = order_features["local_disorder_pressure_mean"]
    collapse = order_features["local_collapse_pressure_mean"]
    boundary_contrast = order_features["segment_boundary_contrast"]
    contact_density = order_features["contact_prior_density"]
    contact_weight = order_features["contact_prior_mean_weight"]
    contact_max = order_features["contact_prior_max_weight"]
    long_range = order_features["long_range_closure_potential"]
    breaker = order_features["breaker_density"]
    interruption = order_features["hydrophobic_interruption_rate"]
    charge_repulsion = order_features["same_charge_repulsion_pressure"]
    opposite_charge = order_features["opposite_charge_pairing_potential"]
    cysteine = order_features["cysteine_bridge_potential"]
    crossing = order_features["contact_prior_crossing_fraction"]
    contact_signal = _rounded(contact_density * 0.45 + contact_weight * 1.2 + contact_max * 0.35)
    hydrophobic_core = _rounded(
        collapse * 0.48
        + order_features["hydrophobic_cluster_periodicity"] * 0.12
        + order_features["hydrophobic_cluster_max_fraction"] * 1.4
        + contact_weight * 0.64
        + cysteine * 0.10
        - interruption * 0.18
    )
    loop_disorder = _rounded(
        disorder * 0.50
        + breaker * 0.38
        + interruption * 0.18
        + boundary_contrast * 0.16
        + charge_repulsion * 0.12
    )
    secondary_balance = _rounded(
        1.0
        - abs(helix - beta) * 0.72
        - disorder * 0.22
        + boundary_contrast * 0.08
    )
    domain_stability = _rounded(
        0.84
        - boundary_contrast * 0.56
        - breaker * 0.28
        - charge_repulsion * 0.22
        + contact_signal * 0.16
        + opposite_charge * 0.18
    )
    signature = FoldingTopologySignature(
        sequence_complexity=float(composition["sequence_complexity"]),
        secondary_structure_balance=secondary_balance,
        contact_map_closure=_rounded(
            contact_signal * 0.72
            + hydrophobic_core * 0.18
            + long_range * 0.20
            - loop_disorder * 0.12
        ),
        hydrophobic_core_closure=hydrophobic_core,
        loop_disorder_pressure=loop_disorder,
        domain_boundary_stability=domain_stability,
        long_range_contact_order=_rounded(
            order_features["contact_prior_average_order"] * 2.4
            + long_range * 0.42
            + contact_density * 0.14
        ),
        conformational_flexibility=_rounded(
            loop_disorder * 0.58
            + boundary_contrast * 0.24
            + charge_repulsion * 0.18
            + (1.0 - contact_signal) * 0.06
        ),
        knot_or_entanglement_signature=_rounded(
            crossing * 0.50
            + long_range * 0.36
            + contact_density * 0.20
            + order_features["contact_prior_average_order"] * 0.42
        ),
        uncertainty_radius=_rounded(
            0.20
            + loop_disorder * 0.20
            + (1.0 - contact_signal) * 0.12
            + boundary_contrast * 0.16
        ),
    )
    order_score = _rounded(
        signature.contact_map_closure * 0.30
        + signature.hydrophobic_core_closure * 0.24
        + signature.long_range_contact_order * 0.22
        + signature.domain_boundary_stability * 0.12
        + (1.0 - signature.loop_disorder_pressure) * 0.12
    )
    return OrderAwarePrediction(
        protein_id=protein_id,
        sequence_length=len(normalized),
        topology_signature=signature,
        fold_class=classify_broad_fold_topology(signature),
        order_features=order_features,
        contact_prior_edges=edges,
        order_topology_score=order_score,
        contact_prior_signal=contact_signal,
    )


def _feature_delta(left: Mapping[str, float], right: Mapping[str, float]) -> float:
    return round(
        _mean(abs(float(left.get(key, 0.0)) - float(right.get(key, 0.0))) for key in ORDER_FEATURE_KEYS),
        6,
    )


def _edge_shape_distance(
    left_edges: Sequence[ContactPriorEdge],
    right_edges: Sequence[ContactPriorEdge],
) -> float:
    left = {
        (edge.segment_i, edge.segment_j): edge.edge_weight
        for edge in left_edges[:20]
    }
    right = {
        (edge.segment_i, edge.segment_j): edge.edge_weight
        for edge in right_edges[:20]
    }
    keys = sorted(set(left).union(right))
    if not keys:
        return 0.0
    return round(_mean(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in keys), 6)


def order_aware_control_separation_rows(
    references: Sequence[FoldingReferenceExample],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for reference in references:
        real_sequence = normalize_sequence(reference.sequence)
        real_prediction = predict_order_aware_topology_signature(
            real_sequence,
            protein_id=reference.protein_id,
        )
        for control_kind in ORDER_CONTROL_KINDS:
            control_sequence = _control_sequence(
                real_sequence,
                control_kind,
                protein_id=reference.protein_id,
            )
            control_prediction = predict_order_aware_topology_signature(
                control_sequence,
                protein_id=f"{reference.protein_id}:{control_kind}",
            )
            signature_delta = _signature_delta_mean(
                real_prediction.topology_signature,
                control_prediction.topology_signature,
            )
            feature_delta = _feature_delta(
                real_prediction.order_features,
                control_prediction.order_features,
            )
            edge_delta = _edge_shape_distance(
                real_prediction.contact_prior_edges,
                control_prediction.contact_prior_edges,
            )
            score_delta = abs(
                real_prediction.order_topology_score
                - control_prediction.order_topology_score
            )
            separation = _rounded(
                signature_delta * 0.75
                + feature_delta * 1.45
                + edge_delta * 1.15
                + score_delta * 1.00
            )
            rows.append(
                {
                    "protein_id": reference.protein_id,
                    "control_kind": control_kind,
                    "same_composition": sorted(real_sequence) == sorted(control_sequence),
                    "real_sequence_score": real_prediction.order_topology_score,
                    "control_sequence_score": control_prediction.order_topology_score,
                    "score_delta": round(score_delta, 6),
                    "signature_delta_mean": signature_delta,
                    "order_feature_delta_mean": feature_delta,
                    "contact_prior_edge_delta_mean": edge_delta,
                    "separation_score": separation,
                    "real_fold_class": real_prediction.fold_class,
                    "control_fold_class": control_prediction.fold_class,
                    "fold_class_changed": real_prediction.fold_class
                    != control_prediction.fold_class,
                    "control_sequence_written": False,
                }
            )
    return rows


def order_aware_benchmark_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    evidence_by_id = {row.protein_id: row for row in evidence_rows}
    rows: list[dict[str, object]] = []
    for reference in references:
        evidence = evidence_by_id[reference.protein_id]
        prediction = predict_order_aware_topology_signature(
            reference.sequence,
            protein_id=reference.protein_id,
        )
        label_class = normalize_fold_class(reference.reference_fold_class)
        structure_similarity = contact_map_proxy_similarity(
            prediction.topology_signature,
            evidence.structure_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        label_similarity = contact_map_proxy_similarity(
            prediction.topology_signature,
            reference.reference_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        structure_class_match = prediction.fold_class == evidence.structure_fold_class
        label_class_match = prediction.fold_class == label_class
        failure_reason = ""
        if not structure_class_match:
            failure_reason = "order_aware_prediction_structure_class_mismatch"
        elif structure_similarity < 0.70:
            failure_reason = "low_order_aware_prediction_structure_similarity"
        rows.append(
            {
                "protein_id": reference.protein_id,
                "sequence_length": prediction.sequence_length,
                "source_accession": reference.source_accession,
                "reference_structure_source": reference.reference_structure_source,
                "predicted_fold_class": prediction.fold_class,
                "structure_fold_class": evidence.structure_fold_class,
                "label_fold_class": label_class,
                "prediction_vs_structure_score": structure_similarity,
                "prediction_vs_label_score": label_similarity,
                "prediction_structure_class_match": structure_class_match,
                "prediction_label_class_match": label_class_match,
                "order_topology_score": prediction.order_topology_score,
                "contact_prior_signal": prediction.contact_prior_signal,
                "predicted_topology_signature": signature_to_dict(
                    prediction.topology_signature
                ),
                "structure_topology_signature": signature_to_dict(
                    evidence.structure_topology_signature
                ),
                "order_features": prediction.order_features,
                "failure_reason": failure_reason,
                "folding_problem_solved": False,
                "folding_solution_claim_created": False,
            }
        )
    return rows


def contact_prior_rows(
    references: Sequence[FoldingReferenceExample],
    *,
    top_edges_per_protein: int = 15,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for reference in references:
        prediction = predict_order_aware_topology_signature(
            reference.sequence,
            protein_id=reference.protein_id,
        )
        for edge in prediction.contact_prior_edges[:top_edges_per_protein]:
            rows.append(edge.to_dict())
    return rows


def _mean_from_rows(rows: Sequence[Mapping[str, object]], key: str) -> float:
    return round(_mean(float(row.get(key, 0.0)) for row in rows), 6)


def build_order_aware_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = order_aware_benchmark_rows(references, evidence_rows)
    control_rows = order_aware_control_separation_rows(references)
    contact_rows = contact_prior_rows(references)
    structure_matches = sum(
        1 for row in rows if bool(row["prediction_structure_class_match"])
    )
    label_matches = sum(1 for row in rows if bool(row["prediction_label_class_match"]))
    sensitivity = _mean_from_rows(control_rows, "separation_score")
    shuffled_separation = _mean_from_rows(control_rows, "separation_score")
    contact_prior_signal = _mean_from_rows(rows, "contact_prior_signal")
    contact_prior_signal_seen = contact_prior_signal > 0.10 and len(contact_rows) > 0
    recipe_order_blind = sensitivity <= 0.05
    composition_only_warning = sensitivity < 0.25 or not contact_prior_signal_seen
    revision_required = (
        composition_only_warning
        or structure_matches < len(rows)
        or any(row["failure_reason"] for row in rows)
    )
    return {
        "benchmark_kind": ORDER_AWARE_BENCHMARK_KIND,
        "source_structure_benchmark_kind": STRUCTURE_BENCHMARK_KIND,
        "source_label_benchmark_kind": LABEL_BENCHMARK_KIND_V0,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "benchmark_size": len(rows),
        "predictor_input_boundary": (
            "sequence_only_no_labels_no_structure_answers"
        ),
        "sequence_order_sensitivity_score": sensitivity,
        "real_vs_shuffled_separation_mean": shuffled_separation,
        "contact_prior_signal_mean": contact_prior_signal,
        "contact_prior_signal_seen": contact_prior_signal_seen,
        "recipe_order_blind": recipe_order_blind,
        "composition_only_warning": composition_only_warning,
        "prediction_vs_structure_match_count": structure_matches,
        "prediction_vs_structure_mismatch_count": len(rows) - structure_matches,
        "prediction_vs_structure_accuracy": round(structure_matches / max(len(rows), 1), 6),
        "prediction_vs_label_match_count": label_matches,
        "prediction_vs_label_accuracy": round(label_matches / max(len(rows), 1), 6),
        "mean_prediction_vs_structure_score": _mean_from_rows(
            rows,
            "prediction_vs_structure_score",
        ),
        "mean_prediction_vs_label_score": _mean_from_rows(
            rows,
            "prediction_vs_label_score",
        ),
        "control_rows_reviewed": len(control_rows),
        "contact_prior_edges_reported": len(contact_rows),
        "revision_required": revision_required,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "next_scaling_gate": (
            "Do not scale to 50 rows until sequence_order_sensitivity_score "
            "is nonzero and real_vs_shuffled_separation_mean remains meaningful."
        ),
        "boundary_statement": (
            "This order-aware recipe derives local windows, segment boundaries, "
            "hydrophobic clusters, charge patterns, breaker placement, and a "
            "predicted contact-prior graph from sequence only. Labels and "
            "structure evidence are used only after prediction for scoring."
        ),
        "rows": rows,
    }


def write_order_aware_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    contact_rows: Sequence[Mapping[str, object]],
    control_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    contact_prior_path: Path,
    control_separation_path: Path,
    dashboard_path: Path,
) -> tuple[Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(
        rows,
        rows_path,
        nested_json_fields={
            "predicted_topology_signature",
            "structure_topology_signature",
            "order_features",
        },
    )
    _write_csv_rows(contact_rows, contact_prior_path, nested_json_fields=set())
    _write_csv_rows(control_rows, control_separation_path, nested_json_fields=set())
    dashboard_path.write_text(
        render_order_aware_dashboard(report, control_rows),
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        contact_prior_path,
        control_separation_path,
        dashboard_path,
    )


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
    *,
    nested_json_fields: set[str],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            file.write("")
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {}
            for key, value in row.items():
                if key in nested_json_fields:
                    output[key] = json.dumps(value, sort_keys=True, separators=(",", ":"))
                else:
                    output[key] = value
            writer.writerow(output)
    return path


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _control_plot_rows(control_rows: Sequence[Mapping[str, object]]) -> str:
    grouped: dict[str, list[Mapping[str, object]]] = {}
    for row in control_rows:
        grouped.setdefault(str(row["protein_id"]), []).append(row)
    html_rows: list[str] = []
    for protein_id, rows in grouped.items():
        real_score = float(rows[0].get("real_sequence_score", 0.0))
        real_width = round(real_score * 100, 3)
        control_bars = [
            "<div class=\"bar-pair\">"
            f"<span>{_escape(row.get('control_kind', ''))}</span>"
            "<div class=\"bar-track\">"
            f"<div class=\"bar control\" style=\"width:{round(float(row.get('control_sequence_score', 0.0)) * 100, 3)}%\"></div>"
            "</div>"
            f"<strong>{_escape(row.get('control_sequence_score', 0.0))}</strong>"
            "</div>"
            for row in rows
        ]
        html_rows.append(
            "<section class=\"protein\">"
            f"<h3>{_escape(protein_id)}</h3>"
            "<div class=\"bar-pair real\">"
            "<span>real sequence</span>"
            "<div class=\"bar-track\">"
            f"<div class=\"bar real\" style=\"width:{real_width}%\"></div>"
            "</div>"
            f"<strong>{real_score}</strong>"
            "</div>"
            + "".join(control_bars)
            + "</section>"
        )
    return "".join(html_rows)


def render_order_aware_dashboard(
    report: Mapping[str, object],
    control_rows: Sequence[Mapping[str, object]],
) -> str:
    metrics = {
        "Order Sensitivity": report.get("sequence_order_sensitivity_score", 0.0),
        "Real vs Controls": report.get("real_vs_shuffled_separation_mean", 0.0),
        "Contact Prior": report.get("contact_prior_signal_seen", False),
        "Order Blind": report.get("recipe_order_blind", True),
        "Composition Warning": report.get("composition_only_warning", True),
        "Claim Allowed": report.get("claim_allowed", False),
    }
    metric_cards = "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(value)}</strong>"
        "</div>"
        for label, value in metrics.items()
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Order-Aware Folding Topology Benchmark</title>
<style>
body {{ margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17211c; background: #fff; }}
header {{ padding: 28px clamp(18px, 4vw, 54px); background: #eef5f2; border-bottom: 1px solid #cad6d0; }}
main {{ padding: 24px clamp(18px, 4vw, 54px) 52px; }}
h1 {{ margin: 0 0 10px; font-size: clamp(28px, 4vw, 44px); letter-spacing: 0; }}
h2 {{ font-size: 22px; letter-spacing: 0; }}
h3 {{ margin: 18px 0 8px; font-size: 16px; }}
p {{ max-width: 980px; line-height: 1.55; }}
.warning {{ display: inline-block; padding: 8px 10px; border: 1px solid #9b3d2e; color: #9b3d2e; background: #fff8f6; font-weight: 700; }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 10px; margin-top: 18px; }}
.metric {{ border: 1px solid #cad6d0; background: #f7faf8; padding: 12px; min-height: 82px; }}
.metric span {{ display: block; color: #5d6762; font-size: 13px; }}
.metric strong {{ display: block; margin-top: 8px; font-size: 22px; }}
.protein {{ border-top: 1px solid #cad6d0; padding: 12px 0; }}
.bar-pair {{ display: grid; grid-template-columns: minmax(170px, 260px) 1fr 70px; gap: 10px; align-items: center; margin: 8px 0; }}
.bar-track {{ height: 16px; border: 1px solid #cad6d0; background: #fff; }}
.bar {{ height: 100%; background: #59708a; }}
.bar.real {{ background: #0d766e; }}
.bar.control {{ background: #9b6a34; }}
@media (max-width: 760px) {{ .bar-pair {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<header>
  <p class="warning">NOT A FOLDING SOLUTION / ORDER-AWARE HYPOTHESIS TEST</p>
  <h1>Order-Aware Folding Topology Benchmark</h1>
  <p>This page visualizes whether the sequence-only recipe reacts when order is
  broken while composition is preserved. Labels and structure evidence are used
  only after prediction for scoring.</p>
  <div class="metrics">{metric_cards}</div>
</header>
<main>
  <section>
    <h2>Real Sequence vs Shuffled Controls</h2>
    <p>Bars show the order-aware topology score for the real sequence and each
    internal perturbation control. Control sequences are not written.</p>
    {_control_plot_rows(control_rows)}
  </section>
</main>
</body>
</html>
"""


def load_order_aware_inputs(
    benchmark_file: Path,
    structure_evidence_file: Path,
) -> tuple[tuple[FoldingReferenceExample, ...], tuple[StructureEvidenceRow, ...]]:
    return load_real10_with_structure_evidence(
        benchmark_file,
        structure_evidence_file,
    )
