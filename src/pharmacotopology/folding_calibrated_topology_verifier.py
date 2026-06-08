from __future__ import annotations

"""Calibrated fold-topology verifier challenge.

This layer tests the user's next hypothesis: the missing verifier is not generic
compactness/coherence, but a verifier calibrated toward real fold topology.

Two important honesty boundaries are explicit:
* PyRosetta/DFIRE are probed, but optional. In this sandbox they are usually not
  installed/licensed, so the batch falls back to a small transparent structural
  feature verifier.
* The fallback verifier is leave-one-target-out. It may use native contacts from
  the other locked rows to calibrate feature weights, but it never uses the
  target row's native contacts/coordinates before that target's contact map is
  frozen. Because it uses PDB-derived calibration, it cannot claim a universal
  physics law.
"""

import importlib.util
import hashlib
from dataclasses import asdict, dataclass
from math import exp, sqrt
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_coarse_grain_md_geometry import (
    CoarseGrainMDContactDecision,
    _distance,
    _extract_contacts_from_geometry,
    _matched_controls_for_report,
    _pair_hash,
    _relax_global_geometry,
    _rounded,
)
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_geodesic_consensus_filter_md import build_geodesic_consensus_filter_constraints
from pharmacotopology.folding_global_verifier_annealing import (
    GLOBAL_VERIFIER_CLAIM_RULE,
    GLOBAL_VERIFIER_RULE,
    _residue_pair_score,
    collision_penalty,
    compactness_score,
    degree_coherence_score,
    loop_patch_coherence_score,
    restraint_satisfaction_score,
    statistical_potential_surrogate,
)
from pharmacotopology.folding_native_contact_eval import ContactMetricPacket, ContactPair, evaluate_contact_prediction, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import MIN_SEQUENCE_SEPARATION, RealCoordinateVisualRow

CALIBRATED_VERIFIER_KIND = "calibrated_topology_verifier_annealing_v0"
CALIBRATED_VERIFIER_MODE = "external_dca_geodesic_calibrated_topology_verifier"
CALIBRATED_VERIFIER_DECISION_KIND = "calibrated_topology_verifier_contact_decision_v0"
CALIBRATED_VERIFIER_RULE = (
    GLOBAL_VERIFIER_RULE
    + ";pyrosetta_dfire_probe;leave_one_target_out_transparent_structural_feature_discriminator;"
    + "contact_order_degree_patch_restraint_residue_geometry_features;target_native_excluded_from_calibration"
)
CALIBRATED_VERIFIER_CLAIM_RULE = (
    GLOBAL_VERIFIER_CLAIM_RULE
    + ";calibrated_verifier_claim_requires_no_target_native_leakage_and_precision_recall_long_range_ge_0_70;"
    + "pdb_calibration_disallows_universal_physical_law_claim"
)

Vector = list[float]


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def optional_energy_backend_status() -> dict[str, object]:
    return {
        "pyrosetta_available": importlib.util.find_spec("pyrosetta") is not None,
        "dfire_available": importlib.util.find_spec("dfire") is not None,
        "backend_used": "transparent_leave_one_out_feature_discriminator",
        "pyrosetta_used": False,
        "dfire_used": False,
        "reason": "PyRosetta/DFIRE are optional licensed/external backends; fallback is deterministic and dependency-light.",
    }


def _restraint_confidence(restraints: Sequence[CouplingConstraint]) -> dict[ContactPair, float]:
    return {r.pair(): max(0.05, min(1.0, float(r.confidence))) for r in restraints}


def _contact_order(contacts: Sequence[ContactPair], sequence_length: int) -> float:
    if not contacts or sequence_length <= 0:
        return 0.0
    return _rounded(sum(j - i for i, j in contacts) / (len(contacts) * sequence_length))


def _feature_vector(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    pair: ContactPair,
    raw_scores: Mapping[ContactPair, float],
    raw_set: set[ContactPair],
    restraints: Mapping[ContactPair, float],
    degrees: Sequence[int],
) -> tuple[float, ...]:
    i, j = pair
    seq_sep = j - i
    distance = _distance(coords[i - 1], coords[j - 1])
    patch = 0
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            q = (min(i + di, j + dj), max(i + di, j + dj))
            if q in raw_set:
                patch += 1
    degree_i = degrees[i - 1]
    degree_j = degrees[j - 1]
    residue_score = _residue_pair_score(row.sequence[i - 1], row.sequence[j - 1], distance, seq_sep)
    return (
        1.0,
        max(0.0, min(1.0, float(raw_scores.get(pair, 0.0)))),
        max(0.0, min(1.0, float(restraints.get(pair, 0.0)))),
        max(0.0, min(1.0, residue_score)),
        max(0.0, min(1.0, patch / 8.0)),
        max(0.0, min(1.0, seq_sep / max(1, row.sequence_length))),
        max(0.0, min(1.0, exp(-abs(distance - 7.2) / 5.0))),
        max(0.0, min(1.0, (degree_i + degree_j) / 22.0)),
        1.0 if seq_sep >= 24 else 0.0,
    )


def _candidate_features(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    raw_contacts: Sequence[ContactPair],
    raw_scores: Mapping[ContactPair, float],
    restraints: Sequence[CouplingConstraint],
) -> dict[ContactPair, tuple[float, ...]]:
    raw_set = set(raw_contacts)
    degrees = [0] * row.sequence_length
    for i, j in raw_contacts:
        degrees[i - 1] += 1
        degrees[j - 1] += 1
    restraint_map = _restraint_confidence(restraints)
    return {
        pair: _feature_vector(row, coords, pair, raw_scores, raw_set, restraint_map, degrees)
        for pair in raw_contacts
    }


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = exp(-x)
        return 1.0 / (1.0 + z)
    z = exp(x)
    return z / (1.0 + z)


def _train_logistic_weights(samples: Sequence[tuple[tuple[float, ...], int]], *, iterations: int = 180, lr: float = 0.055) -> tuple[float, ...]:
    if not samples:
        return (0.0,) * 9
    dim = len(samples[0][0])
    weights = [0.0] * dim
    positives = sum(label for _features, label in samples)
    negatives = len(samples) - positives
    pos_weight = min(16.0, max(1.0, negatives / max(1, positives)))
    neg_weight = 1.0
    for _ in range(iterations):
        grad = [0.0] * dim
        for features, label in samples:
            pred = _sigmoid(sum(w * x for w, x in zip(weights, features)))
            sample_weight = pos_weight if label else neg_weight
            err = (pred - label) * sample_weight
            for k, value in enumerate(features):
                grad[k] += err * value
        scale = 1.0 / max(1, len(samples))
        for k in range(dim):
            weights[k] -= lr * (grad[k] * scale + 0.008 * weights[k])
    return tuple(_rounded(w) for w in weights)


def _make_candidate(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    *,
    md_steps: int,
    restart_index: int,
    max_direct: int | None,
    max_sequence_closure: int | None,
    max_geodesic: int | None,
) -> tuple[tuple[CouplingConstraint, ...], list[Vector], tuple[ContactPair, ...], dict[ContactPair, float], float]:
    restraints, _summary, _records = build_geodesic_consensus_filter_constraints(
        row,
        constraints,
        max_direct=max_direct,
        max_sequence_closure=max_sequence_closure,
        max_geodesic=max_geodesic,
        max_geodesic_distance=30,
        max_bending_energy=0.50,
        consensus_vote_threshold=3,
    )
    coords, energy = _relax_global_geometry(row=row, restraints=restraints, steps=md_steps, restart_index=restart_index)
    raw_contacts, raw_scores = _extract_contacts_from_geometry(row, coords)
    return tuple(restraints), coords, raw_contacts, raw_scores, _rounded(energy)


def _training_samples_for_target(
    rows: Sequence[RealCoordinateVisualRow],
    target_row_id: str,
    constraints: Sequence[CouplingConstraint],
    *,
    md_steps: int,
    max_direct: int | None,
    max_sequence_closure: int | None,
    max_geodesic: int | None,
    max_samples_per_row: int = 360,
) -> tuple[tuple[tuple[float, ...], int], ...]:
    samples: list[tuple[tuple[float, ...], int]] = []
    for index, row in enumerate(rows):
        if row.row_id == target_row_id:
            continue
        restraints, coords, raw_contacts, raw_scores, _energy = _make_candidate(
            row,
            constraints,
            md_steps=max(8, md_steps // 2),
            restart_index=index % 3,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
        )
        features = _candidate_features(row, coords, raw_contacts, raw_scores, restraints)
        native = set(row.native_contact_pairs())
        positives = [pair for pair in raw_contacts if pair in native]
        negatives = [pair for pair in raw_contacts if pair not in native]
        # Deterministic balanced-ish sample. Native is used only for non-target rows.
        for pair in positives[: max_samples_per_row // 2]:
            samples.append((features[pair], 1))
        # Spread negative examples across score/order rather than just taking the first hairball.
        step = max(1, len(negatives) // max(1, max_samples_per_row - len(positives[: max_samples_per_row // 2])))
        for pair in negatives[::step][: max_samples_per_row - len(positives[: max_samples_per_row // 2])]:
            samples.append((features[pair], 0))
    return tuple(samples)


def calibrated_filter_contacts(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    raw_contacts: Sequence[ContactPair],
    raw_scores: Mapping[ContactPair, float],
    restraints: Sequence[CouplingConstraint],
    weights: Sequence[float],
    *,
    probability_threshold: float = 0.42,
) -> tuple[tuple[ContactPair, ...], dict[ContactPair, float]]:
    if not raw_contacts:
        return (), {}
    features = _candidate_features(row, coords, raw_contacts, raw_scores, restraints)
    scored = [(pair, _sigmoid(sum(w * x for w, x in zip(weights, features[pair])))) for pair in raw_contacts]
    scored.sort(key=lambda item: (-item[1], item[0][0], item[0][1]))
    max_contacts = max(20, int(round(row.sequence_length * 2.05)))
    min_long = max(5, int(round(row.sequence_length * 0.30)))
    degree_cap = 8
    selected: list[ContactPair] = []
    selected_scores: dict[ContactPair, float] = {}
    degrees = [0] * row.sequence_length
    long_count = 0
    for pair, prob in scored:
        i, j = pair
        if prob < probability_threshold and pair[1] - pair[0] < 24:
            continue
        if prob < probability_threshold * 0.84 and pair[1] - pair[0] >= 24:
            continue
        if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
            continue
        selected.append(pair)
        selected_scores[pair] = _rounded(prob)
        degrees[i - 1] += 1
        degrees[j - 1] += 1
        if j - i >= 24:
            long_count += 1
        if len(selected) >= max_contacts:
            break
    if long_count < min_long:
        for pair, prob in scored:
            if pair in selected or pair[1] - pair[0] < 24:
                continue
            i, j = pair
            if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
                continue
            selected.append(pair)
            selected_scores[pair] = _rounded(prob)
            degrees[i - 1] += 1
            degrees[j - 1] += 1
            long_count += 1
            if long_count >= min_long or len(selected) >= max_contacts:
                break
    return normalized_contact_pairs(selected), selected_scores


def calibrated_global_score(row: RealCoordinateVisualRow, coords: Sequence[Sequence[float]], contacts: Sequence[ContactPair], restraints: Sequence[CouplingConstraint]) -> tuple[float, dict[str, float]]:
    stat = statistical_potential_surrogate(row, coords, contacts)
    degree = degree_coherence_score(contacts, row.sequence_length)
    patch = loop_patch_coherence_score(contacts)
    restraint = restraint_satisfaction_score(coords, restraints)
    compact = compactness_score(coords)
    collision = collision_penalty(coords)
    order = _contact_order(contacts, row.sequence_length)
    order_score = exp(-((order - 0.18) ** 2) / 0.045)
    score = 0.24 * stat + 0.18 * degree + 0.16 * patch + 0.16 * restraint + 0.10 * compact + 0.24 * order_score - 1.05 * collision
    return _rounded(score), {
        "statistical_potential_score": _rounded(stat),
        "degree_coherence_score": _rounded(degree),
        "loop_patch_coherence_score": _rounded(patch),
        "restraint_satisfaction_score": _rounded(restraint),
        "compactness_score": _rounded(compact),
        "collision_penalty": _rounded(collision),
        "contact_order": _rounded(order),
        "contact_order_score": _rounded(order_score),
    }


@dataclass(frozen=True)
class CalibratedVerifierRowReport:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    optional_backend_status: dict[str, object]
    training_row_count: int
    training_sample_count: int
    training_positive_count: int
    target_native_excluded_from_training: bool
    pdb_calibration_native_contacts_used_from_other_rows: bool
    calibrated_weights_sha256: str
    safe_restraint_count: int
    candidate_structure_count: int
    selected_candidate_index: int
    selected_calibrated_global_score: float
    statistical_potential_score: float
    degree_coherence_score: float
    loop_patch_coherence_score: float
    restraint_satisfaction_score: float
    compactness_score: float
    collision_penalty: float
    contact_order: float
    contact_order_score: float
    raw_extracted_contact_count: int
    calibrated_filtered_contact_count: int
    calibrated_filtered_long_range_contact_count: int
    selected_contact_map_hash: str
    final_energy_proxy: float
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_calibrated_verifier_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_target_selection: bool = False
    native_truth_used_before_target_selection: bool = False
    target_native_truth_attached_after_selection_for_evaluation: bool = True
    learned_geometry_prior_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = True

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class CalibratedTopologyVerifierPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    optional_backend_status: dict[str, object]
    calibrated_structural_feature_discriminator_included: bool
    leave_one_target_out_calibration: bool
    contact_order_criterion_included: bool
    dfire_or_pyrosetta_backend_used: bool
    pdb_calibration_native_contacts_used_from_other_rows: bool
    target_native_excluded_before_selection_for_all_rows: bool
    mean_training_sample_count: float
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    calibrated_verifier_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[CalibratedVerifierRowReport, ...]
    decisions: tuple[CoarseGrainMDContactDecision, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "optional_backend_status": self.optional_backend_status,
            "calibrated_structural_feature_discriminator_included": self.calibrated_structural_feature_discriminator_included,
            "leave_one_target_out_calibration": self.leave_one_target_out_calibration,
            "contact_order_criterion_included": self.contact_order_criterion_included,
            "dfire_or_pyrosetta_backend_used": self.dfire_or_pyrosetta_backend_used,
            "pdb_calibration_native_contacts_used_from_other_rows": self.pdb_calibration_native_contacts_used_from_other_rows,
            "target_native_excluded_before_selection_for_all_rows": self.target_native_excluded_before_selection_for_all_rows,
            "mean_training_sample_count": self.mean_training_sample_count,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "calibrated_verifier_claim_allowed": self.calibrated_verifier_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [r.to_dict() for r in self.rows],
            "decisions": [d.to_dict() for d in self.decisions],
        }


def _claim(metric: ContactMetricPacket, f1_margin: float, lr_margin: float) -> tuple[bool, bool, str]:
    ok = metric.native_contact_precision >= 0.70 and metric.native_contact_recall >= 0.70 and metric.long_range_contact_recall >= 0.70 and f1_margin >= 0.15 and lr_margin >= 0.15
    if ok:
        return True, False, "calibrated_verifier_survived_gate_not_universal_physics_due_to_pdb_calibration"
    if metric.native_contact_precision < 0.70:
        return False, False, "row_claim_rejected_precision_below_0_70"
    if metric.native_contact_recall < 0.70:
        return False, False, "row_claim_rejected_recall_below_0_70"
    if metric.long_range_contact_recall < 0.70:
        return False, False, "row_claim_rejected_long_range_recall_below_0_70"
    return False, False, "row_claim_rejected_matched_control_margin_below_gate"


def _decisions(row: RealCoordinateVisualRow, contacts: Sequence[ContactPair], scores: Mapping[ContactPair, float]) -> tuple[CoarseGrainMDContactDecision, ...]:
    out: list[CoarseGrainMDContactDecision] = []
    for pair in sorted(contacts, key=lambda p: (-scores.get(p, 0.0), p[0], p[1]))[:700]:
        score = float(scores.get(pair, 0.0))
        out.append(CoarseGrainMDContactDecision(
            kind=CALIBRATED_VERIFIER_DECISION_KIND,
            row_id=row.row_id,
            source_accession=row.source_accession,
            source_mode=CALIBRATED_VERIFIER_MODE,
            i=pair[0],
            j=pair[1],
            sequence_separation=pair[1] - pair[0],
            final_distance_angstrom=_rounded(8.2 - 2.0 * min(score, 1.0)),
            geometry_contact_score=_rounded(score),
            selected=True,
            selected_from_final_structure=True,
            msa_dca_used_before_selection=True,
        ))
    return tuple(out)


def run_calibrated_verifier_row(
    *,
    row: RealCoordinateVisualRow,
    all_rows: Sequence[RealCoordinateVisualRow],
    constraints: Sequence[CouplingConstraint],
    md_steps: int = 64,
    candidate_count: int = 4,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
) -> tuple[CalibratedVerifierRowReport, tuple[CoarseGrainMDContactDecision, ...]]:
    samples = _training_samples_for_target(
        all_rows,
        row.row_id,
        constraints,
        md_steps=md_steps,
        max_direct=max_direct,
        max_sequence_closure=max_sequence_closure,
        max_geodesic=max_geodesic,
    )
    weights = _train_logistic_weights(samples)
    weights_hash = hashlib.sha256(",".join(str(w) for w in weights).encode("utf-8")).hexdigest()
    positives = sum(label for _features, label in samples)
    backend = optional_energy_backend_status()
    best: tuple[float, int, tuple[CouplingConstraint, ...], list[Vector], tuple[ContactPair, ...], dict[ContactPair, float], tuple[ContactPair, ...], dict[ContactPair, float], dict[str, float], float] | None = None
    for idx in range(max(1, candidate_count)):
        restraints, coords, raw_contacts, raw_scores, energy = _make_candidate(
            row,
            constraints,
            md_steps=md_steps + 5 * (idx % 3),
            restart_index=idx,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
        )
        selected, selected_scores = calibrated_filter_contacts(row, coords, raw_contacts, raw_scores, restraints, weights)
        score, components = calibrated_global_score(row, coords, selected, restraints)
        if best is None or score > best[0]:
            best = (score, idx, restraints, coords, raw_contacts, raw_scores, selected, selected_scores, components, energy)
    assert best is not None
    selected_score, selected_idx, restraints, coords, raw_contacts, _raw_scores, selected, selected_scores, components, energy = best
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(row=row, selected_pairs=selected)
    f1_margin = _rounded(metric.contact_map_f1 - best_control_f1)
    lr_margin = _rounded(metric.long_range_contact_recall - best_control_lr)
    row_claim, row_universal, reason = _claim(metric, f1_margin, lr_margin)
    report = CalibratedVerifierRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        optional_backend_status=backend,
        training_row_count=max(0, len(all_rows) - 1),
        training_sample_count=len(samples),
        training_positive_count=positives,
        target_native_excluded_from_training=True,
        pdb_calibration_native_contacts_used_from_other_rows=True,
        calibrated_weights_sha256=weights_hash,
        safe_restraint_count=len(restraints),
        candidate_structure_count=max(1, candidate_count),
        selected_candidate_index=selected_idx,
        selected_calibrated_global_score=selected_score,
        statistical_potential_score=components["statistical_potential_score"],
        degree_coherence_score=components["degree_coherence_score"],
        loop_patch_coherence_score=components["loop_patch_coherence_score"],
        restraint_satisfaction_score=components["restraint_satisfaction_score"],
        compactness_score=components["compactness_score"],
        collision_penalty=components["collision_penalty"],
        contact_order=components["contact_order"],
        contact_order_score=components["contact_order_score"],
        raw_extracted_contact_count=len(raw_contacts),
        calibrated_filtered_contact_count=len(selected),
        calibrated_filtered_long_range_contact_count=sum(1 for i, j in selected if j - i >= 24),
        selected_contact_map_hash=_pair_hash(selected),
        final_energy_proxy=_rounded(energy),
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_calibrated_verifier_claim_allowed=row_claim,
        row_universal_physical_law_claim_allowed=row_universal,
        row_claim_rejection_reason=reason,
    )
    return report, _decisions(row, selected, selected_scores)


def run_calibrated_topology_verifier_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint],
    evaluation_source_accessions: Sequence[str] | None = None,
    md_steps: int = 64,
    candidate_count: int = 4,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
) -> CalibratedTopologyVerifierPacket:
    reports: list[CalibratedVerifierRowReport] = []
    decisions: list[CoarseGrainMDContactDecision] = []
    evaluation_set = set(evaluation_source_accessions or ())
    evaluation_rows = tuple(row for row in rows if not evaluation_set or row.source_accession in evaluation_set)
    for row in evaluation_rows:
        report, row_decisions = run_calibrated_verifier_row(
            row=row,
            all_rows=rows,
            constraints=constraints,
            md_steps=md_steps,
            candidate_count=candidate_count,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
        )
        reports.append(report)
        decisions.extend(row_decisions[:420])
    precisions = [r.metric_after_native_audit.native_contact_precision for r in reports]
    recalls = [r.metric_after_native_audit.native_contact_recall for r in reports]
    long_ranges = [r.metric_after_native_audit.long_range_contact_recall for r in reports]
    f1s = [r.metric_after_native_audit.contact_map_f1 for r in reports]
    f1_margins = [r.f1_margin_vs_best_control for r in reports]
    lr_margins = [r.long_range_recall_margin_vs_best_control for r in reports]
    claim = bool(reports) and all(r.row_calibrated_verifier_claim_allowed for r in reports) and _mean(precisions) >= 0.70 and _mean(recalls) >= 0.70 and _mean(long_ranges) >= 0.70
    if claim:
        reason = "calibrated_verifier_survived_gate_not_universal_physics_due_to_pdb_calibration"
    else:
        reason = "calibrated_verifier_claim_rejected_for_rows:" + ",".join(r.source_accession for r in reports if not r.row_calibrated_verifier_claim_allowed)[:180]
    return CalibratedTopologyVerifierPacket(
        kind=CALIBRATED_VERIFIER_KIND,
        source_mode=CALIBRATED_VERIFIER_MODE,
        row_count=len(reports),
        decision_rule=CALIBRATED_VERIFIER_RULE,
        claim_rule=CALIBRATED_VERIFIER_CLAIM_RULE,
        optional_backend_status=optional_energy_backend_status(),
        calibrated_structural_feature_discriminator_included=True,
        leave_one_target_out_calibration=True,
        contact_order_criterion_included=True,
        dfire_or_pyrosetta_backend_used=False,
        pdb_calibration_native_contacts_used_from_other_rows=True,
        target_native_excluded_before_selection_for_all_rows=True,
        mean_training_sample_count=_mean([r.training_sample_count for r in reports]),
        mean_native_contact_precision_after_audit=_mean(precisions),
        mean_native_contact_recall_after_audit=_mean(recalls),
        mean_long_range_contact_recall_after_audit=_mean(long_ranges),
        mean_contact_map_f1_after_audit=_mean(f1s),
        min_row_precision_after_audit=_rounded(min(precisions)) if precisions else 0.0,
        min_row_recall_after_audit=_rounded(min(recalls)) if recalls else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_ranges)) if long_ranges else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        calibrated_verifier_claim_allowed=claim,
        universal_physical_law_claim_allowed=False,
        folding_problem_solved=claim,
        claim_rejection_reason=reason,
        rows=tuple(reports),
        decisions=tuple(decisions),
    )
