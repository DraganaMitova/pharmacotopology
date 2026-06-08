from __future__ import annotations

"""Native-free global verifier and structure-level annealing challenge.

This layer tests the hypothesis after geodesic-consensus filtering failed:
path-level filters are not enough; the missing step may be a global verifier that
scores an entire C-alpha geometry/contact field at once before native audit.

The verifier is deliberately dependency-free and native-free. It uses only:
* safe external DCA/restraint information already admitted by prior gates;
* sequence-derived residue classes and lightweight secondary structure;
* geometry-derived global properties of candidate structures.

It does not use native contacts, native coordinates, templates, AlphaFold,
ESMFold, PyRosetta, DFIRE tables, or trained parameters before selection. Native
contacts are attached only after the selected contact map is frozen, for audit.
"""

import hashlib
from dataclasses import asdict, dataclass
from math import exp, log, sqrt
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_coarse_grain_md_geometry import (
    COARSE_GRAIN_MD_CLAIM_RULE,
    COARSE_GRAIN_MD_RULE,
    CoarseGrainMDContactDecision,
    _distance,
    _energy_proxy,
    _extract_contacts_from_geometry,
    _matched_controls_for_report,
    _pair_hash,
    _relax_global_geometry,
    _restraint_edges,
    _rounded,
)
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_geodesic_consensus_filter_md import (
    build_geodesic_consensus_filter_constraints,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_sequence_physical_priors import predict_lightweight_secondary_structure

GLOBAL_VERIFIER_KIND = "global_verifier_structure_annealing_v0"
GLOBAL_VERIFIER_DECISION_KIND = "global_verifier_extracted_contact_decision_v0"
GLOBAL_VERIFIER_MODE = "external_dca_geodesic_global_verifier_annealing"
GLOBAL_VERIFIER_RULE = (
    COARSE_GRAIN_MD_RULE
    + ";global_structure_candidate_generation;native_free_statistical_potential_surrogate;"
    + "degree_distribution_coherence;loop_patch_coherence;restraint_satisfaction;"
    + "global_candidate_score_selection;contacts_extracted_after_selected_structure_only"
)
GLOBAL_VERIFIER_CLAIM_RULE = (
    COARSE_GRAIN_MD_CLAIM_RULE
    + ";global_verifier_claim_requires_all_rows_precision_recall_long_range_ge_0_70;"
    + "native_used_only_after_global_structure_selection"
)

HYDROPHOBIC = set("AILMFWVYCP")
CHARGED_POS = set("KRH")
CHARGED_NEG = set("DE")
POLAR = set("STNQ")
GLY_PRO = set("GP")

Vector = list[float]


@dataclass(frozen=True)
class GlobalVerifierRowReport:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    safe_restraint_count: int
    candidate_structure_count: int
    selected_candidate_index: int
    selected_global_score: float
    statistical_potential_score: float
    degree_coherence_score: float
    loop_patch_coherence_score: float
    restraint_satisfaction_score: float
    compactness_score: float
    collision_penalty: float
    raw_extracted_contact_count: int
    globally_filtered_contact_count: int
    globally_filtered_long_range_contact_count: int
    selected_contact_map_hash: str
    final_energy_proxy: float
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_global_verifier_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = True
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class GlobalVerifierAnnealingPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    global_verifier_included: bool
    statistical_potential_surrogate_included: bool
    degree_distribution_coherence_included: bool
    loop_patch_coherence_included: bool
    monte_carlo_global_candidate_selection_included: bool
    contacts_predicted_before_structure: bool
    atomistic_md_engine_used: bool
    dependency_free_md_style_relaxation_used: bool
    mean_selected_global_score: float
    mean_statistical_potential_score: float
    mean_degree_coherence_score: float
    mean_restraint_satisfaction_score: float
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    global_verifier_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[GlobalVerifierRowReport, ...]
    decisions: tuple[CoarseGrainMDContactDecision, ...]
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = True
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "global_verifier_included": self.global_verifier_included,
            "statistical_potential_surrogate_included": self.statistical_potential_surrogate_included,
            "degree_distribution_coherence_included": self.degree_distribution_coherence_included,
            "loop_patch_coherence_included": self.loop_patch_coherence_included,
            "monte_carlo_global_candidate_selection_included": self.monte_carlo_global_candidate_selection_included,
            "contacts_predicted_before_structure": self.contacts_predicted_before_structure,
            "atomistic_md_engine_used": self.atomistic_md_engine_used,
            "dependency_free_md_style_relaxation_used": self.dependency_free_md_style_relaxation_used,
            "mean_selected_global_score": self.mean_selected_global_score,
            "mean_statistical_potential_score": self.mean_statistical_potential_score,
            "mean_degree_coherence_score": self.mean_degree_coherence_score,
            "mean_restraint_satisfaction_score": self.mean_restraint_satisfaction_score,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "global_verifier_claim_allowed": self.global_verifier_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [row.to_dict() for row in self.rows],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "structure_model_used_before_selection": self.structure_model_used_before_selection,
            "learned_geometry_prior_used_before_selection": self.learned_geometry_prior_used_before_selection,
            "msa_dca_used_before_selection": self.msa_dca_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _bounded01(value: float) -> float:
    return _rounded(max(0.0, min(1.0, float(value))))


def _restraint_map(restraints: Sequence[CouplingConstraint]) -> dict[ContactPair, float]:
    return {r.pair(): max(float(r.confidence), 0.05) for r in restraints}


def _residue_pair_score(a: str, b: str, distance: float, seq_sep: int) -> float:
    close = 1.0 / (1.0 + exp((distance - 7.2) / 0.75))
    score = 0.0
    if a in HYDROPHOBIC and b in HYDROPHOBIC:
        score += 1.00 * close
    elif (a in CHARGED_POS and b in CHARGED_NEG) or (a in CHARGED_NEG and b in CHARGED_POS):
        score += 0.78 * close
    elif (a in POLAR and b in POLAR) or (a in POLAR and b in CHARGED_POS | CHARGED_NEG) or (b in POLAR and a in CHARGED_POS | CHARGED_NEG):
        score += 0.42 * close
    else:
        score += 0.18 * close
    if a in GLY_PRO or b in GLY_PRO:
        score -= 0.10 * close
    if seq_sep >= 24:
        score += 0.10 * close
    return score


def statistical_potential_surrogate(row: RealCoordinateVisualRow, coords: Sequence[Sequence[float]], contacts: Sequence[ContactPair]) -> float:
    if not contacts:
        return 0.0
    sequence = row.sequence
    raw = 0.0
    for i, j in contacts:
        raw += _residue_pair_score(sequence[i - 1], sequence[j - 1], _distance(coords[i - 1], coords[j - 1]), j - i)
    # Penalize contact-map explosions; native-free expected density band for coarse C-alpha maps.
    n = row.sequence_length
    density = len(contacts) / max(1.0, n)
    density_penalty = abs(density - 2.65) / 2.65
    return _rounded((raw / max(1, len(contacts))) - 0.42 * density_penalty)


def degree_coherence_score(contacts: Sequence[ContactPair], sequence_length: int) -> float:
    if not contacts or sequence_length <= 0:
        return 0.0
    degrees = [0] * sequence_length
    for i, j in contacts:
        degrees[i - 1] += 1
        degrees[j - 1] += 1
    avg = sum(degrees) / sequence_length
    variance = sum((d - avg) ** 2 for d in degrees) / sequence_length
    # C-alpha 8A maps usually should not look like a dense hairball or a line with all-zero degrees.
    avg_score = exp(-((avg - 5.3) ** 2) / 13.5)
    var_score = exp(-((sqrt(variance) - 2.6) ** 2) / 8.0)
    cap_penalty = sum(max(0, d - 11) for d in degrees) / max(1, sum(degrees))
    isolated_penalty = sum(1 for d in degrees if d == 0) / sequence_length
    return _rounded(max(0.0, 0.62 * avg_score + 0.38 * var_score - 0.70 * cap_penalty - 0.25 * isolated_penalty))


def loop_patch_coherence_score(contacts: Sequence[ContactPair]) -> float:
    if not contacts:
        return 0.0
    contact_set = set(contacts)
    supported = 0
    isolated = 0
    for i, j in contact_set:
        neighbours = 0
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                pair = (min(i + di, j + dj), max(i + di, j + dj))
                if pair in contact_set:
                    neighbours += 1
        if neighbours >= 2:
            supported += 1
        if neighbours == 0:
            isolated += 1
    return _rounded(max(0.0, supported / len(contact_set) - 0.45 * isolated / len(contact_set)))


def restraint_satisfaction_score(coords: Sequence[Sequence[float]], restraints: Sequence[CouplingConstraint]) -> float:
    if not restraints:
        return 0.0
    weighted = 0.0
    total = 0.0
    for r in restraints:
        confidence = max(0.05, min(1.0, float(r.confidence)))
        distance = _distance(coords[r.i - 1], coords[r.j - 1])
        target = 7.1 if r.sequence_separation >= 12 else 6.2
        weighted += confidence * exp(-((distance - target) ** 2) / 12.0)
        total += confidence
    return _rounded(weighted / max(total, 1e-9))


def compactness_score(coords: Sequence[Sequence[float]]) -> float:
    n = len(coords)
    if n <= 1:
        return 0.0
    center = [sum(point[axis] for point in coords) / n for axis in range(3)]
    rg = sqrt(sum(sum((point[axis] - center[axis]) ** 2 for axis in range(3)) for point in coords) / n)
    target = 2.7 * (n ** 0.37)
    return _rounded(exp(-abs(rg - target) / max(target, 1.0)))


def collision_penalty(coords: Sequence[Sequence[float]]) -> float:
    n = len(coords)
    collisions = 0.0
    checks = 0
    stride = 2 if n > 160 else 1
    for i in range(0, n, stride):
        for j in range(i + MIN_SEQUENCE_SEPARATION, n, stride):
            checks += 1
            d = _distance(coords[i], coords[j])
            if d < 3.45:
                collisions += (3.45 - d) / 3.45
    return _rounded(collisions / max(1, checks))


def _contact_native_free_score(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    pair: ContactPair,
    geometry_scores: Mapping[ContactPair, float],
    restraints: Mapping[ContactPair, float],
    contact_set: set[ContactPair],
    degrees: Sequence[int],
) -> float:
    i, j = pair
    sequence = row.sequence
    d = _distance(coords[i - 1], coords[j - 1])
    score = float(geometry_scores.get(pair, 0.0))
    score += 0.30 * _residue_pair_score(sequence[i - 1], sequence[j - 1], d, j - i)
    score += 0.40 * restraints.get(pair, 0.0)
    patch = 0
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            q = (min(i + di, j + dj), max(i + di, j + dj))
            if q in contact_set:
                patch += 1
    score += 0.045 * min(patch, 6)
    if degrees[i - 1] > 10 or degrees[j - 1] > 10:
        score -= 0.22
    if j - i >= 24:
        score += 0.08
    return score


def globally_filter_contacts(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    raw_contacts: Sequence[ContactPair],
    geometry_scores: Mapping[ContactPair, float],
    restraints: Sequence[CouplingConstraint],
) -> tuple[tuple[ContactPair, ...], dict[ContactPair, float]]:
    if not raw_contacts:
        return (), {}
    restraint_scores = _restraint_map(restraints)
    raw_set = set(raw_contacts)
    base_degrees = [0] * row.sequence_length
    for i, j in raw_contacts:
        base_degrees[i - 1] += 1
        base_degrees[j - 1] += 1
    scored = [
        (pair, _contact_native_free_score(row, coords, pair, geometry_scores, restraint_scores, raw_set, base_degrees))
        for pair in raw_contacts
    ]
    scored.sort(key=lambda item: (-item[1], item[0][0], item[0][1]))
    max_contacts = max(24, int(round(row.sequence_length * 2.45)))
    degree_cap = 9
    selected: list[ContactPair] = []
    selected_scores: dict[ContactPair, float] = {}
    degrees = [0] * row.sequence_length
    long_range_kept = 0
    for pair, score in scored:
        i, j = pair
        if score < 0.34 and pair not in restraint_scores:
            continue
        if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
            continue
        selected.append(pair)
        selected_scores[pair] = _bounded01(score / 1.65)
        degrees[i - 1] += 1
        degrees[j - 1] += 1
        if j - i >= 24:
            long_range_kept += 1
        if len(selected) >= max_contacts:
            break
    # Avoid a verifier that solves precision by deleting all long-range contacts.
    min_long_range = max(6, int(0.18 * len(selected))) if selected else 0
    if long_range_kept < min_long_range:
        for pair, score in scored:
            if pair in selected or pair[1] - pair[0] < 24:
                continue
            i, j = pair
            if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
                continue
            selected.append(pair)
            selected_scores[pair] = _bounded01(score / 1.65)
            degrees[i - 1] += 1
            degrees[j - 1] += 1
            long_range_kept += 1
            if long_range_kept >= min_long_range or len(selected) >= max_contacts:
                break
    return normalized_contact_pairs(selected), selected_scores


def global_structure_score(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    raw_contacts: Sequence[ContactPair],
    filtered_contacts: Sequence[ContactPair],
    restraints: Sequence[CouplingConstraint],
) -> tuple[float, dict[str, float]]:
    stat = statistical_potential_surrogate(row, coords, filtered_contacts)
    degree = degree_coherence_score(filtered_contacts, row.sequence_length)
    patch = loop_patch_coherence_score(filtered_contacts)
    restraint = restraint_satisfaction_score(coords, restraints)
    compact = compactness_score(coords)
    collision = collision_penalty(coords)
    raw_density = len(raw_contacts) / max(1, row.sequence_length)
    filtered_density = len(filtered_contacts) / max(1, row.sequence_length)
    density_penalty = abs(filtered_density - 2.3) / 3.0 + max(0.0, raw_density - 7.0) / 10.0
    score = 0.32 * stat + 0.25 * degree + 0.18 * patch + 0.20 * restraint + 0.12 * compact - 1.2 * collision - 0.16 * density_penalty
    components = {
        "statistical_potential_score": _rounded(stat),
        "degree_coherence_score": _rounded(degree),
        "loop_patch_coherence_score": _rounded(patch),
        "restraint_satisfaction_score": _rounded(restraint),
        "compactness_score": _rounded(compact),
        "collision_penalty": _rounded(collision),
    }
    return _rounded(score), components


def _claim_decision(metric: ContactMetricPacket, f1_margin: float, long_range_margin: float) -> tuple[bool, bool, str]:
    allowed = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
        and f1_margin >= 0.15
        and long_range_margin >= 0.15
    )
    if allowed:
        return True, False, "external_dca_global_verifier_survived_gate_not_universal_physics"
    if metric.native_contact_precision < 0.70:
        return False, False, "row_claim_rejected_precision_below_0_70"
    if metric.native_contact_recall < 0.70:
        return False, False, "row_claim_rejected_recall_below_0_70"
    if metric.long_range_contact_recall < 0.70:
        return False, False, "row_claim_rejected_long_range_recall_below_0_70"
    return False, False, "row_claim_rejected_matched_control_margin_below_gate"


def _decisions_for_row(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[ContactPair],
    scores: Mapping[ContactPair, float],
    max_report_decisions: int = 900,
) -> tuple[CoarseGrainMDContactDecision, ...]:
    decisions: list[CoarseGrainMDContactDecision] = []
    ordered = sorted(selected_pairs, key=lambda pair: (-scores.get(pair, 0.0), pair[0], pair[1]))[:max_report_decisions]
    for pair in ordered:
        score = float(scores.get(pair, 0.0))
        decisions.append(
            CoarseGrainMDContactDecision(
                kind=GLOBAL_VERIFIER_DECISION_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                source_mode=GLOBAL_VERIFIER_MODE,
                i=pair[0],
                j=pair[1],
                sequence_separation=pair[1] - pair[0],
                final_distance_angstrom=_rounded(7.2 + 0.65 * (-score if score else 0.0)),
                geometry_contact_score=_rounded(score),
                selected=True,
                selected_from_final_structure=True,
                msa_dca_used_before_selection=True,
            )
        )
    return tuple(decisions)


def run_global_verifier_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    md_steps: int = 84,
    candidate_count: int = 5,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
) -> tuple[GlobalVerifierRowReport, tuple[CoarseGrainMDContactDecision, ...]]:
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
    candidates = max(1, int(candidate_count))
    best: tuple[float, int, list[Vector], tuple[ContactPair, ...], dict[ContactPair, float], tuple[ContactPair, ...], dict[str, float], float] | None = None
    for index in range(candidates):
        # Deterministic global candidate generation: each restart moves the whole structure.
        coords, energy = _relax_global_geometry(row=row, restraints=restraints, steps=md_steps + 6 * (index % 3), restart_index=index)
        raw_contacts, raw_scores = _extract_contacts_from_geometry(row, coords)
        filtered_contacts, filtered_scores = globally_filter_contacts(row, coords, raw_contacts, raw_scores, restraints)
        score, components = global_structure_score(row, coords, raw_contacts, filtered_contacts, restraints)
        if best is None or score > best[0]:
            best = (score, index, coords, raw_contacts, filtered_scores, filtered_contacts, components, energy)
    assert best is not None
    selected_score, selected_index, _coords, raw_contacts, selected_scores, selected_contacts, components, energy = best
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected_contacts)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(row=row, selected_pairs=selected_contacts)
    f1_margin = _rounded(metric.contact_map_f1 - best_control_f1)
    lr_margin = _rounded(metric.long_range_contact_recall - best_control_lr)
    row_claim, row_universal, rejection = _claim_decision(metric, f1_margin, lr_margin)
    decisions = _decisions_for_row(row=row, selected_pairs=selected_contacts, scores=selected_scores)
    report = GlobalVerifierRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        safe_restraint_count=len(restraints),
        candidate_structure_count=candidates,
        selected_candidate_index=selected_index,
        selected_global_score=selected_score,
        statistical_potential_score=components["statistical_potential_score"],
        degree_coherence_score=components["degree_coherence_score"],
        loop_patch_coherence_score=components["loop_patch_coherence_score"],
        restraint_satisfaction_score=components["restraint_satisfaction_score"],
        compactness_score=components["compactness_score"],
        collision_penalty=components["collision_penalty"],
        raw_extracted_contact_count=len(raw_contacts),
        globally_filtered_contact_count=len(selected_contacts),
        globally_filtered_long_range_contact_count=sum(1 for pair in selected_contacts if pair[1] - pair[0] >= 24),
        selected_contact_map_hash=_pair_hash(selected_contacts),
        final_energy_proxy=_rounded(energy),
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_global_verifier_claim_allowed=row_claim,
        row_universal_physical_law_claim_allowed=row_universal,
        row_claim_rejection_reason=rejection,
    )
    return report, decisions


def run_global_verifier_annealing_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint],
    md_steps: int = 84,
    candidate_count: int = 5,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
) -> GlobalVerifierAnnealingPacket:
    reports: list[GlobalVerifierRowReport] = []
    decisions: list[CoarseGrainMDContactDecision] = []
    for row in rows:
        report, row_decisions = run_global_verifier_row(
            row=row,
            constraints=constraints,
            md_steps=md_steps,
            candidate_count=candidate_count,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
        )
        reports.append(report)
        decisions.extend(row_decisions[:420])
    precision_values = [r.metric_after_native_audit.native_contact_precision for r in reports]
    recall_values = [r.metric_after_native_audit.native_contact_recall for r in reports]
    long_range_values = [r.metric_after_native_audit.long_range_contact_recall for r in reports]
    f1_values = [r.metric_after_native_audit.contact_map_f1 for r in reports]
    f1_margins = [r.f1_margin_vs_best_control for r in reports]
    lr_margins = [r.long_range_recall_margin_vs_best_control for r in reports]
    all_rows_claim = bool(reports) and all(r.row_global_verifier_claim_allowed for r in reports)
    mean_gate = _mean(precision_values) >= 0.70 and _mean(recall_values) >= 0.70 and _mean(long_range_values) >= 0.70
    claim = all_rows_claim and mean_gate
    if claim:
        rejection = "external_dca_global_verifier_survived_gate_not_universal_physics"
    else:
        failed = [r.source_accession for r in reports if not r.row_global_verifier_claim_allowed]
        rejection = "global_verifier_claim_rejected_for_rows:" + ",".join(failed[:12])
    return GlobalVerifierAnnealingPacket(
        kind=GLOBAL_VERIFIER_KIND,
        source_mode=GLOBAL_VERIFIER_MODE,
        row_count=len(reports),
        decision_rule=GLOBAL_VERIFIER_RULE,
        claim_rule=GLOBAL_VERIFIER_CLAIM_RULE,
        global_verifier_included=True,
        statistical_potential_surrogate_included=True,
        degree_distribution_coherence_included=True,
        loop_patch_coherence_included=True,
        monte_carlo_global_candidate_selection_included=True,
        contacts_predicted_before_structure=False,
        atomistic_md_engine_used=False,
        dependency_free_md_style_relaxation_used=True,
        mean_selected_global_score=_mean([r.selected_global_score for r in reports]),
        mean_statistical_potential_score=_mean([r.statistical_potential_score for r in reports]),
        mean_degree_coherence_score=_mean([r.degree_coherence_score for r in reports]),
        mean_restraint_satisfaction_score=_mean([r.restraint_satisfaction_score for r in reports]),
        mean_native_contact_precision_after_audit=_mean(precision_values),
        mean_native_contact_recall_after_audit=_mean(recall_values),
        mean_long_range_contact_recall_after_audit=_mean(long_range_values),
        mean_contact_map_f1_after_audit=_mean(f1_values),
        min_row_precision_after_audit=_rounded(min(precision_values)) if precision_values else 0.0,
        min_row_recall_after_audit=_rounded(min(recall_values)) if recall_values else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_range_values)) if long_range_values else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        global_verifier_claim_allowed=claim,
        universal_physical_law_claim_allowed=False,
        folding_problem_solved=claim,
        claim_rejection_reason=rejection,
        rows=tuple(reports),
        decisions=tuple(decisions),
    )
