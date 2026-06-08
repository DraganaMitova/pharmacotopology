from __future__ import annotations

"""Long-range calibrated topology verifier challenge.

This layer tests the next hypothesis: the previous calibrated verifier became too
conservative and erased long-range topology, so the verifier must be calibrated
on true long-range contacts from non-target reference structures.

Honesty boundary:
* The target row's native contacts/coordinates are never used before that row's
  contact map is frozen.
* The long-range potential is leave-one-target-out and PDB-derived, so it is a
  calibrated empirical verifier, not a universal physical law.
* Optional PyRosetta/DFIRE/ESMFold-teacher hooks are represented in the packet,
  but the default sandbox path is deterministic and dependency-light.
"""

from collections import defaultdict
from dataclasses import asdict, dataclass
from math import exp, log
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_calibrated_topology_verifier import (
    CALIBRATED_VERIFIER_DECISION_KIND,
    _candidate_features,
    _claim,
    _decisions,
    _make_candidate,
    _mean,
    _sigmoid,
    _train_logistic_weights,
    _training_samples_for_target,
    calibrated_global_score,
    optional_energy_backend_status,
)
from pharmacotopology.folding_coarse_grain_md_geometry import (
    CoarseGrainMDContactDecision,
    _matched_controls_for_report,
    _pair_hash,
    _rounded,
)
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_global_verifier_annealing import _residue_pair_score
from pharmacotopology.folding_native_contact_eval import ContactMetricPacket, ContactPair, evaluate_contact_prediction, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow

LONG_RANGE_CALIBRATED_KIND = "long_range_calibrated_topology_verifier_v0"
LONG_RANGE_CALIBRATED_MODE = "external_dca_geodesic_long_range_calibrated_verifier"
LONG_RANGE_DECISION_KIND = "long_range_calibrated_contact_decision_v0"
LONG_RANGE_RULE = (
    "leave_one_target_out_long_range_contact_potential;aa_pair_sep_bin_contact_probability;"
    "long_range_reserve_selection;structure_level_score_with_lr_potential;target_native_excluded_before_selection"
)
LONG_RANGE_CLAIM_RULE = (
    "claim_requires_precision_recall_long_range_ge_0_70_and_matched_control_margins;"
    "pdb_long_range_calibration_disallows_universal_physical_law_claim"
)

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
LONG_RANGE_THRESHOLD = 24
POTENTIAL_SEP_FLOOR = 30


def _sep_bin(sep: int) -> int:
    if sep < 30:
        return 0
    if sep < 40:
        return 30
    if sep < 60:
        return 40
    if sep < 90:
        return 60
    if sep < 130:
        return 90
    if sep < 200:
        return 130
    return 200


def _aa_pair(a: str, b: str) -> str:
    a = a if a in AA_ALPHABET else "X"
    b = b if b in AA_ALPHABET else "X"
    return "".join(sorted((a, b)))


@dataclass(frozen=True)
class LongRangePotential:
    target_row_id: str
    reference_row_count: int
    possible_long_range_pair_count: int
    native_long_range_contact_count: int
    base_probability: float
    probability_by_key: dict[str, float]
    odds_by_key: dict[str, float]
    max_probability: float
    target_native_excluded_from_potential: bool = True

    def key(self, aa_i: str, aa_j: str, sep: int) -> str:
        return f"{_aa_pair(aa_i, aa_j)}:{_sep_bin(sep)}"

    def probability(self, aa_i: str, aa_j: str, sep: int) -> float:
        if sep < POTENTIAL_SEP_FLOOR:
            return 0.0
        return self.probability_by_key.get(self.key(aa_i, aa_j, sep), self.base_probability)

    def score(self, aa_i: str, aa_j: str, sep: int) -> float:
        if sep < POTENTIAL_SEP_FLOOR:
            return 0.0
        raw = self.odds_by_key.get(self.key(aa_i, aa_j, sep), 0.0)
        return _rounded(1.0 / (1.0 + exp(-raw)))

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        # Keep reports compact but reproducible enough.
        payload["probability_key_count"] = len(self.probability_by_key)
        payload["top_probability_keys"] = sorted(self.probability_by_key.items(), key=lambda kv: (-kv[1], kv[0]))[:12]
        payload.pop("probability_by_key")
        payload.pop("odds_by_key")
        return payload


def build_leave_one_out_long_range_potential(rows: Sequence[RealCoordinateVisualRow], target_row_id: str) -> LongRangePotential:
    possible: dict[str, int] = defaultdict(int)
    contacts: dict[str, int] = defaultdict(int)
    possible_total = 0
    contact_total = 0
    ref_count = 0
    for row in rows:
        if row.row_id == target_row_id:
            continue
        ref_count += 1
        native = set(row.native_contact_pairs())
        seq = row.sequence
        for i in range(1, row.sequence_length + 1):
            for j in range(i + POTENTIAL_SEP_FLOOR, row.sequence_length + 1):
                key = f"{_aa_pair(seq[i - 1], seq[j - 1])}:{_sep_bin(j - i)}"
                possible[key] += 1
                possible_total += 1
                if (i, j) in native:
                    contacts[key] += 1
                    contact_total += 1
    # Conservative Beta prior, so rare keys do not explode.
    base = (contact_total + 1.0) / (possible_total + 2.0) if possible_total else 0.01
    probabilities: dict[str, float] = {}
    odds: dict[str, float] = {}
    for key, n in possible.items():
        c = contacts.get(key, 0)
        # stronger shrinkage for sparse reference set
        p = (c + 0.75 * base) / (n + 0.75)
        probabilities[key] = _rounded(p)
        odds[key] = _rounded(log(max(1e-6, p) / max(1e-6, base)))
    return LongRangePotential(
        target_row_id=target_row_id,
        reference_row_count=ref_count,
        possible_long_range_pair_count=possible_total,
        native_long_range_contact_count=contact_total,
        base_probability=_rounded(base),
        probability_by_key=probabilities,
        odds_by_key=odds,
        max_probability=_rounded(max(probabilities.values()) if probabilities else 0.0),
    )


def long_range_contact_score(pair: ContactPair, sequence: str, potential: LongRangePotential) -> float:
    i, j = pair
    return potential.score(sequence[i - 1], sequence[j - 1], j - i)


def chemical_score(aa_i: str, aa_j: str) -> float:
    hydrophobic = set("AVILMFWY")
    positive = set("KRH")
    negative = set("DE")
    polar = set("STNQ")
    if aa_i not in AA_ALPHABET or aa_j not in AA_ALPHABET:
        return 0.1
    if aa_i in hydrophobic and aa_j in hydrophobic:
        return 1.0
    if (aa_i in positive and aa_j in negative) or (aa_i in negative and aa_j in positive):
        return 0.9
    if aa_i in polar and aa_j in polar:
        return 0.5
    if (aa_i in polar and aa_j in hydrophobic) or (aa_i in hydrophobic and aa_j in polar):
        return 0.3
    return 0.1


def combined_score(i: int, j: int, sequence: str, stat_score: float, *, chem_weight: float = 0.7, stat_weight: float = 0.3) -> float:
    if i <= 0 or j <= 0 or i > len(sequence) or j > len(sequence):
        return stat_score
    chem = chemical_score(sequence[i - 1], sequence[j - 1])
    return chem_weight * chem + stat_weight * stat_score


def _patch_support(pair: ContactPair, raw_set: set[ContactPair]) -> float:
    i, j = pair
    patch = 0
    for di in (-2, -1, 0, 1, 2):
        for dj in (-2, -1, 0, 1, 2):
            if di == 0 and dj == 0:
                continue
            q = (min(i + di, j + dj), max(i + di, j + dj))
            if q in raw_set:
                patch += 1
    return max(0.0, min(1.0, patch / 16.0))


def long_range_aware_filter_contacts(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    raw_contacts: Sequence[ContactPair],
    raw_scores: Mapping[ContactPair, float],
    restraints: Sequence[CouplingConstraint],
    weights: Sequence[float],
    potential: LongRangePotential,
    *,
    probability_threshold: float = 0.34,
    use_chemical_score: bool = False,
) -> tuple[tuple[ContactPair, ...], dict[ContactPair, float]]:
    if not raw_contacts:
        return (), {}
    features = _candidate_features(row, coords, raw_contacts, raw_scores, restraints)
    restraint_map = {r.pair(): max(0.0, min(1.0, float(r.confidence))) for r in restraints}
    raw_set = set(raw_contacts)
    scored: list[tuple[ContactPair, float, float]] = []
    for pair in raw_contacts:
        calibrated = _sigmoid(sum(w * x for w, x in zip(weights, features[pair])))
        lr = long_range_contact_score(pair, row.sequence, potential)
        sep = pair[1] - pair[0]
        residue = _residue_pair_score(row.sequence[pair[0] - 1], row.sequence[pair[1] - 1], 7.2, sep)
        if use_chemical_score:
            residue = combined_score(pair[0], pair[1], row.sequence, residue)
        patch = _patch_support(pair, raw_set)
        restraint = restraint_map.get(pair, 0.0)
        raw = max(0.0, min(1.0, float(raw_scores.get(pair, 0.0))))
        if sep >= POTENTIAL_SEP_FLOOR:
            score = 0.26 * calibrated + 0.31 * lr + 0.17 * patch + 0.12 * raw + 0.10 * restraint + 0.04 * residue
        else:
            score = 0.42 * calibrated + 0.22 * patch + 0.16 * raw + 0.12 * restraint + 0.08 * residue
        scored.append((pair, _rounded(score), _rounded(lr)))
    scored.sort(key=lambda item: (-item[1], -item[2], item[0][0], item[0][1]))
    max_contacts = max(30, int(round(row.sequence_length * 2.35)))
    target_long = max(8, int(round(row.sequence_length * 0.32)))
    degree_cap = 7
    selected: list[ContactPair] = []
    scores: dict[ContactPair, float] = {}
    degrees = [0] * row.sequence_length
    long_count = 0

    def try_add(pair: ContactPair, score: float) -> bool:
        nonlocal long_count
        if pair in scores:
            return False
        i, j = pair
        if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
            return False
        selected.append(pair)
        scores[pair] = _rounded(score)
        degrees[i - 1] += 1
        degrees[j - 1] += 1
        if j - i >= LONG_RANGE_THRESHOLD:
            long_count += 1
        return True

    # First, rescue statistically plausible long-range contacts so calibration does not erase topology.
    for pair, score, lr in scored:
        if pair[1] - pair[0] < POTENTIAL_SEP_FLOOR:
            continue
        if score < probability_threshold * 0.72 and lr < 0.50:
            continue
        try_add(pair, score)
        if long_count >= target_long:
            break
    # Then fill with globally plausible contacts.
    for pair, score, _lr in scored:
        sep = pair[1] - pair[0]
        gate = probability_threshold * (0.86 if sep >= LONG_RANGE_THRESHOLD else 1.0)
        if score < gate:
            continue
        try_add(pair, score)
        if len(selected) >= max_contacts:
            break
    # If too sparse, add top-ranked contacts without relaxing degree cap too much.
    for pair, score, _lr in scored:
        if len(selected) >= max_contacts:
            break
        if score < probability_threshold * 0.62:
            break
        try_add(pair, score)
    return normalized_contact_pairs(selected), scores


def long_range_global_score(
    row: RealCoordinateVisualRow,
    coords: Sequence[Sequence[float]],
    contacts: Sequence[ContactPair],
    restraints: Sequence[CouplingConstraint],
    potential: LongRangePotential,
) -> tuple[float, dict[str, float]]:
    base, components = calibrated_global_score(row, coords, contacts, restraints)
    if contacts:
        lr_contacts = [p for p in contacts if p[1] - p[0] >= POTENTIAL_SEP_FLOOR]
        lr_score = mean([long_range_contact_score(p, row.sequence, potential) for p in lr_contacts]) if lr_contacts else 0.0
        lr_fraction = len(lr_contacts) / len(contacts)
    else:
        lr_score = 0.0
        lr_fraction = 0.0
    # Target moderate long-range presence, not hairball overload.
    lr_balance = exp(-((lr_fraction - 0.34) ** 2) / 0.075)
    score = 0.72 * base + 0.20 * lr_score + 0.08 * lr_balance
    out = dict(components)
    out["long_range_potential_score"] = _rounded(lr_score)
    out["long_range_contact_fraction"] = _rounded(lr_fraction)
    out["long_range_balance_score"] = _rounded(lr_balance)
    return _rounded(score), out


@dataclass(frozen=True)
class LongRangeCalibratedRowReport:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    optional_backend_status: dict[str, object]
    long_range_potential: dict[str, object]
    training_row_count: int
    training_sample_count: int
    training_positive_count: int
    target_native_excluded_from_training: bool
    target_native_excluded_from_long_range_potential: bool
    pdb_long_range_native_contacts_used_from_other_rows: bool
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
    contact_order: float
    contact_order_score: float
    long_range_potential_score: float
    long_range_contact_fraction: float
    long_range_balance_score: float
    raw_extracted_contact_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    final_energy_proxy: float
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_long_range_calibrated_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_target_selection: bool = False
    native_truth_used_before_target_selection: bool = False
    target_native_truth_attached_after_selection_for_evaluation: bool = True
    learned_geometry_prior_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = True
    esmfold_teacher_transfer_used: bool = False

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class LongRangeCalibratedVerifierPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    optional_backend_status: dict[str, object]
    long_range_reference_potential_included: bool
    leave_one_target_out_long_range_calibration: bool
    target_native_excluded_before_selection_for_all_rows: bool
    pdb_long_range_native_contacts_used_from_other_rows: bool
    esmfold_teacher_transfer_used: bool
    universal_physical_law_claim_allowed: bool
    mean_training_sample_count: float
    mean_reference_long_range_base_probability: float
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    long_range_calibrated_verifier_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[LongRangeCalibratedRowReport, ...]
    decisions: tuple[CoarseGrainMDContactDecision, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "optional_backend_status": self.optional_backend_status,
            "long_range_reference_potential_included": self.long_range_reference_potential_included,
            "leave_one_target_out_long_range_calibration": self.leave_one_target_out_long_range_calibration,
            "target_native_excluded_before_selection_for_all_rows": self.target_native_excluded_before_selection_for_all_rows,
            "pdb_long_range_native_contacts_used_from_other_rows": self.pdb_long_range_native_contacts_used_from_other_rows,
            "esmfold_teacher_transfer_used": self.esmfold_teacher_transfer_used,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "mean_training_sample_count": self.mean_training_sample_count,
            "mean_reference_long_range_base_probability": self.mean_reference_long_range_base_probability,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "long_range_calibrated_verifier_claim_allowed": self.long_range_calibrated_verifier_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [r.to_dict() for r in self.rows],
            "decisions": [d.to_dict() for d in self.decisions],
        }


def run_long_range_calibrated_verifier_row(
    *,
    row: RealCoordinateVisualRow,
    all_rows: Sequence[RealCoordinateVisualRow],
    constraints: Sequence[CouplingConstraint],
    md_steps: int = 56,
    candidate_count: int = 4,
    max_direct: int | None = 80,
    max_sequence_closure: int | None = 120,
    max_geodesic: int | None = 260,
    use_chemical_score: bool = False,
) -> tuple[LongRangeCalibratedRowReport, tuple[CoarseGrainMDContactDecision, ...]]:
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
    positives = sum(label for _features, label in samples)
    potential = build_leave_one_out_long_range_potential(all_rows, row.row_id)
    backend = optional_energy_backend_status()
    best = None
    for idx in range(max(1, candidate_count)):
        restraints, coords, raw_contacts, raw_scores, energy = _make_candidate(
            row,
            constraints,
            md_steps=md_steps + 4 * (idx % 3),
            restart_index=idx,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
        )
        selected, selected_scores = long_range_aware_filter_contacts(
            row,
            coords,
            raw_contacts,
            raw_scores,
            restraints,
            weights,
            potential,
            use_chemical_score=use_chemical_score,
        )
        score, components = long_range_global_score(row, coords, selected, restraints, potential)
        if best is None or score > best[0]:
            best = (score, idx, restraints, coords, raw_contacts, selected, selected_scores, components, energy)
    assert best is not None
    selected_score, selected_idx, restraints, _coords, raw_contacts, selected, selected_scores, components, energy = best
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(row=row, selected_pairs=selected)
    f1_margin = _rounded(metric.contact_map_f1 - best_control_f1)
    lr_margin = _rounded(metric.long_range_contact_recall - best_control_lr)
    row_claim, row_universal, reason = _claim(metric, f1_margin, lr_margin)
    report = LongRangeCalibratedRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        optional_backend_status=backend,
        long_range_potential=potential.to_dict(),
        training_row_count=max(0, len(all_rows) - 1),
        training_sample_count=len(samples),
        training_positive_count=positives,
        target_native_excluded_from_training=True,
        target_native_excluded_from_long_range_potential=True,
        pdb_long_range_native_contacts_used_from_other_rows=True,
        safe_restraint_count=len(restraints),
        candidate_structure_count=max(1, candidate_count),
        selected_candidate_index=selected_idx,
        selected_global_score=selected_score,
        statistical_potential_score=components["statistical_potential_score"],
        degree_coherence_score=components["degree_coherence_score"],
        loop_patch_coherence_score=components["loop_patch_coherence_score"],
        restraint_satisfaction_score=components["restraint_satisfaction_score"],
        compactness_score=components["compactness_score"],
        collision_penalty=components["collision_penalty"],
        contact_order=components["contact_order"],
        contact_order_score=components["contact_order_score"],
        long_range_potential_score=components["long_range_potential_score"],
        long_range_contact_fraction=components["long_range_contact_fraction"],
        long_range_balance_score=components["long_range_balance_score"],
        raw_extracted_contact_count=len(raw_contacts),
        selected_contact_count=len(selected),
        selected_long_range_contact_count=sum(1 for i, j in selected if j - i >= LONG_RANGE_THRESHOLD),
        selected_contact_map_hash=_pair_hash(selected),
        final_energy_proxy=_rounded(energy),
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_long_range_calibrated_claim_allowed=row_claim,
        row_universal_physical_law_claim_allowed=row_universal,
        row_claim_rejection_reason=reason,
    )
    decisions = []
    for decision in _decisions(row, selected, selected_scores):
        decisions.append(CoarseGrainMDContactDecision(
            kind=LONG_RANGE_DECISION_KIND,
            row_id=decision.row_id,
            source_accession=decision.source_accession,
            source_mode=LONG_RANGE_CALIBRATED_MODE,
            i=decision.i,
            j=decision.j,
            sequence_separation=decision.sequence_separation,
            final_distance_angstrom=decision.final_distance_angstrom,
            geometry_contact_score=decision.geometry_contact_score,
            selected=decision.selected,
            selected_from_final_structure=True,
            msa_dca_used_before_selection=True,
        ))
    return report, tuple(decisions)


def run_long_range_calibrated_verifier_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint],
    evaluation_source_accessions: Sequence[str] | None = None,
    md_steps: int = 56,
    candidate_count: int = 4,
    max_direct: int | None = 80,
    max_sequence_closure: int | None = 120,
    max_geodesic: int | None = 260,
    use_chemical_score: bool = False,
) -> LongRangeCalibratedVerifierPacket:
    evaluation_set = set(evaluation_source_accessions or ())
    evaluation_rows = tuple(row for row in rows if not evaluation_set or row.source_accession in evaluation_set)
    reports: list[LongRangeCalibratedRowReport] = []
    decisions: list[CoarseGrainMDContactDecision] = []
    for row in evaluation_rows:
        report, row_decisions = run_long_range_calibrated_verifier_row(
            row=row,
            all_rows=rows,
            constraints=constraints,
            md_steps=md_steps,
            candidate_count=candidate_count,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
            use_chemical_score=use_chemical_score,
        )
        reports.append(report)
        decisions.extend(row_decisions[:500])
    precisions = [r.metric_after_native_audit.native_contact_precision for r in reports]
    recalls = [r.metric_after_native_audit.native_contact_recall for r in reports]
    long_ranges = [r.metric_after_native_audit.long_range_contact_recall for r in reports]
    f1s = [r.metric_after_native_audit.contact_map_f1 for r in reports]
    f1_margins = [r.f1_margin_vs_best_control for r in reports]
    lr_margins = [r.long_range_recall_margin_vs_best_control for r in reports]
    claim = bool(reports) and all(r.row_long_range_calibrated_claim_allowed for r in reports) and _mean(precisions) >= 0.70 and _mean(recalls) >= 0.70 and _mean(long_ranges) >= 0.70
    reason = "long_range_calibrated_verifier_survived_gate_not_universal_physics_due_to_pdb_calibration" if claim else "long_range_calibrated_claim_rejected_for_rows:" + ",".join(r.source_accession for r in reports if not r.row_long_range_calibrated_claim_allowed)[:180]
    return LongRangeCalibratedVerifierPacket(
        kind=LONG_RANGE_CALIBRATED_KIND,
        source_mode=LONG_RANGE_CALIBRATED_MODE,
        row_count=len(reports),
        decision_rule=LONG_RANGE_RULE,
        claim_rule=LONG_RANGE_CLAIM_RULE,
        optional_backend_status=optional_energy_backend_status(),
        long_range_reference_potential_included=True,
        leave_one_target_out_long_range_calibration=True,
        target_native_excluded_before_selection_for_all_rows=True,
        pdb_long_range_native_contacts_used_from_other_rows=True,
        esmfold_teacher_transfer_used=False,
        universal_physical_law_claim_allowed=False,
        mean_training_sample_count=_mean([r.training_sample_count for r in reports]),
        mean_reference_long_range_base_probability=_mean([float(r.long_range_potential.get("base_probability", 0.0)) for r in reports]),
        mean_native_contact_precision_after_audit=_mean(precisions),
        mean_native_contact_recall_after_audit=_mean(recalls),
        mean_long_range_contact_recall_after_audit=_mean(long_ranges),
        mean_contact_map_f1_after_audit=_mean(f1s),
        min_row_precision_after_audit=_rounded(min(precisions)) if precisions else 0.0,
        min_row_recall_after_audit=_rounded(min(recalls)) if recalls else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_ranges)) if long_ranges else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        long_range_calibrated_verifier_claim_allowed=claim,
        folding_problem_solved=claim,
        claim_rejection_reason=reason,
        rows=tuple(reports),
        decisions=tuple(decisions),
    )
