from __future__ import annotations

"""Global contact factor-graph + ensemble selector.

This layer addresses the failure mode found by the five-axis challenge: the
physics terms existed, but contact selection was still local/pairwise.  The
module therefore builds explicit all-or-none cooperative contact factors,
residue-degree mutex constraints, and a deterministic top-K ensemble of global
solutions.

Two regimes are intentionally separated:

* sequence-only factor graph: uses only the five-axis sequence scores and keeps
  the universal claim gate strict.
* learned-geometry factor graph: can consume an external non-AlphaFold,
  no-MSA/no-template predicted PDB (for example ESMFold/OmegaFold/SPIRED) as an
  independent global geometry prior.  This can solve a target, but it is not a
  sequence-only physical-law claim.

Native/coordinate truth is only attached after selection for audit metrics.
"""

import hashlib
from collections import defaultdict
from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from statistics import mean
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_five_axis_physics import (
    FiveAxisContactDecision,
    build_five_axis_contact_decisions,
    matched_control_pairs,
)
from pharmacotopology.folding_independent_contact_evidence import (
    IndependentContactEvidencePair,
    contact_evidence_from_predicted_pdb,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    contact_map_hash,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_single_sequence_structure_source import is_alphafold_like_source_id


GLOBAL_FACTOR_GRAPH_KIND = "global_contact_factor_graph_ensemble_v0"
GLOBAL_FACTOR_GRAPH_DECISION_KIND = "global_factor_graph_contact_decision_v0"
GLOBAL_FACTOR_GRAPH_FACTOR_KIND = "global_factor_graph_all_or_none_factor_v0"
GLOBAL_FACTOR_GRAPH_MUTEX_KIND = "global_factor_graph_residue_degree_mutex_v0"
GLOBAL_FACTOR_GRAPH_SOLUTION_KIND = "global_factor_graph_top_k_solution_v0"
SEQUENCE_ONLY_MODE = "sequence_only_five_axis_factor_graph"
LEARNED_GEOMETRY_MODE = "learned_geometry_prior_factor_graph"
GLOBAL_FACTOR_GRAPH_RULE = (
    "build_contact_factors_from_score_neighbourhoods;select_all_or_none_groups;"
    "enforce_residue_degree_mutex;emit_top_k_global_solutions;"
    "consensus_probability_threshold;native_audit_after_selection_only"
)
SEQUENCE_ONLY_CLAIM_RULE = (
    "sequence_only_global_claim_requires_mean_precision_ge_0_70_and_mean_recall_ge_0_70_"
    "and_mean_long_range_recall_ge_0_70_and_no_learned_prior"
)
LEARNED_GEOMETRY_CLAIM_RULE = (
    "target_solved_claim_allowed_only_for_non_alphafold_no_msa_no_template_external_geometry_"
    "with_precision_ge_0_70_and_recall_ge_0_70;not_a_universal_physical_law_claim"
)
DISALLOWED_MSA_TEMPLATE_TOKENS = (
    "jackhmmer",
    "hhblits",
    "hhsearch",
    "hmmer",
    "template",
    "pdb_template",
    "msa_file",
    "a3m",
    "stockholm",
    ".sto",
)


@dataclass(frozen=True)
class GlobalFactorGraphModelSpec:
    source_id: str
    pdb_path: str
    chain_id: str | None = None
    allow_alphafold_source: bool = False
    allow_msa_or_template_source: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class GlobalContactFactor:
    kind: str
    factor_id: str
    row_id: str
    source_accession: str
    source_mode: str
    factor_type: str
    factor_source: str
    contacts: tuple[ContactPair, ...]
    contact_count: int
    long_range_contact_count: int
    mean_pair_score: float
    min_pair_score: float
    max_pair_score: float
    objective_score: float
    anchor_i: int
    anchor_j: int
    all_or_none: bool = True
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["contacts"] = [list(pair) for pair in self.contacts]
        return payload


@dataclass(frozen=True)
class GlobalMutexGroup:
    kind: str
    mutex_id: str
    row_id: str
    source_accession: str
    residue_index: int
    max_selected_degree: int
    candidate_contact_count: int
    candidate_contacts_hash: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class GlobalFactorGraphSolution:
    kind: str
    row_id: str
    source_accession: str
    solution_index: int
    selected_factor_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    objective_score: float
    selected_contact_map_hash: str
    skipped_factor_count_due_to_budget: int
    skipped_factor_count_due_to_mutex: int
    selected_factor_ids: tuple[str, ...]
    metric_after_native_audit: ContactMetricPacket
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def selected_pairs(self) -> tuple[ContactPair, ...]:
        raise RuntimeError("selected pairs are stored on the row-level consensus, not in this compact solution record")

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["selected_factor_ids"] = list(self.selected_factor_ids)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class GlobalFactorGraphContactDecision:
    kind: str
    row_id: str
    source_accession: str
    source_mode: str
    i: int
    j: int
    sequence_separation: int
    base_pair_score: float
    factor_vote_count: int
    selected_solution_count: int
    selection_probability: float
    selected: bool
    selection_reason: str
    supporting_factor_ids: tuple[str, ...]
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["supporting_factor_ids"] = list(self.supporting_factor_ids)
        return payload


@dataclass(frozen=True)
class GlobalFactorGraphRowReport:
    row_id: str
    source_accession: str
    source_mode: str
    sequence_hash: str
    sequence_length: int
    candidate_pair_count: int
    all_or_none_factor_count: int
    mutex_group_count: int
    ensemble_solution_count: int
    consensus_probability_threshold: float
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    metric_after_native_audit: ContactMetricPacket
    best_solution_metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_factor_graph_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    learned_geometry_prior_used_before_selection: bool = False
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "row_id": self.row_id,
            "source_accession": self.source_accession,
            "source_mode": self.source_mode,
            "sequence_hash": self.sequence_hash,
            "sequence_length": self.sequence_length,
            "candidate_pair_count": self.candidate_pair_count,
            "all_or_none_factor_count": self.all_or_none_factor_count,
            "mutex_group_count": self.mutex_group_count,
            "ensemble_solution_count": self.ensemble_solution_count,
            "consensus_probability_threshold": self.consensus_probability_threshold,
            "selected_contact_count": self.selected_contact_count,
            "selected_long_range_contact_count": self.selected_long_range_contact_count,
            "selected_contact_map_hash": self.selected_contact_map_hash,
            "metric_after_native_audit": self.metric_after_native_audit.to_dict(),
            "best_solution_metric_after_native_audit": self.best_solution_metric_after_native_audit.to_dict(),
            "matched_control_count": self.matched_control_count,
            "best_control_f1_after_audit": self.best_control_f1_after_audit,
            "best_control_long_range_recall_after_audit": self.best_control_long_range_recall_after_audit,
            "f1_margin_vs_best_control": self.f1_margin_vs_best_control,
            "long_range_recall_margin_vs_best_control": self.long_range_recall_margin_vs_best_control,
            "row_factor_graph_claim_allowed": self.row_factor_graph_claim_allowed,
            "row_universal_physical_law_claim_allowed": self.row_universal_physical_law_claim_allowed,
            "row_claim_rejection_reason": self.row_claim_rejection_reason,
            "learned_geometry_prior_used_before_selection": self.learned_geometry_prior_used_before_selection,
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "native_truth_attached_after_selection_for_evaluation": self.native_truth_attached_after_selection_for_evaluation,
            "alphafold_used_before_selection": self.alphafold_used_before_selection,
            "msa_used_before_selection": self.msa_used_before_selection,
            "template_used_before_selection": self.template_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


@dataclass(frozen=True)
class GlobalFactorGraphPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    factor_graph_included: bool
    all_or_none_factors_included: bool
    mutex_constraints_included: bool
    top_k_ensemble_included: bool
    learned_geometry_prior_used_before_selection: bool
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    factor_graph_ensemble_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[GlobalFactorGraphRowReport, ...]
    factors: tuple[GlobalContactFactor, ...]
    mutex_groups: tuple[GlobalMutexGroup, ...]
    solutions: tuple[GlobalFactorGraphSolution, ...]
    decisions: tuple[GlobalFactorGraphContactDecision, ...]
    model_specs: tuple[GlobalFactorGraphModelSpec, ...] = ()
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "factor_graph_included": self.factor_graph_included,
            "all_or_none_factors_included": self.all_or_none_factors_included,
            "mutex_constraints_included": self.mutex_constraints_included,
            "top_k_ensemble_included": self.top_k_ensemble_included,
            "learned_geometry_prior_used_before_selection": self.learned_geometry_prior_used_before_selection,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "factor_graph_ensemble_claim_allowed": self.factor_graph_ensemble_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [row.to_dict() for row in self.rows],
            "factors": [factor.to_dict() for factor in self.factors],
            "mutex_groups": [group.to_dict() for group in self.mutex_groups],
            "solutions": [solution.to_dict() for solution in self.solutions],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "model_specs": [spec.to_dict() for spec in self.model_specs],
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "alphafold_used_before_selection": self.alphafold_used_before_selection,
            "msa_used_before_selection": self.msa_used_before_selection,
            "template_used_before_selection": self.template_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, float(value)))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _pair_hash(pairs: Iterable[ContactPair]) -> str:
    return contact_map_hash(normalized_contact_pairs(pairs))


def _stable_float(*parts: object) -> float:
    encoded = ":".join(str(part) for part in parts).encode("utf-8")
    value = int(hashlib.sha256(encoded).hexdigest()[:12], 16)
    return value / float(0xFFFFFFFFFFFF)


def _looks_msa_or_template_like(source_id: str, pdb_path: str) -> bool:
    haystack = f"{source_id} {pdb_path}".lower().replace("\\", "/")
    return any(token in haystack for token in DISALLOWED_MSA_TEMPLATE_TOKENS)


def _candidate_pairs_from_decisions(decisions: Sequence[FiveAxisContactDecision]) -> tuple[ContactPair, ...]:
    return normalized_contact_pairs(decision.pair() for decision in decisions)


def _five_axis_score_map(decisions: Sequence[FiveAxisContactDecision]) -> dict[ContactPair, float]:
    return {decision.pair(): float(decision.final_score) for decision in decisions}


def _evidence_score_map(evidence: Sequence[IndependentContactEvidencePair]) -> dict[ContactPair, float]:
    grouped: dict[ContactPair, list[float]] = defaultdict(list)
    for item in evidence:
        distance_score = 1.0
        if item.distance_angstrom is not None:
            # 8A is the contact cutoff.  Give shorter predicted contacts a higher
            # score while keeping every predicted contact positive.
            distance_score = _clamp(1.0 - (float(item.distance_angstrom) - 3.5) / 6.0, 0.25, 1.0)
        grouped[item.pair()].append(0.74 * float(item.confidence) + 0.26 * distance_score)
    return {pair: _rounded(mean(values)) for pair, values in grouped.items()}


def _merge_scores_for_learned_geometry(
    *,
    five_axis_scores: Mapping[ContactPair, float],
    evidence_scores: Mapping[ContactPair, float],
) -> dict[ContactPair, float]:
    merged: dict[ContactPair, float] = {}
    for pair, evidence_score in evidence_scores.items():
        physics_score = float(five_axis_scores.get(pair, 0.50))
        merged[pair] = _rounded(0.76 * evidence_score + 0.24 * physics_score)
    return merged


def _local_factor_members(
    *,
    anchor: ContactPair,
    pair_pool: Sequence[ContactPair],
    scores: Mapping[ContactPair, float],
    neighbourhood_radius: int,
    min_relative_score: float,
    source_mode: str,
) -> tuple[ContactPair, ...]:
    ai, aj = anchor
    anchor_score = float(scores.get(anchor, 0.0))
    members: list[ContactPair] = []
    for pair in pair_pool:
        pi, pj = pair
        if abs(pi - ai) > neighbourhood_radius or abs(pj - aj) > neighbourhood_radius:
            continue
        if source_mode == SEQUENCE_ONLY_MODE and float(scores.get(pair, 0.0)) < anchor_score * min_relative_score:
            continue
        members.append(pair)
    if anchor not in members:
        members.append(anchor)
    members = sorted(set(members), key=lambda pair: (-scores.get(pair, 0.0), pair[0], pair[1]))
    # Keep factors compact so all-or-none really means a local cooperative patch,
    # not one huge threshold blob.
    return normalized_contact_pairs(members[:18])


def build_contact_factors(
    *,
    row: RealCoordinateVisualRow,
    source_mode: str,
    scored_pairs: Mapping[ContactPair, float],
    pair_pool: Sequence[ContactPair],
    factor_source: str,
    neighbourhood_radius: int = 2,
    min_relative_score: float = 0.86,
    max_factor_count: int | None = None,
    learned_geometry_prior_used_before_selection: bool = False,
) -> tuple[GlobalContactFactor, ...]:
    ordered = sorted(
        set(pair_pool),
        key=lambda pair: (-float(scored_pairs.get(pair, 0.0)), pair[0], pair[1]),
    )
    if max_factor_count is None:
        max_factor_count = max(32, min(220, row.sequence_length * 2))
    factors: list[GlobalContactFactor] = []
    covered_anchors: set[ContactPair] = set()
    for anchor in ordered:
        if len(factors) >= max_factor_count:
            break
        if anchor in covered_anchors:
            continue
        anchor_score = float(scored_pairs.get(anchor, 0.0))
        if anchor_score <= 0.0:
            continue
        members = _local_factor_members(
            anchor=anchor,
            pair_pool=ordered,
            scores=scored_pairs,
            neighbourhood_radius=neighbourhood_radius,
            min_relative_score=min_relative_score,
            source_mode=source_mode,
        )
        if not members:
            continue
        scores = [float(scored_pairs.get(pair, 0.0)) for pair in members]
        long_count = sum(1 for pair in members if pair[1] - pair[0] >= 24)
        mean_score = mean(scores)
        min_score = min(scores)
        max_score = max(scores)
        objective = _score(sum(scores) + 0.035 * sqrt(len(members)) + 0.018 * long_count)
        factor_id = (
            f"{row.row_id}:factor:{source_mode}:"
            f"{len(factors) + 1:04d}:{anchor[0]}_{anchor[1]}:{_pair_hash(members)[:12]}"
        )
        factors.append(
            GlobalContactFactor(
                kind=GLOBAL_FACTOR_GRAPH_FACTOR_KIND,
                factor_id=factor_id,
                row_id=row.row_id,
                source_accession=row.source_accession,
                source_mode=source_mode,
                factor_type="all_or_none",
                factor_source=factor_source,
                contacts=members,
                contact_count=len(members),
                long_range_contact_count=long_count,
                mean_pair_score=_rounded(mean_score),
                min_pair_score=_rounded(min_score),
                max_pair_score=_rounded(max_score),
                objective_score=objective,
                anchor_i=anchor[0],
                anchor_j=anchor[1],
                learned_geometry_prior_used_before_selection=learned_geometry_prior_used_before_selection,
            )
        )
        # Avoid many duplicate factors around the exact same cooperative patch.
        for member in members:
            covered_anchors.add(member)
    return tuple(factors)


def build_residue_mutex_groups(
    *,
    row: RealCoordinateVisualRow,
    pair_pool: Sequence[ContactPair],
    max_selected_degree: int,
) -> tuple[GlobalMutexGroup, ...]:
    by_residue: dict[int, list[ContactPair]] = defaultdict(list)
    for pair in set(pair_pool):
        by_residue[pair[0]].append(pair)
        by_residue[pair[1]].append(pair)
    groups: list[GlobalMutexGroup] = []
    for residue_index, contacts in sorted(by_residue.items()):
        normalized = normalized_contact_pairs(contacts)
        if len(normalized) <= max_selected_degree:
            continue
        groups.append(
            GlobalMutexGroup(
                kind=GLOBAL_FACTOR_GRAPH_MUTEX_KIND,
                mutex_id=f"{row.row_id}:residue_mutex:{residue_index:04d}",
                row_id=row.row_id,
                source_accession=row.source_accession,
                residue_index=residue_index,
                max_selected_degree=max_selected_degree,
                candidate_contact_count=len(normalized),
                candidate_contacts_hash=_pair_hash(normalized),
            )
        )
    return tuple(groups)


def _factor_order_for_solution(
    factors: Sequence[GlobalContactFactor],
    *,
    solution_index: int,
    diversity_temperature: float,
) -> list[GlobalContactFactor]:
    return sorted(
        factors,
        key=lambda factor: (
            -(
                factor.objective_score
                + diversity_temperature
                * (_stable_float(factor.factor_id, "solution", solution_index) - 0.5)
            ),
            factor.factor_id,
        ),
    )


def _solve_one_global_solution(
    *,
    row: RealCoordinateVisualRow,
    factors: Sequence[GlobalContactFactor],
    solution_index: int,
    max_selected_contacts: int,
    max_residue_degree: int,
    diversity_temperature: float,
) -> tuple[tuple[GlobalContactFactor, ...], tuple[ContactPair, ...], float, int, int]:
    selected_factors: list[GlobalContactFactor] = []
    selected_pairs: set[ContactPair] = set()
    residue_degree: dict[int, int] = defaultdict(int)
    skipped_budget = 0
    skipped_mutex = 0
    objective = 0.0
    for factor in _factor_order_for_solution(
        factors,
        solution_index=solution_index,
        diversity_temperature=diversity_temperature,
    ):
        new_pairs = [pair for pair in factor.contacts if pair not in selected_pairs]
        if not new_pairs:
            continue
        if len(selected_pairs) + len(new_pairs) > max_selected_contacts:
            skipped_budget += 1
            continue
        proposed_degree = dict(residue_degree)
        violates = False
        for pair in new_pairs:
            proposed_degree[pair[0]] = proposed_degree.get(pair[0], 0) + 1
            proposed_degree[pair[1]] = proposed_degree.get(pair[1], 0) + 1
            if proposed_degree[pair[0]] > max_residue_degree or proposed_degree[pair[1]] > max_residue_degree:
                violates = True
                break
        if violates:
            skipped_mutex += 1
            continue
        selected_factors.append(factor)
        objective += factor.objective_score
        for pair in new_pairs:
            selected_pairs.add(pair)
            residue_degree[pair[0]] += 1
            residue_degree[pair[1]] += 1
    return (
        tuple(selected_factors),
        normalized_contact_pairs(selected_pairs),
        _score(objective),
        skipped_budget,
        skipped_mutex,
    )


def solve_top_k_factor_graph(
    *,
    row: RealCoordinateVisualRow,
    factors: Sequence[GlobalContactFactor],
    source_mode: str,
    ensemble_size: int = 8,
    max_selected_contacts: int | None = None,
    max_residue_degree: int | None = None,
    diversity_temperature: float | None = None,
) -> tuple[tuple[GlobalFactorGraphSolution, ...], tuple[tuple[ContactPair, ...], ...]]:
    max_contacts = max_selected_contacts or max(64, int(round(row.sequence_length * 3.05)))
    if max_residue_degree is None:
        max_residue_degree = 24 if source_mode == LEARNED_GEOMETRY_MODE else 14
    if diversity_temperature is None:
        diversity_temperature = 0.004 if source_mode == LEARNED_GEOMETRY_MODE else 0.045
    solutions: list[GlobalFactorGraphSolution] = []
    solution_pairs: list[tuple[ContactPair, ...]] = []
    for solution_index in range(1, max(1, ensemble_size) + 1):
        selected_factors, pairs, objective, skipped_budget, skipped_mutex = _solve_one_global_solution(
            row=row,
            factors=factors,
            solution_index=solution_index,
            max_selected_contacts=max_contacts,
            max_residue_degree=max_residue_degree,
            diversity_temperature=diversity_temperature,
        )
        metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=pairs)
        solution_pairs.append(pairs)
        solutions.append(
            GlobalFactorGraphSolution(
                kind=GLOBAL_FACTOR_GRAPH_SOLUTION_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                solution_index=solution_index,
                selected_factor_count=len(selected_factors),
                selected_contact_count=len(pairs),
                selected_long_range_contact_count=sum(1 for pair in pairs if pair[1] - pair[0] >= 24),
                objective_score=objective,
                selected_contact_map_hash=_pair_hash(pairs),
                skipped_factor_count_due_to_budget=skipped_budget,
                skipped_factor_count_due_to_mutex=skipped_mutex,
                selected_factor_ids=tuple(factor.factor_id for factor in selected_factors),
                metric_after_native_audit=metric,
                learned_geometry_prior_used_before_selection=source_mode == LEARNED_GEOMETRY_MODE,
            )
        )
    return tuple(solutions), tuple(solution_pairs)


def _consensus_pairs(
    solution_pairs: Sequence[Sequence[ContactPair]],
    *,
    threshold: float,
) -> tuple[ContactPair, ...]:
    counts: dict[ContactPair, int] = defaultdict(int)
    for pairs in solution_pairs:
        for pair in set(pairs):
            counts[pair] += 1
    denom = max(1, len(solution_pairs))
    return normalized_contact_pairs(pair for pair, count in counts.items() if count / denom >= threshold)


def _contact_decisions(
    *,
    row: RealCoordinateVisualRow,
    source_mode: str,
    base_scores: Mapping[ContactPair, float],
    factors: Sequence[GlobalContactFactor],
    solution_pairs: Sequence[Sequence[ContactPair]],
    consensus_pairs: Sequence[ContactPair],
    consensus_probability_threshold: float,
    max_report_decisions: int = 1200,
) -> tuple[GlobalFactorGraphContactDecision, ...]:
    pair_to_factor_ids: dict[ContactPair, list[str]] = defaultdict(list)
    for factor in factors:
        for pair in factor.contacts:
            pair_to_factor_ids[pair].append(factor.factor_id)
    selected_set = set(consensus_pairs)
    solution_counts: dict[ContactPair, int] = defaultdict(int)
    for pairs in solution_pairs:
        for pair in set(pairs):
            solution_counts[pair] += 1
    denom = max(1, len(solution_pairs))
    all_pairs = set(pair_to_factor_ids) | selected_set
    ordered_pairs = sorted(
        all_pairs,
        key=lambda pair: (
            pair not in selected_set,
            -solution_counts.get(pair, 0) / denom,
            -float(base_scores.get(pair, 0.0)),
            pair[0],
            pair[1],
        ),
    )[:max_report_decisions]
    decisions: list[GlobalFactorGraphContactDecision] = []
    for pair in ordered_pairs:
        probability = _rounded(solution_counts.get(pair, 0) / denom)
        selected = pair in selected_set
        reason = (
            f"selected_by_top_k_consensus_probability_ge_{consensus_probability_threshold:.2f}"
            if selected
            else "below_top_k_consensus_probability_threshold"
        )
        decisions.append(
            GlobalFactorGraphContactDecision(
                kind=GLOBAL_FACTOR_GRAPH_DECISION_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                source_mode=source_mode,
                i=pair[0],
                j=pair[1],
                sequence_separation=pair[1] - pair[0],
                base_pair_score=_rounded(float(base_scores.get(pair, 0.0))),
                factor_vote_count=len(pair_to_factor_ids.get(pair, ())),
                selected_solution_count=solution_counts.get(pair, 0),
                selection_probability=probability,
                selected=selected,
                selection_reason=reason,
                supporting_factor_ids=tuple(pair_to_factor_ids.get(pair, ())[:8]),
                learned_geometry_prior_used_before_selection=source_mode == LEARNED_GEOMETRY_MODE,
            )
        )
    return tuple(decisions)


def _matched_controls_for_report(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[ContactPair],
    candidate_pairs: Sequence[ContactPair],
) -> tuple[float, float, int]:
    controls = []
    for index in range(1, max(2, min(6, int(sqrt(max(1, len(selected_pairs))))) ) + 1):
        control_pairs = matched_control_pairs(
            row=row,
            selected_pairs=selected_pairs,
            candidate_pairs=candidate_pairs,
            control_index=index,
        )
        if not control_pairs:
            continue
        controls.append(evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=control_pairs))
    best_f1 = max((metric.contact_map_f1 for metric in controls), default=0.0)
    best_lr = max((metric.long_range_contact_recall for metric in controls), default=0.0)
    return best_f1, best_lr, len(controls)


def _empty_metric(row: RealCoordinateVisualRow) -> ContactMetricPacket:
    return evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=())


def _load_learned_geometry_evidence(
    *,
    row: RealCoordinateVisualRow,
    model_specs: Sequence[GlobalFactorGraphModelSpec],
) -> tuple[tuple[IndependentContactEvidencePair, ...], bool, bool, bool, tuple[GlobalFactorGraphModelSpec, ...]]:
    evidence: list[IndependentContactEvidencePair] = []
    alphafold_used = False
    msa_used = False
    template_used = False
    accepted_specs: list[GlobalFactorGraphModelSpec] = []
    for spec in model_specs:
        alphafold_like = is_alphafold_like_source_id(f"{spec.source_id} {spec.pdb_path}")
        msa_template_like = _looks_msa_or_template_like(spec.source_id, spec.pdb_path)
        if alphafold_like and not spec.allow_alphafold_source:
            alphafold_used = True
            continue
        if msa_template_like and not spec.allow_msa_or_template_source:
            msa_used = True
            template_used = True
            continue
        try:
            source_evidence = contact_evidence_from_predicted_pdb(
                row=row,
                pdb_path=Path(spec.pdb_path),
                source_id=spec.source_id,
                source_family="learned_geometry_prior",
                source_kind="external_non_alphafold_predicted_structure_contacts_v0",
                chain_id=spec.chain_id,
            )
        except (OSError, ValueError):
            continue
        evidence.extend(source_evidence)
        accepted_specs.append(spec)
    return tuple(evidence), alphafold_used, msa_used, template_used, tuple(accepted_specs)


def run_global_factor_graph_row(
    *,
    row: RealCoordinateVisualRow,
    source_mode: str = SEQUENCE_ONLY_MODE,
    model_specs: Sequence[GlobalFactorGraphModelSpec] = (),
    max_sequence_separation: int | None = 160,
    ensemble_size: int = 8,
    consensus_probability_threshold: float | None = None,
    max_selected_contacts: int | None = None,
) -> tuple[
    GlobalFactorGraphRowReport,
    tuple[GlobalContactFactor, ...],
    tuple[GlobalMutexGroup, ...],
    tuple[GlobalFactorGraphSolution, ...],
    tuple[GlobalFactorGraphContactDecision, ...],
    tuple[GlobalFactorGraphModelSpec, ...],
]:
    five_axis_decisions = build_five_axis_contact_decisions(
        row,
        max_sequence_separation=max_sequence_separation if source_mode == SEQUENCE_ONLY_MODE else None,
    )
    five_axis_scores = _five_axis_score_map(five_axis_decisions)
    learned_evidence: tuple[IndependentContactEvidencePair, ...] = ()
    alphafold_used = False
    msa_used = False
    template_used = False
    accepted_specs: tuple[GlobalFactorGraphModelSpec, ...] = ()
    if source_mode == LEARNED_GEOMETRY_MODE:
        learned_evidence, alphafold_used, msa_used, template_used, accepted_specs = _load_learned_geometry_evidence(
            row=row,
            model_specs=model_specs,
        )
        evidence_scores = _evidence_score_map(learned_evidence)
        scored_pairs = _merge_scores_for_learned_geometry(
            five_axis_scores=five_axis_scores,
            evidence_scores=evidence_scores,
        )
        pair_pool = normalized_contact_pairs(scored_pairs)
        factor_source = "external_non_alphafold_learned_geometry_contact_clusters"
        learned_prior = bool(learned_evidence)
        neighbourhood_radius = 2
        min_relative_score = 0.0
        max_factor_count = max(64, row.sequence_length * 3)
        if consensus_probability_threshold is None:
            consensus_probability_threshold = 0.50
    else:
        scored_pairs = five_axis_scores
        pair_pool = _candidate_pairs_from_decisions(five_axis_decisions)
        factor_source = "sequence_only_five_axis_score_clusters"
        learned_prior = False
        neighbourhood_radius = 2
        min_relative_score = 0.86
        max_factor_count = max(40, min(180, row.sequence_length))
        if consensus_probability_threshold is None:
            consensus_probability_threshold = 0.55

    if not pair_pool or (source_mode == LEARNED_GEOMETRY_MODE and not learned_evidence):
        empty = _empty_metric(row)
        report = GlobalFactorGraphRowReport(
            row_id=row.row_id,
            source_accession=row.source_accession,
            source_mode=source_mode,
            sequence_hash=row.sequence_sha256,
            sequence_length=row.sequence_length,
            candidate_pair_count=0,
            all_or_none_factor_count=0,
            mutex_group_count=0,
            ensemble_solution_count=0,
            consensus_probability_threshold=float(consensus_probability_threshold or 0.0),
            selected_contact_count=0,
            selected_long_range_contact_count=0,
            selected_contact_map_hash=_pair_hash(()),
            metric_after_native_audit=empty,
            best_solution_metric_after_native_audit=empty,
            matched_control_count=0,
            best_control_f1_after_audit=0.0,
            best_control_long_range_recall_after_audit=0.0,
            f1_margin_vs_best_control=0.0,
            long_range_recall_margin_vs_best_control=0.0,
            row_factor_graph_claim_allowed=False,
            row_universal_physical_law_claim_allowed=False,
            row_claim_rejection_reason="missing_usable_contact_source_for_global_factor_graph",
            learned_geometry_prior_used_before_selection=learned_prior,
            alphafold_used_before_selection=alphafold_used,
            msa_used_before_selection=msa_used,
            template_used_before_selection=template_used,
        )
        return report, (), (), (), (), accepted_specs

    factors = build_contact_factors(
        row=row,
        source_mode=source_mode,
        scored_pairs=scored_pairs,
        pair_pool=pair_pool,
        factor_source=factor_source,
        neighbourhood_radius=neighbourhood_radius,
        min_relative_score=min_relative_score,
        max_factor_count=max_factor_count,
        learned_geometry_prior_used_before_selection=learned_prior,
    )
    max_degree = 24 if source_mode == LEARNED_GEOMETRY_MODE else 14
    mutex_groups = build_residue_mutex_groups(
        row=row,
        pair_pool=pair_pool,
        max_selected_degree=max_degree,
    )
    solutions, solution_pairs = solve_top_k_factor_graph(
        row=row,
        factors=factors,
        source_mode=source_mode,
        ensemble_size=ensemble_size,
        max_selected_contacts=max_selected_contacts,
        max_residue_degree=max_degree,
    )
    consensus_pairs = _consensus_pairs(
        solution_pairs,
        threshold=float(consensus_probability_threshold),
    )
    consensus_metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=consensus_pairs,
    )
    best_solution_metric = max(
        (solution.metric_after_native_audit for solution in solutions),
        key=lambda metric: (metric.contact_map_f1, metric.native_contact_recall, metric.native_contact_precision),
        default=_empty_metric(row),
    )
    candidate_for_controls = pair_pool if source_mode == SEQUENCE_ONLY_MODE else _candidate_pairs_from_decisions(five_axis_decisions)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(
        row=row,
        selected_pairs=consensus_pairs,
        candidate_pairs=candidate_for_controls,
    )
    f1_margin = _score(consensus_metric.contact_map_f1 - best_control_f1)
    lr_margin = _score(consensus_metric.long_range_contact_recall - best_control_lr)
    if source_mode == LEARNED_GEOMETRY_MODE:
        row_claim_allowed = (
            learned_prior
            and not alphafold_used
            and not msa_used
            and not template_used
            and consensus_metric.native_contact_precision >= 0.70
            and consensus_metric.native_contact_recall >= 0.70
        )
        universal_allowed = False
        if row_claim_allowed:
            rejection = "target_solved_by_global_factor_graph_over_non_alphafold_learned_geometry_prior"
        elif not learned_prior:
            rejection = "row_claim_rejected_missing_learned_geometry_prior"
        elif consensus_metric.native_contact_precision < 0.70:
            rejection = "row_claim_rejected_precision_below_0_70"
        elif consensus_metric.native_contact_recall < 0.70:
            rejection = "row_claim_rejected_recall_below_0_70"
        else:
            rejection = "row_claim_rejected_global_factor_graph_gate_failed"
    else:
        row_claim_allowed = (
            consensus_metric.native_contact_precision >= 0.70
            and consensus_metric.native_contact_recall >= 0.70
            and consensus_metric.long_range_contact_recall >= 0.70
            and f1_margin >= 0.15
            and lr_margin >= 0.15
        )
        universal_allowed = row_claim_allowed
        if row_claim_allowed:
            rejection = "sequence_only_row_claim_survived_global_factor_graph_gate"
        elif consensus_metric.native_contact_precision < 0.70:
            rejection = "sequence_only_row_claim_rejected_precision_below_0_70"
        elif consensus_metric.native_contact_recall < 0.70:
            rejection = "sequence_only_row_claim_rejected_recall_below_0_70"
        elif consensus_metric.long_range_contact_recall < 0.70:
            rejection = "sequence_only_row_claim_rejected_long_range_recall_below_0_70"
        else:
            rejection = "sequence_only_row_claim_rejected_matched_control_margin_below_gate"

    decisions = _contact_decisions(
        row=row,
        source_mode=source_mode,
        base_scores=scored_pairs,
        factors=factors,
        solution_pairs=solution_pairs,
        consensus_pairs=consensus_pairs,
        consensus_probability_threshold=float(consensus_probability_threshold),
    )
    report = GlobalFactorGraphRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        source_mode=source_mode,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        candidate_pair_count=len(pair_pool),
        all_or_none_factor_count=len(factors),
        mutex_group_count=len(mutex_groups),
        ensemble_solution_count=len(solutions),
        consensus_probability_threshold=float(consensus_probability_threshold),
        selected_contact_count=len(consensus_pairs),
        selected_long_range_contact_count=sum(1 for pair in consensus_pairs if pair[1] - pair[0] >= 24),
        selected_contact_map_hash=_pair_hash(consensus_pairs),
        metric_after_native_audit=consensus_metric,
        best_solution_metric_after_native_audit=best_solution_metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_factor_graph_claim_allowed=row_claim_allowed,
        row_universal_physical_law_claim_allowed=universal_allowed,
        row_claim_rejection_reason=rejection,
        learned_geometry_prior_used_before_selection=learned_prior,
        alphafold_used_before_selection=alphafold_used,
        msa_used_before_selection=msa_used,
        template_used_before_selection=template_used,
    )
    return report, factors, mutex_groups, solutions, decisions, accepted_specs


def run_global_factor_graph_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    source_mode: str = SEQUENCE_ONLY_MODE,
    model_specs_by_accession: Mapping[str, Sequence[GlobalFactorGraphModelSpec]] | None = None,
    max_sequence_separation: int | None = 160,
    ensemble_size: int = 8,
    consensus_probability_threshold: float | None = None,
) -> GlobalFactorGraphPacket:
    row_reports: list[GlobalFactorGraphRowReport] = []
    all_factors: list[GlobalContactFactor] = []
    all_mutex_groups: list[GlobalMutexGroup] = []
    all_solutions: list[GlobalFactorGraphSolution] = []
    all_decisions: list[GlobalFactorGraphContactDecision] = []
    accepted_specs: list[GlobalFactorGraphModelSpec] = []
    specs_by_accession = model_specs_by_accession or {}
    for row in rows:
        specs = tuple(specs_by_accession.get(row.source_accession, ()))
        report, factors, mutex_groups, solutions, decisions, specs_accepted = run_global_factor_graph_row(
            row=row,
            source_mode=source_mode,
            model_specs=specs,
            max_sequence_separation=max_sequence_separation,
            ensemble_size=ensemble_size,
            consensus_probability_threshold=consensus_probability_threshold,
        )
        row_reports.append(report)
        all_factors.extend(factors[:240])
        all_mutex_groups.extend(mutex_groups[:240])
        all_solutions.extend(solutions)
        selected_decisions = [decision for decision in decisions if decision.selected]
        rejected_tail = [decision for decision in decisions if not decision.selected][:80]
        all_decisions.extend(selected_decisions + rejected_tail)
        accepted_specs.extend(specs_accepted)

    precision_values = [row.metric_after_native_audit.native_contact_precision for row in row_reports]
    recall_values = [row.metric_after_native_audit.native_contact_recall for row in row_reports]
    long_recall_values = [row.metric_after_native_audit.long_range_contact_recall for row in row_reports]
    f1_values = [row.metric_after_native_audit.contact_map_f1 for row in row_reports]
    f1_margins = [row.f1_margin_vs_best_control for row in row_reports]
    lr_margins = [row.long_range_recall_margin_vs_best_control for row in row_reports]
    learned_prior = source_mode == LEARNED_GEOMETRY_MODE and any(
        row.learned_geometry_prior_used_before_selection for row in row_reports
    )
    if source_mode == LEARNED_GEOMETRY_MODE:
        factor_claim = bool(row_reports) and all(row.row_factor_graph_claim_allowed for row in row_reports)
        universal_allowed = False
        if factor_claim:
            rejection = "global_factor_graph_ensemble_target_claim_survived_learned_geometry_gate_not_universal_physics"
        else:
            failed = [row.source_accession for row in row_reports if not row.row_factor_graph_claim_allowed]
            rejection = "global_factor_graph_learned_geometry_claim_rejected_for_rows:" + ",".join(failed[:12])
        claim_rule = LEARNED_GEOMETRY_CLAIM_RULE
    else:
        factor_claim = (
            bool(row_reports)
            and all(row.row_factor_graph_claim_allowed for row in row_reports)
            and _mean(precision_values) >= 0.70
            and _mean(recall_values) >= 0.70
            and _mean(long_recall_values) >= 0.70
        )
        universal_allowed = factor_claim and not learned_prior
        if factor_claim:
            rejection = "sequence_only_global_factor_graph_claim_survived_strict_gate"
        else:
            failed = [row.source_accession for row in row_reports if not row.row_factor_graph_claim_allowed]
            rejection = "sequence_only_global_factor_graph_claim_rejected_for_rows:" + ",".join(failed[:12])
        claim_rule = SEQUENCE_ONLY_CLAIM_RULE

    return GlobalFactorGraphPacket(
        kind=GLOBAL_FACTOR_GRAPH_KIND,
        source_mode=source_mode,
        row_count=len(row_reports),
        decision_rule=GLOBAL_FACTOR_GRAPH_RULE,
        claim_rule=claim_rule,
        factor_graph_included=True,
        all_or_none_factors_included=True,
        mutex_constraints_included=True,
        top_k_ensemble_included=True,
        learned_geometry_prior_used_before_selection=learned_prior,
        mean_native_contact_precision_after_audit=_mean(precision_values),
        mean_native_contact_recall_after_audit=_mean(recall_values),
        mean_long_range_contact_recall_after_audit=_mean(long_recall_values),
        mean_contact_map_f1_after_audit=_mean(f1_values),
        min_row_precision_after_audit=_rounded(min(precision_values)) if precision_values else 0.0,
        min_row_recall_after_audit=_rounded(min(recall_values)) if recall_values else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_recall_values)) if long_recall_values else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        factor_graph_ensemble_claim_allowed=factor_claim,
        universal_physical_law_claim_allowed=universal_allowed,
        folding_problem_solved=factor_claim,
        claim_rejection_reason=rejection,
        rows=tuple(row_reports),
        factors=tuple(all_factors),
        mutex_groups=tuple(all_mutex_groups),
        solutions=tuple(all_solutions),
        decisions=tuple(all_decisions),
        model_specs=tuple(accepted_specs),
        coordinate_truth_used_before_selection=False,
        native_truth_used_before_selection=False,
        alphafold_used_before_selection=any(row.alphafold_used_before_selection for row in row_reports),
        msa_used_before_selection=any(row.msa_used_before_selection for row in row_reports),
        template_used_before_selection=any(row.template_used_before_selection for row in row_reports),
        raw_sequence_exposed=False,
    )
