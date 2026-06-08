from __future__ import annotations

"""MSA-free learned-geometry model ensemble and physics refinement.

This module is intentionally a *consumer* of MSA-free structure predictions.  It
never downloads models, never calls AlphaFold, never uses templates/MSAs, and
never reads native coordinates before selecting contacts.  It combines whatever
single-sequence predicted structures are available (ESMFold, OmegaFold, SPIRED,
Chai-1, Boltz, or any compatible PDB output) with sequence-only physical priors
and a small iterative contact-graph refinement pass.
"""

from collections import defaultdict
from dataclasses import asdict, dataclass
from math import log
from pathlib import Path
from statistics import mean
from typing import Iterable, Mapping, Sequence

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
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow
from pharmacotopology.folding_sequence_physical_priors import (
    SequenceContactPhysicalPrior,
    build_sequence_physical_prior_scores,
)
from pharmacotopology.folding_single_sequence_structure_source import is_alphafold_like_source_id


MSA_FREE_MODEL_ENSEMBLE_KIND = "msa_free_learned_global_geometry_ensemble_v0"
MSA_FREE_MODEL_CONSENSUS_KIND = "msa_free_model_consensus_contact_decision_v0"
MSA_FREE_MODEL_ENSEMBLE_REPORT_KIND = "msa_free_model_ensemble_report_v0"

DISALLOWED_PROVENANCE_TOKENS = (
    "alphafold",
    "alpha_fold",
    "afdb",
    "af-",
    "colabfold",
    "jackhmmer",
    "hhblits",
    "hhsearch",
    "hmmer",
    "template",
    "pdb_template",
    "msa_file",
    "a3m",
    "sto",
    "stockholm",
)


@dataclass(frozen=True)
class MSAFreeStructureModelSpec:
    source_id: str
    pdb_path: str
    chain_id: str | None = None
    allow_alphafold_source: bool = False
    allow_msa_or_template_source: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MSAFreeStructureModelContacts:
    source_id: str
    pdb_path: str
    chain_id: str | None
    status: str
    rejected: bool
    parse_error: str
    contact_count: int
    mean_contact_confidence: float
    contact_map_hash: str
    alphafold_like_source: bool
    msa_or_template_like_source: bool
    raw_sequence_exposed: bool = False
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MSAFreeConsensusContactDecision:
    kind: str
    row_id: str
    source_accession: str
    i: int
    j: int
    model_vote_count: int
    model_vote_fraction: float
    model_source_ids: tuple[str, ...]
    mean_model_confidence: float
    physical_prior_score: float
    contact_energy_score: float
    secondary_structure_score: float
    degree_consistency_score: float
    coupling_support_score: float
    candidate_region_support_score: float
    cooperative_context_score: float
    loop_entropy_score: float
    final_score: float
    selected: bool
    selection_reason: str
    iteration_index: int
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["model_source_ids"] = list(self.model_source_ids)
        return payload


@dataclass(frozen=True)
class MSAFreeModelEnsemblePacket:
    kind: str
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    model_count_requested: int
    usable_model_count: int
    rejected_model_count: int
    available_model_contact_count: int
    candidate_pair_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    iteration_count_requested: int
    iteration_count_executed: int
    converged: bool
    selection_threshold: float
    max_selected_contacts: int
    direct_best_model_source_id: str
    direct_best_model_metric: ContactMetricPacket
    consensus_metric: ContactMetricPacket
    folding_problem_solved: bool
    direct_structure_solved: bool
    consensus_contact_collapse_solved: bool
    folding_solution_mode: str
    claim_rejection_reason: str
    models: tuple[MSAFreeStructureModelContacts, ...]
    decisions: tuple[MSAFreeConsensusContactDecision, ...]
    selected_pairs_hash: str
    native_truth_used_before_selection: bool = False
    coordinate_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "row_id": self.row_id,
            "source_accession": self.source_accession,
            "sequence_hash": self.sequence_hash,
            "sequence_length": self.sequence_length,
            "model_count_requested": self.model_count_requested,
            "usable_model_count": self.usable_model_count,
            "rejected_model_count": self.rejected_model_count,
            "available_model_contact_count": self.available_model_contact_count,
            "candidate_pair_count": self.candidate_pair_count,
            "selected_contact_count": self.selected_contact_count,
            "selected_long_range_contact_count": self.selected_long_range_contact_count,
            "iteration_count_requested": self.iteration_count_requested,
            "iteration_count_executed": self.iteration_count_executed,
            "converged": self.converged,
            "selection_threshold": self.selection_threshold,
            "max_selected_contacts": self.max_selected_contacts,
            "direct_best_model_source_id": self.direct_best_model_source_id,
            "direct_best_model_metric": self.direct_best_model_metric.to_dict(),
            "consensus_metric": self.consensus_metric.to_dict(),
            "folding_problem_solved": self.folding_problem_solved,
            "direct_structure_solved": self.direct_structure_solved,
            "consensus_contact_collapse_solved": self.consensus_contact_collapse_solved,
            "folding_solution_mode": self.folding_solution_mode,
            "claim_rejection_reason": self.claim_rejection_reason,
            "models": [model.to_dict() for model in self.models],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "selected_pairs_hash": self.selected_pairs_hash,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "alphafold_used_before_selection": self.alphafold_used_before_selection,
            "msa_used_before_selection": self.msa_used_before_selection,
            "template_used_before_selection": self.template_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def _empty_metric(row: RealCoordinateVisualRow) -> ContactMetricPacket:
    return evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=())


def _looks_msa_or_template_like(source_id: str, path: str) -> bool:
    normalized = f"{source_id} {path}".lower().replace("/", "_").replace("-", "-")
    return any(token in normalized for token in DISALLOWED_PROVENANCE_TOKENS if token not in {"alphafold", "alpha_fold", "afdb", "af-"})


def _mean_confidence(evidence: Sequence[IndependentContactEvidencePair]) -> float:
    if not evidence:
        return 0.0
    return _rounded(mean(item.confidence for item in evidence))


def _load_model_contacts(
    *,
    row: RealCoordinateVisualRow,
    spec: MSAFreeStructureModelSpec,
) -> tuple[MSAFreeStructureModelContacts, tuple[IndependentContactEvidencePair, ...]]:
    source_text = f"{spec.source_id} {spec.pdb_path}"
    alphafold_like = is_alphafold_like_source_id(source_text)
    msa_like = _looks_msa_or_template_like(spec.source_id, spec.pdb_path)
    pdb_path = Path(spec.pdb_path)
    if alphafold_like and not spec.allow_alphafold_source:
        return (
            MSAFreeStructureModelContacts(
                source_id=spec.source_id,
                pdb_path=str(pdb_path),
                chain_id=spec.chain_id,
                status="rejected_alphafold_like_source_for_msa_free_ensemble",
                rejected=True,
                parse_error="",
                contact_count=0,
                mean_contact_confidence=0.0,
                contact_map_hash=contact_map_hash(()),
                alphafold_like_source=True,
                msa_or_template_like_source=msa_like,
            ),
            (),
        )
    if msa_like and not spec.allow_msa_or_template_source:
        return (
            MSAFreeStructureModelContacts(
                source_id=spec.source_id,
                pdb_path=str(pdb_path),
                chain_id=spec.chain_id,
                status="rejected_msa_or_template_like_source_for_msa_free_ensemble",
                rejected=True,
                parse_error="",
                contact_count=0,
                mean_contact_confidence=0.0,
                contact_map_hash=contact_map_hash(()),
                alphafold_like_source=alphafold_like,
                msa_or_template_like_source=True,
            ),
            (),
        )
    if not pdb_path.exists() or pdb_path.stat().st_size <= 0:
        return (
            MSAFreeStructureModelContacts(
                source_id=spec.source_id,
                pdb_path=str(pdb_path),
                chain_id=spec.chain_id,
                status="missing_predicted_structure_file",
                rejected=False,
                parse_error="",
                contact_count=0,
                mean_contact_confidence=0.0,
                contact_map_hash=contact_map_hash(()),
                alphafold_like_source=alphafold_like,
                msa_or_template_like_source=msa_like,
            ),
            (),
        )
    try:
        evidence = contact_evidence_from_predicted_pdb(
            row=row,
            pdb_path=pdb_path,
            source_id=spec.source_id,
            source_kind="msa_free_learned_structure_model_contacts_v0",
            source_family="msa_free_learned_structure_model",
            chain_id=spec.chain_id,
        )
    except Exception as exc:  # noqa: BLE001 - safe parse boundary for external predictor outputs
        return (
            MSAFreeStructureModelContacts(
                source_id=spec.source_id,
                pdb_path=str(pdb_path),
                chain_id=spec.chain_id,
                status="predicted_structure_parse_failed",
                rejected=False,
                parse_error=f"{type(exc).__name__}: {exc}",
                contact_count=0,
                mean_contact_confidence=0.0,
                contact_map_hash=contact_map_hash(()),
                alphafold_like_source=alphafold_like,
                msa_or_template_like_source=msa_like,
            ),
            (),
        )
    pairs = tuple(item.pair() for item in evidence)
    return (
        MSAFreeStructureModelContacts(
            source_id=spec.source_id,
            pdb_path=str(pdb_path),
            chain_id=spec.chain_id,
            status="prediction_available",
            rejected=False,
            parse_error="",
            contact_count=len(evidence),
            mean_contact_confidence=_mean_confidence(evidence),
            contact_map_hash=contact_map_hash(pairs),
            alphafold_like_source=alphafold_like,
            msa_or_template_like_source=msa_like,
        ),
        evidence,
    )


def _pair_to_models(
    evidence_by_model: Mapping[str, Sequence[IndependentContactEvidencePair]],
) -> dict[ContactPair, list[IndependentContactEvidencePair]]:
    grouped: dict[ContactPair, list[IndependentContactEvidencePair]] = defaultdict(list)
    for items in evidence_by_model.values():
        for item in items:
            grouped[item.pair()].append(item)
    return dict(grouped)


def _support_set(pairs: Iterable[Sequence[int]]) -> set[ContactPair]:
    return set(normalized_contact_pairs(pairs))


def _contact_graph_degrees(pairs: Iterable[ContactPair]) -> dict[int, int]:
    degrees: dict[int, int] = defaultdict(int)
    for left, right in pairs:
        degrees[left] += 1
        degrees[right] += 1
    return dict(degrees)


def _cooperative_context_score(pair: ContactPair, current_pairs: set[ContactPair]) -> float:
    """Score whether a pair sits inside a local predicted contact patch.

    This encodes a small piece of cooperativity: isolated contacts are less
    trusted than contacts embedded in a local cluster.  It uses only predicted
    contacts from prior refinement iterations, never native truth.
    """

    if not current_pairs:
        return 0.50
    left, right = pair
    near = 0
    for delta_left in range(-2, 3):
        for delta_right in range(-2, 3):
            if delta_left == 0 and delta_right == 0:
                continue
            candidate = (left + delta_left, right + delta_right)
            if candidate[0] < 1 or candidate[1] <= candidate[0]:
                continue
            if candidate in current_pairs:
                near += 1
    return _rounded(min(1.0, 0.30 + near / 12.0))


def _loop_entropy_score(pair: ContactPair, sequence_length: int) -> float:
    """Bounded proxy for loop-closure entropy.

    Long-range contacts are allowed, but extremely long closures receive a mild
    penalty unless model consensus/physics support overcomes it.
    """

    separation = pair[1] - pair[0]
    if separation <= 24:
        return 0.88
    ratio = separation / max(1, sequence_length)
    penalty = min(0.42, log(1.0 + separation / 24.0) / 5.0 + 0.18 * ratio)
    return _rounded(1.0 - penalty)


def _candidate_support_score(pair: ContactPair, candidate_region_pairs: set[ContactPair]) -> float:
    return 1.0 if pair in candidate_region_pairs else 0.35


def _coupling_support_score(pair: ContactPair, coupling_pairs: set[ContactPair]) -> float:
    return 1.0 if pair in coupling_pairs else 0.25


def _build_decisions(
    *,
    row: RealCoordinateVisualRow,
    pair_model_items: Mapping[ContactPair, Sequence[IndependentContactEvidencePair]],
    usable_model_count: int,
    candidate_pairs: set[ContactPair],
    coupling_pairs: set[ContactPair],
    candidate_region_pairs: set[ContactPair],
    current_pairs: set[ContactPair],
    iteration_index: int,
    selection_threshold: float,
    max_selected_contacts: int,
) -> tuple[MSAFreeConsensusContactDecision, ...]:
    physical_priors = build_sequence_physical_prior_scores(
        row=row,
        candidate_pairs=candidate_pairs,
        current_pairs=current_pairs,
    )
    decisions: list[MSAFreeConsensusContactDecision] = []
    for pair in sorted(candidate_pairs):
        items = tuple(pair_model_items.get(pair, ()))
        model_ids = tuple(sorted({item.source_id for item in items}))
        vote_count = len(model_ids)
        vote_fraction = _rounded(vote_count / max(1, usable_model_count))
        raw_mean_conf = _mean_confidence(items)
        # Ordinary PDB B-factors are not calibrated model confidence.  When a
        # non-AlphaFold MSA-free predictor leaves this field near zero, treat
        # confidence as unknown rather than as evidence against the contact.
        mean_conf = max(raw_mean_conf, 0.72) if items else 0.0
        prior = physical_priors.get(pair)
        if prior is None:
            # Should be rare, but keep the score safe and native-free.
            physical_score = 0.0
            energy_score = 0.0
            ss_score = 0.0
            degree_score = 0.0
        else:
            physical_score = prior.physical_prior_score
            energy_score = prior.contact_energy_score
            ss_score = prior.secondary_structure_score
            degree_score = prior.degree_consistency_score
        coop_score = _cooperative_context_score(pair, current_pairs)
        entropy_score = _loop_entropy_score(pair, row.sequence_length)
        coupling_score = _coupling_support_score(pair, coupling_pairs)
        candidate_score = _candidate_support_score(pair, candidate_region_pairs)
        final_score = _rounded(
            0.46 * vote_fraction
            + 0.18 * mean_conf
            + 0.16 * physical_score
            + 0.08 * coop_score
            + 0.06 * entropy_score
            + 0.04 * coupling_score
            + 0.02 * candidate_score
        )
        has_model = vote_count > 0
        selected = bool(has_model and final_score >= selection_threshold)
        reason = "selected_model_consensus_physics_refined" if selected else "below_threshold_or_no_model_vote"
        decisions.append(
            MSAFreeConsensusContactDecision(
                kind=MSA_FREE_MODEL_CONSENSUS_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                i=pair[0],
                j=pair[1],
                model_vote_count=vote_count,
                model_vote_fraction=vote_fraction,
                model_source_ids=model_ids,
                mean_model_confidence=mean_conf,
                physical_prior_score=physical_score,
                contact_energy_score=energy_score,
                secondary_structure_score=ss_score,
                degree_consistency_score=degree_score,
                coupling_support_score=coupling_score,
                candidate_region_support_score=candidate_score,
                cooperative_context_score=coop_score,
                loop_entropy_score=entropy_score,
                final_score=final_score,
                selected=selected,
                selection_reason=reason,
                iteration_index=iteration_index,
            )
        )
    selected_decisions = [item for item in decisions if item.selected]
    if len(selected_decisions) > max_selected_contacts:
        keep = {
            item.pair()
            for item in sorted(
                selected_decisions,
                key=lambda item: (
                    item.final_score,
                    item.model_vote_count,
                    item.mean_model_confidence,
                    item.sequence_separation if hasattr(item, "sequence_separation") else item.j - item.i,
                ),
                reverse=True,
            )[:max_selected_contacts]
        }
        decisions = [
            item
            if item.pair() in keep or not item.selected
            else MSAFreeConsensusContactDecision(
                **{**item.to_dict(), "model_source_ids": tuple(item.model_source_ids), "selected": False, "selection_reason": "trimmed_by_contact_budget"}
            )
            for item in decisions
        ]
    return tuple(decisions)


def run_msa_free_model_ensemble(
    *,
    row: RealCoordinateVisualRow,
    model_specs: Sequence[MSAFreeStructureModelSpec],
    coupling_pairs: Iterable[Sequence[int]] = (),
    candidate_region_pairs: Iterable[Sequence[int]] = (),
    iteration_count: int = 4,
    selection_threshold: float = 0.62,
    solved_precision_threshold: float = 0.70,
    solved_recall_threshold: float = 0.70,
    max_selected_contacts: int | None = None,
) -> MSAFreeModelEnsemblePacket:
    """Run learned-geometry consensus + physics refinement and audit afterward.

    Native contacts are used only at the final evaluation step through
    ``evaluate_contact_prediction``.  They are not part of model loading,
    scoring, consensus, iteration, or selection.
    """

    max_contacts = max_selected_contacts or max(64, int(round(row.sequence_length * 3.05)))
    iterations = max(1, int(iteration_count))
    coupling_support = _support_set(coupling_pairs)
    candidate_region_support = _support_set(candidate_region_pairs)

    model_records: list[MSAFreeStructureModelContacts] = []
    evidence_by_model: dict[str, tuple[IndependentContactEvidencePair, ...]] = {}
    for spec in model_specs:
        record, evidence = _load_model_contacts(row=row, spec=spec)
        model_records.append(record)
        if record.status == "prediction_available" and not record.rejected and evidence:
            evidence_by_model[record.source_id] = evidence

    usable_model_count = len(evidence_by_model)
    pair_model_items = _pair_to_models(evidence_by_model)
    model_pairs = set(pair_model_items)
    all_candidate_pairs = set(model_pairs) | coupling_support | candidate_region_support
    if not all_candidate_pairs:
        empty_metric = _empty_metric(row)
        return MSAFreeModelEnsemblePacket(
            kind=MSA_FREE_MODEL_ENSEMBLE_REPORT_KIND,
            row_id=row.row_id,
            source_accession=row.source_accession,
            sequence_hash=row.sequence_sha256,
            sequence_length=row.sequence_length,
            model_count_requested=len(model_specs),
            usable_model_count=0,
            rejected_model_count=sum(1 for item in model_records if item.rejected),
            available_model_contact_count=0,
            candidate_pair_count=0,
            selected_contact_count=0,
            selected_long_range_contact_count=0,
            iteration_count_requested=iterations,
            iteration_count_executed=0,
            converged=True,
            selection_threshold=selection_threshold,
            max_selected_contacts=max_contacts,
            direct_best_model_source_id="none",
            direct_best_model_metric=empty_metric,
            consensus_metric=empty_metric,
            folding_problem_solved=False,
            direct_structure_solved=False,
            consensus_contact_collapse_solved=False,
            folding_solution_mode="abstain_no_usable_msa_free_model_output",
            claim_rejection_reason="missing_usable_msa_free_model_output",
            models=tuple(model_records),
            decisions=(),
            selected_pairs_hash=contact_map_hash(()),
        )

    current_pairs = set(model_pairs)
    final_decisions: tuple[MSAFreeConsensusContactDecision, ...] = ()
    converged = False
    executed = 0
    for iteration_index in range(1, iterations + 1):
        executed = iteration_index
        final_decisions = _build_decisions(
            row=row,
            pair_model_items=pair_model_items,
            usable_model_count=usable_model_count,
            candidate_pairs=all_candidate_pairs,
            coupling_pairs=coupling_support,
            candidate_region_pairs=candidate_region_support,
            current_pairs=current_pairs,
            iteration_index=iteration_index,
            selection_threshold=selection_threshold,
            max_selected_contacts=max_contacts,
        )
        next_pairs = {decision.pair() for decision in final_decisions if decision.selected}
        if next_pairs == current_pairs:
            converged = True
            current_pairs = next_pairs
            break
        current_pairs = next_pairs

    selected_pairs = tuple(sorted(current_pairs))
    consensus_metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=selected_pairs,
    )

    best_source_id = "none"
    best_metric = _empty_metric(row)
    for source_id, items in sorted(evidence_by_model.items()):
        metric = evaluate_contact_prediction(
            native_pairs=row.native_contact_pairs(),
            predicted_pairs=[item.pair() for item in items],
        )
        if (metric.contact_map_f1, metric.native_contact_precision, metric.native_contact_recall) > (
            best_metric.contact_map_f1,
            best_metric.native_contact_precision,
            best_metric.native_contact_recall,
        ):
            best_source_id = source_id
            best_metric = metric

    direct_solved = bool(
        best_metric.native_contact_precision >= solved_precision_threshold
        and best_metric.native_contact_recall >= solved_recall_threshold
        and usable_model_count > 0
    )
    consensus_solved = bool(
        consensus_metric.native_contact_precision >= solved_precision_threshold
        and consensus_metric.native_contact_recall >= solved_recall_threshold
        and usable_model_count > 0
    )
    solved = bool(direct_solved or consensus_solved)
    if consensus_solved and usable_model_count >= 2:
        mode = "multi_model_consensus_physics_refined"
    elif consensus_solved:
        mode = "single_model_consensus_physics_refined"
    elif direct_solved:
        mode = "direct_msa_free_single_sequence_structure"
    elif usable_model_count > 0:
        mode = "abstain_msa_free_models_available_but_threshold_not_met"
    else:
        mode = "abstain_no_usable_msa_free_model_output"

    if solved:
        rejection = "none"
    elif usable_model_count <= 0:
        rejection = "missing_usable_msa_free_model_output"
    else:
        rejection = "msa_free_model_ensemble_did_not_reach_precision_recall_thresholds"

    return MSAFreeModelEnsemblePacket(
        kind=MSA_FREE_MODEL_ENSEMBLE_REPORT_KIND,
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        model_count_requested=len(model_specs),
        usable_model_count=usable_model_count,
        rejected_model_count=sum(1 for item in model_records if item.rejected),
        available_model_contact_count=sum(item.contact_count for item in model_records),
        candidate_pair_count=len(all_candidate_pairs),
        selected_contact_count=len(selected_pairs),
        selected_long_range_contact_count=sum(1 for left, right in selected_pairs if right - left >= 24),
        iteration_count_requested=iterations,
        iteration_count_executed=executed,
        converged=converged,
        selection_threshold=selection_threshold,
        max_selected_contacts=max_contacts,
        direct_best_model_source_id=best_source_id,
        direct_best_model_metric=best_metric,
        consensus_metric=consensus_metric,
        folding_problem_solved=solved,
        direct_structure_solved=direct_solved,
        consensus_contact_collapse_solved=consensus_solved,
        folding_solution_mode=mode,
        claim_rejection_reason=rejection,
        models=tuple(model_records),
        decisions=final_decisions,
        selected_pairs_hash=contact_map_hash(selected_pairs),
    )


def scan_model_pdb_specs(
    *,
    root: Path,
    default_chain_id: str | None = None,
    allow_alphafold_source: bool = False,
    allow_msa_or_template_source: bool = False,
) -> tuple[MSAFreeStructureModelSpec, ...]:
    """Find PDB files under a directory and turn them into safe model specs."""

    specs: list[MSAFreeStructureModelSpec] = []
    for path in sorted(root.rglob("*.pdb")):
        if path.name.startswith("._"):
            continue
        source_id = path.stem.replace(" ", "_")
        specs.append(
            MSAFreeStructureModelSpec(
                source_id=source_id,
                pdb_path=str(path),
                chain_id=default_chain_id,
                allow_alphafold_source=allow_alphafold_source,
                allow_msa_or_template_source=allow_msa_or_template_source,
            )
        )
    return tuple(specs)
