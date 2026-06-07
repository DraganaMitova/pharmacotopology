from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from statistics import mean
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactMetricPacket, evaluate_contact_prediction
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    CONTACT_CUTOFF_ANGSTROM,
    MIN_SEQUENCE_SEPARATION,
    CoordinatePoint,
    RealCoordinateVisualRow,
)


INDEPENDENT_CONTACT_EVIDENCE_KIND = "independent_contact_evidence_v0"
INDEPENDENT_CONTACT_ENSEMBLE_KIND = "independent_contact_ensemble_probe_v0"
EXTERNAL_STRUCTURE_CONTACT_SOURCE_KIND = "external_predicted_structure_contacts_v0"
NATIVE_COORDINATE_LEAKAGE_POSITIVE_CONTROL_KIND = (
    "native_coordinate_leakage_positive_control_v0"
)
CANDIDATE_REGION_SOURCE_KIND = "candidate_region_sequence_closure_source_v0"
EXTERNAL_COUPLING_CONTACT_SOURCE_KIND = "external_evolutionary_coupling_contact_source_v0"


@dataclass(frozen=True)
class ExternalCAPoint:
    sequence_index: int
    residue_number: int
    chain_id: str
    x: float
    y: float
    z: float
    confidence: float = 1.0

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class IndependentContactEvidencePair:
    row_id: str
    source_accession: str
    source_id: str
    source_kind: str
    source_family: str
    i: int
    j: int
    confidence: float
    distance_angstrom: float | None = None
    mean_coordinate_confidence: float | None = None
    sequence_separation: int = 0
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> tuple[int, int]:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class EnsembleContactDecision:
    row_id: str
    source_accession: str
    i: int
    j: int
    vote_count: int
    vote_fraction: float
    source_ids: tuple[str, ...]
    source_kinds: tuple[str, ...]
    source_families: tuple[str, ...]
    mean_confidence: float
    max_confidence: float
    selected: bool
    selection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> tuple[int, int]:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["source_ids"] = list(self.source_ids)
        data["source_kinds"] = list(self.source_kinds)
        data["source_families"] = list(self.source_families)
        return data


@dataclass(frozen=True)
class EnsembleContactReport:
    row_id: str
    source_accession: str
    evidence_kind: str
    source_family_count: int
    evidence_pair_count: int
    candidate_region_pair_count: int
    external_coupling_pair_count: int
    independent_structure_pair_count: int
    final_pair_count: int
    final_long_range_pair_count: int
    min_votes_required: int
    require_candidate_region_support: bool
    require_independent_structure_support: bool
    contact_precision: float
    contact_recall: float
    long_range_precision: float
    long_range_recall: float
    native_contact_count: int
    native_long_range_contact_count: int
    true_positive_contacts: int
    true_positive_long_range_contacts: int
    benchmark_claim_allowed: bool
    claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def _distance(left: CoordinatePoint | ExternalCAPoint, right: CoordinatePoint | ExternalCAPoint) -> float:
    return sqrt(
        (left.x - right.x) ** 2
        + (left.y - right.y) ** 2
        + (left.z - right.z) ** 2
    )


def parse_ca_pdb_points(
    pdb_path: Path,
    *,
    chain_id: str | None = None,
) -> tuple[ExternalCAPoint, ...]:
    """Parse C-alpha coordinates from a PDB file without exposing sequence labels.

    AlphaFold DB PDB files store pLDDT-like residue confidence in the B-factor
    column.  For ordinary PDB files the same field is just B-factor, so callers
    decide whether to treat it as model confidence.
    """

    points: list[ExternalCAPoint] = []
    seen_residues: set[tuple[str, int, str]] = set()
    sequence_index = 0
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        parsed_chain = line[21].strip() or "A"
        if chain_id is not None and parsed_chain != chain_id:
            continue
        try:
            residue_number = int(line[22:26].strip())
            insertion_code = line[26].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            confidence = float(line[60:66].strip()) / 100.0
        except ValueError:
            continue
        residue_key = (parsed_chain, residue_number, insertion_code)
        if residue_key in seen_residues:
            continue
        seen_residues.add(residue_key)
        sequence_index += 1
        points.append(
            ExternalCAPoint(
                sequence_index=sequence_index,
                residue_number=residue_number,
                chain_id=parsed_chain,
                x=x,
                y=y,
                z=z,
                confidence=_rounded(confidence),
            )
        )
    return tuple(points)


def contact_evidence_from_points(
    *,
    row: RealCoordinateVisualRow,
    points: Sequence[CoordinatePoint | ExternalCAPoint],
    source_id: str,
    source_kind: str,
    source_family: str,
    contact_cutoff_angstrom: float = CONTACT_CUTOFF_ANGSTROM,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
    coordinate_truth_used_before_selection: bool = False,
    native_truth_used_before_selection: bool = False,
) -> tuple[IndependentContactEvidencePair, ...]:
    if len(points) != row.sequence_length:
        raise ValueError(
            "external coordinate point count does not match benchmark sequence length: "
            f"{len(points)} != {row.sequence_length} for {row.source_accession}"
        )
    ordered = tuple(sorted(points, key=lambda point: point.sequence_index))
    evidence: list[IndependentContactEvidencePair] = []
    for left_index, left in enumerate(ordered):
        for right in ordered[left_index + 1 :]:
            separation = right.sequence_index - left.sequence_index
            if separation < minimum_sequence_separation:
                continue
            distance = _distance(left, right)
            if distance > contact_cutoff_angstrom:
                continue
            left_conf = getattr(left, "confidence", 1.0)
            right_conf = getattr(right, "confidence", 1.0)
            mean_confidence = _rounded(mean([float(left_conf), float(right_conf)]))
            evidence.append(
                IndependentContactEvidencePair(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    source_id=source_id,
                    source_kind=source_kind,
                    source_family=source_family,
                    i=left.sequence_index,
                    j=right.sequence_index,
                    confidence=mean_confidence,
                    distance_angstrom=_score(distance),
                    mean_coordinate_confidence=mean_confidence,
                    sequence_separation=separation,
                    coordinate_truth_used_before_selection=coordinate_truth_used_before_selection,
                    native_truth_used_before_selection=native_truth_used_before_selection,
                )
            )
    return tuple(evidence)


def contact_evidence_from_predicted_pdb(
    *,
    row: RealCoordinateVisualRow,
    pdb_path: Path,
    source_id: str,
    source_family: str = "independent_structure",
    source_kind: str = EXTERNAL_STRUCTURE_CONTACT_SOURCE_KIND,
    chain_id: str | None = None,
    contact_cutoff_angstrom: float = CONTACT_CUTOFF_ANGSTROM,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> tuple[IndependentContactEvidencePair, ...]:
    points = parse_ca_pdb_points(pdb_path, chain_id=chain_id)
    return contact_evidence_from_points(
        row=row,
        points=points,
        source_id=source_id,
        source_kind=source_kind,
        source_family=source_family,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        minimum_sequence_separation=minimum_sequence_separation,
        coordinate_truth_used_before_selection=False,
        native_truth_used_before_selection=False,
    )


def native_coordinate_positive_control_evidence(
    row: RealCoordinateVisualRow,
    *,
    source_id: str = "native_coordinate_positive_control",
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> tuple[IndependentContactEvidencePair, ...]:
    return contact_evidence_from_points(
        row=row,
        points=row.coordinate_points,
        source_id=source_id,
        source_kind=NATIVE_COORDINATE_LEAKAGE_POSITIVE_CONTROL_KIND,
        source_family="independent_structure",
        minimum_sequence_separation=minimum_sequence_separation,
        coordinate_truth_used_before_selection=True,
        native_truth_used_before_selection=True,
    )


def candidate_region_evidence_from_events(
    *,
    row: RealCoordinateVisualRow,
    events: Sequence[NucleusClosureEvent],
    confidence: float = 0.50,
) -> tuple[IndependentContactEvidencePair, ...]:
    evidence_by_pair: dict[tuple[int, int], IndependentContactEvidencePair] = {}
    for event in events:
        for pair in event.candidate_region_pairs():
            if pair[0] < 1 or pair[1] > row.sequence_length:
                continue
            if pair not in evidence_by_pair:
                evidence_by_pair[pair] = IndependentContactEvidencePair(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    source_id="candidate_region_pool",
                    source_kind=CANDIDATE_REGION_SOURCE_KIND,
                    source_family="candidate_region",
                    i=pair[0],
                    j=pair[1],
                    confidence=_rounded(confidence),
                    sequence_separation=pair[1] - pair[0],
                )
    return tuple(sorted(evidence_by_pair.values(), key=lambda item: item.pair()))


def coupling_contact_evidence(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
) -> tuple[IndependentContactEvidencePair, ...]:
    evidence: list[IndependentContactEvidencePair] = []
    for constraint in constraints:
        pair = constraint.pair()
        if pair[0] < 1 or pair[1] > row.sequence_length:
            continue
        evidence.append(
            IndependentContactEvidencePair(
                row_id=row.row_id,
                source_accession=row.source_accession,
                source_id="external_dca_coupling",
                source_kind=EXTERNAL_COUPLING_CONTACT_SOURCE_KIND,
                source_family="external_coupling",
                i=pair[0],
                j=pair[1],
                confidence=_rounded(constraint.confidence),
                sequence_separation=pair[1] - pair[0],
            )
        )
    return tuple(evidence)


def load_contact_evidence_json(path: Path) -> tuple[IndependentContactEvidencePair, ...]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(parsed, Mapping):
        rows = parsed.get("contacts", parsed.get("pairs", ()))
    else:
        rows = parsed
    if not isinstance(rows, list):
        raise ValueError(f"independent contact evidence must be a list or object with contacts: {path}")
    evidence: list[IndependentContactEvidencePair] = []
    for index, item in enumerate(rows):
        if not isinstance(item, Mapping):
            raise ValueError(f"contact evidence row {index} is not an object")
        try:
            i = int(item["i"])
            j = int(item["j"])
            row_id = str(item["row_id"])
            source_accession = str(item["source_accession"])
        except KeyError as exc:
            raise ValueError(f"contact evidence row {index} missing field {exc.args[0]}") from exc
        if j <= i:
            i, j = j, i
        evidence.append(
            IndependentContactEvidencePair(
                row_id=row_id,
                source_accession=source_accession,
                source_id=str(item.get("source_id", path.stem)),
                source_kind=str(item.get("source_kind", EXTERNAL_STRUCTURE_CONTACT_SOURCE_KIND)),
                source_family=str(item.get("source_family", "independent_structure")),
                i=i,
                j=j,
                confidence=_rounded(float(item.get("confidence", 1.0))),
                distance_angstrom=(
                    _score(float(item["distance_angstrom"]))
                    if item.get("distance_angstrom") is not None
                    else None
                ),
                mean_coordinate_confidence=(
                    _rounded(float(item["mean_coordinate_confidence"]))
                    if item.get("mean_coordinate_confidence") is not None
                    else None
                ),
                sequence_separation=int(item.get("sequence_separation", j - i)),
                coordinate_truth_used_before_selection=bool(
                    item.get("coordinate_truth_used_before_selection", False)
                ),
                native_truth_used_before_selection=bool(
                    item.get("native_truth_used_before_selection", False)
                ),
                raw_sequence_exposed=bool(item.get("raw_sequence_exposed", False)),
            )
        )
    return tuple(evidence)


def write_contact_evidence_json(
    path: Path,
    *,
    evidence: Sequence[IndependentContactEvidencePair],
    source_kind: str,
) -> None:
    payload = {
        "kind": INDEPENDENT_CONTACT_EVIDENCE_KIND,
        "source_kind": source_kind,
        "contact_count": len(evidence),
        "contacts": [item.to_dict() for item in evidence],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def ensemble_contact_decisions(
    *,
    row: RealCoordinateVisualRow,
    evidence_sources: Sequence[IndependentContactEvidencePair],
    min_votes_required: int = 2,
    require_candidate_region_support: bool = True,
    require_independent_structure_support: bool = True,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> tuple[EnsembleContactDecision, ...]:
    grouped: dict[tuple[int, int], list[IndependentContactEvidencePair]] = {}
    for item in evidence_sources:
        if item.row_id != row.row_id:
            continue
        pair = item.pair()
        if pair[1] - pair[0] < minimum_sequence_separation:
            continue
        grouped.setdefault(pair, []).append(item)

    decisions: list[EnsembleContactDecision] = []
    for pair, items in sorted(grouped.items()):
        source_families = tuple(sorted({item.source_family for item in items}))
        source_ids = tuple(sorted({item.source_id for item in items}))
        source_kinds = tuple(sorted({item.source_kind for item in items}))
        confidences = [item.confidence for item in items]
        vote_count = len(source_families)
        has_candidate = "candidate_region" in source_families
        has_independent_structure = "independent_structure" in source_families
        selected = (
            vote_count >= min_votes_required
            and (has_candidate or not require_candidate_region_support)
            and (has_independent_structure or not require_independent_structure_support)
        )
        if selected:
            if require_candidate_region_support and require_independent_structure_support:
                reason = "min_votes_with_candidate_and_independent_structure"
            elif require_candidate_region_support:
                reason = "min_votes_with_candidate_region"
            elif require_independent_structure_support:
                reason = "min_votes_with_independent_structure"
            else:
                reason = "min_votes"
        elif require_candidate_region_support and not has_candidate:
            reason = "missing_candidate_region_support"
        elif require_independent_structure_support and not has_independent_structure:
            reason = "missing_independent_structure_support"
        else:
            reason = "insufficient_independent_votes"
        decisions.append(
            EnsembleContactDecision(
                row_id=row.row_id,
                source_accession=row.source_accession,
                i=pair[0],
                j=pair[1],
                vote_count=vote_count,
                vote_fraction=_rounded(vote_count / max(1, len(source_families))),
                source_ids=source_ids,
                source_kinds=source_kinds,
                source_families=source_families,
                mean_confidence=_score(mean(confidences) if confidences else 0.0),
                max_confidence=_score(max(confidences) if confidences else 0.0),
                selected=selected,
                selection_reason=reason,
                coordinate_truth_used_before_selection=any(
                    item.coordinate_truth_used_before_selection for item in items
                ),
                native_truth_used_before_selection=any(
                    item.native_truth_used_before_selection for item in items
                ),
                raw_sequence_exposed=any(item.raw_sequence_exposed for item in items),
            )
        )
    return tuple(decisions)


def evaluate_ensemble_contacts(
    *,
    row: RealCoordinateVisualRow,
    evidence_sources: Sequence[IndependentContactEvidencePair],
    min_votes_required: int = 2,
    require_candidate_region_support: bool = True,
    require_independent_structure_support: bool = True,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> tuple[EnsembleContactReport, tuple[EnsembleContactDecision, ...], ContactMetricPacket]:
    decisions = ensemble_contact_decisions(
        row=row,
        evidence_sources=evidence_sources,
        min_votes_required=min_votes_required,
        require_candidate_region_support=require_candidate_region_support,
        require_independent_structure_support=require_independent_structure_support,
        minimum_sequence_separation=minimum_sequence_separation,
    )
    selected_pairs = {decision.pair() for decision in decisions if decision.selected}
    native_pairs = set(row.native_contact_pairs())
    native_long_pairs = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    selected_long_pairs = {pair for pair in selected_pairs if pair[1] - pair[0] >= 24}
    long_tp = selected_long_pairs & native_long_pairs
    metric = evaluate_contact_prediction(
        native_pairs=native_pairs,
        predicted_pairs=selected_pairs,
    )
    source_families = {item.source_family for item in evidence_sources if item.row_id == row.row_id}
    coordinate_truth_used = any(
        item.coordinate_truth_used_before_selection for item in evidence_sources if item.row_id == row.row_id
    ) or any(decision.coordinate_truth_used_before_selection for decision in decisions if decision.selected)
    native_truth_used = any(
        item.native_truth_used_before_selection for item in evidence_sources if item.row_id == row.row_id
    ) or any(decision.native_truth_used_before_selection for decision in decisions if decision.selected)
    raw_sequence_exposed = any(
        item.raw_sequence_exposed for item in evidence_sources if item.row_id == row.row_id
    ) or any(decision.raw_sequence_exposed for decision in decisions if decision.selected)
    independent_structure_pair_count = len(
        {
            item.pair()
            for item in evidence_sources
            if item.row_id == row.row_id and item.source_family == "independent_structure"
        }
    )
    benchmark_claim_allowed = bool(
        independent_structure_pair_count > 0
        and not coordinate_truth_used
        and not native_truth_used
        and not raw_sequence_exposed
    )
    if independent_structure_pair_count <= 0:
        rejection = "missing_independent_structure_source"
    elif coordinate_truth_used or native_truth_used:
        rejection = "independent_source_is_native_coordinate_leakage_positive_control"
    elif raw_sequence_exposed:
        rejection = "raw_sequence_exposed"
    else:
        rejection = "none"
    report = EnsembleContactReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        evidence_kind=INDEPENDENT_CONTACT_ENSEMBLE_KIND,
        source_family_count=len(source_families),
        evidence_pair_count=len({item.pair() for item in evidence_sources if item.row_id == row.row_id}),
        candidate_region_pair_count=len(
            {
                item.pair()
                for item in evidence_sources
                if item.row_id == row.row_id and item.source_family == "candidate_region"
            }
        ),
        external_coupling_pair_count=len(
            {
                item.pair()
                for item in evidence_sources
                if item.row_id == row.row_id and item.source_family == "external_coupling"
            }
        ),
        independent_structure_pair_count=independent_structure_pair_count,
        final_pair_count=len(selected_pairs),
        final_long_range_pair_count=len(selected_long_pairs),
        min_votes_required=min_votes_required,
        require_candidate_region_support=require_candidate_region_support,
        require_independent_structure_support=require_independent_structure_support,
        contact_precision=metric.native_contact_precision,
        contact_recall=metric.native_contact_recall,
        long_range_precision=_rounded(len(long_tp) / len(selected_long_pairs)) if selected_long_pairs else 0.0,
        long_range_recall=_rounded(len(long_tp) / len(native_long_pairs)) if native_long_pairs else 0.0,
        native_contact_count=len(native_pairs),
        native_long_range_contact_count=len(native_long_pairs),
        true_positive_contacts=metric.true_positive_contacts,
        true_positive_long_range_contacts=len(long_tp),
        benchmark_claim_allowed=benchmark_claim_allowed,
        claim_rejection_reason=rejection,
        coordinate_truth_used_before_selection=coordinate_truth_used,
        native_truth_used_before_selection=native_truth_used,
        raw_sequence_exposed=raw_sequence_exposed,
    )
    return report, decisions, metric
