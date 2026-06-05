from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from pharmacotopology.folding_evolutionary_constraints import (
    COUPLING_CONSTRAINT_KIND,
    EVOLUTIONARY_COUPLING_LAYER_KIND,
    CouplingConstraint,
    CouplingDataset,
)
from pharmacotopology.folding_external_coupling_sources import (
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
    ExternalCouplingQualityPolicy,
    SERIOUS_EXTERNAL_COUPLING_POLICY,
    accepted_external_source_kind,
    rejected_external_source_kind,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


EXTERNAL_COUPLING_IMPORT_KIND = "external_evolutionary_coupling_import_v0"

REQUIRED_EXTERNAL_CONSTRAINT_FIELDS = (
    "row_id",
    "source_accession",
    "i",
    "j",
    "sequence_separation",
    "normalized_separation",
    "confidence",
    "raw_score",
    "apc_corrected_score",
    "rank",
    "rank_fraction",
    "source_kind",
    "msa_source_kind",
    "msa_sha256",
    "msa_depth",
    "effective_sequence_count",
    "effective_sequence_count_over_length",
    "target_coverage",
    "focus_sequence_mapping_confidence",
    "coordinate_truth_used_to_build_constraint",
    "native_truth_used_before_coupling_selection",
    "structure_model_used",
    "raw_sequence_exposed",
)


@dataclass(frozen=True)
class ExternalCouplingRowStatus:
    row_id: str
    source_accession: str
    row_external_status: str
    rejection_reason: str
    raw_constraint_count: int
    accepted_constraint_count: int
    target_coverage: float
    focus_sequence_mapping_confidence: float
    effective_sequence_count_over_length: float
    top_l_couplings_available: bool
    coordinate_truth_used_to_build_constraints: bool
    native_truth_used_before_coupling_selection: bool
    structure_model_used: bool
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalCouplingConstraintAudit:
    row_id: str
    source_accession: str
    constraint_id: str
    i: int
    j: int
    sequence_separation: int
    normalized_separation: float
    confidence: float
    raw_score: float
    apc_corrected_score: float
    rank: int
    rank_fraction: float
    source_kind: str
    msa_source_kind: str
    msa_sha256: str
    msa_depth: int
    effective_sequence_count: float
    effective_sequence_count_over_length: float
    target_coverage: float
    focus_sequence_mapping_confidence: float
    coordinate_truth_used_to_build_constraint: bool
    native_truth_used_before_coupling_selection: bool
    structure_model_used: bool
    native_contact_supported: bool
    monomer_coordinate_unsupported: bool
    possible_interdomain_or_allosteric_signal: bool
    possible_oligomer_interface_signal: bool
    benchmark_counts_as_false_positive: bool
    row_external_status: str
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalCouplingImportResult:
    dataset: CouplingDataset
    row_statuses: tuple[ExternalCouplingRowStatus, ...]
    constraint_audits: tuple[ExternalCouplingConstraintAudit, ...]
    batch_id: str
    policy: ExternalCouplingQualityPolicy
    any_coordinate_taint: bool
    any_native_truth_taint: bool
    any_structure_model_used: bool


def _score(value: float) -> float:
    return round(value, 6)


def _constraint_id(row_id: str, i: int, j: int, rank: int) -> str:
    digest = hashlib.sha256(f"{row_id}:{i}:{j}:{rank}".encode("utf-8")).hexdigest()
    return f"external_{digest[:16]}"


def _require_constraint_fields(raw: Mapping[str, Any]) -> None:
    missing = [key for key in REQUIRED_EXTERNAL_CONSTRAINT_FIELDS if key not in raw]
    if missing:
        raise ValueError(
            "external constraint missing required field(s): "
            + ", ".join(sorted(missing))
        )


def _as_float(raw: Mapping[str, Any], key: str) -> float:
    return float(raw[key])


def _as_int(raw: Mapping[str, Any], key: str) -> int:
    return int(raw[key])


def _source_kind_failure(source_kind: str) -> Optional[str]:
    if rejected_external_source_kind(source_kind):
        return f"rejected external coupling source kind: {source_kind}"
    if not accepted_external_source_kind(source_kind):
        return f"unsupported external coupling source kind: {source_kind}"
    return None


def _row_status(
    *,
    row: RealCoordinateVisualRow,
    raw_constraints: Sequence[Mapping[str, Any]],
    policy: ExternalCouplingQualityPolicy,
    coordinate_taint: bool,
    native_truth_taint: bool,
    structure_model_used: bool,
) -> tuple[str, str, float, float, float, bool]:
    if coordinate_taint or native_truth_taint or structure_model_used:
        if coordinate_taint:
            reason = "coordinate_truth_used_to_build_constraint=true"
        elif native_truth_taint:
            reason = "native_truth_used_before_coupling_selection=true"
        else:
            reason = "structure_model_used=true"
        return (
            "external_couplings_rejected_coordinate_taint",
            reason,
            0.0,
            0.0,
            0.0,
            False,
        )
    if not raw_constraints:
        return (
            "external_couplings_rejected_low_depth",
            "no_external_couplings_for_preregistered_row",
            0.0,
            0.0,
            0.0,
            False,
        )

    target_coverage = min(_as_float(raw, "target_coverage") for raw in raw_constraints)
    mapping_confidence = min(
        _as_float(raw, "focus_sequence_mapping_confidence")
        for raw in raw_constraints
    )
    depth_over_length = min(
        _as_float(raw, "effective_sequence_count_over_length")
        for raw in raw_constraints
    )
    top_l_available = len(raw_constraints) >= row.sequence_length
    explicit_top_l = [
        bool(raw.get("top_L_couplings_available"))
        for raw in raw_constraints
        if "top_L_couplings_available" in raw
    ]
    if explicit_top_l:
        top_l_available = all(explicit_top_l)

    if target_coverage < policy.target_coverage_min:
        return (
            "external_couplings_rejected_low_coverage",
            "target_coverage_below_threshold",
            _score(target_coverage),
            _score(mapping_confidence),
            _score(depth_over_length),
            top_l_available,
        )
    if mapping_confidence < policy.focus_sequence_mapping_confidence_min:
        return (
            "external_couplings_rejected_mapping_ambiguous",
            "focus_sequence_mapping_confidence_below_threshold",
            _score(target_coverage),
            _score(mapping_confidence),
            _score(depth_over_length),
            top_l_available,
        )
    if depth_over_length < policy.effective_sequence_count_over_length_min:
        return (
            "external_couplings_rejected_low_depth",
            "effective_sequence_count_over_length_below_threshold",
            _score(target_coverage),
            _score(mapping_confidence),
            _score(depth_over_length),
            top_l_available,
        )
    if policy.require_top_l_couplings and not top_l_available:
        return (
            "external_couplings_rejected_low_depth",
            "top_L_couplings_not_available",
            _score(target_coverage),
            _score(mapping_confidence),
            _score(depth_over_length),
            top_l_available,
        )
    return (
        "external_couplings_available",
        "",
        _score(target_coverage),
        _score(mapping_confidence),
        _score(depth_over_length),
        top_l_available,
    )


def import_external_coupling_dataset(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    external_coupling_file: Path,
    policy: ExternalCouplingQualityPolicy = SERIOUS_EXTERNAL_COUPLING_POLICY,
) -> ExternalCouplingImportResult:
    parsed = json.loads(external_coupling_file.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError("external coupling file must contain a JSON object")

    batch_id = str(parsed.get("batch_id", ""))
    if batch_id != EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID:
        raise ValueError(f"unsupported external coupling batch_id: {batch_id}")
    if str(parsed.get("layer_kind", "")) != EVOLUTIONARY_COUPLING_LAYER_KIND:
        raise ValueError("external coupling file has unsupported layer_kind")
    if str(parsed.get("constraint_kind", "")) != COUPLING_CONSTRAINT_KIND:
        raise ValueError("external coupling file has unsupported constraint_kind")
    if bool(parsed.get("raw_sequence_exposed", False)):
        raise ValueError("external coupling file must not expose raw sequence text")

    coupling_source_kind = str(parsed.get("coupling_source_kind", ""))
    source_failure = _source_kind_failure(coupling_source_kind)
    if source_failure is not None:
        raise ValueError(source_failure)

    expected_row_ids = tuple(row.row_id for row in rows)
    preregistered = tuple(parsed.get("benchmark_row_ids_preregistered", ()))
    if preregistered != expected_row_ids:
        raise ValueError("benchmark_row_ids_preregistered must match frozen rows")

    raw_constraints_payload = parsed.get("constraints", [])
    if not isinstance(raw_constraints_payload, list):
        raise ValueError("external coupling constraints must be a list")
    if not all(isinstance(raw, Mapping) for raw in raw_constraints_payload):
        raise ValueError("external coupling constraints must be JSON objects")
    raw_constraints = tuple(
        raw for raw in raw_constraints_payload if isinstance(raw, Mapping)
    )

    row_by_id = {row.row_id: row for row in rows}
    raw_by_row: dict[str, list[Mapping[str, Any]]] = {
        row.row_id: [] for row in rows
    }
    seen_pairs: set[tuple[str, int, int]] = set()
    for raw in raw_constraints:
        _require_constraint_fields(raw)
        row_id = str(raw.get("row_id", ""))
        row = row_by_id.get(row_id)
        if row is None:
            raise ValueError(f"external constraint row not preregistered: {row_id}")
        if str(raw.get("source_accession", "")) != row.source_accession:
            raise ValueError(f"source accession mismatch for {row_id}")
        constraint_source_kind = str(raw.get("source_kind", coupling_source_kind))
        source_failure = _source_kind_failure(constraint_source_kind)
        if source_failure is not None:
            raise ValueError(source_failure)
        if bool(raw.get("raw_sequence_exposed", False)):
            raise ValueError("external constraint must not expose raw sequence text")
        pair_key = (row_id, int(raw["i"]), int(raw["j"]))
        if pair_key in seen_pairs:
            raise ValueError(f"duplicate external coupling pair rejected: {pair_key}")
        seen_pairs.add(pair_key)
        raw_by_row[row_id].append(raw)

    top_level_external_used = bool(
        parsed.get("external_evolutionary_couplings_used", False)
    )
    if top_level_external_used and not raw_constraints:
        raise ValueError(
            "external_evolutionary_couplings_used=true with zero external constraints"
        )

    any_coordinate_taint = bool(
        parsed.get("coordinate_truth_used_to_build_constraints", False)
    )
    any_native_truth_taint = bool(
        parsed.get("native_truth_used_before_coupling_selection", False)
    )
    any_structure_model_used = bool(parsed.get("structure_model_used", False))
    accepted_constraints: list[CouplingConstraint] = []
    audits: list[ExternalCouplingConstraintAudit] = []
    statuses: list[ExternalCouplingRowStatus] = []

    status_by_row: dict[str, str] = {}
    for row in rows:
        row_raw = tuple(raw_by_row[row.row_id])
        row_coordinate_taint = any(
            bool(raw.get("coordinate_truth_used_to_build_constraint", False))
            for raw in row_raw
        )
        row_native_taint = any(
            bool(raw.get("native_truth_used_before_coupling_selection", False))
            for raw in row_raw
        )
        row_structure_used = any(
            bool(raw.get("structure_model_used", False)) for raw in row_raw
        )
        any_coordinate_taint = any_coordinate_taint or row_coordinate_taint
        any_native_truth_taint = any_native_truth_taint or row_native_taint
        any_structure_model_used = any_structure_model_used or row_structure_used
        (
            status,
            reason,
            target_coverage,
            mapping_confidence,
            depth_over_length,
            top_l_available,
        ) = _row_status(
            row=row,
            raw_constraints=row_raw,
            policy=policy,
            coordinate_taint=row_coordinate_taint,
            native_truth_taint=row_native_taint,
            structure_model_used=row_structure_used,
        )
        status_by_row[row.row_id] = status
        statuses.append(
            ExternalCouplingRowStatus(
                row_id=row.row_id,
                source_accession=row.source_accession,
                row_external_status=status,
                rejection_reason=reason,
                raw_constraint_count=len(row_raw),
                accepted_constraint_count=(
                    len(row_raw) if status == "external_couplings_available" else 0
                ),
                target_coverage=target_coverage,
                focus_sequence_mapping_confidence=mapping_confidence,
                effective_sequence_count_over_length=depth_over_length,
                top_l_couplings_available=top_l_available,
                coordinate_truth_used_to_build_constraints=row_coordinate_taint,
                native_truth_used_before_coupling_selection=row_native_taint,
                structure_model_used=row_structure_used,
            )
        )

    native_pairs_by_row = {
        row.row_id: set(row.native_contact_pairs())
        for row in rows
    }
    for raw in raw_constraints:
        row = row_by_id[str(raw["row_id"])]
        i = _as_int(raw, "i")
        j = _as_int(raw, "j")
        if j <= i:
            raise ValueError("external coupling constraints must use i < j")
        if i < 1 or j > row.sequence_length:
            raise ValueError(f"external constraint outside row bounds: {row.row_id}")
        sequence_separation = j - i
        if _as_int(raw, "sequence_separation") != sequence_separation:
            raise ValueError(f"invalid external sequence separation: {row.row_id}")
        normalized = _score(sequence_separation / row.sequence_length)
        if abs(_as_float(raw, "normalized_separation") - normalized) > 0.000001:
            raise ValueError(f"invalid external normalized separation: {row.row_id}")
        rank = _as_int(raw, "rank")
        native_supported = (i, j) in native_pairs_by_row[row.row_id]
        status = status_by_row[row.row_id]
        audit = ExternalCouplingConstraintAudit(
            row_id=row.row_id,
            source_accession=row.source_accession,
            constraint_id=str(
                raw.get("constraint_id") or _constraint_id(row.row_id, i, j, rank)
            ),
            i=i,
            j=j,
            sequence_separation=sequence_separation,
            normalized_separation=normalized,
            confidence=_score(_as_float(raw, "confidence")),
            raw_score=_score(_as_float(raw, "raw_score")),
            apc_corrected_score=_score(_as_float(raw, "apc_corrected_score")),
            rank=rank,
            rank_fraction=_score(_as_float(raw, "rank_fraction")),
            source_kind=str(raw.get("source_kind", coupling_source_kind)),
            msa_source_kind=str(raw.get("msa_source_kind", "")),
            msa_sha256=str(raw.get("msa_sha256", "")),
            msa_depth=_as_int(raw, "msa_depth"),
            effective_sequence_count=_score(
                _as_float(raw, "effective_sequence_count")
            ),
            effective_sequence_count_over_length=_score(
                _as_float(raw, "effective_sequence_count_over_length")
            ),
            target_coverage=_score(_as_float(raw, "target_coverage")),
            focus_sequence_mapping_confidence=_score(
                _as_float(raw, "focus_sequence_mapping_confidence")
            ),
            coordinate_truth_used_to_build_constraint=bool(
                raw.get("coordinate_truth_used_to_build_constraint", False)
            ),
            native_truth_used_before_coupling_selection=bool(
                raw.get("native_truth_used_before_coupling_selection", False)
            ),
            structure_model_used=bool(raw.get("structure_model_used", False)),
            native_contact_supported=native_supported,
            monomer_coordinate_unsupported=not native_supported,
            possible_interdomain_or_allosteric_signal=not native_supported,
            possible_oligomer_interface_signal=not native_supported,
            benchmark_counts_as_false_positive=not native_supported,
            row_external_status=status,
        )
        audits.append(audit)
        if status != "external_couplings_available":
            continue
        accepted_constraints.append(
            CouplingConstraint(
                row_id=row.row_id,
                source_accession=row.source_accession,
                constraint_id=audit.constraint_id,
                i=i,
                j=j,
                sequence_separation=sequence_separation,
                normalized_separation=normalized,
                confidence=audit.confidence,
                constraint_class=str(
                    raw.get("constraint_class", "external_dca_coupling")
                ),
                source_kind=audit.source_kind,
                coordinate_truth_used_to_build_constraint=(
                    audit.coordinate_truth_used_to_build_constraint
                ),
                native_truth_used_before_coupling_selection=(
                    audit.native_truth_used_before_coupling_selection
                ),
                structure_model_used=audit.structure_model_used,
                raw_score=audit.raw_score,
                apc_corrected_score=audit.apc_corrected_score,
                rank=audit.rank,
                rank_fraction=audit.rank_fraction,
            )
        )

    accepted_external_couplings_used = bool(accepted_constraints)
    if top_level_external_used and not accepted_external_couplings_used:
        raise ValueError(
            "external_evolutionary_couplings_used=true but no constraints passed "
            "external quality gates"
        )

    dataset = CouplingDataset(
        layer_kind=EVOLUTIONARY_COUPLING_LAYER_KIND,
        constraint_kind=COUPLING_CONSTRAINT_KIND,
        source_benchmark_file=str(parsed.get("source_benchmark_file", "")),
        source_constraint_kind=str(parsed.get("source_constraint_kind", "")),
        coupling_source_kind=coupling_source_kind,
        coordinate_truth_used_to_build_constraints=any_coordinate_taint,
        native_truth_used_before_coupling_selection=any_native_truth_taint,
        external_evolutionary_couplings_used=accepted_external_couplings_used,
        raw_sequence_exposed=False,
        constraints=tuple(accepted_constraints),
        structure_model_used_before_coupling_selection=any_structure_model_used,
    )
    return ExternalCouplingImportResult(
        dataset=dataset,
        row_statuses=tuple(statuses),
        constraint_audits=tuple(audits),
        batch_id=batch_id,
        policy=policy,
        any_coordinate_taint=any_coordinate_taint,
        any_native_truth_taint=any_native_truth_taint,
        any_structure_model_used=any_structure_model_used,
    )


def write_imported_external_coupling_dataset(
    result: ExternalCouplingImportResult,
    path: Path,
) -> Path:
    payload = {
        "layer_kind": result.dataset.layer_kind,
        "batch_id": result.batch_id,
        "constraint_kind": result.dataset.constraint_kind,
        "coupling_source_kind": result.dataset.coupling_source_kind,
        "external_evolutionary_couplings_used": (
            result.dataset.external_evolutionary_couplings_used
        ),
        "coordinate_truth_used_to_build_constraints": (
            result.dataset.coordinate_truth_tainted
        ),
        "native_truth_used_before_coupling_selection": (
            result.dataset.native_truth_tainted
        ),
        "oracle_constraint_control": result.dataset.oracle_constraint_control,
        "raw_sequence_exposed": result.dataset.raw_sequence_exposed,
        "constraints": [constraint.to_dict() for constraint in result.dataset.constraints],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
