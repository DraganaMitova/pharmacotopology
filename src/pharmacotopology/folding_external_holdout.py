from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from pharmacotopology.folding_architecture_axis import (
    ARCHITECTURE_AXIS_BENCHMARK_KIND,
    ARCHITECTURE_AXIS_SIGNATURE_KIND,
    architecture_evidence_packet_from_sequence,
)
from pharmacotopology.folding_axis_adjudication import AXIS_NAMES, UNKNOWN_BY_AXIS
from pharmacotopology.folding_axis_profile import (
    AXIS_PROFILE_SIGNATURE_KIND,
    _profile_axes_for_row,
)
from pharmacotopology.folding_regime_analysis import predict_regime_routed_gate
from pharmacotopology.folding_topology import normalize_sequence


EXTERNAL_HOLDOUT_BENCHMARK_KIND = "external_fold_family_holdout_benchmark"
EXTERNAL_HOLDOUT_SPLIT = "external_fold_family_100"
DEFAULT_DEVELOPMENT_BENCHMARK_FILE = Path(
    "data/folding_benchmarks_real_50.locked.json"
)
REQUIRED_TRUTH_AXES = (
    "secondary_structure_axis",
    "architecture_axis",
    "order_axis",
    "environment_axis",
)


@dataclass(frozen=True)
class ExternalHoldoutRow:
    row_id: str
    source_id: str
    source_kind: str
    sequence: str
    sequence_sha256: str
    length: int
    external_family_id: str
    external_family_name: str
    external_family_group: str
    holdout_split: str
    truth_axes: dict[str, str]
    truth_scope: str
    evidence_notes: tuple[str, ...]

    def to_safe_dict(self) -> dict[str, object]:
        data = asdict(self)
        data.pop("sequence", None)
        return data


def _sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _short_hash(sequence: str) -> str:
    return _sha256(sequence)[:16]


def _bool_mean(values: Sequence[bool]) -> float:
    if not values:
        return 0.0
    return round(sum(1 for value in values if value) / len(values), 6)


def _known_axis(axis: str, value: object) -> bool:
    return str(value) != UNKNOWN_BY_AXIS[axis]


def _axis_match(axis: str, predicted: object, truth: object) -> bool:
    return _known_axis(axis, predicted) and _known_axis(axis, truth) and predicted == truth


def _load_json(path: Path) -> Mapping[str, Any]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{path} must contain a JSON object")
    return parsed


def _raw_rows(data: Mapping[str, Any]) -> list[Mapping[str, object]]:
    rows = data.get("references")
    if not isinstance(rows, list):
        raise ValueError("Holdout dataset must include a references list")
    output: list[Mapping[str, object]] = []
    for index, row in enumerate(rows, start=1):
        if not isinstance(row, Mapping):
            raise ValueError(f"Holdout row {index} must be an object")
        output.append(row)
    return output


def load_external_holdout_rows(path: Path) -> tuple[ExternalHoldoutRow, ...]:
    rows: list[ExternalHoldoutRow] = []
    for index, row in enumerate(_raw_rows(_load_json(path)), start=1):
        sequence = normalize_sequence(str(row.get("sequence", "")))
        truth_axes = row.get("truth_axes", {})
        if not isinstance(truth_axes, Mapping):
            raise ValueError(f"row[{index}].truth_axes must be an object")
        missing_axes = [axis for axis in REQUIRED_TRUTH_AXES if axis not in truth_axes]
        if missing_axes:
            raise ValueError(f"row[{index}].truth_axes missing {','.join(missing_axes)}")
        expected_sha = str(row.get("sequence_sha256", ""))
        actual_sha = _sha256(sequence)
        if expected_sha != actual_sha:
            raise ValueError(f"row[{index}].sequence_sha256_mismatch")
        length = int(row.get("length", 0))
        if length != len(sequence):
            raise ValueError(f"row[{index}].length_mismatch")
        evidence_notes = row.get("evidence_notes", ())
        if isinstance(evidence_notes, list):
            notes = tuple(str(note) for note in evidence_notes)
        else:
            notes = (str(evidence_notes),) if evidence_notes else ()
        rows.append(
            ExternalHoldoutRow(
                row_id=str(row["row_id"]),
                source_id=str(row["source_id"]),
                source_kind=str(row["source_kind"]),
                sequence=sequence,
                sequence_sha256=actual_sha,
                length=length,
                external_family_id=str(row["external_family_id"]),
                external_family_name=str(row["external_family_name"]),
                external_family_group=str(row["external_family_group"]),
                holdout_split=str(row["holdout_split"]),
                truth_axes={axis: str(truth_axes[axis]) for axis in REQUIRED_TRUTH_AXES},
                truth_scope=str(row["truth_scope"]),
                evidence_notes=notes,
            )
        )
    return tuple(rows)


def _development_rows(path: Path) -> list[Mapping[str, object]]:
    data = _load_json(path)
    return _raw_rows(data)


def development_overlap_summary(
    holdout_rows: Sequence[ExternalHoldoutRow],
    *,
    development_benchmark_file: Path,
) -> dict[str, object]:
    development_rows = _development_rows(development_benchmark_file)
    dev_row_ids = {str(row.get("protein_id", row.get("row_id", ""))) for row in development_rows}
    dev_source_ids = {
        str(row.get("source_id", row.get("source_accession", "")))
        for row in development_rows
    }
    dev_sequences = {
        _sha256(normalize_sequence(str(row.get("sequence", ""))))
        for row in development_rows
        if row.get("sequence")
    }
    dev_family_ids = {
        str(row.get("external_family_id", ""))
        for row in development_rows
        if row.get("external_family_id")
    }
    row_overlap = [row.row_id for row in holdout_rows if row.row_id in dev_row_ids]
    source_overlap = [
        row.source_id for row in holdout_rows if row.source_id in dev_source_ids
    ]
    sequence_overlap = [
        row.row_id for row in holdout_rows if row.sequence_sha256 in dev_sequences
    ]
    family_overlap = [
        row.external_family_id
        for row in holdout_rows
        if row.external_family_id in dev_family_ids
    ]
    development_overlap_count = len(set(row_overlap)) + len(set(source_overlap))
    return {
        "development_overlap_count": development_overlap_count,
        "development_row_id_overlap_count": len(set(row_overlap)),
        "development_source_id_overlap_count": len(set(source_overlap)),
        "development_sequence_overlap_count": len(set(sequence_overlap)),
        "development_family_overlap_count": len(set(family_overlap)),
        "development_overlap_rows": tuple(sorted(set(row_overlap))),
        "development_sequence_overlap_rows": tuple(sorted(set(sequence_overlap))),
        "development_family_overlap_ids": tuple(sorted(set(family_overlap))),
        "holdout_non_overlap_valid": (
            development_overlap_count == 0
            and not sequence_overlap
            and not family_overlap
        ),
    }


def validate_holdout_lock(
    rows: Sequence[ExternalHoldoutRow],
    *,
    development_benchmark_file: Path,
) -> dict[str, object]:
    row_ids = [row.row_id for row in rows]
    source_ids = [row.source_id for row in rows]
    sequence_hashes = [row.sequence_sha256 for row in rows]
    family_ids = [row.external_family_id for row in rows]
    violations = []
    if len(rows) != 100:
        violations.append("holdout_row_count_not_100")
    if len(set(row_ids)) != len(row_ids):
        violations.append("duplicate_row_id")
    if len(set(source_ids)) != len(source_ids):
        violations.append("duplicate_source_id")
    if len(set(sequence_hashes)) != len(sequence_hashes):
        violations.append("duplicate_sequence_sha256")
    for index, row in enumerate(rows, start=1):
        if row.holdout_split != EXTERNAL_HOLDOUT_SPLIT:
            violations.append(f"row[{index}].holdout_split_mismatch")
        if row.truth_scope not in {
            "global_chain",
            "solved_fragment",
            "domain_fragment",
            "uncertain",
        }:
            violations.append(f"row[{index}].truth_scope_unknown")
    overlap = development_overlap_summary(
        rows,
        development_benchmark_file=development_benchmark_file,
    )
    if not overlap["holdout_non_overlap_valid"]:
        violations.append("development_overlap_detected")
    return {
        "holdout_row_count": len(rows),
        "holdout_unique_sequence_count": len(set(sequence_hashes)),
        "holdout_unique_family_count": len(set(family_ids)),
        "holdout_lock_valid": not violations,
        "holdout_lock_violations": tuple(violations),
        **overlap,
    }


def _source_row_for_profile(row: ExternalHoldoutRow) -> dict[str, object]:
    routed = predict_regime_routed_gate(row.sequence, protein_id=row.row_id)
    return {
        "protein_id": row.row_id,
        "sequence_length": row.length,
        "protein_regime": routed.regime_prediction.protein_regime,
        "predicted_fold_class": routed.predicted_fold_class,
        "confidence": routed.confidence,
        "forced_prediction": routed.forced_prediction,
        "abstained": routed.abstained,
        "gate_path": " | ".join(routed.gate_path),
        "gate_decision_reason": routed.gate_decision_reason,
    }


def _axis_conflict_axes(
    *,
    predicted_axes: Mapping[str, str],
    truth_axes: Mapping[str, str],
) -> tuple[str, ...]:
    conflicts = []
    for axis in AXIS_NAMES:
        predicted = predicted_axes[axis]
        truth = truth_axes[axis]
        if _known_axis(axis, predicted) and _known_axis(axis, truth) and predicted != truth:
            conflicts.append(axis)
    return tuple(conflicts)


def _failure_cohort(row: Mapping[str, object]) -> str:
    conflict_axes = str(row["conflict_axes"])
    if conflict_axes:
        if "environment_axis" in conflict_axes:
            return "membrane_topology_conflict"
        if "architecture_axis" in conflict_axes:
            return "architecture_axis_conflict"
        if "secondary_structure_axis" in conflict_axes:
            return "secondary_axis_conflict"
        if "order_axis" in conflict_axes:
            return "order_axis_conflict"
    if row["truth_architecture_axis"] == "repeat_like" and row[
        "architecture_axis_prediction"
    ] == "unknown":
        return "repeat_like_under_claimed"
    if row["truth_architecture_axis"] == "fragment_scope" and row[
        "architecture_axis_prediction"
    ] == "unknown":
        return "fragment_scope_over_abstention"
    if row["truth_environment_axis"] == "membrane_like" and row[
        "profile_environment_axis"
    ] == "unknown":
        return "membrane_topology_over_abstention"
    if row["truth_architecture_axis"] == "multidomain_or_segmented" and row[
        "architecture_axis_prediction"
    ] == "unknown":
        return "large_chain_multidomain_under_claimed"
    if row["truth_secondary_structure_axis"] == "alpha_beta_mixed" and row[
        "profile_secondary_structure_axis"
    ] == "weak_or_unknown":
        return "alpha_beta_boundary_ambiguity"
    if row["truth_secondary_structure_axis"] == "beta_rich" and row[
        "profile_secondary_structure_axis"
    ] == "weak_or_unknown":
        return "beta_rich_low_complexity_confusion"
    if row["truth_order_axis"] == "ordered" and row["profile_order_axis"] in {
        "mixed_or_uncertain",
        "disordered_flexible",
    }:
        return "disorder_folded_domain_mimic"
    if bool(row["any_axis_claimed"]):
        return "safe_axis_claimed"
    return "safe_abstention_unresolved"


def external_holdout_rows(
    holdout_rows: Sequence[ExternalHoldoutRow],
) -> list[dict[str, object]]:
    rows = []
    for holdout_row in holdout_rows:
        source_row = _source_row_for_profile(holdout_row)
        profile_axes, _profile_reasons = _profile_axes_for_row(source_row)
        architecture_packet = architecture_evidence_packet_from_sequence(
            holdout_row.sequence,
            protein_id=holdout_row.row_id,
        )
        combined_axes = dict(profile_axes)
        combined_axes["architecture_axis"] = (
            architecture_packet.architecture_axis_prediction
        )
        truth_axes = holdout_row.truth_axes
        profile_conflicts = _axis_conflict_axes(
            predicted_axes=profile_axes,
            truth_axes=truth_axes,
        )
        combined_conflicts = _axis_conflict_axes(
            predicted_axes=combined_axes,
            truth_axes=truth_axes,
        )
        architecture_conflict = (
            _known_axis(
                "architecture_axis",
                architecture_packet.architecture_axis_prediction,
            )
            and _known_axis("architecture_axis", truth_axes["architecture_axis"])
            and architecture_packet.architecture_axis_prediction
            != truth_axes["architecture_axis"]
        )
        safe_axis_claim_count = 0
        unsafe_axis_claim_count = 0
        for axis in AXIS_NAMES:
            predicted = combined_axes[axis]
            truth = truth_axes[axis]
            if not _known_axis(axis, predicted):
                continue
            if _known_axis(axis, truth) and predicted == truth:
                safe_axis_claim_count += 1
            elif _known_axis(axis, truth) and predicted != truth:
                unsafe_axis_claim_count += 1
        safe_axis_recovered_axes = [
            axis
            for axis in AXIS_NAMES
            if _known_axis(axis, combined_axes[axis])
            and not bool(source_row["forced_prediction"])
        ]
        row = {
            "row_id": holdout_row.row_id,
            "source_id": holdout_row.source_id,
            "source_kind": holdout_row.source_kind,
            "sequence_sha256": holdout_row.sequence_sha256,
            "sequence_hash_short": _short_hash(holdout_row.sequence),
            "length": holdout_row.length,
            "external_family_id": holdout_row.external_family_id,
            "external_family_name": holdout_row.external_family_name,
            "external_family_group": holdout_row.external_family_group,
            "holdout_split": holdout_row.holdout_split,
            "truth_scope": holdout_row.truth_scope,
            "evidence_notes": ";".join(holdout_row.evidence_notes),
            "source_predicted_fold_class": source_row["predicted_fold_class"],
            "source_forced_prediction": source_row["forced_prediction"],
            "source_abstained": source_row["abstained"],
            "source_confidence": source_row["confidence"],
            "protein_regime": source_row["protein_regime"],
            "gate_path": source_row["gate_path"],
            "gate_decision_reason": source_row["gate_decision_reason"],
            "profile_secondary_structure_axis": profile_axes[
                "secondary_structure_axis"
            ],
            "profile_order_axis": profile_axes["order_axis"],
            "profile_environment_axis": profile_axes["environment_axis"],
            "profile_architecture_axis": profile_axes["architecture_axis"],
            "architecture_axis_prediction": architecture_packet.architecture_axis_prediction,
            "architecture_axis_confidence": architecture_packet.architecture_axis_confidence,
            "architecture_axis_claim_allowed": architecture_packet.architecture_axis_claim_allowed,
            "architecture_axis_abstention_reason": architecture_packet.architecture_axis_abstention_reason,
            "truth_secondary_structure_axis": truth_axes[
                "secondary_structure_axis"
            ],
            "truth_architecture_axis": truth_axes["architecture_axis"],
            "truth_order_axis": truth_axes["order_axis"],
            "truth_environment_axis": truth_axes["environment_axis"],
            "axis_profile_same_axis_conflict": bool(profile_conflicts),
            "axis_profile_conflict_axes": ";".join(profile_conflicts),
            "architecture_axis_same_axis_conflict": architecture_conflict,
            "combined_same_axis_conflict": bool(combined_conflicts),
            "conflict_axes": ";".join(combined_conflicts),
            "safe_axis_claim_count": safe_axis_claim_count,
            "unsafe_axis_claim_count": unsafe_axis_claim_count,
            "safe_axis_recovered_axes": ";".join(safe_axis_recovered_axes),
            "safe_axis_recovered_count": len(safe_axis_recovered_axes),
            "unsafe_class_recovery": False,
            "guard_override": False,
            "any_axis_claimed": any(
                _known_axis(axis, combined_axes[axis]) for axis in AXIS_NAMES
            ),
            "global_fold_class_claim_allowed": False,
            "axis_profile_claim_allowed": True,
            "architecture_axis_claim_allowed": architecture_packet.architecture_axis_claim_allowed,
            "folding_problem_solved": False,
            "claim_allowed": False,
            "prediction_source_layer": "current_safe_axis_stack",
        }
        row["failure_cohort"] = _failure_cohort(row)
        rows.append(row)
    return rows


def _coverage(rows: Sequence[Mapping[str, object]], key: str, axis: str) -> float:
    return _bool_mean([_known_axis(axis, row[key]) for row in rows])


def _safe_axis_claim_total(rows: Sequence[Mapping[str, object]]) -> int:
    return sum(int(row["safe_axis_claim_count"]) for row in rows)


def _unsafe_axis_claim_total(rows: Sequence[Mapping[str, object]]) -> int:
    return sum(int(row["unsafe_axis_claim_count"]) for row in rows)


def family_summary_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["external_family_id"])].append(row)
    output = []
    for family_id, family_rows in sorted(grouped.items()):
        conflicts = [row for row in family_rows if bool(row["combined_same_axis_conflict"])]
        abstentions = [
            row for row in family_rows if not bool(row["any_axis_claimed"])
        ]
        output.append(
            {
                "external_family_id": family_id,
                "external_family_name": family_rows[0]["external_family_name"],
                "external_family_group": family_rows[0]["external_family_group"],
                "row_count": len(family_rows),
                "axis_profile_coverage": _bool_mean(
                    [bool(row["any_axis_claimed"]) for row in family_rows]
                ),
                "architecture_axis_coverage": _bool_mean(
                    [
                        bool(row["architecture_axis_claim_allowed"])
                        for row in family_rows
                    ]
                ),
                "family_level_abstention_count": len(abstentions),
                "family_level_conflict_count": len(conflicts),
                "dominant_failure_cohort": Counter(
                    str(row["failure_cohort"]) for row in family_rows
                ).most_common(1)[0][0],
            }
        )
    return output


def axis_conflict_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    keys = (
        "row_id",
        "external_family_id",
        "external_family_group",
        "truth_scope",
        "conflict_axes",
        "source_predicted_fold_class",
        "profile_secondary_structure_axis",
        "profile_order_axis",
        "profile_environment_axis",
        "architecture_axis_prediction",
        "truth_secondary_structure_axis",
        "truth_architecture_axis",
        "truth_order_axis",
        "truth_environment_axis",
        "failure_cohort",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["combined_same_axis_conflict"])
    ]


def abstention_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    keys = (
        "row_id",
        "external_family_id",
        "external_family_group",
        "length",
        "truth_scope",
        "protein_regime",
        "source_predicted_fold_class",
        "architecture_axis_abstention_reason",
        "truth_secondary_structure_axis",
        "truth_architecture_axis",
        "truth_order_axis",
        "truth_environment_axis",
        "failure_cohort",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if not bool(row["any_axis_claimed"])
        or not bool(row["architecture_axis_claim_allowed"])
    ]


def failure_cohort_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    counts: dict[tuple[str, str, str, str, str, str, str, str, str], dict[str, object]] = {}
    for row in rows:
        key = (
            str(row["failure_cohort"]),
            str(row["external_family_group"]),
            _length_band(int(row["length"])),
            str(row["truth_scope"]),
            str(row["truth_environment_axis"]),
            str(row["truth_architecture_axis"]),
            str(row["truth_secondary_structure_axis"]),
            str(row["truth_order_axis"]),
            str(row["architecture_axis_abstention_reason"]),
        )
        if key not in counts:
            counts[key] = {
                "failure_cohort": key[0],
                "family_group": key[1],
                "length_band": key[2],
                "truth_scope": key[3],
                "environment_axis": key[4],
                "architecture_axis": key[5],
                "secondary_structure_axis": key[6],
                "order_axis": key[7],
                "abstention_reason": key[8],
                "row_count": 0,
                "conflict_count": 0,
                "abstention_count": 0,
                "example_rows": [],
                "prediction_source_layer": "current_safe_axis_stack",
            }
        item = counts[key]
        item["row_count"] = int(item["row_count"]) + 1
        item["conflict_count"] = int(item["conflict_count"]) + int(
            bool(row["combined_same_axis_conflict"])
        )
        item["abstention_count"] = int(item["abstention_count"]) + int(
            not bool(row["any_axis_claimed"])
        )
        examples = item["example_rows"]
        if isinstance(examples, list) and len(examples) < 6:
            examples.append(row["row_id"])
    output = list(counts.values())
    for item in output:
        item["example_rows"] = ";".join(item["example_rows"])
    return sorted(
        output,
        key=lambda item: (-int(item["row_count"]), str(item["failure_cohort"])),
    )


def _length_band(length: int) -> str:
    if length < 70:
        return "short_fragment_under_70"
    if length < 140:
        return "small_70_139"
    if length < 260:
        return "medium_140_259"
    if length < 420:
        return "large_260_419"
    return "very_large_420_plus"


def build_external_holdout_report(
    holdout_rows: Sequence[ExternalHoldoutRow],
    *,
    holdout_file: Path,
    development_benchmark_file: Path,
) -> dict[str, object]:
    rows = external_holdout_rows(holdout_rows)
    lock = validate_holdout_lock(
        holdout_rows,
        development_benchmark_file=development_benchmark_file,
    )
    family_rows = family_summary_rows(rows)
    axis_profile_same_axis_conflicts = [
        row for row in rows if bool(row["axis_profile_same_axis_conflict"])
    ]
    architecture_same_axis_conflicts = [
        row for row in rows if bool(row["architecture_axis_same_axis_conflict"])
    ]
    forced_same_axis_conflicts = [
        row
        for row in rows
        if bool(row["source_forced_prediction"])
        and bool(row["combined_same_axis_conflict"])
    ]
    high_confidence_wrong = [
        row
        for row in rows
        if bool(row["source_forced_prediction"])
        and float(row["source_confidence"]) >= 0.58
        and bool(row["combined_same_axis_conflict"])
    ]
    family_failure_count = sum(
        1 for row in family_rows if str(row["dominant_failure_cohort"]) != "safe_axis_claimed"
    )
    family_abstention_count = sum(
        int(row["family_level_abstention_count"]) for row in family_rows
    )
    family_conflict_count = sum(
        1 for row in family_rows if int(row["family_level_conflict_count"]) > 0
    )
    unsafe_axis_claim_count = _unsafe_axis_claim_total(rows)
    architecture_axis_claim_allowed = not architecture_same_axis_conflicts
    report = {
        "benchmark_kind": EXTERNAL_HOLDOUT_BENCHMARK_KIND,
        "holdout_file": str(holdout_file),
        "development_benchmark_file": str(development_benchmark_file),
        "holdout_split": EXTERNAL_HOLDOUT_SPLIT,
        "source_axis_profile_signature_kind": AXIS_PROFILE_SIGNATURE_KIND,
        "source_architecture_axis_signature_kind": ARCHITECTURE_AXIS_SIGNATURE_KIND,
        "source_architecture_axis_benchmark_kind": ARCHITECTURE_AXIS_BENCHMARK_KIND,
        "prediction_stack": (
            "regime_routing;fold_axis_safety_guards;axis_profile_recovery;"
            "architecture_axis_adjudication"
        ),
        "prediction_logic_changed_in_this_batch": False,
        "thresholds_changed_in_this_batch": False,
        **lock,
        "collapsed_class_coverage": _bool_mean(
            [bool(row["source_forced_prediction"]) for row in rows]
        ),
        "axis_profile_coverage": _bool_mean(
            [bool(row["any_axis_claimed"]) for row in rows]
        ),
        "secondary_axis_coverage": _coverage(
            rows,
            "profile_secondary_structure_axis",
            "secondary_structure_axis",
        ),
        "architecture_axis_coverage": _bool_mean(
            [bool(row["architecture_axis_claim_allowed"]) for row in rows]
        ),
        "order_axis_coverage": _coverage(rows, "profile_order_axis", "order_axis"),
        "environment_axis_coverage": _coverage(
            rows,
            "profile_environment_axis",
            "environment_axis",
        ),
        "axis_profile_same_axis_conflict_count": len(
            axis_profile_same_axis_conflicts
        ),
        "architecture_axis_same_axis_conflict_count": len(
            architecture_same_axis_conflicts
        ),
        "forced_same_axis_conflict_count": len(forced_same_axis_conflicts),
        "high_confidence_wrong_count_after_axis_scoring": len(high_confidence_wrong),
        "safe_axis_claim_count": _safe_axis_claim_total(rows),
        "unsafe_axis_claim_count": unsafe_axis_claim_count,
        "safe_axis_recovered_count": sum(
            int(row["safe_axis_recovered_count"]) for row in rows
        ),
        "unsafe_class_recovery_count": sum(
            1 for row in rows if bool(row["unsafe_class_recovery"])
        ),
        "guard_override_count": sum(1 for row in rows if bool(row["guard_override"])),
        "family_level_failure_count": family_failure_count,
        "family_level_abstention_count": family_abstention_count,
        "family_level_conflict_count": family_conflict_count,
        "family_generalization_status": (
            "conflicts_detected"
            if unsafe_axis_claim_count
            or axis_profile_same_axis_conflicts
            or architecture_same_axis_conflicts
            else "safe_with_abstentions"
        ),
        "global_fold_class_claim_allowed": False,
        "axis_profile_claim_allowed": True,
        "architecture_axis_claim_allowed": architecture_axis_claim_allowed,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "artifact_reproducible": True,
        "boundary_statement": (
            "This falsification layer runs the current safe axis stack over a "
            "locked non-overlapping holdout. It does not tune thresholds, change "
            "guards, recover global fold classes, export raw sequences in "
            "artifacts, or claim protein folding is solved."
        ),
        "rows": rows,
    }
    return report


def build_external_holdout_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": "external_fold_family_holdout_safety_certificate",
        "external_fold_family_holdout_complete": True,
        "holdout_row_count": report["holdout_row_count"],
        "development_overlap_count": report["development_overlap_count"],
        "development_sequence_overlap_count": report[
            "development_sequence_overlap_count"
        ],
        "holdout_non_overlap_valid": report["holdout_non_overlap_valid"],
        "holdout_lock_valid": report["holdout_lock_valid"],
        "prediction_logic_changed_in_this_batch": False,
        "thresholds_changed_in_this_batch": False,
        "guard_overrides": report["guard_override_count"],
        "unsafe_class_recovery_count": report["unsafe_class_recovery_count"],
        "unsafe_axis_claim_count": report["unsafe_axis_claim_count"],
        "global_fold_class_claim_allowed": report[
            "global_fold_class_claim_allowed"
        ],
        "axis_profile_claim_allowed": report["axis_profile_claim_allowed"],
        "architecture_axis_claim_allowed": report["architecture_axis_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "raw_sequences_exported": False,
        "output_artifacts": tuple(output_names),
    }


def write_external_holdout_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    family_rows: Sequence[Mapping[str, object]],
    conflicts: Sequence[Mapping[str, object]],
    abstentions: Sequence[Mapping[str, object]],
    failure_cohorts: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    family_summary_path: Path,
    axis_conflicts_path: Path,
    abstentions_path: Path,
    failure_cohorts_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    output_names = (
        report_path.name,
        rows_path.name,
        family_summary_path.name,
        axis_conflicts_path.name,
        abstentions_path.name,
        failure_cohorts_path.name,
        dashboard_path.name,
        certificate_path.name,
    )
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(rows, rows_path)
    _write_csv_rows(family_rows, family_summary_path)
    _write_csv_rows(conflicts, axis_conflicts_path)
    _write_csv_rows(abstentions, abstentions_path)
    _write_csv_rows(failure_cohorts, failure_cohorts_path)
    dashboard_path.write_text(render_external_holdout_dashboard(report), encoding="utf-8")
    certificate = build_external_holdout_certificate(
        report,
        output_names=output_names,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        family_summary_path,
        axis_conflicts_path,
        abstentions_path,
        failure_cohorts_path,
        dashboard_path,
        certificate_path,
    )


def _write_csv_rows(rows: Sequence[Mapping[str, object]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _metric_cards(report: Mapping[str, object]) -> str:
    labels = (
        "holdout_row_count",
        "holdout_non_overlap_valid",
        "axis_profile_coverage",
        "architecture_axis_coverage",
        "axis_profile_same_axis_conflict_count",
        "architecture_axis_same_axis_conflict_count",
        "high_confidence_wrong_count_after_axis_scoring",
        "unsafe_class_recovery_count",
        "family_generalization_status",
        "claim_allowed",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _mapping_table(title: str, mapping: object) -> str:
    if not isinstance(mapping, Mapping) or not mapping:
        return ""
    rows = []
    for key, value in mapping.items():
        rows.append(
            "<tr>"
            f"<td>{_escape(key)}</td>"
            f"<td>{_escape(value)}</td>"
            "</tr>"
        )
    return (
        f"<section><h2>{_escape(title)}</h2>"
        "<table><thead><tr><th>key</th><th>value</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def _dominant_failure_table(report: Mapping[str, object]) -> str:
    rows = failure_cohort_rows(report.get("rows", []))[:12]
    body = "".join(
        "<tr>"
        f"<td>{_escape(row['failure_cohort'])}</td>"
        f"<td>{_escape(row['family_group'])}</td>"
        f"<td>{_escape(row['row_count'])}</td>"
        f"<td>{_escape(row['conflict_count'])}</td>"
        f"<td>{_escape(row['abstention_reason'])}</td>"
        f"<td>{_escape(row['example_rows'])}</td>"
        "</tr>"
        for row in rows
    )
    return (
        "<section><h2>Dominant Failure Cohorts</h2>"
        "<table><thead><tr>"
        "<th>cohort</th><th>family group</th><th>rows</th><th>conflicts</th>"
        "<th>abstention reason</th><th>examples</th>"
        "</tr></thead><tbody>"
        + body
        + "</tbody></table></section>"
    )


def render_external_holdout_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>External Fold-Family Holdout</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f6f7f2;
      color: #1f2523;
    }}
    header {{
      padding: 32px;
      background: #21312d;
      color: #f6f7f2;
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      margin-top: 20px;
    }}
    .metric {{
      background: #ffffff;
      border: 1px solid #d5ddd5;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #58635e;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 22px;
    }}
    .rule-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 12px;
    }}
    .rule {{
      background: #fff;
      border: 1px solid #d5ddd5;
      border-radius: 6px;
      padding: 14px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #fff;
      border: 1px solid #d5ddd5;
      border-radius: 6px;
      overflow: hidden;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e8eee8;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
    th {{
      background: #e8eee8;
      color: #34423d;
    }}
  </style>
</head>
<body>
  <header>
    <h1>External Holdout, Not Development Benchmark</h1>
    <p>No threshold tuning in this batch. The current safe axis stack is evaluated as-is on non-overlapping fold-family rows.</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Safety Sections</h2>
      <div class="rule-grid">
        <div class="rule"><strong>No Threshold Tuning In This Batch</strong><br>Prediction logic is called unchanged.</div>
        <div class="rule"><strong>Family Non-Overlap Check</strong><br>Row IDs, source IDs, and sequence hashes are checked against the development 50.</div>
        <div class="rule"><strong>Axis Profile Generalization</strong><br>Partial axis claims are scored by truth axis after prediction.</div>
        <div class="rule"><strong>Architecture Axis Generalization</strong><br>Architecture packets are evaluated without rule repair.</div>
        <div class="rule"><strong>Abstention Is Allowed</strong><br>Coverage drops are valid falsification signals.</div>
        <div class="rule"><strong>Unsafe Class Recovery Remains Forbidden</strong><br>Global fold class stays locked.</div>
        <div class="rule"><strong>Global Fold Class Still Locked</strong><br>claim_allowed remains false.</div>
        <div class="rule"><strong>Next Repair Candidates</strong><br>Use cohorts below; do not infer fixes from this dashboard alone.</div>
      </div>
    </section>
    {_mapping_table("Family Non-Overlap Check", {
        "development_overlap_count": report.get("development_overlap_count"),
        "development_sequence_overlap_count": report.get("development_sequence_overlap_count"),
        "development_family_overlap_count": report.get("development_family_overlap_count"),
        "holdout_non_overlap_valid": report.get("holdout_non_overlap_valid"),
    })}
    {_mapping_table("Axis Coverage", {
        "collapsed_class_coverage": report.get("collapsed_class_coverage"),
        "axis_profile_coverage": report.get("axis_profile_coverage"),
        "secondary_axis_coverage": report.get("secondary_axis_coverage"),
        "architecture_axis_coverage": report.get("architecture_axis_coverage"),
        "order_axis_coverage": report.get("order_axis_coverage"),
        "environment_axis_coverage": report.get("environment_axis_coverage"),
    })}
    {_dominant_failure_table(report)}
  </main>
</body>
</html>
"""
