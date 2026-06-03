from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

from pharmacotopology.folding_hierarchical_gates import load_hierarchical_gate_inputs
from pharmacotopology.folding_regime_analysis import (
    REGIME_ANALYSIS_BENCHMARK_KIND,
    regime_analysis_rows,
)
from pharmacotopology.folding_structure_benchmark import StructureEvidenceRow
from pharmacotopology.folding_topology import (
    HYDROPHOBIC_AMINO_ACIDS,
    FoldingReferenceExample,
    normalize_fold_class,
    normalize_sequence,
    sequence_features,
)


FOLD_AXIS_ADJUDICATION_BENCHMARK_KIND = "orthogonal_fold_axis_truth_adjudication"
AXIS_SIGNATURE_KIND = "orthogonal_fold_axis_projection"

SECONDARY_AXIS_VALUES = (
    "alpha_rich",
    "beta_rich",
    "alpha_beta_mixed",
    "weak_or_unknown",
)
ARCHITECTURE_AXIS_VALUES = (
    "compact_single_domain",
    "multidomain_or_segmented",
    "repeat_like",
    "fragment_scope",
    "unknown",
)
ORDER_AXIS_VALUES = (
    "ordered",
    "disordered_flexible",
    "mixed_or_uncertain",
)
ENVIRONMENT_AXIS_VALUES = (
    "soluble_like",
    "membrane_like",
    "unknown",
)
AXIS_NAMES = (
    "secondary_structure_axis",
    "architecture_axis",
    "order_axis",
    "environment_axis",
)
UNKNOWN_BY_AXIS = {
    "secondary_structure_axis": "weak_or_unknown",
    "architecture_axis": "unknown",
    "order_axis": "mixed_or_uncertain",
    "environment_axis": "unknown",
}


def _mean(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _bool_mean(values: Iterable[bool]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return round(sum(1 for value in values if value) / len(values), 6)


def _fraction(sequence: str, alphabet: frozenset[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in alphabet) / len(sequence)


def _secondary_axis_from_fold_class(fold_class: str) -> str:
    if fold_class in {"alpha_rich", "beta_rich", "alpha_beta_mixed"}:
        return fold_class
    return "weak_or_unknown"


def _label_axes(label_class: str) -> dict[str, str]:
    label_class = normalize_fold_class(label_class)
    return {
        "secondary_structure_axis": _secondary_axis_from_fold_class(label_class),
        "architecture_axis": (
            "multidomain_or_segmented"
            if label_class == "multidomain_boundary"
            else "unknown"
        ),
        "order_axis": (
            "disordered_flexible"
            if label_class == "disordered_flexible"
            else ("ordered" if label_class != "disordered_flexible" else "mixed_or_uncertain")
        ),
        "environment_axis": "unknown",
    }


def _structure_axes(
    *,
    structure_class: str,
    structure: StructureEvidenceRow,
    sequence_length: int,
) -> dict[str, str]:
    features = structure.structure_features
    residue_count = float(features.get("residue_count", sequence_length))
    domain_boundary = float(features.get("domain_boundary_signal", 0.0))
    return {
        "secondary_structure_axis": _secondary_axis_from_fold_class(structure_class),
        "architecture_axis": _structure_architecture_axis(
            structure_class=structure_class,
            residue_count=residue_count,
            domain_boundary=domain_boundary,
        ),
        "order_axis": (
            "disordered_flexible"
            if structure.evidence_kind == "disorder_reference"
            or structure_class == "disordered_flexible"
            else "ordered"
        ),
        "environment_axis": "unknown",
    }


def _structure_architecture_axis(
    *,
    structure_class: str,
    residue_count: float,
    domain_boundary: float,
) -> str:
    if residue_count < 70:
        return "fragment_scope"
    if structure_class == "multidomain_boundary" or (
        residue_count >= 180 and domain_boundary >= 0.55
    ):
        return "multidomain_or_segmented"
    if structure_class in {"alpha_rich", "beta_rich", "alpha_beta_mixed"}:
        return "compact_single_domain"
    return "unknown"


def _predicted_axes(row: Mapping[str, object]) -> dict[str, str]:
    predicted_class = str(row["predicted_fold_class"])
    protein_regime = str(row["protein_regime"])
    forced = bool(row["forced_prediction"])
    abstained = bool(row["abstained"])
    return {
        "secondary_structure_axis": _secondary_axis_from_fold_class(predicted_class),
        "architecture_axis": _predicted_architecture_axis(
            predicted_class=predicted_class,
            protein_regime=protein_regime,
        ),
        "order_axis": _predicted_order_axis(
            predicted_class=predicted_class,
            protein_regime=protein_regime,
            forced=forced,
            abstained=abstained,
        ),
        "environment_axis": (
            "membrane_like"
            if protein_regime == "membrane_like"
            else ("unknown" if protein_regime == "ambiguous_regime" else "soluble_like")
        ),
    }


def _predicted_architecture_axis(
    *,
    predicted_class: str,
    protein_regime: str,
) -> str:
    if protein_regime == "small_peptide_or_fragment":
        return "fragment_scope"
    if protein_regime == "repeat_like":
        return "repeat_like"
    if protein_regime == "multidomain_modular" or predicted_class == "multidomain_boundary":
        return "multidomain_or_segmented"
    if protein_regime == "compact_single_domain":
        return "compact_single_domain"
    return "unknown"


def _predicted_order_axis(
    *,
    predicted_class: str,
    protein_regime: str,
    forced: bool,
    abstained: bool,
) -> str:
    if predicted_class == "disordered_flexible":
        return "disordered_flexible"
    if protein_regime == "intrinsically_disordered" and forced:
        return "disordered_flexible"
    if predicted_class in {
        "alpha_rich",
        "beta_rich",
        "alpha_beta_mixed",
        "multidomain_boundary",
    }:
        return "ordered"
    if abstained:
        return "mixed_or_uncertain"
    return "mixed_or_uncertain"


def _environment_axis_from_sequence(sequence: str) -> str:
    normalized = normalize_sequence(sequence)
    features = sequence_features(normalized)
    hydrophobic = float(features["hydrophobic_fraction"])
    hydrophobic_run = _max_run_fraction(normalized, HYDROPHOBIC_AMINO_ACIDS)
    if len(normalized) >= 180 and hydrophobic >= 0.50:
        return "membrane_like"
    if hydrophobic >= 0.48 and hydrophobic_run >= 0.05:
        return "membrane_like"
    return "unknown"


def _max_run_fraction(sequence: str, alphabet: frozenset[str]) -> float:
    current = 0
    best = 0
    for residue in sequence:
        if residue in alphabet:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return round(best / max(len(sequence), 1), 6)


def _adjudicated_truth_axes(
    *,
    label_axes: Mapping[str, str],
    structure_axes: Mapping[str, str],
    sequence_environment_axis: str,
) -> dict[str, str]:
    truth = {}
    for axis in AXIS_NAMES:
        structure_value = structure_axes[axis]
        label_value = label_axes[axis]
        if axis == "environment_axis":
            truth[axis] = sequence_environment_axis
        elif structure_value != UNKNOWN_BY_AXIS[axis]:
            truth[axis] = structure_value
        elif label_value != UNKNOWN_BY_AXIS[axis]:
            truth[axis] = label_value
        else:
            truth[axis] = UNKNOWN_BY_AXIS[axis]
    return truth


def _axis_scorable(axis: str, predicted_value: str, truth_value: str) -> bool:
    return (
        predicted_value != UNKNOWN_BY_AXIS[axis]
        and truth_value != UNKNOWN_BY_AXIS[axis]
    )


def _axis_match(axis: str, predicted_value: str, truth_value: str) -> bool:
    if not _axis_scorable(axis, predicted_value, truth_value):
        return False
    return predicted_value == truth_value


def _same_axis_structure_label_conflicts(
    label_axes: Mapping[str, str],
    structure_axes: Mapping[str, str],
) -> tuple[str, ...]:
    conflicts = []
    for axis in AXIS_NAMES:
        label_value = label_axes[axis]
        structure_value = structure_axes[axis]
        if (
            label_value != UNKNOWN_BY_AXIS[axis]
            and structure_value != UNKNOWN_BY_AXIS[axis]
            and label_value != structure_value
        ):
            conflicts.append(axis)
    return tuple(conflicts)


def _populated_axis_names(axis_values: Mapping[str, str]) -> tuple[str, ...]:
    return tuple(
        axis for axis in AXIS_NAMES if axis_values[axis] != UNKNOWN_BY_AXIS[axis]
    )


def _taxonomy_collapse_reason(
    *,
    label_class: str,
    structure_class: str,
    label_axes: Mapping[str, str],
    structure_axes: Mapping[str, str],
    same_axis_conflicts: Sequence[str],
) -> str:
    if label_class == structure_class:
        return ""
    if same_axis_conflicts:
        return "single_class_disagreement_contains_same_axis_conflict"
    label_populated = ",".join(_populated_axis_names(label_axes))
    structure_populated = ",".join(_populated_axis_names(structure_axes))
    return (
        "single_class_disagreement_is_orthogonal_axis_projection:"
        f"label_axes={label_populated or 'none'};"
        f"structure_axes={structure_populated or 'none'}"
    )


def _manual_review_reasons(row: Mapping[str, object]) -> tuple[str, ...]:
    reasons = []
    if bool(row["orthogonal_axis_disagreement"]):
        reasons.append("single_class_taxonomy_collapse")
    if bool(row["true_same_axis_conflict"]):
        reasons.append("true_same_axis_conflict")
    if int(row["axis_unscorable_count"]) >= len(AXIS_NAMES):
        reasons.append("axis_unscorable")
    if bool(row["high_confidence_wrong_after_axis_scoring"]):
        reasons.append("high_confidence_wrong_after_axis_scoring")
    if str(row["protein_regime"]) == "membrane_like":
        reasons.append("membrane_regime_not_fold_class")
    if str(row["predicted_architecture_axis"]) == "fragment_scope":
        reasons.append("fragment_evidence_not_global_fold_truth")
    return tuple(reasons)


def axis_adjudication_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    evidence_by_id = {row.protein_id: row for row in evidence_rows}
    regime_rows_by_id = {
        str(row["protein_id"]): row
        for row in regime_analysis_rows(references, evidence_rows)
    }
    rows = []
    for reference in references:
        regime_row = regime_rows_by_id[reference.protein_id]
        structure = evidence_by_id[reference.protein_id]
        label_class = normalize_fold_class(reference.reference_fold_class)
        structure_class = structure.structure_fold_class
        sequence_length = len(normalize_sequence(reference.sequence))
        label_axis_values = _label_axes(label_class)
        structure_axis_values = _structure_axes(
            structure_class=structure_class,
            structure=structure,
            sequence_length=sequence_length,
        )
        sequence_environment_axis = _environment_axis_from_sequence(reference.sequence)
        truth_axis_values = _adjudicated_truth_axes(
            label_axes=label_axis_values,
            structure_axes=structure_axis_values,
            sequence_environment_axis=sequence_environment_axis,
        )
        predicted_axis_values = _predicted_axes(regime_row)
        axis_scores = {}
        axis_unscorable_count = 0
        axis_conflicts = []
        for axis in AXIS_NAMES:
            predicted_value = predicted_axis_values[axis]
            truth_value = truth_axis_values[axis]
            scorable = _axis_scorable(axis, predicted_value, truth_value)
            match = _axis_match(axis, predicted_value, truth_value)
            if not scorable:
                axis_unscorable_count += 1
            elif not match:
                axis_conflicts.append(axis)
            axis_scores[axis] = {
                "predicted": predicted_value,
                "truth": truth_value,
                "scorable": scorable,
                "match": match,
            }
        structure_label_same_axis_conflicts = _same_axis_structure_label_conflicts(
            label_axis_values,
            structure_axis_values,
        )
        structure_label_disagreement = label_class != structure_class
        orthogonal_axis_disagreement = (
            structure_label_disagreement and not structure_label_same_axis_conflicts
        )
        high_confidence_wrong_after_axis = (
            bool(regime_row["forced_prediction"])
            and float(regime_row["confidence"]) >= 0.58
            and bool(axis_conflicts)
        )
        row: dict[str, object] = {
            "protein_id": reference.protein_id,
            "sequence_length": sequence_length,
            "axis_signature_kind": AXIS_SIGNATURE_KIND,
            "protein_regime": regime_row["protein_regime"],
            "predicted_fold_class": regime_row["predicted_fold_class"],
            "structure_fold_class": structure_class,
            "label_fold_class": label_class,
            "confidence": regime_row["confidence"],
            "forced_prediction": regime_row["forced_prediction"],
            "abstained": regime_row["abstained"],
            "gate_path": regime_row["gate_path"],
            "gate_decision_reason": regime_row["gate_decision_reason"],
            "single_class_label_structure_disagreement": structure_label_disagreement,
            "orthogonal_axis_disagreement": orthogonal_axis_disagreement,
            "structure_label_same_axis_conflict": bool(
                structure_label_same_axis_conflicts
            ),
            "structure_label_same_axis_conflict_axes": ";".join(
                structure_label_same_axis_conflicts
            ),
            "true_same_axis_conflict": bool(axis_conflicts),
            "same_axis_conflict_axes": ";".join(axis_conflicts),
            "axis_conflict_axes": ";".join(axis_conflicts),
            "axis_unscorable_count": axis_unscorable_count,
            "high_confidence_wrong_after_axis_scoring": high_confidence_wrong_after_axis,
            "taxonomy_collapse_reason": _taxonomy_collapse_reason(
                label_class=label_class,
                structure_class=structure_class,
                label_axes=label_axis_values,
                structure_axes=structure_axis_values,
                same_axis_conflicts=structure_label_same_axis_conflicts,
            ),
            "folding_problem_solved": False,
            "folding_solution_claim_created": False,
        }
        for axis in AXIS_NAMES:
            short = axis.replace("_axis", "")
            row[f"predicted_{axis}"] = predicted_axis_values[axis]
            row[f"structure_{axis}"] = structure_axis_values[axis]
            row[f"label_{axis}"] = label_axis_values[axis]
            row[f"adjudicated_truth_{axis}"] = truth_axis_values[axis]
            row[f"{short}_axis_scorable"] = axis_scores[axis]["scorable"]
            row[f"{short}_axis_match"] = axis_scores[axis]["match"]
        manual_reasons = _manual_review_reasons(row)
        row["manual_review_required"] = bool(manual_reasons)
        row["manual_review_reasons"] = ";".join(manual_reasons)
        rows.append(row)
    return rows


def axis_conflict_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "predicted_fold_class",
        "structure_fold_class",
        "label_fold_class",
        "single_class_label_structure_disagreement",
        "orthogonal_axis_disagreement",
        "structure_label_same_axis_conflict",
        "structure_label_same_axis_conflict_axes",
        "true_same_axis_conflict",
        "same_axis_conflict_axes",
        "axis_conflict_axes",
        "forced_prediction",
        "abstained",
        "taxonomy_collapse_reason",
        "manual_review_reasons",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["single_class_label_structure_disagreement"])
        or bool(row["true_same_axis_conflict"])
        or bool(row["axis_conflict_axes"])
    ]


def axis_manual_review_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "predicted_fold_class",
        "structure_fold_class",
        "label_fold_class",
        "orthogonal_axis_disagreement",
        "true_same_axis_conflict",
        "axis_unscorable_count",
        "high_confidence_wrong_after_axis_scoring",
        "manual_review_reasons",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["manual_review_required"])
    ]


def axis_confusion_matrix_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    counts: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in rows:
        for axis in AXIS_NAMES:
            predicted = str(row[f"predicted_{axis}"])
            truth = str(row[f"adjudicated_truth_{axis}"])
            key = (axis, predicted, truth)
            if key not in counts:
                counts[key] = {
                    "axis_name": axis,
                    "predicted_axis_value": predicted,
                    "truth_axis_value": truth,
                    "count": 0,
                    "match_count": 0,
                    "conflict_count": 0,
                    "unscorable_count": 0,
                }
            item = counts[key]
            item["count"] = int(item["count"]) + 1
            scorable = bool(
                row[f"{axis.replace('_axis', '')}_axis_scorable"]
            )
            match = bool(row[f"{axis.replace('_axis', '')}_axis_match"])
            item["match_count"] = int(item["match_count"]) + int(match)
            item["conflict_count"] = int(item["conflict_count"]) + int(
                scorable and not match
            )
            item["unscorable_count"] = int(item["unscorable_count"]) + int(
                not scorable
            )
    return sorted(
        counts.values(),
        key=lambda item: (
            str(item["axis_name"]),
            str(item["predicted_axis_value"]),
            str(item["truth_axis_value"]),
        ),
    )


def _axis_accuracy(
    rows: Sequence[Mapping[str, object]],
    *,
    axis: str,
) -> dict[str, object]:
    prefix = axis.replace("_axis", "")
    scorable_rows = [
        row for row in rows if bool(row[f"{prefix}_axis_scorable"])
    ]
    return {
        "scorable_count": len(scorable_rows),
        "unscorable_count": len(rows) - len(scorable_rows),
        "accuracy": _bool_mean(
            bool(row[f"{prefix}_axis_match"]) for row in scorable_rows
        ),
    }


def _taxonomy_note_counts(rows: Sequence[Mapping[str, object]]) -> dict[str, int]:
    return {
        "architecture_not_secondary_structure": sum(
            1
            for row in rows
            if row["label_architecture_axis"] != "unknown"
            and row["structure_secondary_structure_axis"] != "weak_or_unknown"
        ),
        "membrane_regime_not_fold_class": sum(
            1 for row in rows if row["predicted_environment_axis"] == "membrane_like"
        ),
        "disorder_not_beta_alpha_absence": sum(
            1
            for row in rows
            if row["label_order_axis"] == "disordered_flexible"
            or row["structure_order_axis"] == "disordered_flexible"
        ),
        "fragment_evidence_not_global_fold_truth": sum(
            1
            for row in rows
            if row["predicted_architecture_axis"] == "fragment_scope"
            or row["structure_architecture_axis"] == "fragment_scope"
        ),
    }


def _axis_conflict_count(
    rows: Sequence[Mapping[str, object]],
    *,
    axis_name: Optional[str] = None,
    forced: Optional[bool] = None,
    abstained: Optional[bool] = None,
) -> int:
    count = 0
    for row in rows:
        conflict_axes = {
            axis for axis in str(row["axis_conflict_axes"]).split(";") if axis
        }
        if axis_name is not None and axis_name not in conflict_axes:
            continue
        if axis_name is None and not conflict_axes:
            continue
        if forced is not None and bool(row["forced_prediction"]) is not forced:
            continue
        if abstained is not None and bool(row["abstained"]) is not abstained:
            continue
        count += 1
    return count


def _guard_abstention_count(rows: Sequence[Mapping[str, object]], marker: str) -> int:
    return sum(
        1
        for row in rows
        if bool(row["abstained"]) and marker in str(row["gate_path"])
    )


def build_axis_adjudication_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = axis_adjudication_rows(references, evidence_rows)
    structure_label_disagreements = [
        row for row in rows if bool(row["single_class_label_structure_disagreement"])
    ]
    orthogonal_axis_disagreements = [
        row for row in rows if bool(row["orthogonal_axis_disagreement"])
    ]
    true_same_axis_conflicts = [
        row for row in rows if bool(row["true_same_axis_conflict"])
    ]
    structure_label_same_axis_conflicts = [
        row for row in rows if bool(row["structure_label_same_axis_conflict"])
    ]
    high_confidence_wrong_after_axis = [
        row for row in rows if bool(row["high_confidence_wrong_after_axis_scoring"])
    ]
    axis_unscorable_count = sum(int(row["axis_unscorable_count"]) for row in rows)
    folded_domain_mimic_abstained_count = _guard_abstention_count(
        rows,
        "abstained_folded_domain_mimic_disorder_conflict",
    )
    secondary_axis_ambiguity_abstained_count = _guard_abstention_count(
        rows,
        "secondary_structure_gate:abstained_alpha_mixed_ambiguity",
    )
    return {
        "benchmark_kind": FOLD_AXIS_ADJUDICATION_BENCHMARK_KIND,
        "source_regime_analysis_benchmark_kind": REGIME_ANALYSIS_BENCHMARK_KIND,
        "axis_signature_kind": AXIS_SIGNATURE_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "predictor_input_boundary": "sequence_only_no_labels_no_structure_answers",
        "truth_adjudication_boundary": (
            "labels_structure_sources_and_reference_axes_used_only_after_prediction"
        ),
        "benchmark_size": len(rows),
        "single_class_taxonomy_collapse_detected": bool(
            orthogonal_axis_disagreements
        ),
        "structure_label_disagreement_count": len(structure_label_disagreements),
        "orthogonal_axis_disagreement_count": len(orthogonal_axis_disagreements),
        "true_same_axis_conflict_count": len(true_same_axis_conflicts),
        "structure_label_same_axis_conflict_count": len(
            structure_label_same_axis_conflicts
        ),
        "axis_unscorable_count": axis_unscorable_count,
        "axis_unscorable_row_count": sum(
            1 for row in rows if int(row["axis_unscorable_count"]) > 0
        ),
        "high_confidence_wrong_count_after_axis_scoring": len(
            high_confidence_wrong_after_axis
        ),
        "forced_same_axis_conflict_count": _axis_conflict_count(
            rows,
            forced=True,
        ),
        "forced_order_axis_conflict_count": _axis_conflict_count(
            rows,
            axis_name="order_axis",
            forced=True,
        ),
        "forced_secondary_axis_conflict_count": _axis_conflict_count(
            rows,
            axis_name="secondary_structure_axis",
            forced=True,
        ),
        "abstained_axis_conflict_count": _axis_conflict_count(
            rows,
            abstained=True,
        ),
        "regime_axis_conflict_count": (
            _axis_conflict_count(rows, axis_name="architecture_axis")
            + _axis_conflict_count(rows, axis_name="environment_axis")
        ),
        "folded_domain_mimic_abstained_count": (
            folded_domain_mimic_abstained_count
        ),
        "secondary_axis_ambiguity_abstained_count": (
            secondary_axis_ambiguity_abstained_count
        ),
        "coverage_loss_from_safety_guards": (
            folded_domain_mimic_abstained_count
            + secondary_axis_ambiguity_abstained_count
        ),
        "artifact_reproducible": True,
        "forced_prediction_count": sum(
            1 for row in rows if bool(row["forced_prediction"])
        ),
        "abstained_prediction_count": sum(1 for row in rows if bool(row["abstained"])),
        "axis_accuracy": {
            axis: _axis_accuracy(rows, axis=axis) for axis in AXIS_NAMES
        },
        "taxonomy_note_counts": _taxonomy_note_counts(rows),
        "manual_review_row_count": sum(
            1 for row in rows if bool(row["manual_review_required"])
        ),
        "manual_review_rows": axis_manual_review_rows(rows),
        "revision_required": True,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "boundary_statement": (
            "This layer adjudicates the existing sequence-only regime-routed "
            "predictions against orthogonal truth axes after prediction. It does "
            "not retune routing, generate new protein predictions, export raw "
            "sequences, or claim that folding is solved."
        ),
        "rows": rows,
    }


def write_axis_adjudication_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    conflicts: Sequence[Mapping[str, object]],
    manual_review_rows: Sequence[Mapping[str, object]],
    confusion_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    conflicts_path: Path,
    manual_review_path: Path,
    confusion_matrices_path: Path,
    dashboard_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(rows, rows_path)
    _write_csv_rows(conflicts, conflicts_path)
    _write_csv_rows(manual_review_rows, manual_review_path)
    _write_csv_rows(confusion_rows, confusion_matrices_path)
    dashboard_path.write_text(render_axis_dashboard(report), encoding="utf-8")
    return (
        report_path,
        rows_path,
        conflicts_path,
        manual_review_path,
        confusion_matrices_path,
        dashboard_path,
    )


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
) -> Path:
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
        "single_class_taxonomy_collapse_detected",
        "structure_label_disagreement_count",
        "orthogonal_axis_disagreement_count",
        "true_same_axis_conflict_count",
        "forced_same_axis_conflict_count",
        "axis_unscorable_count",
        "high_confidence_wrong_count_after_axis_scoring",
        "claim_allowed",
        "folding_problem_solved",
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
        if isinstance(value, Mapping):
            rendered = ", ".join(
                f"{_escape(nested_key)}: {_escape(nested_value)}"
                for nested_key, nested_value in value.items()
            )
        else:
            rendered = _escape(value)
        rows.append(
            "<tr>"
            f"<td>{_escape(key)}</td>"
            f"<td>{rendered}</td>"
            "</tr>"
        )
    return (
        f"<section><h2>{_escape(title)}</h2>"
        "<table><thead><tr><th>key</th><th>value</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def _axis_confusion_preview(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        for axis in AXIS_NAMES:
            predicted = str(row[f"predicted_{axis}"])
            truth = str(row[f"adjudicated_truth_{axis}"])
            counts[axis][f"{predicted} -> {truth}"] += 1
    sections = []
    for axis in AXIS_NAMES:
        body = "".join(
            "<tr>"
            f"<td>{_escape(pair)}</td><td>{_escape(count)}</td>"
            "</tr>"
            for pair, count in counts[axis].most_common(10)
        )
        sections.append(
            f"<h3>{_escape(axis)}</h3>"
            "<table><thead><tr><th>predicted -> truth</th><th>count</th>"
            "</tr></thead><tbody>"
            + body
            + "</tbody></table>"
        )
    return "<section><h2>Axis Confusion Matrices</h2>" + "".join(sections) + "</section>"


def _manual_review_table(report: Mapping[str, object]) -> str:
    rows = report.get("manual_review_rows", [])
    if not isinstance(rows, Sequence) or not rows:
        return "<section><h2>Manual Review Rows</h2><p>No rows.</p></section>"
    body = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        body.append(
            "<tr>"
            f"<td>{_escape(row.get('protein_id', ''))}</td>"
            f"<td>{_escape(row.get('protein_regime', ''))}</td>"
            f"<td>{_escape(row.get('predicted_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('structure_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('label_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('manual_review_reasons', ''))}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Manual Review Rows</h2>"
        "<table><thead><tr><th>protein</th><th>regime</th><th>predicted</th>"
        "<th>structure</th><th>label</th><th>reasons</th></tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def render_axis_dashboard(report: Mapping[str, object]) -> str:
    title = "Fold Axis Truth Adjudication"
    callouts = (
        "Single-Class Benchmark Is Lossy",
        "Architecture ≠ Secondary Structure",
        "Membrane Regime ≠ Fold Class",
        "Disorder ≠ Beta/Alpha Absence",
        "Fragment Evidence ≠ Global Fold Truth",
    )
    return (
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        f"<title>{_escape(title)}</title>"
        "<style>"
        "body{font-family:Arial,sans-serif;margin:0;background:#f7f8f5;color:#1f2933;}"
        "header{background:#263238;color:white;padding:28px 32px;}"
        "main{padding:24px 32px;}"
        ".metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:18px 0;}"
        ".metric,.callout{background:white;border:1px solid #d7ddd8;border-radius:6px;padding:14px;}"
        ".metric span{display:block;font-size:12px;color:#52616b;margin-bottom:6px;}"
        ".metric strong{font-size:20px;}"
        ".callouts{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:10px;margin:18px 0;}"
        "section{margin:24px 0;}"
        "table{border-collapse:collapse;width:100%;background:white;border:1px solid #d7ddd8;}"
        "th,td{border:1px solid #d7ddd8;padding:8px;text-align:left;vertical-align:top;font-size:13px;}"
        "th{background:#e9efec;}"
        "</style></head><body>"
        f"<header><h1>{_escape(title)}</h1>"
        "<p>Collapsed fold classes are split into secondary-structure, "
        "architecture, order, and environment axes after prediction.</p></header>"
        "<main>"
        f"<div class=\"metrics\">{_metric_cards(report)}</div>"
        "<div class=\"callouts\">"
        + "".join(f"<div class=\"callout\">{_escape(callout)}</div>" for callout in callouts)
        + "</div>"
        + _mapping_table("Axis Accuracy", report.get("axis_accuracy"))
        + _mapping_table("Taxonomy Notes", report.get("taxonomy_note_counts"))
        + _axis_confusion_preview(report)
        + _manual_review_table(report)
        + "</main></body></html>"
    )
