from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_axis_adjudication import (
    AXIS_NAMES,
    AXIS_SIGNATURE_KIND,
    UNKNOWN_BY_AXIS,
    axis_adjudication_rows,
)
from pharmacotopology.folding_regime_analysis import REGIME_ANALYSIS_BENCHMARK_KIND
from pharmacotopology.folding_structure_benchmark import StructureEvidenceRow
from pharmacotopology.folding_topology import FoldingReferenceExample


FOLD_AXIS_PROFILE_BENCHMARK_KIND = "fold_axis_profile_coverage_recovery"
AXIS_PROFILE_SIGNATURE_KIND = "axis_safe_partial_fold_profile"
GLOBAL_PROFILE_CLASS = "insufficient_topology_evidence"


def _bool_mean(values: Sequence[bool]) -> float:
    if not values:
        return 0.0
    return round(sum(1 for value in values if value) / len(values), 6)


def _known_axis_value(axis: str, value: object) -> bool:
    return str(value) != UNKNOWN_BY_AXIS[axis]


def _axis_scorable(axis: str, predicted_value: object, truth_value: object) -> bool:
    return _known_axis_value(axis, predicted_value) and _known_axis_value(
        axis,
        truth_value,
    )


def _axis_match(axis: str, predicted_value: object, truth_value: object) -> bool:
    if not _axis_scorable(axis, predicted_value, truth_value):
        return False
    return str(predicted_value) == str(truth_value)


def _empty_profile_axes() -> dict[str, str]:
    return {axis: UNKNOWN_BY_AXIS[axis] for axis in AXIS_NAMES}


def _source_class_axis_profile(
    row: Mapping[str, object],
) -> tuple[dict[str, str], dict[str, str]]:
    axes = _empty_profile_axes()
    reasons = {axis: "" for axis in AXIS_NAMES}
    source_class = str(row["predicted_fold_class"])
    if not bool(row["forced_prediction"]):
        return axes, reasons

    if source_class in {"alpha_rich", "beta_rich", "alpha_beta_mixed"}:
        axes["secondary_structure_axis"] = source_class
        axes["order_axis"] = "ordered"
        reasons["secondary_structure_axis"] = (
            "source_forced_class_projected_to_secondary_axis"
        )
        reasons["order_axis"] = "source_forced_folded_class_projected_to_order_axis"
    elif source_class == "disordered_flexible":
        axes["order_axis"] = "disordered_flexible"
        reasons["order_axis"] = "source_forced_disorder_class_projected_to_order_axis"
    elif source_class == "multidomain_boundary":
        axes["order_axis"] = "ordered"
        reasons["order_axis"] = (
            "source_forced_architecture_class_projected_to_order_axis"
        )
    return axes, reasons


def _recover_secondary_axis(
    row: Mapping[str, object],
    axes: dict[str, str],
    reasons: dict[str, str],
) -> None:
    if _known_axis_value("secondary_structure_axis", axes["secondary_structure_axis"]):
        return
    gate_path = str(row["gate_path"])
    if "secondary_structure_gate:alpha_periodic_compact" in gate_path:
        axes["secondary_structure_axis"] = "alpha_rich"
        reasons["secondary_structure_axis"] = (
            "recovered_from_specific_alpha_periodic_secondary_gate"
        )
    elif "secondary_structure_gate:beta_pairing_supported" in gate_path:
        axes["secondary_structure_axis"] = "beta_rich"
        reasons["secondary_structure_axis"] = (
            "recovered_from_specific_beta_pairing_secondary_gate"
        )


def _recover_order_axis(
    row: Mapping[str, object],
    axes: dict[str, str],
    reasons: dict[str, str],
) -> None:
    if _known_axis_value("order_axis", axes["order_axis"]):
        return
    gate_path = str(row["gate_path"])
    protein_regime = str(row["protein_regime"])
    if "regime_router:abstained_folded_domain_mimic_disorder_conflict" in gate_path:
        axes["order_axis"] = "ordered"
        reasons["order_axis"] = "recovered_from_folded_domain_mimic_guard"
    elif "secondary_structure_gate:abstained_alpha_mixed_ambiguity" in gate_path:
        axes["order_axis"] = "ordered"
        reasons["order_axis"] = "recovered_from_alpha_mixed_ambiguity_guard"
    elif (
        protein_regime not in {"intrinsically_disordered", "ambiguous_regime"}
        and "disorder_gate:foldable_candidate" in gate_path
        and (
            "compactness_gate:compact_or_borderline_supported" in gate_path
            or "compactness_gate:compact_closure_supported" in gate_path
        )
    ):
        axes["order_axis"] = "ordered"
        reasons["order_axis"] = "recovered_from_foldable_compact_axis_evidence"


def _recover_architecture_axis(
    row: Mapping[str, object],
    axes: dict[str, str],
    reasons: dict[str, str],
) -> None:
    if str(row["protein_regime"]) == "small_peptide_or_fragment":
        axes["architecture_axis"] = "fragment_scope"
        reasons["architecture_axis"] = "recovered_fragment_scope_from_regime"


def _recover_environment_axis(
    row: Mapping[str, object],
    axes: dict[str, str],
    reasons: dict[str, str],
) -> None:
    if str(row["protein_regime"]) == "membrane_like":
        axes["environment_axis"] = "membrane_like"
        reasons["environment_axis"] = "recovered_membrane_axis_from_regime"


def _profile_axes_for_row(
    row: Mapping[str, object],
) -> tuple[dict[str, str], dict[str, str]]:
    axes, reasons = _source_class_axis_profile(row)
    _recover_secondary_axis(row, axes, reasons)
    _recover_order_axis(row, axes, reasons)
    _recover_architecture_axis(row, axes, reasons)
    _recover_environment_axis(row, axes, reasons)
    return axes, reasons


def _guard_override_created(
    row: Mapping[str, object],
    profile_axes: Mapping[str, str],
) -> bool:
    gate_path = str(row["gate_path"])
    if "secondary_structure_gate:abstained_alpha_mixed_ambiguity" in gate_path:
        return _known_axis_value(
            "secondary_structure_axis",
            profile_axes["secondary_structure_axis"],
        )
    if "regime_router:abstained_folded_domain_mimic_disorder_conflict" in gate_path:
        return (
            profile_axes["order_axis"] == "disordered_flexible"
            or _known_axis_value(
                "secondary_structure_axis",
                profile_axes["secondary_structure_axis"],
            )
        )
    return False


def axis_profile_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    source_rows = axis_adjudication_rows(references, evidence_rows)
    rows: list[dict[str, object]] = []
    for source_row in source_rows:
        profile_axes, profile_reasons = _profile_axes_for_row(source_row)
        conflict_axes = []
        scorable_count = 0
        match_count = 0
        safe_recovered_axes = []
        row: dict[str, object] = {
            "protein_id": source_row["protein_id"],
            "sequence_length": source_row["sequence_length"],
            "axis_profile_signature_kind": AXIS_PROFILE_SIGNATURE_KIND,
            "source_axis_signature_kind": AXIS_SIGNATURE_KIND,
            "protein_regime": source_row["protein_regime"],
            "source_predicted_fold_class": source_row["predicted_fold_class"],
            "profile_global_fold_class": GLOBAL_PROFILE_CLASS,
            "source_forced_prediction": source_row["forced_prediction"],
            "source_abstained": source_row["abstained"],
            "confidence": source_row["confidence"],
            "gate_path": source_row["gate_path"],
            "gate_decision_reason": source_row["gate_decision_reason"],
            "global_class_claim_allowed": False,
            "collapsed_class_recovered": False,
            "unsafe_class_recovery": False,
        }
        for axis in AXIS_NAMES:
            profile_value = profile_axes[axis]
            truth_value = source_row[f"adjudicated_truth_{axis}"]
            scorable = _axis_scorable(axis, profile_value, truth_value)
            match = _axis_match(axis, profile_value, truth_value)
            if scorable:
                scorable_count += 1
                match_count += int(match)
                if not match:
                    conflict_axes.append(axis)
            if bool(source_row["abstained"]) and _known_axis_value(
                axis,
                profile_value,
            ):
                safe_recovered_axes.append(axis)
            prefix = axis.replace("_axis", "")
            row[f"profile_{axis}"] = profile_value
            row[f"profile_{prefix}_claim_allowed"] = _known_axis_value(
                axis,
                profile_value,
            )
            row[f"profile_{prefix}_claim_reason"] = profile_reasons[axis]
            row[f"adjudicated_truth_{axis}"] = truth_value
            row[f"{prefix}_axis_scorable"] = scorable
            row[f"{prefix}_axis_match"] = match
        row["axis_profile_claim_count"] = sum(
            1 for axis in AXIS_NAMES if _known_axis_value(axis, profile_axes[axis])
        )
        row["axis_profile_has_claim"] = int(row["axis_profile_claim_count"]) > 0
        row["axis_profile_same_axis_conflict"] = bool(conflict_axes)
        row["axis_profile_conflict_axes"] = ";".join(conflict_axes)
        row["axis_profile_scorable_count"] = scorable_count
        row["axis_profile_match_count"] = match_count
        row["safe_axis_recovered_axes"] = ";".join(safe_recovered_axes)
        row["safe_axis_recovered_count"] = len(safe_recovered_axes)
        row["guard_override"] = _guard_override_created(source_row, profile_axes)
        row["folding_problem_solved"] = False
        row["folding_solution_claim_created"] = False
        row["drug_design_created"] = False
        row["molecule_generated"] = False
        row["protein_sequence_design_created"] = False
        rows.append(row)
    return rows


def axis_profile_abstention_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "source_predicted_fold_class",
        "profile_global_fold_class",
        "profile_secondary_structure_axis",
        "profile_architecture_axis",
        "profile_order_axis",
        "profile_environment_axis",
        "safe_axis_recovered_axes",
        "safe_axis_recovered_count",
        "axis_profile_same_axis_conflict",
        "axis_profile_conflict_axes",
        "gate_path",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["source_abstained"])
    ]


def axis_profile_recovery_candidate_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "profile_secondary_structure_axis",
        "profile_architecture_axis",
        "profile_order_axis",
        "profile_environment_axis",
        "safe_axis_recovered_axes",
        "safe_axis_recovered_count",
        "axis_profile_same_axis_conflict",
        "profile_secondary_structure_claim_reason",
        "profile_architecture_claim_reason",
        "profile_order_claim_reason",
        "profile_environment_claim_reason",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if int(row["safe_axis_recovered_count"]) > 0
    ]


def _coverage(rows: Sequence[Mapping[str, object]], axis: str) -> float:
    return _bool_mean(
        [
            _known_axis_value(axis, row[f"profile_{axis}"])
            for row in rows
        ]
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
        "accuracy": _bool_mean(
            [bool(row[f"{prefix}_axis_match"]) for row in scorable_rows]
        ),
        "scorable_count": len(scorable_rows),
        "claimed_count": sum(
            1 for row in rows if _known_axis_value(axis, row[f"profile_{axis}"])
        ),
        "unscorable_count": len(rows) - len(scorable_rows),
    }


def _recovery_distribution(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        for axis in str(row["safe_axis_recovered_axes"]).split(";"):
            if axis:
                counts[axis] += 1
    return {axis: counts[axis] for axis in AXIS_NAMES}


def build_axis_profile_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = axis_profile_rows(references, evidence_rows)
    forced_count = sum(1 for row in rows if bool(row["source_forced_prediction"]))
    unsafe_class_recovery_count = sum(
        1 for row in rows if bool(row["unsafe_class_recovery"])
    )
    guard_override_count = sum(1 for row in rows if bool(row["guard_override"]))
    axis_profile_same_axis_conflict_count = sum(
        1 for row in rows if bool(row["axis_profile_same_axis_conflict"])
    )
    high_confidence_wrong_after_axis = sum(
        1
        for row in rows
        if bool(row["source_forced_prediction"])
        and float(row["confidence"]) >= 0.58
        and bool(row["axis_profile_same_axis_conflict"])
    )
    safe_axis_recovered_count = sum(
        int(row["safe_axis_recovered_count"]) for row in rows
    )
    return {
        "benchmark_kind": FOLD_AXIS_PROFILE_BENCHMARK_KIND,
        "source_regime_analysis_benchmark_kind": REGIME_ANALYSIS_BENCHMARK_KIND,
        "axis_profile_signature_kind": AXIS_PROFILE_SIGNATURE_KIND,
        "source_axis_signature_kind": AXIS_SIGNATURE_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "predictor_input_boundary": "sequence_only_no_labels_no_structure_answers",
        "truth_scoring_boundary": (
            "labels_structure_sources_and_reference_axes_used_only_after_axis_profile"
        ),
        "benchmark_size": len(rows),
        "collapsed_class_coverage": _bool_mean(
            [bool(row["source_forced_prediction"]) for row in rows]
        ),
        "axis_profile_coverage": _bool_mean(
            [bool(row["axis_profile_has_claim"]) for row in rows]
        ),
        "secondary_axis_coverage": _coverage(rows, "secondary_structure_axis"),
        "architecture_axis_coverage": _coverage(rows, "architecture_axis"),
        "order_axis_coverage": _coverage(rows, "order_axis"),
        "environment_axis_coverage": _coverage(rows, "environment_axis"),
        "safe_axis_recovered_count": safe_axis_recovered_count,
        "safe_axis_recovered_row_count": sum(
            1 for row in rows if int(row["safe_axis_recovered_count"]) > 0
        ),
        "safe_axis_recovered_distribution": _recovery_distribution(rows),
        "unsafe_class_recovery_count": unsafe_class_recovery_count,
        "guard_override_count": guard_override_count,
        "forced_same_axis_conflict_count": sum(
            1
            for row in rows
            if bool(row["source_forced_prediction"])
            and bool(row["axis_profile_same_axis_conflict"])
        ),
        "axis_profile_same_axis_conflict_count": (
            axis_profile_same_axis_conflict_count
        ),
        "high_confidence_wrong_count_after_axis_scoring": (
            high_confidence_wrong_after_axis
        ),
        "forced_prediction_count": forced_count,
        "abstained_prediction_count": sum(
            1 for row in rows if bool(row["source_abstained"])
        ),
        "axis_accuracy": {
            axis: _axis_accuracy(rows, axis=axis) for axis in AXIS_NAMES
        },
        "global_fold_class_claim_allowed": False,
        "axis_profile_claim_allowed": (
            axis_profile_same_axis_conflict_count == 0
            and high_confidence_wrong_after_axis == 0
            and unsafe_class_recovery_count == 0
            and guard_override_count == 0
            and safe_axis_recovered_count > 0
        ),
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "artifact_reproducible": True,
        "boundary_statement": (
            "This layer recovers only axis-level profile claims from existing "
            "sequence-only gates. It does not recover collapsed fold-class "
            "coverage, override abstention guards, export raw sequences, or "
            "claim that folding is solved."
        ),
        "rows": rows,
    }


def build_axis_profile_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": "fold_axis_profile_safety_certificate",
        "benchmark_kind": report["benchmark_kind"],
        "axis_profile_signature_kind": report["axis_profile_signature_kind"],
        "source_benchmark_file": report["source_benchmark_file"],
        "structure_evidence_file": report["structure_evidence_file"],
        "benchmark_size": report["benchmark_size"],
        "collapsed_class_coverage": report["collapsed_class_coverage"],
        "axis_profile_coverage": report["axis_profile_coverage"],
        "safe_axis_recovered_count": report["safe_axis_recovered_count"],
        "unsafe_class_recovery_count": report["unsafe_class_recovery_count"],
        "guard_override_count": report["guard_override_count"],
        "axis_profile_same_axis_conflict_count": report[
            "axis_profile_same_axis_conflict_count"
        ],
        "high_confidence_wrong_count_after_axis_scoring": report[
            "high_confidence_wrong_count_after_axis_scoring"
        ],
        "global_fold_class_claim_allowed": report[
            "global_fold_class_claim_allowed"
        ],
        "axis_profile_claim_allowed": report["axis_profile_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "raw_sequences_exported": False,
        "output_artifacts": tuple(output_names),
    }


def write_axis_profile_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    abstention_rows: Sequence[Mapping[str, object]],
    recovery_candidate_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    abstentions_path: Path,
    recovery_candidates_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    output_names = (
        report_path.name,
        rows_path.name,
        abstentions_path.name,
        recovery_candidates_path.name,
        dashboard_path.name,
        certificate_path.name,
    )
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(rows, rows_path)
    _write_csv_rows(abstention_rows, abstentions_path)
    _write_csv_rows(recovery_candidate_rows, recovery_candidates_path)
    dashboard_path.write_text(render_axis_profile_dashboard(report), encoding="utf-8")
    certificate = build_axis_profile_certificate(report, output_names=output_names)
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        abstentions_path,
        recovery_candidates_path,
        dashboard_path,
        certificate_path,
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
        "collapsed_class_coverage",
        "axis_profile_coverage",
        "safe_axis_recovered_count",
        "unsafe_class_recovery_count",
        "guard_override_count",
        "axis_profile_same_axis_conflict_count",
        "high_confidence_wrong_count_after_axis_scoring",
        "axis_profile_claim_allowed",
        "global_fold_class_claim_allowed",
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


def _recovery_preview(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    preview_rows = [
        row
        for row in rows
        if isinstance(row, Mapping) and int(row["safe_axis_recovered_count"]) > 0
    ][:20]
    body = "".join(
        "<tr>"
        f"<td>{_escape(row['protein_id'])}</td>"
        f"<td>{_escape(row['protein_regime'])}</td>"
        f"<td>{_escape(row['profile_secondary_structure_axis'])}</td>"
        f"<td>{_escape(row['profile_architecture_axis'])}</td>"
        f"<td>{_escape(row['profile_order_axis'])}</td>"
        f"<td>{_escape(row['profile_environment_axis'])}</td>"
        f"<td>{_escape(row['safe_axis_recovered_axes'])}</td>"
        "</tr>"
        for row in preview_rows
    )
    return (
        "<section><h2>Axis-Safe Recovery Candidates</h2>"
        "<table><thead><tr>"
        "<th>protein_id</th><th>regime</th><th>secondary</th>"
        "<th>architecture</th><th>order</th><th>environment</th><th>recovered axes</th>"
        "</tr></thead><tbody>"
        + body
        + "</tbody></table></section>"
    )


def render_axis_profile_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Fold Axis Profile Coverage Recovery</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f7f4ee;
      color: #202326;
    }}
    header {{
      padding: 32px;
      background: #24322f;
      color: #f8f4ec;
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
      margin: 20px 0 0;
    }}
    .metric {{
      background: #ffffff;
      border: 1px solid #d8d1c5;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #5f665f;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 24px;
    }}
    .rule-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}
    .rule {{
      background: #fff;
      border: 1px solid #d8d1c5;
      border-radius: 6px;
      padding: 14px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #fff;
      border: 1px solid #d8d1c5;
      border-radius: 6px;
      overflow: hidden;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #ebe5da;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
    th {{
      background: #eee7dc;
      color: #39423d;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Fold Axis Profile Coverage Recovery</h1>
    <p>Partial axis claims are allowed only where existing sequence-only gates support them. Global fold-class recovery remains refused.</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Safety Rules</h2>
      <div class="rule-grid">
        <div class="rule"><strong>No Collapsed Class Recovery</strong><br>Abstained rows stay globally unknown.</div>
        <div class="rule"><strong>I Know This Axis</strong><br>Specific secondary, order, scope, or membrane evidence can populate one axis.</div>
        <div class="rule"><strong>I Do Not Know That Axis</strong><br>Architecture stays mostly unknown until its own evidence layer exists.</div>
        <div class="rule"><strong>No Guard Overrides</strong><br>Safety abstention guards are not turned back into class predictions.</div>
      </div>
    </section>
    {_mapping_table("Axis Coverage", {
        "secondary_axis_coverage": report.get("secondary_axis_coverage"),
        "architecture_axis_coverage": report.get("architecture_axis_coverage"),
        "order_axis_coverage": report.get("order_axis_coverage"),
        "environment_axis_coverage": report.get("environment_axis_coverage"),
    })}
    {_mapping_table("Axis Accuracy", report.get("axis_accuracy"))}
    {_mapping_table("Recovery Distribution", report.get("safe_axis_recovered_distribution"))}
    {_recovery_preview(report)}
  </main>
</body>
</html>
"""
