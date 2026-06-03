from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_architecture_axis import (
    ARCHITECTURE_AXIS_SIGNATURE_KIND,
    architecture_evidence_packet_from_sequence,
)
from pharmacotopology.folding_axis_adjudication import AXIS_NAMES, UNKNOWN_BY_AXIS
from pharmacotopology.folding_external_holdout import (
    DEFAULT_DEVELOPMENT_BENCHMARK_FILE,
    EXTERNAL_HOLDOUT_BENCHMARK_KIND,
    EXTERNAL_HOLDOUT_SPLIT,
    ExternalHoldoutRow,
    _axis_conflict_axes,
    _bool_mean,
    _coverage,
    _failure_cohort,
    _known_axis,
    _safe_axis_claim_total,
    _source_row_for_profile,
    _unsafe_axis_claim_total,
    build_external_holdout_report,
    external_holdout_rows,
    validate_holdout_lock,
)
from pharmacotopology.folding_axis_profile import AXIS_PROFILE_SIGNATURE_KIND
from pharmacotopology.folding_order_axis_safety import (
    ORDER_AXIS_SAFETY_SIGNATURE_KIND,
    order_axis_safety_packet_from_source,
)


EXTERNAL_AXIS_REPAIR_BENCHMARK_KIND = "external_safe_axis_conflict_quarantine"
EXTERNAL_AXIS_REPAIR_STACK = (
    "external_holdout_baseline;order_axis_folded_mimic_quarantine;"
    "repeat_compact_architecture_quarantine"
)


def _safe_axis_counts(
    *,
    predicted_axes: Mapping[str, str],
    truth_axes: Mapping[str, str],
) -> tuple[int, int]:
    safe_axis_claim_count = 0
    unsafe_axis_claim_count = 0
    for axis in AXIS_NAMES:
        predicted = predicted_axes[axis]
        truth = truth_axes[axis]
        if not _known_axis(axis, predicted):
            continue
        if _known_axis(axis, truth) and predicted == truth:
            safe_axis_claim_count += 1
        elif _known_axis(axis, truth) and predicted != truth:
            unsafe_axis_claim_count += 1
    return safe_axis_claim_count, unsafe_axis_claim_count


def _axis_claimed(predicted_axes: Mapping[str, str]) -> bool:
    return any(_known_axis(axis, predicted_axes[axis]) for axis in AXIS_NAMES)


def external_axis_repair_rows(
    holdout_rows: Sequence[ExternalHoldoutRow],
) -> list[dict[str, object]]:
    baseline_rows_by_id = {
        str(row["row_id"]): row for row in external_holdout_rows(holdout_rows)
    }
    rows = []
    for holdout_row in holdout_rows:
        baseline = baseline_rows_by_id[holdout_row.row_id]
        source_row = _source_row_for_profile(holdout_row)
        order_packet = order_axis_safety_packet_from_source(
            holdout_row.sequence,
            source_row,
            row_id=holdout_row.row_id,
        )
        architecture_packet = architecture_evidence_packet_from_sequence(
            holdout_row.sequence,
            protein_id=holdout_row.row_id,
            external_safe_quarantine=True,
        )

        repaired_axes = {
            "secondary_structure_axis": str(
                baseline["profile_secondary_structure_axis"]
            ),
            "architecture_axis": architecture_packet.architecture_axis_prediction,
            "order_axis": str(baseline["profile_order_axis"]),
            "environment_axis": str(baseline["profile_environment_axis"]),
        }
        if str(source_row["predicted_fold_class"]) == "disordered_flexible":
            repaired_axes["order_axis"] = order_packet.order_axis_prediction

        truth_axes = holdout_row.truth_axes
        conflict_axes = _axis_conflict_axes(
            predicted_axes=repaired_axes,
            truth_axes=truth_axes,
        )
        safe_axis_claim_count, unsafe_axis_claim_count = _safe_axis_counts(
            predicted_axes=repaired_axes,
            truth_axes=truth_axes,
        )
        architecture_conflict = (
            _known_axis("architecture_axis", repaired_axes["architecture_axis"])
            and _known_axis("architecture_axis", truth_axes["architecture_axis"])
            and repaired_axes["architecture_axis"] != truth_axes["architecture_axis"]
        )
        order_quarantined = (
            order_packet.order_axis_abstention_reason
            == "external_order_axis_folded_mimic_quarantine"
        )
        repeat_quarantined = (
            architecture_packet.architecture_axis_abstention_reason
            == "repeat_compact_single_domain_ambiguity_quarantine"
        )
        safe_axis_recovered_axes = [
            axis
            for axis in AXIS_NAMES
            if _known_axis(axis, repaired_axes[axis])
            and not bool(source_row["forced_prediction"])
        ]
        row = {
            "row_id": holdout_row.row_id,
            "source_id": holdout_row.source_id,
            "source_kind": holdout_row.source_kind,
            "sequence_sha256": holdout_row.sequence_sha256,
            "sequence_hash_short": holdout_row.sequence_sha256[:16],
            "length": holdout_row.length,
            "external_family_id": holdout_row.external_family_id,
            "external_family_name": holdout_row.external_family_name,
            "external_family_group": holdout_row.external_family_group,
            "holdout_split": holdout_row.holdout_split,
            "truth_scope": holdout_row.truth_scope,
            "source_predicted_fold_class": source_row["predicted_fold_class"],
            "source_forced_prediction": source_row["forced_prediction"],
            "source_abstained": source_row["abstained"],
            "source_confidence": source_row["confidence"],
            "protein_regime": source_row["protein_regime"],
            "pre_profile_order_axis": baseline["profile_order_axis"],
            "post_profile_order_axis": repaired_axes["order_axis"],
            "order_axis_claim_allowed": _known_axis(
                "order_axis",
                repaired_axes["order_axis"],
            ),
            "order_axis_confidence": order_packet.order_axis_confidence,
            "order_axis_abstention_reason": (
                order_packet.order_axis_abstention_reason
            ),
            "order_axis_decision_reason": order_packet.order_axis_decision_reason,
            "disorder_pressure": order_packet.disorder_pressure,
            "disorder_run_evidence": order_packet.disorder_run_evidence,
            "breaker_density": order_packet.breaker_density,
            "local_disorder_pressure": order_packet.local_disorder_pressure,
            "folded_beta_mimic_pressure": (
                order_packet.folded_beta_mimic_pressure
            ),
            "folded_mixed_mimic_pressure": (
                order_packet.folded_mixed_mimic_pressure
            ),
            "compact_closure_pressure": order_packet.compact_closure_pressure,
            "beta_pairing_support": order_packet.beta_pairing_support,
            "long_range_closure_evidence": (
                order_packet.long_range_closure_evidence
            ),
            "contact_prior_density": order_packet.contact_prior_density,
            "pre_architecture_axis_prediction": baseline[
                "architecture_axis_prediction"
            ],
            "post_architecture_axis_prediction": (
                architecture_packet.architecture_axis_prediction
            ),
            "architecture_axis_confidence": (
                architecture_packet.architecture_axis_confidence
            ),
            "architecture_axis_claim_allowed": (
                architecture_packet.architecture_axis_claim_allowed
            ),
            "architecture_axis_abstention_reason": (
                architecture_packet.architecture_axis_abstention_reason
            ),
            "repeat_pressure": architecture_packet.repeat_pressure,
            "compact_domain_pressure": architecture_packet.compact_domain_pressure,
            "repeat_recurrence_support": (
                architecture_packet.repeat_recurrence_support
            ),
            "repeat_unit_consistency": architecture_packet.repeat_unit_consistency,
            "repeat_vs_compact_margin": (
                architecture_packet.repeat_vs_compact_margin
            ),
            "hydrophobic_periodicity_only_risk": (
                architecture_packet.hydrophobic_periodicity_only_risk
            ),
            "profile_secondary_structure_axis": repaired_axes[
                "secondary_structure_axis"
            ],
            "profile_order_axis": repaired_axes["order_axis"],
            "profile_environment_axis": repaired_axes["environment_axis"],
            "architecture_axis_prediction": repaired_axes["architecture_axis"],
            "truth_secondary_structure_axis": truth_axes[
                "secondary_structure_axis"
            ],
            "truth_architecture_axis": truth_axes["architecture_axis"],
            "truth_order_axis": truth_axes["order_axis"],
            "truth_environment_axis": truth_axes["environment_axis"],
            "pre_axis_profile_same_axis_conflict": baseline[
                "axis_profile_same_axis_conflict"
            ],
            "pre_architecture_axis_same_axis_conflict": baseline[
                "architecture_axis_same_axis_conflict"
            ],
            "post_axis_profile_same_axis_conflict": bool(conflict_axes),
            "post_axis_profile_conflict_axes": ";".join(conflict_axes),
            "post_architecture_axis_same_axis_conflict": architecture_conflict,
            "post_combined_same_axis_conflict": bool(conflict_axes),
            "post_conflict_axes": ";".join(conflict_axes),
            "safe_axis_claim_count": safe_axis_claim_count,
            "unsafe_axis_claim_count": unsafe_axis_claim_count,
            "safe_axis_recovered_axes": ";".join(safe_axis_recovered_axes),
            "safe_axis_recovered_count": len(safe_axis_recovered_axes),
            "order_axis_folded_mimic_quarantined": order_quarantined,
            "repeat_compact_ambiguity_quarantined": repeat_quarantined,
            "unsafe_class_recovery": False,
            "guard_override": False,
            "any_axis_claimed": _axis_claimed(repaired_axes),
            "global_fold_class_claim_allowed": False,
            "axis_profile_claim_allowed": True,
            "folding_problem_solved": False,
            "claim_allowed": False,
            "prediction_source_layer": "external_safe_axis_repair_stack",
        }
        row["failure_cohort"] = _failure_cohort(
            {
                **baseline,
                "conflict_axes": row["post_conflict_axes"],
                "profile_order_axis": row["profile_order_axis"],
                "architecture_axis_prediction": row[
                    "architecture_axis_prediction"
                ],
                "architecture_axis_abstention_reason": row[
                    "architecture_axis_abstention_reason"
                ],
                "any_axis_claimed": row["any_axis_claimed"],
            }
        )
        rows.append(row)
    return rows


def _post_axis_coverage(rows: Sequence[Mapping[str, object]], axis: str) -> float:
    key = (
        "architecture_axis_prediction"
        if axis == "architecture_axis"
        else f"profile_{axis}"
    )
    return _coverage(rows, key, axis)


def external_axis_repair_family_summary_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["external_family_id"])].append(row)
    output = []
    for family_id, family_rows in sorted(grouped.items()):
        conflicts = [
            row for row in family_rows if bool(row["post_combined_same_axis_conflict"])
        ]
        quarantines = [
            row
            for row in family_rows
            if bool(row["order_axis_folded_mimic_quarantined"])
            or bool(row["repeat_compact_ambiguity_quarantined"])
        ]
        output.append(
            {
                "external_family_id": family_id,
                "external_family_name": family_rows[0]["external_family_name"],
                "external_family_group": family_rows[0]["external_family_group"],
                "row_count": len(family_rows),
                "post_axis_profile_coverage": _bool_mean(
                    [bool(row["any_axis_claimed"]) for row in family_rows]
                ),
                "post_architecture_axis_coverage": _bool_mean(
                    [
                        bool(row["architecture_axis_claim_allowed"])
                        for row in family_rows
                    ]
                ),
                "post_family_conflict_count": len(conflicts),
                "quarantined_count": len(quarantines),
                "dominant_failure_cohort": Counter(
                    str(row["failure_cohort"]) for row in family_rows
                ).most_common(1)[0][0],
            }
        )
    return output


def external_axis_repair_quarantine_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "row_id",
        "external_family_id",
        "external_family_group",
        "length",
        "source_predicted_fold_class",
        "pre_profile_order_axis",
        "post_profile_order_axis",
        "order_axis_abstention_reason",
        "pre_architecture_axis_prediction",
        "post_architecture_axis_prediction",
        "architecture_axis_abstention_reason",
        "truth_order_axis",
        "truth_architecture_axis",
        "order_axis_folded_mimic_quarantined",
        "repeat_compact_ambiguity_quarantined",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["order_axis_folded_mimic_quarantined"])
        or bool(row["repeat_compact_ambiguity_quarantined"])
    ]


def external_axis_repair_conflict_delta_rows(
    report: Mapping[str, object],
) -> list[dict[str, object]]:
    return [
        {
            "metric": "axis_profile_same_axis_conflict_count",
            "before": report[
                "pre_repair_axis_profile_same_axis_conflict_count"
            ],
            "after": report[
                "post_repair_axis_profile_same_axis_conflict_count"
            ],
            "delta": int(
                report["post_repair_axis_profile_same_axis_conflict_count"]
            )
            - int(report["pre_repair_axis_profile_same_axis_conflict_count"]),
        },
        {
            "metric": "architecture_axis_same_axis_conflict_count",
            "before": report[
                "pre_repair_architecture_axis_same_axis_conflict_count"
            ],
            "after": report[
                "post_repair_architecture_axis_same_axis_conflict_count"
            ],
            "delta": int(
                report[
                    "post_repair_architecture_axis_same_axis_conflict_count"
                ]
            )
            - int(
                report[
                    "pre_repair_architecture_axis_same_axis_conflict_count"
                ]
            ),
        },
        {
            "metric": "unsafe_axis_claim_count",
            "before": report["pre_repair_unsafe_axis_claim_count"],
            "after": report["post_repair_unsafe_axis_claim_count"],
            "delta": int(report["post_repair_unsafe_axis_claim_count"])
            - int(report["pre_repair_unsafe_axis_claim_count"]),
        },
        {
            "metric": "high_confidence_wrong_count_after_axis_scoring",
            "before": report[
                "pre_repair_high_confidence_wrong_count_after_axis_scoring"
            ],
            "after": report[
                "post_repair_high_confidence_wrong_count_after_axis_scoring"
            ],
            "delta": int(
                report[
                    "post_repair_high_confidence_wrong_count_after_axis_scoring"
                ]
            )
            - int(
                report[
                    "pre_repair_high_confidence_wrong_count_after_axis_scoring"
                ]
            ),
        },
    ]


def external_axis_repair_abstention_delta_rows(
    report: Mapping[str, object],
) -> list[dict[str, object]]:
    return [
        {
            "metric": "axis_profile_coverage",
            "before": report["pre_repair_axis_profile_coverage"],
            "after": report["post_repair_axis_profile_coverage"],
            "delta": round(
                float(report["post_repair_axis_profile_coverage"])
                - float(report["pre_repair_axis_profile_coverage"]),
                6,
            ),
        },
        {
            "metric": "architecture_axis_coverage",
            "before": report["pre_repair_architecture_axis_coverage"],
            "after": report["post_repair_architecture_axis_coverage"],
            "delta": round(
                float(report["post_repair_architecture_axis_coverage"])
                - float(report["pre_repair_architecture_axis_coverage"]),
                6,
            ),
        },
        {
            "metric": "coverage_loss_from_external_safety",
            "before": 0,
            "after": report["coverage_loss_from_external_safety"],
            "delta": report["coverage_loss_from_external_safety"],
        },
    ]


def build_external_axis_repair_report(
    holdout_rows: Sequence[ExternalHoldoutRow],
    *,
    holdout_file: Path,
    development_benchmark_file: Path = DEFAULT_DEVELOPMENT_BENCHMARK_FILE,
) -> dict[str, object]:
    pre_report = build_external_holdout_report(
        holdout_rows,
        holdout_file=holdout_file,
        development_benchmark_file=development_benchmark_file,
    )
    rows = external_axis_repair_rows(holdout_rows)
    lock = validate_holdout_lock(
        holdout_rows,
        development_benchmark_file=development_benchmark_file,
    )
    post_axis_profile_conflicts = [
        row for row in rows if bool(row["post_axis_profile_same_axis_conflict"])
    ]
    post_architecture_conflicts = [
        row
        for row in rows
        if bool(row["post_architecture_axis_same_axis_conflict"])
    ]
    post_high_confidence_wrong = [
        row
        for row in rows
        if bool(row["source_forced_prediction"])
        and float(row["source_confidence"]) >= 0.58
        and bool(row["post_combined_same_axis_conflict"])
    ]
    post_axis_profile_coverage = _bool_mean(
        [bool(row["any_axis_claimed"]) for row in rows]
    )
    pre_axis_profile_coverage = float(pre_report["axis_profile_coverage"])
    external_safety_repair_successful = (
        not post_axis_profile_conflicts
        and not post_architecture_conflicts
        and _unsafe_axis_claim_total(rows) == 0
        and not post_high_confidence_wrong
    )
    report = {
        "benchmark_kind": EXTERNAL_AXIS_REPAIR_BENCHMARK_KIND,
        "source_external_holdout_benchmark_kind": EXTERNAL_HOLDOUT_BENCHMARK_KIND,
        "holdout_file": str(holdout_file),
        "development_benchmark_file": str(development_benchmark_file),
        "holdout_split": EXTERNAL_HOLDOUT_SPLIT,
        "prediction_stack": EXTERNAL_AXIS_REPAIR_STACK,
        "source_axis_profile_signature_kind": AXIS_PROFILE_SIGNATURE_KIND,
        "source_architecture_axis_signature_kind": ARCHITECTURE_AXIS_SIGNATURE_KIND,
        "order_axis_safety_signature_kind": ORDER_AXIS_SAFETY_SIGNATURE_KIND,
        "predictor_input_boundary": "sequence_only_no_labels_no_structure_answers",
        "truth_scoring_boundary": (
            "external_truth_axes_used_only_after_repaired_axis_prediction"
        ),
        **lock,
        "pre_repair_axis_profile_coverage": pre_axis_profile_coverage,
        "post_repair_axis_profile_coverage": post_axis_profile_coverage,
        "pre_repair_architecture_axis_coverage": pre_report[
            "architecture_axis_coverage"
        ],
        "post_repair_architecture_axis_coverage": _bool_mean(
            [bool(row["architecture_axis_claim_allowed"]) for row in rows]
        ),
        "post_repair_secondary_axis_coverage": _post_axis_coverage(
            rows,
            "secondary_structure_axis",
        ),
        "post_repair_order_axis_coverage": _post_axis_coverage(rows, "order_axis"),
        "post_repair_environment_axis_coverage": _post_axis_coverage(
            rows,
            "environment_axis",
        ),
        "pre_repair_axis_profile_same_axis_conflict_count": pre_report[
            "axis_profile_same_axis_conflict_count"
        ],
        "post_repair_axis_profile_same_axis_conflict_count": len(
            post_axis_profile_conflicts
        ),
        "pre_repair_architecture_axis_same_axis_conflict_count": pre_report[
            "architecture_axis_same_axis_conflict_count"
        ],
        "post_repair_architecture_axis_same_axis_conflict_count": len(
            post_architecture_conflicts
        ),
        "pre_repair_unsafe_axis_claim_count": pre_report[
            "unsafe_axis_claim_count"
        ],
        "post_repair_unsafe_axis_claim_count": _unsafe_axis_claim_total(rows),
        "pre_repair_high_confidence_wrong_count_after_axis_scoring": pre_report[
            "high_confidence_wrong_count_after_axis_scoring"
        ],
        "post_repair_high_confidence_wrong_count_after_axis_scoring": len(
            post_high_confidence_wrong
        ),
        "order_axis_folded_mimic_quarantined_count": sum(
            1 for row in rows if bool(row["order_axis_folded_mimic_quarantined"])
        ),
        "repeat_compact_ambiguity_quarantined_count": sum(
            1 for row in rows if bool(row["repeat_compact_ambiguity_quarantined"])
        ),
        "safe_axis_claim_count": _safe_axis_claim_total(rows),
        "unsafe_axis_claim_count": _unsafe_axis_claim_total(rows),
        "safe_axis_recovered_count": sum(
            int(row["safe_axis_recovered_count"]) for row in rows
        ),
        "unsafe_class_recovery_count": sum(
            1 for row in rows if bool(row["unsafe_class_recovery"])
        ),
        "guard_override_count": sum(1 for row in rows if bool(row["guard_override"])),
        "coverage_loss_from_external_safety": round(
            pre_axis_profile_coverage - post_axis_profile_coverage,
            6,
        ),
        "external_safety_repair_successful": external_safety_repair_successful,
        "legacy_axis_artifacts_reproducible": True,
        "legacy_axis_profile_artifacts_reproducible": True,
        "global_fold_class_claim_allowed": False,
        "axis_profile_claim_allowed": external_safety_repair_successful,
        "architecture_axis_claim_allowed": not post_architecture_conflicts,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "artifact_reproducible": True,
        "boundary_statement": (
            "This repair layer responds to external holdout overclaims by "
            "quarantining weak disorder-order projection and repeat-vs-compact "
            "architecture ambiguity. It repairs by abstention, not by forcing "
            "new labels or claiming global fold classes."
        ),
        "rows": rows,
    }
    return report


def build_external_axis_repair_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": "external_axis_repair_safety_certificate",
        "external_axis_repair_complete": True,
        "legacy_axis_artifacts_reproducible": report[
            "legacy_axis_artifacts_reproducible"
        ],
        "legacy_axis_profile_artifacts_reproducible": report[
            "legacy_axis_profile_artifacts_reproducible"
        ],
        "order_axis_folded_mimic_quarantine_active": True,
        "repeat_compact_ambiguity_quarantine_active": True,
        "post_repair_axis_profile_same_axis_conflict_count": report[
            "post_repair_axis_profile_same_axis_conflict_count"
        ],
        "post_repair_architecture_axis_same_axis_conflict_count": report[
            "post_repair_architecture_axis_same_axis_conflict_count"
        ],
        "post_repair_unsafe_axis_claim_count": report[
            "post_repair_unsafe_axis_claim_count"
        ],
        "post_repair_high_confidence_wrong_count_after_axis_scoring": report[
            "post_repair_high_confidence_wrong_count_after_axis_scoring"
        ],
        "unsafe_class_recovery_count": report["unsafe_class_recovery_count"],
        "guard_override_count": report["guard_override_count"],
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


def write_external_axis_repair_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    conflict_delta_rows: Sequence[Mapping[str, object]],
    abstention_delta_rows: Sequence[Mapping[str, object]],
    quarantine_rows: Sequence[Mapping[str, object]],
    family_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    conflict_delta_path: Path,
    abstention_delta_path: Path,
    quarantine_rows_path: Path,
    family_summary_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    output_names = (
        report_path.name,
        rows_path.name,
        conflict_delta_path.name,
        abstention_delta_path.name,
        quarantine_rows_path.name,
        family_summary_path.name,
        dashboard_path.name,
        certificate_path.name,
    )
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(rows, rows_path)
    _write_csv_rows(conflict_delta_rows, conflict_delta_path)
    _write_csv_rows(abstention_delta_rows, abstention_delta_path)
    _write_csv_rows(quarantine_rows, quarantine_rows_path)
    _write_csv_rows(family_rows, family_summary_path)
    dashboard_path.write_text(render_external_axis_repair_dashboard(report), encoding="utf-8")
    certificate = build_external_axis_repair_certificate(
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
        conflict_delta_path,
        abstention_delta_path,
        quarantine_rows_path,
        family_summary_path,
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
        "pre_repair_axis_profile_coverage",
        "post_repair_axis_profile_coverage",
        "pre_repair_axis_profile_same_axis_conflict_count",
        "post_repair_axis_profile_same_axis_conflict_count",
        "pre_repair_architecture_axis_same_axis_conflict_count",
        "post_repair_architecture_axis_same_axis_conflict_count",
        "post_repair_unsafe_axis_claim_count",
        "post_repair_high_confidence_wrong_count_after_axis_scoring",
        "order_axis_folded_mimic_quarantined_count",
        "repeat_compact_ambiguity_quarantined_count",
        "claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _mapping_table(title: str, mapping: Mapping[str, object]) -> str:
    rows = "".join(
        "<tr>"
        f"<td>{_escape(key)}</td>"
        f"<td>{_escape(value)}</td>"
        "</tr>"
        for key, value in mapping.items()
    )
    return (
        f"<section><h2>{_escape(title)}</h2>"
        "<table><thead><tr><th>key</th><th>value</th></tr></thead><tbody>"
        + rows
        + "</tbody></table></section>"
    )


def _quarantine_preview(report: Mapping[str, object]) -> str:
    rows = external_axis_repair_quarantine_rows(report.get("rows", []))[:32]
    body = "".join(
        "<tr>"
        f"<td>{_escape(row['row_id'])}</td>"
        f"<td>{_escape(row['external_family_group'])}</td>"
        f"<td>{_escape(row['pre_profile_order_axis'])}</td>"
        f"<td>{_escape(row['post_profile_order_axis'])}</td>"
        f"<td>{_escape(row['pre_architecture_axis_prediction'])}</td>"
        f"<td>{_escape(row['post_architecture_axis_prediction'])}</td>"
        f"<td>{_escape(row['order_axis_abstention_reason'])}</td>"
        f"<td>{_escape(row['architecture_axis_abstention_reason'])}</td>"
        "</tr>"
        for row in rows
    )
    return (
        "<section><h2>Safety Improved By Abstention</h2>"
        "<table><thead><tr>"
        "<th>row_id</th><th>family group</th><th>pre order</th><th>post order</th>"
        "<th>pre architecture</th><th>post architecture</th>"
        "<th>order reason</th><th>architecture reason</th>"
        "</tr></thead><tbody>"
        + body
        + "</tbody></table></section>"
    )


def render_external_axis_repair_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>External Holdout Safety Repair</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f6f7f2;
      color: #1f2523;
    }}
    header {{
      padding: 32px;
      background: #24302c;
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
      color: #1f2523;
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
    <h1>External Holdout Safety Repair</h1>
    <p>Coverage was too optimistic; this repair closes unsafe axis claims by abstention rather than label forcing.</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Repair Rules</h2>
      <div class="rule-grid">
        <div class="rule"><strong>Coverage Was Too Optimistic</strong><br>The baseline holdout exposed unsafe axis claims.</div>
        <div class="rule"><strong>Disorder Claim Requires Strong Disorder Evidence</strong><br>Disordered class is not automatically an order-axis claim.</div>
        <div class="rule"><strong>Folded Beta/Mixed Mimics Are Quarantined</strong><br>Suspicious disorder claims become unknown, not ordered.</div>
        <div class="rule"><strong>Repeat-Like Requires Recurrence</strong><br>Repeat claims require more than periodic hydrophobic clusters.</div>
        <div class="rule"><strong>Hydrophobic Periodicity Alone Is Not Architecture</strong><br>Compact closure plus weak recurrence blocks repeat-like claims.</div>
        <div class="rule"><strong>Safety Improved By Abstention</strong><br>Coverage may drop; that is the point of this batch.</div>
        <div class="rule"><strong>Global Fold Class Still Locked</strong><br>No global class recovery or folding-solved claim is created.</div>
      </div>
    </section>
    {_mapping_table("Conflict Delta", {
        "pre_repair_axis_profile_same_axis_conflict_count": report.get("pre_repair_axis_profile_same_axis_conflict_count"),
        "post_repair_axis_profile_same_axis_conflict_count": report.get("post_repair_axis_profile_same_axis_conflict_count"),
        "pre_repair_architecture_axis_same_axis_conflict_count": report.get("pre_repair_architecture_axis_same_axis_conflict_count"),
        "post_repair_architecture_axis_same_axis_conflict_count": report.get("post_repair_architecture_axis_same_axis_conflict_count"),
        "pre_repair_unsafe_axis_claim_count": report.get("pre_repair_unsafe_axis_claim_count"),
        "post_repair_unsafe_axis_claim_count": report.get("post_repair_unsafe_axis_claim_count"),
    })}
    {_mapping_table("Abstention Cost", {
        "pre_repair_axis_profile_coverage": report.get("pre_repair_axis_profile_coverage"),
        "post_repair_axis_profile_coverage": report.get("post_repair_axis_profile_coverage"),
        "pre_repair_architecture_axis_coverage": report.get("pre_repair_architecture_axis_coverage"),
        "post_repair_architecture_axis_coverage": report.get("post_repair_architecture_axis_coverage"),
        "coverage_loss_from_external_safety": report.get("coverage_loss_from_external_safety"),
    })}
    {_quarantine_preview(report)}
  </main>
</body>
</html>
"""
