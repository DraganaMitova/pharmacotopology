from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    build_hierarchical_gate_report,
    gate_failure_rows,
    gate_path_rows,
    hierarchical_gate_rows,
    load_hierarchical_gate_inputs,
    write_hierarchical_gate_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_10_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_rows.csv"
)
DEFAULT_GATE_PATHS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_gate_paths.csv"
)
DEFAULT_GATE_FAILURES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_gate_failures.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_dashboard.html"
)
DEFAULT_BASELINE_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_report.json"
)


def _load_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return dict(data) if isinstance(data, Mapping) else {}


def _class_distribution(rows: Sequence[Mapping[str, object]]) -> dict[str, int]:
    return dict(sorted(Counter(str(row["structure_fold_class"]) for row in rows).items()))


def _accuracy_by_class(rows: Sequence[Mapping[str, object]]) -> dict[str, float]:
    totals = Counter(str(row["structure_fold_class"]) for row in rows)
    matches = Counter(
        str(row["structure_fold_class"])
        for row in rows
        if bool(row["prediction_structure_class_match"])
    )
    return {
        fold_class: round(matches[fold_class] / total, 6)
        for fold_class, total in sorted(totals.items())
    }


def _abstention_rate_by_class(rows: Sequence[Mapping[str, object]]) -> dict[str, float]:
    totals = Counter(str(row["structure_fold_class"]) for row in rows)
    abstentions = Counter(
        str(row["structure_fold_class"])
        for row in rows
        if bool(row["abstained"])
    )
    return {
        fold_class: round(abstentions[fold_class] / total, 6)
        for fold_class, total in sorted(totals.items())
    }


def _gate_path_distribution_by_class(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, dict[str, int]]:
    distribution: dict[str, Counter[str]] = {}
    for row in rows:
        fold_class = str(row["structure_fold_class"])
        path = str(row["gate_path"]).split(" | ")[0]
        distribution.setdefault(fold_class, Counter())[path] += 1
    return {
        fold_class: dict(sorted(paths.items()))
        for fold_class, paths in sorted(distribution.items())
    }


def _confusion_matrix(rows: Sequence[Mapping[str, object]]) -> dict[str, dict[str, int]]:
    actual_classes = sorted({str(row["structure_fold_class"]) for row in rows})
    predicted_classes = sorted({str(row["predicted_fold_class"]) for row in rows})
    matrix: dict[str, dict[str, int]] = {
        actual: {predicted: 0 for predicted in predicted_classes}
        for actual in actual_classes
    }
    for row in rows:
        matrix[str(row["structure_fold_class"])][str(row["predicted_fold_class"])] += 1
    return matrix


def _per_class_stability_status(rows: Sequence[Mapping[str, object]]) -> dict[str, str]:
    accuracy = _accuracy_by_class(rows)
    abstention = _abstention_rate_by_class(rows)
    status = {}
    for fold_class in sorted(accuracy):
        if accuracy[fold_class] >= 0.50:
            status[fold_class] = "stable_enough_for_review"
        elif abstention.get(fold_class, 0.0) >= 0.50:
            status[fold_class] = "mostly_abstained"
        else:
            status[fold_class] = "unstable_failure_review_required"
    return status


def _stability_status(report: Mapping[str, object]) -> str:
    if (
        float(report.get("prediction_vs_structure_accuracy", 0.0)) >= 0.50
        and int(report.get("high_confidence_wrong_count", 0)) <= 2
        and int(report.get("false_beta_from_disorder_count", 0)) == 0
        and int(report.get("flexible_segmentation_false_multidomain_count", 0)) <= 2
    ):
        return "stable_on_50_row_stress_test"
    if float(report.get("prediction_vs_structure_accuracy", 0.0)) < 0.35:
        return "unstable_accuracy_drop"
    return "mixed_stability_review_required"


def _failure_mode_delta(current: Mapping[str, object], baseline: Mapping[str, object]) -> dict[str, object]:
    keys = (
        "false_beta_from_disorder_count",
        "false_mixed_from_alpha_count",
        "flexible_segmentation_false_multidomain_count",
        "high_confidence_wrong_count",
    )
    result: dict[str, object] = {}
    for key in keys:
        current_value = int(current.get(key, 0))
        baseline_value = int(baseline.get(key, 0))
        if current_value == baseline_value == 0:
            status = "preserved_zero"
        elif current_value <= baseline_value:
            status = "not_worse"
        else:
            status = "new_or_increased_failure_mode"
        result[key] = {
            "baseline": baseline_value,
            "current": current_value,
            "delta": current_value - baseline_value,
            "status": status,
        }
    return result


def augment_report_for_stability(
    report: dict[str, object],
    *,
    rows: Sequence[Mapping[str, object]],
    benchmark_file: Path,
    baseline_report_file: Path,
) -> dict[str, object]:
    benchmark_metadata = _load_json(benchmark_file)
    baseline = _load_json(baseline_report_file)
    external_rows = int(benchmark_metadata.get("external_rows", len(rows)))
    report["external_rows"] = external_rows
    report["class_distribution"] = benchmark_metadata.get(
        "class_distribution",
        _class_distribution(rows),
    )
    report["structure_class_distribution"] = _class_distribution(rows)
    report["accuracy_by_class"] = _accuracy_by_class(rows)
    report["abstention_rate_by_class"] = _abstention_rate_by_class(rows)
    report["gate_path_distribution_by_true_class"] = (
        _gate_path_distribution_by_class(rows)
    )
    report["confusion_matrix"] = _confusion_matrix(rows)
    report["high_confidence_wrong_cases"] = [
        {
            "protein_id": row["protein_id"],
            "predicted_fold_class": row["predicted_fold_class"],
            "structure_fold_class": row["structure_fold_class"],
            "confidence": row["confidence"],
            "gate_path": row["gate_path"],
        }
        for row in rows
        if bool(row["high_confidence_wrong"])
    ]
    report["per_class_stability_status"] = _per_class_stability_status(rows)
    report["locked_after_generation"] = benchmark_metadata.get(
        "locked_after_generation",
        False,
    )
    report["no_retuning_flag"] = benchmark_metadata.get("no_retuning_flag", False)
    report["lock_certificate"] = benchmark_metadata.get("lock_certificate", {})
    report["stability_test_of_commit"] = benchmark_metadata.get(
        "stability_test_of_commit",
        "",
    )
    report["result_compared_to_10_row_benchmark"] = {
        "baseline_report": str(baseline_report_file),
        "baseline_benchmark_size": baseline.get("benchmark_size", 0),
        "current_benchmark_size": report.get("benchmark_size", 0),
        "recipe_frozen_before_50_row_run": True,
    }
    report["accuracy_delta_from_10"] = round(
        float(report.get("prediction_vs_structure_accuracy", 0.0))
        - float(baseline.get("prediction_vs_structure_accuracy", 0.0)),
        6,
    )
    report["abstention_delta_from_10"] = int(report.get("abstained_prediction_count", 0)) - int(
        baseline.get("abstained_prediction_count", 0)
    )
    report["high_confidence_wrong_delta_from_10"] = int(
        report.get("high_confidence_wrong_count", 0)
    ) - int(baseline.get("high_confidence_wrong_count", 0))
    report["failure_modes_preserved_or_new"] = _failure_mode_delta(report, baseline)
    report["stability_status"] = _stability_status(report)
    return report


def write_hierarchical_gate_certificate(
    report: Mapping[str, object],
    output_path: Path,
) -> Path:
    payload = {
        "benchmark_kind": report.get("benchmark_kind"),
        "benchmark_size": report.get("benchmark_size"),
        "external_rows": report.get("external_rows"),
        "class_distribution": report.get("class_distribution"),
        "locked_after_generation": report.get("locked_after_generation"),
        "no_retuning_flag": report.get("no_retuning_flag"),
        "stability_status": report.get("stability_status"),
        "prediction_vs_structure_accuracy": report.get(
            "prediction_vs_structure_accuracy"
        ),
        "high_confidence_wrong_count": report.get("high_confidence_wrong_count"),
        "claim_allowed": report.get("claim_allowed"),
        "folding_problem_solved": report.get("folding_problem_solved"),
        "lock_certificate": report.get("lock_certificate", {}),
        "result_compared_to_10_row_benchmark": report.get(
            "result_compared_to_10_row_benchmark",
            {},
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def write_hierarchical_gate_confusion_matrix(
    rows: Sequence[Mapping[str, object]],
    output_path: Path,
) -> Path:
    actual_classes = sorted({str(row["structure_fold_class"]) for row in rows})
    predicted_classes = sorted({str(row["predicted_fold_class"]) for row in rows})
    matrix: dict[str, dict[str, int]] = {
        actual: {predicted: 0 for predicted in predicted_classes}
        for actual in actual_classes
    }
    for row in rows:
        matrix[str(row["structure_fold_class"])][str(row["predicted_fold_class"])] += 1
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=["structure_fold_class"] + predicted_classes,
            lineterminator="\n",
        )
        writer.writeheader()
        for actual in actual_classes:
            output = {"structure_fold_class": actual}
            output.update(matrix[actual])
            writer.writerow(output)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the hierarchical folding decision gates against the locked "
            "10-row folding benchmark and emit gate path diagnostics."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--gate-paths-output", default=str(DEFAULT_GATE_PATHS_PATH))
    parser.add_argument(
        "--gate-failures-output",
        default=str(DEFAULT_GATE_FAILURES_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument(
        "--baseline-report",
        default=str(DEFAULT_BASELINE_REPORT_PATH),
        help="10-row hierarchical gate report used for stability deltas.",
    )
    parser.add_argument(
        "--certificate-output",
        default="",
        help="Optional stability certificate JSON path.",
    )
    parser.add_argument(
        "--confusion-output",
        default="",
        help="Optional structure-vs-prediction confusion matrix CSV path.",
    )
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_hierarchical_gate_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = hierarchical_gate_rows(references, evidence_rows)
    paths = gate_path_rows(rows)
    failures = gate_failure_rows(rows)
    report = build_hierarchical_gate_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    report = augment_report_for_stability(
        report,
        rows=rows,
        benchmark_file=benchmark_file,
        baseline_report_file=Path(args.baseline_report),
    )
    outputs = write_hierarchical_gate_outputs(
        report=report,
        rows=rows,
        path_rows=paths,
        failure_rows=failures,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        gate_paths_path=Path(args.gate_paths_output),
        gate_failures_path=Path(args.gate_failures_output),
        dashboard_path=Path(args.dashboard_output),
    )
    extra_outputs = []
    if args.certificate_output:
        extra_outputs.append(
            write_hierarchical_gate_certificate(
                report,
                Path(args.certificate_output),
            )
        )
    if args.confusion_output:
        extra_outputs.append(
            write_hierarchical_gate_confusion_matrix(
                rows,
                Path(args.confusion_output),
            )
        )
    for output in outputs:
        print(output)
    for output in extra_outputs:
        print(output)


if __name__ == "__main__":
    main()
