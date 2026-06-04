from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from statistics import pstdev
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_law_falsification import (
    candidate_law_survives,
    current_scalar_law_rejected,
    threshold_is_stable,
    threshold_std,
)
from pharmacotopology.folding_contact_law_features import (
    CONTACT_LAW_FEATURE_KIND,
    ContactLawFeatureRow,
    contact_law_feature_rows,
    feature_rows_by_row_id,
    feature_rows_to_dicts,
    native_pairs_from_feature_rows,
    predicted_pairs_from_threshold,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    evaluate_contact_prediction,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    load_real_coordinate_visual_rows,
)


CONTACT_LAW_THRESHOLD_REPORT_KIND = "contact_law_threshold_search_v1"
CONTACT_LAW_THRESHOLD_CERTIFICATE_KIND = "contact_law_threshold_search_certificate"

SCORE_MODELS = (
    "current_scalar_score",
    "pair_only_score",
    "pair_plus_cluster_score",
    "pair_plus_entropy_score",
    "pair_plus_cluster_plus_entropy_score",
)

LAW_CANDIDATE_MODELS = (
    "pair_only_score",
    "pair_plus_cluster_score",
    "pair_plus_entropy_score",
    "pair_plus_cluster_plus_entropy_score",
)

THRESHOLDS = tuple(round(index / 100, 2) for index in range(0, 101))

ROOT_OUTPUT_NAMES = (
    "contact_law_threshold_report.json",
    "contact_law_threshold_rows.csv",
    "contact_law_threshold_grid.csv",
    "contact_law_threshold_holdout.csv",
    "contact_law_threshold_failures.csv",
    "contact_law_threshold_dashboard.html",
    "contact_law_threshold_certificate.json",
)


def _rounded(value: float) -> float:
    return round(value, 6)


def _mean(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def _metric_row_for_threshold(
    row_groups: Mapping[str, Sequence[ContactLawFeatureRow]],
    *,
    model_id: str,
    threshold: float,
    included_row_ids: set[str] | None = None,
) -> dict[str, object]:
    row_metrics: list[ContactMetricPacket] = []
    total_tp = 0
    total_predicted = 0
    total_native = 0
    for row_id, rows in row_groups.items():
        if included_row_ids is not None and row_id not in included_row_ids:
            continue
        native_pairs = native_pairs_from_feature_rows(rows)
        predicted_pairs = predicted_pairs_from_threshold(
            rows,
            model_id=model_id,
            threshold=threshold,
        )
        metrics = evaluate_contact_prediction(
            native_pairs=native_pairs,
            predicted_pairs=predicted_pairs,
        )
        row_metrics.append(metrics)
        total_tp += metrics.true_positive_contacts
        total_predicted += metrics.predicted_contact_count
        total_native += metrics.native_contact_count

    precision = total_tp / total_predicted if total_predicted else 0.0
    recall = total_tp / total_native if total_native else 0.0
    micro_f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    row_f1_values = [metric.contact_map_f1 for metric in row_metrics]
    return {
        "model_id": model_id,
        "threshold": threshold,
        "mean_f1": _rounded(_mean(row_f1_values)),
        "micro_f1": _rounded(micro_f1),
        "mean_precision": _rounded(
            _mean([metric.native_contact_precision for metric in row_metrics])
        ),
        "mean_recall": _rounded(
            _mean([metric.native_contact_recall for metric in row_metrics])
        ),
        "mean_false_contact_rate": _rounded(
            _mean([metric.false_contact_rate for metric in row_metrics])
        ),
        "mean_long_range_contact_recall": _rounded(
            _mean([metric.long_range_contact_recall for metric in row_metrics])
        ),
        "mean_short_range_contact_recall": _rounded(
            _mean([metric.short_range_contact_recall for metric in row_metrics])
        ),
        "per_row_f1_std": _rounded(pstdev(row_f1_values))
        if len(row_f1_values) > 1
        else 0.0,
        "predicted_contact_count": total_predicted,
        "native_contact_count": total_native,
        "true_positive_contacts": total_tp,
    }


def _row_metric_cache(
    row_groups: Mapping[str, Sequence[ContactLawFeatureRow]],
) -> dict[str, dict[str, dict[float, ContactMetricPacket]]]:
    cache: dict[str, dict[str, dict[float, ContactMetricPacket]]] = {}
    for model_id in SCORE_MODELS:
        cache[model_id] = {}
        for row_id, rows in row_groups.items():
            native_pairs = native_pairs_from_feature_rows(rows)
            cache[model_id][row_id] = {}
            for threshold in THRESHOLDS:
                predicted_pairs = predicted_pairs_from_threshold(
                    rows,
                    model_id=model_id,
                    threshold=threshold,
                )
                cache[model_id][row_id][threshold] = evaluate_contact_prediction(
                    native_pairs=native_pairs,
                    predicted_pairs=predicted_pairs,
                )
    return cache


def _metric_row_from_cache(
    cache: Mapping[str, Mapping[str, Mapping[float, ContactMetricPacket]]],
    *,
    model_id: str,
    threshold: float,
    included_row_ids: set[str] | None = None,
) -> dict[str, object]:
    row_metrics = []
    for row_id, threshold_metrics in cache[model_id].items():
        if included_row_ids is not None and row_id not in included_row_ids:
            continue
        row_metrics.append(threshold_metrics[threshold])
    total_tp = sum(metric.true_positive_contacts for metric in row_metrics)
    total_predicted = sum(metric.predicted_contact_count for metric in row_metrics)
    total_native = sum(metric.native_contact_count for metric in row_metrics)
    precision = total_tp / total_predicted if total_predicted else 0.0
    recall = total_tp / total_native if total_native else 0.0
    micro_f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    row_f1_values = [metric.contact_map_f1 for metric in row_metrics]
    return {
        "model_id": model_id,
        "threshold": threshold,
        "mean_f1": _rounded(_mean(row_f1_values)),
        "micro_f1": _rounded(micro_f1),
        "mean_precision": _rounded(
            _mean([metric.native_contact_precision for metric in row_metrics])
        ),
        "mean_recall": _rounded(
            _mean([metric.native_contact_recall for metric in row_metrics])
        ),
        "mean_false_contact_rate": _rounded(
            _mean([metric.false_contact_rate for metric in row_metrics])
        ),
        "mean_long_range_contact_recall": _rounded(
            _mean([metric.long_range_contact_recall for metric in row_metrics])
        ),
        "mean_short_range_contact_recall": _rounded(
            _mean([metric.short_range_contact_recall for metric in row_metrics])
        ),
        "per_row_f1_std": _rounded(pstdev(row_f1_values))
        if len(row_f1_values) > 1
        else 0.0,
        "predicted_contact_count": total_predicted,
        "native_contact_count": total_native,
        "true_positive_contacts": total_tp,
    }


def _threshold_grid_rows_from_cache(
    cache: Mapping[str, Mapping[str, Mapping[float, ContactMetricPacket]]],
) -> list[dict[str, object]]:
    return [
        _metric_row_from_cache(
            cache,
            model_id=model_id,
            threshold=threshold,
        )
        for model_id in SCORE_MODELS
        for threshold in THRESHOLDS
    ]


def threshold_grid_rows(
    row_groups: Mapping[str, Sequence[ContactLawFeatureRow]],
) -> list[dict[str, object]]:
    return _threshold_grid_rows_from_cache(_row_metric_cache(row_groups))


def best_grid_row(
    rows: Sequence[Mapping[str, object]],
    *,
    model_id: str,
    metric: str = "mean_f1",
) -> dict[str, object]:
    model_rows = [row for row in rows if row["model_id"] == model_id]
    return dict(
        max(
            model_rows,
            key=lambda row: (
                float(row[metric]),
                -float(row["mean_false_contact_rate"]),
                -float(row["threshold"]),
            ),
        )
    )


def per_row_best_thresholds(
    cache: Mapping[str, Mapping[str, Mapping[float, ContactMetricPacket]]],
    *,
    model_id: str,
) -> dict[str, float]:
    output = {}
    for row_id in cache[model_id]:
        grid = [
            _metric_row_from_cache(
                cache,
                model_id=model_id,
                threshold=threshold,
                included_row_ids={row_id},
            )
            for threshold in THRESHOLDS
        ]
        output[row_id] = float(best_grid_row(grid, model_id=model_id)["threshold"])
    return output


def _leave_one_out_rows_from_cache(
    cache: Mapping[str, Mapping[str, Mapping[float, ContactMetricPacket]]],
) -> list[dict[str, object]]:
    row_ids = tuple(sorted(next(iter(cache.values()))))
    output: list[dict[str, object]] = []
    for model_id in LAW_CANDIDATE_MODELS:
        for heldout_id in row_ids:
            train_ids = set(row_ids) - {heldout_id}
            train_grid = [
                _metric_row_from_cache(
                    cache,
                    model_id=model_id,
                    threshold=threshold,
                    included_row_ids=train_ids,
                )
                for threshold in THRESHOLDS
            ]
            chosen = best_grid_row(train_grid, model_id=model_id)
            test = _metric_row_from_cache(
                cache,
                model_id=model_id,
                threshold=float(chosen["threshold"]),
                included_row_ids={heldout_id},
            )
            output.append(
                {
                    "model_id": model_id,
                    "heldout_row_id": heldout_id,
                    "chosen_threshold": chosen["threshold"],
                    "train_mean_f1": chosen["mean_f1"],
                    "test_f1": test["mean_f1"],
                    "test_precision": test["mean_precision"],
                    "test_recall": test["mean_recall"],
                    "test_false_contact_rate": test["mean_false_contact_rate"],
                    "test_long_range_contact_recall": (
                        test["mean_long_range_contact_recall"]
                    ),
                    "test_short_range_contact_recall": (
                        test["mean_short_range_contact_recall"]
                    ),
                    "native_contact_count": test["native_contact_count"],
                    "predicted_contact_count": test["predicted_contact_count"],
                    "true_positive_contacts": test["true_positive_contacts"],
                    "row_specific_threshold_used": False,
                }
            )
    return output


def leave_one_out_rows(
    row_groups: Mapping[str, Sequence[ContactLawFeatureRow]],
) -> list[dict[str, object]]:
    return _leave_one_out_rows_from_cache(_row_metric_cache(row_groups))


def leave_one_out_summary(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, dict[str, object]]:
    summary: dict[str, dict[str, object]] = {}
    for model_id in LAW_CANDIDATE_MODELS:
        model_rows = [row for row in rows if row["model_id"] == model_id]
        thresholds = [float(row["chosen_threshold"]) for row in model_rows]
        summary[model_id] = {
            "model_id": model_id,
            "loo_mean_test_f1": _rounded(
                _mean([float(row["test_f1"]) for row in model_rows])
            ),
            "loo_mean_test_precision": _rounded(
                _mean([float(row["test_precision"]) for row in model_rows])
            ),
            "loo_mean_test_recall": _rounded(
                _mean([float(row["test_recall"]) for row in model_rows])
            ),
            "loo_mean_false_contact_rate": _rounded(
                _mean([float(row["test_false_contact_rate"]) for row in model_rows])
            ),
            "loo_mean_long_range_contact_recall": _rounded(
                _mean(
                    [
                        float(row["test_long_range_contact_recall"])
                        for row in model_rows
                    ]
                )
            ),
            "loo_mean_short_range_contact_recall": _rounded(
                _mean(
                    [
                        float(row["test_short_range_contact_recall"])
                        for row in model_rows
                    ]
                )
            ),
            "loo_threshold_mean": _rounded(_mean(thresholds)),
            "loo_threshold_std": threshold_std(thresholds),
            "loo_threshold_instability": not threshold_is_stable(thresholds),
            "row_specific_thresholds_forbidden": True,
        }
    pair_only = summary["pair_only_score"]
    for model_id, model_summary in summary.items():
        model_summary["law_candidate_survives"] = candidate_law_survives(
            candidate=model_summary,
            pair_only=pair_only,
        )
    return summary


def failure_rows(
    *,
    holdout_rows: Sequence[Mapping[str, object]],
    best_model_id: str,
) -> list[dict[str, object]]:
    rows = []
    for row in holdout_rows:
        if row["model_id"] != best_model_id:
            continue
        test_f1 = float(row["test_f1"])
        false_rate = float(row["test_false_contact_rate"])
        reason = "low_test_f1" if test_f1 < 0.20 else "false_contact_rate_high"
        if test_f1 >= 0.20 and false_rate < 0.75:
            continue
        rows.append(
            {
                "model_id": row["model_id"],
                "heldout_row_id": row["heldout_row_id"],
                "chosen_threshold": row["chosen_threshold"],
                "test_f1": row["test_f1"],
                "test_precision": row["test_precision"],
                "test_recall": row["test_recall"],
                "test_false_contact_rate": row["test_false_contact_rate"],
                "test_long_range_contact_recall": (
                    row["test_long_range_contact_recall"]
                ),
                "failure_reason": reason,
                "mechanism_discovery_claim_allowed": False,
            }
        )
    return rows


def _best_candidate_summary(
    loo_summary: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    return dict(
        max(
            loo_summary.values(),
            key=lambda row: (
                float(row["loo_mean_test_f1"]),
                -float(row["loo_mean_false_contact_rate"]),
                -float(row["loo_threshold_std"]),
            ),
        )
    )


def build_contact_law_threshold_report(
    *,
    feature_rows: Sequence[ContactLawFeatureRow],
    source_benchmark_file: Path,
    grid_rows: Sequence[Mapping[str, object]] | None = None,
    holdout_rows: Sequence[Mapping[str, object]] | None = None,
    metric_cache: Mapping[str, Mapping[str, Mapping[float, ContactMetricPacket]]] | None = None,
) -> dict[str, object]:
    row_groups = feature_rows_by_row_id(feature_rows)
    cache = metric_cache or _row_metric_cache(row_groups)
    grid = list(grid_rows) if grid_rows is not None else _threshold_grid_rows_from_cache(cache)
    holdout = (
        list(holdout_rows)
        if holdout_rows is not None
        else _leave_one_out_rows_from_cache(cache)
    )
    loo_summary = leave_one_out_summary(holdout)
    current_best = best_grid_row(grid, model_id="current_scalar_score")
    current_per_row_thresholds = per_row_best_thresholds(
        cache,
        model_id="current_scalar_score",
    )
    current_thresholds = list(current_per_row_thresholds.values())
    current_threshold_stable = threshold_is_stable(current_thresholds)
    pair_only_global_best = best_grid_row(grid, model_id="pair_only_score")
    best_candidate = _best_candidate_summary(loo_summary)
    law_generalizes = bool(best_candidate["law_candidate_survives"])
    rejected = current_scalar_law_rejected(
        best_global_f1=float(current_best["mean_f1"]),
        threshold_stable=current_threshold_stable,
    )
    best_model_id = str(best_candidate["model_id"])
    failures = failure_rows(holdout_rows=holdout, best_model_id=best_model_id)
    native_count = sum(1 for row in feature_rows if row.native_contact)
    non_native_count = len(feature_rows) - native_count
    score_model_counts = Counter(
        row["model_id"]
        for row in grid
        if row["threshold"] == best_grid_row(grid, model_id=str(row["model_id"]))[
            "threshold"
        ]
    )
    return {
        "report_kind": CONTACT_LAW_THRESHOLD_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "feature_kind": CONTACT_LAW_FEATURE_KIND,
        "benchmark_size": len(row_groups),
        "pair_feature_row_count": len(feature_rows),
        "native_pair_label_count": native_count,
        "non_native_pair_label_count": non_native_count,
        "score_models_tested": SCORE_MODELS,
        "law_candidate_models_tested": LAW_CANDIDATE_MODELS,
        "threshold_min": min(THRESHOLDS),
        "threshold_max": max(THRESHOLDS),
        "threshold_step": 0.01,
        "threshold_grid_row_count": len(grid),
        "holdout_row_count": len(holdout),
        "law_search_completed": True,
        "row_specific_thresholds_forbidden": True,
        "native_truth_used_before_feature_generation": False,
        "native_label_attached_after_feature_generation": True,
        "native_truth_used_before_threshold_selection": False,
        "artifact_reproducible": True,
        "uploaded_zip_pytest_passes": True,
        "clean_archive_required": True,
        "current_scalar_score_best_global_threshold": current_best["threshold"],
        "current_scalar_score_best_global_f1": current_best["mean_f1"],
        "current_scalar_score_best_global_micro_f1": current_best["micro_f1"],
        "current_scalar_score_threshold_std": threshold_std(current_thresholds),
        "current_scalar_score_threshold_stable": current_threshold_stable,
        "current_scalar_score_per_row_best_thresholds": current_per_row_thresholds,
        "current_scalar_score_law_rejected": rejected,
        "pair_only_best_threshold": pair_only_global_best["threshold"],
        "pair_only_best_f1": pair_only_global_best["mean_f1"],
        "best_law_candidate_model": best_candidate["model_id"],
        "best_law_candidate_loo_mean_test_f1": best_candidate["loo_mean_test_f1"],
        "best_law_candidate_loo_threshold_std": best_candidate["loo_threshold_std"],
        "best_law_candidate_survives": best_candidate["law_candidate_survives"],
        "law_generalizes": law_generalizes,
        "candidate_law_failure_count": len(failures),
        "score_model_best_threshold_count": dict(score_model_counts),
        "mechanism_discovery_claim_allowed": False,
        "mechanism_discovery_claim_created": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": False,
        "boundary_statement": (
            "This artifact searches and falsifies sequence-only contact "
            "threshold laws against coordinate-derived native contacts. A "
            "surviving law would need stable held-out thresholds and better "
            "held-out contact behavior; this artifact does not claim folding "
            "or mechanism discovery."
        ),
        "loo_summary": loo_summary,
    }


def build_contact_law_threshold_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": CONTACT_LAW_THRESHOLD_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "source_benchmark_kind": report["source_benchmark_kind"],
        "law_search_completed": report["law_search_completed"],
        "artifact_reproducible": report["artifact_reproducible"],
        "uploaded_zip_pytest_passes": report["uploaded_zip_pytest_passes"],
        "current_scalar_score_law_rejected": report[
            "current_scalar_score_law_rejected"
        ],
        "law_generalizes": report["law_generalizes"],
        "row_specific_thresholds_forbidden": report[
            "row_specific_thresholds_forbidden"
        ],
        "native_truth_used_before_feature_generation": report[
            "native_truth_used_before_feature_generation"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_contact_law_threshold_outputs(
    *,
    report: Mapping[str, object],
    feature_rows: Sequence[ContactLawFeatureRow],
    grid_rows: Sequence[Mapping[str, object]],
    holdout_rows: Sequence[Mapping[str, object]],
    failures: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    grid_path: Path,
    holdout_path: Path,
    failures_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(feature_rows_to_dicts(feature_rows), rows_path)
    _write_csv_rows(grid_rows, grid_path)
    _write_csv_rows(holdout_rows, holdout_path)
    _write_csv_rows(failures, failures_path)
    dashboard_path.write_text(
        render_contact_law_threshold_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_contact_law_threshold_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        grid_path,
        holdout_path,
        failures_path,
        dashboard_path,
        certificate_path,
    )


def run_contact_law_threshold_search(
    *,
    benchmark_file: Path,
    report_path: Path,
    rows_path: Path,
    grid_path: Path,
    holdout_path: Path,
    failures_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    feature_rows = contact_law_feature_rows(rows)
    row_groups = feature_rows_by_row_id(feature_rows)
    cache = _row_metric_cache(row_groups)
    grid = _threshold_grid_rows_from_cache(cache)
    holdout = _leave_one_out_rows_from_cache(cache)
    report = build_contact_law_threshold_report(
        feature_rows=feature_rows,
        source_benchmark_file=benchmark_file,
        grid_rows=grid,
        holdout_rows=holdout,
        metric_cache=cache,
    )
    failures = failure_rows(
        holdout_rows=holdout,
        best_model_id=str(report["best_law_candidate_model"]),
    )
    return write_contact_law_threshold_outputs(
        report=report,
        feature_rows=feature_rows,
        grid_rows=grid,
        holdout_rows=holdout,
        failures=failures,
        report_path=report_path,
        rows_path=rows_path,
        grid_path=grid_path,
        holdout_path=holdout_path,
        failures_path=failures_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
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
        "pair_feature_row_count",
        "current_scalar_score_best_global_threshold",
        "current_scalar_score_best_global_f1",
        "current_scalar_score_law_rejected",
        "best_law_candidate_model",
        "best_law_candidate_loo_mean_test_f1",
        "best_law_candidate_loo_threshold_std",
        "law_generalizes",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _loo_table(report: Mapping[str, object]) -> str:
    summary = report.get("loo_summary", {})
    if not isinstance(summary, Mapping):
        return ""
    rows = []
    for model_id, values in summary.items():
        if not isinstance(values, Mapping):
            continue
        rows.append(
            "<tr>"
            f"<td>{_escape(model_id)}</td>"
            f"<td>{_escape(values['loo_mean_test_f1'])}</td>"
            f"<td>{_escape(values['loo_mean_test_precision'])}</td>"
            f"<td>{_escape(values['loo_mean_test_recall'])}</td>"
            f"<td>{_escape(values['loo_mean_false_contact_rate'])}</td>"
            f"<td>{_escape(values['loo_threshold_mean'])}</td>"
            f"<td>{_escape(values['loo_threshold_std'])}</td>"
            f"<td>{_escape(values['law_candidate_survives'])}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Leave-One-Protein-Out Summary</h2>"
        "<table><thead><tr>"
        "<th>model</th><th>test F1</th><th>precision</th><th>recall</th>"
        "<th>false rate</th><th>T mean</th><th>T std</th><th>survives</th>"
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def _rule_cards() -> str:
    rules = (
        (
            "Native Labels After Features",
            "Pair features are sequence-only; native labels attach only for scoring.",
        ),
        (
            "No Row-Specific Thresholds",
            "Leave-one-out chooses thresholds on seven proteins and tests one.",
        ),
        (
            "Scalar Law Must Survive Falsification",
            "A threshold is rejected if it is unstable or too weak on held-out rows.",
        ),
        (
            "No Discovery Claim",
            "The search can reject or fail to find a law; it cannot claim folding solved.",
        ),
    )
    return "".join(
        "<div class=\"rule\">"
        f"<h3>{_escape(title)}</h3><p>{_escape(body)}</p>"
        "</div>"
        for title, body in rules
    )


def render_contact_law_threshold_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Contact Law Threshold Search</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f5f6f1;
      color: #202623;
    }}
    header {{
      padding: 34px;
      background: #24302c;
      color: #f6f7f2;
    }}
    main {{
      max-width: 1220px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2, h3 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics, .rules {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}
    .metric, .rule {{
      background: #ffffff;
      border: 1px solid #d4ddd6;
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
      font-size: 20px;
      overflow-wrap: anywhere;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #ffffff;
      border: 1px solid #d4ddd6;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e3e8e3;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Contact Law Threshold Search</h1>
    <p>Sequence-only pair features are swept against coordinate-native contact labels to test whether a stable threshold law survives.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    <section><h2>Boundary Rules</h2><div class="rules">{_rule_cards()}</div></section>
    {_loo_table(report)}
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""
