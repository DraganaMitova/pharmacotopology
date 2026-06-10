from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_coupling_negative_controls import (  # noqa: E402
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    import_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    _build_multiscale_physical_contexts,
    _run_multiscale_phase_aligned_external_novelty_boundary_selector,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_holdout_1ubq.locked.json")
DEFAULT_EXTERNAL_COUPLING_FILE = Path(
    "data/folding_real_coordinate_holdout_1ubq_external_couplings.v0.locked.json"
)
DEFAULT_REPORT = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "external_phase_holdout_1ubq_smoke_report.json"
)
DEFAULT_SELECTED_EVENTS = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "external_phase_holdout_1ubq_selected_events.csv"
)


def _metric_payload(run) -> dict[str, object]:
    metric = run.metric
    false_event_count = round(metric.false_nucleus_rate * metric.selected_event_count)
    return {
        "selector_name": run.selector_name,
        "control_kind": run.control_kind,
        "constraint_count": run.constraint_count,
        "selected_event_count": metric.selected_event_count,
        "false_event_count": false_event_count,
        "false_nucleus_rate": metric.false_nucleus_rate,
        "cluster_precision": metric.contact_cluster_precision,
        "long_range_contact_recall": metric.long_range_contact_recall,
        "coupling_constraint_recall": metric.coupling_constraint_recall,
        "real_vs_decoy_coupling_enrichment_ratio": (
            metric.real_vs_decoy_coupling_enrichment_ratio
        ),
    }


def run_external_phase_holdout_smoke(
    *,
    benchmark_file: Path,
    external_coupling_file: Path,
    report_output: Path,
    selected_events_output: Path,
    control_name: str,
) -> tuple[Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    import_result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_coupling_file,
    )
    physical_contexts = _build_multiscale_physical_contexts(rows)
    real_run = _run_multiscale_phase_aligned_external_novelty_boundary_selector(
        rows=rows,
        dataset=import_result.dataset,
        selector_name="external_multiscale_phase_aligned_external_novelty_boundary",
        control_kind="external_real_holdout_smoke",
        physical_contexts=physical_contexts,
    )
    controls = generate_external_coupling_negative_controls(
        rows=rows,
        dataset=import_result.dataset,
    )
    control_dataset = controls[control_name]
    control_run = _run_multiscale_phase_aligned_external_novelty_boundary_selector(
        rows=rows,
        dataset=control_dataset.dataset,
        selector_name=control_name,
        control_kind=control_dataset.control_kind,
        physical_contexts=physical_contexts,
    )
    real_metric = _metric_payload(real_run)
    control_metric = _metric_payload(control_run)
    report = {
        "report_kind": (
            "external_phase_aligned_external_novelty_single_holdout_smoke_v0"
        ),
        "benchmark_file": str(benchmark_file),
        "external_coupling_file": str(external_coupling_file),
        "row_count": len(rows),
        "row_statuses": [status.to_dict() for status in import_result.row_statuses],
        "accepted_external_constraint_count": len(import_result.dataset.constraints),
        "coordinate_truth_used_to_build_constraints": (
            import_result.dataset.coordinate_truth_tainted
        ),
        "native_truth_used_before_coupling_selection": (
            import_result.dataset.native_truth_tainted
        ),
        "oracle_constraint_control": import_result.dataset.oracle_constraint_control,
        "external_evolutionary_couplings_used": (
            import_result.dataset.external_evolutionary_couplings_used
        ),
        "selector": real_metric,
        "matched_control": control_metric,
        "beats_random_long_range_same_count_control": (
            real_metric["false_nucleus_rate"] <= control_metric["false_nucleus_rate"]
            and real_metric["long_range_contact_recall"]
            >= control_metric["long_range_contact_recall"]
            and real_metric["cluster_precision"] > control_metric["cluster_precision"]
        ),
        "folding_problem_solved": False,
        "claim_allowed": False,
    }
    report_output.parent.mkdir(parents=True, exist_ok=True)
    report_output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv_rows(real_run.selected_rows, selected_events_output)
    return report_output, selected_events_output


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the current phase-aligned external novelty selector on one "
            "holdout protein plus one matched control."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--external-coupling-file",
        default=str(DEFAULT_EXTERNAL_COUPLING_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT))
    parser.add_argument("--selected-events-output", default=str(DEFAULT_SELECTED_EVENTS))
    parser.add_argument(
        "--control-name",
        default="external_random_long_range_same_count",
    )
    args = parser.parse_args()

    for output in run_external_phase_holdout_smoke(
        benchmark_file=Path(args.benchmark_file),
        external_coupling_file=Path(args.external_coupling_file),
        report_output=Path(args.report_output),
        selected_events_output=Path(args.selected_events_output),
        control_name=args.control_name,
    ):
        print(output)


if __name__ == "__main__":
    main()
