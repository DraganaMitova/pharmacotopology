from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    load_hierarchical_gate_inputs,
)
from pharmacotopology.folding_axis_adjudication import (  # noqa: E402
    axis_adjudication_rows,
    axis_conflict_rows,
    axis_confusion_matrix_rows,
    axis_manual_review_rows,
    build_axis_adjudication_report,
    write_axis_adjudication_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_50.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_50_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_adjudication_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_rows.csv"
)
DEFAULT_CONFLICTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_conflicts.csv"
)
DEFAULT_MANUAL_REVIEW_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_manual_review.csv"
)
DEFAULT_CONFUSION_MATRICES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_confusion_matrices.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_dashboard.html"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Adjudicate collapsed fold labels into orthogonal truth axes over "
            "the locked 50-row regime-routed benchmark."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--conflicts-output", default=str(DEFAULT_CONFLICTS_PATH))
    parser.add_argument(
        "--manual-review-output",
        default=str(DEFAULT_MANUAL_REVIEW_PATH),
    )
    parser.add_argument(
        "--confusion-matrices-output",
        default=str(DEFAULT_CONFUSION_MATRICES_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_hierarchical_gate_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = axis_adjudication_rows(references, evidence_rows)
    conflicts = axis_conflict_rows(rows)
    manual_review_rows = axis_manual_review_rows(rows)
    confusion_rows = axis_confusion_matrix_rows(rows)
    report = build_axis_adjudication_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_axis_adjudication_outputs(
        report=report,
        rows=rows,
        conflicts=conflicts,
        manual_review_rows=manual_review_rows,
        confusion_rows=confusion_rows,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        conflicts_path=Path(args.conflicts_output),
        manual_review_path=Path(args.manual_review_output),
        confusion_matrices_path=Path(args.confusion_matrices_output),
        dashboard_path=Path(args.dashboard_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
