from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_external_holdout import (  # noqa: E402
    DEFAULT_DEVELOPMENT_BENCHMARK_FILE,
    abstention_rows,
    axis_conflict_rows,
    build_external_holdout_report,
    external_holdout_rows,
    failure_cohort_rows,
    family_summary_rows,
    load_external_holdout_rows,
    write_external_holdout_outputs,
)


DEFAULT_HOLDOUT_FILE = Path(
    "data/folding_benchmarks_external_fold_family_100.locked.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_rows.csv"
)
DEFAULT_FAMILY_SUMMARY_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_family_summary.csv"
)
DEFAULT_CONFLICTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_axis_conflicts.csv"
)
DEFAULT_ABSTENTIONS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_abstentions.csv"
)
DEFAULT_FAILURE_COHORTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_failure_cohorts.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate the current safe folding-axis stack on the locked "
            "external fold-family holdout without tuning thresholds."
        )
    )
    parser.add_argument("--holdout-file", default=str(DEFAULT_HOLDOUT_FILE))
    parser.add_argument(
        "--development-benchmark-file",
        default=str(DEFAULT_DEVELOPMENT_BENCHMARK_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument(
        "--family-summary-output",
        default=str(DEFAULT_FAMILY_SUMMARY_PATH),
    )
    parser.add_argument("--axis-conflicts-output", default=str(DEFAULT_CONFLICTS_PATH))
    parser.add_argument("--abstentions-output", default=str(DEFAULT_ABSTENTIONS_PATH))
    parser.add_argument(
        "--failure-cohorts-output",
        default=str(DEFAULT_FAILURE_COHORTS_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    holdout_file = Path(args.holdout_file)
    development_benchmark_file = Path(args.development_benchmark_file)
    holdout_rows = load_external_holdout_rows(holdout_file)
    rows = external_holdout_rows(holdout_rows)
    report = build_external_holdout_report(
        holdout_rows,
        holdout_file=holdout_file,
        development_benchmark_file=development_benchmark_file,
    )
    outputs = write_external_holdout_outputs(
        report=report,
        rows=rows,
        family_rows=family_summary_rows(rows),
        conflicts=axis_conflict_rows(rows),
        abstentions=abstention_rows(rows),
        failure_cohorts=failure_cohort_rows(rows),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        family_summary_path=Path(args.family_summary_output),
        axis_conflicts_path=Path(args.axis_conflicts_output),
        abstentions_path=Path(args.abstentions_output),
        failure_cohorts_path=Path(args.failure_cohorts_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
