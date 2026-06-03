from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_external_axis_repair import (  # noqa: E402
    build_external_axis_repair_report,
    external_axis_repair_abstention_delta_rows,
    external_axis_repair_conflict_delta_rows,
    external_axis_repair_family_summary_rows,
    external_axis_repair_quarantine_rows,
    external_axis_repair_rows,
    write_external_axis_repair_outputs,
)
from pharmacotopology.folding_external_holdout import (  # noqa: E402
    DEFAULT_DEVELOPMENT_BENCHMARK_FILE,
    load_external_holdout_rows,
)


DEFAULT_HOLDOUT_FILE = Path(
    "data/folding_benchmarks_external_fold_family_100.locked.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_rows.csv"
)
DEFAULT_CONFLICT_DELTA_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_conflict_delta.csv"
)
DEFAULT_ABSTENTION_DELTA_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_abstention_delta.csv"
)
DEFAULT_QUARANTINE_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_quarantine_rows.csv"
)
DEFAULT_FAMILY_SUMMARY_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_family_summary.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_axis_repair_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Apply external-safe order and architecture quarantine to the "
            "locked external fold-family holdout."
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
        "--conflict-delta-output",
        default=str(DEFAULT_CONFLICT_DELTA_PATH),
    )
    parser.add_argument(
        "--abstention-delta-output",
        default=str(DEFAULT_ABSTENTION_DELTA_PATH),
    )
    parser.add_argument(
        "--quarantine-rows-output",
        default=str(DEFAULT_QUARANTINE_ROWS_PATH),
    )
    parser.add_argument(
        "--family-summary-output",
        default=str(DEFAULT_FAMILY_SUMMARY_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    holdout_file = Path(args.holdout_file)
    development_benchmark_file = Path(args.development_benchmark_file)
    holdout_rows = load_external_holdout_rows(holdout_file)
    rows = external_axis_repair_rows(holdout_rows)
    report = build_external_axis_repair_report(
        holdout_rows,
        holdout_file=holdout_file,
        development_benchmark_file=development_benchmark_file,
    )
    outputs = write_external_axis_repair_outputs(
        report=report,
        rows=rows,
        conflict_delta_rows=external_axis_repair_conflict_delta_rows(report),
        abstention_delta_rows=external_axis_repair_abstention_delta_rows(report),
        quarantine_rows=external_axis_repair_quarantine_rows(rows),
        family_rows=external_axis_repair_family_summary_rows(rows),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        conflict_delta_path=Path(args.conflict_delta_output),
        abstention_delta_path=Path(args.abstention_delta_output),
        quarantine_rows_path=Path(args.quarantine_rows_output),
        family_summary_path=Path(args.family_summary_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
