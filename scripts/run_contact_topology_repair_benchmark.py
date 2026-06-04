from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_contact_topology_repair_benchmark import (  # noqa: E402
    run_contact_topology_repair_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_mechanism_visual_12.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_rows.csv"
)
DEFAULT_GAP_ANALYSIS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_gap_analysis.csv"
)
DEFAULT_FAILURE_COHORTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_failure_cohorts.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_certificate.json"
)
DEFAULT_VISUALS_ROOT = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_repair_visuals"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run sequence-only contact-topology repair against the locked "
            "12-row visual mechanism benchmark, then score native gaps after "
            "prediction."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument(
        "--gap-analysis-output",
        default=str(DEFAULT_GAP_ANALYSIS_PATH),
    )
    parser.add_argument(
        "--failure-cohorts-output",
        default=str(DEFAULT_FAILURE_COHORTS_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    parser.add_argument("--visuals-root", default=str(DEFAULT_VISUALS_ROOT))
    args = parser.parse_args()

    outputs = run_contact_topology_repair_benchmark(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        gap_analysis_path=Path(args.gap_analysis_output),
        failure_cohorts_path=Path(args.failure_cohorts_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
        visuals_root=Path(args.visuals_root),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
