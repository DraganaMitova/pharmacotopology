from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    run_folding_nucleus_closure_search,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_report.json"
)
DEFAULT_EVENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_events.csv"
)
DEFAULT_TRAJECTORY_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_trajectory.csv"
)
DEFAULT_METRICS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_metrics.csv"
)
DEFAULT_FAILURES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_failures.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_nucleus_closure_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Search sequence-only cooperative segment closure events against "
            "the real-coordinate visual benchmark. This tests nucleus-like "
            "closure behavior, not a solved folding law."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--events-output", default=str(DEFAULT_EVENTS_PATH))
    parser.add_argument("--trajectory-output", default=str(DEFAULT_TRAJECTORY_PATH))
    parser.add_argument("--metrics-output", default=str(DEFAULT_METRICS_PATH))
    parser.add_argument("--failures-output", default=str(DEFAULT_FAILURES_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_folding_nucleus_closure_search(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        events_path=Path(args.events_output),
        trajectory_path=Path(args.trajectory_output),
        metrics_path=Path(args.metrics_output),
        failures_path=Path(args.failures_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
