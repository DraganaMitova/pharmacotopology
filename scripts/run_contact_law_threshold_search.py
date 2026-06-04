from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_contact_threshold_search import (  # noqa: E402
    run_contact_law_threshold_search,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_rows.csv"
)
DEFAULT_GRID_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_grid.csv"
)
DEFAULT_HOLDOUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_holdout.csv"
)
DEFAULT_FAILURES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_failures.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_law_threshold_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Search sequence-only contact-law thresholds against the locked "
            "real-coordinate visual benchmark. This is a falsification bench, "
            "not a folding solver."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--grid-output", default=str(DEFAULT_GRID_PATH))
    parser.add_argument("--holdout-output", default=str(DEFAULT_HOLDOUT_PATH))
    parser.add_argument("--failures-output", default=str(DEFAULT_FAILURES_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_contact_law_threshold_search(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        grid_path=Path(args.grid_output),
        holdout_path=Path(args.holdout_output),
        failures_path=Path(args.failures_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
