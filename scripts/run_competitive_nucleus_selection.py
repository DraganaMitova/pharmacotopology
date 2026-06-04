from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_nucleus_competition import (  # noqa: E402
    run_competitive_nucleus_selection,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_report.json"
)
DEFAULT_SELECTED_EVENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_selected_events.csv"
)
DEFAULT_REJECTIONS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_rejections.csv"
)
DEFAULT_COMPATIBILITY_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_compatibility.csv"
)
DEFAULT_TRAJECTORY_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_trajectory.csv"
)
DEFAULT_METRICS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_metrics.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/competitive_nucleus_selection_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Apply competitive closure selection and frustration filtering to "
            "the coordinate-native nucleus benchmark. This is an audit layer, "
            "not a protein-folding solution."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument(
        "--selected-events-output",
        default=str(DEFAULT_SELECTED_EVENTS_PATH),
    )
    parser.add_argument("--rejections-output", default=str(DEFAULT_REJECTIONS_PATH))
    parser.add_argument(
        "--compatibility-output",
        default=str(DEFAULT_COMPATIBILITY_PATH),
    )
    parser.add_argument("--trajectory-output", default=str(DEFAULT_TRAJECTORY_PATH))
    parser.add_argument("--metrics-output", default=str(DEFAULT_METRICS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_competitive_nucleus_selection(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        selected_events_path=Path(args.selected_events_output),
        rejections_path=Path(args.rejections_output),
        compatibility_path=Path(args.compatibility_output),
        trajectory_path=Path(args.trajectory_output),
        metrics_path=Path(args.metrics_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()

