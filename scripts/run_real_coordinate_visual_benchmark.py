from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    run_real_coordinate_visual_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visual_8_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visual_8_rows.csv"
)
DEFAULT_CONTACT_METRICS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visual_8_contact_metrics.csv"
)
DEFAULT_NATIVE_CONTACT_SUMMARY_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visual_8_native_contact_summary.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visual_8_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visual_8_certificate.json"
)
DEFAULT_VISUALS_ROOT = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_coordinate_visuals"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the real-coordinate visual contact benchmark. Native maps are "
            "derived from locked C-alpha coordinate traces after sequence-only "
            "prediction."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument(
        "--contact-metrics-output",
        default=str(DEFAULT_CONTACT_METRICS_PATH),
    )
    parser.add_argument(
        "--native-contact-summary-output",
        default=str(DEFAULT_NATIVE_CONTACT_SUMMARY_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    parser.add_argument("--visuals-root", default=str(DEFAULT_VISUALS_ROOT))
    args = parser.parse_args()

    outputs = run_real_coordinate_visual_benchmark(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        contact_metrics_path=Path(args.contact_metrics_output),
        native_contact_summary_path=Path(args.native_contact_summary_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
        visuals_root=Path(args.visuals_root),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
