from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_physical_selection import (  # noqa: E402
    run_active_physical_selection_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/active_physical_selection_report.json"
)
DEFAULT_SELECTORS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/active_physical_selection_selectors.csv"
)
DEFAULT_SELECTED_EVENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/active_physical_selection_selected_events.csv"
)
DEFAULT_ABLATION_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/active_physical_selection_ablation.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/active_physical_selection_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/active_physical_selection_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Use coarse physical closure-state terms as active selection "
            "controls, then ablate each term. This is a bounded falsification "
            "benchmark, not a folding mechanism discovery claim."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--selectors-output", default=str(DEFAULT_SELECTORS_PATH))
    parser.add_argument(
        "--selected-events-output",
        default=str(DEFAULT_SELECTED_EVENTS_PATH),
    )
    parser.add_argument("--ablation-output", default=str(DEFAULT_ABLATION_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_active_physical_selection_benchmark(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        selectors_path=Path(args.selectors_output),
        selected_events_path=Path(args.selected_events_output),
        ablation_path=Path(args.ablation_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
