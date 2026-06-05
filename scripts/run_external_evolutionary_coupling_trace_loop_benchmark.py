from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    run_external_evolutionary_coupling_trace_loop_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_ORACLE_COUPLING_FILE = Path(
    "data/folding_real_coordinate_visual_8_couplings.locked.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_report.json"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_certificate.json"
)
DEFAULT_SELECTORS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_selectors.csv"
)
DEFAULT_SELECTED_EVENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_selected_events.csv"
)
DEFAULT_FRONTIER_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_frontier.csv"
)
DEFAULT_CONTROLS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_controls.csv"
)
DEFAULT_ROW_STATUS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_row_status.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_dashboard.html"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_V0 against a "
            "provenance-locked external MSA/DCA coupling file."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", required=True)
    parser.add_argument(
        "--oracle-coupling-file",
        default=str(DEFAULT_ORACLE_COUPLING_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    parser.add_argument("--selectors-output", default=str(DEFAULT_SELECTORS_PATH))
    parser.add_argument(
        "--selected-events-output",
        default=str(DEFAULT_SELECTED_EVENTS_PATH),
    )
    parser.add_argument("--frontier-output", default=str(DEFAULT_FRONTIER_PATH))
    parser.add_argument("--controls-output", default=str(DEFAULT_CONTROLS_PATH))
    parser.add_argument("--row-status-output", default=str(DEFAULT_ROW_STATUS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    args = parser.parse_args()

    outputs = run_external_evolutionary_coupling_trace_loop_benchmark(
        benchmark_file=Path(args.benchmark_file),
        external_coupling_file=Path(args.external_coupling_file),
        oracle_coupling_file=Path(args.oracle_coupling_file),
        report_path=Path(args.report_output),
        certificate_path=Path(args.certificate_output),
        selectors_path=Path(args.selectors_output),
        selected_events_path=Path(args.selected_events_output),
        frontier_path=Path(args.frontier_output),
        controls_path=Path(args.controls_output),
        row_status_path=Path(args.row_status_output),
        dashboard_path=Path(args.dashboard_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
