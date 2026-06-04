from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_coupling_nucleus_selector import (  # noqa: E402
    run_coupling_nucleus_selector_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_COUPLING_FILE = Path(
    "data/folding_real_coordinate_visual_8_couplings.locked.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_report.json"
)
DEFAULT_SELECTORS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_selectors.csv"
)
DEFAULT_SELECTED_EVENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_selected_events.csv"
)
DEFAULT_ASSESSMENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_assessments.csv"
)
DEFAULT_DECOYS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_decoys.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Test whether a locked coupling-constraint channel can select "
            "native-compatible folding nuclei while preserving future closure "
            "paths. The checked-in coupling file is an oracle control unless "
            "replaced with external MSA/DCA couplings."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--selectors-output", default=str(DEFAULT_SELECTORS_PATH))
    parser.add_argument(
        "--selected-events-output",
        default=str(DEFAULT_SELECTED_EVENTS_PATH),
    )
    parser.add_argument("--assessments-output", default=str(DEFAULT_ASSESSMENTS_PATH))
    parser.add_argument("--decoys-output", default=str(DEFAULT_DECOYS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_coupling_nucleus_selector_benchmark(
        benchmark_file=Path(args.benchmark_file),
        coupling_file=Path(args.coupling_file),
        report_path=Path(args.report_output),
        selectors_path=Path(args.selectors_output),
        selected_events_path=Path(args.selected_events_output),
        assessments_path=Path(args.assessments_output),
        decoys_path=Path(args.decoys_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
