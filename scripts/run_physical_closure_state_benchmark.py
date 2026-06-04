from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_closure_state_builder import (  # noqa: E402
    run_physical_closure_state_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_report.json"
)
DEFAULT_STATES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_states.csv"
)
DEFAULT_DECOYS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_decoys.csv"
)
DEFAULT_RANK_ENRICHMENT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_rank_enrichment.csv"
)
DEFAULT_METRICS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_metrics.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/physical_closure_state_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Instantiate graph-selected closure events as coarse physical "
            "states and compare them with matched decoys. This is not an "
            "atomistic folding simulation."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--states-output", default=str(DEFAULT_STATES_PATH))
    parser.add_argument("--decoys-output", default=str(DEFAULT_DECOYS_PATH))
    parser.add_argument(
        "--rank-enrichment-output",
        default=str(DEFAULT_RANK_ENRICHMENT_PATH),
    )
    parser.add_argument("--metrics-output", default=str(DEFAULT_METRICS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_physical_closure_state_benchmark(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        states_path=Path(args.states_output),
        decoys_path=Path(args.decoys_output),
        rank_enrichment_path=Path(args.rank_enrichment_output),
        metrics_path=Path(args.metrics_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()

