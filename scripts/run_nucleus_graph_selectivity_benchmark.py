from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_nucleus_graph_selectivity import (  # noqa: E402
    run_nucleus_graph_selectivity_benchmark,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_report.json"
)
DEFAULT_GRAPHS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_graphs.csv"
)
DEFAULT_SELECTED_EVENTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_selected_events.csv"
)
DEFAULT_REJECTIONS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_rejections.csv"
)
DEFAULT_DECOYS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_decoys.csv"
)
DEFAULT_RANK_ENRICHMENT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_rank_enrichment.csv"
)
DEFAULT_METRICS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_metrics.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/nucleus_graph_selectivity_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Score nucleus closure graphs and falsify them against matched "
            "decoys. This is a graph-selectivity audit, not a folding solution."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--graphs-output", default=str(DEFAULT_GRAPHS_PATH))
    parser.add_argument(
        "--selected-events-output",
        default=str(DEFAULT_SELECTED_EVENTS_PATH),
    )
    parser.add_argument("--rejections-output", default=str(DEFAULT_REJECTIONS_PATH))
    parser.add_argument("--decoys-output", default=str(DEFAULT_DECOYS_PATH))
    parser.add_argument(
        "--rank-enrichment-output",
        default=str(DEFAULT_RANK_ENRICHMENT_PATH),
    )
    parser.add_argument("--metrics-output", default=str(DEFAULT_METRICS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    outputs = run_nucleus_graph_selectivity_benchmark(
        benchmark_file=Path(args.benchmark_file),
        report_path=Path(args.report_output),
        graphs_path=Path(args.graphs_output),
        selected_events_path=Path(args.selected_events_output),
        rejections_path=Path(args.rejections_output),
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

