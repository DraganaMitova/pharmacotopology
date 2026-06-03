from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_structure_benchmark import (  # noqa: E402
    build_structure_benchmark_report,
    load_real10_with_structure_evidence,
    sequence_order_control_rows,
    write_structure_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_10_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_structure_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_structure_rows.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_structure_dashboard.html"
)
DEFAULT_ORDER_CONTROLS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_order_controls.csv"
)
DEFAULT_FALSIFICATION_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_falsification_report.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the structure-derived folding topology benchmark and "
            "sequence-order falsification controls."
        )
    )
    parser.add_argument(
        "--benchmark-file",
        default=str(DEFAULT_BENCHMARK_FILE),
        help="Locked external benchmark file.",
    )
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
        help="Structure evidence JSON produced by extract_structure_topology_signatures.py.",
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument(
        "--order-controls-output",
        default=str(DEFAULT_ORDER_CONTROLS_PATH),
    )
    parser.add_argument(
        "--falsification-output",
        default=str(DEFAULT_FALSIFICATION_REPORT_PATH),
    )
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_real10_with_structure_evidence(
        benchmark_file,
        structure_evidence_file,
    )
    report = build_structure_benchmark_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    control_rows = sequence_order_control_rows(references)
    outputs = write_structure_outputs(
        report=report,
        control_rows=control_rows,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        dashboard_path=Path(args.dashboard_output),
        order_controls_path=Path(args.order_controls_output),
        falsification_report_path=Path(args.falsification_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
