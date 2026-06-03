from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    build_hierarchical_gate_report,
    gate_failure_rows,
    gate_path_rows,
    hierarchical_gate_rows,
    load_hierarchical_gate_inputs,
    write_hierarchical_gate_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_10_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_rows.csv"
)
DEFAULT_GATE_PATHS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_gate_paths.csv"
)
DEFAULT_GATE_FAILURES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_gate_failures.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_hierarchical_gate_dashboard.html"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the hierarchical folding decision gates against the locked "
            "10-row folding benchmark and emit gate path diagnostics."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--gate-paths-output", default=str(DEFAULT_GATE_PATHS_PATH))
    parser.add_argument(
        "--gate-failures-output",
        default=str(DEFAULT_GATE_FAILURES_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_hierarchical_gate_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = hierarchical_gate_rows(references, evidence_rows)
    paths = gate_path_rows(rows)
    failures = gate_failure_rows(rows)
    report = build_hierarchical_gate_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_hierarchical_gate_outputs(
        report=report,
        rows=rows,
        path_rows=paths,
        failure_rows=failures,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        gate_paths_path=Path(args.gate_paths_output),
        gate_failures_path=Path(args.gate_failures_output),
        dashboard_path=Path(args.dashboard_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
