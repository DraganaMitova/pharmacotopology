from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_architecture_axis import (  # noqa: E402
    architecture_axis_abstention_rows,
    architecture_axis_conflict_rows,
    architecture_axis_rows,
    build_architecture_axis_report,
    write_architecture_axis_outputs,
)
from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    load_hierarchical_gate_inputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_50.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_50_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_rows.csv"
)
DEFAULT_CONFLICTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_conflicts.csv"
)
DEFAULT_ABSTENTIONS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_abstentions.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Adjudicate sequence-only architecture-axis evidence over the "
            "locked 50-row folding benchmark."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--conflicts-output", default=str(DEFAULT_CONFLICTS_PATH))
    parser.add_argument(
        "--abstentions-output",
        default=str(DEFAULT_ABSTENTIONS_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_hierarchical_gate_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = architecture_axis_rows(references, evidence_rows)
    report = build_architecture_axis_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_architecture_axis_outputs(
        report=report,
        rows=rows,
        conflicts=architecture_axis_conflict_rows(rows),
        abstentions=architecture_axis_abstention_rows(rows),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        conflicts_path=Path(args.conflicts_output),
        abstentions_path=Path(args.abstentions_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
