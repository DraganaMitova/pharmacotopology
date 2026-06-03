from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_motif_alignment import (  # noqa: E402
    build_motif_alignment_report,
    evidence_conflict_rows,
    failure_diagnosis_rows,
    load_motif_alignment_inputs,
    motif_alignment_rows,
    write_motif_alignment_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_10_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_motif_alignment_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_motif_alignment_rows.csv"
)
DEFAULT_FAILURE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_failure_diagnosis.csv"
)
DEFAULT_CONFLICTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_evidence_conflicts.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_motif_alignment_dashboard.html"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the motif-to-structure alignment layer against the locked "
            "10-row folding benchmark and emit failure-diagnosis artifacts."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument("--failure-output", default=str(DEFAULT_FAILURE_PATH))
    parser.add_argument("--conflicts-output", default=str(DEFAULT_CONFLICTS_PATH))
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_motif_alignment_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = motif_alignment_rows(references, evidence_rows)
    failure_rows = failure_diagnosis_rows(rows)
    conflict_rows = evidence_conflict_rows(rows)
    report = build_motif_alignment_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_motif_alignment_outputs(
        report=report,
        rows=rows,
        failure_rows=failure_rows,
        conflict_rows=conflict_rows,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        failure_path=Path(args.failure_output),
        conflicts_path=Path(args.conflicts_output),
        dashboard_path=Path(args.dashboard_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
