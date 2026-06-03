from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_axis_profile import (  # noqa: E402
    axis_profile_abstention_rows,
    axis_profile_recovery_candidate_rows,
    axis_profile_rows,
    build_axis_profile_report,
    write_axis_profile_outputs,
)
from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    load_hierarchical_gate_inputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_50.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_50_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_rows.csv"
)
DEFAULT_ABSTENTIONS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_abstentions.csv"
)
DEFAULT_RECOVERY_CANDIDATES_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_recovery_candidates.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Recover safe partial fold-axis coverage from the locked 50-row "
            "regime-routed benchmark without recovering collapsed fold classes."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--structure-evidence-file",
        default=str(DEFAULT_STRUCTURE_EVIDENCE_FILE),
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument(
        "--abstentions-output",
        default=str(DEFAULT_ABSTENTIONS_PATH),
    )
    parser.add_argument(
        "--recovery-candidates-output",
        default=str(DEFAULT_RECOVERY_CANDIDATES_PATH),
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
    rows = axis_profile_rows(references, evidence_rows)
    report = build_axis_profile_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_axis_profile_outputs(
        report=report,
        rows=rows,
        abstention_rows=axis_profile_abstention_rows(rows),
        recovery_candidate_rows=axis_profile_recovery_candidate_rows(rows),
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        abstentions_path=Path(args.abstentions_output),
        recovery_candidates_path=Path(args.recovery_candidates_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
