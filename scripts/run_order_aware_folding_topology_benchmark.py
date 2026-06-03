from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_order_aware import (  # noqa: E402
    build_order_aware_report,
    contact_prior_rows,
    load_order_aware_inputs,
    order_aware_benchmark_rows,
    order_aware_control_separation_rows,
    write_order_aware_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_10_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_order_aware_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_order_aware_rows.csv"
)
DEFAULT_CONTACT_PRIOR_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_contact_prior.csv"
)
DEFAULT_CONTROL_SEPARATION_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_control_separation.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_10_order_aware_dashboard.html"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the order-aware folding topology recipe and sequence-order "
            "falsification controls against the locked 10-row benchmark."
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
        "--contact-prior-output",
        default=str(DEFAULT_CONTACT_PRIOR_PATH),
    )
    parser.add_argument(
        "--control-separation-output",
        default=str(DEFAULT_CONTROL_SEPARATION_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_order_aware_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = order_aware_benchmark_rows(references, evidence_rows)
    controls = order_aware_control_separation_rows(references)
    contacts = contact_prior_rows(references)
    report = build_order_aware_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_order_aware_outputs(
        report=report,
        rows=rows,
        contact_rows=contacts,
        control_rows=controls,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        contact_prior_path=Path(args.contact_prior_output),
        control_separation_path=Path(args.control_separation_output),
        dashboard_path=Path(args.dashboard_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
