from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_visual_mechanism_audit import (  # noqa: E402
    build_visual_mechanism_audit_report,
    write_visual_mechanism_audit_outputs,
)


DEFAULT_BASELINE_REPORT = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_12_report.json"
)
DEFAULT_REPAIR_REPORT = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_report.json"
)
DEFAULT_SOURCE_BENCHMARK = Path("data/folding_mechanism_visual_12.locked.json")
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_audit_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_audit_rows.csv"
)
DEFAULT_OVERFIT_RISKS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_audit_overfit_risks.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_audit_dashboard.html"
)
DEFAULT_CERTIFICATE_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_audit_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Audit visual folding mechanism claims, mark the 12-row benchmark "
            "as toy/coarse/internal, and report contact-repair overfit risk."
        )
    )
    parser.add_argument("--baseline-report", default=str(DEFAULT_BASELINE_REPORT))
    parser.add_argument("--repair-report", default=str(DEFAULT_REPAIR_REPORT))
    parser.add_argument("--source-benchmark", default=str(DEFAULT_SOURCE_BENCHMARK))
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_PATH))
    parser.add_argument("--rows-output", default=str(DEFAULT_ROWS_PATH))
    parser.add_argument(
        "--overfit-risks-output",
        default=str(DEFAULT_OVERFIT_RISKS_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    parser.add_argument("--certificate-output", default=str(DEFAULT_CERTIFICATE_PATH))
    args = parser.parse_args()

    report = build_visual_mechanism_audit_report(
        baseline_report_path=Path(args.baseline_report),
        repair_report_path=Path(args.repair_report),
        source_benchmark_file=Path(args.source_benchmark),
    )
    outputs = write_visual_mechanism_audit_outputs(
        report=report,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        overfit_risks_path=Path(args.overfit_risks_output),
        dashboard_path=Path(args.dashboard_output),
        certificate_path=Path(args.certificate_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
