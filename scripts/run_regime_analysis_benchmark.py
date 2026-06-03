from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    load_hierarchical_gate_inputs,
)
from pharmacotopology.folding_regime_analysis import (  # noqa: E402
    abstention_analysis_rows,
    build_regime_analysis_report,
    failure_cohort_rows,
    high_confidence_wrong_rows,
    regime_analysis_rows,
    write_regime_analysis_outputs,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_50.locked.json")
DEFAULT_STRUCTURE_EVIDENCE_FILE = Path(
    "data/folding_benchmarks_real_50_structure_evidence.json"
)
DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_regime_analysis_report.json"
)
DEFAULT_ROWS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_regime_rows.csv"
)
DEFAULT_FAILURE_COHORTS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_failure_cohorts.csv"
)
DEFAULT_HIGH_CONFIDENCE_WRONG_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_high_confidence_wrong.csv"
)
DEFAULT_ABSTENTION_ANALYSIS_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_abstention_analysis.csv"
)
DEFAULT_DASHBOARD_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_50_regime_dashboard.html"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run sequence-only protein regime routing and failure-cohort "
            "analysis over the locked 50-row folding benchmark."
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
        "--failure-cohorts-output",
        default=str(DEFAULT_FAILURE_COHORTS_PATH),
    )
    parser.add_argument(
        "--high-confidence-wrong-output",
        default=str(DEFAULT_HIGH_CONFIDENCE_WRONG_PATH),
    )
    parser.add_argument(
        "--abstention-analysis-output",
        default=str(DEFAULT_ABSTENTION_ANALYSIS_PATH),
    )
    parser.add_argument("--dashboard-output", default=str(DEFAULT_DASHBOARD_PATH))
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    structure_evidence_file = Path(args.structure_evidence_file)
    references, evidence_rows = load_hierarchical_gate_inputs(
        benchmark_file,
        structure_evidence_file,
    )
    rows = regime_analysis_rows(references, evidence_rows)
    cohorts = failure_cohort_rows(rows)
    high_confidence_wrong = high_confidence_wrong_rows(rows)
    abstentions = abstention_analysis_rows(rows)
    report = build_regime_analysis_report(
        references,
        evidence_rows,
        source_benchmark_file=benchmark_file,
        structure_evidence_file=structure_evidence_file,
    )
    outputs = write_regime_analysis_outputs(
        report=report,
        rows=rows,
        cohort_rows=cohorts,
        high_confidence_wrong=high_confidence_wrong,
        abstention_rows=abstentions,
        report_path=Path(args.report_output),
        rows_path=Path(args.rows_output),
        failure_cohorts_path=Path(args.failure_cohorts_output),
        high_confidence_wrong_path=Path(args.high_confidence_wrong_output),
        abstention_analysis_path=Path(args.abstention_analysis_output),
        dashboard_path=Path(args.dashboard_output),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
