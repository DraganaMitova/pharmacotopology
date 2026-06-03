from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from pharmacotopology.folding_metrics import (  # noqa: E402
    evidence_readiness_summary,
    fold_class_match_rate,
    mean_contact_map_similarity,
    summarize_benchmark,
)
from pharmacotopology.folding_topology import (  # noqa: E402
    run_folding_topology_benchmark,
)
from run_folding_topology_benchmark import (  # noqa: E402
    write_folding_benchmark_outputs,
)


def test_folding_metrics_stay_bounded() -> None:
    comparisons = run_folding_topology_benchmark()
    summary = summarize_benchmark(comparisons)

    assert 0.0 <= mean_contact_map_similarity(comparisons) <= 1.0
    assert 0.0 <= fold_class_match_rate(comparisons) <= 1.0
    assert 0.0 <= summary["mean_uncertainty_radius"] <= 1.0
    assert evidence_readiness_summary(comparisons) == {
        "benchmark_shell_only": len(comparisons)
    }


def test_folding_benchmark_writes_json_and_csv(tmp_path: Path) -> None:
    comparisons = run_folding_topology_benchmark()
    report_path, csv_path = write_folding_benchmark_outputs(
        comparisons,
        tmp_path / "folding_report.json",
        tmp_path / "folding.csv",
    )

    report_text = report_path.read_text(encoding="utf-8")
    csv_text = csv_path.read_text(encoding="utf-8")

    assert "protein_folding_topology_hypothesis_benchmark" in report_text
    assert "folding_problem_solved" in report_text
    assert "predicted_topology_signature" in csv_text
    assert "reference_topology_signature" in csv_text
    assert "benchmark_shell_only" in csv_text
