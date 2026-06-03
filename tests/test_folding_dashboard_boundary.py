from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from pharmacotopology.folding_topology import (  # noqa: E402
    run_folding_topology_benchmark,
)
from render_folding_benchmark_dashboard import (  # noqa: E402
    render_folding_benchmark_dashboard,
)
from run_folding_topology_benchmark import (  # noqa: E402
    write_folding_benchmark_outputs,
)


def test_folding_benchmark_dashboard_contains_proof_sections(
    tmp_path: Path,
) -> None:
    comparisons = run_folding_topology_benchmark()
    report_path, csv_path = write_folding_benchmark_outputs(
        comparisons,
        tmp_path / "real_folding_500_report.json",
        tmp_path / "real_folding_500_rows.csv",
    )
    output_path = render_folding_benchmark_dashboard(
        report_path,
        csv_path,
        tmp_path / "real_folding_500_dashboard.html",
    )
    html = output_path.read_text(encoding="utf-8")

    assert "Folding Benchmark Dashboard" in html
    assert "NOT A FOLDING SOLUTION" in html
    assert "Confusion Matrix" in html
    assert "Similarity Distribution" in html
    assert "Per-Class Accuracy" in html
    assert "Failure Table" in html
    assert "Radar Overlay" in html
    assert "Locked Benchmark Certificate" in html
    assert "NOT A DRUG-DESIGN TOOL" in html
