from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from pharmacotopology.layer import run_clean_pharmacotopology_layer
from render_pharmacotopology_dashboard import render_dashboard


def test_pharmacotopology_dashboard_creates_html(tmp_path: Path) -> None:
    run_clean_pharmacotopology_layer(tmp_path)

    output_path = render_dashboard(
        tmp_path / "memory.jsonl",
        tmp_path / "pharmacotopology_dashboard.html",
    )

    assert output_path.exists()
    assert output_path.suffix == ".html"


def test_pharmacotopology_dashboard_contains_required_labels(tmp_path: Path) -> None:
    run_clean_pharmacotopology_layer(tmp_path)
    output_path = render_dashboard(
        tmp_path / "memory.jsonl",
        tmp_path / "pharmacotopology_dashboard.html",
    )
    html = output_path.read_text(encoding="utf-8")

    assert "Φ.review" in html
    assert "SIMULATION ONLY" in html
    assert "NOT MEDICAL ADVICE" in html
    assert "pathology_reduction_score" in html
    assert "collapse_cost_score" in html
    assert "net_topology_health_score" in html


def test_pharmacotopology_dashboard_avoids_prescribing_language(
    tmp_path: Path,
) -> None:
    run_clean_pharmacotopology_layer(tmp_path)
    output_path = render_dashboard(
        tmp_path / "memory.jsonl",
        tmp_path / "pharmacotopology_dashboard.html",
    )
    html = output_path.read_text(encoding="utf-8").lower()

    forbidden_phrases = (
        "dosage",
        "dose",
        "prescription",
        "prescribe",
        "patient use",
        "real-world prescribing",
        "treatment recommendation",
        "medication guidance",
        "pill works",
        "pill does not work",
    )
    for phrase in forbidden_phrases:
        assert phrase not in html
