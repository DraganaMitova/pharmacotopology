from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from export_pharmacotopology_csv import export_csv, read_phi_packet
from pharmacotopology.layer import get_topology_profile, run_clean_pharmacotopology_layer
from run_sensitivity_analysis import (
    run_sensitivity_analysis,
    write_sensitivity_outputs,
)


def test_csv_export_writes_rankings_and_deltas(tmp_path: Path) -> None:
    run_clean_pharmacotopology_layer(tmp_path)
    packet = read_phi_packet(tmp_path / "memory.jsonl")

    rankings_path, deltas_path = export_csv(
        packet["Φ.review"],
        tmp_path / "rankings.csv",
        tmp_path / "deltas.csv",
    )

    rankings = rankings_path.read_text(encoding="utf-8")
    deltas = deltas_path.read_text(encoding="utf-8")

    assert "mechanism_id" in rankings
    assert "evidence_readiness_label" in rankings
    assert "net_topology_health_lower" in rankings
    assert "dimension" in deltas
    assert "collapse_cost" in deltas


def test_sensitivity_analysis_writes_json_and_csv(tmp_path: Path) -> None:
    report = run_sensitivity_analysis(
        profile=get_topology_profile("anxiety_like"),
        pressure_step=0.05,
        dimensions=("threat_propagation",),
    )

    assert report["clinical_use_allowed"] is False
    assert report["practical_use"] == "ranking_robustness_review"
    assert len(report["variants"]) == 2
    assert report["variants"][0]["dimension"] == "threat_propagation"

    report_path, csv_path = write_sensitivity_outputs(
        report,
        tmp_path / "sensitivity.json",
        tmp_path / "sensitivity.csv",
    )

    assert report_path.exists()
    assert csv_path.exists()
    assert "top_mechanism_id" in csv_path.read_text(encoding="utf-8")
