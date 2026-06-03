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
from explore_sensitivity import explore_sensitivity, write_explorer_outputs
from pharmacotopology.layer import get_topology_profile, run_clean_pharmacotopology_layer
from render_profile_comparison_dashboard import render_profile_comparison_dashboard
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
    assert "primary_evidence_sources" in rankings
    assert "calibration_blockers" in rankings
    assert "net_topology_health_lower" in rankings
    assert "dimension" in deltas
    assert "collapse_cost" in deltas


def test_profile_comparison_dashboard_contains_cross_profile_sections(
    tmp_path: Path,
) -> None:
    output_path = render_profile_comparison_dashboard(tmp_path / "multi.html")
    html = output_path.read_text(encoding="utf-8")

    assert "Multi-Profile Pharmacotopology Dashboard" in html
    assert "Baseline Topology Maps" in html
    assert "Combined Ranking Table" in html
    assert "Mechanism/Profile Heatmap" in html
    assert "Cross-Profile Summary" in html
    assert "mixed_state_like" in html
    assert "NOT MEDICAL ADVICE" in html


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


def test_sensitivity_explorer_writes_distribution_and_samples(
    tmp_path: Path,
) -> None:
    report = explore_sensitivity(
        profile_key="schizophrenia_like",
        mechanism_id="nmda_support_like",
        vary_field="collapse_cost",
        vary_range=(0.05, 0.25),
        samples=12,
        noise=0.02,
        seed=7,
    )

    assert report["clinical_use_allowed"] is False
    assert report["practical_use"] == "assumption_sensitivity_exploration"
    assert report["target_distribution"]["net_score_p05"] <= (
        report["target_distribution"]["net_score_p95"]
    )
    assert len(report["samples_table"]) == 12
    assert report["robustness"]

    report_path, samples_path = write_explorer_outputs(
        report,
        tmp_path / "explorer.json",
        tmp_path / "samples.csv",
    )

    assert report_path.exists()
    assert samples_path.exists()
    assert "varied_value" in samples_path.read_text(encoding="utf-8")
