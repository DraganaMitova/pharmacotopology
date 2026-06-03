import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_architecture_axis import (  # noqa: E402
    ARCHITECTURE_AXIS_BENCHMARK_KIND,
    ARCHITECTURE_AXIS_SIGNATURE_KIND,
    architecture_axis_abstention_rows,
    architecture_axis_conflict_rows,
    architecture_axis_rows,
    build_architecture_axis_report,
    write_architecture_axis_outputs,
)
from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    load_hierarchical_gate_inputs,
)


BENCHMARK_50 = ROOT / "data" / "folding_benchmarks_real_50.locked.json"
REL_BENCHMARK_50 = Path("data/folding_benchmarks_real_50.locked.json")
STRUCTURE_EVIDENCE_50 = (
    ROOT / "data" / "folding_benchmarks_real_50_structure_evidence.json"
)
REL_STRUCTURE_EVIDENCE_50 = Path(
    "data/folding_benchmarks_real_50_structure_evidence.json"
)
REPORT_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_architecture_axis_report.json"
)
ROWS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_architecture_axis_rows.csv"
)
CONFLICTS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_architecture_axis_conflicts.csv"
)
ABSTENTIONS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_architecture_axis_abstentions.csv"
)
DASHBOARD_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_architecture_axis_dashboard.html"
)
CERTIFICATE_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_architecture_axis_certificate.json"
)


def test_architecture_axis_report_recovers_safe_coverage_only() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    report = build_architecture_axis_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_50,
        structure_evidence_file=STRUCTURE_EVIDENCE_50,
    )

    assert report["benchmark_kind"] == ARCHITECTURE_AXIS_BENCHMARK_KIND
    assert report["architecture_axis_signature_kind"] == (
        ARCHITECTURE_AXIS_SIGNATURE_KIND
    )
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_labels_no_structure_answers"
    )
    assert report["benchmark_size"] == 50
    assert report["previous_profile_architecture_axis_coverage"] == 0.14
    assert report["architecture_axis_coverage"] == 0.32
    assert report["architecture_axis_coverage"] > report[
        "previous_profile_architecture_axis_coverage"
    ]
    assert report["architecture_axis_claim_allowed_count"] == 16
    assert report["architecture_axis_abstained_count"] == 34
    assert report["architecture_axis_conflict_count"] == 0
    assert report["architecture_axis_safe_claim_count"] == 16
    assert report["architecture_axis_same_axis_conflict_count"] == 0
    assert report["architecture_axis_accuracy"] == 1.0
    assert report["fragment_scope_detected_count"] == 7
    assert report["compact_single_domain_claim_count"] == 6
    assert report["multidomain_claim_count"] == 3
    assert report["multidomain_abstained_count"] == 9
    assert report["repeat_like_claim_count"] == 0
    assert report["architecture_claim_distribution"] == {
        "compact_single_domain": 6,
        "fragment_scope": 7,
        "multidomain_or_segmented": 3,
        "repeat_like": 0,
        "unknown": 0,
    }
    assert report["architecture_claim_without_secondary_leakage_count"] == 0
    assert report["architecture_claim_without_label_leakage_count"] == 0
    assert report["architecture_axis_recovered_from_profile_count"] == 7
    assert report["global_fold_class_claim_allowed"] is False
    assert report["axis_profile_claim_allowed"] is True
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False
    assert report["artifact_reproducible"] is True


def test_architecture_axis_rows_separate_architecture_from_other_axes() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    rows = {
        row["row_id"]: row
        for row in architecture_axis_rows(references, evidence)
    }

    fragment = rows["pdb_1R69_A_lambda_repressor_ntd"]
    assert fragment["sequence_hash"] == "4181a9e1c39ff426"
    assert fragment["architecture_axis_prediction"] == "fragment_scope"
    assert fragment["architecture_axis_claim_allowed"] is True
    assert fragment["architecture_axis_recovered_from_profile"] is True
    assert fragment["adjudicated_truth_architecture_axis"] == "fragment_scope"
    assert fragment["architecture_axis_same_axis_conflict"] is False

    compact = rows["pdb_1MBN_A_myoglobin"]
    assert compact["architecture_axis_prediction"] == "compact_single_domain"
    assert compact["architecture_axis_claim_allowed"] is True
    assert compact["profile_architecture_axis"] == "unknown"
    assert compact["architecture_axis_recovered_from_profile"] is False
    assert compact["architecture_secondary_leakage_used"] is False
    assert compact["architecture_label_leakage_used"] is False

    modular = rows["pdb_1IGT_B_igg2a_heavy_chain"]
    assert modular["architecture_axis_prediction"] == "multidomain_or_segmented"
    assert modular["architecture_axis_claim_allowed"] is True
    assert modular["adjudicated_truth_architecture_axis"] == (
        "multidomain_or_segmented"
    )
    assert modular["architecture_axis_same_axis_conflict"] is False

    rhodopsin = rows["pdb_1F88_A_rhodopsin"]
    assert rhodopsin["protein_regime"] == "membrane_like"
    assert rhodopsin["architecture_axis_prediction"] == "unknown"
    assert rhodopsin["architecture_axis_claim_allowed"] is False
    assert rhodopsin["architecture_axis_abstention_reason"] == (
        "membrane segmentation can reflect helices rather than global domains"
    )
    assert rhodopsin["architecture_axis_same_axis_conflict"] is False

    groel = rows["pdb_1AON_A_groel"]
    assert groel["architecture_axis_prediction"] == "unknown"
    assert groel["architecture_axis_claim_allowed"] is False
    assert groel["architecture_axis_abstention_reason"] == (
        "repeat/segmentation signals are mixed"
    )
    assert groel["adjudicated_truth_architecture_axis"] == "compact_single_domain"


def test_tracked_architecture_axis_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_50.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE_50.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_50.read_text(encoding="utf-8").splitlines()))
    abstentions = list(
        csv.DictReader(ABSTENTIONS_50.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD_50.read_text(encoding="utf-8")

    assert report["benchmark_kind"] == ARCHITECTURE_AXIS_BENCHMARK_KIND
    assert certificate["certificate_kind"] == (
        "architecture_axis_evidence_safety_certificate"
    )
    assert len(rows) == 50
    assert len(abstentions) == 34
    assert CONFLICTS_50.read_text(encoding="utf-8") == ""
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert rows[0]["architecture_axis_signature_kind"] == (
        ARCHITECTURE_AXIS_SIGNATURE_KIND
    )
    assert certificate["raw_sequences_exported"] is False
    assert certificate["global_fold_class_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Architecture Axis Evidence Adjudication" in dashboard
    assert "Architecture Is Not Secondary Structure" in dashboard
    assert "Length Is Not Enough" in dashboard
    assert "Fragments Stay Scoped" in dashboard
    assert "No Global Class Recovery" in dashboard


def test_architecture_axis_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    references, evidence = load_hierarchical_gate_inputs(
        REL_BENCHMARK_50,
        REL_STRUCTURE_EVIDENCE_50,
    )
    rows = architecture_axis_rows(references, evidence)
    report = build_architecture_axis_report(
        references,
        evidence,
        source_benchmark_file=REL_BENCHMARK_50,
        structure_evidence_file=REL_STRUCTURE_EVIDENCE_50,
    )
    outputs = {
        "report": tmp_path / REPORT_50.name,
        "rows": tmp_path / ROWS_50.name,
        "conflicts": tmp_path / CONFLICTS_50.name,
        "abstentions": tmp_path / ABSTENTIONS_50.name,
        "dashboard": tmp_path / DASHBOARD_50.name,
        "certificate": tmp_path / CERTIFICATE_50.name,
    }

    write_architecture_axis_outputs(
        report=report,
        rows=rows,
        conflicts=architecture_axis_conflict_rows(rows),
        abstentions=architecture_axis_abstention_rows(rows),
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        conflicts_path=outputs["conflicts"],
        abstentions_path=outputs["abstentions"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT_50,
        "rows": ROWS_50,
        "conflicts": CONFLICTS_50,
        "abstentions": ABSTENTIONS_50,
        "dashboard": DASHBOARD_50,
        "certificate": CERTIFICATE_50,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
