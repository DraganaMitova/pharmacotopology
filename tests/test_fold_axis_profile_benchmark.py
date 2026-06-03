import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_axis_profile import (  # noqa: E402
    AXIS_PROFILE_SIGNATURE_KIND,
    FOLD_AXIS_PROFILE_BENCHMARK_KIND,
    axis_profile_abstention_rows,
    axis_profile_recovery_candidate_rows,
    axis_profile_rows,
    build_axis_profile_report,
    write_axis_profile_outputs,
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
    / "real_folding_50_axis_profile_report.json"
)
ROWS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_profile_rows.csv"
)
ABSTENTIONS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_profile_abstentions.csv"
)
RECOVERY_CANDIDATES_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_profile_recovery_candidates.csv"
)
DASHBOARD_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_profile_dashboard.html"
)
CERTIFICATE_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_profile_certificate.json"
)


def test_axis_profile_report_recovers_axis_coverage_without_class_recovery() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    report = build_axis_profile_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_50,
        structure_evidence_file=STRUCTURE_EVIDENCE_50,
    )

    assert report["benchmark_kind"] == FOLD_AXIS_PROFILE_BENCHMARK_KIND
    assert report["axis_profile_signature_kind"] == AXIS_PROFILE_SIGNATURE_KIND
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_labels_no_structure_answers"
    )
    assert report["benchmark_size"] == 50
    assert report["collapsed_class_coverage"] == 0.28
    assert report["axis_profile_coverage"] == 0.86
    assert report["axis_profile_coverage"] > report["collapsed_class_coverage"]
    assert report["secondary_axis_coverage"] == 0.16
    assert report["architecture_axis_coverage"] == 0.14
    assert report["order_axis_coverage"] == 0.84
    assert report["environment_axis_coverage"] == 0.04
    assert report["safe_axis_recovered_count"] == 39
    assert report["safe_axis_recovered_row_count"] == 29
    assert report["safe_axis_recovered_distribution"] == {
        "architecture_axis": 7,
        "environment_axis": 2,
        "order_axis": 28,
        "secondary_structure_axis": 2,
    }
    assert report["unsafe_class_recovery_count"] == 0
    assert report["guard_override_count"] == 0
    assert report["forced_same_axis_conflict_count"] == 0
    assert report["axis_profile_same_axis_conflict_count"] == 0
    assert report["high_confidence_wrong_count_after_axis_scoring"] == 0
    assert report["global_fold_class_claim_allowed"] is False
    assert report["axis_profile_claim_allowed"] is True
    assert report["folding_problem_solved"] is False
    assert report["artifact_reproducible"] is True
    assert report["axis_accuracy"]["order_axis"] == {
        "accuracy": 1.0,
        "claimed_count": 42,
        "scorable_count": 42,
        "unscorable_count": 8,
    }
    assert report["axis_accuracy"]["secondary_structure_axis"] == {
        "accuracy": 1.0,
        "claimed_count": 8,
        "scorable_count": 8,
        "unscorable_count": 42,
    }


def test_axis_profile_rows_keep_guards_and_recover_only_axes() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    rows = {
        row["protein_id"]: row
        for row in axis_profile_rows(references, evidence)
    }

    basic_fgf = rows["pdb_1BFG_A_basic_fgf"]
    assert basic_fgf["source_abstained"] is True
    assert basic_fgf["profile_global_fold_class"] == "insufficient_topology_evidence"
    assert basic_fgf["global_class_claim_allowed"] is False
    assert basic_fgf["profile_secondary_structure_axis"] == "weak_or_unknown"
    assert basic_fgf["profile_order_axis"] == "ordered"
    assert basic_fgf["safe_axis_recovered_axes"] == "order_axis"
    assert basic_fgf["guard_override"] is False
    assert basic_fgf["axis_profile_same_axis_conflict"] is False

    barstar = rows["pdb_1BTA_A_barstar"]
    assert barstar["source_abstained"] is True
    assert barstar["profile_secondary_structure_axis"] == "weak_or_unknown"
    assert barstar["profile_order_axis"] == "ordered"
    assert barstar["safe_axis_recovered_axes"] == "order_axis"
    assert barstar["guard_override"] is False
    assert barstar["axis_profile_same_axis_conflict"] is False

    rop = rows["pdb_1ROP_A_rop_protein"]
    assert rop["source_abstained"] is True
    assert rop["profile_secondary_structure_axis"] == "alpha_rich"
    assert rop["profile_architecture_axis"] == "fragment_scope"
    assert rop["profile_order_axis"] == "ordered"
    assert rop["safe_axis_recovered_axes"] == (
        "secondary_structure_axis;architecture_axis;order_axis"
    )
    assert rop["global_class_claim_allowed"] is False

    rhodopsin = rows["pdb_1F88_A_rhodopsin"]
    assert rhodopsin["profile_secondary_structure_axis"] == "weak_or_unknown"
    assert rhodopsin["profile_architecture_axis"] == "unknown"
    assert rhodopsin["profile_order_axis"] == "ordered"
    assert rhodopsin["profile_environment_axis"] == "membrane_like"
    assert rhodopsin["safe_axis_recovered_axes"] == "order_axis;environment_axis"
    assert rhodopsin["axis_profile_same_axis_conflict"] is False

    stathmin = rows["disprot_P16949_stathmin"]
    assert stathmin["source_abstained"] is True
    assert stathmin["profile_order_axis"] == "mixed_or_uncertain"
    assert stathmin["safe_axis_recovered_axes"] == ""
    assert stathmin["guard_override"] is False


def test_tracked_axis_profile_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_50.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE_50.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_50.read_text(encoding="utf-8").splitlines()))
    abstentions = list(
        csv.DictReader(ABSTENTIONS_50.read_text(encoding="utf-8").splitlines())
    )
    recovery_candidates = list(
        csv.DictReader(
            RECOVERY_CANDIDATES_50.read_text(encoding="utf-8").splitlines()
        )
    )
    dashboard = DASHBOARD_50.read_text(encoding="utf-8")

    assert report["benchmark_kind"] == FOLD_AXIS_PROFILE_BENCHMARK_KIND
    assert certificate["certificate_kind"] == "fold_axis_profile_safety_certificate"
    assert len(rows) == 50
    assert len(abstentions) == 36
    assert len(recovery_candidates) == 29
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert rows[0]["axis_profile_signature_kind"] == AXIS_PROFILE_SIGNATURE_KIND
    assert certificate["raw_sequences_exported"] is False
    assert certificate["global_fold_class_claim_allowed"] is False
    assert certificate["axis_profile_claim_allowed"] is True
    assert "Fold Axis Profile Coverage Recovery" in dashboard
    assert "No Collapsed Class Recovery" in dashboard
    assert "I Know This Axis" in dashboard
    assert "I Do Not Know That Axis" in dashboard
    assert "No Guard Overrides" in dashboard


def test_axis_profile_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    references, evidence = load_hierarchical_gate_inputs(
        REL_BENCHMARK_50,
        REL_STRUCTURE_EVIDENCE_50,
    )
    rows = axis_profile_rows(references, evidence)
    report = build_axis_profile_report(
        references,
        evidence,
        source_benchmark_file=REL_BENCHMARK_50,
        structure_evidence_file=REL_STRUCTURE_EVIDENCE_50,
    )
    outputs = {
        "report": tmp_path / REPORT_50.name,
        "rows": tmp_path / ROWS_50.name,
        "abstentions": tmp_path / ABSTENTIONS_50.name,
        "recovery": tmp_path / RECOVERY_CANDIDATES_50.name,
        "dashboard": tmp_path / DASHBOARD_50.name,
        "certificate": tmp_path / CERTIFICATE_50.name,
    }

    write_axis_profile_outputs(
        report=report,
        rows=rows,
        abstention_rows=axis_profile_abstention_rows(rows),
        recovery_candidate_rows=axis_profile_recovery_candidate_rows(rows),
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        abstentions_path=outputs["abstentions"],
        recovery_candidates_path=outputs["recovery"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT_50,
        "rows": ROWS_50,
        "abstentions": ABSTENTIONS_50,
        "recovery": RECOVERY_CANDIDATES_50,
        "dashboard": DASHBOARD_50,
        "certificate": CERTIFICATE_50,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
