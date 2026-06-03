import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_regime_analysis import (  # noqa: E402
    REGIME_ANALYSIS_BENCHMARK_KIND,
    REGIME_SIGNATURE_KIND,
    build_regime_analysis_report,
    detect_protein_regime,
    load_hierarchical_gate_inputs,
    regime_analysis_rows,
)


BENCHMARK_50 = ROOT / "data" / "folding_benchmarks_real_50.locked.json"
STRUCTURE_EVIDENCE_50 = (
    ROOT / "data" / "folding_benchmarks_real_50_structure_evidence.json"
)
LEGACY_ARTIFACT_DIR = (
    ROOT / "first_contact_clean_pharmacotopology_layer_run" / "archived_legacy"
)
REPORT_50 = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_50_regime_analysis_report.json"
)
ROWS_50 = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_50_regime_rows.csv"
)
FAILURE_COHORTS_50 = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_50_failure_cohorts.csv"
)
HIGH_CONFIDENCE_WRONG_50 = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_50_high_confidence_wrong.csv"
)
ABSTENTION_ANALYSIS_50 = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_50_abstention_analysis.csv"
)
DASHBOARD_50 = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_50_regime_dashboard.html"
)


def test_regime_detection_uses_sequence_only_features() -> None:
    references, _evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    rhodopsin = next(
        reference
        for reference in references
        if reference.protein_id == "pdb_1F88_A_rhodopsin"
    )
    nupr1 = next(
        reference
        for reference in references
        if reference.protein_id == "disprot_O60356_nupr1"
    )

    rhodopsin_regime = detect_protein_regime(
        rhodopsin.sequence,
        protein_id=rhodopsin.protein_id,
    )
    nupr1_regime = detect_protein_regime(
        nupr1.sequence,
        protein_id=nupr1.protein_id,
    )

    assert rhodopsin_regime.regime_detection_used is True
    assert rhodopsin_regime.protein_regime == "membrane_like"
    assert rhodopsin_regime.regime_allowed_gate_path == (
        "abstain_until_membrane_gate_exists"
    )
    assert nupr1_regime.protein_regime == "intrinsically_disordered"
    assert nupr1_regime.regime_allowed_gate_path == "disorder_gate_or_abstain"


def test_regime_analysis_report_records_required_diagnostic_fields() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    report = build_regime_analysis_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_50,
        structure_evidence_file=STRUCTURE_EVIDENCE_50,
    )

    required_fields = {
        "benchmark_size",
        "prediction_vs_structure_accuracy",
        "prediction_vs_label_accuracy",
        "forced_prediction_count",
        "abstained_prediction_count",
        "high_confidence_wrong_count",
        "regime_detection_used",
        "regime_accuracy",
        "regime_confidence_mean",
        "accuracy_by_regime",
        "abstention_by_regime",
        "high_confidence_wrong_by_regime",
        "dominant_failure_cohort",
        "new_failure_modes_detected",
        "old_failure_modes_still_closed",
        "structure_label_disagreement_count",
        "ambiguous_reference_count",
        "possible_bad_rows_count",
        "rows_requiring_manual_review",
        "revision_required",
        "claim_allowed",
        "folding_problem_solved",
    }

    assert required_fields.issubset(report)
    assert report["benchmark_kind"] == REGIME_ANALYSIS_BENCHMARK_KIND
    assert report["protein_regime_signature_kind"] == REGIME_SIGNATURE_KIND
    assert report["regime_detection_boundary"] == (
        "sequence_only_no_labels_no_pdb_no_cath_no_disprot_truth"
    )
    assert report["benchmark_size"] == 50
    assert report["prediction_vs_structure_accuracy"] == 0.28
    assert report["prediction_vs_label_accuracy"] == 0.28
    assert report["forced_prediction_count"] == 14
    assert report["abstained_prediction_count"] == 36
    assert report["high_confidence_wrong_count"] == 0
    assert report["regime_detection_used"] is True
    assert report["regime_accuracy"] == 0.74
    assert report["regime_confidence_mean"] == 0.738114
    assert report["structure_label_disagreement_count"] == 17
    assert report["ambiguous_reference_count"] == 22
    assert report["possible_bad_rows_count"] == 8
    assert report["manual_review_row_count"] == 23
    assert report["hierarchical_high_confidence_wrong_count_before_regime_routing"] == 6
    assert report["hierarchical_high_confidence_wrong_prevented_by_regime_routing"] == 6
    assert report["dominant_failure_cohort"]["failure_cohort"] == (
        "abstained_unresolved"
    )
    assert report["dominant_failure_cohort"]["count"] == 12
    assert report["failure_cohort_distribution"] == {
        "abstained_on_regime_risk": 6,
        "abstained_on_structure_label_disagreement": 12,
        "abstained_unresolved": 12,
        "correct_forced": 14,
        "regime_prevented_high_confidence_wrong": 6,
    }
    assert report["old_failure_modes_still_closed"] == {
        "false_beta_from_disorder_count": {"count": 0, "status": "closed"},
        "false_mixed_from_alpha_count": {"count": 0, "status": "closed"},
        "flexible_segmentation_false_multidomain_count": {
            "count": 0,
            "status": "closed",
        },
    }
    assert "structure_label_disagreement_requires_manual_review" in report[
        "new_failure_modes_detected"
    ]
    assert report["revision_required"] is True
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_regime_guard_abstains_folded_domain_mimic_disorder_case() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    rows = {
        row["protein_id"]: row
        for row in regime_analysis_rows(references, evidence)
    }
    basic_fgf = rows["pdb_1BFG_A_basic_fgf"]

    assert basic_fgf["protein_regime"] == "intrinsically_disordered"
    assert basic_fgf["predicted_fold_class"] == "insufficient_topology_evidence"
    assert basic_fgf["forced_prediction"] is False
    assert basic_fgf["abstained"] is True
    assert "abstained_folded_domain_mimic_disorder_conflict" in basic_fgf[
        "gate_path"
    ]


def test_tracked_regime_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_50.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_50.read_text(encoding="utf-8").splitlines()))
    cohort_rows = list(
        csv.DictReader(FAILURE_COHORTS_50.read_text(encoding="utf-8").splitlines())
    )
    high_confidence_wrong_rows = list(
        csv.DictReader(
            HIGH_CONFIDENCE_WRONG_50.read_text(encoding="utf-8").splitlines()
        )
    )
    abstention_rows = list(
        csv.DictReader(
            ABSTENTION_ANALYSIS_50.read_text(encoding="utf-8").splitlines()
        )
    )
    dashboard = DASHBOARD_50.read_text(encoding="utf-8")

    assert report["benchmark_kind"] == REGIME_ANALYSIS_BENCHMARK_KIND
    assert len(rows) == 50
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert rows[0]["protein_regime_signature_kind"] == REGIME_SIGNATURE_KIND
    assert any(row["protein_regime"] == "membrane_like" for row in rows)
    assert any(row["protein_regime"] == "intrinsically_disordered" for row in rows)
    assert cohort_rows[0]["failure_cohort"]
    assert len(high_confidence_wrong_rows) == 0
    assert len(abstention_rows) == 36
    assert "Failure Cohort Heatmap" in dashboard
    assert "Old Failure Counters Stayed Closed" in dashboard
