import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_axis_adjudication import (  # noqa: E402
    AXIS_SIGNATURE_KIND,
    FOLD_AXIS_ADJUDICATION_BENCHMARK_KIND,
    axis_adjudication_rows,
    build_axis_adjudication_report,
)
from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    load_hierarchical_gate_inputs,
)


BENCHMARK_50 = ROOT / "data" / "folding_benchmarks_real_50.locked.json"
STRUCTURE_EVIDENCE_50 = (
    ROOT / "data" / "folding_benchmarks_real_50_structure_evidence.json"
)
REPORT_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_adjudication_report.json"
)
ROWS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_rows.csv"
)
CONFLICTS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_conflicts.csv"
)
MANUAL_REVIEW_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_manual_review.csv"
)
CONFUSION_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_confusion_matrices.csv"
)
DASHBOARD_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_axis_dashboard.html"
)


def test_axis_adjudication_report_splits_collapsed_truth_axes() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    report = build_axis_adjudication_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_50,
        structure_evidence_file=STRUCTURE_EVIDENCE_50,
    )

    assert report["benchmark_kind"] == FOLD_AXIS_ADJUDICATION_BENCHMARK_KIND
    assert report["axis_signature_kind"] == AXIS_SIGNATURE_KIND
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_labels_no_structure_answers"
    )
    assert report["truth_adjudication_boundary"] == (
        "labels_structure_sources_and_reference_axes_used_only_after_prediction"
    )
    assert report["benchmark_size"] == 50
    assert report["single_class_taxonomy_collapse_detected"] is True
    assert report["structure_label_disagreement_count"] == 17
    assert report["orthogonal_axis_disagreement_count"] == 10
    assert report["true_same_axis_conflict_count"] == 10
    assert report["structure_label_same_axis_conflict_count"] == 7
    assert report["axis_unscorable_count"] == 144
    assert report["axis_unscorable_row_count"] == 50
    assert report["high_confidence_wrong_count_after_axis_scoring"] == 0
    assert report["forced_prediction_count"] == 16
    assert report["abstained_prediction_count"] == 34
    assert report["axis_accuracy"]["order_axis"] == {
        "accuracy": 0.9375,
        "scorable_count": 16,
        "unscorable_count": 34,
    }
    assert report["taxonomy_note_counts"]["membrane_regime_not_fold_class"] == 2
    assert report["manual_review_row_count"] == 30
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_axis_rows_separate_taxonomy_collapse_from_model_conflicts() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_50,
        STRUCTURE_EVIDENCE_50,
    )
    rows = {
        row["protein_id"]: row
        for row in axis_adjudication_rows(references, evidence)
    }

    rhodopsin = rows["pdb_1F88_A_rhodopsin"]
    assert rhodopsin["protein_regime"] == "membrane_like"
    assert rhodopsin["structure_fold_class"] == "alpha_rich"
    assert rhodopsin["label_fold_class"] == "multidomain_boundary"
    assert rhodopsin["predicted_environment_axis"] == "membrane_like"
    assert rhodopsin["adjudicated_truth_environment_axis"] == "membrane_like"
    assert rhodopsin["orthogonal_axis_disagreement"] is True
    assert rhodopsin["true_same_axis_conflict"] is False
    assert "membrane_regime_not_fold_class" in rhodopsin["manual_review_reasons"]

    basic_fgf = rows["pdb_1BFG_A_basic_fgf"]
    assert basic_fgf["predicted_order_axis"] == "disordered_flexible"
    assert basic_fgf["adjudicated_truth_order_axis"] == "ordered"
    assert basic_fgf["true_same_axis_conflict"] is True
    assert basic_fgf["axis_conflict_axes"] == "order_axis"

    barstar = rows["pdb_1BTA_A_barstar"]
    assert barstar["predicted_secondary_structure_axis"] == "alpha_rich"
    assert barstar["adjudicated_truth_secondary_structure_axis"] == (
        "alpha_beta_mixed"
    )
    assert barstar["true_same_axis_conflict"] is True
    assert barstar["axis_conflict_axes"] == "secondary_structure_axis"


def test_tracked_axis_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_50.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_50.read_text(encoding="utf-8").splitlines()))
    conflicts = list(
        csv.DictReader(CONFLICTS_50.read_text(encoding="utf-8").splitlines())
    )
    manual_review = list(
        csv.DictReader(MANUAL_REVIEW_50.read_text(encoding="utf-8").splitlines())
    )
    confusion = list(
        csv.DictReader(CONFUSION_50.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD_50.read_text(encoding="utf-8")

    assert report["benchmark_kind"] == FOLD_AXIS_ADJUDICATION_BENCHMARK_KIND
    assert len(rows) == 50
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert rows[0]["axis_signature_kind"] == AXIS_SIGNATURE_KIND
    assert len(conflicts) == 22
    assert len(manual_review) == 30
    assert len(confusion) >= 20
    assert "Single-Class Benchmark Is Lossy" in dashboard
    assert "Architecture ≠ Secondary Structure" in dashboard
    assert "Membrane Regime ≠ Fold Class" in dashboard
    assert "Disorder ≠ Beta/Alpha Absence" in dashboard
    assert "Fragment Evidence ≠ Global Fold Truth" in dashboard
