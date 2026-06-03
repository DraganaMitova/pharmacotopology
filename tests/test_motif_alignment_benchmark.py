import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_motif_alignment import (  # noqa: E402
    ABSTAINED_CLASS,
    MOTIF_ALIGNMENT_BENCHMARK_KIND,
    MOTIF_SIGNATURE_KIND,
    build_motif_alignment_report,
    evidence_conflict_rows,
    failure_diagnosis_rows,
    load_motif_alignment_inputs,
    motif_alignment_rows,
    predict_motif_alignment,
)


BENCHMARK_FILE = ROOT / "data" / "folding_benchmarks_real_10.locked.json"
STRUCTURE_EVIDENCE_FILE = (
    ROOT / "data" / "folding_benchmarks_real_10_structure_evidence.json"
)
MOTIF_REPORT = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_10_motif_alignment_report.json"
)
MOTIF_ROWS = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_10_motif_alignment_rows.csv"
)
FAILURE_DIAGNOSIS = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_10_failure_diagnosis.csv"
)
EVIDENCE_CONFLICTS = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_10_evidence_conflicts.csv"
)


def test_motif_prediction_returns_evidence_before_class_claim() -> None:
    references, _ = load_motif_alignment_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    prediction = predict_motif_alignment(
        references[0].sequence,
        protein_id=references[0].protein_id,
    )
    evidence_keys = tuple(prediction.motif_evidence.to_dict())

    assert evidence_keys[:8] == (
        "alpha_periodicity_evidence",
        "beta_alternation_evidence",
        "compact_core_evidence",
        "disorder_run_evidence",
        "domain_boundary_evidence",
        "long_range_closure_evidence",
        "breaker_turn_evidence",
        "charge_frustration_evidence",
    )
    assert prediction.motif_signal_seen is True
    assert prediction.predicted_fold_class in {
        ABSTAINED_CLASS,
        "alpha_rich",
        "beta_rich",
        "alpha_beta_mixed",
        "multidomain_boundary",
        "disordered_flexible",
    }


def test_motif_alignment_report_uses_uncertainty_gate() -> None:
    references, evidence = load_motif_alignment_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    report = build_motif_alignment_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_FILE,
        structure_evidence_file=STRUCTURE_EVIDENCE_FILE,
    )

    assert report["benchmark_kind"] == MOTIF_ALIGNMENT_BENCHMARK_KIND
    assert report["topology_evidence_vector_kind"] == MOTIF_SIGNATURE_KIND
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_labels_no_structure_answers"
    )
    assert report["truth_channels_used_only_after_prediction"] is True
    assert report["prediction_vs_structure_accuracy"] == 0.1
    assert report["prediction_vs_label_accuracy"] == 0.1
    assert report["sequence_order_sensitivity_score"] == 0.278079
    assert report["real_vs_shuffled_separation_mean"] == 0.278079
    assert report["contact_prior_signal_seen"] is True
    assert report["motif_signal_seen"] is True
    assert report["evidence_conflict_mean"] == 0.895054
    assert report["uncertainty_gating_used"] is True
    assert report["forced_prediction_count"] == 2
    assert report["abstained_prediction_count"] == 8
    assert report["high_confidence_wrong_count"] == 0
    assert report["revision_required"] is True
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_motif_rows_keep_evidence_columns_first() -> None:
    references, evidence = load_motif_alignment_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    rows = motif_alignment_rows(references, evidence)

    assert list(rows[0])[:16] == [
        "protein_id",
        "sequence_length",
        "topology_evidence_vector_kind",
        "motif_signal_seen",
        "alpha_periodicity_evidence",
        "beta_alternation_evidence",
        "compact_core_evidence",
        "disorder_run_evidence",
        "domain_boundary_evidence",
        "long_range_closure_evidence",
        "breaker_turn_evidence",
        "charge_frustration_evidence",
        "local_alpha_pressure_evidence",
        "local_beta_pressure_evidence",
        "local_disorder_pressure_evidence",
        "mixed_motif_evidence",
    ]
    assert all("sequence" not in row for row in rows)
    assert all("control_sequence" not in row for row in rows)
    assert failure_diagnosis_rows(rows)
    assert evidence_conflict_rows(rows)


def test_tracked_motif_alignment_outputs_have_diagnosis_fields() -> None:
    report = json.loads(MOTIF_REPORT.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(MOTIF_ROWS.read_text(encoding="utf-8").splitlines()))
    failure_rows = list(
        csv.DictReader(FAILURE_DIAGNOSIS.read_text(encoding="utf-8").splitlines())
    )
    conflict_rows = list(
        csv.DictReader(EVIDENCE_CONFLICTS.read_text(encoding="utf-8").splitlines())
    )

    assert report["benchmark_kind"] == MOTIF_ALIGNMENT_BENCHMARK_KIND
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["high_confidence_wrong_count"] == 0
    assert len(rows) == 10
    assert rows[0]["topology_evidence_vector_kind"] == MOTIF_SIGNATURE_KIND
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert len(failure_rows) == 9
    assert "likely_failure_reason" in failure_rows[0]
    assert "confidence" in failure_rows[0]
    assert len(conflict_rows) == 10
    assert "evidence_conflict_score" in conflict_rows[0]
    assert "uncertainty_radius" in conflict_rows[0]
    assert "claim_strength" in conflict_rows[0]
