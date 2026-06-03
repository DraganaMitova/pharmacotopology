import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_hierarchical_gates import (  # noqa: E402
    HIERARCHICAL_GATE_BENCHMARK_KIND,
    HIERARCHICAL_GATE_SIGNATURE_KIND,
    build_hierarchical_gate_report,
    gate_failure_rows,
    gate_path_rows,
    hierarchical_gate_rows,
    load_hierarchical_gate_inputs,
    predict_hierarchical_gate,
)
from pharmacotopology.folding_motif_alignment import (  # noqa: E402
    ABSTAINED_CLASS,
    MOTIF_SIGNATURE_KIND,
)


BENCHMARK_FILE = ROOT / "data" / "folding_benchmarks_real_10.locked.json"
STRUCTURE_EVIDENCE_FILE = (
    ROOT / "data" / "folding_benchmarks_real_10_structure_evidence.json"
)
LEGACY_ARTIFACT_DIR = (
    ROOT / "first_contact_clean_pharmacotopology_layer_run" / "archived_legacy"
)
GATE_REPORT = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_hierarchical_gate_report.json"
)
GATE_ROWS = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_hierarchical_gate_rows.csv"
)
GATE_PATHS = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_gate_paths.csv"
)
GATE_FAILURES = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_gate_failures.csv"
)


def test_hierarchical_prediction_uses_staged_gate_path() -> None:
    references, _ = load_hierarchical_gate_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    prediction = predict_hierarchical_gate(
        references[0].sequence,
        protein_id=references[0].protein_id,
    )

    assert prediction.hierarchy_prediction_used is True
    assert prediction.beta_evidence_requires_pairing is True
    assert prediction.alpha_evidence_requires_periodicity is True
    assert prediction.gate_path[0].startswith("disorder_gate:")
    assert prediction.gate_path[1].startswith("compactness_gate:")
    assert prediction.gate_path[2].startswith("segmentation_gate:")
    assert prediction.gate_path[3].startswith("secondary_structure_gate:")
    assert prediction.predicted_fold_class in {
        ABSTAINED_CLASS,
        "alpha_rich",
        "beta_rich",
        "alpha_beta_mixed",
        "multidomain_boundary",
        "disordered_flexible",
    }


def test_hierarchical_gate_report_reduces_specific_failure_modes() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    report = build_hierarchical_gate_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_FILE,
        structure_evidence_file=STRUCTURE_EVIDENCE_FILE,
    )

    assert report["benchmark_kind"] == HIERARCHICAL_GATE_BENCHMARK_KIND
    assert report["topology_evidence_vector_kind"] == MOTIF_SIGNATURE_KIND
    assert report["hierarchical_gate_signature_kind"] == (
        HIERARCHICAL_GATE_SIGNATURE_KIND
    )
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_labels_no_structure_answers"
    )
    assert report["truth_channels_used_only_after_prediction"] is True
    assert report["prediction_vs_structure_accuracy"] == 0.6
    assert report["prediction_vs_label_accuracy"] == 0.6
    assert report["forced_prediction_count"] == 6
    assert report["abstained_prediction_count"] == 4
    assert report["high_confidence_wrong_count"] == 0
    assert report["evidence_conflict_mean"] == 0.549305
    assert report["disorder_gate_accuracy"] == 1.0
    assert report["compactness_gate_accuracy"] == 1.0
    assert report["segmentation_gate_accuracy"] == 1.0
    assert report["secondary_structure_gate_accuracy"] == 0.428571
    assert report["flexible_segmentation_false_multidomain_count"] == 0
    assert report["false_beta_from_disorder_count"] == 0
    assert report["false_mixed_from_alpha_count"] == 0
    assert report["hierarchy_changed_raw_prediction_count"] == 9
    assert report["hierarchy_prevented_false_multidomain_count"] == 2
    assert report["hierarchy_prevented_false_beta_count"] == 1
    assert report["hierarchy_prevented_false_mixed_count"] == 2
    assert report["revision_required"] is True
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_hierarchical_rows_expose_required_gate_fields() -> None:
    references, evidence = load_hierarchical_gate_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    rows = hierarchical_gate_rows(references, evidence)
    paths = gate_path_rows(rows)
    failures = gate_failure_rows(rows)

    assert list(rows[0])[:16] == [
        "protein_id",
        "sequence_length",
        "topology_evidence_vector_kind",
        "hierarchical_gate_signature_kind",
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
    assert len(paths) == 10
    assert len(failures) == 4
    assert any(row["flexible_segmentation_warning"] for row in rows)
    assert any(row["hierarchy_prevented_false_beta"] for row in rows)


def test_tracked_hierarchical_gate_outputs_have_expected_fields() -> None:
    report = json.loads(GATE_REPORT.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(GATE_ROWS.read_text(encoding="utf-8").splitlines()))
    path_rows = list(
        csv.DictReader(GATE_PATHS.read_text(encoding="utf-8").splitlines())
    )
    failure_rows = list(
        csv.DictReader(GATE_FAILURES.read_text(encoding="utf-8").splitlines())
    )

    assert report["benchmark_kind"] == HIERARCHICAL_GATE_BENCHMARK_KIND
    assert report["high_confidence_wrong_count"] == 0
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert len(rows) == 10
    assert rows[0]["hierarchical_gate_signature_kind"] == (
        HIERARCHICAL_GATE_SIGNATURE_KIND
    )
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert len(path_rows) == 10
    assert "gate_path" in path_rows[0]
    assert len(failure_rows) == 4
    assert "gate_decision_reason" in failure_rows[0]
