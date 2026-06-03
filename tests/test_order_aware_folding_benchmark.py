import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_order_aware import (  # noqa: E402
    ORDER_AWARE_BENCHMARK_KIND,
    build_order_aware_report,
    contact_prior_rows,
    load_order_aware_inputs,
    order_aware_control_separation_rows,
    predict_order_aware_topology_signature,
)
from pharmacotopology.folding_structure_benchmark import (  # noqa: E402
    _control_sequence,
    _signature_delta_mean,
)
from pharmacotopology.folding_topology import normalize_sequence  # noqa: E402


BENCHMARK_FILE = ROOT / "data" / "folding_benchmarks_real_10.locked.json"
STRUCTURE_EVIDENCE_FILE = (
    ROOT / "data" / "folding_benchmarks_real_10_structure_evidence.json"
)
LEGACY_ARTIFACT_DIR = (
    ROOT / "first_contact_clean_pharmacotopology_layer_run" / "archived_legacy"
)
ORDER_AWARE_REPORT = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_order_aware_report.json"
)
CONTROL_SEPARATION = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_control_separation.csv"
)
CONTACT_PRIOR = (
    LEGACY_ARTIFACT_DIR
    / "real_folding_10_contact_prior.csv"
)


def test_order_aware_prediction_changes_under_order_controls() -> None:
    references, _ = load_order_aware_inputs(BENCHMARK_FILE, STRUCTURE_EVIDENCE_FILE)
    reference = references[0]
    sequence = normalize_sequence(reference.sequence)
    real = predict_order_aware_topology_signature(
        sequence,
        protein_id=reference.protein_id,
    )
    control_sequence = _control_sequence(
        sequence,
        "hydrophobic_cluster_destroyed",
        protein_id=reference.protein_id,
    )
    control = predict_order_aware_topology_signature(
        control_sequence,
        protein_id=f"{reference.protein_id}:control",
    )

    assert sorted(sequence) == sorted(control_sequence)
    assert _signature_delta_mean(real.topology_signature, control.topology_signature) > 0
    assert real.order_features["hydrophobic_interruption_rate"] != (
        control.order_features["hydrophobic_interruption_rate"]
    )
    assert real.contact_prior_edges


def test_order_aware_report_exposes_nonzero_order_signal() -> None:
    references, evidence = load_order_aware_inputs(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    report = build_order_aware_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_FILE,
        structure_evidence_file=STRUCTURE_EVIDENCE_FILE,
    )
    controls = order_aware_control_separation_rows(references)
    contacts = contact_prior_rows(references)

    assert report["benchmark_kind"] == ORDER_AWARE_BENCHMARK_KIND
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_labels_no_structure_answers"
    )
    assert report["sequence_order_sensitivity_score"] == 0.278079
    assert report["real_vs_shuffled_separation_mean"] == 0.278079
    assert report["contact_prior_signal_seen"] is True
    assert report["recipe_order_blind"] is False
    assert report["composition_only_warning"] is False
    assert report["prediction_vs_structure_accuracy"] == 0.2
    assert report["prediction_vs_label_accuracy"] == 0.1
    assert report["revision_required"] is True
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert len(controls) == 50
    assert contacts
    assert all(row["control_sequence_written"] is False for row in controls)
    assert all("control_sequence" not in row for row in controls)


def test_tracked_order_aware_outputs_have_expected_boundary_fields() -> None:
    report = json.loads(ORDER_AWARE_REPORT.read_text(encoding="utf-8"))
    control_rows = list(
        csv.DictReader(CONTROL_SEPARATION.read_text(encoding="utf-8").splitlines())
    )
    contact_rows = list(
        csv.DictReader(CONTACT_PRIOR.read_text(encoding="utf-8").splitlines())
    )

    assert report["benchmark_kind"] == ORDER_AWARE_BENCHMARK_KIND
    assert report["sequence_order_sensitivity_score"] == 0.278079
    assert report["composition_only_warning"] is False
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert len(control_rows) == 50
    assert "control_sequence" not in control_rows[0]
    assert {row["same_composition"] for row in control_rows} == {"True"}
    assert contact_rows
    assert "edge_weight" in contact_rows[0]
    assert "label_fold_class" not in contact_rows[0]
