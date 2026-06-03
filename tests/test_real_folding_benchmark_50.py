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
)
from pharmacotopology.folding_reference_loader import (  # noqa: E402
    load_folding_reference_dataset,
)
from pharmacotopology.folding_structure_benchmark import (  # noqa: E402
    load_structure_evidence_rows,
)


BENCHMARK_50 = ROOT / "data" / "folding_benchmarks_real_50.locked.json"
STRUCTURE_EVIDENCE_50 = (
    ROOT / "data" / "folding_benchmarks_real_50_structure_evidence.json"
)
REPORT_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_hierarchical_gate_report.json"
)
ROWS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_hierarchical_gate_rows.csv"
)
PATHS_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_gate_paths.csv"
)
FAILURES_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_gate_failures.csv"
)
CERTIFICATE_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_certificate.json"
)
CONFUSION_50 = (
    ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "real_folding_50_confusion_matrix.csv"
)


def test_locked_real_50_dataset_is_stratified_external_data() -> None:
    dataset = load_folding_reference_dataset(BENCHMARK_50, require_external=True)
    payload = json.loads(BENCHMARK_50.read_text(encoding="utf-8"))

    assert dataset.validation.valid is True
    assert dataset.validation.references_loaded == 50
    assert dataset.validation.external_reference_count == 50
    assert payload["locked_after_generation"] is True
    assert payload["no_retuning_flag"] is True
    assert payload["stability_test_of_commit"] == "17dc6e2"
    assert payload["class_distribution"] == {
        "alpha_beta_mixed": 10,
        "alpha_rich": 10,
        "beta_rich": 10,
        "disordered_flexible": 10,
        "multidomain_boundary": 10,
    }
    assert payload["folding_problem_solved"] is False


def test_real_50_structure_evidence_keeps_coordinate_and_disorder_channels() -> None:
    rows = load_structure_evidence_rows(STRUCTURE_EVIDENCE_50)

    assert len(rows) == 50
    assert sum(1 for row in rows if row.evidence_kind == "coordinate_contact_graph") == 40
    assert sum(1 for row in rows if row.evidence_kind == "disorder_reference") == 10
    assert all(row.folding_problem_solved is False for row in rows)
    assert all(row.folding_solution_claim_created is False for row in rows)


def test_real_50_hierarchical_report_records_stability_result() -> None:
    report = json.loads(REPORT_50.read_text(encoding="utf-8"))

    assert report["benchmark_kind"] == HIERARCHICAL_GATE_BENCHMARK_KIND
    assert report["benchmark_size"] == 50
    assert report["external_rows"] == 50
    assert report["class_distribution"] == {
        "alpha_beta_mixed": 10,
        "alpha_rich": 10,
        "beta_rich": 10,
        "disordered_flexible": 10,
        "multidomain_boundary": 10,
    }
    assert report["prediction_vs_structure_accuracy"] == 0.34
    assert report["prediction_vs_label_accuracy"] == 0.38
    assert report["forced_prediction_count"] == 27
    assert report["abstained_prediction_count"] == 23
    assert report["high_confidence_wrong_count"] == 6
    assert report["false_beta_from_disorder_count"] == 0
    assert report["false_mixed_from_alpha_count"] == 0
    assert report["flexible_segmentation_false_multidomain_count"] == 0
    assert report["disorder_gate_accuracy"] == 0.9
    assert report["compactness_gate_accuracy"] == 0.82
    assert report["segmentation_gate_accuracy"] == 0.9
    assert report["secondary_structure_gate_accuracy"] == 0.194444
    assert report["hierarchy_changed_raw_prediction_count"] == 38
    assert report["accuracy_delta_from_10"] == -0.26
    assert report["abstention_delta_from_10"] == 19
    assert report["high_confidence_wrong_delta_from_10"] == 6
    assert report["stability_status"] == "unstable_accuracy_drop"
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["revision_required"] is True


def test_real_50_tracked_outputs_have_required_surfaces() -> None:
    rows = list(csv.DictReader(ROWS_50.read_text(encoding="utf-8").splitlines()))
    path_rows = list(csv.DictReader(PATHS_50.read_text(encoding="utf-8").splitlines()))
    failure_rows = list(
        csv.DictReader(FAILURES_50.read_text(encoding="utf-8").splitlines())
    )
    certificate = json.loads(CERTIFICATE_50.read_text(encoding="utf-8"))
    confusion_rows = list(
        csv.DictReader(CONFUSION_50.read_text(encoding="utf-8").splitlines())
    )

    assert len(rows) == 50
    assert "sequence" not in rows[0]
    assert "control_sequence" not in rows[0]
    assert len(path_rows) == 50
    assert "gate_path" in path_rows[0]
    assert len(failure_rows) == 33
    assert certificate["benchmark_size"] == 50
    assert certificate["external_rows"] == 50
    assert certificate["folding_problem_solved"] is False
    assert len(confusion_rows) == 5
    assert "insufficient_topology_evidence" in confusion_rows[0]
