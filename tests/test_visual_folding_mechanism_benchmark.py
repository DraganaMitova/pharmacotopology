import csv
import json
from dataclasses import replace
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_contact_topology import (  # noqa: E402
    CONTACT_TOPOLOGY_SIGNATURE_KIND,
    PREDICTOR_INPUT_BOUNDARY,
    VISUAL_MECHANISM_BENCHMARK_KIND,
    load_visual_mechanism_rows,
    predict_contact_topology,
    validate_visual_mechanism_lock,
)
from pharmacotopology.folding_visual_mechanism_benchmark import (  # noqa: E402
    PER_ROW_VISUAL_NAMES,
    ROOT_OUTPUT_NAMES,
    VISUAL_MECHANISM_CERTIFICATE_KIND,
    VISUAL_MECHANISM_SIGNATURE_KIND,
    build_visual_mechanism_report,
    contact_metric_rows,
    failure_cohort_rows,
    safe_visual_mechanism_rows,
    visual_mechanism_packets,
    write_visual_mechanism_outputs,
)


BENCHMARK_12 = ROOT / "data" / "folding_mechanism_visual_12.locked.json"
REL_BENCHMARK_12 = Path("data/folding_mechanism_visual_12.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT_12 = RUN_DIR / "visual_mechanism_12_report.json"
ROWS_12 = RUN_DIR / "visual_mechanism_12_rows.csv"
CONTACT_METRICS_12 = RUN_DIR / "visual_mechanism_12_contact_metrics.csv"
FAILURE_COHORTS_12 = RUN_DIR / "visual_mechanism_12_failure_cohorts.csv"
DASHBOARD_12 = RUN_DIR / "visual_mechanism_12_dashboard.html"
CERTIFICATE_12 = RUN_DIR / "visual_mechanism_12_certificate.json"
VISUALS_ROOT = RUN_DIR / "visuals"


def _report_from_locked_data(source_path: Path = REL_BENCHMARK_12):
    rows = load_visual_mechanism_rows(source_path)
    packets = visual_mechanism_packets(rows)
    report = build_visual_mechanism_report(
        packets=packets,
        source_benchmark_file=source_path,
        lock_validation=validate_visual_mechanism_lock(rows),
    )
    return rows, packets, report


def test_visual_mechanism_report_meets_safety_and_visibility_targets() -> None:
    rows, packets, report = _report_from_locked_data()

    assert len(rows) == 12
    assert report["benchmark_kind"] == VISUAL_MECHANISM_BENCHMARK_KIND
    assert report["visual_mechanism_signature_kind"] == (
        VISUAL_MECHANISM_SIGNATURE_KIND
    )
    assert report["predictor_input_boundary"] == (
        "sequence_only_no_native_contacts_no_truth_axes"
    )
    assert report["benchmark_size"] == 12
    assert report["visual_artifacts_generated_for_rows"] == 12
    assert report["visual_files_per_row"] == 7
    assert report["visual_artifacts_generated_count"] == 90
    assert report["contact_map_f1_computed_count"] == 12
    assert report["visible_partial_success_count"] == 5
    assert report["visible_failure_count"] == 7
    assert report["failures_visualized"] is True
    assert report["mean_contact_map_f1"] == 0.104031
    assert report["max_contact_map_f1"] == 0.363636
    assert report["failure_cohorts"] == {
        "beta_long_range_pairing_failure": 3,
        "disorder_control_overclosure_failure": 1,
        "false_contact_overprediction": 2,
        "membrane_contact_topology_failure": 1,
        "visible_partial_success": 5,
    }
    assert report["native_truth_used_before_prediction"] is False
    assert report["raw_sequence_exposed"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False
    assert all(not packet.prediction.raw_sequence_exposed for packet in packets)
    assert all(
        not packet.prediction.native_truth_used_before_prediction
        for packet in packets
    )


def test_contact_prediction_is_sequence_only_before_native_scoring() -> None:
    row = load_visual_mechanism_rows(BENCHMARK_12)[0]
    original = predict_contact_topology(row.sequence, row_id=row.row_id)
    changed_truth = replace(
        row,
        native_contact_pairs=((1, row.length),),
        native_contact_map_hash="changed-after-prediction",
        truth_axes={
            "secondary_structure_axis": "beta_rich",
            "architecture_axis": "multidomain_or_segmented",
            "order_axis": "disordered_flexible",
            "environment_axis": "membrane_like",
        },
    )
    repeated = predict_contact_topology(
        changed_truth.sequence,
        row_id=changed_truth.row_id,
    )

    assert original.contact_topology_signature_kind == CONTACT_TOPOLOGY_SIGNATURE_KIND
    assert original.predictor_input_boundary == PREDICTOR_INPUT_BOUNDARY
    assert original.predicted_contact_map_hash == repeated.predicted_contact_map_hash
    assert original.predicted_contact_pairs == repeated.predicted_contact_pairs
    assert original.native_truth_used_before_prediction is False
    assert repeated.native_truth_used_before_prediction is False


def test_checked_in_visual_mechanism_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_12.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE_12.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_12.read_text(encoding="utf-8").splitlines()))
    contact_metrics = list(
        csv.DictReader(CONTACT_METRICS_12.read_text(encoding="utf-8").splitlines())
    )
    failure_cohorts = list(
        csv.DictReader(FAILURE_COHORTS_12.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD_12.read_text(encoding="utf-8")
    visual_files = sorted(path for path in VISUALS_ROOT.glob("*/*") if path.is_file())

    assert report["benchmark_kind"] == VISUAL_MECHANISM_BENCHMARK_KIND
    assert certificate["certificate_kind"] == VISUAL_MECHANISM_CERTIFICATE_KIND
    assert len(rows) == 12
    assert len(contact_metrics) == 12
    assert len(failure_cohorts) == 5
    assert len(visual_files) == 84
    assert "sequence" not in rows[0]
    assert "raw_sequence" not in rows[0]
    assert certificate["raw_sequence_exposed"] is False
    assert certificate["native_truth_used_before_prediction"] is False
    assert certificate["global_folding_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "This Is A Mechanism Visualization, Not A Solved Folding Engine" in dashboard
    assert "Contact Maps Are The First Proof Target" in dashboard
    assert "Trajectory Is Coarse-Grained" in dashboard
    assert "Native Contacts Are Used Only After Prediction" in dashboard
    assert "Failures Are Visual Evidence" in dashboard
    assert "Global Folding Claim Remains Locked" in dashboard
    for row in rows:
        row_dir = VISUALS_ROOT / row["row_id"]
        for name in PER_ROW_VISUAL_NAMES:
            assert (row_dir / name).exists()


def test_visual_mechanism_outputs_do_not_export_raw_sequences() -> None:
    locked_rows = load_visual_mechanism_rows(BENCHMARK_12)
    generated_paths = list(RUN_DIR.glob("visual_mechanism_12_*"))
    generated_paths.extend(path for path in VISUALS_ROOT.glob("*/*") if path.is_file())

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in locked_rows:
            assert row.sequence not in text


def test_visual_mechanism_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    rows, packets, report = _report_from_locked_data()
    outputs = {
        "report": tmp_path / REPORT_12.name,
        "rows": tmp_path / ROWS_12.name,
        "contact_metrics": tmp_path / CONTACT_METRICS_12.name,
        "failure_cohorts": tmp_path / FAILURE_COHORTS_12.name,
        "dashboard": tmp_path / DASHBOARD_12.name,
        "certificate": tmp_path / CERTIFICATE_12.name,
    }
    write_visual_mechanism_outputs(
        report=report,
        packets=packets,
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        contact_metrics_path=outputs["contact_metrics"],
        failure_cohorts_path=outputs["failure_cohorts"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
        visuals_root=tmp_path / "visuals",
    )

    checked_in = {
        "report": REPORT_12,
        "rows": ROWS_12,
        "contact_metrics": CONTACT_METRICS_12,
        "failure_cohorts": FAILURE_COHORTS_12,
        "dashboard": DASHBOARD_12,
        "certificate": CERTIFICATE_12,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
    for row in rows:
        for name in PER_ROW_VISUAL_NAMES:
            assert (tmp_path / "visuals" / row.row_id / name).read_text(
                encoding="utf-8"
            ) == (VISUALS_ROOT / row.row_id / name).read_text(encoding="utf-8")


def test_visual_mechanism_row_helpers_keep_safe_shapes() -> None:
    _, packets, _ = _report_from_locked_data()
    rows = safe_visual_mechanism_rows(packets)
    metrics = contact_metric_rows(packets)
    cohorts = failure_cohort_rows(packets)

    assert len(rows) == 12
    assert len(metrics) == 12
    assert sum(int(row["row_count"]) for row in cohorts) == 12
    assert set(ROOT_OUTPUT_NAMES) == {
        "visual_mechanism_12_report.json",
        "visual_mechanism_12_rows.csv",
        "visual_mechanism_12_contact_metrics.csv",
        "visual_mechanism_12_failure_cohorts.csv",
        "visual_mechanism_12_dashboard.html",
        "visual_mechanism_12_certificate.json",
    }
    assert all("sequence" not in row for row in rows)
    assert all(row["raw_sequence_exposed"] is False for row in rows)
    assert all(row["global_folding_claim_allowed"] is False for row in rows)
