import csv
import json
from dataclasses import replace
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_contact_topology import predict_contact_topology  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    PER_ROW_VISUAL_NAMES,
    REAL_COORDINATE_NATIVE_KIND,
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    REAL_COORDINATE_VISUAL_CERTIFICATE_KIND,
    REAL_COORDINATE_VISUAL_SIGNATURE_KIND,
    ROOT_OUTPUT_NAMES,
    build_real_coordinate_visual_report,
    contact_metric_rows,
    load_real_coordinate_visual_rows,
    native_contact_summary_rows,
    parse_pdb_ca_coordinate_points,
    real_coordinate_visual_packets,
    safe_real_coordinate_visual_rows,
    validate_real_coordinate_visual_lock,
    write_real_coordinate_visual_outputs,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT_8 = RUN_DIR / "real_coordinate_visual_8_report.json"
ROWS_8 = RUN_DIR / "real_coordinate_visual_8_rows.csv"
CONTACT_METRICS_8 = RUN_DIR / "real_coordinate_visual_8_contact_metrics.csv"
NATIVE_SUMMARY_8 = RUN_DIR / "real_coordinate_visual_8_native_contact_summary.csv"
DASHBOARD_8 = RUN_DIR / "real_coordinate_visual_8_dashboard.html"
CERTIFICATE_8 = RUN_DIR / "real_coordinate_visual_8_certificate.json"
VISUALS_ROOT = RUN_DIR / "real_coordinate_visuals"


PDB_FIXTURE = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  CA  LEU A   2       3.800   0.000   0.000  1.00 10.00           C
ATOM      3  CA  GLY A   3       7.600   0.000   0.000  1.00 10.00           C
ATOM      4  CA  VAL A   4       0.000   4.000   0.000  1.00 10.00           C
ATOM      5  CA  ILE A   5       3.800   4.000   0.000  1.00 10.00           C
ATOM      6  CA  THR A   6       7.600   4.000   0.000  1.00 10.00           C
END
"""


def _report_from_locked_data(source_path: Path = REL_BENCHMARK_8):
    rows = load_real_coordinate_visual_rows(source_path)
    packets = real_coordinate_visual_packets(rows)
    report = build_real_coordinate_visual_report(
        packets=packets,
        source_benchmark_file=source_path,
        lock_validation=validate_real_coordinate_visual_lock(rows),
    )
    return rows, packets, report


def test_pdb_ca_points_can_derive_coordinate_contacts() -> None:
    points = parse_pdb_ca_coordinate_points(PDB_FIXTURE, chain_id="A")

    assert len(points) == 6
    assert points[0].sequence_index == 1
    assert points[-1].sequence_index == 6
    assert points[0].x == 0.0
    assert points[-1].y == 4.0


def test_real_coordinate_visual_report_meets_boundary_targets() -> None:
    rows, packets, report = _report_from_locked_data()

    assert len(rows) == 8
    assert report["benchmark_kind"] == REAL_COORDINATE_VISUAL_BENCHMARK_KIND
    assert report["real_coordinate_visual_signature_kind"] == (
        REAL_COORDINATE_VISUAL_SIGNATURE_KIND
    )
    assert report["prediction_input_boundary"] == (
        "sequence_only_no_native_contacts_no_truth_axes"
    )
    assert report["truth_scoring_boundary"] == (
        "coordinates_and_coordinate_native_contacts_used_only_after_sequence_prediction"
    )
    assert report["native_contact_derivation_kind"] == REAL_COORDINATE_NATIVE_KIND
    assert report["real_coordinate_native_contacts_extracted"] is True
    assert report["toy_locked_contact_targets_used"] is False
    assert report["coarse_ca_only"] is True
    assert report["full_atomic_folding_available"] is False
    assert report["benchmark_size"] == 8
    assert report["coordinate_backed_row_count"] == 8
    assert report["visual_artifacts_generated_for_rows"] == 8
    assert report["visual_files_per_row"] == 8
    assert report["visual_artifacts_generated_count"] == 70
    assert report["contact_map_f1_computed_count"] == 8
    assert report["mean_contact_map_f1"] == 0.051198
    assert report["max_contact_map_f1"] == 0.089776
    assert report["mean_native_contact_precision"] == 0.158448
    assert report["mean_native_contact_recall"] == 0.030616
    assert report["visible_partial_success_count"] == 3
    assert report["visible_failure_count"] == 5
    assert report["failure_cohorts"] == {
        "coordinate_architecture_gap": 1,
        "coordinate_beta_long_range_gap": 2,
        "coordinate_false_contact_overprediction": 1,
        "coordinate_low_recall_gap": 1,
        "coordinate_visible_partial_success": 3,
    }
    assert report["native_truth_used_before_prediction"] is False
    assert report["coordinate_truth_used_before_prediction"] is False
    assert report["raw_sequence_exposed"] is False
    assert report["repair_heuristic_applied"] is False
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert all(not packet.prediction.raw_sequence_exposed for packet in packets)
    assert all(
        not packet.prediction.native_truth_used_before_prediction
        for packet in packets
    )


def test_coordinate_truth_is_blind_to_prediction() -> None:
    row = load_real_coordinate_visual_rows(BENCHMARK_8)[0]
    original = predict_contact_topology(row.sequence, row_id=row.row_id)
    changed_coordinate_truth = replace(
        row,
        coordinate_points=tuple(reversed(row.coordinate_points)),
        coordinate_trace_hash="changed-after-prediction",
        native_contact_map_hash="changed-after-prediction",
    )
    repeated = predict_contact_topology(
        changed_coordinate_truth.sequence,
        row_id=changed_coordinate_truth.row_id,
    )

    assert original.predicted_contact_map_hash == repeated.predicted_contact_map_hash
    assert original.predicted_contact_pairs == repeated.predicted_contact_pairs
    assert original.native_truth_used_before_prediction is False
    assert repeated.native_truth_used_before_prediction is False


def test_checked_in_real_coordinate_visual_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_8.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE_8.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_8.read_text(encoding="utf-8").splitlines()))
    contact_metrics = list(
        csv.DictReader(CONTACT_METRICS_8.read_text(encoding="utf-8").splitlines())
    )
    native_summary = list(
        csv.DictReader(NATIVE_SUMMARY_8.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD_8.read_text(encoding="utf-8")
    visual_files = sorted(path for path in VISUALS_ROOT.glob("*/*") if path.is_file())

    assert report["benchmark_kind"] == REAL_COORDINATE_VISUAL_BENCHMARK_KIND
    assert certificate["certificate_kind"] == REAL_COORDINATE_VISUAL_CERTIFICATE_KIND
    assert len(rows) == 8
    assert len(contact_metrics) == 8
    assert len(native_summary) == 8
    assert len(visual_files) == 64
    assert "sequence" not in rows[0]
    assert "raw_sequence" not in rows[0]
    assert certificate["real_coordinate_native_contacts_extracted"] is True
    assert certificate["toy_locked_contact_targets_used"] is False
    assert certificate["coarse_ca_only"] is True
    assert certificate["native_truth_used_before_prediction"] is False
    assert certificate["coordinate_truth_used_before_prediction"] is False
    assert certificate["raw_sequence_exposed"] is False
    assert certificate["repair_heuristic_applied"] is False
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["global_folding_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Real Coordinate Native Contacts" in dashboard
    assert "Prediction Is Blind To Coordinates" in dashboard
    assert "C-alpha Coarse, Not Full Atomic Folding" in dashboard
    assert "No Repair Heuristic Applied" in dashboard
    assert "Global Folding Claim Remains Locked" in dashboard
    for row in rows:
        row_dir = VISUALS_ROOT / row["row_id"]
        for name in PER_ROW_VISUAL_NAMES:
            assert (row_dir / name).exists()


def test_real_coordinate_visual_outputs_do_not_export_raw_sequences() -> None:
    locked_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("real_coordinate_visual_8_*"))
    generated_paths.extend(path for path in VISUALS_ROOT.glob("*/*") if path.is_file())

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in locked_rows:
            assert row.sequence not in text


def test_real_coordinate_visual_artifacts_are_reproducible(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(ROOT)
    rows, packets, report = _report_from_locked_data()
    outputs = {
        "report": tmp_path / REPORT_8.name,
        "rows": tmp_path / ROWS_8.name,
        "contact_metrics": tmp_path / CONTACT_METRICS_8.name,
        "native_summary": tmp_path / NATIVE_SUMMARY_8.name,
        "dashboard": tmp_path / DASHBOARD_8.name,
        "certificate": tmp_path / CERTIFICATE_8.name,
    }
    write_real_coordinate_visual_outputs(
        report=report,
        packets=packets,
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        contact_metrics_path=outputs["contact_metrics"],
        native_contact_summary_path=outputs["native_summary"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
        visuals_root=tmp_path / "real_coordinate_visuals",
    )

    checked_in = {
        "report": REPORT_8,
        "rows": ROWS_8,
        "contact_metrics": CONTACT_METRICS_8,
        "native_summary": NATIVE_SUMMARY_8,
        "dashboard": DASHBOARD_8,
        "certificate": CERTIFICATE_8,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
    for row in rows:
        for name in PER_ROW_VISUAL_NAMES:
            assert (tmp_path / "real_coordinate_visuals" / row.row_id / name).read_text(
                encoding="utf-8"
            ) == (VISUALS_ROOT / row.row_id / name).read_text(encoding="utf-8")


def test_real_coordinate_row_helpers_keep_safe_shapes() -> None:
    _, packets, _ = _report_from_locked_data()
    rows = safe_real_coordinate_visual_rows(packets)
    metrics = contact_metric_rows(packets)
    summaries = native_contact_summary_rows(packets)

    assert len(rows) == 8
    assert len(metrics) == 8
    assert len(summaries) == 8
    assert set(ROOT_OUTPUT_NAMES) == {
        "real_coordinate_visual_8_report.json",
        "real_coordinate_visual_8_rows.csv",
        "real_coordinate_visual_8_contact_metrics.csv",
        "real_coordinate_visual_8_native_contact_summary.csv",
        "real_coordinate_visual_8_dashboard.html",
        "real_coordinate_visual_8_certificate.json",
    }
    assert all("sequence" not in row for row in rows)
    assert all(row["raw_sequence_exposed"] is False for row in rows)
    assert all(row["repair_heuristic_applied"] is False for row in rows)
    assert all(row["global_folding_claim_allowed"] is False for row in rows)
