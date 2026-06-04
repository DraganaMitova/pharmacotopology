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
    load_visual_mechanism_rows,
    predict_contact_topology,
    validate_visual_mechanism_lock,
)
from pharmacotopology.folding_contact_topology_repair_benchmark import (  # noqa: E402
    CONTACT_TOPOLOGY_REPAIR_BENCHMARK_KIND,
    CONTACT_TOPOLOGY_REPAIR_CERTIFICATE_KIND,
    PER_ROW_REPAIR_VISUAL_NAMES,
    build_contact_topology_repair_report,
    contact_topology_repair_packets,
    repair_failure_cohort_rows,
    repair_gap_rows,
    safe_repair_rows,
    write_contact_topology_repair_outputs,
)
from pharmacotopology.folding_long_range_contact_repair import (  # noqa: E402
    REPAIR_INPUT_BOUNDARY,
    REPAIRED_CONTACT_TOPOLOGY_SIGNATURE_KIND,
    repair_contact_topology,
)


BENCHMARK_12 = ROOT / "data" / "folding_mechanism_visual_12.locked.json"
REL_BENCHMARK_12 = Path("data/folding_mechanism_visual_12.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT_12 = RUN_DIR / "contact_topology_repair_12_report.json"
ROWS_12 = RUN_DIR / "contact_topology_repair_12_rows.csv"
GAP_ANALYSIS_12 = RUN_DIR / "contact_topology_repair_12_gap_analysis.csv"
FAILURE_COHORTS_12 = RUN_DIR / "contact_topology_repair_12_failure_cohorts.csv"
DASHBOARD_12 = RUN_DIR / "contact_topology_repair_12_dashboard.html"
CERTIFICATE_12 = RUN_DIR / "contact_topology_repair_12_certificate.json"
VISUALS_ROOT = RUN_DIR / "contact_repair_visuals"


def _report_from_locked_data(source_path: Path = REL_BENCHMARK_12):
    rows = load_visual_mechanism_rows(source_path)
    packets = contact_topology_repair_packets(rows)
    report = build_contact_topology_repair_report(
        packets=packets,
        source_benchmark_file=source_path,
        lock_validation=validate_visual_mechanism_lock(rows),
    )
    return rows, packets, report


def test_contact_topology_repair_improves_visible_and_long_range_metrics() -> None:
    rows, packets, report = _report_from_locked_data()

    assert len(rows) == 12
    assert report["benchmark_kind"] == CONTACT_TOPOLOGY_REPAIR_BENCHMARK_KIND
    assert report["benchmark_size"] == 12
    assert report["visual_artifacts_generated_for_rows"] == 12
    assert report["visual_artifacts_generated_count"] == 90
    assert report["contact_map_f1_computed_count"] == 12
    assert report["baseline_visible_partial_success_count"] == 5
    assert report["visible_partial_success_count"] == 8
    assert report["visible_partial_success_delta"] == 3
    assert report["baseline_visible_failure_count"] == 7
    assert report["visible_failure_count"] == 4
    assert report["baseline_mean_contact_map_f1"] == 0.104031
    assert report["repaired_mean_contact_map_f1"] == 0.263729
    assert report["mean_contact_map_f1_delta"] == 0.159698
    assert report["baseline_mean_long_range_contact_recall"] == 0.333333
    assert report["repaired_mean_long_range_contact_recall"] == 0.805556
    assert report["long_range_contact_recall_delta"] == 0.472223
    assert report["baseline_mean_beta_pairing_contact_recall"] == 0.018519
    assert report["repaired_mean_beta_pairing_contact_recall"] == 0.672619
    assert report["beta_pairing_contact_recall_delta"] == 0.6541
    assert report["premature_compaction_count"] == 5
    assert report["visual_failure_cohort_count"] == 4
    assert report["failure_cohorts"] == {
        "bad_beta_pairing": 1,
        "disorder_over_collapse": 1,
        "membrane_mis_topology": 1,
        "premature_compaction": 1,
        "visible_partial_success": 8,
    }
    assert report["native_truth_used_before_prediction"] is False
    assert report["native_truth_used_before_repair"] is False
    assert report["raw_sequence_exposed"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert all(not packet.repair.native_truth_used_before_repair for packet in packets)
    assert all(not packet.repair.raw_sequence_exposed for packet in packets)


def test_repair_prediction_is_sequence_only_before_native_scoring() -> None:
    row = load_visual_mechanism_rows(BENCHMARK_12)[3]
    baseline = predict_contact_topology(row.sequence, row_id=row.row_id)
    original = repair_contact_topology(row.sequence, baseline_prediction=baseline)
    changed_truth = replace(
        row,
        native_contact_pairs=((1, row.length),),
        native_contact_map_hash="changed-after-prediction",
        truth_axes={
            "secondary_structure_axis": "alpha_rich",
            "architecture_axis": "multidomain_or_segmented",
            "order_axis": "disordered_flexible",
            "environment_axis": "membrane_like",
        },
    )
    repeated_baseline = predict_contact_topology(
        changed_truth.sequence,
        row_id=changed_truth.row_id,
    )
    repeated = repair_contact_topology(
        changed_truth.sequence,
        baseline_prediction=repeated_baseline,
    )

    assert original.repaired_prediction.contact_topology_signature_kind == (
        REPAIRED_CONTACT_TOPOLOGY_SIGNATURE_KIND
    )
    assert original.repair_input_boundary == REPAIR_INPUT_BOUNDARY
    assert original.repaired_prediction.predicted_contact_map_hash == (
        repeated.repaired_prediction.predicted_contact_map_hash
    )
    assert original.repaired_prediction.predicted_contact_pairs == (
        repeated.repaired_prediction.predicted_contact_pairs
    )
    assert original.native_truth_used_before_repair is False
    assert repeated.native_truth_used_before_repair is False


def test_checked_in_contact_repair_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT_12.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE_12.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS_12.read_text(encoding="utf-8").splitlines()))
    gap_rows = list(
        csv.DictReader(GAP_ANALYSIS_12.read_text(encoding="utf-8").splitlines())
    )
    failure_rows = list(
        csv.DictReader(FAILURE_COHORTS_12.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD_12.read_text(encoding="utf-8")
    visual_files = sorted(path for path in VISUALS_ROOT.glob("*/*") if path.is_file())

    assert report["benchmark_kind"] == CONTACT_TOPOLOGY_REPAIR_BENCHMARK_KIND
    assert certificate["certificate_kind"] == CONTACT_TOPOLOGY_REPAIR_CERTIFICATE_KIND
    assert len(rows) == 12
    assert len(gap_rows) == 12
    assert len(failure_rows) == 5
    assert len(visual_files) == 84
    assert "sequence" not in rows[0]
    assert "raw_sequence" not in rows[0]
    assert certificate["raw_sequence_exposed"] is False
    assert certificate["native_truth_used_before_prediction"] is False
    assert certificate["native_truth_used_before_repair"] is False
    assert certificate["global_folding_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Contact Topology Repair Workbench" in dashboard
    assert "Native Gap Analysis Happens After Prediction" in dashboard
    assert "Long-Range Repair Is Sequence-Only" in dashboard
    assert "Failures Stay Visible" in dashboard
    assert "Global Folding Claim Remains Locked" in dashboard
    for row in rows:
        row_dir = VISUALS_ROOT / row["row_id"]
        for name in PER_ROW_REPAIR_VISUAL_NAMES:
            assert (row_dir / name).exists()


def test_contact_repair_outputs_do_not_export_raw_sequences() -> None:
    locked_rows = load_visual_mechanism_rows(BENCHMARK_12)
    generated_paths = list(RUN_DIR.glob("contact_topology_repair_12_*"))
    generated_paths.extend(path for path in VISUALS_ROOT.glob("*/*") if path.is_file())

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in locked_rows:
            assert row.sequence not in text


def test_contact_repair_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    rows, packets, report = _report_from_locked_data()
    outputs = {
        "report": tmp_path / REPORT_12.name,
        "rows": tmp_path / ROWS_12.name,
        "gap": tmp_path / GAP_ANALYSIS_12.name,
        "failure": tmp_path / FAILURE_COHORTS_12.name,
        "dashboard": tmp_path / DASHBOARD_12.name,
        "certificate": tmp_path / CERTIFICATE_12.name,
    }
    write_contact_topology_repair_outputs(
        report=report,
        packets=packets,
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        gap_analysis_path=outputs["gap"],
        failure_cohorts_path=outputs["failure"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
        visuals_root=tmp_path / "contact_repair_visuals",
    )

    checked_in = {
        "report": REPORT_12,
        "rows": ROWS_12,
        "gap": GAP_ANALYSIS_12,
        "failure": FAILURE_COHORTS_12,
        "dashboard": DASHBOARD_12,
        "certificate": CERTIFICATE_12,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
    for row in rows:
        for name in PER_ROW_REPAIR_VISUAL_NAMES:
            assert (
                tmp_path / "contact_repair_visuals" / row.row_id / name
            ).read_text(encoding="utf-8") == (
                VISUALS_ROOT / row.row_id / name
            ).read_text(encoding="utf-8")


def test_contact_repair_row_helpers_keep_safe_shapes() -> None:
    _, packets, _ = _report_from_locked_data()
    rows = safe_repair_rows(packets)
    gap_rows = repair_gap_rows(packets)
    cohorts = repair_failure_cohort_rows(packets)

    assert len(rows) == 12
    assert len(gap_rows) == 12
    assert sum(int(row["repaired_row_count"]) for row in cohorts) == 12
    assert all("sequence" not in row for row in rows)
    assert all(row["raw_sequence_exposed"] is False for row in rows)
    assert all(row["global_folding_claim_allowed"] is False for row in rows)
    assert any(
        row["repaired_failure_mechanism"] == "premature_compaction"
        for row in rows
    )
