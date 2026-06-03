import csv
import json
from dataclasses import replace
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_external_holdout import (  # noqa: E402
    EXTERNAL_HOLDOUT_BENCHMARK_KIND,
    EXTERNAL_HOLDOUT_SPLIT,
    DEFAULT_DEVELOPMENT_BENCHMARK_FILE,
    abstention_rows,
    axis_conflict_rows,
    build_external_holdout_report,
    external_holdout_rows,
    failure_cohort_rows,
    family_summary_rows,
    load_external_holdout_rows,
    validate_holdout_lock,
    write_external_holdout_outputs,
)


HOLDOUT_FILE = ROOT / "data" / "folding_benchmarks_external_fold_family_100.locked.json"
REL_HOLDOUT_FILE = Path("data/folding_benchmarks_external_fold_family_100.locked.json")
DEVELOPMENT_FILE = ROOT / DEFAULT_DEVELOPMENT_BENCHMARK_FILE
REL_DEVELOPMENT_FILE = DEFAULT_DEVELOPMENT_BENCHMARK_FILE
REPORT = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_report.json"
ROWS = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_rows.csv"
FAMILY_SUMMARY = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_family_summary.csv"
CONFLICTS = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_axis_conflicts.csv"
ABSTENTIONS = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_abstentions.csv"
FAILURE_COHORTS = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_failure_cohorts.csv"
DASHBOARD = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_dashboard.html"
CERTIFICATE = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_fold_family_100_certificate.json"


def test_holdout_lock_has_100_non_overlapping_rows() -> None:
    rows = load_external_holdout_rows(HOLDOUT_FILE)
    lock = validate_holdout_lock(
        rows,
        development_benchmark_file=DEVELOPMENT_FILE,
    )

    assert len(rows) == 100
    assert all(row.sequence_sha256 for row in rows)
    assert all(row.sequence_sha256 == row.sequence_sha256.lower() for row in rows)
    assert all(row.holdout_split == EXTERNAL_HOLDOUT_SPLIT for row in rows)
    assert lock["holdout_row_count"] == 100
    assert lock["holdout_unique_sequence_count"] == 100
    assert lock["holdout_unique_family_count"] == 23
    assert lock["development_overlap_count"] == 0
    assert lock["development_sequence_overlap_count"] == 0
    assert lock["development_family_overlap_count"] == 0
    assert lock["holdout_lock_valid"] is True
    assert lock["holdout_non_overlap_valid"] is True


def test_external_holdout_report_tracks_generalization_failures() -> None:
    rows = load_external_holdout_rows(HOLDOUT_FILE)
    report = build_external_holdout_report(
        rows,
        holdout_file=HOLDOUT_FILE,
        development_benchmark_file=DEVELOPMENT_FILE,
    )

    assert report["benchmark_kind"] == EXTERNAL_HOLDOUT_BENCHMARK_KIND
    assert report["holdout_row_count"] == 100
    assert report["holdout_lock_valid"] is True
    assert report["holdout_non_overlap_valid"] is True
    assert report["collapsed_class_coverage"] == 0.37
    assert report["axis_profile_coverage"] == 0.78
    assert report["secondary_axis_coverage"] == 0.05
    assert report["architecture_axis_coverage"] == 0.17
    assert report["order_axis_coverage"] == 0.76
    assert report["environment_axis_coverage"] == 0.17
    assert report["axis_profile_same_axis_conflict_count"] == 22
    assert report["architecture_axis_same_axis_conflict_count"] == 14
    assert report["forced_same_axis_conflict_count"] == 22
    assert report["high_confidence_wrong_count_after_axis_scoring"] == 17
    assert report["safe_axis_claim_count"] == 79
    assert report["unsafe_axis_claim_count"] == 36
    assert report["safe_axis_recovered_count"] == 78
    assert report["unsafe_class_recovery_count"] == 0
    assert report["guard_override_count"] == 0
    assert report["family_level_failure_count"] == 14
    assert report["family_level_abstention_count"] == 22
    assert report["family_level_conflict_count"] == 13
    assert report["family_generalization_status"] == "conflicts_detected"
    assert report["global_fold_class_claim_allowed"] is False
    assert report["axis_profile_claim_allowed"] is True
    assert report["architecture_axis_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False
    assert report["prediction_logic_changed_in_this_batch"] is False
    assert report["thresholds_changed_in_this_batch"] is False


def test_prediction_outputs_do_not_read_truth_axes_before_prediction() -> None:
    original = load_external_holdout_rows(HOLDOUT_FILE)[0]
    changed_truth = replace(
        original,
        truth_axes={
            "secondary_structure_axis": "beta_rich",
            "architecture_axis": "repeat_like",
            "order_axis": "disordered_flexible",
            "environment_axis": "membrane_like",
        },
    )

    original_prediction = external_holdout_rows([original])[0]
    changed_prediction = external_holdout_rows([changed_truth])[0]
    prediction_keys = (
        "source_predicted_fold_class",
        "source_forced_prediction",
        "source_abstained",
        "source_confidence",
        "protein_regime",
        "gate_path",
        "profile_secondary_structure_axis",
        "profile_order_axis",
        "profile_environment_axis",
        "profile_architecture_axis",
        "architecture_axis_prediction",
        "architecture_axis_confidence",
        "architecture_axis_claim_allowed",
        "architecture_axis_abstention_reason",
    )
    for key in prediction_keys:
        assert original_prediction[key] == changed_prediction[key]
    assert original_prediction["truth_secondary_structure_axis"] != changed_prediction[
        "truth_secondary_structure_axis"
    ]


def test_external_holdout_artifacts_do_not_expose_raw_sequences() -> None:
    holdout_rows = load_external_holdout_rows(HOLDOUT_FILE)
    artifact_paths = (
        REPORT,
        ROWS,
        FAMILY_SUMMARY,
        CONFLICTS,
        ABSTENTIONS,
        FAILURE_COHORTS,
        DASHBOARD,
        CERTIFICATE,
    )
    artifact_text = "\n".join(
        path.read_text(encoding="utf-8") for path in artifact_paths
    )
    for row in holdout_rows:
        assert row.sequence not in artifact_text
    assert "sequence_sha256" in ROWS.read_text(encoding="utf-8")
    assert "sequence" not in next(
        csv.DictReader(ROWS.read_text(encoding="utf-8").splitlines())
    )


def test_tracked_external_holdout_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS.read_text(encoding="utf-8").splitlines()))
    families = list(
        csv.DictReader(FAMILY_SUMMARY.read_text(encoding="utf-8").splitlines())
    )
    conflicts = list(csv.DictReader(CONFLICTS.read_text(encoding="utf-8").splitlines()))
    abstentions = list(
        csv.DictReader(ABSTENTIONS.read_text(encoding="utf-8").splitlines())
    )
    cohorts = list(
        csv.DictReader(FAILURE_COHORTS.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["benchmark_kind"] == EXTERNAL_HOLDOUT_BENCHMARK_KIND
    assert certificate["external_fold_family_holdout_complete"] is True
    assert certificate["holdout_row_count"] == report["holdout_row_count"] == 100
    assert certificate["development_overlap_count"] == report[
        "development_overlap_count"
    ] == 0
    assert certificate["development_sequence_overlap_count"] == report[
        "development_sequence_overlap_count"
    ] == 0
    assert certificate["holdout_non_overlap_valid"] is True
    assert certificate["prediction_logic_changed_in_this_batch"] is False
    assert certificate["thresholds_changed_in_this_batch"] is False
    assert certificate["guard_overrides"] == 0
    assert certificate["unsafe_class_recovery_count"] == 0
    assert certificate["global_fold_class_claim_allowed"] is False
    assert certificate["axis_profile_claim_allowed"] is True
    assert certificate["architecture_axis_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert certificate["claim_allowed"] is False
    assert len(rows) == 100
    assert len(families) == 23
    assert len(conflicts) == 36
    assert len(abstentions) == 83
    assert len(cohorts) == 24
    assert "External Holdout, Not Development Benchmark" in dashboard
    assert "No Threshold Tuning In This Batch" in dashboard
    assert "Family Non-Overlap Check" in dashboard
    assert "Axis Profile Generalization" in dashboard
    assert "Architecture Axis Generalization" in dashboard
    assert "Abstention Is Allowed" in dashboard
    assert "Unsafe Class Recovery Remains Forbidden" in dashboard
    assert "Global Fold Class Still Locked" in dashboard
    assert "Dominant Failure Cohorts" in dashboard
    assert "Next Repair Candidates" in dashboard


def test_external_holdout_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    holdout_rows = load_external_holdout_rows(REL_HOLDOUT_FILE)
    rows = external_holdout_rows(holdout_rows)
    report = build_external_holdout_report(
        holdout_rows,
        holdout_file=REL_HOLDOUT_FILE,
        development_benchmark_file=REL_DEVELOPMENT_FILE,
    )
    outputs = {
        "report": tmp_path / REPORT.name,
        "rows": tmp_path / ROWS.name,
        "family": tmp_path / FAMILY_SUMMARY.name,
        "conflicts": tmp_path / CONFLICTS.name,
        "abstentions": tmp_path / ABSTENTIONS.name,
        "cohorts": tmp_path / FAILURE_COHORTS.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }

    write_external_holdout_outputs(
        report=report,
        rows=rows,
        family_rows=family_summary_rows(rows),
        conflicts=axis_conflict_rows(rows),
        abstentions=abstention_rows(rows),
        failure_cohorts=failure_cohort_rows(rows),
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        family_summary_path=outputs["family"],
        axis_conflicts_path=outputs["conflicts"],
        abstentions_path=outputs["abstentions"],
        failure_cohorts_path=outputs["cohorts"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "rows": ROWS,
        "family": FAMILY_SUMMARY,
        "conflicts": CONFLICTS,
        "abstentions": ABSTENTIONS,
        "cohorts": FAILURE_COHORTS,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
