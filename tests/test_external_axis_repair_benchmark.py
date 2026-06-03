import csv
import json
from dataclasses import replace
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_architecture_axis import (  # noqa: E402
    architecture_evidence_packet_from_sequence,
)
from pharmacotopology.folding_external_axis_repair import (  # noqa: E402
    EXTERNAL_AXIS_REPAIR_BENCHMARK_KIND,
    build_external_axis_repair_report,
    external_axis_repair_abstention_delta_rows,
    external_axis_repair_conflict_delta_rows,
    external_axis_repair_family_summary_rows,
    external_axis_repair_quarantine_rows,
    external_axis_repair_rows,
    write_external_axis_repair_outputs,
)
from pharmacotopology.folding_external_holdout import (  # noqa: E402
    DEFAULT_DEVELOPMENT_BENCHMARK_FILE,
    load_external_holdout_rows,
)
from pharmacotopology.folding_order_axis_safety import (  # noqa: E402
    order_axis_safety_packet_from_source,
)


HOLDOUT_FILE = ROOT / "data" / "folding_benchmarks_external_fold_family_100.locked.json"
REL_HOLDOUT_FILE = Path("data/folding_benchmarks_external_fold_family_100.locked.json")
DEVELOPMENT_FILE = ROOT / DEFAULT_DEVELOPMENT_BENCHMARK_FILE
REL_DEVELOPMENT_FILE = DEFAULT_DEVELOPMENT_BENCHMARK_FILE
REPORT = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_report.json"
ROWS = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_rows.csv"
CONFLICT_DELTA = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_conflict_delta.csv"
ABSTENTION_DELTA = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_abstention_delta.csv"
QUARANTINE_ROWS = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_quarantine_rows.csv"
FAMILY_SUMMARY = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_family_summary.csv"
DASHBOARD = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_dashboard.html"
CERTIFICATE = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_axis_repair_certificate.json"


def test_external_axis_repair_closes_unsafe_axis_claims_by_abstention() -> None:
    holdout_rows = load_external_holdout_rows(HOLDOUT_FILE)
    report = build_external_axis_repair_report(
        holdout_rows,
        holdout_file=HOLDOUT_FILE,
        development_benchmark_file=DEVELOPMENT_FILE,
    )

    assert report["benchmark_kind"] == EXTERNAL_AXIS_REPAIR_BENCHMARK_KIND
    assert report["holdout_row_count"] == 100
    assert report["pre_repair_axis_profile_coverage"] == 0.78
    assert report["post_repair_axis_profile_coverage"] == 0.55
    assert report["post_repair_axis_profile_coverage"] < report[
        "pre_repair_axis_profile_coverage"
    ]
    assert report["pre_repair_axis_profile_same_axis_conflict_count"] == 22
    assert report["post_repair_axis_profile_same_axis_conflict_count"] == 0
    assert report["pre_repair_architecture_axis_same_axis_conflict_count"] == 14
    assert report["post_repair_architecture_axis_same_axis_conflict_count"] == 0
    assert report["pre_repair_unsafe_axis_claim_count"] == 36
    assert report["post_repair_unsafe_axis_claim_count"] == 0
    assert report["pre_repair_high_confidence_wrong_count_after_axis_scoring"] == 17
    assert report["post_repair_high_confidence_wrong_count_after_axis_scoring"] == 0
    assert report["order_axis_folded_mimic_quarantined_count"] == 22
    assert report["repeat_compact_ambiguity_quarantined_count"] == 14
    assert report["coverage_loss_from_external_safety"] == 0.23
    assert report["external_safety_repair_successful"] is True
    assert report["legacy_axis_artifacts_reproducible"] is True
    assert report["legacy_axis_profile_artifacts_reproducible"] is True
    assert report["unsafe_class_recovery_count"] == 0
    assert report["guard_override_count"] == 0
    assert report["global_fold_class_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False


def test_order_axis_safety_quarantines_folded_mimics_without_flipping_to_ordered() -> None:
    rows = {row.row_id: row for row in load_external_holdout_rows(HOLDOUT_FILE)}
    folded_mimic = rows["holdout_023_beta_rich_01"]
    true_disorder = rows["holdout_076_disordered_or_mixed_order_01"]
    source_row = {
        "predicted_fold_class": "disordered_flexible",
        "confidence": 0.60,
        "protein_regime": "intrinsically_disordered",
        "forced_prediction": True,
    }

    mimic_packet = order_axis_safety_packet_from_source(
        folded_mimic.sequence,
        source_row,
        row_id=folded_mimic.row_id,
    )
    disorder_packet = order_axis_safety_packet_from_source(
        true_disorder.sequence,
        source_row,
        row_id=true_disorder.row_id,
    )

    assert mimic_packet.order_axis_prediction == "mixed_or_uncertain"
    assert mimic_packet.order_axis_claim_allowed is False
    assert mimic_packet.order_axis_abstention_reason == (
        "external_order_axis_folded_mimic_quarantine"
    )
    assert disorder_packet.order_axis_prediction == "disordered_flexible"
    assert disorder_packet.order_axis_claim_allowed is True
    assert disorder_packet.order_axis_abstention_reason == ""


def test_architecture_repeat_requires_recurrence_not_hydrophobic_periodicity_only() -> None:
    rows = {row.row_id: row for row in load_external_holdout_rows(HOLDOUT_FILE)}
    compact_periodic = rows["holdout_003_alpha_rich_soluble_01"]

    baseline = architecture_evidence_packet_from_sequence(
        compact_periodic.sequence,
        protein_id=compact_periodic.row_id,
    )
    repaired = architecture_evidence_packet_from_sequence(
        compact_periodic.sequence,
        protein_id=compact_periodic.row_id,
        external_safe_quarantine=True,
    )

    assert baseline.architecture_axis_prediction == "repeat_like"
    assert baseline.architecture_axis_claim_allowed is True
    assert repaired.architecture_axis_prediction == "unknown"
    assert repaired.architecture_axis_claim_allowed is False
    assert repaired.architecture_axis_abstention_reason == (
        "repeat_compact_single_domain_ambiguity_quarantine"
    )
    assert repaired.repeat_recurrence_support < 0.70
    assert repaired.hydrophobic_periodicity_only_risk >= 0.40


def test_external_axis_repair_prediction_does_not_read_truth_axes() -> None:
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
    original_prediction = external_axis_repair_rows([original])[0]
    changed_prediction = external_axis_repair_rows([changed_truth])[0]
    prediction_keys = (
        "source_predicted_fold_class",
        "source_forced_prediction",
        "source_confidence",
        "protein_regime",
        "post_profile_order_axis",
        "order_axis_abstention_reason",
        "post_architecture_axis_prediction",
        "architecture_axis_abstention_reason",
        "profile_secondary_structure_axis",
        "profile_environment_axis",
    )
    for key in prediction_keys:
        assert original_prediction[key] == changed_prediction[key]
    assert original_prediction["truth_order_axis"] != changed_prediction[
        "truth_order_axis"
    ]


def test_external_axis_repair_artifacts_do_not_expose_raw_sequences() -> None:
    holdout_rows = load_external_holdout_rows(HOLDOUT_FILE)
    artifact_paths = (
        REPORT,
        ROWS,
        CONFLICT_DELTA,
        ABSTENTION_DELTA,
        QUARANTINE_ROWS,
        FAMILY_SUMMARY,
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


def test_tracked_external_axis_repair_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS.read_text(encoding="utf-8").splitlines()))
    conflict_delta = list(
        csv.DictReader(CONFLICT_DELTA.read_text(encoding="utf-8").splitlines())
    )
    abstention_delta = list(
        csv.DictReader(ABSTENTION_DELTA.read_text(encoding="utf-8").splitlines())
    )
    quarantine_rows = list(
        csv.DictReader(QUARANTINE_ROWS.read_text(encoding="utf-8").splitlines())
    )
    families = list(
        csv.DictReader(FAMILY_SUMMARY.read_text(encoding="utf-8").splitlines())
    )
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["benchmark_kind"] == EXTERNAL_AXIS_REPAIR_BENCHMARK_KIND
    assert certificate["external_axis_repair_complete"] is True
    assert certificate["legacy_axis_artifacts_reproducible"] is True
    assert certificate["legacy_axis_profile_artifacts_reproducible"] is True
    assert certificate["order_axis_folded_mimic_quarantine_active"] is True
    assert certificate["repeat_compact_ambiguity_quarantine_active"] is True
    assert certificate["post_repair_axis_profile_same_axis_conflict_count"] == 0
    assert certificate["post_repair_architecture_axis_same_axis_conflict_count"] == 0
    assert certificate["post_repair_unsafe_axis_claim_count"] == 0
    assert certificate["post_repair_high_confidence_wrong_count_after_axis_scoring"] == 0
    assert certificate["unsafe_class_recovery_count"] == 0
    assert certificate["guard_override_count"] == 0
    assert certificate["global_fold_class_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert certificate["claim_allowed"] is False
    assert len(rows) == 100
    assert len(conflict_delta) == 4
    assert len(abstention_delta) == 3
    assert len(quarantine_rows) == 36
    assert len(families) == 23
    assert "External Holdout Safety Repair" in dashboard
    assert "Coverage Was Too Optimistic" in dashboard
    assert "Disorder Claim Requires Strong Disorder Evidence" in dashboard
    assert "Folded Beta/Mixed Mimics Are Quarantined" in dashboard
    assert "Repeat-Like Requires Recurrence" in dashboard
    assert "Hydrophobic Periodicity Alone Is Not Architecture" in dashboard
    assert "Safety Improved By Abstention" in dashboard
    assert "Global Fold Class Still Locked" in dashboard


def test_external_axis_repair_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    holdout_rows = load_external_holdout_rows(REL_HOLDOUT_FILE)
    rows = external_axis_repair_rows(holdout_rows)
    report = build_external_axis_repair_report(
        holdout_rows,
        holdout_file=REL_HOLDOUT_FILE,
        development_benchmark_file=REL_DEVELOPMENT_FILE,
    )
    outputs = {
        "report": tmp_path / REPORT.name,
        "rows": tmp_path / ROWS.name,
        "conflict_delta": tmp_path / CONFLICT_DELTA.name,
        "abstention_delta": tmp_path / ABSTENTION_DELTA.name,
        "quarantine": tmp_path / QUARANTINE_ROWS.name,
        "family": tmp_path / FAMILY_SUMMARY.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }

    write_external_axis_repair_outputs(
        report=report,
        rows=rows,
        conflict_delta_rows=external_axis_repair_conflict_delta_rows(report),
        abstention_delta_rows=external_axis_repair_abstention_delta_rows(report),
        quarantine_rows=external_axis_repair_quarantine_rows(rows),
        family_rows=external_axis_repair_family_summary_rows(rows),
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        conflict_delta_path=outputs["conflict_delta"],
        abstention_delta_path=outputs["abstention_delta"],
        quarantine_rows_path=outputs["quarantine"],
        family_summary_path=outputs["family"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "rows": ROWS,
        "conflict_delta": CONFLICT_DELTA,
        "abstention_delta": ABSTENTION_DELTA,
        "quarantine": QUARANTINE_ROWS,
        "family": FAMILY_SUMMARY,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
