import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_contact_topology import (  # noqa: E402
    load_visual_mechanism_rows,
)
from pharmacotopology.folding_visual_mechanism_audit import (  # noqa: E402
    VISUAL_MECHANISM_AUDIT_CERTIFICATE_KIND,
    VISUAL_MECHANISM_AUDIT_KIND,
    build_visual_mechanism_audit_report,
    write_visual_mechanism_audit_outputs,
)


BENCHMARK_12 = ROOT / "data" / "folding_mechanism_visual_12.locked.json"
REL_BENCHMARK_12 = Path("data/folding_mechanism_visual_12.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
BASELINE_REPORT = RUN_DIR / "visual_mechanism_12_report.json"
REPAIR_REPORT = RUN_DIR / "contact_topology_repair_12_report.json"
REL_BASELINE_REPORT = Path(
    "first_contact_clean_pharmacotopology_layer_run/visual_mechanism_12_report.json"
)
REL_REPAIR_REPORT = Path(
    "first_contact_clean_pharmacotopology_layer_run/contact_topology_repair_12_report.json"
)
AUDIT_REPORT = RUN_DIR / "visual_mechanism_audit_report.json"
AUDIT_ROWS = RUN_DIR / "visual_mechanism_audit_rows.csv"
AUDIT_OVERFIT = RUN_DIR / "visual_mechanism_audit_overfit_risks.csv"
AUDIT_DASHBOARD = RUN_DIR / "visual_mechanism_audit_dashboard.html"
AUDIT_CERTIFICATE = RUN_DIR / "visual_mechanism_audit_certificate.json"


def _audit_report_from_checked_artifacts():
    return build_visual_mechanism_audit_report(
        baseline_report_path=REL_BASELINE_REPORT,
        repair_report_path=REL_REPAIR_REPORT,
        source_benchmark_file=REL_BENCHMARK_12,
    )


def test_locked_visual_12_is_marked_toy_and_non_discovery() -> None:
    data = json.loads(BENCHMARK_12.read_text(encoding="utf-8"))

    assert data["benchmark_scope"] == "toy_coarse_internal_contact_map_benchmark"
    assert data["visual_12_is_toy_benchmark"] is True
    assert data["coarse_native_contacts_only"] is True
    assert data["contact_repair_overfit_risk_expected"] is True
    assert data["mechanism_discovery_claim_allowed"] is False
    assert data["folding_problem_solved"] is False
    assert data["global_folding_claim_allowed"] is False


def test_visual_mechanism_audit_freezes_claim_boundary() -> None:
    report = _audit_report_from_checked_artifacts()

    assert report["audit_kind"] == VISUAL_MECHANISM_AUDIT_KIND
    assert report["visual_12_is_toy_benchmark"] is True
    assert report["coarse_native_contacts_only"] is True
    assert report["full_physical_native_contacts_available"] is False
    assert report["baseline_visible_partial_success_count"] == 5
    assert report["repaired_visible_partial_success_count"] == 8
    assert report["visible_partial_success_delta"] == 3
    assert report["baseline_visible_failure_count"] == 7
    assert report["repaired_visible_failure_count"] == 4
    assert report["baseline_mean_contact_map_f1"] == 0.104031
    assert report["repaired_mean_contact_map_f1"] == 0.263729
    assert report["hardcoded_beta_registry_pair_templates_detected"] is True
    assert report["hardcoded_beta_registry_pattern_families"] == (
        "compact_beta_registry_centers",
        "fixed_cross_sheet_anchors",
    )
    assert report["contact_repair_overfit_risk_reported"] is True
    assert report["contact_repair_overfit_risk_level"] == (
        "high_for_visual_12_benchmark"
    )
    assert report["overfit_risk_row_count"] == 5
    assert report["beta_template_success_gain_row_count"] == 3
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["artifact_reproducible"] is True
    assert report["clean_archive_required"] is True
    assert report["finder_zip_allowed"] is False


def test_checked_in_visual_mechanism_audit_outputs_have_expected_surfaces() -> None:
    report = json.loads(AUDIT_REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(AUDIT_CERTIFICATE.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(AUDIT_ROWS.read_text(encoding="utf-8").splitlines()))
    overfit_rows = list(
        csv.DictReader(AUDIT_OVERFIT.read_text(encoding="utf-8").splitlines())
    )
    dashboard = AUDIT_DASHBOARD.read_text(encoding="utf-8")

    assert report["audit_kind"] == VISUAL_MECHANISM_AUDIT_KIND
    assert certificate["certificate_kind"] == VISUAL_MECHANISM_AUDIT_CERTIFICATE_KIND
    assert len(rows) == 12
    assert len(overfit_rows) == 5
    assert "sequence" not in rows[0]
    assert "raw_sequence" not in rows[0]
    assert certificate["visual_12_is_toy_benchmark"] is True
    assert certificate["contact_repair_overfit_risk_reported"] is True
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert certificate["artifact_reproducible"] is True
    assert certificate["clean_archive_required"] is True
    assert certificate["finder_zip_allowed"] is False
    assert "Visual Mechanism Audit" in dashboard
    assert "Toy Benchmark" in dashboard
    assert "Overfit Risk Reported" in dashboard
    assert "No Discovery Claim" in dashboard
    assert "Archive Rule" in dashboard


def test_visual_mechanism_audit_reports_beta_template_risk_rows() -> None:
    overfit_rows = list(
        csv.DictReader(AUDIT_OVERFIT.read_text(encoding="utf-8").splitlines())
    )
    by_id = {row["row_id"]: row for row in overfit_rows}

    assert set(by_id) == {
        "visual_004_beta_csp",
        "visual_005_beta_sh3",
        "visual_006_beta_protein_l",
        "visual_007_mixed_protein_g",
        "visual_008_mixed_ubiquitin",
    }
    assert by_id["visual_004_beta_csp"]["overfit_risk_reason"] == (
        "beta_registry_template_created_success_gain_on_toy_benchmark"
    )
    assert by_id["visual_008_mixed_ubiquitin"]["overfit_risk_reason"] == (
        "beta_registry_template_still_drives_failure_or_compaction"
    )


def test_visual_mechanism_audit_outputs_do_not_export_raw_sequences() -> None:
    locked_rows = load_visual_mechanism_rows(BENCHMARK_12)
    generated_paths = [
        AUDIT_REPORT,
        AUDIT_ROWS,
        AUDIT_OVERFIT,
        AUDIT_DASHBOARD,
        AUDIT_CERTIFICATE,
    ]

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in locked_rows:
            assert row.sequence not in text


def test_visual_mechanism_audit_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    report = _audit_report_from_checked_artifacts()
    outputs = {
        "report": tmp_path / AUDIT_REPORT.name,
        "rows": tmp_path / AUDIT_ROWS.name,
        "overfit": tmp_path / AUDIT_OVERFIT.name,
        "dashboard": tmp_path / AUDIT_DASHBOARD.name,
        "certificate": tmp_path / AUDIT_CERTIFICATE.name,
    }
    write_visual_mechanism_audit_outputs(
        report=report,
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        overfit_risks_path=outputs["overfit"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": AUDIT_REPORT,
        "rows": AUDIT_ROWS,
        "overfit": AUDIT_OVERFIT,
        "dashboard": AUDIT_DASHBOARD,
        "certificate": AUDIT_CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")
