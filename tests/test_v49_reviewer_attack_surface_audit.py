from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_v49_reviewer_attack_surface_audit_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v49_audit_passes_and_writes_certificate() -> None:
    runner = _load(RUNNER, "v49_runner_pass")
    paths = runner.run_v49()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["status"] == "V49_REVIEWER_ATTACK_SURFACE_HARDENING_PASSED"
    assert cert["evidence_taxonomy_present"] is True
    assert cert["spatial_proxy_boundary_present"] is True
    assert cert["operator_scoring_rubric_present"] is True
    assert cert["latest_completed_cycle_documented"] is True
    assert cert["negative_control_inventory_present"] is True
    assert cert["v46_full_cycle_visible"] is True
    assert cert["readme_claim_boundary_preserved"] is True
    assert cert["passed_control_count"] == cert["control_count"]
    assert cert["negative_null_control_count"] >= 8
    assert cert["folding_problem_solved"] is False
    assert cert["claim_allowed"] is True
    assert paths["report"].exists()


def test_v49_docs_freeze_taxonomy_source_fields_and_scoring_labels() -> None:
    runner = _load(RUNNER, "v49_runner_docs")
    taxonomy = (ROOT / "docs" / "EVIDENCE_TAXONOMY.md").read_text(encoding="utf-8")
    rubric = (ROOT / "docs" / "OPERATOR_SCORING_RUBRIC.md").read_text(encoding="utf-8")
    for category in runner.EVIDENCE_CLASSES:
        assert category in taxonomy
    for field in runner.REQUIRED_SOURCE_ROLE_FIELDS:
        assert field in taxonomy
    for field in runner.REQUIRED_SCORE_FIELDS:
        assert field in rubric
    for label in runner.REQUIRED_SCORE_LABELS:
        assert label in rubric
    assert "never biological prediction evidence" in taxonomy
    assert "never biological claim evidence" in taxonomy


def test_v49_readme_answers_the_four_reviewer_questions() -> None:
    runner = _load(RUNNER, "v49_runner_readme")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "## Evidence Boundary" in readme
    assert "## Operator Scoring" in readme
    assert "## Latest Completed Full Cycle" in readme
    assert "## Negative / Null Controls" in readme
    assert "pure_non_coordinate" in readme
    assert "spatial_proxy_non_coordinate" in readme
    assert "spatial-proxy evidence can encode spatial information" in readme
    assert "coordinate-derived sources are blocked before sealing" in readme.lower()
    assert runner.V46_CERT_REL in readme
    assert runner.V46_REPORT_REL in readme
    assert "V46 completed seal -> holdout -> validation" in readme
    assert len(set(runner.negative_control_ids_from_readme(readme))) >= 8


def test_v49_audit_fails_when_required_readme_section_is_missing() -> None:
    runner = _load(RUNNER, "v49_runner_failure")
    readme = (ROOT / "README.md").read_text(encoding="utf-8").replace("## Evidence Boundary", "## Evidence Boundary Removed")
    taxonomy = (ROOT / "docs" / "EVIDENCE_TAXONOMY.md").read_text(encoding="utf-8")
    rubric = (ROOT / "docs" / "OPERATOR_SCORING_RUBRIC.md").read_text(encoding="utf-8")
    inventory = (ROOT / "docs" / "NEGATIVE_CONTROL_INVENTORY.md").read_text(encoding="utf-8")
    cert = runner.audit_texts(readme, taxonomy, rubric, inventory)
    assert cert["status"] == "V49_REVIEWER_ATTACK_SURFACE_HARDENING_FAILED"
    assert "readme_evidence_boundary_section_present" in cert["failed_checks"]


def test_v49_future_source_rows_fail_when_role_fields_are_missing() -> None:
    runner = _load(RUNNER, "v49_runner_schema")
    complete = {
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "spatial_proxy": False,
        "coordinate_derived": False,
        "internal_runtime": False,
        "allowed_for_prediction": True,
        "allowed_for_holdout": False,
        "allowed_for_claim": True,
        "leakage_risk": "low",
        "rationale": "Sequence-only input.",
    }
    assert runner.missing_source_role_fields(complete) == []
    incomplete = {"source_class": "pure_non_coordinate"}
    assert set(runner.missing_source_role_fields(incomplete)) == set(runner.REQUIRED_SOURCE_ROLE_FIELDS) - {"source_class"}
