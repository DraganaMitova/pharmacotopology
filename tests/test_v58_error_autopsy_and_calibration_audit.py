from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
V58_RUNNER = ROOT / "scripts" / "run_v58_real_sequence_time_blind_folding_replication_gate_v0.py"
AUDIT_RUNNER = ROOT / "scripts" / "run_v58_error_autopsy_and_calibration_audit_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def audit_outputs() -> dict[str, object]:
    v58 = _load(V58_RUNNER, "v58_runner_for_autopsy")
    v58.run_v58()
    audit = _load(AUDIT_RUNNER, "v58_autopsy_runner")
    paths = audit.run_v58_autopsy()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    table = json.loads(paths["table"].read_text(encoding="utf-8"))
    return {"paths": paths, "certificate": cert, "table": table}


def test_v58_autopsy_certificate_preserves_failures_without_engine_changes(audit_outputs: dict[str, object]) -> None:
    cert = audit_outputs["certificate"]
    assert isinstance(cert, dict)
    assert cert["status"] == "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_COMPLETED_REVIEW_REQUIRED"
    assert cert["target_count"] == 20
    assert cert["raw_supported_count"] == 0
    assert cert["raw_failure_count"] == 20
    assert cert["accepted_count"] == 0
    assert cert["strict_accepted_count"] == 0
    assert cert["abstain_recommended_count"] == 20
    assert cert["accepted_accuracy"] is None
    assert cert["strict_accepted_accuracy"] is None
    assert cert["abstention_rate"] == 1.0
    assert cert["overall_accuracy"] == 0.0
    assert cert["failure_preserved"] is True
    assert cert["could_failures_have_been_avoided_by_abstention"] is True
    assert cert["engine_biology_modified"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["atomistic_md_executed"] is False
    assert cert["readme_touched"] is False


def test_v58_autopsy_table_has_required_self_decision_fields(audit_outputs: dict[str, object]) -> None:
    table = audit_outputs["table"]
    assert isinstance(table, dict)
    rows = table["rows"]
    assert len(rows) == 20
    required_fields = {
        "target_id",
        "selected_regime",
        "expected_regime_postseal",
        "sequence_only_selected_regime",
        "confidence_signal",
        "pass_level_1",
        "pass_level_2",
        "pass_level_3",
        "failure_bucket",
        "would_abstention_have_been_correct",
        "post_autopsy_decision",
    }
    for row in rows:
        assert required_fields <= row.keys()
        assert row["post_autopsy_decision"] in {"accepted", "accepted_with_caution", "abstain_recommended"}
        assert row["confidence_signal"] in {
            "self_consistent",
            "contested_top_signal",
            "selected_against_stronger_internal_signal",
            "self_abstained",
        }


def test_v58_autopsy_abstains_only_on_preserved_raw_failures(audit_outputs: dict[str, object]) -> None:
    rows = audit_outputs["table"]["rows"]
    raw_failures = [row for row in rows if not row["pass_level_1"] or not row["pass_level_3"]]
    abstained = [row for row in rows if row["post_autopsy_decision"] == "abstain_recommended"]
    accepted = [row for row in rows if row["post_autopsy_decision"] in {"accepted", "accepted_with_caution"}]
    assert len(raw_failures) == 20
    assert {row["target_id"] for row in raw_failures} == {row["target_id"] for row in abstained}
    assert all(row["would_abstention_have_been_correct"] is True for row in abstained)
    assert accepted == []


def test_v58_autopsy_failure_buckets_match_observed_v58_misses(audit_outputs: dict[str, object]) -> None:
    cert = audit_outputs["certificate"]
    assert cert["failure_bucket_counts"] == {
        "globular_vs_interface_or_oligomer_ambiguity": 2,
        "insufficient_evidence_or_missing_context": 16,
        "right_regime_wrong_topology_or_observable": 2,
    }


def test_v58_autopsy_outputs_are_written_and_do_not_touch_readme(audit_outputs: dict[str, object]) -> None:
    paths = audit_outputs["paths"]
    for path in paths.values():
        assert path.exists()
    assert "README" not in "\n".join(str(path) for path in paths.values())
    report = paths["report"].read_text(encoding="utf-8")
    assert "Engine biology modified: `False`" in report
