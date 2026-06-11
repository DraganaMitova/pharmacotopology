from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v67_builds_requested_mixed_fast_discovery_shard() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V67"
        / "v67_mixed_fast_discovery_target_manifest.json"
    )
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
        / "v67_e64_mixed_fast_discovery_200_certificate.json"
    )
    assert manifest["batch_id"] == "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
    assert manifest["engine_version_used"] == "E64"
    assert manifest["baseline_engine_version"] == "E63"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == {
        "NEW_RCSB_NONREDUNDANT_E64_DISCOVERY": 100,
        "V66_FAILED_ACCEPTED_REPLAY": 70,
        "SENTINEL_STABILITY_REPLAY": 30,
    }
    assert cert["batch_id"] == "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
    assert cert["batch_mode"] == "mixed_fast_discovery_shard"
    assert cert["engine_version_used"] == "E64"
    assert cert["targets_total"] == 200
    assert cert["controls_passed"] is True
    assert cert["passed_control_count"] == cert["control_count"] == 18


def test_v67_records_e64_lineage_revision() -> None:
    e64_cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "E64"
        / "e64_v66_tm_complex_lineage_revision_certificate.json"
    )
    assert e64_cert["engine_version"] == "E64"
    assert e64_cert["lineage"] == ["E60", "E61", "E62", "E63", "E64"]
    assert e64_cert["source_batch_trigger"] == "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200"
    assert e64_cert["next_required_batch"] == "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
    assert e64_cert["claim_allowed"] is False


def test_v67_answers_old_failure_and_sentinel_questions() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V67"
        / "v67_mixed_fast_discovery_certificate.json"
    )
    repair = cert["old_v66_failure_repair"]
    assert repair == {
        "old_v66_failure_count": 70,
        "old_v66_failures_repaired_by_e64": 1,
        "old_v66_failures_remaining": 69,
        "question_answer": "E64 repaired 1 of the 70 V66 failed-accepted targets; 69 remain for failure grammar mining.",
    }
    assert cert["sentinel_regression_count"] == 0
    assert cert["sentinel_regressions"] == []
    assert cert["v67_questions_answered"] == {
        "did_e64_repair_old_v66_failures": repair["question_answer"],
        "new_dominant_failure_mode": "assembly_required_core_vs_membrane_topology",
        "did_e64_create_regressions": "0 sentinel regressions.",
        "is_engine_over_accepting_or_abstaining": "over_accepting_relative_to_abstention",
    }


def test_v67_mines_failure_grammar_without_hiding_autopsy_fields() -> None:
    failure_report = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V67"
        / "v67_mixed_fast_discovery_failure_report.json"
    )
    assert failure_report["dominant_failure_mode"] == "assembly_required_core_vs_membrane_topology"
    assert failure_report["dominant_failure_count"] == 48
    assert failure_report["failure_modes"] == {
        "assembly_required_core_vs_membrane_topology": 48,
        "cofactor_locked_basin_vs_membrane_topology": 33,
        "short_interface_context_missing": 12,
        "disorder_or_low_complexity_context": 12,
        "assembly_required_core": 1,
        "over_abstain_membrane_multidomain_folding_proteostasis": 1,
    }
    rows = failure_report["failure_grammar_rows"]
    assert len(rows) == 107
    assert all(row["engine_thought"] for row in rows)
    assert all(row["reality_showed"] for row in rows)
    assert all(row["missing_esperanto_word"] for row in rows)
    assert Counter(row["failure_mode"] for row in rows).most_common(1) == [
        ("assembly_required_core_vs_membrane_topology", 48)
    ]


def test_v67_metrics_show_over_acceptance_mining_posture() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
        / "v67_e64_mixed_fast_discovery_200_certificate.json"
    )
    assert cert["status"] == "V67_E64_MIXED_FAST_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED"
    assert cert["supported_count"] == 93
    assert cert["failed_accepted_count"] == 106
    assert cert["abstain_count"] == 1
    assert cert["accepted_count"] == 199
    assert cert["accepted_accuracy"] == pytest.approx(0.46733668341708545)
    assert cert["raw_accuracy"] == pytest.approx(0.465)
    assert cert["coverage"] == pytest.approx(0.995)
    assert cert["accept_abstain_posture"] == "over_accepting_relative_to_abstention"
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False
