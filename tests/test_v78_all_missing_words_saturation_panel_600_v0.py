from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
V78_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V78"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v78_uses_requested_explicit_composition() -> None:
    manifest = _read(V78_ROOT / "v78_all_missing_words_target_manifest.json")
    expected = {
        "COILED_COIL_POSITIVE": 100,
        "COILED_COIL_NEAR_NEGATIVE_ASSEMBLY": 60,
        "COILED_COIL_NEAR_NEGATIVE_GLOBULAR": 40,
        "REPEAT_SOLENOID_POSITIVE": 100,
        "REPEAT_NEAR_NEGATIVE_BETA": 60,
        "REPEAT_NEAR_NEGATIVE_MULTIDOMAIN": 40,
        "KNOTTED_OR_SLIPKNOT_POSITIVE": 50,
        "KNOT_NEAR_NEGATIVE_GLOBULAR": 40,
        "KNOT_NEAR_NEGATIVE_BETA_OR_REPEAT": 30,
        "V77_SIGNAL_TM_SENTINEL_REPLAY": 40,
        "V76_SECRETORY_DISULFIDE_SENTINEL_REPLAY": 40,
        "V75_MULTIDOMAIN_ASSEMBLY_MEMBRANE_METAL_DISORDER_SENTINEL_REPLAY": 40,
        "RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS": 40,
    }

    assert manifest["batch_id"] == "V78_ALL_MISSING_WORDS_SATURATION_PANEL_600"
    assert manifest["engine_version_used"] == "E72"
    assert manifest["target_count_selected"] == 680
    assert manifest["composition_rule"] == expected
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected
    assert manifest["uses_static_observable_thresholds"] is False


def test_v78_passes_zero_failed_accepted_and_claim_gate() -> None:
    cert = _read(V78_ROOT / "v78_all_missing_words_certificate.json")

    assert cert["status"] == "V78_ALL_MISSING_WORDS_SATURATION_PANEL_PASSED"
    assert cert["requested_nominal_target_count"] == 600
    assert cert["actual_target_count_from_requested_composition"] == 680
    assert cert["targets_total"] == 680
    assert cert["accepted_count"] == 640
    assert cert["accepted_supported"] == 640
    assert cert["clean_abstain"] == 40
    assert cert["clean_abstain_supported"] == 40
    assert cert["failed_accepted_count"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(1.0)
    assert cert["controls_passed"] is True
    assert cert["sentinel_regressions"] == 0
    assert cert["candidate_grammars_remaining"] == []
    assert cert["protein_esperanto_language_saturation_status"] == "current_known_missing_words_implemented"
    assert cert["uses_static_observable_thresholds"] is False
    assert cert["static_thresholds_removed_from_engine"] is True
    assert cert["withheld_context_leakage_detected"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["next_required_batch"] == "V79_BLIND_RCSB_DISCOVERY_500_NO_KNOWN_MISSING_WORD_QUEUE"


def test_v78_rows_pass_matched_controls_and_promoted_word_checks() -> None:
    scoring = _read(V78_ROOT / "v78_all_missing_words_scoring_report.json")
    rows = scoring["rows"]
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]

    assert len(rows) == 680
    assert len(accepted) == 640
    assert len(abstained) == 40
    assert all(row["accepted_supported"] for row in accepted)
    assert all(row["matched_control_dominance_passed"] for row in accepted)
    assert all(row["clean_abstain_supported"] for row in abstained)
    assert all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in rows)
    assert all(row["missing_word_candidate"] is None for row in rows)
    assert all(row["known_coiled_coil_word"] == "coiled_coil_register" for row in rows if row["panel_group"] == "COILED_COIL_POSITIVE")
    assert all(row["known_repeat_solenoid_word"] == "repeat_solenoid_topology" for row in rows if row["panel_group"] == "REPEAT_SOLENOID_POSITIVE")
    assert all(row["known_knotted_topology_word"] == "knotted_topology" for row in rows if row["panel_group"] == "KNOTTED_OR_SLIPKNOT_POSITIVE")
