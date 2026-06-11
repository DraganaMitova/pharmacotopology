from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def v61_outputs() -> dict[str, object]:
    paths = {
        "raw_candidate_cache": ROOT / "data" / "protein_esperanto_engine" / "V61" / "intake" / "raw_rcsb_30pct_representative_entities.json",
        "target_manifest": ROOT / "data" / "protein_esperanto_engine" / "V61" / "v61_rcsb_nonredundant_100_target_manifest.json",
        "engine_declaration": ROOT / "data" / "protein_esperanto_engine" / "V61" / "v61_rcsb_nonredundant_100_engine_declaration.json",
        "scoring_report": ROOT / "data" / "protein_esperanto_engine" / "V61" / "v61_rcsb_nonredundant_100_scoring_report.json",
        "failure_report": ROOT / "data" / "protein_esperanto_engine" / "V61" / "v61_rcsb_nonredundant_100_failure_report.json",
        "certificate": ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V61_RCSB_NONREDUNDANT_100_BATCH" / "v61_rcsb_nonredundant_100_batch_certificate.json",
        "report": ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V61_RCSB_NONREDUNDANT_100_BATCH" / "V61_RCSB_NONREDUNDANT_100_BATCH_REPORT.md",
        "protein_universe_ledger": ROOT / "data" / "protein_esperanto_engine" / "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN" / "ledgers" / "protein_universe_ledger_v0.json",
        "engine_version_ledger": ROOT / "data" / "protein_esperanto_engine" / "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN" / "ledgers" / "engine_version_ledger_v0.json",
        "failure_grammar_ledger": ROOT / "data" / "protein_esperanto_engine" / "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN" / "ledgers" / "failure_grammar_ledger_v0.json",
        "claim_ledger": ROOT / "data" / "protein_esperanto_engine" / "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN" / "ledgers" / "claim_ledger_v0.json",
    }
    for path in paths.values():
        assert path.exists()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["target_manifest"].read_text(encoding="utf-8"))
    scoring = json.loads(paths["scoring_report"].read_text(encoding="utf-8"))
    failure_ledger = json.loads(paths["failure_grammar_ledger"].read_text(encoding="utf-8"))
    return {
        "paths": paths,
        "certificate": cert,
        "manifest": manifest,
        "scoring": scoring,
        "failure_ledger": failure_ledger,
    }


def test_v61_runs_100_target_nonredundant_batch(v61_outputs: dict[str, object]) -> None:
    cert = v61_outputs["certificate"]
    assert isinstance(cert, dict)
    assert cert["batch_id"] == "V61_RCSB_NONREDUNDANT_100_BATCH"
    assert cert["targets_total"] == 100
    assert cert["target_selection_manual"] is False
    assert cert["accepted_count"] + cert["abstain_count"] == 100
    assert cert["supported_count"] <= cert["accepted_count"]
    assert cert["failed_accepted_count"] == cert["accepted_count"] - cert["supported_count"]
    assert cert["coverage"] == cert["accepted_count"] / 100
    assert cert["raw_accuracy"] == cert["supported_count"] / 100
    assert cert["engine_modified_during_batch"] is False
    assert cert["readme_modified"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["contact_truth_used_before_seal"] is False
    assert cert["alphafold_used_before_seal"] is False
    assert cert["controls_passed"] is True


def test_v61_manifest_is_automatic_rcsb_30pct_representatives(v61_outputs: dict[str, object]) -> None:
    manifest = v61_outputs["manifest"]
    assert isinstance(manifest, dict)
    assert manifest["kind"] == "V61_RCSB_NONREDUNDANT_100_TARGET_MANIFEST_v0"
    assert manifest["target_selection_manual"] is False
    assert manifest["target_count_selected"] == 100
    assert manifest["sequence_cluster_representative_selection"] is True
    assert manifest["sequence_cluster_identity_cutoff"] == 30
    rows = manifest["selected_targets"]
    clusters = [row["sequence_cluster_30_id"] for row in rows]
    assert len(set(clusters)) == 100
    assert all(row["source_database"] == "RCSB_PDB" for row in rows)
    assert all(row["structure_determination_methodology"].lower() == "experimental" for row in rows)
    assert all("protein" in row["polymer_type"].lower() for row in rows)
    assert all(40 <= row["sequence_length"] <= 800 for row in rows)


def test_v61_frozen_engine_ledger_starts_from_d927781(v61_outputs: dict[str, object]) -> None:
    paths = v61_outputs["paths"]
    assert isinstance(paths, dict)
    engine = json.loads(paths["engine_declaration"].read_text(encoding="utf-8"))
    ledger = json.loads(paths["engine_version_ledger"].read_text(encoding="utf-8"))
    assert engine["engine_version_used"] == "E60"
    assert engine["engine_start_commit_required"] == "d927781"
    assert engine["engine_source_last_commit"].startswith("d927781")
    assert engine["engine_modified_during_batch"] is False
    assert ledger["versions"][0]["engine_modified_during_batch"] is False
    assert ledger["versions"][0]["commit"].startswith("d927781")


def test_v61_protein_universe_ledger_has_required_schema(v61_outputs: dict[str, object]) -> None:
    paths = v61_outputs["paths"]
    assert isinstance(paths, dict)
    ledger = json.loads(paths["protein_universe_ledger"].read_text(encoding="utf-8"))
    assert ledger["kind"] == "V61_PROTEIN_UNIVERSE_LEDGER_v0"
    assert ledger["row_count"] == 100
    required = {
        "protein_id",
        "sequence",
        "source_database",
        "length",
        "organism",
        "experimental_structure_available",
        "fold_class_available",
        "disorder_data_available",
        "membrane_data_available",
        "kinetics_data_available",
        "process_data_available",
        "eligible_for_prediction",
        "eligible_for_holdout",
        "batch_id",
    }
    for row in ledger["rows"]:
        assert required <= row.keys()
        assert row["batch_id"] == "V61_RCSB_NONREDUNDANT_100_BATCH"
        assert row["experimental_structure_available"] is True
        assert row["eligible_for_prediction"] is True
        assert row["eligible_for_holdout"] is True


def test_v61_failures_are_reported_as_missing_esperanto(v61_outputs: dict[str, object]) -> None:
    cert = v61_outputs["certificate"]
    failure_ledger = v61_outputs["failure_ledger"]
    assert isinstance(cert, dict)
    assert isinstance(failure_ledger, dict)
    expected_failure_rows = cert["failed_accepted_count"] + cert["abstain_count"]
    assert failure_ledger["failure_count"] == expected_failure_rows
    assert sum(cert["failure_modes"].values()) == expected_failure_rows
    assert len(cert["missing_esperanto_candidates"]) == expected_failure_rows
    for row in failure_ledger["rows"]:
        assert row["failure_type"]
        assert row["missing_esperanto_rule"]
        assert "later engine version" in row["control_that_prevents_overfitting"]


def test_v61_seals_every_packet_before_postseal_holdout(v61_outputs: dict[str, object]) -> None:
    scoring = v61_outputs["scoring"]
    assert isinstance(scoring, dict)
    rows = scoring["rows"]
    assert len(rows) == 100
    packet_root = ROOT / "data" / "protein_esperanto_engine" / "V61" / "sealed_predictions"
    holdout_root = ROOT / "data" / "protein_esperanto_engine" / "V61" / "holdouts_postseal"
    packets = list(packet_root.glob("*/sealed_simulation_packet.json"))
    assert len(packets) == 100
    for packet_path in packets:
        packet = json.loads(packet_path.read_text(encoding="utf-8"))
        holdout = json.loads((holdout_root / packet["target_id"] / "postseal_holdout_manifest.json").read_text(encoding="utf-8"))
        assert packet["sealed_before_holdout"] is True
        assert packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"] == 0
        assert packet["evidence_manifest"]["internal_runtime_source_count_for_prediction"] == 0
        assert holdout["holdout_opened_after_prediction_hash"] == packet["prediction_hash"]
        assert holdout["postseal_sources"][0]["used_before_prediction"] is False


def test_v61_controls_include_readme_and_leakage_guards(v61_outputs: dict[str, object]) -> None:
    cert = v61_outputs["certificate"]
    assert isinstance(cert, dict)
    controls = {row["control_id"]: row for row in cert["controls"]}
    for control_id in [
        "target_selection_manual_false",
        "targets_total_100",
        "rcsb_experimental_protein_entities_only",
        "sequence_cluster_representative_selection",
        "engine_starts_as_d927781",
        "engine_modified_during_batch_false",
        "all_sealed_before_holdout",
        "coordinate_leakage_control",
        "alphafold_leakage_control",
        "holdout_opened_before_seal_control",
        "internal_runtime_leakage_control",
        "shuffled_sequence_controls_reported",
        "readme_modified_false",
    ]:
        assert controls[control_id]["passed"] is True
