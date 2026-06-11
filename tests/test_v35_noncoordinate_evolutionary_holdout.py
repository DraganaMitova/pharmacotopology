from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v35_noncoordinate_evolutionary_holdout_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v35_holdout", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v34_passed() -> dict:
    return {
        "control_status": "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED",
        "claim_allowed": False,
        "new_MD_allowed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "next_action": "move_to_non_coordinate_external_evolutionary_or_blind_holdout_tests",
    }


def _valid_row(file_path: str) -> dict:
    return {
        "target": "KcsA",
        "file_path": file_path,
        "evidence_type": "sequence_conservation_filter_signature",
        "source_name": "External potassium channel family conservation source",
        "source_url_or_citation": "UniProt/Pfam-style external family conservation source, accessioned local holdout fixture",
        "source_date_or_version": "fixture_2026-06-11",
        "source_boundary": "external_noncoordinate_sequence_family_evidence_no_pdb_coordinates",
        "allowed_use": "claim_disabled_V35_noncoordinate_holdout_gate_only",
        "provenance_notes": "Sequence-family conservation context only; no PDB coordinates, no native metrics, no MD.",
        "coordinate_derived": False,
        "native_coordinate_used_before_selection": False,
        "claim_allowed": False,
    }


def _valid_content() -> str:
    return "\n".join([
        "target,family_context,filter_signature,ion_specificity,evolutionary_basis,notes",
        "KcsA,KcsA-like potassium channel family,TVGYG,K+,external MSA conservation across homolog family,sequence-only evidence no coordinates",
        "",
    ])


def _configure_tmp_repo(mod, tmp_path: Path) -> None:
    mod.REPO_ROOT = tmp_path
    mod.RUN_ROOT = tmp_path / "first_contact_clean_pharmacotopology_layer_run"
    mod.KCSA_NONCOORDINATE_ROOT = tmp_path / "data" / "external_constraints" / "KcsA" / "noncoordinate_evolutionary"


def test_empty_manifest_clean_abstains_with_claims_closed(tmp_path: Path) -> None:
    mod = _load_module()
    _configure_tmp_repo(mod, tmp_path)
    cert = mod.build_v35({"rows": []}, _v34_passed())
    assert cert["control_status"] == "V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE"
    assert cert["evidence_row_count"] == 0
    assert cert["noncoordinate_external_source_count"] == 0
    assert cert["operator_candidate_found"] is False
    assert cert["passed_control_count"] == 9
    assert cert["control_count"] == 9
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["positive_folding_evidence_found"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["fixed_threshold_policy"] == "forbidden"
    assert cert["native_metrics_used_for_selection"] is False


def test_v35_required_negative_controls_all_pass() -> None:
    mod = _load_module()
    controls = mod.run_v35_controls()
    assert len(controls) == 9
    assert [row["passed"] for row in controls] == [True] * 9
    assert {row["control_id"] for row in controls} == {
        "missing_noncoordinate_source_clean_abstains",
        "coordinate_derived_kcsa_csv_is_blocked",
        "internal_runtime_source_is_blocked",
        "placeholder_source_name_or_citation_is_invalid",
        "annotation_only_source_does_not_open_v35",
        "wrong_target_in_kcsa_v35_path_does_not_open",
        "renamed_potassium_filter_signature_fails_content_validation",
        "ion_specificity_relabel_away_from_potassium_fails",
        "generic_channel_annotation_is_not_enough",
    }


def test_coordinate_derived_source_is_blocked(tmp_path: Path) -> None:
    mod = _load_module()
    _configure_tmp_repo(mod, tmp_path)
    row = _valid_row("data/external_constraints/KcsA/pore_filter/kcsa_1bl8_pore_filter_external_contacts.csv")
    row["evidence_type"] = "coordinate_derived_contact"
    row["source_name"] = "RCSB PDB 1BL8 coordinate contact source"
    row["source_url_or_citation"] = "RCSB PDB 1BL8 coordinate-derived contacts"
    row["coordinate_derived"] = True
    cert = mod.build_v35({"rows": [row]}, _v34_passed())
    assert cert["control_status"] == "V35_BLOCKED_COORDINATE_DERIVED_SOURCE_SUPPLIED"
    assert cert["coordinate_derived_source_count"] == 1
    assert cert["noncoordinate_external_source_count"] == 0
    assert cert["claim_allowed"] is False


def test_internal_runtime_source_is_blocked(tmp_path: Path) -> None:
    mod = _load_module()
    _configure_tmp_repo(mod, tmp_path)
    row = _valid_row("first_contact_clean_pharmacotopology_layer_run/V33_NEGATIVE_CONTROLS/internal_poison_source.csv")
    row["evidence_type"] = "internal_runtime_report"
    row["source_name"] = "Internal runtime report"
    row["source_url_or_citation"] = "first_contact_clean_pharmacotopology_layer_run generated report"
    cert = mod.build_v35({"rows": [row]}, _v34_passed())
    assert cert["control_status"] == "V35_BLOCKED_INTERNAL_RUNTIME_SOURCE_SUPPLIED"
    assert cert["internal_runtime_source_count"] == 1
    assert cert["noncoordinate_external_source_count"] == 0


def test_invalid_rows_do_not_open_v35(tmp_path: Path) -> None:
    mod = _load_module()
    _configure_tmp_repo(mod, tmp_path)
    path = "data/external_constraints/KcsA/noncoordinate_evolutionary/source.csv"
    row = _valid_row(path)
    row["source_name"] = "TODO placeholder"
    row["source_url_or_citation"] = "citation needed"
    out = mod.evaluate_manifest_rows([row], {path: _valid_content()})
    assert out["source_status"] == "V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE"
    assert "row_0:placeholder_or_missing_source_name" in out["failed_checks"]
    assert "row_0:placeholder_or_missing_source_url_or_citation" in out["failed_checks"]

    annotation = _valid_row(path)
    annotation["evidence_type"] = "annotation_only_claim"
    out = mod.evaluate_manifest_rows([annotation], {path: _valid_content()})
    assert out["source_status"] == "V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE"
    assert "row_0:excluded_evidence_type:annotation_only_claim" in out["failed_checks"]

    wrong_target = _valid_row(path)
    wrong_target["target"] = "XCL1_lymphotactin"
    out = mod.evaluate_manifest_rows([wrong_target], {path: _valid_content()})
    assert out["source_status"] == "V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE"
    assert "row_0:target_not_KcsA" in out["failed_checks"]


def test_content_validation_rejects_damaged_signatures() -> None:
    mod = _load_module()
    no_signature = mod.validate_noncoordinate_content(_valid_content().replace("TVGYG", "GYAVT"), "sequence_conservation_filter_signature")
    assert no_signature["valid"] is False
    assert "missing_potassium_channel_filter_signature" in no_signature["failures"]

    no_potassium = mod.validate_noncoordinate_content(
        _valid_content().replace("K+", "Na+").replace("potassium channel", "sodium channel"),
        "sequence_conservation_filter_signature",
    )
    assert no_potassium["valid"] is False
    assert "missing_potassium_ion_specificity" in no_potassium["failures"]

    generic = mod.validate_noncoordinate_content("KcsA channel protein generic annotation only", "sequence_conservation_filter_signature")
    assert generic["valid"] is False
    assert "missing_potassium_channel_filter_signature" in generic["failures"]
    assert "missing_external_evolutionary_or_family_basis" in generic["failures"]


def test_positive_noncoordinate_source_passes_claim_disabled(tmp_path: Path) -> None:
    mod = _load_module()
    _configure_tmp_repo(mod, tmp_path)
    source = tmp_path / "data" / "external_constraints" / "KcsA" / "noncoordinate_evolutionary" / "kcsa_family_signature.csv"
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text(_valid_content(), encoding="utf-8")
    rel = "data/external_constraints/KcsA/noncoordinate_evolutionary/kcsa_family_signature.csv"
    cert = mod.build_v35({"rows": [_valid_row(rel)]}, _v34_passed())
    assert cert["source_status"] == "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT_READY_CLAIM_DISABLED"
    assert cert["control_status"] == "V35_NONCOORDINATE_EVOLUTIONARY_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED"
    assert cert["evidence_row_count"] == 1
    assert cert["noncoordinate_external_source_count"] == 1
    assert cert["coordinate_derived_source_count"] == 0
    assert cert["internal_runtime_source_count"] == 0
    assert cert["operator_candidate_found"] is True
    assert cert["selected_source_paths"] == [rel]
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_writer_outputs_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    _configure_tmp_repo(mod, tmp_path)
    cert = mod.build_v35({"rows": []}, _v34_passed())
    paths = mod.write_outputs(tmp_path / "out", cert)
    for path in paths.values():
        assert path.exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["next_decision"]["claim_allowed"] is False
    assert written["artifacts"]["certificate"] == str(paths["certificate"])
