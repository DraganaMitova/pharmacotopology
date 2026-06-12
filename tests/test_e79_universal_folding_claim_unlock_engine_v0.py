from __future__ import annotations

from pathlib import Path

from pharmacotopology.protein_esperanto_engine import (
    E79_ENGINE_REVISION,
    E79_HARD_REGIME_FAMILIES,
    E79_HARD_REGIME_FAMILY_ALIASES,
    E79_REQUIRED_SEALED_OUTPUTS,
    E79_UNLOCK_COMPONENTS,
    atomistic_or_validated_physical_executor,
    e79_normalize_hard_regime_family,
    e79_unlock_engine_manifest,
    e79_universal_claim_firewall,
    external_blind_benchmark_export,
    family_generalization_evaluator,
    fresh_target_resolver,
    real_fold_holdout_loader,
    target_fold_evaluator,
)


def test_e79_manifest_declares_unlock_components_without_field_setting() -> None:
    manifest = e79_unlock_engine_manifest()

    assert manifest["kind"] == "E79_UNIVERSAL_FOLDING_CLAIM_UNLOCK_ENGINE_MANIFEST_v0"
    assert manifest["engine_revision"] == E79_ENGINE_REVISION
    assert manifest["components"] == E79_UNLOCK_COMPONENTS
    assert manifest["hard_regime_families"] == E79_HARD_REGIME_FAMILIES
    assert manifest["hard_regime_family_aliases"] == E79_HARD_REGIME_FAMILY_ALIASES
    assert manifest["required_sealed_outputs"] == E79_REQUIRED_SEALED_OUTPUTS
    assert manifest["fresh_targets_required_for_universal_claim"] is True
    assert manifest["deterministic_variants_for_universal_claim_allowed"] is False
    assert manifest["proxy_physical_execution_for_claim_allowed"] is False
    assert manifest["universal_folding_solution_claim_allowed_by_field_setting"] is False
    assert manifest["protein_folding_solved_by_field_setting"] is False
    assert e79_normalize_hard_regime_family("knotted_topology") == "knot_slipknot"
    assert e79_normalize_hard_regime_family("coiled_coil_register") == "coiled_coil"
    assert e79_normalize_hard_regime_family("repeat_solenoid_topology") == "repeat"


def test_e79_firewall_blocks_empty_unlock_evidence(tmp_path: Path) -> None:
    fresh = fresh_target_resolver(required_families=E79_HARD_REGIME_FAMILIES, candidate_targets=[])
    holdouts = real_fold_holdout_loader(holdout_rows=[])
    physical = atomistic_or_validated_physical_executor(execution_rows=[])
    target_fold = target_fold_evaluator(
        fresh_resolution=fresh,
        sealed_predictions=[],
        fold_holdouts=holdouts,
        physical_execution=physical,
    )
    family = family_generalization_evaluator(
        required_families=E79_HARD_REGIME_FAMILIES,
        target_fold_evaluation=target_fold,
        sentinels_preserved=True,
        failed_accepted_count=0,
        unsupported_claims=0,
    )
    export = external_blind_benchmark_export(
        benchmark_rows=[],
        export_path=str(tmp_path / "external_blind_export.json"),
    )
    firewall = e79_universal_claim_firewall(
        fresh_resolution=fresh,
        physical_execution=physical,
        target_fold_evaluation=target_fold,
        family_generalization=family,
        external_benchmark=export,
        unresolved_classes=family["missing_families"],
    )

    assert fresh["fresh_target_shortage"] is True
    assert physical["proxy_physical_execution_used_for_claim"] is False
    assert target_fold["target_fold_claim_count"] == 0
    assert family["general_solution_candidate_claim_allowed"] is False
    assert export["external_blind_benchmark_exported"] is True
    assert export["external_blind_benchmark_passed"] is False
    assert firewall["claim_5_target_fold_supported"] is False
    assert firewall["general_solution_candidate_claim_allowed"] is False
    assert firewall["universal_folding_solution_claim_allowed"] is False
    assert firewall["protein_folding_solved"] is False
    assert firewall["unsupported_fold_claims"] == 0
    assert firewall["unsupported_physical_claims"] == 0
    assert "fresh_target_shortage" in firewall["blocked_reasons"]
    assert "no_real_or_validated_physical_execution" in firewall["blocked_reasons"]
    assert "target_fold_claim_count_zero" in firewall["blocked_reasons"]


def _raw_family_for_test(family: str) -> str:
    aliases = {
        "cofactor_metal": "metal_ligand_basin",
        "coiled_coil": "coiled_coil_register",
        "knot_slipknot": "knotted_topology",
        "repeat": "repeat_solenoid_topology",
    }
    return aliases.get(family, family)


def _positive_unlock_rows() -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    candidates = []
    predictions = []
    holdouts = []
    executions = []
    benchmark_rows = []
    for index, family in enumerate(E79_HARD_REGIME_FAMILIES, start=1):
        raw_family = _raw_family_for_test(family)
        target_id = f"E79_UNLOCK_{index}_{family}"
        candidates.append({
            "target_id": target_id,
            "family": raw_family,
            "fresh_blind_target": True,
            "deterministic_variant": False,
            "nonredundant": True,
        })
        prediction_hash = f"sealed_prediction_hash_{index}"
        holdout_hash = f"postseal_holdout_hash_{index}"
        predictions.append({
            "target_id": target_id,
            "family": raw_family,
            "protein_sentence": f"sealed sentence for {raw_family}",
            "topology_prediction": f"sealed topology for {raw_family}",
            "contact_order_prediction": f"sealed contact order for {raw_family}",
            "fold_family_prediction": raw_family,
            "physical_execution_plan": f"sealed physical plan for {raw_family}",
            "prediction_hash": prediction_hash,
            "preseal_allowed_sources": {
                "coordinate_truth_used_before_prediction": False,
                "native_contacts_used_before_prediction": False,
                "structure_model_used_as_prediction_input": False,
            },
        })
        holdouts.append({
            "target_id": target_id,
            "opened_after_prediction_hash": True,
            "used_before_prediction": False,
            "coordinate_holdout_available": True,
            "native_contact_map_available": True,
            "topology_class_available": True,
            "assembly_state_available": True,
            "observable_family": raw_family,
            "source_hash": holdout_hash,
            "supports_selected_prediction": True,
            "wrong_grammar_supports": False,
            "bag_of_words_supports": False,
            "masked_sentence_supports": False,
            "wrong_target_supports": False,
        })
        executions.append({
            "target_id": target_id,
            "execution_backend": "validated_coarse_protocol",
            "real_execution_performed": True,
            "validated_coarse_or_atomistic": True,
            "proxy_physical_execution": False,
            "target_specific_topology_environment": True,
            "selected_sentence_execution": {"supports_selected": True},
            "wrong_sentence_execution": {"supports_selected": False},
            "bag_of_words_execution": {"supports_selected": False},
            "masked_sentence_execution": {"supports_selected": False},
            "wrong_target_observable_execution": {"supports_selected": False},
            "native_truth_used_before_execution": False,
        })
        benchmark_rows.append({
            "target_id": target_id,
            "exported": True,
            "sealed_prediction_hash": prediction_hash,
            "postseal_holdout_hash": holdout_hash,
            "target_fold_claim_allowed": True,
            "coordinate_native_leakage": False,
        })
    return candidates, predictions, holdouts, executions, benchmark_rows


def test_e79_can_unlock_only_when_fresh_families_holdouts_physics_and_external_export_all_pass() -> None:
    candidates, predictions, holdout_rows, execution_rows, benchmark_rows = _positive_unlock_rows()
    fresh = fresh_target_resolver(required_families=E79_HARD_REGIME_FAMILIES, candidate_targets=candidates)
    holdouts = real_fold_holdout_loader(holdout_rows=holdout_rows)
    physical = atomistic_or_validated_physical_executor(execution_rows=execution_rows)
    target_fold = target_fold_evaluator(
        fresh_resolution=fresh,
        sealed_predictions=predictions,
        fold_holdouts=holdouts,
        physical_execution=physical,
    )
    family = family_generalization_evaluator(
        required_families=E79_HARD_REGIME_FAMILIES,
        target_fold_evaluation=target_fold,
        sentinels_preserved=True,
        failed_accepted_count=0,
        unsupported_claims=target_fold["unsupported_fold_claims"] + physical["unsupported_physical_claims"],
    )
    export = external_blind_benchmark_export(
        benchmark_rows=benchmark_rows,
        export_path="data/protein_esperanto_engine/V85/unit_external_export.json",
    )
    firewall = e79_universal_claim_firewall(
        fresh_resolution=fresh,
        physical_execution=physical,
        target_fold_evaluation=target_fold,
        family_generalization=family,
        external_benchmark=export,
        unresolved_classes=[],
    )

    assert fresh["fresh_target_shortage"] is False
    assert fresh["missing_required_families"] == []
    assert physical["real_or_validated_execution_count"] == len(E79_HARD_REGIME_FAMILIES)
    assert physical["proxy_physical_execution_used_for_claim"] is False
    assert target_fold["target_fold_claim_count"] == len(E79_HARD_REGIME_FAMILIES)
    assert family["all_major_grammar_families_represented"] is True
    assert family["general_solution_candidate_claim_allowed"] is True
    assert export["external_blind_benchmark_passed"] is True
    assert firewall["blocked_reasons"] == []
    assert firewall["claim_5_target_fold_supported"] is True
    assert firewall["general_solution_candidate_claim_allowed"] is True
    assert firewall["universal_folding_solution_claim_allowed"] is True
    assert firewall["protein_folding_solved"] is True
