from __future__ import annotations

from pharmacotopology.protein_esperanto_engine import (
    BackboneCoordinateEmitter,
    ConstraintToDistanceMapCompiler,
    E80_ENGINE_COMPONENTS,
    E80_ENGINE_REVISION,
    E80_REQUIRED_HARD_FAMILIES,
    FoldGeometryCompiler,
    FoldQualityEvaluator,
    PhysicalRelaxationExecutor,
    RealHoldoutCoordinateLoader,
    TopologyConstraintCompiler,
    UniversalSolutionUnlockFirewall,
    e80_engine_manifest,
    e80_fresh_target_resolver,
    e80_normalize_hard_family,
    external_blind_benchmark_export,
)


def test_e80_manifest_and_family_aliases_are_geometry_not_field_setting() -> None:
    manifest = e80_engine_manifest()

    assert manifest["kind"] == "E80_REAL_FOLD_GEOMETRY_AND_PHYSICAL_CALIBRATION_ENGINE_MANIFEST_v0"
    assert manifest["engine_revision"] == E80_ENGINE_REVISION
    assert manifest["components"] == E80_ENGINE_COMPONENTS
    assert manifest["required_hard_families"] == E80_REQUIRED_HARD_FAMILIES
    assert manifest["native_coordinate_truth_allowed_before_prediction"] is False
    assert manifest["native_contacts_allowed_before_prediction"] is False
    assert manifest["external_fold_model_prediction_input_allowed"] is False
    assert manifest["proxy_physical_execution_for_target_fold_claim_allowed"] is False
    assert manifest["protein_folding_solved_by_field_setting"] is False
    assert e80_normalize_hard_family("knotted_topology") == "knot_slipknot"
    assert e80_normalize_hard_family("coiled_coil_register") == "coiled_coil"
    assert e80_normalize_hard_family("repeat_solenoid_topology") == "repeat_solenoid"
    assert e80_normalize_hard_family("closed_beta_topology") == "beta_closure"


def test_e80_empty_evidence_blocks_universal_solution() -> None:
    fresh = e80_fresh_target_resolver(required_families=E80_REQUIRED_HARD_FAMILIES, candidate_targets=[])
    holdout = RealHoldoutCoordinateLoader().load(holdout_rows=[])
    physical = PhysicalRelaxationExecutor().execute(predicted_fold_packets=[], execution_rows=[])
    quality = FoldQualityEvaluator().evaluate(
        predicted_fold_packets=[],
        holdout_packet=holdout,
        physical_relaxation_packet=physical,
    )
    external = external_blind_benchmark_export(
        benchmark_rows=[],
        export_path="data/protein_esperanto_engine/V86/unit_empty_export.json",
    )
    firewall = UniversalSolutionUnlockFirewall().evaluate(
        fresh_resolution=fresh,
        holdout_packet=holdout,
        physical_relaxation_packet=physical,
        fold_quality_packet=quality,
        external_benchmark=external,
        token_only_acceptance_count=0,
        required_families=E80_REQUIRED_HARD_FAMILIES,
    )

    assert fresh["fresh_target_shortage"] is True
    assert holdout["real_fold_holdout_count"] == 0
    assert physical["real_or_validated_physical_execution_count"] == 0
    assert quality["target_fold_claim_count"] == 0
    assert firewall["universal_folding_solution_claim_allowed"] is False
    assert firewall["protein_folding_solved"] is False
    assert "fresh_target_shortage" in firewall["blocked_reasons"]
    assert "real_fold_holdout_count_zero" in firewall["blocked_reasons"]
    assert "real_or_validated_physical_execution_count_zero" in firewall["blocked_reasons"]
    assert "target_fold_claim_count_zero" in firewall["blocked_reasons"]


def _sequence_for_family(family: str) -> str:
    if family == "secretory_disulfide":
        return "ACDECGHIKCLMNCQRSTC"
    if family == "intrinsic_disorder_no_single_fold":
        return "GSGSGSDEKRGSGSDEKR"
    return "ACDEFGHIKLMNPQRSTVWY"


def test_e80_can_unlock_when_geometry_holdout_physics_and_external_export_all_support() -> None:
    candidates = []
    predicted_packets = []
    holdout_rows = []
    execution_rows = []
    benchmark_rows = []
    geometry = FoldGeometryCompiler()
    distance = ConstraintToDistanceMapCompiler()
    topology = TopologyConstraintCompiler()
    emitter = BackboneCoordinateEmitter()
    for index, family in enumerate(E80_REQUIRED_HARD_FAMILIES, start=1):
        target_id = f"E80_TARGET_{index}_{family}"
        sequence = _sequence_for_family(family)
        candidates.append({
            "target_id": target_id,
            "family": family,
            "sequence": sequence,
            "fresh_blind_target": True,
            "deterministic_variant": False,
            "nonredundant": True,
        })
        sentence = {
            "target_id": target_id,
            "hard_family": family,
            "protein_sentence": f"sealed sentence for {family}",
            "predicted_contact_pairs": [{"residue_i": 1, "residue_j": len(sequence)}],
        }
        fold = geometry.compile(
            protein_sentence_packet=sentence,
            operator_state_propagation_summary={"target_id": target_id, "coordinate_truth_used": False},
            hypothesized_interaction_language_map={"target_id": target_id},
            sequence=sequence,
            allowed_preseal_evidence={
                "coordinate_truth_used_before_prediction": False,
                "native_coordinates_used_before_prediction": False,
                "native_contacts_used_before_prediction": False,
                "native_topology_labels_used_before_prediction": False,
                "postseal_annotations_used_before_prediction": False,
                "structure_model_used_as_prediction_input": False,
                "alphafold_or_external_fold_model_used_as_prediction_input": False,
            },
        )
        distance_map = distance.compile(fold_constraint_packet=fold)
        topology_packet = topology.compile(fold_constraint_packet=fold)
        prediction = emitter.emit(
            sequence=sequence,
            fold_constraint_packet=fold,
            distance_map_packet=distance_map,
            topology_constraint_packet=topology_packet,
        )
        predicted_packets.append(prediction)
        holdout_hash = f"holdout_hash_{index}"
        holdout_rows.append({
            "target_id": target_id,
            "opened_after_prediction_hash": True,
            "used_before_prediction": False,
            "coordinate_holdout_available": True,
            "contact_map_available": True,
            "topology_holdout_available": True,
            "fold_family": family,
            "topology_class": prediction["predicted_topology"],
            "holdout_hash": holdout_hash,
            "native_contact_map": prediction["predicted_contact_map"],
            "selected_contact_supports": True,
            "selected_topology_supports": True,
            "selected_family_supports": True,
            "contact_order_support": True,
            "long_range_contact_enrichment_supported": True,
            "disulfide_tm_assembly_ligand_knot_repeat_support": True,
            "family_specific_correctness": True,
            "wrong_grammar_supports": False,
            "bag_of_words_supports": False,
            "masked_sentence_supports": False,
            "wrong_target_supports": False,
        })
        execution_rows.append({
            "target_id": target_id,
            "execution_backend": "validated_coarse_protocol",
            "validated_coarse_execution": True,
            "real_openmm_execution": False,
            "proxy_only": False,
            "energy_before": index,
            "energy_after": index,
            "constraint_violation_before": index,
            "constraint_violation_after": index,
        })
        benchmark_rows.append({
            "target_id": target_id,
            "exported": True,
            "sealed_prediction_hash": prediction["prediction_hash"],
            "postseal_holdout_hash": holdout_hash,
            "target_fold_claim_allowed": True,
            "coordinate_native_leakage": False,
        })
    fresh = e80_fresh_target_resolver(required_families=E80_REQUIRED_HARD_FAMILIES, candidate_targets=candidates)
    holdout = RealHoldoutCoordinateLoader().load(holdout_rows=holdout_rows)
    physical = PhysicalRelaxationExecutor().execute(
        predicted_fold_packets=predicted_packets,
        execution_rows=execution_rows,
    )
    quality = FoldQualityEvaluator().evaluate(
        predicted_fold_packets=predicted_packets,
        holdout_packet=holdout,
        physical_relaxation_packet=physical,
    )
    external = external_blind_benchmark_export(
        benchmark_rows=benchmark_rows,
        export_path="data/protein_esperanto_engine/V86/unit_positive_export.json",
    )
    firewall = UniversalSolutionUnlockFirewall().evaluate(
        fresh_resolution=fresh,
        holdout_packet=holdout,
        physical_relaxation_packet=physical,
        fold_quality_packet=quality,
        external_benchmark=external,
        token_only_acceptance_count=0,
        required_families=E80_REQUIRED_HARD_FAMILIES,
    )

    assert fresh["fresh_target_shortage"] is False
    assert holdout["real_fold_holdout_count"] == len(E80_REQUIRED_HARD_FAMILIES)
    assert physical["real_or_validated_physical_execution_count"] == len(E80_REQUIRED_HARD_FAMILIES)
    assert physical["proxy_physical_execution_used_for_claim"] is False
    assert quality["target_fold_claim_count"] == len(E80_REQUIRED_HARD_FAMILIES)
    assert firewall["every_required_hard_family_has_supported_target_fold_claim"] is True
    assert firewall["unsupported_fold_claims"] == 0
    assert firewall["unsupported_physical_claims"] == 0
    assert firewall["external_blind_benchmark_passed"] is True
    assert firewall["universal_folding_solution_claim_allowed"] is True
    assert firewall["protein_folding_solved"] is True
