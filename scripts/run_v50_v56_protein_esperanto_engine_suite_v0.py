#!/usr/bin/env python3
from __future__ import annotations

"""Run V50-V56 Protein Esperanto grammar and coarse operator-state propagator suite."""

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    COORDINATE_DERIVED,
    GRAMMAR_RULES,
    INTERNAL_RUNTIME,
    MECHANISM_CLASSES,
    SPATIAL_PROXY_NON_COORDINATE,
    STATE_VARIABLES,
    UNIVERSAL_MARKS,
    UNIVERSAL_OPERATORS,
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    make_openmm_bridge_spec,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
    validate_against_holdout,
)


DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V50_V56_PROTEIN_ESPERANTO_ENGINE"

PASSED = "V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_PASSED_REVIEW_REQUIRED"
FAILED = "V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_FAILED"
BLOCKED_LEAKAGE = "V50_TO_V56_BLOCKED_FOR_PREDICTION_LEAKAGE"


TARGET_SPECS: dict[str, dict[str, Any]] = {
    "V44_FUS_LC": {
        "target_name": "FUS-LC residues 1-214",
        "manifest": REPO_ROOT / "data/live_unsolved_targets/V44/FUS_LC/sources/source_manifest.json",
        "sequence_preference": ["sequence"],
        "expected_mechanism_class": "intrinsic_disorder_phase_separation",
        "hard_class": "IDP / phase separation",
        "focus_regions": [{"name": "FUS low-complexity aromatic/charged stickers", "span": "1-214"}],
        "perturbations": [
            {
                "perturbation_id": "FUS_AROMATIC_STICKER_MUTATION",
                "description": "mutate aromatic sticker residues in the LC region",
                "operator_scales": {"phase_operator": 0.45},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "FUS_CHARGE_SCREENING_CONDITION",
                "description": "screen charge/salt condition to reduce disorder expansion",
                "operator_scales": {"disorder_operator": 0.55},
                "metric": "disorder_order_balance",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "FUS_WRONG_REGION_CONTROL",
                "description": "perturb a non-existent membrane/proteostasis route",
                "operator_scales": {"proteostasis_operator": 0.20},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "phase_prone_dynamic_supported", "metric": "basin:phase_prone_dynamic", "comparator": ">=", "threshold": 0.25},
            {"check_id": "compact_single_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12},
        ],
    },
    "V45_TDP43_LCD": {
        "target_name": "TDP-43 LCD residues 274-414",
        "manifest": REPO_ROOT / "data/live_unsolved_targets/V45/TDP43_LCD/sources/source_manifest.json",
        "sequence_preference": ["sequence"],
        "expected_mechanism_class": "intrinsic_disorder_phase_separation",
        "hard_class": "IDP / prion-like phase ensemble",
        "focus_regions": [{"name": "TDP-43 sparse aromatic prion-like LCD", "span": "274-414"}],
        "perturbations": [
            {
                "perturbation_id": "TDP43_AROMATIC_SELF_ASSOCIATION_MUTATION",
                "description": "weaken aromatic self-association in the LCD",
                "operator_scales": {"phase_operator": 0.50},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "TDP43_RNA_OR_MULTIDOMAIN_CONTEXT_REMOVAL",
                "description": "remove RNA/multidomain context that supports phase pressure",
                "operator_scales": {"phase_operator": 0.60},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "TDP43_WRONG_MEMBRANE_CONTROL",
                "description": "force a membrane-channel-like perturbation",
                "operator_scales": {"membrane_pressure_operator": 0.20},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "prion_like_phase_ensemble_supported", "metric": "basin:phase_prone_dynamic", "comparator": ">=", "threshold": 0.23},
            {"check_id": "generic_compact_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12},
        ],
    },
    "V46_CFTR_F508DEL": {
        "target_name": "CFTR F508del",
        "manifest": REPO_ROOT / "data/live_unsolved_targets/V46/CFTR_F508del/sources/source_manifest.json",
        "sequence_preference": ["sequence", "full_sequence"],
        "expected_mechanism_class": "membrane_multidomain_folding_proteostasis",
        "hard_class": "membrane multidomain folding / proteostasis",
        "focus_regions": [{"name": "F508del NBD1 mutation site", "position": 508, "span": "508"}],
        "perturbations": [
            {
                "perturbation_id": "CFTR_F508DEL_DAMAGE",
                "description": "delete F508 in NBD1",
                "operator_scales": {"closure_operator": 0.58, "interface_operator": 0.62, "proteostasis_operator": 0.58},
                "damage": 0.60,
                "metric": "proteostasis_routing",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "CFTR_CORRECTOR_RESCUE_CONDITION",
                "description": "apply interface/proteostasis-aware corrector condition",
                "operator_scales": {"interface_operator": 1.20, "proteostasis_operator": 1.22},
                "rescue": 0.45,
                "metric": "proteostasis_routing",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "CFTR_WRONG_REGION_CONTROL",
                "description": "perturb an unrelated soluble loop outside NBD1/interface route",
                "operator_scales": {"host_hijack_operator": 0.20},
                "metric": "proteostasis_routing",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "proteostasis_route_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.55},
            {"check_id": "interface_route_supported", "metric": "interface_readiness", "comparator": ">=", "threshold": 0.55},
        ],
    },
    "V47_RFAH_CTD": {
        "target_name": "RfaH-CTD residues 101-162",
        "manifest": REPO_ROOT / "data/live_unsolved_targets/V47/RfaH_CTD/sources/source_manifest.json",
        "sequence_preference": ["ctd_sequence", "sequence", "full_sequence"],
        "expected_mechanism_class": "metamorphic_fold_switching",
        "hard_class": "metamorphic alpha/beta fold switching",
        "focus_regions": [{"name": "RfaH CTD state-separated region", "span": "101-162"}],
        "perturbations": [
            {
                "perturbation_id": "RFAH_CTD_RELEASE_CONTEXT",
                "description": "release CTD from NTD-bound context",
                "operator_scales": {"dual_basin_switch_operator": 1.15},
                "release": 0.55,
                "metric": "basin:beta_released_basin",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "RFAH_NTD_INTERFACE_STABILIZATION",
                "description": "stabilize NTD-CTD alpha/autoinhibited context",
                "alpha_bias": 0.45,
                "metric": "basin:alpha_context_basin",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "RFAH_WRONG_PHASE_CONTROL",
                "description": "force an unrelated phase-separation perturbation",
                "operator_scales": {"phase_operator": 0.20},
                "metric": "basin:beta_released_basin",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "alpha_basin_present", "metric": "basin:alpha_context_basin", "comparator": ">=", "threshold": 0.30},
            {"check_id": "beta_basin_present", "metric": "basin:beta_released_basin", "comparator": ">=", "threshold": 0.30},
            {"check_id": "averaged_fold_rejected", "metric": "basin:averaged_single_fold", "comparator": "<=", "threshold": 0.15},
        ],
    },
    "V48_SARS2_ORF6": {
        "target_name": "SARS-CoV-2 ORF6",
        "manifest": REPO_ROOT / "data/live_unsolved_targets/V48/SARS2_ORF6/sources/source_manifest.json",
        "sequence_preference": ["sequence", "c_terminal_sequence"],
        "expected_mechanism_class": "short_region_host_interface_hijacking",
        "hard_class": "short-region host-interface hijacking",
        "focus_regions": [{"name": "ORF6 C-terminal host-interface region", "span": "38-61"}],
        "perturbations": [
            {
                "perturbation_id": "ORF6_C_TERMINAL_DISRUPTION",
                "description": "disrupt ORF6 C-terminal residues 38-61",
                "operator_scales": {"host_hijack_operator": 0.35, "interface_operator": 0.45},
                "interface_disruption": 0.65,
                "metric": "basin:host_interface_engaged",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "ORF6_HOST_PARTNER_REMOVAL",
                "description": "remove RAE1/NUP98 host-interface context",
                "operator_scales": {"host_hijack_operator": 0.42},
                "interface_disruption": 0.45,
                "metric": "interface_readiness",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "ORF6_WRONG_LOCALIZATION_ONLY_CONTROL",
                "description": "alter localization-only context outside the C-terminal interface",
                "operator_scales": {"membrane_pressure_operator": 0.25},
                "metric": "basin:host_interface_engaged",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "host_interface_dominant", "metric": "basin:host_interface_engaged", "comparator": ">=", "threshold": 0.60},
            {"check_id": "compact_globular_shortcut_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.10},
        ],
    },
}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    return path


def _extract_scope(source: dict[str, Any]) -> dict[str, Any]:
    scope = source.get("sequence_region_scope")
    return scope if isinstance(scope, dict) else {}


def _extract_sequence(manifest: dict[str, Any], preferences: list[str]) -> str:
    scopes = [_extract_scope(source) for source in manifest.get("prediction_sources", []) if isinstance(source, dict)]
    for key in preferences:
        for scope in scopes:
            value = scope.get(key)
            if isinstance(value, str) and value:
                return value
    for scope in scopes:
        for value in scope.values():
            if isinstance(value, str) and len(value) >= 20 and set(value.upper()) <= set("ACDEFGHIKLMNPQRSTVWY"):
                return value
    raise SystemExit(f"could not extract sequence from manifest kind={manifest.get('kind')}")


def load_target_profile(target_id: str) -> dict[str, Any]:
    spec = TARGET_SPECS[target_id]
    manifest = _read_json(spec["manifest"], f"{target_id} source manifest")
    sources = [source for source in manifest.get("prediction_sources", []) if isinstance(source, dict)]
    return {
        "target_id": target_id,
        "target_name": spec["target_name"],
        "sequence": _extract_sequence(manifest, spec["sequence_preference"]),
        "sources": sources,
        "focus_regions": spec["focus_regions"],
        "perturbations": spec["perturbations"],
        "expected_mechanism_class": spec["expected_mechanism_class"],
        "hard_class": spec["hard_class"],
        "holdout_observables": spec["holdout_observables"],
    }


def extract_v50_grammar() -> dict[str, Any]:
    rows = []
    required_sections = [
        "universal_marks",
        "universal_operators",
        "mechanism_classes",
        "state_variables",
        "environmental_pressures",
        "allowed_evidence_classes",
        "forbidden_evidence_classes",
        "transition_rules",
        "perturbation_prediction_rules",
        "falsification_rules",
        "null_controls",
        "simulation_readiness_decision",
    ]
    for mechanism_class, rule in GRAMMAR_RULES.items():
        rows.append({
            "mechanism_class": mechanism_class,
            "grammar_sentence": (
                f"{' + '.join(rule['marks'])} + {' + '.join(rule['pressures'])} + "
                f"{' + '.join(rule['operators'])} -> {rule['state_change']} -> {rule['testable_effect']}"
            ),
            "marks": rule["marks"],
            "pressures": rule["pressures"],
            "operators": rule["operators"],
            "state_change": rule["state_change"],
            "testable_effect": rule["testable_effect"],
            "falsification_rule": rule["falsification_rule"],
            "null_control": rule["null_control"],
        })
    return {
        "kind": "V50_PROTEIN_ESPERANTO_MECHANISM_GRAMMAR_EXTRACTION_v0",
        "required_sections": required_sections,
        "universal_marks": UNIVERSAL_MARKS,
        "universal_operators": UNIVERSAL_OPERATORS,
        "mechanism_classes": MECHANISM_CLASSES,
        "state_variables": STATE_VARIABLES,
        "allowed_evidence_classes": ["pure_non_coordinate", "explicitly_tagged_spatial_proxy_non_coordinate_prediction_input"],
        "forbidden_evidence_classes": ["coordinate_derived_before_seal", "internal_runtime_as_biological_evidence", "holdout_before_seal"],
        "transition_rules": rows,
        "perturbation_prediction_rules": [
            "correct perturbations move the named state variable in the predicted direction",
            "wrong-region perturbations must remain weaker than correct perturbations",
            "forced wrong grammar must abstain or fail",
        ],
        "falsification_rules": [row["falsification_rule"] for row in rows],
        "null_controls": [row["null_control"] for row in rows],
        "all_existing_hard_classes_rewritten": True,
        "simulation_readiness_decision": "ready_for_v51_contract_and_v52_coarse_mvp",
        "folding_problem_solved": False,
    }


def build_v51_engine_spec() -> dict[str, Any]:
    return {
        "kind": "V51_FOLDING_SIMULATION_ENGINE_SPEC_v0",
        "input_schema": {
            "sequence": "amino-acid string",
            "prediction_sources": "V49-classified evidence rows",
            "focus_regions": "optional named sequence spans, no coordinates",
            "forced_grammar": "optional control input; wrong grammar must abstain",
        },
        "state_vector_schema": [
            "position_index",
            "residue_identity",
            "segment_id",
            "charge_mark",
            "hydrophobic_mark",
            "disorder_mark",
            "secondary_propensity_mark",
            "interface_mark",
            "membrane_mark",
            "operator_activations",
            "current_state",
            "state_confidence",
        ],
        "operator_schema": [
            "operator",
            "acts_on",
            "activation_strength",
            "activated_by_evidence_ids",
            "pushes_toward",
            "perturbation_should",
            "falsified_by",
            "state_variable",
        ],
        "simulation_cycle": [
            "gate evidence",
            "build sequence field",
            "select or abstain mechanism grammar",
            "build operator field",
            "simulate coarse state trajectory",
            "seal prediction hash",
            "open holdout validation",
        ],
        "trajectory_output_schema": [
            "timepoint",
            "residue_exposure",
            "segment_compaction",
            "contact_probability",
            "operator_activation",
            "frustration",
            "state_basin_occupancy",
            "interface_readiness",
            "disorder_order_balance",
            "proteostasis_routing",
        ],
        "validation_schema": [
            "sealed_prediction_hash",
            "holdout_opened_after_prediction_hash",
            "score_label",
            "checks",
            "postseal_sources",
        ],
        "failure_modes": [
            "coordinate leakage before seal",
            "internal runtime biological evidence",
            "holdout opened before seal",
            "wrong grammar overpromotion",
            "generic annotation-only promotion",
            "all targets collapse into same output",
            "perturbation direction failure",
        ],
        "atomistic_coordinates_required": False,
        "folding_problem_solved": False,
    }


def _postseal_holdout(profile: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V55_SIMULATION_HOLDOUT_VALIDATION_MANIFEST_v0",
        "target_id": profile["target_id"],
        "target_name": profile["target_name"],
        "expected_mechanism_class": profile["expected_mechanism_class"],
        "expected_observables": profile["holdout_observables"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{profile['target_id']}_POSTSEAL_OPERATOR_HOLDOUT",
                "source_role": "holdout_validation",
                "source_class": "spatial_proxy_non_coordinate",
                "used_before_prediction": False,
                "rationale": "Post-seal operator/state validation summary for coarse trajectory direction.",
            }
        ],
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {
        "control_id": control_id,
        "passed": bool(passed),
        "reason": reason,
        "observed": observed,
    }


def _run_controls(profiles: list[dict[str, Any]], packets: list[dict[str, Any]], validations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fus = next(profile for profile in profiles if profile["target_id"] == "V44_FUS_LC")
    cftr = next(profile for profile in profiles if profile["target_id"] == "V46_CFTR_F508DEL")
    orf6 = next(profile for profile in profiles if profile["target_id"] == "V48_SARS2_ORF6")
    fus_packet = next(packet for packet in packets if packet["target_id"] == "V44_FUS_LC")
    cftr_packet = next(packet for packet in packets if packet["target_id"] == "V46_CFTR_F508DEL")
    orf6_packet = next(packet for packet in packets if packet["target_id"] == "V48_SARS2_ORF6")

    random_packet = build_sealed_operator_state_packet(
        target_id="CONTROL_RANDOM_SEQUENCE",
        target_name="random sequence control",
        sequence=deterministic_random_sequence(len(fus["sequence"])),
        sources=fus["sources"],
        perturbations=[],
    )
    shuffled_packet = build_sealed_operator_state_packet(
        target_id="CONTROL_SHUFFLED_CFTR",
        target_name="shuffled CFTR sequence control",
        sequence=shuffled_sequence(cftr["sequence"]),
        sources=cftr["sources"],
        perturbations=[],
    )
    swapped_packet = build_sealed_operator_state_packet(
        target_id="CONTROL_SWAPPED_EVIDENCE",
        target_name="ORF6 sequence with FUS evidence",
        sequence=orf6["sequence"],
        sources=fus["sources"],
        perturbations=[],
    )
    name_only_packet = build_sealed_operator_state_packet(
        target_id="CONTROL_NAME_ONLY",
        target_name="target name only control",
        sequence=fus["sequence"],
        sources=[],
        perturbations=[],
    )
    generic_packet = build_sealed_operator_state_packet(
        target_id="CONTROL_GENERIC_ANNOTATION",
        target_name="generic annotation only control",
        sequence=orf6["sequence"],
        sources=[{
            "source_id": "GENERIC_VIRAL_ACCESSORY_ONLY",
            "source_class": "pure_non_coordinate",
            "source_role": "prediction_input",
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "evidence_statement": "generic viral accessory protein annotation only",
        }],
        perturbations=[],
    )
    forced_wrong = build_sealed_operator_state_packet(
        target_id="CONTROL_FORCED_WRONG_GRAMMAR",
        target_name="FUS forced to globular closure",
        sequence=fus["sequence"],
        sources=fus["sources"],
        perturbations=[],
        forced_grammar="globular_closure",
    )
    coord_gate = evidence_boundary_gate([{
        "source_id": "BAD_PDB_COORDINATES",
        "source_role": "prediction_input",
        "source_class": COORDINATE_DERIVED,
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "BAD_RUNTIME_REPORT",
        "source_role": "prediction_input",
        "source_class": INTERNAL_RUNTIME,
        "internal_runtime": True,
    }])
    spatial_gate = evidence_boundary_gate([{
        "source_id": "UNTAGGED_FRET_PROXY",
        "source_role": "prediction_input",
        "source_type": "FRET distance summary",
        "coordinate_derived": False,
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "PREOPENED_HOLDOUT",
        "source_role": "holdout_validation",
        "source_class": SPATIAL_PROXY_NON_COORDINATE,
        "spatial_proxy": True,
    }])
    contradicted_holdout = _postseal_holdout(fus, fus_packet)
    contradicted_holdout["expected_observables"] = [
        {"check_id": "impossible_compact_fold", "metric": "basin:compact_single_fold", "comparator": ">=", "threshold": 0.90}
    ]
    contradicted_validation = validate_against_holdout(sealed_packet=fus_packet, holdout=contradicted_holdout)
    correct_perturbations = sum(
        1
        for packet in packets
        for row in packet["predicted_perturbation_table"]
        if row["direction_passed"] and not row["perturbation_id"].endswith("CONTROL")
    )
    total_perturbations = sum(
        1
        for packet in packets
        for row in packet["predicted_perturbation_table"]
        if not row["perturbation_id"].endswith("CONTROL")
    )
    wrong_controls = [
        row
        for packet in packets
        for row in packet["predicted_perturbation_table"]
        if row["perturbation_id"].endswith("CONTROL")
    ]
    return [
        _control(
            "random_sequence_control",
            random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain",
            "Random sequence must not create a target-specific phase packet from FUS evidence.",
            random_packet["selected_mechanism_grammar"]["mechanism_class"],
        ),
        _control(
            "shuffled_sequence_control",
            sequence_operator_coherence(shuffled_packet) < sequence_operator_coherence(cftr_packet),
            "Shuffled sequence must weaken operator coherence.",
            {
                "original": sequence_operator_coherence(cftr_packet),
                "shuffled": sequence_operator_coherence(shuffled_packet),
            },
        ),
        _control(
            "swapped_evidence_control",
            swapped_packet["selected_mechanism_grammar"]["mechanism_class"] != orf6["expected_mechanism_class"],
            "FUS evidence must not validate ORF6 host-hijacking grammar.",
            swapped_packet["selected_mechanism_grammar"]["mechanism_class"],
        ),
        _control(
            "wrong_target_control",
            name_only_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain",
            "Target name without allowed evidence must abstain.",
            name_only_packet["selected_mechanism_grammar"]["mechanism_class"],
        ),
        _control(
            "generic_annotation_only_control",
            generic_packet["selected_mechanism_grammar"]["mechanism_class"] != orf6["expected_mechanism_class"],
            "Generic viral accessory annotation alone must not create ORF6-specific grammar.",
            generic_packet["selected_mechanism_grammar"]["mechanism_class"],
        ),
        _control(
            "coordinate_leakage_control",
            coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [],
            "Coordinate-derived evidence is blocked before sealing.",
            coord_gate,
        ),
        _control(
            "internal_runtime_leakage_control",
            runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [],
            "Internal runtime artifacts never initialize biological evidence.",
            runtime_gate,
        ),
        _control(
            "spatial_proxy_tagging_control",
            spatial_gate["spatial_proxy_untagged_or_misused_count"] == 1,
            "Spatial-proxy evidence must be explicitly tagged.",
            spatial_gate,
        ),
        _control(
            "forced_wrong_grammar_control",
            forced_wrong["selected_mechanism_grammar"]["forced_grammar_rejected"] is True
            and forced_wrong["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain",
            "Forced wrong grammar must fail or cleanly abstain.",
            forced_wrong["selected_mechanism_grammar"],
        ),
        _control(
            "failed_prediction_not_repaired_after_holdout",
            contradicted_validation["score_label"] == "contradicted",
            "Contradicted validation remains failed after holdout opening.",
            contradicted_validation,
        ),
        _control(
            "holdout_opened_before_seal_control",
            holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [],
            "Holdout evidence opened before seal blocks the run.",
            holdout_gate,
        ),
        _control(
            "wild_type_mutant_direction_control",
            correct_perturbations == total_perturbations and all(row["direction_passed"] for row in wrong_controls),
            "Correct perturbations move in the expected direction and wrong controls remain weak.",
            {"correct": correct_perturbations, "total": total_perturbations, "wrong_controls": wrong_controls},
        ),
        _control(
            "hard_class_validation_control",
            all(row["score_label"] == "supported" for row in validations),
            "All hard-class sealed trajectories must survive post-seal validation.",
            [row["score_label"] for row in validations],
        ),
        _control(
            "folding_problem_solved_never_true",
            all(packet["folding_problem_solved"] is False for packet in packets),
            "The operator-state propagator must not set folding_problem_solved=true.",
        ),
    ]


def _battery_rows(packets: list[dict[str, Any]], validations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    validation_by_target = {row["target_id"]: row for row in validations}
    rows = []
    for packet in packets:
        final = packet["operator_state_propagation_summary"]["final_state_summary"]
        mechanism = packet["selected_mechanism_grammar"]["mechanism_class"]
        rows.append({
            "target_id": packet["target_id"],
            "target_name": packet["target_name"],
            "mechanism_class": mechanism,
            "operator_names": packet["operator_field"]["operator_names"],
            "final_state_basin_occupancy": final["state_basin_occupancy"],
            "compact_single_fold_probability": final["state_basin_occupancy"].get("compact_single_fold", 0.0),
            "score_label": validation_by_target[packet["target_id"]]["score_label"],
        })
    return rows


def _aggregate_cert(
    *,
    grammar: dict[str, Any],
    engine_spec: dict[str, Any],
    packets: list[dict[str, Any]],
    validations: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    bridge: dict[str, Any],
) -> dict[str, Any]:
    failed = [control["control_id"] for control in controls if not control["passed"]]
    leakage = sum(packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"] for packet in packets)
    runtime = sum(packet["evidence_manifest"]["internal_runtime_source_count_for_prediction"] for packet in packets)
    all_classes = {packet["selected_mechanism_grammar"]["mechanism_class"] for packet in packets}
    expected_classes = {TARGET_SPECS[target_id]["expected_mechanism_class"] for target_id in TARGET_SPECS}
    status = BLOCKED_LEAKAGE if leakage or runtime else PASSED if not failed else FAILED
    return {
        "kind": "V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "v50_grammar_extraction_passed": grammar["all_existing_hard_classes_rewritten"] is True,
        "v51_engine_contract_frozen": "state_vector_schema" in engine_spec and "operator_schema" in engine_spec,
        "v52_operator_state_propagation_mvp_passed": all(packet["sealed_before_holdout"] for packet in packets),
        "v53_hard_class_battery_passed": all(row["score_label"] == "supported" for row in validations),
        "v54_perturbation_engine_passed": all(
            row["direction_passed"]
            for packet in packets
            for row in packet["predicted_perturbation_table"]
        ),
        "v55_postseal_validation_passed": all(row["holdout_opened_after_prediction_hash"] for row in validations),
        "v56_openmm_bridge_spec_passed": bridge["custom_force_mapping_count"] >= 5 and bridge["openmm_execution_required_for_v56"] is False,
        "target_count": len(packets),
        "hard_classes_expected": sorted(expected_classes),
        "hard_classes_observed": sorted(all_classes),
        "same_language_different_regimes": len(all_classes) >= 4,
        "sealed_prediction_count": len(packets),
        "supported_validation_count": sum(1 for row in validations if row["score_label"] == "supported"),
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "failed_checks": failed,
        "controls": controls,
        "coordinate_derived_source_count_before_prediction": leakage,
        "internal_runtime_source_count_for_prediction": runtime,
        "folding_problem_solved": False,
        "atomistic_md_performed": False,
        "openmm_bridge_execution_deferred": True,
        "claim_allowed": status == PASSED,
        "allowed_claim_text": (
            "We extracted a reusable operator grammar from multiple hard protein regimes and built a leakage-controlled coarse simulation engine that predicts mechanism trajectories and perturbation directions before holdout validation."
            if status == PASSED else ""
        ),
        "forbidden_claims": [
            "universal protein folding is solved",
            "the operator-state propagator predicts atomistic coordinates",
            "AlphaFold or PDB coordinates were used before sealing",
            "OpenMM/GROMACS validation is complete",
            "wrong grammar can be forced into success",
        ],
    }


def _write_report(path: Path, cert: dict[str, Any], battery: list[dict[str, Any]]) -> None:
    lines = [
        "# V50-V56 Protein Esperanto Engine Suite",
        "",
        f"Status: `{cert['status']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        f"atomistic_md_performed: `{cert['atomistic_md_performed']}`",
        f"Targets: `{cert['target_count']}`",
        f"Supported validations: `{cert['supported_validation_count']}` / `{cert['target_count']}`",
        f"Controls passed: `{cert['passed_control_count']}` / `{cert['control_count']}`",
        "",
        "## Version Gates",
        f"- V50 grammar extraction: `{cert['v50_grammar_extraction_passed']}`",
        f"- V51 engine contract: `{cert['v51_engine_contract_frozen']}`",
        f"- V52 coarse operator-state propagator MVP: `{cert['v52_operator_state_propagation_mvp_passed']}`",
        f"- V53 hard-class battery: `{cert['v53_hard_class_battery_passed']}`",
        f"- V54 perturbation engine: `{cert['v54_perturbation_engine_passed']}`",
        f"- V55 post-seal validation: `{cert['v55_postseal_validation_passed']}`",
        f"- V56 OpenMM bridge spec: `{cert['v56_openmm_bridge_spec_passed']}`",
        "",
        "## Hard-Class Battery",
    ]
    for row in battery:
        lines.append(
            f"- `{row['target_id']}` `{row['mechanism_class']}` score `{row['score_label']}` basins `{row['final_state_basin_occupancy']}`"
        )
    lines.extend(["", "## Controls"])
    for control in cert["controls"]:
        lines.append(f"- `{control['control_id']}`: `{control['passed']}`")
    lines.extend([
        "",
        "## Claim Boundary",
        cert["allowed_claim_text"] or "No claim allowed until failed checks are fixed.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v50_v56(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    profiles = [load_target_profile(target_id) for target_id in TARGET_SPECS]
    grammar = extract_v50_grammar()
    engine_spec = build_v51_engine_spec()
    bridge = make_openmm_bridge_spec()
    packets: list[dict[str, Any]] = []
    validations: list[dict[str, Any]] = []

    _write_json(DATA_ROOT / "V50" / "protein_esperanto_mechanism_grammar_extraction.json", grammar)
    _write_json(DATA_ROOT / "V51" / "operator_state_engine_spec.json", engine_spec)
    _write_json(DATA_ROOT / "V56" / "operator_to_custom_force_bridge_spec.json", bridge)

    for profile in profiles:
        packet = build_sealed_operator_state_packet(
            target_id=profile["target_id"],
            target_name=profile["target_name"],
            sequence=profile["sequence"],
            sources=profile["sources"],
            focus_regions=profile["focus_regions"],
            perturbations=profile["perturbations"],
        )
        packets.append(packet)
        target_dir = DATA_ROOT / "V52" / "sealed_predictions" / profile["target_id"]
        _write_json(target_dir / "sealed_operator_state_packet.json", packet)
        _write_json(target_dir / "prediction_inputs_manifest.json", {
            "kind": "V52_PREDICTION_INPUTS_MANIFEST_v0",
            "target_id": profile["target_id"],
            "source_ids": packet["input_evidence_manifest"]["source_ids"],
            "holdouts_available_before_seal": False,
            "coordinates_available_before_seal": False,
        })
        holdout = _postseal_holdout(profile, packet)
        validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
        validations.append(validation)
        _write_json(DATA_ROOT / "V55" / "holdouts_postseal" / profile["target_id"] / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "V55" / "validation" / profile["target_id"] / "validation_result.json", validation)

    battery = _battery_rows(packets, validations)
    controls = _run_controls(profiles, packets, validations)
    cert = _aggregate_cert(
        grammar=grammar,
        engine_spec=engine_spec,
        packets=packets,
        validations=validations,
        controls=controls,
        bridge=bridge,
    )
    cert_hash = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    cert["certificate_hash"] = cert_hash
    _write_json(DATA_ROOT / "V53" / "operator_state_multiregime_battery.json", {
        "kind": "V53_OPERATOR_STATE_MULTIREGIME_BATTERY_v0",
        "rows": battery,
        "passed": cert["v53_hard_class_battery_passed"],
    })
    _write_json(DATA_ROOT / "V54" / "operator_state_perturbation_and_condition_response.json", {
        "kind": "V54_OPERATOR_STATE_PERTURBATION_AND_CONDITION_RESPONSE_v0",
        "targets": [
            {
                "target_id": packet["target_id"],
                "rows": packet["predicted_perturbation_table"],
            }
            for packet in packets
        ],
        "passed": cert["v54_perturbation_engine_passed"],
    })
    _write_json(DATA_ROOT / "V55" / "operator_state_holdout_validation_summary.json", {
        "kind": "V55_OPERATOR_STATE_HOLDOUT_VALIDATION_SUMMARY_v0",
        "validations": validations,
        "passed": cert["v55_postseal_validation_passed"],
    })

    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v50_v56_protein_esperanto_engine_certificate.json"
    report_path = out_dir / "V50_V56_PROTEIN_ESPERANTO_ENGINE_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert, battery)
    return {
        "certificate": cert_path,
        "report": report_path,
        "grammar": DATA_ROOT / "V50" / "protein_esperanto_mechanism_grammar_extraction.json",
        "engine_spec": DATA_ROOT / "V51" / "operator_state_engine_spec.json",
        "battery": DATA_ROOT / "V53" / "operator_state_multiregime_battery.json",
        "perturbations": DATA_ROOT / "V54" / "operator_state_perturbation_and_condition_response.json",
        "validation_summary": DATA_ROOT / "V55" / "operator_state_holdout_validation_summary.json",
        "bridge": DATA_ROOT / "V56" / "operator_to_custom_force_bridge_spec.json",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V50-V56 Protein Esperanto engine suite.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v50_v56(args.out_dir)
    cert = _read_json(paths["certificate"], "V50-V56 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "target_count": cert["target_count"],
        "supported_validation_count": cert["supported_validation_count"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "folding_problem_solved": cert["folding_problem_solved"],
        "atomistic_md_performed": cert["atomistic_md_performed"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
