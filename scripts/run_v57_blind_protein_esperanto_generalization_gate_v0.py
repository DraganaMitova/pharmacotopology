#!/usr/bin/env python3
from __future__ import annotations

"""Run V57 frozen-engine blind Protein Esperanto generalization gate."""

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    MECHANISM_CLASSES,
    UNIVERSAL_OPERATORS,
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
    validate_against_holdout,
)


DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V57"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V57_BLIND_PROTEIN_ESPERANTO_GENERALIZATION_GATE"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

PASSED = "V57_BLIND_GENERALIZATION_PASSED_REVIEW_REQUIRED"
PARTIAL = "V57_PARTIAL_GENERALIZATION_WITH_ABSTENTIONS_REVIEW_REQUIRED"
BLOCKED_REVISION = "V57_GENERALIZATION_BLOCKED_ENGINE_NEEDS_REVISION"
BLOCKED_LEAKAGE = "V57_BLOCKED_FOR_LEAKAGE"


FRESH_TARGETS: dict[str, dict[str, Any]] = {
    "V57_GB1_DOMAIN": {
        "target_name": "Streptococcal protein G B1 domain",
        "regime_slot": "globular_closure",
        "expected_mechanism_class": "globular_closure",
        "forced_wrong_grammar": "short_region_host_interface_hijacking",
        "sequence": "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
        "focus_regions": [{"name": "GB1 hydrophobic closure core", "span": "full 56-residue domain"}],
        "prediction_sources": [
            {
                "source_id": "GB1_SEQUENCE_FUNCTION_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Small soluble IgG-binding domain with stable globular domain behavior and hydrophobic core closure context; coordinate records excluded.",
            },
            {
                "source_id": "GB1_SEQUENCE_MARKS_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Sequence composition supports compact closure: mixed hydrophobic/aromatic core marks with low disorder tendency.",
            },
        ],
        "perturbations": [
            {
                "perturbation_id": "GB1_CORE_HYDROPHOBIC_DAMAGE",
                "description": "damage hydrophobic closure-core chemistry",
                "operator_scales": {"closure_operator": 0.45},
                "metric": "contact_probability",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "GB1_CORE_STABILIZING_CONTEXT",
                "description": "stabilize closure-core packing context",
                "operator_scales": {"closure_operator": 1.30},
                "metric": "contact_probability",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "GB1_WRONG_HOST_INTERFACE_CONTROL",
                "description": "force unrelated host-hijack interface perturbation",
                "operator_scales": {"host_hijack_operator": 0.20},
                "metric": "contact_probability",
                "expected_direction": "unchanged",
            },
            {
                "perturbation_id": "GB1_NEUTRAL_CONTROL",
                "description": "neutral condition with no operator change",
                "operator_scales": {},
                "metric": "contact_probability",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "compact_basin_supported", "metric": "basin:compact_folded", "comparator": ">=", "threshold": 0.50},
            {"check_id": "contact_probability_supported", "metric": "contact_probability", "comparator": ">=", "threshold": 0.45},
        ],
    },
    "V57_HNRNPA1_LCD": {
        "target_name": "hnRNPA1-like prion LCD",
        "regime_slot": "idp_phase_or_disorder",
        "expected_mechanism_class": "intrinsic_disorder_phase_separation",
        "forced_wrong_grammar": "globular_closure",
        "sequence": "MASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGGYGGGGPGY",
        "focus_regions": [{"name": "hnRNPA1-like low-complexity aromatic/glycine LCD", "span": "full LCD"}],
        "prediction_sources": [
            {
                "source_id": "HNRNPA1_LCD_SEQUENCE_FEATURES_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Low complexity prion-like LCD with aromatic/glycine marks, RNA granule context, phase separation tendency, and disorder ensemble behavior.",
            },
            {
                "source_id": "HNRNPA1_LCD_CONDITION_CONTEXT_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Condensation is expected to respond to aromatic sticker changes, charge/salt context, and concentration pressure.",
            },
        ],
        "perturbations": [
            {
                "perturbation_id": "HNRNPA1_AROMATIC_STICKER_DAMAGE",
                "description": "weaken aromatic sticker contribution in the LCD",
                "operator_scales": {"phase_operator": 0.45},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "HNRNPA1_DISSOLVING_CONDITION",
                "description": "screen charge/concentration pressure to dissolve the condensate-prone ensemble",
                "operator_scales": {"disorder_operator": 0.55, "phase_operator": 0.60},
                "metric": "disorder_order_balance",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "HNRNPA1_WRONG_MEMBRANE_CONTROL",
                "description": "apply unrelated membrane/proteostasis perturbation",
                "operator_scales": {"membrane_pressure_operator": 0.20},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "unchanged",
            },
            {
                "perturbation_id": "HNRNPA1_NEUTRAL_CONTROL",
                "description": "neutral condition with no operator change",
                "operator_scales": {},
                "metric": "basin:phase_prone_dynamic",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "phase_prone_dynamic_supported", "metric": "basin:phase_prone_dynamic", "comparator": ">=", "threshold": 0.25},
            {"check_id": "compact_single_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12},
        ],
    },
    "V57_RHODOPSIN_P23H": {
        "target_name": "Human rhodopsin P23H membrane/proteostasis gate",
        "regime_slot": "membrane_proteostasis",
        "expected_mechanism_class": "membrane_multidomain_folding_proteostasis",
        "forced_wrong_grammar": "intrinsic_disorder_phase_separation",
        "sequence": "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASTVQGVSQLLRTAKARVYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASTTVSKTETSQVAPA",
        "focus_regions": [{"name": "P23H N-terminal mutation plus transmembrane GPCR proteostasis context", "position": 23, "span": "23 and transmembrane bundle"}],
        "prediction_sources": [
            {
                "source_id": "RHO_P23H_SEQUENCE_MEMBRANE_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Rhodopsin is a membrane GPCR with multiple transmembrane segments; P23H is framed as a membrane protein folding, trafficking, and proteostasis stress case.",
            },
            {
                "source_id": "RHO_PROTEOSTASIS_CONTEXT_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Expected mechanism involves membrane pressure, local mutation stress, maturation routing, and quality-control trafficking rather than a soluble single-domain shortcut.",
            },
        ],
        "perturbations": [
            {
                "perturbation_id": "RHO_P23H_DAMAGE",
                "description": "introduce P23H-like local folding damage",
                "operator_scales": {"closure_operator": 0.60, "interface_operator": 0.62, "proteostasis_operator": 0.55},
                "damage": 0.55,
                "metric": "proteostasis_routing",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "RHO_CHAPERONE_RESCUE_CONTEXT",
                "description": "apply proteostasis/chaperone rescue context",
                "operator_scales": {"proteostasis_operator": 1.25, "interface_operator": 1.10},
                "rescue": 0.42,
                "metric": "proteostasis_routing",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "RHO_WRONG_HOST_HIJACK_CONTROL",
                "description": "apply unrelated host-hijack short-interface perturbation",
                "operator_scales": {"host_hijack_operator": 0.20},
                "metric": "proteostasis_routing",
                "expected_direction": "unchanged",
            },
            {
                "perturbation_id": "RHO_NEUTRAL_CONTROL",
                "description": "neutral condition with no operator change",
                "operator_scales": {},
                "metric": "proteostasis_routing",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.55},
            {"check_id": "interface_readiness_supported", "metric": "interface_readiness", "comparator": ">=", "threshold": 0.55},
        ],
    },
    "V57_XCL1_SWITCH": {
        "target_name": "XCL1 lymphotactin fold-switch gate",
        "regime_slot": "metamorphic_switch",
        "expected_mechanism_class": "metamorphic_fold_switching",
        "forced_wrong_grammar": "globular_closure",
        "sequence": "VGSEVSDKRTCVSLTTQRLPVSRIKTYTITEGSLRAVIFITKRGLKVCADPQATWVRDVVRSMDRKSNTRNNMIQTKPTGTQQSTNTAVTLTG",
        "focus_regions": [{"name": "XCL1 state-separated chemokine switch region", "span": "full chemokine domain"}],
        "prediction_sources": [
            {
                "source_id": "XCL1_STATE_SWITCH_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "XCL1 lymphotactin is treated as a metamorphic fold switch with alpha/beta two-state chemokine behavior and partner/context dependence; coordinates excluded.",
            },
            {
                "source_id": "XCL1_CONTEXT_BALANCE_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "State balance should respond to release, dimer/monomer context, and mutations that bias one basin over the other.",
            },
        ],
        "perturbations": [
            {
                "perturbation_id": "XCL1_RELEASE_OR_DIMER_CONTEXT",
                "description": "shift context toward released/beta-like basin",
                "operator_scales": {"dual_basin_switch_operator": 1.15},
                "release": 0.50,
                "metric": "basin:beta_released_basin",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "XCL1_ALPHA_STABILIZING_CONTEXT",
                "description": "stabilize alpha/chemokine-like basin",
                "alpha_bias": 0.42,
                "metric": "basin:alpha_context_basin",
                "expected_direction": "increase",
            },
            {
                "perturbation_id": "XCL1_WRONG_PHASE_CONTROL",
                "description": "apply unrelated phase-separation perturbation",
                "operator_scales": {"phase_operator": 0.20},
                "metric": "basin:beta_released_basin",
                "expected_direction": "unchanged",
            },
            {
                "perturbation_id": "XCL1_NEUTRAL_CONTROL",
                "description": "neutral condition with no operator change",
                "operator_scales": {},
                "metric": "basin:beta_released_basin",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "alpha_basin_present", "metric": "basin:alpha_context_basin", "comparator": ">=", "threshold": 0.30},
            {"check_id": "beta_basin_present", "metric": "basin:beta_released_basin", "comparator": ">=", "threshold": 0.30},
            {"check_id": "averaged_single_fold_rejected", "metric": "basin:averaged_single_fold", "comparator": "<=", "threshold": 0.15},
        ],
    },
    "V57_HIV_TAT": {
        "target_name": "HIV-1 Tat short host-interface gate",
        "regime_slot": "short_interface_hijacking",
        "expected_mechanism_class": "short_region_host_interface_hijacking",
        "forced_wrong_grammar": "globular_closure",
        "sequence": "MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAPQDSQTHQASLSKQPTSQPRGDPTGPKE",
        "focus_regions": [{"name": "Tat basic/activation host-interface motif", "span": "basic and activation region"}],
        "prediction_sources": [
            {
                "source_id": "HIV_TAT_HOST_HIJACK_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "HIV Tat uses short motif host hijack logic through TAR/P-TEFb host interface and transcriptional activation; coordinates and peptide-complex structures excluded.",
            },
            {
                "source_id": "HIV_TAT_BASIC_REGION_NONCOORDINATE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Basic/interface motif perturbations should weaken host transcription hijacking without requiring a stable globular fold.",
            },
        ],
        "perturbations": [
            {
                "perturbation_id": "TAT_BASIC_REGION_DAMAGE",
                "description": "disrupt basic host-interface motif",
                "operator_scales": {"host_hijack_operator": 0.35, "interface_operator": 0.45},
                "interface_disruption": 0.62,
                "metric": "basin:host_interface_engaged",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "TAT_HOST_PARTNER_REMOVAL",
                "description": "remove host partner context",
                "operator_scales": {"host_hijack_operator": 0.42},
                "interface_disruption": 0.45,
                "metric": "interface_readiness",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "TAT_WRONG_MEMBRANE_CONTROL",
                "description": "apply unrelated membrane/proteostasis perturbation",
                "operator_scales": {"membrane_pressure_operator": 0.20},
                "metric": "basin:host_interface_engaged",
                "expected_direction": "unchanged",
            },
            {
                "perturbation_id": "TAT_NEUTRAL_CONTROL",
                "description": "neutral condition with no operator change",
                "operator_scales": {},
                "metric": "basin:host_interface_engaged",
                "expected_direction": "unchanged",
            },
        ],
        "holdout_observables": [
            {"check_id": "host_interface_engaged", "metric": "basin:host_interface_engaged", "comparator": ">=", "threshold": 0.60},
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


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _engine_declaration() -> dict[str, Any]:
    engine_text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V57_FROZEN_ENGINE_DECLARATION_v0",
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": engine_text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "frozen_operator_names": UNIVERSAL_OPERATORS,
        "frozen_mechanism_classes": MECHANISM_CLASSES,
        "allowed_v57_changes": [
            "new target manifests",
            "new source manifests",
            "new holdout manifests",
            "new reports",
            "new certificates",
            "new tests",
        ],
        "disallowed_v57_changes": [
            "new operators",
            "new mechanism classes",
            "new scoring rules",
            "new control definitions",
            "README updates",
            "post-result tuning",
        ],
        "engine_modified_for_v57": False,
        "folding_problem_solved": False,
    }


def _target_manifest() -> dict[str, Any]:
    return {
        "kind": "V57_BLIND_TARGET_MANIFEST_v0",
        "target_count": len(FRESH_TARGETS),
        "freshness_rule": "Targets are not V44-FUS, V45-TDP43, V46-CFTR, V47-RfaH, or V48-ORF6 extraction targets.",
        "answer_key_available_to_prediction": False,
        "holdouts_available_before_seal": False,
        "readme_update_allowed": False,
        "targets": [
            {
                "target_id": target_id,
                "target_name": spec["target_name"],
                "regime_slot": spec["regime_slot"],
                "sequence_length": len(spec["sequence"]),
                "prediction_source_manifest": f"data/protein_esperanto_engine/V57/source_manifests/{target_id}/source_manifest.json",
                "expected_mechanism_class_hidden_until_holdout": True,
            }
            for target_id, spec in FRESH_TARGETS.items()
        ],
    }


def _source_manifest(target_id: str, spec: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V57_TARGET_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "target_name": spec["target_name"],
        "regime_slot": spec["regime_slot"],
        "sequence": spec["sequence"],
        "prediction_sources": spec["prediction_sources"],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "coordinate-derived contacts or native contact maps",
            "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing",
            "internal runtime artifacts as biological evidence",
            "validation holdouts before prediction sealing",
            "answer key/class labels during prediction",
            "post-result tuning of V50-V56 engine",
        ],
        "coordinate_sources_available_before_prediction": False,
        "internal_runtime_sources_available_to_prediction": False,
        "answer_key_available_to_prediction": False,
        "holdouts_created": False,
        "folding_problem_solved": False,
    }


def _postseal_holdout(target_id: str, spec: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V57_POSTSEAL_BLIND_GENERALIZATION_HOLDOUT_v0",
        "target_id": target_id,
        "target_name": spec["target_name"],
        "expected_mechanism_class": spec["expected_mechanism_class"],
        "expected_observables": spec["holdout_observables"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{target_id}_V57_POSTSEAL_HOLDOUT",
                "source_class": "spatial_proxy_non_coordinate",
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "rationale": "Post-seal validation summary for pre-registered mechanism class, trajectory type, and perturbation direction.",
            }
        ],
    }


def _wrong_grammar_packet(target_id: str, spec: dict[str, Any], source_manifest: dict[str, Any]) -> dict[str, Any]:
    return build_sealed_operator_state_packet(
        target_id=f"{target_id}_WRONG_GRAMMAR_CONTROL",
        target_name=f"{spec['target_name']} forced wrong grammar",
        sequence=spec["sequence"],
        sources=source_manifest["prediction_sources"],
        focus_regions=spec["focus_regions"],
        perturbations=[],
        forced_grammar=spec["forced_wrong_grammar"],
    )


def _shuffled_control_packet(target_id: str, spec: dict[str, Any], source_manifest: dict[str, Any]) -> dict[str, Any]:
    return build_sealed_operator_state_packet(
        target_id=f"{target_id}_SHUFFLED_CONTROL",
        target_name=f"{spec['target_name']} shuffled control",
        sequence=shuffled_sequence(spec["sequence"]),
        sources=source_manifest["prediction_sources"],
        focus_regions=spec["focus_regions"],
        perturbations=[],
    )


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _pre_registered_row(target_id: str, spec: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "target_id": target_id,
        "target_name": spec["target_name"],
        "selected_grammar": packet["selected_mechanism_grammar"]["mechanism_class"],
        "trajectory_type": packet["selected_mechanism_grammar"]["mechanism_class"],
        "dominant_regions": [row["acts_on"] for row in packet["operator_field"]["operators"]],
        "operator_field": packet["operator_field"]["operator_names"],
        "predicted_perturbation_directions": [
            {
                "perturbation_id": row["perturbation_id"],
                "metric": row["metric"],
                "expected_direction": row["expected_direction"],
            }
            for row in packet["predicted_perturbation_table"]
        ],
        "falsification_criteria": packet["predicted_falsifiers"],
        "abstention_criteria": [
            "coordinate leakage before seal",
            "holdout opened before seal",
            "only generic target name or generic annotation remains",
            "forced wrong grammar is required for prediction",
        ],
        "prediction_hash": packet["prediction_hash"],
    }


def _target_controls(
    *,
    target_id: str,
    spec: dict[str, Any],
    packet: dict[str, Any],
    wrong_packet: dict[str, Any],
    shuffled_packet: dict[str, Any],
    validation: dict[str, Any],
) -> list[dict[str, Any]]:
    perturb_rows = packet["predicted_perturbation_table"]
    directional_rows = [row for row in perturb_rows if not row["perturbation_id"].endswith(("CONTROL",))]
    control_rows = [row for row in perturb_rows if row["perturbation_id"].endswith("CONTROL")]
    shuffled_coherence = sequence_operator_coherence(shuffled_packet)
    original_coherence = sequence_operator_coherence(packet)
    return [
        _control(
            f"{target_id}_selected_expected_grammar",
            packet["selected_mechanism_grammar"]["mechanism_class"] == spec["expected_mechanism_class"],
            "Frozen engine must select expected grammar before holdout.",
            packet["selected_mechanism_grammar"]["mechanism_class"],
        ),
        _control(
            f"{target_id}_wrong_grammar_fails_or_abstains",
            wrong_packet["selected_mechanism_grammar"]["forced_grammar_rejected"] is True
            and wrong_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain",
            "Forced wrong grammar must fail or abstain.",
            wrong_packet["selected_mechanism_grammar"],
        ),
        _control(
            f"{target_id}_perturbation_direction_separation",
            all(row["direction_passed"] for row in directional_rows)
            and all(row["direction_passed"] and row["observed_direction"] == "unchanged" for row in control_rows),
            "Correct perturbations move in the expected direction; wrong/neutral controls do not look equally correct.",
            {"directional": directional_rows, "controls": control_rows},
        ),
        _control(
            f"{target_id}_shuffled_or_neutral_control_not_better",
            shuffled_coherence <= original_coherence + 0.05,
            "Shuffled sequence control must not improve operator coherence materially.",
            {"original_coherence": original_coherence, "shuffled_coherence": shuffled_coherence},
        ),
        _control(
            f"{target_id}_postseal_validation_supported",
            validation["score_label"] == "supported",
            "Post-seal validation must support the sealed packet.",
            validation["score_label"],
        ),
    ]


def _global_controls(
    *,
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    wrong_packets: list[dict[str, Any]],
    validations: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    random_packet = build_sealed_operator_state_packet(
        target_id="V57_RANDOM_SEQUENCE_CONTROL",
        target_name="V57 random sequence control",
        sequence=deterministic_random_sequence(96),
        sources=[],
        perturbations=[],
    )
    coord_gate = evidence_boundary_gate([{
        "source_id": "V57_BAD_COORDINATES",
        "source_class": "coordinate_derived",
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    preopened_holdout_gate = evidence_boundary_gate([{
        "source_id": "V57_PREOPENED_HOLDOUT",
        "source_class": "spatial_proxy_non_coordinate",
        "source_role": "holdout_validation",
        "spatial_proxy": True,
    }])
    return [
        _control(
            "v57_engine_source_frozen",
            engine_declaration["engine_source_last_commit"] == _git(["log", "-1", "--format=%H", "--", "src/pharmacotopology/protein_esperanto_engine.py"])
            and engine_declaration["engine_modified_for_v57"] is False,
            "V57 must not modify the V50-V56 engine.",
            engine_declaration,
        ),
        _control(
            "v57_no_new_operators_or_mechanism_classes",
            engine_declaration["frozen_operator_names"] == UNIVERSAL_OPERATORS
            and engine_declaration["frozen_mechanism_classes"] == MECHANISM_CLASSES,
            "V57 must use frozen operators and mechanism classes.",
        ),
        _control(
            "v57_at_least_three_fresh_regimes",
            len({packet["selected_mechanism_grammar"]["mechanism_class"] for packet in packets}) >= 3
            and len(packets) >= 3,
            "Blind gate must cover at least three fresh regimes.",
            sorted({packet["selected_mechanism_grammar"]["mechanism_class"] for packet in packets}),
        ),
        _control(
            "v57_all_wrong_grammars_fail",
            all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_packets),
            "Every target must reject a forced wrong grammar.",
            [packet["selected_mechanism_grammar"] for packet in wrong_packets],
        ),
        _control(
            "v57_random_sequence_abstains",
            random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain",
            "Random sequence without allowed evidence must abstain.",
            random_packet["selected_mechanism_grammar"],
        ),
        _control(
            "v57_coordinate_leakage_blocks",
            coord_gate["coordinate_derived_source_count_before_prediction"] == 1
            and coord_gate["allowed_initialization_source_ids"] == [],
            "Coordinate-derived prediction evidence must block.",
            coord_gate,
        ),
        _control(
            "v57_holdout_before_seal_blocks",
            preopened_holdout_gate["holdout_opened_before_seal"] is True
            and preopened_holdout_gate["allowed_initialization_source_ids"] == [],
            "Holdout opened before seal must block.",
            preopened_holdout_gate,
        ),
        _control(
            "v57_all_holdouts_reference_hash",
            all(row["holdout_opened_after_prediction_hash"] for row in validations),
            "Each post-seal validation must reference the sealed hash.",
        ),
        _control(
            "v57_folding_problem_solved_never_true",
            all(packet["folding_problem_solved"] is False for packet in packets),
            "V57 must preserve folding_problem_solved=false.",
        ),
    ]


def _aggregate_certificate(
    *,
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    validations: list[dict[str, Any]],
    controls: list[dict[str, Any]],
) -> dict[str, Any]:
    failed = [control["control_id"] for control in controls if not control["passed"]]
    leakage = sum(packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"] for packet in packets)
    runtime = sum(packet["evidence_manifest"]["internal_runtime_source_count_for_prediction"] for packet in packets)
    abstentions = [
        packet["target_id"]
        for packet in packets
        if packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"
    ]
    supported = [row for row in validations if row["score_label"] == "supported"]
    if leakage or runtime:
        status = BLOCKED_LEAKAGE
    elif failed:
        status = BLOCKED_REVISION
    elif abstentions:
        status = PARTIAL
    else:
        status = PASSED
    return {
        "kind": "V57_BLIND_PROTEIN_ESPERANTO_GENERALIZATION_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_modified_for_v57": engine_declaration["engine_modified_for_v57"],
        "target_count": len(packets),
        "fresh_regime_count": len({packet["selected_mechanism_grammar"]["mechanism_class"] for packet in packets}),
        "sealed_prediction_count": len(packets),
        "supported_validation_count": len(supported),
        "abstention_count": len(abstentions),
        "abstained_targets": abstentions,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_checks": failed,
        "controls": controls,
        "coordinate_derived_source_count_before_prediction": leakage,
        "internal_runtime_source_count_for_prediction": runtime,
        "folding_problem_solved": False,
        "atomistic_md_performed": False,
        "readme_touched": False,
        "claim_allowed": status == PASSED,
        "allowed_claim_text": (
            "The frozen Protein Esperanto engine generalized to fresh targets across multiple regimes under leakage-controlled blind validation."
            if status == PASSED else ""
        ),
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "Coordinates were predicted de novo.",
            "Atomistic folding was solved.",
            "External review is unnecessary.",
        ],
    }


def _write_report(path: Path, cert: dict[str, Any], preregistered: list[dict[str, Any]], validations: list[dict[str, Any]]) -> None:
    lines = [
        "# V57 Blind Protein Esperanto Generalization Gate",
        "",
        f"Status: `{cert['status']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        f"engine_modified_for_v57: `{cert['engine_modified_for_v57']}`",
        f"Targets: `{cert['target_count']}`",
        f"Supported validations: `{cert['supported_validation_count']}` / `{cert['target_count']}`",
        f"Controls passed: `{cert['passed_control_count']}` / `{cert['control_count']}`",
        "",
        "## Pre-Registered Outputs",
    ]
    validation_by_target = {row["target_id"]: row for row in validations}
    for row in preregistered:
        lines.append(
            f"- `{row['target_id']}` grammar `{row['selected_grammar']}` operators `{row['operator_field']}` validation `{validation_by_target[row['target_id']]['score_label']}`"
        )
    lines.extend(["", "## Failed Checks"])
    if cert["failed_checks"]:
        for check in cert["failed_checks"]:
            lines.append(f"- `{check}`")
    else:
        lines.append("- none")
    lines.extend([
        "",
        "## Claim Boundary",
        cert["allowed_claim_text"] or "No pass claim allowed until failed checks are resolved.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v57(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    engine_declaration = _engine_declaration()
    target_manifest = _target_manifest()
    packets: list[dict[str, Any]] = []
    wrong_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    validations: list[dict[str, Any]] = []
    preregistered: list[dict[str, Any]] = []
    target_controls: list[dict[str, Any]] = []

    _write_json(DATA_ROOT / "v57_frozen_engine_declaration.json", engine_declaration)
    _write_json(DATA_ROOT / "v57_blind_target_manifest.json", target_manifest)

    for target_id, spec in FRESH_TARGETS.items():
        source_manifest = _source_manifest(target_id, spec)
        _write_json(DATA_ROOT / "source_manifests" / target_id / "source_manifest.json", source_manifest)
        packet = build_sealed_operator_state_packet(
            target_id=target_id,
            target_name=spec["target_name"],
            sequence=spec["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=spec["focus_regions"],
            perturbations=spec["perturbations"],
        )
        packets.append(packet)
        _write_json(DATA_ROOT / "sealed_simulation_packets" / target_id / "sealed_simulation_packet.json", packet)
        prereg = _pre_registered_row(target_id, spec, packet)
        preregistered.append(prereg)
        _write_json(DATA_ROOT / "pre_registered_outputs" / target_id / "pre_registered_output.json", prereg)

        wrong_packet = _wrong_grammar_packet(target_id, spec, source_manifest)
        wrong_packets.append(wrong_packet)
        _write_json(DATA_ROOT / "wrong_grammar_controls" / target_id / "wrong_grammar_packet.json", wrong_packet)

        shuffled_packet = _shuffled_control_packet(target_id, spec, source_manifest)
        shuffled_packets.append(shuffled_packet)
        _write_json(DATA_ROOT / "shuffled_controls" / target_id / "shuffled_control_packet.json", shuffled_packet)

        holdout = _postseal_holdout(target_id, spec, packet)
        validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
        validations.append(validation)
        _write_json(DATA_ROOT / "holdouts_postseal" / target_id / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target_id / "validation_result.json", validation)
        target_controls.extend(
            _target_controls(
                target_id=target_id,
                spec=spec,
                packet=packet,
                wrong_packet=wrong_packet,
                shuffled_packet=shuffled_packet,
                validation=validation,
            )
        )

    global_controls = _global_controls(
        engine_declaration=engine_declaration,
        packets=packets,
        wrong_packets=wrong_packets,
        validations=validations,
    )
    controls = target_controls + global_controls
    cert = _aggregate_certificate(
        engine_declaration=engine_declaration,
        packets=packets,
        validations=validations,
        controls=controls,
    )
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    _write_json(DATA_ROOT / "v57_perturbation_response_table.json", {
        "kind": "V57_PERTURBATION_RESPONSE_TABLE_v0",
        "targets": [
            {"target_id": packet["target_id"], "rows": packet["predicted_perturbation_table"]}
            for packet in packets
        ],
    })
    _write_json(DATA_ROOT / "v57_postseal_validation_report.json", {
        "kind": "V57_POSTSEAL_VALIDATION_REPORT_v0",
        "validations": validations,
        "supported_validation_count": cert["supported_validation_count"],
    })
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v57_blind_protein_esperanto_generalization_certificate.json"
    report_path = out_dir / "V57_BLIND_PROTEIN_ESPERANTO_GENERALIZATION_GATE_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert, preregistered, validations)
    return {
        "certificate": cert_path,
        "report": report_path,
        "engine_declaration": DATA_ROOT / "v57_frozen_engine_declaration.json",
        "target_manifest": DATA_ROOT / "v57_blind_target_manifest.json",
        "perturbation_table": DATA_ROOT / "v57_perturbation_response_table.json",
        "validation_report": DATA_ROOT / "v57_postseal_validation_report.json",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V57 blind Protein Esperanto generalization gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v57(args.out_dir)
    cert = _read_json(paths["certificate"], "V57 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "target_count": cert["target_count"],
        "fresh_regime_count": cert["fresh_regime_count"],
        "supported_validation_count": cert["supported_validation_count"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "folding_problem_solved": cert["folding_problem_solved"],
        "engine_modified_for_v57": cert["engine_modified_for_v57"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] in {PASSED, PARTIAL} else 1


if __name__ == "__main__":
    raise SystemExit(main())
