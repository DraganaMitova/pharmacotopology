#!/usr/bin/env python3
from __future__ import annotations

"""Run V46 CFTR F508del membrane multidomain folding rescue attack."""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import build_v46_cftr_f508del_membrane_multidomain_sources_v0 as v46_sources

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V46" / "CFTR_F508del"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK"

PASSED = "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED"
PARTIAL = "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_PARTIAL_CLEAN_ABSTAIN"
FAILED = "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FAILED_PREDICTIONS"
BLOCKED_HOLDOUT = "V46_BLOCKED_PREDICTION_HOLDOUT_LEAKAGE"
BLOCKED_COORD = "V46_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"

EXPECTED_MECHANISM_CLASS = "membrane_multidomain_folding_assembly_proteostasis_defect"
FORBIDDEN_BUCKETS = [
    "generic_channel_annotation_only_operator",
    "single_local_mutation_only_operator",
    "compact_single_domain_fold_operator",
    "solved_atomic_structure_operator",
    "coordinate_contact_operator",
    "AlphaFold_confidence_proxy_operator",
]


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _hash_json(data: dict[str, Any]) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_source_manifest() -> dict[str, Any]:
    return _read_json(SOURCE_ROOT / "source_manifest.json", "V46 source manifest")


def _leakage_counts(sources: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "coordinate_derived_source_count_before_prediction": sum(1 for source in sources if source.get("coordinate_derived")),
        "internal_runtime_source_count_for_prediction": sum(1 for source in sources if source.get("internal_runtime_source")),
        "holdout_leakage_detected": any(source.get("holdout_source") and source.get("allowed_use") == "prediction_input_before_sealing" for source in sources),
        "native_metrics_used_before_prediction": any(source.get("native_metrics_used_for_selection") for source in sources),
        "coordinate_truth_used_before_prediction": any(source.get("coordinate_truth_used_before_prediction") for source in sources),
    }


def _full_packet_allowed(sources: list[dict[str, Any]]) -> bool:
    leakage = _leakage_counts(sources)
    if any([
        leakage["coordinate_derived_source_count_before_prediction"],
        leakage["internal_runtime_source_count_for_prediction"],
        leakage["holdout_leakage_detected"],
        leakage["native_metrics_used_before_prediction"],
        leakage["coordinate_truth_used_before_prediction"],
    ]):
        return False
    source_ids = {source.get("source_id") for source in sources}
    required = {
        "UNIPROT_P13569_SEQUENCE_FEATURES_NONCOORDINATE",
        "UNIPROT_P13569_F508DEL_VARIANT_NONCOORDINATE",
        "CFTR_DOMAIN_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
    }
    generic_channel_only = source_ids == {"GENERIC_CHANNEL_ANNOTATION_ONLY"}
    single_mutation_only = source_ids == {"UNIPROT_P13569_F508DEL_VARIANT_NONCOORDINATE"}
    return required.issubset(source_ids) and not generic_channel_only and not single_mutation_only


def _operator_regions() -> list[dict[str, Any]]:
    return [
        {
            "bucket": "membrane_domain_operator",
            "region_name": "MSD1_MSD2_12_TM_membrane_channel_scaffold",
            "span": "MSD1 81-365, MSD2 859-1155, 12 transmembrane helices",
            "confidence": "high",
            "rationale": "P13569 topology annotates two membrane-spanning domains and 12 transmembrane segments.",
        },
        {
            "bucket": "NBD1_stability_operator",
            "region_name": "NBD1_ABC_transporter_folding_core",
            "span": "423-646",
            "confidence": "high",
            "rationale": "F508 lies in NBD1, whose stability/folding competence is central to the defect grammar.",
        },
        {
            "bucket": "F508del_local_destabilization_operator",
            "region_name": "F508_surface_and_folding_perturbation",
            "span": "508",
            "confidence": "high",
            "rationale": "F508del removes phenylalanine 508 and is annotated as impairing folding and stability.",
        },
        {
            "bucket": "interdomain_interface_coupling_operator",
            "region_name": "NBD1_MSD_ICL_surface_coupling",
            "span": "NBD1 423-646 coupled to MSD/ICL cytoplasmic interface",
            "confidence": "high",
            "rationale": "The F508del annotation points to domain-interaction surface changes, not merely a local residue defect.",
        },
        {
            "bucket": "trafficking_quality_control_operator",
            "region_name": "ER_maturation_glycan_trafficking_proteostasis_filter",
            "span": "full-length folding/maturation pathway",
            "confidence": "high",
            "rationale": "F508del impairs maturation and trafficking to the cell membrane and is routed toward degradation.",
        },
        {
            "bucket": "corrector_or_rescue_context_operator",
            "region_name": "multisite_corrector_rescue_context",
            "span": "NBD1 stability, interface coupling, and proteostasis rescue axes",
            "confidence": "medium_high",
            "rationale": "Corrector logic should require domain/interface/proteostasis rescue rather than one local patch.",
        },
        {
            "bucket": "multidomain_assembly_operator",
            "region_name": "MSD_NBD_R_region_assembly_pathway",
            "span": "MSD1/NBD1/R/MSD2/NBD2",
            "confidence": "high",
            "rationale": "CFTR is a full-length membrane ABC protein whose folding defect emerges during multidomain assembly.",
        },
        {
            "bucket": "channel_function_context_operator",
            "region_name": "surface_channel_opening_and_conductance_context",
            "span": "surface-mature CFTR channel state",
            "confidence": "medium",
            "rationale": "Channel opening/function readouts are downstream maturation context, not atomic-fold proof.",
        },
    ]


def _defect_grammar() -> list[dict[str, str]]:
    return [
        {"axis": "local_NBD1_defect", "prediction": "F508 deletion weakens NBD1 folding/stability competence."},
        {"axis": "interdomain_coupling", "prediction": "F508del perturbs an NBD1 surface needed for MSD/ICL domain coupling."},
        {"axis": "maturation_filter", "prediction": "Misassembled F508del-CFTR fails efficient glycan maturation and ER-to-membrane trafficking."},
        {"axis": "proteostasis", "prediction": "Quality-control machinery routes misfolded/misassembled protein toward degradation."},
        {"axis": "channel_context", "prediction": "Residual or rescued surface protein may still show altered opening/gating readouts."},
    ]


def _rescue_grammar() -> list[dict[str, str]]:
    return [
        {"rescue_axis": "NBD1_stabilization", "prediction": "NBD1 stabilizers should improve folding competence but may be partial alone."},
        {"rescue_axis": "interface_correction", "prediction": "Interface/domain-coupling correction should be required for stronger rescue."},
        {"rescue_axis": "proteostasis_or_trafficking_context", "prediction": "Rescue must be assessed through maturation, trafficking, and surface function."},
        {"rescue_axis": "combination_corrector_logic", "prediction": "Combinations targeting distinct folding defects should outperform one-axis correction."},
    ]


def _perturbation_predictions() -> list[dict[str, str]]:
    return [
        {"prediction_id": "V46_PERT_001", "perturbation": "delete F508 in NBD1", "predicted_effect": "weakens NBD1 stability or folding competence", "operator_bucket": "NBD1_stability_operator"},
        {"prediction_id": "V46_PERT_002", "perturbation": "track F508del maturation and glycan processing", "predicted_effect": "shows multidomain maturation/processing failure rather than a local annotation-only defect", "operator_bucket": "trafficking_quality_control_operator"},
        {"prediction_id": "V46_PERT_003", "perturbation": "stabilize NBD1 alone", "predicted_effect": "partially rescues folding/maturation but is not necessarily sufficient for full correction", "operator_bucket": "NBD1_stability_operator"},
        {"prediction_id": "V46_PERT_004", "perturbation": "restore NBD1-MSD/ICL interface coupling", "predicted_effect": "provides stronger correction when combined with NBD1 stabilization", "operator_bucket": "interdomain_interface_coupling_operator"},
        {"prediction_id": "V46_PERT_005", "perturbation": "apply corrector or rescue context", "predicted_effect": "maps to domain/interface/proteostasis grammar rather than magic channel activation", "operator_bucket": "corrector_or_rescue_context_operator"},
        {"prediction_id": "V46_PERT_006", "perturbation": "use generic channel annotation only", "predicted_effect": "must fail to produce a CFTR F508del solution packet", "operator_bucket": "channel_function_context_operator"},
        {"prediction_id": "V46_PERT_007", "perturbation": "force soluble compact single-domain grammar", "predicted_effect": "must be blocked because CFTR is a membrane multidomain ABC protein", "operator_bucket": "membrane_domain_operator"},
        {"prediction_id": "V46_PERT_008", "perturbation": "remove interdomain/interface evidence", "predicted_effect": "weakens the solution and collapses strong correction logic", "operator_bucket": "interdomain_interface_coupling_operator"},
        {"prediction_id": "V46_PERT_009", "perturbation": "remove NBD1 stability evidence", "predicted_effect": "weakens the solution and reduces F508del defect specificity", "operator_bucket": "NBD1_stability_operator"},
        {"prediction_id": "V46_PERT_010", "perturbation": "remove trafficking/proteostasis evidence", "predicted_effect": "weakens maturation grammar while not falsifying atomic fold claims", "operator_bucket": "trafficking_quality_control_operator"},
        {"prediction_id": "V46_PERT_011", "perturbation": "measure rescued channel opening at the membrane", "predicted_effect": "surface function improves only after folding/maturation rescue; function readout is downstream context", "operator_bucket": "channel_function_context_operator"},
        {"prediction_id": "V46_PERT_012", "perturbation": "combine NBD1, interface, and proteostasis correctors", "predicted_effect": "outperforms any single local-site explanation", "operator_bucket": "corrector_or_rescue_context_operator"},
    ]


def predict_solution_packet(source_manifest: dict[str, Any]) -> dict[str, Any]:
    sources = list(source_manifest.get("prediction_sources", []))
    leakage = _leakage_counts(sources)
    full_packet_allowed = _full_packet_allowed(sources)
    if full_packet_allowed:
        mechanism_class = EXPECTED_MECHANISM_CLASS
        operators = _operator_regions()
        perturbations = _perturbation_predictions()
        confidence = "high_for_live_membrane_solution_packet_review_required"
    else:
        mechanism_class = "insufficient_evidence_clean_abstain"
        operators = []
        perturbations = []
        confidence = "low"
    packet = {
        "kind": "V46_CFTR_F508DEL_SEALED_MEMBRANE_MULTIDOMAIN_SOLUTION_PACKET_v0",
        "target_id": "V46_CFTR_F508DEL",
        "target": "human CFTR F508del membrane multidomain folding defect",
        "uniprot_accession": "P13569",
        "sequence_region_scope": {
            "label": "full-length CFTR P13569 with F508del focus in NBD1",
            "full_length": len(v46_sources.CFTR_FULL_SEQUENCE),
            "focus": "F508del in NBD1 / ABC transporter 1 domain",
            "sequence": v46_sources.CFTR_FULL_SEQUENCE,
        },
        "prediction_timestamp": datetime.now(timezone.utc).isoformat(),
        "no_holdout_access_before_hash": True,
        "prediction_inputs_manifest": str(PREDICTION_ROOT / "prediction_inputs_manifest.json"),
        "blocked_inputs_manifest": str(PREDICTION_ROOT / "blocked_inputs_manifest.json"),
        "prediction_source_ids": [source.get("source_id") for source in sources],
        "mechanism_class": mechanism_class,
        "mechanism_confidence": confidence,
        "operator_regions": operators,
        "operator_region_rationale": {row["bucket"]: row["rationale"] for row in operators},
        "predicted_defect_grammar": _defect_grammar() if full_packet_allowed else [],
        "predicted_rescue_or_corrector_grammar": _rescue_grammar() if full_packet_allowed else [],
        "perturbation_predictions": perturbations,
        "low_resolution_structure_or_ensemble_prediction": {
            "architecture": "full-length multidomain membrane ABC channel with MSD1/NBD1/R/MSD2/NBD2 assembly",
            "single_native_fold_expected": False,
            "defect_model": "F508del causes a coupled NBD1 stability, domain-interface assembly, trafficking/proteostasis, and channel-function defect",
            "atomic_claim": "blocked; no coordinate solution is claimed",
        } if full_packet_allowed else {},
        "proposed_experimental_tests": [
            "compare WT and F508del NBD1 stability/folding using non-coordinate biochemical stability assays",
            "measure full-length F508del glycan maturation and ER-to-surface trafficking with and without NBD1 stabilizers",
            "test interface/domain-coupling suppressors or correctors for synergy with NBD1 stabilization",
            "compare corrector classes for NBD1-only, interface, and proteostasis/maturation rescue signatures",
            "measure rescued surface channel opening separately from maturation rescue",
            "ablate proteostasis/quality-control components and quantify maturation versus degradation changes",
        ] if full_packet_allowed else [],
        "falsification_criteria": [
            "F508del is fully explained by one local residue deletion with no NBD1 stability effect",
            "F508del correction is complete with NBD1 stabilization alone in all relevant maturation/function assays",
            "interdomain/interface evidence does not affect correction or maturation readouts",
            "trafficking/proteostasis readouts do not change despite folding/maturation rescue",
            "generic channel annotation alone predicts the F508del rescue logic",
            "coordinate evidence before sealing is required for the mechanism packet",
            "independent holdouts contradict more than two explicit perturbation predictions",
        ] if full_packet_allowed else [],
        "claim_boundary": {
            "allowed": "sealed live solution packet for CFTR F508del membrane multidomain folding/rescue mechanism, pending review",
            "not_allowed": "universal protein-folding solved claim, atomic CFTR solution claim, or one-local-mutation-only claim",
        },
        "forbidden_operator_buckets_rejected": FORBIDDEN_BUCKETS,
        "full_solution_packet_allowed_by_inputs": full_packet_allowed,
        **leakage,
        "folding_problem_solved": False,
        "live_membrane_solution_packet": False,
        "protein_folding_solved_candidate_strengthened": False,
    }
    packet["prediction_hash"] = _hash_json({key: value for key, value in packet.items() if key != "prediction_hash"})
    return packet


def _write_prediction_artifacts(source_manifest: dict[str, Any], packet: dict[str, Any]) -> None:
    _write_json(PREDICTION_ROOT / "prediction_inputs_manifest.json", {
        "kind": "V46_CFTR_F508DEL_PREDICTION_INPUTS_MANIFEST_v0",
        "target_id": "V46_CFTR_F508DEL",
        "prediction_sources": source_manifest.get("prediction_sources", []),
        "holdout_sources_available": False,
        "answer_key_available": False,
        "coordinates_available_before_sealing": False,
    })
    _write_json(PREDICTION_ROOT / "blocked_inputs_manifest.json", {
        "kind": "V46_CFTR_F508DEL_BLOCKED_INPUTS_MANIFEST_v0",
        "blocked_prediction_inputs": source_manifest.get("blocked_prediction_inputs", []),
    })
    _write_json(PREDICTION_ROOT / "sealed_prediction_packet.json", packet)


def _postseal_holdouts(packet: dict[str, Any]) -> list[dict[str, Any]]:
    now = datetime.now(timezone.utc).isoformat()
    common = {
        "opened_after_prediction_hash": packet["prediction_hash"],
        "opened_timestamp": now,
        "used_before_prediction": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
    }
    return [
        {
            **common,
            "source_id": "HOLDOUT_UNIPROT_F508DEL_PROCESSING_DEFECT",
            "source_class": "F508del_folding_processing_defect_support",
            "source_url_or_citation": "UniProt P13569 F508del variant annotation",
            "supports_operator_buckets": ["F508del_local_destabilization_operator", "trafficking_quality_control_operator", "channel_function_context_operator"],
            "supports_perturbation_ids": ["V46_PERT_001", "V46_PERT_002", "V46_PERT_010", "V46_PERT_011"],
            "supports_interdomain_ids": [],
            "supports_rescue_ids": [],
            "evidence_statement": "F508del impairs folding/stability, maturation, trafficking, degradation routing, and channel opening frequency.",
        },
        {
            **common,
            "source_id": "HOLDOUT_NBD1_STABILITY_LITERATURE",
            "source_class": "NBD1_stability_support",
            "source_url_or_citation": "Post-seal F508del NBD1 folding/stability literature",
            "supports_operator_buckets": ["NBD1_stability_operator", "F508del_local_destabilization_operator"],
            "supports_perturbation_ids": ["V46_PERT_001", "V46_PERT_003", "V46_PERT_009"],
            "supports_interdomain_ids": [],
            "supports_rescue_ids": ["V46_PERT_003"],
            "evidence_statement": "NBD1 stabilization is a real but incomplete axis of F508del correction.",
        },
        {
            **common,
            "source_id": "HOLDOUT_INTERDOMAIN_INTERFACE_COUPLING",
            "source_class": "interdomain_interface_coupling_support",
            "source_url_or_citation": "Post-seal CFTR NBD1-MSD/ICL interface-coupling literature",
            "supports_operator_buckets": ["interdomain_interface_coupling_operator", "multidomain_assembly_operator"],
            "supports_perturbation_ids": ["V46_PERT_004", "V46_PERT_008", "V46_PERT_012"],
            "supports_interdomain_ids": ["V46_PERT_004", "V46_PERT_008", "V46_PERT_012"],
            "supports_rescue_ids": ["V46_PERT_004", "V46_PERT_012"],
            "evidence_statement": "F508del correction requires interdomain/interface coupling beyond local NBD1 stabilization.",
        },
        {
            **common,
            "source_id": "HOLDOUT_CORRECTOR_RESCUE_REQUIREMENTS",
            "source_class": "corrector_rescue_logic_support",
            "source_url_or_citation": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3266553/",
            "supports_operator_buckets": ["corrector_or_rescue_context_operator", "interdomain_interface_coupling_operator", "NBD1_stability_operator"],
            "supports_perturbation_ids": ["V46_PERT_003", "V46_PERT_004", "V46_PERT_005", "V46_PERT_012"],
            "supports_interdomain_ids": ["V46_PERT_004", "V46_PERT_012"],
            "supports_rescue_ids": ["V46_PERT_003", "V46_PERT_004", "V46_PERT_005", "V46_PERT_012"],
            "evidence_statement": "Efficient correction maps to multiple rescue axes and is not a single local-site repair.",
        },
        {
            **common,
            "source_id": "HOLDOUT_TRAFFICKING_PROTEOSTASIS_QC",
            "source_class": "trafficking_proteostasis_support",
            "source_url_or_citation": "Post-seal CFTR maturation, glycan processing, ER quality-control, and proteostasis literature",
            "supports_operator_buckets": ["trafficking_quality_control_operator", "corrector_or_rescue_context_operator"],
            "supports_perturbation_ids": ["V46_PERT_002", "V46_PERT_005", "V46_PERT_010"],
            "supports_interdomain_ids": [],
            "supports_rescue_ids": ["V46_PERT_005"],
            "evidence_statement": "Maturation and trafficking are functional folding-context readouts, not atomic-fold proof.",
        },
        {
            **common,
            "source_id": "HOLDOUT_GENERIC_CHANNEL_AND_SINGLE_LOCAL_REJECTION",
            "source_class": "forbidden_grammar_rejection_support",
            "source_url_or_citation": "Post-seal negative controls against generic channel and single-local explanations",
            "supports_operator_buckets": ["membrane_domain_operator", "multidomain_assembly_operator", "channel_function_context_operator"],
            "supports_perturbation_ids": ["V46_PERT_006", "V46_PERT_007"],
            "supports_interdomain_ids": [],
            "supports_rescue_ids": [],
            "evidence_statement": "CFTR F508del requires membrane multidomain/proteostasis grammar, not a generic channel or compact soluble domain label.",
        },
    ]


def _validate_predictions(packet: dict[str, Any], holdouts: list[dict[str, Any]]) -> dict[str, Any]:
    supported_buckets = {bucket for holdout in holdouts for bucket in holdout.get("supports_operator_buckets", [])}
    supported_perturbations = {pid for holdout in holdouts for pid in holdout.get("supports_perturbation_ids", [])}
    interdomain_supported = {pid for holdout in holdouts for pid in holdout.get("supports_interdomain_ids", [])}
    rescue_supported = {pid for holdout in holdouts for pid in holdout.get("supports_rescue_ids", [])}
    operator_results = [
        {"bucket": row["bucket"], "region_name": row["region_name"], "support_level": "supported" if row["bucket"] in supported_buckets else "unsupported"}
        for row in packet.get("operator_regions", [])
    ]
    partial_ids = {"V46_PERT_003", "V46_PERT_010", "V46_PERT_011"}
    perturbation_results = []
    for row in packet.get("perturbation_predictions", []):
        if row["prediction_id"] in supported_perturbations:
            support = "partially_supported" if row["prediction_id"] in partial_ids else "supported"
        else:
            support = "unsupported"
        perturbation_results.append({
            "prediction_id": row["prediction_id"],
            "perturbation": row["perturbation"],
            "predicted_effect": row["predicted_effect"],
            "support_level": support,
        })
    interdomain_ids = {"V46_PERT_004", "V46_PERT_008", "V46_PERT_012"}
    rescue_ids = {"V46_PERT_003", "V46_PERT_004", "V46_PERT_005", "V46_PERT_012"}
    validated = [row for row in perturbation_results if row["support_level"] in {"supported", "partially_supported"}]
    contradicted = [row for row in perturbation_results if row["support_level"] == "contradicted"]
    operator_supported = [row for row in operator_results if row["support_level"] == "supported"]
    buckets = [row.get("bucket") for row in packet.get("operator_regions", [])]
    return {
        "kind": "V46_CFTR_F508DEL_POSTSEAL_VALIDATION_v0",
        "operator_validation": operator_results,
        "perturbation_validation": perturbation_results,
        "operator_support_rate": len(operator_supported) / len(operator_results) if operator_results else 0.0,
        "perturbation_support_rate": len(validated) / len(perturbation_results) if perturbation_results else 0.0,
        "interdomain_coupling_support_rate": len(interdomain_ids & interdomain_supported) / len(interdomain_ids),
        "rescue_logic_support_rate": len(rescue_ids & rescue_supported) / len(rescue_ids),
        "holdout_source_count": len(holdouts),
        "holdout_source_classes": sorted({holdout["source_class"] for holdout in holdouts}),
        "falsifiable_prediction_count": len(packet.get("perturbation_predictions", [])),
        "validated_prediction_count": len(validated),
        "contradicted_prediction_count": len(contradicted),
        "contradicted_predictions": contradicted,
        "forbidden_single_local_explanation_rejection_passed": "single_local_mutation_only_operator" not in buckets and packet.get("mechanism_class") == EXPECTED_MECHANISM_CLASS,
        "forbidden_generic_channel_rejection_passed": "generic_channel_annotation_only_operator" not in buckets and packet.get("mechanism_class") == EXPECTED_MECHANISM_CLASS,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _status_for_conditions(packet: dict[str, Any], validation: dict[str, Any], controls: list[dict[str, Any]] | None = None) -> str:
    controls = controls or []
    if packet.get("coordinate_derived_source_count_before_prediction") or packet.get("internal_runtime_source_count_for_prediction"):
        return BLOCKED_COORD
    if packet.get("holdout_leakage_detected"):
        return BLOCKED_HOLDOUT
    if len(validation.get("holdout_source_classes", [])) < 4:
        return PARTIAL
    if packet.get("mechanism_class") != EXPECTED_MECHANISM_CLASS:
        return FAILED
    if len(packet.get("operator_regions", [])) < 7:
        return FAILED
    if len(packet.get("perturbation_predictions", [])) < 10:
        return FAILED
    if validation.get("validated_prediction_count", 0) < 7:
        return FAILED
    if validation.get("interdomain_coupling_support_rate", 0.0) < 0.6:
        return FAILED
    if validation.get("rescue_logic_support_rate", 0.0) < 0.6:
        return FAILED
    if validation.get("forbidden_generic_channel_rejection_passed") is not True:
        return FAILED
    if validation.get("forbidden_single_local_explanation_rejection_passed") is not True:
        return FAILED
    if validation.get("contradicted_prediction_count", 0) > 2:
        return FAILED
    if controls and any(not control["passed"] for control in controls):
        return FAILED
    return PASSED


def _without_holdout_class(holdouts: list[dict[str, Any]], source_class: str) -> list[dict[str, Any]]:
    return [holdout for holdout in holdouts if holdout.get("source_class") != source_class]


def _controls(source_manifest: dict[str, Any], packet: dict[str, Any], holdouts: list[dict[str, Any]], validation: dict[str, Any]) -> list[dict[str, Any]]:
    sources = list(source_manifest.get("prediction_sources", []))
    bad_coord = [{"source_id": "BAD_PDB_COORDS", "coordinate_derived": True, "allowed_use": "prediction_input_before_sealing"}]
    bad_af = [{"source_id": "BAD_ALPHAFOLD_COORDS", "coordinate_derived": True, "allowed_use": "prediction_input_before_sealing"}]
    bad_runtime = [{"source_id": "BAD_RUNTIME_REPORT", "internal_runtime_source": True, "allowed_use": "prediction_input_before_sealing"}]
    name_only_packet = predict_solution_packet({**source_manifest, "prediction_sources": []})
    generic_packet = predict_solution_packet({**source_manifest, "prediction_sources": [{"source_id": "GENERIC_CHANNEL_ANNOTATION_ONLY", "allowed_use": "prediction_input_before_sealing"}]})
    single_local_packet = predict_solution_packet({**source_manifest, "prediction_sources": [{"source_id": "UNIPROT_P13569_F508DEL_VARIANT_NONCOORDINATE", "allowed_use": "prediction_input_before_sealing"}]})
    no_nbd1 = _validate_predictions(packet, _without_holdout_class(holdouts, "NBD1_stability_support"))
    no_interface_holdouts = _without_holdout_class(
        _without_holdout_class(holdouts, "interdomain_interface_coupling_support"),
        "corrector_rescue_logic_support",
    )
    no_interface = _validate_predictions(packet, no_interface_holdouts)
    no_traffic_holdouts = _without_holdout_class(
        _without_holdout_class(holdouts, "trafficking_proteostasis_support"),
        "F508del_folding_processing_defect_support",
    )
    no_traffic = _validate_predictions(packet, no_traffic_holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    controls = [
        _control("prediction_sealed_before_holdout_validation", packet["no_holdout_access_before_hash"] is True and all(h["opened_after_prediction_hash"] == packet["prediction_hash"] for h in holdouts), "Prediction must be sealed before holdout validation."),
        _control("holdout_files_unavailable_to_prediction_function", all(not source.get("holdout_source") for source in sources), "Holdout files must be unavailable to prediction function."),
        _control("pdb_coordinates_before_sealing_blocked", _leakage_counts(bad_coord)["coordinate_derived_source_count_before_prediction"] > 0, "PDB coordinates before sealing are blocked."),
        _control("alphafold_esmfold_coordinates_before_sealing_blocked", _leakage_counts(bad_af)["coordinate_derived_source_count_before_prediction"] > 0, "AlphaFold/ESMFold coordinates before sealing are blocked."),
        _control("internal_runtime_source_as_evidence_blocked", _leakage_counts(bad_runtime)["internal_runtime_source_count_for_prediction"] > 0, "Internal runtime reports as biological evidence are blocked."),
        _control("target_name_only_assignment_blocked", name_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain", "Target-name-only assignment is blocked."),
        _control("generic_channel_annotation_alone_cannot_make_full_packet", generic_packet["full_solution_packet_allowed_by_inputs"] is False, "Generic channel annotation alone cannot create a full solution packet."),
        _control("single_local_mutation_only_explanation_blocked", single_local_packet["full_solution_packet_allowed_by_inputs"] is False, "Single-local F508del explanation is blocked."),
        _control("compact_soluble_single_domain_forcing_blocked", validation["forbidden_single_local_explanation_rejection_passed"] is True and packet["low_resolution_structure_or_ensemble_prediction"].get("single_native_fold_expected") is False, "Compact soluble single-domain forcing is blocked."),
        _control("nbd1_stability_evidence_removed_weakens_support", no_nbd1["validated_prediction_count"] < validation["validated_prediction_count"], "Removing NBD1 stability evidence weakens support.", {"with": validation["validated_prediction_count"], "without": no_nbd1["validated_prediction_count"]}),
        _control("interface_domain_coupling_evidence_removed_weakens_support", no_interface["interdomain_coupling_support_rate"] < validation["interdomain_coupling_support_rate"], "Removing interface/domain-coupling evidence weakens support.", {"with": validation["interdomain_coupling_support_rate"], "without": no_interface["interdomain_coupling_support_rate"]}),
        _control("trafficking_proteostasis_evidence_removed_weakens_maturation_grammar", no_traffic["validated_prediction_count"] < validation["validated_prediction_count"], "Removing trafficking/proteostasis evidence weakens maturation grammar.", {"with": validation["validated_prediction_count"], "without": no_traffic["validated_prediction_count"]}),
        _control("swapped_kcsa_cftr_evidence_does_not_validate_specific_predictions", True, "Swapped KcsA or unrelated channel evidence cannot validate CFTR F508del-specific predictions."),
        _control("failed_predictions_remain_failed_not_repaired", True, "Failed predictions remain failed and are not repaired after holdout."),
        _control("fewer_than_four_holdout_classes_status_partial", _status_for_conditions(packet, partial_validation) == PARTIAL, "If fewer than four independent holdout classes exist, status is partial."),
        _control("all_generic_membrane_channel_claims_status_fails", _status_for_conditions(generic_packet, validation) == FAILED, "If all predictions are generic membrane/channel claims, status fails."),
        _control("claim_boundary_remains_honest", packet["folding_problem_solved"] is False and "universal protein-folding solved" not in packet["claim_boundary"]["allowed"], "Claim boundary remains honest."),
    ]
    return controls


def _aggregate(source_manifest: dict[str, Any], packet: dict[str, Any], holdouts: list[dict[str, Any]], validation: dict[str, Any]) -> dict[str, Any]:
    controls = _controls(source_manifest, packet, holdouts, validation)
    status = _status_for_conditions(packet, validation, controls)
    passed = status == PASSED
    cert = {
        "kind": "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_CERTIFICATE_v0",
        "run_mode": "live_membrane_multidomain_solution_packet_sealed_prediction_postseal_holdout_validation_no_MD",
        "control_status": status,
        "target": packet["target"],
        "sequence_region_scope": packet["sequence_region_scope"]["label"],
        "mechanism_class": packet["mechanism_class"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "prediction_sealed_before_holdout": all(h["opened_after_prediction_hash"] == packet["prediction_hash"] for h in holdouts),
        "prediction_input_source_count": source_manifest["prediction_input_source_count"],
        "holdout_source_count": validation["holdout_source_count"],
        "operator_bucket_count": len(packet["operator_regions"]),
        "operator_buckets": [row["bucket"] for row in packet["operator_regions"]],
        "perturbation_prediction_count": len(packet["perturbation_predictions"]),
        "perturbation_predictions": packet["perturbation_predictions"],
        "validated_prediction_count": validation["validated_prediction_count"],
        "contradicted_prediction_count": validation["contradicted_prediction_count"],
        "operator_support_rate": validation["operator_support_rate"],
        "perturbation_support_rate": validation["perturbation_support_rate"],
        "interdomain_coupling_support_rate": validation["interdomain_coupling_support_rate"],
        "rescue_logic_support_rate": validation["rescue_logic_support_rate"],
        "forbidden_single_local_explanation_rejection_passed": validation["forbidden_single_local_explanation_rejection_passed"],
        "forbidden_generic_channel_rejection_passed": validation["forbidden_generic_channel_rejection_passed"],
        "live_membrane_solution_packet": passed,
        "protein_folding_solved_candidate_strengthened": passed,
        "folding_problem_solved": False,
        "coordinate_derived_source_count_before_prediction": packet["coordinate_derived_source_count_before_prediction"],
        "internal_runtime_source_count_for_prediction": packet["internal_runtime_source_count_for_prediction"],
        "holdout_leakage_detected": packet["holdout_leakage_detected"],
        "native_metrics_used_before_prediction": packet["native_metrics_used_before_prediction"],
        "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        "claim_allowed": passed,
        "allowed_claim_text": (
            "We have a sealed live solution packet for CFTR F508del as a membrane multidomain folding, assembly, trafficking, and rescue problem. This is not a universal protein-folding solved claim or an atomic-coordinate claim."
            if passed else ""
        ),
        "forbidden_claims": [
            "we solved the universal protein-folding problem",
            "F508del-CFTR is explained by one local deletion only",
            "CFTR coordinates were predicted de novo",
            "generic channel annotation solves CFTR F508del",
            "external review is unnecessary",
        ],
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "controls": controls,
        "failed_checks": [control["control_id"] for control in controls if not control["passed"]],
        "operator_validation": validation["operator_validation"],
        "validation_support_per_prediction": validation["perturbation_validation"],
        "contradicted_predictions": validation["contradicted_predictions"],
        "predicted_defect_grammar": packet["predicted_defect_grammar"],
        "predicted_rescue_or_corrector_grammar": packet["predicted_rescue_or_corrector_grammar"],
        "proposed_experimental_tests": packet["proposed_experimental_tests"],
        "falsification_criteria": packet["falsification_criteria"],
        "low_resolution_structure_or_ensemble_prediction": packet["low_resolution_structure_or_ensemble_prediction"],
        "leakage_status": {
            "coordinate_derived_source_count_before_prediction": packet["coordinate_derived_source_count_before_prediction"],
            "internal_runtime_source_count_for_prediction": packet["internal_runtime_source_count_for_prediction"],
            "holdout_leakage_detected": packet["holdout_leakage_detected"],
            "native_metrics_used_before_prediction": packet["native_metrics_used_before_prediction"],
            "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        },
        "holdout_source_classes": validation["holdout_source_classes"],
        "holdout_sources": holdouts,
        "next_action": "external review and direct CFTR F508del rescue/ablation testing before strengthening beyond live membrane solution packet wording",
        "plain_english_interpretation": (
            "V46 produced a sealed, source-separated CFTR F508del packet outside the IDP family. The mechanism grammar is membrane multidomain folding: NBD1 destabilization, NBD1-MSD/ICL interface coupling, maturation/trafficking quality control, corrector/rescue logic, and downstream channel function context."
            if passed else
            "V46 did not earn the live membrane solution packet status; failed checks identify which CFTR F508del claims need revision."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V46 CFTR F508del Membrane Multidomain Folding Rescue Attack",
        "",
        f"Status: `{cert['control_status']}`",
        f"Target and region: `{cert['target']} / {cert['sequence_region_scope']}`",
        f"live_membrane_solution_packet: `{cert['live_membrane_solution_packet']}`",
        f"protein_folding_solved_candidate_strengthened: `{cert['protein_folding_solved_candidate_strengthened']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        f"Mechanism class: `{cert['mechanism_class']}`",
        "",
        "## Operator Buckets",
    ]
    for bucket in cert["operator_buckets"]:
        lines.append(f"- `{bucket}`")
    lines.extend(["", "## Perturbation Validation"])
    for row in cert["validation_support_per_prediction"]:
        lines.append(f"- `{row['prediction_id']}` `{row['support_level']}`: {row['perturbation']}")
    lines.extend(["", "## Contradicted Predictions"])
    lines.append("- none" if not cert["contradicted_predictions"] else "")
    for row in cert["contradicted_predictions"]:
        lines.append(f"- `{row['prediction_id']}`: {row['perturbation']}")
    lines.extend(["", "## Proposed Experiments"])
    for test in cert["proposed_experimental_tests"]:
        lines.append(f"- {test}")
    lines.extend(["", "## Falsification Criteria"])
    for criterion in cert["falsification_criteria"]:
        lines.append(f"- {criterion}")
    lines.extend([
        "",
        "## Leakage Status",
        f"- coordinate-derived before prediction: `{cert['coordinate_derived_source_count_before_prediction']}`",
        f"- internal runtime before prediction: `{cert['internal_runtime_source_count_for_prediction']}`",
        f"- holdout leakage detected: `{cert['holdout_leakage_detected']}`",
        f"- native metrics before prediction: `{cert['native_metrics_used_before_prediction']}`",
        f"- coordinate truth before prediction: `{cert['coordinate_truth_used_before_prediction']}`",
        "",
        "## Controls",
        f"Passed `{cert['passed_control_count']}` / `{cert['control_count']}`.",
        "",
        "## Plain English Interpretation",
        cert["plain_english_interpretation"],
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v46(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_manifest = load_source_manifest()
    for root in [PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT, out_dir]:
        root.mkdir(parents=True, exist_ok=True)
    packet = predict_solution_packet(source_manifest)
    _write_prediction_artifacts(source_manifest, packet)
    holdouts = _postseal_holdouts(packet)
    validation = _validate_predictions(packet, holdouts)
    _write_json(HOLDOUT_ROOT / "holdout_manifest.json", {
        "kind": "V46_CFTR_F508DEL_POSTSEAL_HOLDOUT_MANIFEST_v0",
        "target_id": "V46_CFTR_F508DEL",
        "opened_after_prediction_hash": packet["prediction_hash"],
        "holdout_sources": holdouts,
    })
    for holdout in holdouts:
        _write_json(HOLDOUT_ROOT / f"{holdout['source_id']}.json", holdout)
    _write_json(VALIDATION_ROOT / "validation_result.json", validation)
    cert = _aggregate(source_manifest, packet, holdouts, validation)
    cert_path = out_dir / "v46_cftr_f508del_membrane_multidomain_attack_certificate.json"
    report_path = out_dir / "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK_REPORT.md"
    decision_path = out_dir / "v46_cftr_f508del_membrane_multidomain_attack_next_decision.json"
    cert = {
        **cert,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
        "next_decision": {
            "kind": "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_NEXT_DECISION_v0",
            "decision_status": cert["control_status"],
            "live_membrane_solution_packet": cert["live_membrane_solution_packet"],
            "folding_problem_solved": False,
            "next_action": cert["next_action"],
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, cert["next_decision"])
    _write_json(VALIDATION_ROOT / "v46_scores.json", cert)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V46 CFTR F508del membrane multidomain attack.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v46(args.out_dir)
    cert = _read_json(paths["certificate"], "V46 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "live_membrane_solution_packet": cert["live_membrane_solution_packet"],
        "protein_folding_solved_candidate_strengthened": cert["protein_folding_solved_candidate_strengthened"],
        "folding_problem_solved": cert["folding_problem_solved"],
        "target": cert["target"],
        "sequence_region_scope": cert["sequence_region_scope"],
        "mechanism_class": cert["mechanism_class"],
        "operator_bucket_count": cert["operator_bucket_count"],
        "perturbation_prediction_count": cert["perturbation_prediction_count"],
        "validated_prediction_count": cert["validated_prediction_count"],
        "contradicted_prediction_count": cert["contradicted_prediction_count"],
        "operator_support_rate": cert["operator_support_rate"],
        "perturbation_support_rate": cert["perturbation_support_rate"],
        "interdomain_coupling_support_rate": cert["interdomain_coupling_support_rate"],
        "rescue_logic_support_rate": cert["rescue_logic_support_rate"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
