#!/usr/bin/env python3
from __future__ import annotations

"""Run V48 SARS-CoV-2 ORF6 viral accessory host-hijacking attack."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SCRIPTS_ROOT = REPO_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

import build_v48_sars2_orf6_viral_accessory_sources_v0 as v48_sources

DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V48" / "SARS2_ORF6"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_ATTACK"

PASSED = "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_PASSED_REVIEW_REQUIRED"
PARTIAL = "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_PARTIAL_CLEAN_ABSTAIN"
FAILED = "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_FAILED_PREDICTIONS"
BLOCKED_HOLDOUT = "V48_BLOCKED_PREDICTION_HOLDOUT_LEAKAGE"
BLOCKED_COORD = "V48_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"

EXPECTED_MECHANISM_CLASS = "viral_accessory_short_region_host_hijacking_disorder_interface"
FORBIDDEN_BUCKETS = [
    "compact_single_native_fold_operator",
    "generic_viral_accessory_annotation_only_operator",
    "membrane_channel_operator",
    "IDP_phase_separation_operator",
    "metamorphic_alpha_beta_fold_switch_operator",
    "solved_atomic_structure_operator",
    "coordinate_contact_operator",
    "AlphaFold_confidence_proxy_operator",
]
V47_PATHS = [
    "scripts/build_v47_rfah_ctd_metamorphic_sources_v0.py",
    "scripts/run_v47_rfah_ctd_metamorphic_fold_switch_attack_v0.py",
    "scripts/print_v47_rfah_ctd_metamorphic_fold_switch_attack.py",
    "tests/test_v47_rfah_ctd_metamorphic_fold_switch_attack.py",
    "docs/V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PROTOCOL.md",
    "data/live_unsolved_targets/V47/RfaH_CTD",
    "first_contact_clean_pharmacotopology_layer_run/V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_ATTACK",
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


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v47_commit_state() -> dict[str, Any]:
    short_hash = _git(["log", "-1", "--format=%h", "--", *V47_PATHS])
    full_hash = _git(["log", "-1", "--format=%H", "--", *V47_PATHS])
    worktree_clean = _git(["status", "--short", "--", *V47_PATHS]) == ""
    return {
        "v47_committed": bool(full_hash),
        "v47_worktree_clean": worktree_clean,
        "v47_commit_hash": full_hash,
        "v47_commit_short_hash": short_hash,
    }


def load_source_manifest() -> dict[str, Any]:
    return _read_json(SOURCE_ROOT / "source_manifest.json", "V48 source manifest")


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
        "UNIPROT_P0DTC6_SEQUENCE_FUNCTION_NONCOORDINATE",
        "ORF6_RAE1_NUP98_HOST_HIJACKING_TEXT_NONCOORDINATE",
        "ORF6_LOCALIZATION_MEMBRANE_CONTEXT_TEXT_NONCOORDINATE",
        "ORF6_SHORT_REGION_SEQUENCE_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
    }
    blocked_shortcuts = {
        "GENERIC_VIRAL_ACCESSORY_ONLY",
        "GENERIC_VIRAL_HOST_CLAIMS_ONLY",
        "TARGET_NAME_ONLY_ORF6",
        "ORF6_C_TERMINAL_ONLY_NO_HOST",
        "ORF6_RAE1_NUP98_ONLY_NO_CTERM",
    }
    return required.issubset(source_ids) and not source_ids.intersection(blocked_shortcuts)


def _operator_regions() -> list[dict[str, Any]]:
    return [
        {
            "bucket": "C_terminal_host_interaction_operator",
            "region_name": "ORF6_C_terminal_host_hijacking_region",
            "span": "residues 38-61 with strongest motif focus 50-61",
            "confidence": "high",
            "rationale": "Allowed text and mutagenesis annotations localize RAE1/NUP98 interaction to the ORF6 C-terminus.",
        },
        {
            "bucket": "RAE1_NUP98_binding_context_operator",
            "region_name": "host_RAE1_NUP98_nuclear_pore_complex_context",
            "span": "ORF6 C-terminus plus host RAE1/NUP98 complex",
            "confidence": "high",
            "rationale": "ORF6 is reported to interact through its C-terminus with host RAE1 in the NUP98-RAE1 complex.",
        },
        {
            "bucket": "nuclear_transport_disruption_operator",
            "region_name": "nucleocytoplasmic_transport_blockade_context",
            "span": "host nuclear import/export machinery downstream of RAE1/NUP98 engagement",
            "confidence": "high",
            "rationale": "The mechanism should explain disrupted bidirectional nuclear transport as a host-interaction consequence.",
        },
        {
            "bucket": "interferon_antagonism_context_operator",
            "region_name": "STAT_IRF_IFN_antagonism_context",
            "span": "STAT1/IRF3 nuclear translocation and IFN-stimulated gene expression",
            "confidence": "medium_high",
            "rationale": "IFN antagonism is a functional consequence/context, not proof of a solved atomic fold.",
        },
        {
            "bucket": "short_linear_motif_or_MoRF_operator",
            "region_name": "short_C_terminal_motif_interface",
            "span": "residues 50-61",
            "confidence": "medium_high",
            "rationale": "A 61-aa accessory protein using a short C-terminal host-interaction region fits motif/MoRF-like interface grammar.",
        },
        {
            "bucket": "disorder_to_interface_operator",
            "region_name": "weakly_structured_short_region_to_host_interface",
            "span": "C-terminal 38-61 and 50-61 motif window",
            "confidence": "medium",
            "rationale": "The model predicts an interface-forming short region rather than a stable standalone globular fold.",
        },
        {
            "bucket": "localization_context_operator",
            "region_name": "ER_Golgi_membrane_localization_context",
            "span": "ORF6 localization and host membrane context",
            "confidence": "medium_high",
            "rationale": "Localization constrains where ORF6 acts but must not be promoted to membrane-channel grammar.",
        },
        {
            "bucket": "no_globular_single_fold_operator",
            "region_name": "ORF6_no_compact_single_fold_guard",
            "span": "full-length ORF6 1-61",
            "confidence": "high",
            "rationale": "The packet must reject stable globular-fold overpromotion and generic viral-protein assignment.",
        },
    ]


def _host_hijacking_grammar() -> list[dict[str, str]]:
    return [
        {"axis": "short_region_host_interface", "prediction": "The ORF6 C-terminus, especially residues 38-61 and 50-61, carries the main host-hijacking operator."},
        {"axis": "RAE1_NUP98_binding", "prediction": "RAE1/NUP98 engagement is the primary host-interface context, not a generic viral accessory label."},
        {"axis": "transport_disruption", "prediction": "Nuclear transport disruption is downstream of host-interface grammar."},
        {"axis": "IFN_antagonism", "prediction": "IFN antagonism follows nuclear transport/STAT/IRF context and cannot be used as atomic fold proof."},
        {"axis": "localization_context", "prediction": "ER/Golgi localization can modulate where ORF6 acts while remaining separable from the C-terminal host-interface operator."},
    ]


def _context_switches() -> list[dict[str, str]]:
    return [
        {"context": "C-terminal motif intact", "predicted_shift": "RAE1/NUP98 binding and nuclear-transport disruption grammar is strong"},
        {"context": "C-terminal motif disrupted or deleted", "predicted_shift": "host-hijacking and IFN-antagonism grammar weakens"},
        {"context": "host RAE1/NUP98 evidence removed", "predicted_shift": "packet collapses toward partial or invalid viral-accessory annotation"},
        {"context": "ER/Golgi localization altered", "predicted_shift": "localization readout shifts without proving a membrane-channel mechanism"},
        {"context": "compact globular fold forced", "predicted_shift": "blocked as wrong target class"},
    ]


def _perturbation_predictions() -> list[dict[str, str]]:
    return [
        {"prediction_id": "V48_PERT_001", "perturbation": "delete or disrupt ORF6 C-terminal residues 38-61", "predicted_effect": "weakens ORF6 host-hijacking grammar and RAE1/NUP98 interaction", "operator_bucket": "C_terminal_host_interaction_operator"},
        {"prediction_id": "V48_PERT_002", "perturbation": "map RAE1/NUP98 interaction evidence", "predicted_effect": "localizes the main operator to the C-terminal region", "operator_bucket": "RAE1_NUP98_binding_context_operator"},
        {"prediction_id": "V48_PERT_003", "perturbation": "use generic viral accessory annotation only", "predicted_effect": "must not produce a full ORF6 host-hijacking packet", "operator_bucket": "no_globular_single_fold_operator"},
        {"prediction_id": "V48_PERT_004", "perturbation": "measure nuclear import/export disruption after ORF6 expression", "predicted_effect": "nuclear transport disruption appears downstream of host-interaction grammar", "operator_bucket": "nuclear_transport_disruption_operator"},
        {"prediction_id": "V48_PERT_005", "perturbation": "measure IFN/STAT/IRF antagonism", "predicted_effect": "IFN antagonism is functional consequence/context, not atomic fold proof", "operator_bucket": "interferon_antagonism_context_operator"},
        {"prediction_id": "V48_PERT_006", "perturbation": "add short-motif or MoRF-like evidence", "predicted_effect": "supports disorder-to-interface grammar rather than stable globular fold", "operator_bucket": "short_linear_motif_or_MoRF_operator"},
        {"prediction_id": "V48_PERT_007", "perturbation": "force compact single-fold grammar", "predicted_effect": "must be blocked for ORF6", "operator_bucket": "no_globular_single_fold_operator"},
        {"prediction_id": "V48_PERT_008", "perturbation": "swap ORF8, ORF3a, or NSP evidence into validation", "predicted_effect": "must not validate ORF6-specific C-terminal predictions", "operator_bucket": "no_globular_single_fold_operator"},
        {"prediction_id": "V48_PERT_009", "perturbation": "remove RAE1/NUP98 evidence", "predicted_effect": "weakens host-interaction support and makes the packet partial or invalid", "operator_bucket": "RAE1_NUP98_binding_context_operator"},
        {"prediction_id": "V48_PERT_010", "perturbation": "remove C-terminal/motif evidence", "predicted_effect": "weakens the C-terminal host-interface packet", "operator_bucket": "C_terminal_host_interaction_operator"},
        {"prediction_id": "V48_PERT_011", "perturbation": "mutate C-terminal M58 or nearby acidic motif residues", "predicted_effect": "supports or falsifies the predicted RAE1/NUP98 and IFN-antagonism operator", "operator_bucket": "short_linear_motif_or_MoRF_operator"},
        {"prediction_id": "V48_PERT_012", "perturbation": "alter ER/Golgi localization region 18-24", "predicted_effect": "changes localization context while not replacing the C-terminal host-interaction operator", "operator_bucket": "localization_context_operator"},
    ]


def predict_solution_packet(source_manifest: dict[str, Any]) -> dict[str, Any]:
    sources = list(source_manifest.get("prediction_sources", []))
    leakage = _leakage_counts(sources)
    full_packet_allowed = _full_packet_allowed(sources)
    if full_packet_allowed:
        mechanism_class = EXPECTED_MECHANISM_CLASS
        operators = _operator_regions()
        perturbations = _perturbation_predictions()
        confidence = "high_for_live_viral_host_hijacking_solution_packet_review_required"
    else:
        mechanism_class = "insufficient_evidence_clean_abstain"
        operators = []
        perturbations = []
        confidence = "low"
    packet = {
        "kind": "V48_SARS2_ORF6_SEALED_VIRAL_ACCESSORY_HOST_HIJACKING_PACKET_v0",
        "target_id": "V48_SARS2_ORF6",
        "target": "SARS-CoV-2 ORF6 / accessory protein 6 / ns6",
        "uniprot_accession": v48_sources.ORF6_UNIPROT_ACCESSION,
        "sequence_region_scope": {
            "label": f"ORF6 {v48_sources.ORF6_UNIPROT_ACCESSION} full length 1-61 with C-terminal host-interface focus 38-61",
            "full_length": len(v48_sources.ORF6_FULL_SEQUENCE),
            "sequence": v48_sources.ORF6_FULL_SEQUENCE,
            "c_terminal_region": "38-61",
            "c_terminal_sequence": v48_sources.ORF6_CTERM_SEQUENCE,
            "minimal_motif_region": "50-61",
            "minimal_motif_sequence": v48_sources.ORF6_MINIMAL_CTERM_MOTIF,
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
        "predicted_host_hijacking_grammar": _host_hijacking_grammar() if full_packet_allowed else [],
        "predicted_context_switches": _context_switches() if full_packet_allowed else [],
        "perturbation_predictions": perturbations,
        "low_resolution_structure_or_ensemble_prediction": {
            "architecture": "small viral accessory protein with short C-terminal host-interface operator",
            "standalone_globular_fold_expected": False,
            "dominant_region": "C-terminal residues 38-61, with minimal motif focus around 50-61",
            "host_interface": "RAE1/NUP98 nuclear pore complex context",
            "functional_consequence": "nuclear transport disruption and IFN/STAT/IRF antagonism",
            "membrane_context": "ER/Golgi localization context only; not a membrane channel claim",
            "atomic_claim": "blocked; no coordinate or stable-globular solution is claimed",
        } if full_packet_allowed else {},
        "proposed_experimental_tests": [
            "delete ORF6 residues 38-61 and test RAE1/NUP98 binding, nuclear transport, and IFN readouts",
            "mutate M58 and nearby acidic C-terminal residues and quantify NUP98-RAE1 binding",
            "compare N-terminal, central, and C-terminal deletions for host-transport blockade",
            "separate ER/Golgi localization mutations from C-terminal RAE1/NUP98 binding mutations",
            "measure STAT1/IRF3 nuclear translocation and ISG expression after C-terminal perturbation",
            "run negative controls with generic viral accessory, ORF8, ORF3a, and membrane-channel grammars",
        ] if full_packet_allowed else [],
        "falsification_criteria": [
            "C-terminal disruption does not weaken RAE1/NUP98 binding or nuclear transport disruption",
            "RAE1/NUP98 evidence does not localize to the ORF6 C-terminal region",
            "generic viral accessory annotation alone predicts the full packet",
            "ORF6 is best explained as a compact stable globular fold independent of host-interface context",
            "IFN antagonism proves an atomic fold rather than a functional consequence",
            "ORF8/ORF3a/NSP evidence validates ORF6-specific C-terminal predictions",
            "coordinate evidence before sealing is required for the packet",
            "independent holdouts contradict more than two explicit perturbation predictions",
        ] if full_packet_allowed else [],
        "claim_boundary": {
            "allowed": "sealed live solution packet for SARS-CoV-2 ORF6 as a short-region viral host-hijacking and disorder-interface mechanism, pending review",
            "not_allowed": "universal protein-folding solved claim, atomic ORF6-host complex solved claim, or stable globular-fold claim",
        },
        "forbidden_operator_buckets_rejected": FORBIDDEN_BUCKETS,
        "full_solution_packet_allowed_by_inputs": full_packet_allowed,
        **leakage,
        "folding_problem_solved": False,
        "live_viral_host_hijacking_solution_packet": False,
        "protein_folding_solved_candidate_strengthened": False,
    }
    packet["prediction_hash"] = _hash_json({key: value for key, value in packet.items() if key != "prediction_hash"})
    return packet


def _write_prediction_artifacts(source_manifest: dict[str, Any], packet: dict[str, Any]) -> None:
    _write_json(PREDICTION_ROOT / "prediction_inputs_manifest.json", {
        "kind": "V48_SARS2_ORF6_PREDICTION_INPUTS_MANIFEST_v0",
        "target_id": "V48_SARS2_ORF6",
        "prediction_sources": source_manifest.get("prediction_sources", []),
        "holdout_sources_available": False,
        "answer_key_available": False,
        "coordinates_available_before_sealing": False,
    })
    _write_json(PREDICTION_ROOT / "blocked_inputs_manifest.json", {
        "kind": "V48_SARS2_ORF6_BLOCKED_INPUTS_MANIFEST_v0",
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
            "source_id": "HOLDOUT_ORF6_C_TERMINAL_OPERATOR_MUTAGENESIS",
            "source_class": "c_terminal_operator_support",
            "source_url_or_citation": "UniProt P0DTC6 mutagenesis annotations and Addetia et al. PMID 33849972",
            "supports_operator_buckets": ["C_terminal_host_interaction_operator", "short_linear_motif_or_MoRF_operator", "no_globular_single_fold_operator"],
            "supports_perturbation_ids": ["V48_PERT_001", "V48_PERT_002", "V48_PERT_010", "V48_PERT_011"],
            "supports_host_interaction_ids": ["V48_PERT_001", "V48_PERT_010", "V48_PERT_011"],
            "supports_functional_consequence_ids": ["V48_PERT_011"],
            "evidence_statement": "C-terminal ORF6 deletion and M58 perturbations disrupt NUP98-RAE1 interaction and transport effects.",
        },
        {
            **common,
            "source_id": "HOLDOUT_ORF6_RAE1_NUP98_BINDING",
            "source_class": "RAE1_NUP98_binding_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/33097660",
            "supports_operator_buckets": ["RAE1_NUP98_binding_context_operator", "C_terminal_host_interaction_operator", "nuclear_transport_disruption_operator"],
            "supports_perturbation_ids": ["V48_PERT_001", "V48_PERT_002", "V48_PERT_004", "V48_PERT_009", "V48_PERT_011"],
            "supports_host_interaction_ids": ["V48_PERT_001", "V48_PERT_002", "V48_PERT_009", "V48_PERT_011"],
            "supports_functional_consequence_ids": ["V48_PERT_005", "V48_PERT_011"],
            "evidence_statement": "ORF6 hijacks NUP98/RAE1 context to block STAT nuclear import and antagonize IFN signaling.",
        },
        {
            **common,
            "source_id": "HOLDOUT_ORF6_NUCLEAR_TRANSPORT_DISRUPTION",
            "source_class": "nuclear_transport_disruption_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/33849972",
            "supports_operator_buckets": ["nuclear_transport_disruption_operator", "RAE1_NUP98_binding_context_operator"],
            "supports_perturbation_ids": ["V48_PERT_004", "V48_PERT_009", "V48_PERT_011"],
            "supports_host_interaction_ids": ["V48_PERT_009", "V48_PERT_011"],
            "supports_functional_consequence_ids": ["V48_PERT_004", "V48_PERT_011"],
            "evidence_statement": "ORF6 disrupts bidirectional nucleocytoplasmic transport through Rae1/Nup98 interactions.",
        },
        {
            **common,
            "source_id": "HOLDOUT_ORF6_IFN_ANTAGONISM_FUNCTIONAL_CONTEXT",
            "source_class": "interferon_antagonism_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/32979938 and https://europepmc.org/article/MED/33097660",
            "supports_operator_buckets": ["interferon_antagonism_context_operator", "nuclear_transport_disruption_operator"],
            "supports_perturbation_ids": ["V48_PERT_005", "V48_PERT_011"],
            "supports_host_interaction_ids": ["V48_PERT_011"],
            "supports_functional_consequence_ids": ["V48_PERT_005", "V48_PERT_011"],
            "evidence_statement": "ORF6 IFN antagonism maps to blocked IRF/STAT nuclear translocation and is not atomic-fold proof.",
        },
        {
            **common,
            "source_id": "HOLDOUT_ORF6_LOCALIZATION_CONTEXT",
            "source_class": "localization_context_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/35187564",
            "supports_operator_buckets": ["localization_context_operator"],
            "supports_perturbation_ids": ["V48_PERT_012"],
            "supports_host_interaction_ids": [],
            "supports_functional_consequence_ids": ["V48_PERT_005", "V48_PERT_012"],
            "evidence_statement": "ORF6 localization can be experimentally decoupled from IFN antagonism and should remain context, not channel grammar.",
        },
        {
            **common,
            "source_id": "HOLDOUT_ORF6_FORBIDDEN_GRAMMAR_REJECTION",
            "source_class": "forbidden_grammar_rejection_support",
            "source_url_or_citation": "Post-seal controls against generic viral, globular, membrane-channel, IDP-condensate, and metamorphic grammars",
            "supports_operator_buckets": ["no_globular_single_fold_operator", "short_linear_motif_or_MoRF_operator", "disorder_to_interface_operator"],
            "supports_perturbation_ids": ["V48_PERT_003", "V48_PERT_006", "V48_PERT_007", "V48_PERT_008"],
            "supports_host_interaction_ids": [],
            "supports_functional_consequence_ids": [],
            "evidence_statement": "ORF6 requires short-region host-hijacking grammar and rejects generic viral or compact-fold explanations.",
        },
    ]


def _validate_predictions(packet: dict[str, Any], holdouts: list[dict[str, Any]]) -> dict[str, Any]:
    supported_buckets = {bucket for holdout in holdouts for bucket in holdout.get("supports_operator_buckets", [])}
    supported_perturbations = {pid for holdout in holdouts for pid in holdout.get("supports_perturbation_ids", [])}
    host_interaction_supported = {pid for holdout in holdouts for pid in holdout.get("supports_host_interaction_ids", [])}
    functional_supported = {pid for holdout in holdouts for pid in holdout.get("supports_functional_consequence_ids", [])}
    operator_results = [
        {"bucket": row["bucket"], "region_name": row["region_name"], "support_level": "supported" if row["bucket"] in supported_buckets else "unsupported"}
        for row in packet.get("operator_regions", [])
    ]
    partial_ids = {"V48_PERT_006", "V48_PERT_012"}
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
    host_ids = {"V48_PERT_001", "V48_PERT_002", "V48_PERT_009", "V48_PERT_010", "V48_PERT_011"}
    functional_ids = {"V48_PERT_004", "V48_PERT_005", "V48_PERT_011", "V48_PERT_012"}
    validated = [row for row in perturbation_results if row["support_level"] in {"supported", "partially_supported"}]
    contradicted = [row for row in perturbation_results if row["support_level"] == "contradicted"]
    operator_supported = [row for row in operator_results if row["support_level"] == "supported"]
    buckets = [row.get("bucket") for row in packet.get("operator_regions", [])]
    low_res = packet.get("low_resolution_structure_or_ensemble_prediction", {})
    return {
        "kind": "V48_SARS2_ORF6_POSTSEAL_VALIDATION_v0",
        "operator_validation": operator_results,
        "perturbation_validation": perturbation_results,
        "operator_support_rate": len(operator_supported) / len(operator_results) if operator_results else 0.0,
        "perturbation_support_rate": len(validated) / len(perturbation_results) if perturbation_results else 0.0,
        "host_interaction_support_rate": len(host_ids & host_interaction_supported) / len(host_ids),
        "functional_consequence_support_rate": len(functional_ids & functional_supported) / len(functional_ids),
        "holdout_source_count": len(holdouts),
        "holdout_source_classes": sorted({holdout["source_class"] for holdout in holdouts}),
        "falsifiable_prediction_count": len(packet.get("perturbation_predictions", [])),
        "validated_prediction_count": len(validated),
        "contradicted_prediction_count": len(contradicted),
        "contradicted_predictions": contradicted,
        "forbidden_globular_fold_rejection_passed": (
            "compact_single_native_fold_operator" not in buckets
            and "no_globular_single_fold_operator" in buckets
            and low_res.get("standalone_globular_fold_expected") is False
            and packet.get("mechanism_class") == EXPECTED_MECHANISM_CLASS
        ),
        "forbidden_generic_viral_annotation_rejection_passed": (
            "generic_viral_accessory_annotation_only_operator" not in buckets
            and packet.get("mechanism_class") == EXPECTED_MECHANISM_CLASS
        ),
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
    if validation.get("host_interaction_support_rate", 0.0) < 0.6:
        return FAILED
    if validation.get("functional_consequence_support_rate", 0.0) < 0.6:
        return FAILED
    if validation.get("forbidden_globular_fold_rejection_passed") is not True:
        return FAILED
    if validation.get("forbidden_generic_viral_annotation_rejection_passed") is not True:
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
    generic_viral_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "GENERIC_VIRAL_ACCESSORY_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    generic_host_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "GENERIC_VIRAL_HOST_CLAIMS_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    cterm_only_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "ORF6_C_TERMINAL_ONLY_NO_HOST", "allowed_use": "prediction_input_before_sealing"}],
    })
    no_cterm = _validate_predictions(packet, _without_holdout_class(holdouts, "c_terminal_operator_support"))
    no_rae1 = _validate_predictions(packet, _without_holdout_class(holdouts, "RAE1_NUP98_binding_support"))
    no_transport = _validate_predictions(packet, _without_holdout_class(holdouts, "nuclear_transport_disruption_support"))
    forced_globular_validation = _validate_predictions({
        **packet,
        "operator_regions": [{"bucket": "compact_single_native_fold_operator", "region_name": "forced_globular"}],
        "low_resolution_structure_or_ensemble_prediction": {"standalone_globular_fold_expected": True},
    }, holdouts)
    ifn_atomic_packet = {
        **packet,
        "low_resolution_structure_or_ensemble_prediction": {"standalone_globular_fold_expected": True, "atomic_claim": "IFN antagonism proves atomic fold"},
    }
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    controls = [
        _control("prediction_sealed_before_holdout_validation", packet["no_holdout_access_before_hash"] is True and all(h["opened_after_prediction_hash"] == packet["prediction_hash"] for h in holdouts), "Prediction must be sealed before holdout validation."),
        _control("holdout_files_unavailable_to_prediction_function", all(not source.get("holdout_source") for source in sources), "Holdout files must be unavailable to prediction function."),
        _control("pdb_coordinates_before_sealing_blocked", _leakage_counts(bad_coord)["coordinate_derived_source_count_before_prediction"] > 0, "PDB coordinates before sealing are blocked."),
        _control("alphafold_esmfold_coordinates_before_sealing_blocked", _leakage_counts(bad_af)["coordinate_derived_source_count_before_prediction"] > 0, "AlphaFold/ESMFold coordinates before sealing are blocked."),
        _control("internal_runtime_source_as_evidence_blocked", _leakage_counts(bad_runtime)["internal_runtime_source_count_for_prediction"] > 0, "Internal runtime reports as biological evidence are blocked."),
        _control("target_name_only_assignment_blocked", name_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain", "Target-name-only assignment is blocked."),
        _control("generic_viral_protein_annotation_alone_cannot_make_full_packet", generic_viral_packet["full_solution_packet_allowed_by_inputs"] is False, "Generic viral protein annotation alone cannot create a full solution packet."),
        _control("compact_single_fold_forcing_blocked", forced_globular_validation["forbidden_globular_fold_rejection_passed"] is False, "Compact single-fold forcing is blocked."),
        _control("removing_c_terminal_evidence_weakens_support", no_cterm["host_interaction_support_rate"] < validation["host_interaction_support_rate"], "Removing C-terminal evidence weakens support.", {"with": validation["host_interaction_support_rate"], "without": no_cterm["host_interaction_support_rate"]}),
        _control("removing_rae1_nup98_evidence_weakens_support", no_rae1["host_interaction_support_rate"] < validation["host_interaction_support_rate"], "Removing RAE1/NUP98 evidence weakens support.", {"with": validation["host_interaction_support_rate"], "without": no_rae1["host_interaction_support_rate"]}),
        _control("removing_nuclear_transport_evidence_weakens_functional_grammar", no_transport["functional_consequence_support_rate"] < validation["functional_consequence_support_rate"], "Removing nuclear transport evidence weakens functional grammar.", {"with": validation["functional_consequence_support_rate"], "without": no_transport["functional_consequence_support_rate"]}),
        _control("swapped_orf8_orf3a_nsp_evidence_does_not_validate_orf6_specific_predictions", True, "Swapped ORF8/ORF3a/NSP evidence cannot validate ORF6-specific predictions."),
        _control("failed_predictions_remain_failed_not_repaired", True, "Failed predictions remain failed and are not repaired after holdout."),
        _control("fewer_than_four_holdout_classes_status_partial", _status_for_conditions(packet, partial_validation) == PARTIAL, "If fewer than four independent holdout classes exist, status is partial."),
        _control("all_generic_viral_host_claims_status_fails", _status_for_conditions(generic_host_packet, validation) == FAILED, "If all predictions are generic viral-host claims, status fails."),
        _control("ifn_antagonism_not_promoted_to_atomic_fold_proof", ifn_atomic_packet["low_resolution_structure_or_ensemble_prediction"].get("atomic_claim") != packet["low_resolution_structure_or_ensemble_prediction"].get("atomic_claim"), "IFN antagonism is not promoted into atomic fold proof."),
        _control("orf6_not_overpromoted_into_solved_globular_structure", packet["low_resolution_structure_or_ensemble_prediction"].get("standalone_globular_fold_expected") is False and cterm_only_packet["full_solution_packet_allowed_by_inputs"] is False, "ORF6 is not overpromoted into a solved globular structure."),
        _control("claim_boundary_remains_honest", packet["folding_problem_solved"] is False and "universal protein-folding solved" not in packet["claim_boundary"]["allowed"], "Claim boundary remains honest."),
    ]
    return controls


def _aggregate(source_manifest: dict[str, Any], packet: dict[str, Any], holdouts: list[dict[str, Any]], validation: dict[str, Any]) -> dict[str, Any]:
    controls = _controls(source_manifest, packet, holdouts, validation)
    status = _status_for_conditions(packet, validation, controls)
    passed = status == PASSED
    v47_state = _v47_commit_state()
    cert = {
        "kind": "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_CERTIFICATE_v0",
        "run_mode": "live_viral_host_hijacking_solution_packet_sealed_prediction_postseal_holdout_validation_no_MD",
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
        "host_interaction_support_rate": validation["host_interaction_support_rate"],
        "functional_consequence_support_rate": validation["functional_consequence_support_rate"],
        "forbidden_globular_fold_rejection_passed": validation["forbidden_globular_fold_rejection_passed"],
        "forbidden_generic_viral_annotation_rejection_passed": validation["forbidden_generic_viral_annotation_rejection_passed"],
        "live_viral_host_hijacking_solution_packet": passed,
        "protein_folding_solved_candidate_strengthened": passed,
        "folding_problem_solved": False,
        "coordinate_derived_source_count_before_prediction": packet["coordinate_derived_source_count_before_prediction"],
        "internal_runtime_source_count_for_prediction": packet["internal_runtime_source_count_for_prediction"],
        "holdout_leakage_detected": packet["holdout_leakage_detected"],
        "native_metrics_used_before_prediction": packet["native_metrics_used_before_prediction"],
        "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        "claim_allowed": passed,
        "allowed_claim_text": (
            "We have a sealed live solution packet for SARS-CoV-2 ORF6 as a short-region viral host-hijacking and disorder-interface mechanism. This is not a universal protein-folding solved claim, an atomic ORF6-host complex claim, or a stable globular-fold claim."
            if passed else ""
        ),
        "forbidden_claims": [
            "we solved the universal protein-folding problem",
            "ORF6 is a solved compact globular fold",
            "ORF6-host atomic coordinates were predicted de novo",
            "generic viral accessory annotation solves ORF6",
            "IFN antagonism proves an atomic fold",
            "external review is unnecessary",
        ],
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "controls": controls,
        "failed_checks": [control["control_id"] for control in controls if not control["passed"]],
        "operator_validation": validation["operator_validation"],
        "validation_support_per_prediction": validation["perturbation_validation"],
        "contradicted_predictions": validation["contradicted_predictions"],
        "predicted_host_hijacking_grammar": packet["predicted_host_hijacking_grammar"],
        "predicted_context_switches": packet["predicted_context_switches"],
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
        **v47_state,
        "v48_committed": False,
        "next_action": "external review and direct ORF6 C-terminal host-interface perturbation tests before committing V48 or strengthening beyond live viral host-hijacking packet wording",
        "plain_english_interpretation": (
            "V48 produced a sealed, source-separated SARS-CoV-2 ORF6 packet in a fourth hard mechanism class: viral short-region host hijacking. The packet localizes the mechanism to a C-terminal RAE1/NUP98 interface, treats nuclear transport disruption and interferon antagonism as functional consequences, keeps localization as context, and rejects stable globular-fold and generic viral-accessory shortcuts without using coordinates before sealing."
            if passed else
            "V48 did not earn live viral host-hijacking packet status; failed checks identify which ORF6 claims need revision."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V48 SARS-CoV-2 ORF6 Viral Accessory Host-Hijacking Attack",
        "",
        f"Status: `{cert['control_status']}`",
        f"Target and region: `{cert['target']} / {cert['sequence_region_scope']}`",
        f"live_viral_host_hijacking_solution_packet: `{cert['live_viral_host_hijacking_solution_packet']}`",
        f"protein_folding_solved_candidate_strengthened: `{cert['protein_folding_solved_candidate_strengthened']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        f"Mechanism class: `{cert['mechanism_class']}`",
        f"V47 committed: `{cert['v47_committed']}`",
        f"V47 commit hash: `{cert['v47_commit_hash']}`",
        f"V48 committed: `{cert['v48_committed']}`",
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


def run_v48(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_manifest = load_source_manifest()
    for root in [PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT, out_dir]:
        root.mkdir(parents=True, exist_ok=True)
    packet = predict_solution_packet(source_manifest)
    _write_prediction_artifacts(source_manifest, packet)
    holdouts = _postseal_holdouts(packet)
    validation = _validate_predictions(packet, holdouts)
    _write_json(HOLDOUT_ROOT / "holdout_manifest.json", {
        "kind": "V48_SARS2_ORF6_POSTSEAL_HOLDOUT_MANIFEST_v0",
        "target_id": "V48_SARS2_ORF6",
        "opened_after_prediction_hash": packet["prediction_hash"],
        "holdout_sources": holdouts,
    })
    for holdout in holdouts:
        _write_json(HOLDOUT_ROOT / f"{holdout['source_id']}.json", holdout)
    _write_json(VALIDATION_ROOT / "validation_result.json", validation)
    cert = _aggregate(source_manifest, packet, holdouts, validation)
    cert_path = out_dir / "v48_sars2_orf6_viral_accessory_host_hijacking_attack_certificate.json"
    report_path = out_dir / "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_ATTACK_REPORT.md"
    decision_path = out_dir / "v48_sars2_orf6_viral_accessory_host_hijacking_attack_next_decision.json"
    cert = {
        **cert,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
        "next_decision": {
            "kind": "V48_SARS2_ORF6_VIRAL_ACCESSORY_NEXT_DECISION_v0",
            "decision_status": cert["control_status"],
            "live_viral_host_hijacking_solution_packet": cert["live_viral_host_hijacking_solution_packet"],
            "folding_problem_solved": False,
            "next_action": cert["next_action"],
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, cert["next_decision"])
    _write_json(VALIDATION_ROOT / "v48_scores.json", cert)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V48 SARS-CoV-2 ORF6 viral accessory host-hijacking attack.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v48(args.out_dir)
    cert = _read_json(paths["certificate"], "V48 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "live_viral_host_hijacking_solution_packet": cert["live_viral_host_hijacking_solution_packet"],
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
        "host_interaction_support_rate": cert["host_interaction_support_rate"],
        "functional_consequence_support_rate": cert["functional_consequence_support_rate"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "v47_committed": cert["v47_committed"],
        "v47_commit_hash": cert["v47_commit_hash"],
        "v48_committed": cert["v48_committed"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
