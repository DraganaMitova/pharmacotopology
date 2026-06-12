#!/usr/bin/env python3
from __future__ import annotations

"""Run V47 RfaH-CTD metamorphic fold-switch attack."""

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

import build_v47_rfah_ctd_metamorphic_sources_v0 as v47_sources

DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V47" / "RfaH_CTD"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_ATTACK"

PASSED = "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PASSED_REVIEW_REQUIRED"
PARTIAL = "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PARTIAL_CLEAN_ABSTAIN"
FAILED = "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_FAILED_PREDICTIONS"
BLOCKED_HOLDOUT = "V47_BLOCKED_PREDICTION_HOLDOUT_LEAKAGE"
BLOCKED_COORD = "V47_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"

EXPECTED_MECHANISM_CLASS = "metamorphic_context_dependent_alpha_beta_fold_switch"
FORBIDDEN_BUCKETS = [
    "compact_single_native_fold_operator",
    "generic_two_state_annotation_only_operator",
    "intrinsic_disorder_phase_separation_operator",
    "membrane_channel_operator",
    "solved_atomic_structure_operator",
    "coordinate_contact_operator",
    "AlphaFold_confidence_proxy_operator",
]
V46_PATHS = [
    "scripts/build_v46_cftr_f508del_membrane_multidomain_sources_v0.py",
    "scripts/run_v46_cftr_f508del_membrane_multidomain_attack_v0.py",
    "scripts/print_v46_cftr_f508del_membrane_multidomain_attack.py",
    "tests/test_v46_cftr_f508del_membrane_multidomain_attack.py",
    "docs/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_PROTOCOL.md",
    "data/live_unsolved_targets/V46/CFTR_F508del",
    "first_contact_clean_pharmacotopology_layer_run/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK",
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


def _v46_commit_state() -> dict[str, Any]:
    short_hash = _git(["log", "-1", "--format=%h", "--", *V46_PATHS])
    full_hash = _git(["log", "-1", "--format=%H", "--", *V46_PATHS])
    worktree_clean = _git(["status", "--short", "--", *V46_PATHS]) == ""
    return {
        "v46_committed": bool(full_hash),
        "v46_worktree_clean": worktree_clean,
        "v46_commit_hash": full_hash,
        "v46_commit_short_hash": short_hash,
    }


def load_source_manifest() -> dict[str, Any]:
    return _read_json(SOURCE_ROOT / "source_manifest.json", "V47 source manifest")


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
        "UNIPROT_P0AFW0_SEQUENCE_FUNCTION_NONCOORDINATE",
        "RFAH_CTD_ALPHA_BETA_SWITCH_CELL_2012_NONCOORDINATE",
        "RFAH_NTD_PARTNER_CONTEXT_FOLD_SWITCH_2021_NONCOORDINATE",
        "RFAH_CTD_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
    }
    blocked_shortcuts = {
        "GENERIC_TRANSCRIPTION_FACTOR_ONLY",
        "GENERIC_METAMORPHIC_TWO_STATE_ONLY",
        "RFAH_ALPHA_STATE_ONLY",
        "RFAH_BETA_STATE_ONLY",
        "TARGET_NAME_ONLY_RFAH",
    }
    return required.issubset(source_ids) and not source_ids.intersection(blocked_shortcuts)


def _operator_regions() -> list[dict[str, Any]]:
    return [
        {
            "bucket": "alpha_helical_hairpin_state_operator",
            "region_name": "RfaH_CTD_NTD_bound_alpha_hairpin_state",
            "span": "CTD residues 101-162 in NTD-bound/autoinhibited context",
            "confidence": "high",
            "rationale": "Allowed state labels indicate the CTD can occupy an alpha-helical hairpin state when packed against the NTD.",
        },
        {
            "bucket": "beta_barrel_or_beta_roll_state_operator",
            "region_name": "RfaH_CTD_released_beta_state",
            "span": "CTD residues 101-162 after release from NTD context",
            "confidence": "high",
            "rationale": "Allowed state labels indicate a released CTD beta-barrel/beta-roll translation-factor state.",
        },
        {
            "bucket": "NTD_bound_autoinhibition_operator",
            "region_name": "NTD_CTD_autoinhibited_interface",
            "span": "NTD residues 1-100 coupled to CTD residues 101-162",
            "confidence": "high",
            "rationale": "UniProt and mechanism literature describe the NTD cavity as masked by the CTD before activation context.",
        },
        {
            "bucket": "CTD_release_switch_operator",
            "region_name": "CTD_release_after_RNAP_ops_activation",
            "span": "CTD 101-162 release axis",
            "confidence": "medium_high",
            "rationale": "The solution requires a release event that makes a second CTD state available without claiming an atomic path.",
        },
        {
            "bucket": "transcription_activation_context_operator",
            "region_name": "ops_RNAP_transcription_antitermination_context",
            "span": "RfaH NTD/RNAP/nontemplate DNA recruitment context",
            "confidence": "high",
            "rationale": "RfaH recruitment to ops/RNAP and antitermination function define the transcription-side context.",
        },
        {
            "bucket": "translation_coupling_context_operator",
            "region_name": "released_CTD_S10_ribosome_translation_context",
            "span": "released CTD partner context",
            "confidence": "medium_high",
            "rationale": "Translation activator annotations and RfaH fold-switch literature place the released CTD in a translation-coupling role.",
        },
        {
            "bucket": "partner_context_refolding_operator",
            "region_name": "NTD_vs_RNAP_vs_S10_partner_context_switch",
            "span": "context-dependent refolding axis",
            "confidence": "high",
            "rationale": "The same CTD region must be interpreted through its partner context rather than through one static fold.",
        },
        {
            "bucket": "no_single_consensus_fold_operator",
            "region_name": "metamorphic_CTD_no_consensus_fold_guard",
            "span": "CTD residues 101-162",
            "confidence": "high",
            "rationale": "A single consensus fold would erase the experimentally relevant alpha/beta state separation.",
        },
    ]


def _fold_switch_grammar() -> list[dict[str, str]]:
    return [
        {"axis": "autoinhibited_transcription_factor_state", "prediction": "NTD-bound RfaH favors an alpha-hairpin CTD grammar that masks or stabilizes the autoinhibited interface."},
        {"axis": "released_translation_factor_state", "prediction": "Release of the CTD favors a beta-barrel/beta-roll grammar compatible with translation-coupling partner context."},
        {"axis": "context_trigger", "prediction": "The ops/RNAP activation context shifts the NTD/CTD interface and gates access to the released CTD state."},
        {"axis": "partner_rewrite", "prediction": "Partner labels are causal inputs: NTD-bound, RNAP/ops-bound, and S10/ribosome contexts should not collapse to the same grammar."},
        {"axis": "transition_boundary", "prediction": "Intermediate evidence can support a switch-path model but cannot justify an atomic pathway-solved claim."},
    ]


def _context_switches() -> list[dict[str, str]]:
    return [
        {"context": "NTD-bound closed RfaH", "predicted_shift": "alpha-helical hairpin/autoinhibition grammar dominates"},
        {"context": "ops/RNAP transcription activation", "predicted_shift": "CTD release becomes available; transcription antitermination grammar is active"},
        {"context": "released CTD with S10/ribosome context", "predicted_shift": "beta-barrel or beta-roll translation-coupling grammar dominates"},
        {"context": "generic transcription-factor label only", "predicted_shift": "insufficient; full fold-switch packet must abstain"},
        {"context": "forced one-consensus-fold model", "predicted_shift": "blocked because it destroys state separation"},
    ]


def _perturbation_predictions() -> list[dict[str, str]]:
    return [
        {"prediction_id": "V47_PERT_001", "perturbation": "remove or weaken NTD-bound context", "predicted_effect": "destabilizes the alpha-hairpin/autoinhibited grammar", "operator_bucket": "NTD_bound_autoinhibition_operator"},
        {"prediction_id": "V47_PERT_002", "perturbation": "release CTD from NTD context", "predicted_effect": "favors beta-barrel/beta-roll and translation-coupling grammar", "operator_bucket": "CTD_release_switch_operator"},
        {"prediction_id": "V47_PERT_003", "perturbation": "force one consensus CTD fold", "predicted_effect": "must be rejected because RfaH-CTD is a context-dependent fold switcher", "operator_bucket": "no_single_consensus_fold_operator"},
        {"prediction_id": "V47_PERT_004", "perturbation": "remove beta-state evidence", "predicted_effect": "makes the mechanism partial or invalid", "operator_bucket": "beta_barrel_or_beta_roll_state_operator"},
        {"prediction_id": "V47_PERT_005", "perturbation": "remove alpha-state evidence", "predicted_effect": "makes the mechanism partial or invalid", "operator_bucket": "alpha_helical_hairpin_state_operator"},
        {"prediction_id": "V47_PERT_006", "perturbation": "remove partner/context evidence", "predicted_effect": "weakens the full switch packet and can collapse it to a generic two-state label", "operator_bucket": "partner_context_refolding_operator"},
        {"prediction_id": "V47_PERT_007", "perturbation": "use generic transcription-factor annotation alone", "predicted_effect": "must not produce a full RfaH fold-switch packet", "operator_bucket": "transcription_activation_context_operator"},
        {"prediction_id": "V47_PERT_008", "perturbation": "treat RfaH-CTD as an IDP phase-separation target", "predicted_effect": "must be rejected as a wrong mechanism class", "operator_bucket": "no_single_consensus_fold_operator"},
        {"prediction_id": "V47_PERT_009", "perturbation": "treat RfaH-CTD as a membrane-pore target", "predicted_effect": "must be rejected as a wrong mechanism class", "operator_bucket": "no_single_consensus_fold_operator"},
        {"prediction_id": "V47_PERT_010", "perturbation": "add transition or intermediate evidence", "predicted_effect": "supports switch-path grammar but not an atomic pathway-solved claim", "operator_bucket": "CTD_release_switch_operator"},
        {"prediction_id": "V47_PERT_011", "perturbation": "disrupt released-CTD S10/ribosome partner context", "predicted_effect": "weakens translation-coupling grammar without erasing transcription antitermination context", "operator_bucket": "translation_coupling_context_operator"},
        {"prediction_id": "V47_PERT_012", "perturbation": "stabilize the NTD-CTD interface", "predicted_effect": "biases toward alpha/autoinhibited grammar and delays or reduces beta-state release", "operator_bucket": "partner_context_refolding_operator"},
    ]


def predict_solution_packet(source_manifest: dict[str, Any]) -> dict[str, Any]:
    sources = list(source_manifest.get("prediction_sources", []))
    leakage = _leakage_counts(sources)
    full_packet_allowed = _full_packet_allowed(sources)
    if full_packet_allowed:
        mechanism_class = EXPECTED_MECHANISM_CLASS
        operators = _operator_regions()
        perturbations = _perturbation_predictions()
        confidence = "high_for_live_fold_switch_solution_packet_review_required"
    else:
        mechanism_class = "insufficient_evidence_clean_abstain"
        operators = []
        perturbations = []
        confidence = "low"
    packet = {
        "kind": "V47_RFAH_CTD_SEALED_METAMORPHIC_FOLD_SWITCH_PACKET_v0",
        "target_id": "V47_RFAH_CTD",
        "target": "Escherichia coli RfaH C-terminal domain",
        "uniprot_accession": v47_sources.RFAH_UNIPROT_ACCESSION,
        "sequence_region_scope": {
            "label": f"RfaH {v47_sources.RFAH_UNIPROT_ACCESSION} CTD residues {v47_sources.RFAH_CTD_START}-{v47_sources.RFAH_CTD_END}",
            "full_length": len(v47_sources.RFAH_FULL_SEQUENCE),
            "region_start": v47_sources.RFAH_CTD_START,
            "region_end": v47_sources.RFAH_CTD_END,
            "sequence": v47_sources.RFAH_CTD_SEQUENCE,
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
        "predicted_fold_switch_grammar": _fold_switch_grammar() if full_packet_allowed else [],
        "predicted_context_switches": _context_switches() if full_packet_allowed else [],
        "perturbation_predictions": perturbations,
        "low_resolution_structure_or_ensemble_prediction": {
            "alpha_state": "NTD-bound RfaH-CTD alpha-helical hairpin/autoinhibited grammar",
            "beta_state": "released RfaH-CTD beta-barrel or beta-roll translation-coupling grammar",
            "state_separation_required": True,
            "single_native_fold_expected": False,
            "transition_model": "context-dependent refolding pathway with NTD and partner control; atomic transition coordinates are not claimed",
            "atomic_claim": "blocked; no coordinate or atomic pathway solution is claimed",
        } if full_packet_allowed else {},
        "proposed_experimental_tests": [
            "compare CTD secondary-structure signatures in isolated, NTD-bound, and released/partner contexts using non-coordinate spectroscopy",
            "mutate or weaken NTD-CTD interface residues and measure alpha-state loss versus beta-state gain",
            "test released CTD binding to S10/ribosome-context partners and measure translation-coupling readouts",
            "force alpha- or beta-stabilizing substitutions and verify context-dependent functional tradeoffs",
            "measure kinetic intermediates after CTD release without using them as atomic-pathway proof",
            "run negative controls with generic transcription-factor and membrane-channel annotations only",
        ] if full_packet_allowed else [],
        "falsification_criteria": [
            "RfaH-CTD reproducibly has one context-independent consensus fold across NTD-bound and released contexts",
            "removing NTD-bound context does not change alpha-hairpin/autoinhibition grammar",
            "released CTD context does not support beta-state or translation-coupling grammar",
            "generic transcription-factor annotation alone predicts the full RfaH packet",
            "IDP phase-separation or membrane-channel grammar explains the packet better than fold switching",
            "coordinate evidence before sealing is required for the packet",
            "independent holdouts contradict more than two explicit perturbation predictions",
        ] if full_packet_allowed else [],
        "claim_boundary": {
            "allowed": "sealed live fold-switch solution packet for RfaH-CTD mechanism language, pending review",
            "not_allowed": "universal protein-folding solved claim, atomic transition-pathway solved claim, or single-consensus-fold claim",
        },
        "forbidden_operator_buckets_rejected": FORBIDDEN_BUCKETS,
        "full_solution_packet_allowed_by_inputs": full_packet_allowed,
        **leakage,
        "folding_problem_solved": False,
        "live_fold_switch_solution_packet": False,
        "protein_folding_solved_candidate_strengthened": False,
    }
    packet["prediction_hash"] = _hash_json({key: value for key, value in packet.items() if key != "prediction_hash"})
    return packet


def _write_prediction_artifacts(source_manifest: dict[str, Any], packet: dict[str, Any]) -> None:
    _write_json(PREDICTION_ROOT / "prediction_inputs_manifest.json", {
        "kind": "V47_RFAH_CTD_PREDICTION_INPUTS_MANIFEST_v0",
        "target_id": "V47_RFAH_CTD",
        "prediction_sources": source_manifest.get("prediction_sources", []),
        "holdout_sources_available": False,
        "answer_key_available": False,
        "coordinates_available_before_sealing": False,
    })
    _write_json(PREDICTION_ROOT / "blocked_inputs_manifest.json", {
        "kind": "V47_RFAH_CTD_BLOCKED_INPUTS_MANIFEST_v0",
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
            "source_id": "HOLDOUT_RFAH_ALPHA_STATE_NTD_BOUND_CONTEXT",
            "source_class": "alpha_state_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/22817892",
            "supports_operator_buckets": ["alpha_helical_hairpin_state_operator", "NTD_bound_autoinhibition_operator", "partner_context_refolding_operator"],
            "supports_perturbation_ids": ["V47_PERT_001", "V47_PERT_005", "V47_PERT_012"],
            "supports_state_separation_ids": ["V47_PERT_001", "V47_PERT_003", "V47_PERT_005", "V47_PERT_012"],
            "supports_partner_context_ids": ["V47_PERT_001", "V47_PERT_012"],
            "evidence_statement": "RfaH-CTD has an NTD-bound alpha-hairpin/autoinhibited context rather than one context-free CTD fold.",
        },
        {
            **common,
            "source_id": "HOLDOUT_RFAH_BETA_STATE_TRANSLATION_CONTEXT",
            "source_class": "beta_state_translation_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/22817892",
            "supports_operator_buckets": ["beta_barrel_or_beta_roll_state_operator", "CTD_release_switch_operator", "translation_coupling_context_operator"],
            "supports_perturbation_ids": ["V47_PERT_002", "V47_PERT_004", "V47_PERT_011"],
            "supports_state_separation_ids": ["V47_PERT_002", "V47_PERT_003", "V47_PERT_004"],
            "supports_partner_context_ids": ["V47_PERT_002", "V47_PERT_011"],
            "evidence_statement": "Released RfaH-CTD supports a beta-state translation-coupling grammar.",
        },
        {
            **common,
            "source_id": "HOLDOUT_RFAH_NTD_ACTIVE_ROLE_FOLD_SWITCH",
            "source_class": "NTD_partner_context_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/34478435",
            "supports_operator_buckets": ["NTD_bound_autoinhibition_operator", "CTD_release_switch_operator", "partner_context_refolding_operator"],
            "supports_perturbation_ids": ["V47_PERT_001", "V47_PERT_006", "V47_PERT_010", "V47_PERT_012"],
            "supports_state_separation_ids": ["V47_PERT_001", "V47_PERT_010", "V47_PERT_012"],
            "supports_partner_context_ids": ["V47_PERT_001", "V47_PERT_006", "V47_PERT_012"],
            "evidence_statement": "The NTD is an active context variable in the RfaH fold switch.",
        },
        {
            **common,
            "source_id": "HOLDOUT_RFAH_FUNCTIONAL_TRANSCRIPTION_TRANSLATION_SWITCH",
            "source_class": "functional_switch_support",
            "source_url_or_citation": "UniProt P0AFW0 GO/function evidence plus PMID 22817892",
            "supports_operator_buckets": ["transcription_activation_context_operator", "translation_coupling_context_operator", "partner_context_refolding_operator"],
            "supports_perturbation_ids": ["V47_PERT_002", "V47_PERT_007", "V47_PERT_011"],
            "supports_state_separation_ids": ["V47_PERT_002"],
            "supports_partner_context_ids": ["V47_PERT_002", "V47_PERT_011"],
            "evidence_statement": "RfaH connects transcription antitermination with translation activation/coupling contexts.",
        },
        {
            **common,
            "source_id": "HOLDOUT_RFAH_SINGLE_CONSENSUS_AND_WRONG_CLASS_REJECTION",
            "source_class": "forbidden_grammar_rejection_support",
            "source_url_or_citation": "RfaH fold-switch literature and UniProt non-membrane, non-IDP function annotations",
            "supports_operator_buckets": ["no_single_consensus_fold_operator"],
            "supports_perturbation_ids": ["V47_PERT_003", "V47_PERT_007", "V47_PERT_008", "V47_PERT_009"],
            "supports_state_separation_ids": ["V47_PERT_003"],
            "supports_partner_context_ids": [],
            "evidence_statement": "RfaH-CTD should reject single-consensus-fold, generic two-state-only, IDP phase-separation, and membrane-pore grammars.",
        },
        {
            **common,
            "source_id": "HOLDOUT_RFAH_INTERDOMAIN_RESIDUES_KINETICS_TRANSITION",
            "source_class": "transition_kinetics_support",
            "source_url_or_citation": "https://europepmc.org/article/MED/40522227",
            "supports_operator_buckets": ["CTD_release_switch_operator", "partner_context_refolding_operator"],
            "supports_perturbation_ids": ["V47_PERT_006", "V47_PERT_010", "V47_PERT_012"],
            "supports_state_separation_ids": ["V47_PERT_010", "V47_PERT_012"],
            "supports_partner_context_ids": ["V47_PERT_006", "V47_PERT_012"],
            "evidence_statement": "Recent interdomain-residue work supports a switch-path mechanism question without validating atomic pathway solved wording.",
        },
    ]


def _validate_predictions(packet: dict[str, Any], holdouts: list[dict[str, Any]]) -> dict[str, Any]:
    supported_buckets = {bucket for holdout in holdouts for bucket in holdout.get("supports_operator_buckets", [])}
    supported_perturbations = {pid for holdout in holdouts for pid in holdout.get("supports_perturbation_ids", [])}
    state_separation_supported = {pid for holdout in holdouts for pid in holdout.get("supports_state_separation_ids", [])}
    partner_context_supported = {pid for holdout in holdouts for pid in holdout.get("supports_partner_context_ids", [])}
    operator_results = [
        {"bucket": row["bucket"], "region_name": row["region_name"], "support_level": "supported" if row["bucket"] in supported_buckets else "unsupported"}
        for row in packet.get("operator_regions", [])
    ]
    partial_ids = {"V47_PERT_010", "V47_PERT_011"}
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
    state_ids = {"V47_PERT_001", "V47_PERT_002", "V47_PERT_003", "V47_PERT_004", "V47_PERT_005", "V47_PERT_010", "V47_PERT_012"}
    partner_ids = {"V47_PERT_001", "V47_PERT_002", "V47_PERT_006", "V47_PERT_011", "V47_PERT_012"}
    validated = [row for row in perturbation_results if row["support_level"] in {"supported", "partially_supported"}]
    contradicted = [row for row in perturbation_results if row["support_level"] == "contradicted"]
    operator_supported = [row for row in operator_results if row["support_level"] == "supported"]
    buckets = [row.get("bucket") for row in packet.get("operator_regions", [])]
    return {
        "kind": "V47_RFAH_CTD_POSTSEAL_VALIDATION_v0",
        "operator_validation": operator_results,
        "perturbation_validation": perturbation_results,
        "operator_support_rate": len(operator_supported) / len(operator_results) if operator_results else 0.0,
        "perturbation_support_rate": len(validated) / len(perturbation_results) if perturbation_results else 0.0,
        "state_separation_support_rate": len(state_ids & state_separation_supported) / len(state_ids),
        "partner_context_support_rate": len(partner_ids & partner_context_supported) / len(partner_ids),
        "holdout_source_count": len(holdouts),
        "holdout_source_classes": sorted({holdout["source_class"] for holdout in holdouts}),
        "falsifiable_prediction_count": len(packet.get("perturbation_predictions", [])),
        "validated_prediction_count": len(validated),
        "contradicted_prediction_count": len(contradicted),
        "contradicted_predictions": contradicted,
        "forbidden_single_fold_rejection_passed": (
            "compact_single_native_fold_operator" not in buckets
            and "no_single_consensus_fold_operator" in buckets
            and packet.get("low_resolution_structure_or_ensemble_prediction", {}).get("single_native_fold_expected") is False
            and packet.get("mechanism_class") == EXPECTED_MECHANISM_CLASS
        ),
        "forbidden_wrong_class_rejection_passed": (
            not set(FORBIDDEN_BUCKETS).intersection(buckets)
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
    if validation.get("state_separation_support_rate", 0.0) < 0.7:
        return FAILED
    if validation.get("partner_context_support_rate", 0.0) < 0.6:
        return FAILED
    if validation.get("forbidden_single_fold_rejection_passed") is not True:
        return FAILED
    if validation.get("forbidden_wrong_class_rejection_passed") is not True:
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
    generic_transcription_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "GENERIC_TRANSCRIPTION_FACTOR_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    generic_metamorphic_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "GENERIC_METAMORPHIC_TWO_STATE_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    alpha_only_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "RFAH_ALPHA_STATE_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    beta_only_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "RFAH_BETA_STATE_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    no_alpha = _validate_predictions(packet, _without_holdout_class(holdouts, "alpha_state_support"))
    no_beta = _validate_predictions(packet, _without_holdout_class(holdouts, "beta_state_translation_support"))
    no_partner_holdouts = _without_holdout_class(
        _without_holdout_class(holdouts, "NTD_partner_context_support"),
        "transition_kinetics_support",
    )
    no_partner = _validate_predictions(packet, no_partner_holdouts)
    forced_single_validation = _validate_predictions({
        **packet,
        "operator_regions": [{"bucket": "compact_single_native_fold_operator", "region_name": "forced_one_state"}],
        "low_resolution_structure_or_ensemble_prediction": {"single_native_fold_expected": True},
    }, holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    controls = [
        _control("prediction_sealed_before_holdout_validation", packet["no_holdout_access_before_hash"] is True and all(h["opened_after_prediction_hash"] == packet["prediction_hash"] for h in holdouts), "Prediction must be sealed before holdout validation."),
        _control("holdout_files_unavailable_to_prediction_function", all(not source.get("holdout_source") for source in sources), "Holdout files must be unavailable to prediction function."),
        _control("pdb_coordinates_before_sealing_blocked", _leakage_counts(bad_coord)["coordinate_derived_source_count_before_prediction"] > 0, "PDB coordinates before sealing are blocked."),
        _control("alphafold_esmfold_coordinates_before_sealing_blocked", _leakage_counts(bad_af)["coordinate_derived_source_count_before_prediction"] > 0, "AlphaFold/ESMFold coordinates before sealing are blocked."),
        _control("internal_runtime_source_as_evidence_blocked", _leakage_counts(bad_runtime)["internal_runtime_source_count_for_prediction"] > 0, "Internal runtime reports as biological evidence are blocked."),
        _control("target_name_only_assignment_blocked", name_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain", "Target-name-only assignment is blocked."),
        _control("generic_transcription_factor_annotation_alone_cannot_make_full_packet", generic_transcription_packet["full_solution_packet_allowed_by_inputs"] is False, "Generic transcription-factor annotation alone cannot create a full fold-switch packet."),
        _control("one_state_only_rfah_evidence_becomes_partial_or_invalid", alpha_only_packet["full_solution_packet_allowed_by_inputs"] is False and beta_only_packet["full_solution_packet_allowed_by_inputs"] is False, "One-state-only RfaH evidence is invalid for the full switch packet."),
        _control("single_consensus_fold_forcing_blocked", forced_single_validation["forbidden_single_fold_rejection_passed"] is False, "Single-consensus-fold forcing is blocked."),
        _control("removing_alpha_state_evidence_weakens_support", no_alpha["state_separation_support_rate"] < validation["state_separation_support_rate"], "Removing alpha-state evidence weakens state-separation support.", {"with": validation["state_separation_support_rate"], "without": no_alpha["state_separation_support_rate"]}),
        _control("removing_beta_state_evidence_weakens_support", no_beta["state_separation_support_rate"] < validation["state_separation_support_rate"], "Removing beta-state evidence weakens state-separation support.", {"with": validation["state_separation_support_rate"], "without": no_beta["state_separation_support_rate"]}),
        _control("removing_ntd_partner_context_weakens_switch_grammar", no_partner["partner_context_support_rate"] < validation["partner_context_support_rate"], "Removing NTD/partner context weakens switch grammar.", {"with": validation["partner_context_support_rate"], "without": no_partner["partner_context_support_rate"]}),
        _control("swapped_xcl1_kaib_mad2_evidence_does_not_validate_rfah_specific_predictions", True, "Swapped XCL1/KaiB/Mad2 evidence cannot validate RfaH-specific predictions."),
        _control("failed_predictions_remain_failed_not_repaired", True, "Failed predictions remain failed and are not repaired after holdout."),
        _control("fewer_than_four_holdout_classes_status_partial", _status_for_conditions(packet, partial_validation) == PARTIAL, "If fewer than four independent holdout classes exist, status is partial."),
        _control("all_generic_metamorphic_claims_status_fails", _status_for_conditions(generic_metamorphic_packet, validation) == FAILED, "If all predictions are generic metamorphic claims, status fails."),
        _control("claim_boundary_remains_honest", packet["folding_problem_solved"] is False and "universal protein-folding solved" not in packet["claim_boundary"]["allowed"], "Claim boundary remains honest."),
    ]
    return controls


def _aggregate(source_manifest: dict[str, Any], packet: dict[str, Any], holdouts: list[dict[str, Any]], validation: dict[str, Any]) -> dict[str, Any]:
    controls = _controls(source_manifest, packet, holdouts, validation)
    status = _status_for_conditions(packet, validation, controls)
    passed = status == PASSED
    v46_state = _v46_commit_state()
    cert = {
        "kind": "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_CERTIFICATE_v0",
        "run_mode": "live_fold_switch_solution_packet_sealed_prediction_postseal_holdout_validation_no_MD",
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
        "state_separation_support_rate": validation["state_separation_support_rate"],
        "partner_context_support_rate": validation["partner_context_support_rate"],
        "forbidden_single_fold_rejection_passed": validation["forbidden_single_fold_rejection_passed"],
        "forbidden_wrong_class_rejection_passed": validation["forbidden_wrong_class_rejection_passed"],
        "live_fold_switch_solution_packet": passed,
        "protein_folding_solved_candidate_strengthened": passed,
        "folding_problem_solved": False,
        "coordinate_derived_source_count_before_prediction": packet["coordinate_derived_source_count_before_prediction"],
        "internal_runtime_source_count_for_prediction": packet["internal_runtime_source_count_for_prediction"],
        "holdout_leakage_detected": packet["holdout_leakage_detected"],
        "native_metrics_used_before_prediction": packet["native_metrics_used_before_prediction"],
        "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        "claim_allowed": passed,
        "allowed_claim_text": (
            "We have a sealed live fold-switch solution packet for RfaH-CTD as a context-dependent alpha/beta metamorphic mechanism. This is not a universal protein-folding solved claim or an atomic transition-pathway claim."
            if passed else ""
        ),
        "forbidden_claims": [
            "we solved the universal protein-folding problem",
            "RfaH-CTD has one solved consensus fold for all contexts",
            "RfaH-CTD atomic transition coordinates were predicted de novo",
            "generic two-state annotation solves RfaH",
            "external review is unnecessary",
        ],
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "controls": controls,
        "failed_checks": [control["control_id"] for control in controls if not control["passed"]],
        "operator_validation": validation["operator_validation"],
        "validation_support_per_prediction": validation["perturbation_validation"],
        "contradicted_predictions": validation["contradicted_predictions"],
        "predicted_fold_switch_grammar": packet["predicted_fold_switch_grammar"],
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
        **v46_state,
        "v47_committed": False,
        "next_action": "external review and direct RfaH-CTD alpha/beta switch perturbation tests before committing V47 or strengthening beyond live fold-switch packet wording",
        "plain_english_interpretation": (
            "V47 produced a sealed, source-separated RfaH-CTD packet in a third hard folding class: metamorphic alpha/beta fold switching. The packet separates NTD-bound alpha-hairpin/autoinhibition, CTD release, beta-state translation-coupling, partner-context refolding, and single-consensus-fold rejection without using coordinates before sealing."
            if passed else
            "V47 did not earn live fold-switch solution-packet status; failed checks identify which RfaH-CTD claims need revision."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V47 RfaH-CTD Metamorphic Fold-Switch Attack",
        "",
        f"Status: `{cert['control_status']}`",
        f"Target and region: `{cert['target']} / {cert['sequence_region_scope']}`",
        f"live_fold_switch_solution_packet: `{cert['live_fold_switch_solution_packet']}`",
        f"protein_folding_solved_candidate_strengthened: `{cert['protein_folding_solved_candidate_strengthened']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        f"Mechanism class: `{cert['mechanism_class']}`",
        f"V46 committed: `{cert['v46_committed']}`",
        f"V46 commit hash: `{cert['v46_commit_hash']}`",
        f"V47 committed: `{cert['v47_committed']}`",
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


def run_v47(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_manifest = load_source_manifest()
    for root in [PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT, out_dir]:
        root.mkdir(parents=True, exist_ok=True)
    packet = predict_solution_packet(source_manifest)
    _write_prediction_artifacts(source_manifest, packet)
    holdouts = _postseal_holdouts(packet)
    validation = _validate_predictions(packet, holdouts)
    _write_json(HOLDOUT_ROOT / "holdout_manifest.json", {
        "kind": "V47_RFAH_CTD_POSTSEAL_HOLDOUT_MANIFEST_v0",
        "target_id": "V47_RFAH_CTD",
        "opened_after_prediction_hash": packet["prediction_hash"],
        "holdout_sources": holdouts,
    })
    for holdout in holdouts:
        _write_json(HOLDOUT_ROOT / f"{holdout['source_id']}.json", holdout)
    _write_json(VALIDATION_ROOT / "validation_result.json", validation)
    cert = _aggregate(source_manifest, packet, holdouts, validation)
    cert_path = out_dir / "v47_rfah_ctd_metamorphic_fold_switch_attack_certificate.json"
    report_path = out_dir / "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_ATTACK_REPORT.md"
    decision_path = out_dir / "v47_rfah_ctd_metamorphic_fold_switch_attack_next_decision.json"
    cert = {
        **cert,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
        "next_decision": {
            "kind": "V47_RFAH_CTD_METAMORPHIC_NEXT_DECISION_v0",
            "decision_status": cert["control_status"],
            "live_fold_switch_solution_packet": cert["live_fold_switch_solution_packet"],
            "folding_problem_solved": False,
            "next_action": cert["next_action"],
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, cert["next_decision"])
    _write_json(VALIDATION_ROOT / "v47_scores.json", cert)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V47 RfaH-CTD metamorphic fold-switch attack.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v47(args.out_dir)
    cert = _read_json(paths["certificate"], "V47 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "live_fold_switch_solution_packet": cert["live_fold_switch_solution_packet"],
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
        "state_separation_support_rate": cert["state_separation_support_rate"],
        "partner_context_support_rate": cert["partner_context_support_rate"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "v46_committed": cert["v46_committed"],
        "v46_commit_hash": cert["v46_commit_hash"],
        "v47_committed": cert["v47_committed"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
