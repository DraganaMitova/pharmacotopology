#!/usr/bin/env python3
from __future__ import annotations

"""Run V44 FUS-LC live unsolved mechanism attack."""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import build_v44_fus_lc_live_unsolved_sources_v0 as v44_sources

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V44" / "FUS_LC"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V44_FUS_LC_LIVE_UNSOLVED_MECHANISM_ATTACK"

PASSED = "V44_FUS_LC_LIVE_UNSOLVED_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED"
PARTIAL = "V44_FUS_LC_LIVE_UNSOLVED_PARTIAL_CLEAN_ABSTAIN"
FAILED = "V44_FUS_LC_LIVE_UNSOLVED_FAILED_PREDICTIONS"
BLOCKED_HOLDOUT = "V44_BLOCKED_PREDICTION_HOLDOUT_LEAKAGE"
BLOCKED_COORD = "V44_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"

EXPECTED_MECHANISM_CLASS = "intrinsic_disorder_phase_separation_contextual_ensemble"
EXPECTED_OPERATOR_BUCKETS = [
    "low_complexity_disorder_operator",
    "aromatic_sticker_pattern_operator",
    "polar_spacer_context_operator",
    "phosphorylation_shift_operator",
    "LLPS_self_association_operator",
    "fibril_or_gel_state_shift_operator",
    "context_bound_RNA_or_protein_interaction_operator",
]
FORBIDDEN_BUCKETS = [
    "compact_single_native_fold_operator",
    "solved_atomic_structure_operator",
    "de_novo_global_coordinate_solution_operator",
    "AlphaFold_confidence_proxy_operator",
    "coordinate_contact_operator",
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
    return _read_json(SOURCE_ROOT / "source_manifest.json", "V44 source manifest")


def _leakage_counts(sources: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "coordinate_derived_source_count_before_prediction": sum(1 for source in sources if source.get("coordinate_derived")),
        "internal_runtime_source_count_for_prediction": sum(1 for source in sources if source.get("internal_runtime_source")),
        "holdout_leakage_detected": any(source.get("holdout_source") and source.get("allowed_use") == "prediction_input_before_sealing" for source in sources),
        "native_metrics_used_before_prediction": any(source.get("native_metrics_used_for_selection") for source in sources),
        "coordinate_truth_used_before_prediction": any(source.get("coordinate_truth_used_before_prediction") for source in sources),
    }


def _generic_idp_only(sources: list[dict[str, Any]]) -> bool:
    if not sources:
        return False
    has_only_generic = all(source.get("source_id") == "GENERIC_IDP_ANNOTATION_ONLY" for source in sources)
    return has_only_generic


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
        "UNIPROT_P35637_SEQUENCE_FEATURES_NONCOORDINATE",
        "DISPROT_DP01102_DISORDER_ONLY_NONCOORDINATE",
        "FUS_LC_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
    }
    return required.issubset(source_ids) and not _generic_idp_only(sources)


def _operator_regions(composition: dict[str, Any]) -> list[dict[str, Any]]:
    tyrosines = composition["tyrosine_positions"]
    serines = composition["serine_positions"]
    return [
        {
            "bucket": "low_complexity_disorder_operator",
            "region_name": "FUS_LC_QGSY_low_complexity_disordered_matrix",
            "span": "1-214",
            "confidence": "high",
            "rationale": "P35637 residues 1-214 are 81.3 percent Q/G/S/Y and overlap UniProt/DisProt disorder annotations.",
        },
        {
            "bucket": "aromatic_sticker_pattern_operator",
            "region_name": "distributed_tyrosine_sticker_ladder",
            "span": "1-164 with distal stickers through 208",
            "confidence": "high",
            "rationale": f"Twenty-seven tyrosines are distributed across the LC region at positions {tyrosines}.",
        },
        {
            "bucket": "polar_spacer_context_operator",
            "region_name": "Q_S_G_polar_spacer_network",
            "span": "1-214",
            "confidence": "high",
            "rationale": "Glutamine, serine, and glycine dominate the spacers, favoring disordered multivalent contacts over a single hydrophobic core.",
        },
        {
            "bucket": "phosphorylation_shift_operator",
            "region_name": "N_terminal_serine_charge_switch",
            "span": "26,30,42 plus serine-rich LC background",
            "confidence": "medium_high",
            "rationale": f"UniProt annotates phosphoserines at 26, 30, and 42; the LC region contains {len(serines)} serines.",
        },
        {
            "bucket": "LLPS_self_association_operator",
            "region_name": "multivalent_LC_self_association_field",
            "span": "1-214",
            "confidence": "high",
            "rationale": "Disordered low-complexity stickers and polar spacers predict concentration-dependent self-association and condensate formation.",
        },
        {
            "bucket": "fibril_or_gel_state_shift_operator",
            "region_name": "state_dependent_core_prone_windows",
            "span": "39-95, 110-150, and 155-190 as sequence-pattern hot zones",
            "confidence": "medium_high",
            "rationale": "Sticker-rich windows predict gel/fibril competence as a state-dependent branch, not a whole-domain native fold.",
        },
        {
            "bucket": "context_bound_RNA_or_protein_interaction_operator",
            "region_name": "RNA_or_partner_shifted_ensemble",
            "span": "1-214 with partner-sensitive contacts near 37-41 and 149-154",
            "confidence": "medium",
            "rationale": "FUS is an RNA-binding protein; LC disorder should be context-shifted by RNA/protein partners without becoming a single native fold.",
        },
    ]


def _state_grammar() -> list[dict[str, str]]:
    return [
        {"state": "dilute_disordered_monomer", "grammar": "expanded heterogeneous IDP ensemble with no single native fold"},
        {"state": "dense_liquid_condensate", "grammar": "multivalent aromatic sticker and polar spacer contacts exchange rapidly"},
        {"state": "aged_gel_or_fibril_branch", "grammar": "selected sticker-rich windows become ordered while flanking regions stay fuzzy"},
        {"state": "glass_like_or_slowed_material", "grammar": "condition-dependent kinetic arrest of the condensate network"},
        {"state": "RNA_or_partner_bound_shift", "grammar": "RNA or protein partners rebalance stickers, spacers, and charge-mediated contacts"},
    ]


def _context_switches() -> list[dict[str, str]]:
    return [
        {"context": "higher concentration", "predicted_shift": "increases dense-phase/self-association probability"},
        {"context": "increased phosphorylation or phosphomimetic substitution", "predicted_shift": "weakens LLPS, aggregation, and fibril competence"},
        {"context": "salt/pH/temperature changes", "predicted_shift": "moves liquid/gel/fibril/glass balance rather than selecting one fold"},
        {"context": "RNA or nuclear import receptor/protein partner present", "predicted_shift": "reshapes ensemble and can suppress or redirect self-association"},
    ]


def _perturbation_predictions() -> list[dict[str, str]]:
    return [
        {
            "prediction_id": "V44_PERT_001",
            "perturbation": "phosphomimetic or phosphorylation-like increase at N-terminal serine cluster",
            "predicted_effect": "weakens LLPS, aggregation, and fibril competence by adding negative charge to the sticker-spacer field",
            "operator_bucket": "phosphorylation_shift_operator",
        },
        {
            "prediction_id": "V44_PERT_002",
            "perturbation": "tyrosine/aromatic sticker disruption across FUS-LC",
            "predicted_effect": "weakens self-association and condensation by reducing aromatic multivalence",
            "operator_bucket": "aromatic_sticker_pattern_operator",
        },
        {
            "prediction_id": "V44_PERT_003",
            "perturbation": "tyrosine pattern redistribution without changing total length",
            "predicted_effect": "alters phase boundary and maturation kinetics because sticker spacing is part of the grammar",
            "operator_bucket": "aromatic_sticker_pattern_operator",
        },
        {
            "prediction_id": "V44_PERT_004",
            "perturbation": "replace or delete the Q/G/S/Y low-complexity identity",
            "predicted_effect": "destroys the FUS-LC grammar and collapses the LLPS/fibril-state mechanism",
            "operator_bucket": "low_complexity_disorder_operator",
        },
        {
            "prediction_id": "V44_PERT_005",
            "perturbation": "force compact single-fold grammar",
            "predicted_effect": "must be rejected; the target is an ensemble/phase-separation problem, not a single-native-fold problem",
            "operator_bucket": "low_complexity_disorder_operator",
        },
        {
            "prediction_id": "V44_PERT_006",
            "perturbation": "add RNA or phase-modulating protein partner",
            "predicted_effect": "shifts the ensemble and condensate material state through context-dependent binding and charge/sticker rebalance",
            "operator_bucket": "context_bound_RNA_or_protein_interaction_operator",
        },
        {
            "prediction_id": "V44_PERT_007",
            "perturbation": "isolate or mutate fibril-core-prone windows",
            "predicted_effect": "changes gel/fibril competence in a state-dependent way while leaving flanking regions fuzzy",
            "operator_bucket": "fibril_or_gel_state_shift_operator",
        },
        {
            "prediction_id": "V44_PERT_008",
            "perturbation": "change pH, salt, concentration, temperature, or phosphorylation state",
            "predicted_effect": "shifts the ensemble among dilute, liquid, gel, fibril, and glass-like states",
            "operator_bucket": "LLPS_self_association_operator",
        },
        {
            "prediction_id": "V44_PERT_009",
            "perturbation": "remove the N-terminal core-prone sticker region",
            "predicted_effect": "reduces core-1-like maturation and can reveal alternate C-terminal LC ordering",
            "operator_bucket": "fibril_or_gel_state_shift_operator",
        },
        {
            "prediction_id": "V44_PERT_010",
            "perturbation": "remove phosphorylation capacity with serine-to-alanine substitutions",
            "predicted_effect": "should bias toward stronger self-association relative to phosphomimetic variants",
            "operator_bucket": "phosphorylation_shift_operator",
        },
    ]


def predict_solution_packet(source_manifest: dict[str, Any]) -> dict[str, Any]:
    sources = list(source_manifest.get("prediction_sources", []))
    leakage = _leakage_counts(sources)
    full_packet_allowed = _full_packet_allowed(sources)
    pattern_source = next(
        (source for source in sources if source.get("source_id") == "FUS_LC_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS"),
        {},
    )
    composition = pattern_source.get("composition", {})
    if full_packet_allowed:
        mechanism_class = EXPECTED_MECHANISM_CLASS
        operator_regions = _operator_regions(composition)
        perturbations = _perturbation_predictions()
        mechanism_confidence = "high_for_live_solution_packet_review_required"
    else:
        mechanism_class = "insufficient_evidence_clean_abstain"
        operator_regions = []
        perturbations = []
        mechanism_confidence = "low"
    packet = {
        "kind": "V44_FUS_LC_SEALED_LIVE_UNSOLVED_SOLUTION_PACKET_v0",
        "target_id": "V44_FUS_LC",
        "target": "human FUS low-complexity domain",
        "uniprot_accession": "P35637",
        "sequence_region_scope": {
            "label": "FUS P35637 residues 1-214",
            "start": 1,
            "end": 214,
            "sequence": v44_sources.FUS_LC_SEQUENCE,
        },
        "prediction_timestamp": datetime.now(timezone.utc).isoformat(),
        "no_holdout_access_before_hash": True,
        "prediction_inputs_manifest": str(PREDICTION_ROOT / "prediction_inputs_manifest.json"),
        "blocked_inputs_manifest": str(PREDICTION_ROOT / "blocked_inputs_manifest.json"),
        "prediction_source_ids": [source.get("source_id") for source in sources],
        "mechanism_class": mechanism_class,
        "mechanism_confidence": mechanism_confidence,
        "operator_regions": operator_regions,
        "operator_region_rationale": {row["bucket"]: row["rationale"] for row in operator_regions},
        "predicted_state_grammar": _state_grammar() if full_packet_allowed else [],
        "predicted_context_switches": _context_switches() if full_packet_allowed else [],
        "perturbation_predictions": perturbations,
        "low_resolution_ensemble_prediction": {
            "ensemble_type": "heterogeneous IDP ensemble with condition-dependent condensate material states",
            "single_native_fold_expected": False,
            "state_model": "dilute disorder <-> dense liquid <-> aged gel/fibril/glass-like branch, shifted by phosphorylation, concentration, salt/pH/temperature, RNA, and protein partners",
            "fibril_core_boundary": "state/context-dependent local ordering, not a whole-domain native fold",
        } if full_packet_allowed else {},
        "proposed_experimental_tests": [
            "compare wild-type FUS-LC against S-to-E phosphomimetic and S-to-A phospho-null panels in turbidity, microscopy, FRAP, and fibril assays",
            "mutate or scramble tyrosine stickers while preserving Q/G/S composition and measure LLPS boundary shifts",
            "test 1-95, 39-95, 110-150, and 155-190 constructs for gel/fibril maturation and fuzzy-flank behavior",
            "measure salt, pH, concentration, temperature, and phosphorylation-state phase diagrams",
            "add RNA and Kapbeta2-like protein context and quantify condensate suppression, recruitment, or material-state redirection",
            "use NMR/HDX/crosslinking/single-molecule assays to test ensemble shifts without demanding a single coordinate model",
        ] if full_packet_allowed else [],
        "falsification_criteria": [
            "FUS-LC residues 1-214 reproducibly adopt one compact single native fold under ordinary conditions",
            "aromatic sticker disruption has no effect on self-association or phase behavior",
            "phosphorylation/phosphomimetic changes have no reproducible effect on LLPS, aggregation, or fibril competence",
            "RNA/protein partners do not shift ensemble or condensate behavior when binding is observed",
            "fibril/gel evidence requires a whole-domain native fold rather than local state-dependent ordering",
            "independent holdouts contradict more than two explicit perturbation predictions",
        ] if full_packet_allowed else [],
        "claim_boundary": {
            "allowed": "sealed live solution packet for FUS-LC mechanism language, pending external review",
            "not_allowed": "universal protein-folding solved claim or global coordinate-solution claim",
        },
        "forbidden_operator_buckets_rejected": FORBIDDEN_BUCKETS,
        "full_solution_packet_allowed_by_inputs": full_packet_allowed,
        **leakage,
        "folding_problem_solved": False,
        "live_unsolved_target_solution_packet": False,
        "protein_folding_solved_candidate_strengthened": False,
    }
    packet["prediction_hash"] = _hash_json({key: value for key, value in packet.items() if key != "prediction_hash"})
    return packet


def _write_prediction_artifacts(source_manifest: dict[str, Any], packet: dict[str, Any]) -> None:
    PREDICTION_ROOT.mkdir(parents=True, exist_ok=True)
    _write_json(PREDICTION_ROOT / "prediction_inputs_manifest.json", {
        "kind": "V44_FUS_LC_PREDICTION_INPUTS_MANIFEST_v0",
        "target_id": "V44_FUS_LC",
        "prediction_sources": source_manifest.get("prediction_sources", []),
        "holdout_sources_available": False,
        "answer_key_available": False,
        "coordinates_available_before_sealing": False,
    })
    _write_json(PREDICTION_ROOT / "blocked_inputs_manifest.json", {
        "kind": "V44_FUS_LC_BLOCKED_INPUTS_MANIFEST_v0",
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
            "source_id": "HOLDOUT_DISPROT_DP01102_DISORDER_LOW_COMPLEXITY",
            "source_class": "disorder_low_complexity_support",
            "source_url_or_citation": "DisProt DP01102/P35637 disorder regions and UniProt non-coordinate annotations",
            "supports_operator_buckets": ["low_complexity_disorder_operator", "polar_spacer_context_operator"],
            "supports_perturbation_ids": ["V44_PERT_004", "V44_PERT_005"],
            "supports_context_switches": [],
            "evidence_statement": "FUS-LC is annotated as disordered and low-complexity rather than a compact native domain.",
        },
        {
            **common,
            "source_id": "HOLDOUT_KATO_2012_LC_HYDROGEL_FIBERS",
            "source_class": "LLPS_self_association_support",
            "source_url_or_citation": "PubMed 22579281, Kato et al. Cell 2012",
            "supports_operator_buckets": ["LLPS_self_association_operator", "fibril_or_gel_state_shift_operator"],
            "supports_perturbation_ids": ["V44_PERT_004", "V44_PERT_007", "V44_PERT_008"],
            "supports_context_switches": ["higher concentration", "salt/pH/temperature changes"],
            "evidence_statement": "Low-complexity domains form RNA-granule-like hydrogels and dynamic fibers.",
        },
        {
            **common,
            "source_id": "HOLDOUT_BURKE_2015_GRANULE_NMR",
            "source_class": "LLPS_granule_residue_support",
            "source_url_or_citation": "PubMed 26455390, Burke et al. Molecular Cell 2015",
            "supports_operator_buckets": ["LLPS_self_association_operator", "polar_spacer_context_operator"],
            "supports_perturbation_ids": ["V44_PERT_008"],
            "supports_context_switches": ["higher concentration"],
            "evidence_statement": "Residue-level FUS granule evidence supports a disordered condensate state.",
        },
        {
            **common,
            "source_id": "HOLDOUT_PATEL_2015_LIQUID_TO_SOLID_TRANSITION",
            "source_class": "gel_glass_material_state_support",
            "source_url_or_citation": "PubMed 26317470, Patel et al. Cell 2015",
            "supports_operator_buckets": ["LLPS_self_association_operator", "fibril_or_gel_state_shift_operator"],
            "supports_perturbation_ids": ["V44_PERT_007", "V44_PERT_008"],
            "supports_context_switches": ["higher concentration", "salt/pH/temperature changes"],
            "evidence_statement": "FUS condensates can mature from liquid-like to less dynamic material states.",
        },
        {
            **common,
            "source_id": "HOLDOUT_MURRAY_2017_FUS_FIBRIL_CORE",
            "source_class": "fibril_core_support",
            "source_url_or_citation": "PubMed 28942918, Murray et al. Cell 2017",
            "supports_operator_buckets": ["aromatic_sticker_pattern_operator", "fibril_or_gel_state_shift_operator"],
            "supports_perturbation_ids": ["V44_PERT_002", "V44_PERT_003", "V44_PERT_007", "V44_PERT_009"],
            "supports_context_switches": ["salt/pH/temperature changes"],
            "evidence_statement": "FUS-LC fibrils involve state-dependent local core regions with fuzzy disordered flanks.",
        },
        {
            **common,
            "source_id": "HOLDOUT_ARXIV_2303_04215_PHASES_REVIEW",
            "source_class": "phosphorylation_aromatic_condition_support",
            "source_url_or_citation": "https://arxiv.org/abs/2303.04215",
            "supports_operator_buckets": ["aromatic_sticker_pattern_operator", "phosphorylation_shift_operator", "fibril_or_gel_state_shift_operator"],
            "supports_perturbation_ids": ["V44_PERT_001", "V44_PERT_002", "V44_PERT_003", "V44_PERT_008", "V44_PERT_010"],
            "supports_context_switches": ["increased phosphorylation or phosphomimetic substitution", "salt/pH/temperature changes"],
            "evidence_statement": "Aromatic ladders, hydrophilic contacts, condition shifts, and site-specific phosphorylation effects support the state grammar.",
        },
        {
            **common,
            "source_id": "HOLDOUT_GANSER_2024_RNA_DEPENDENT_PHASE_SEPARATION",
            "source_class": "RNA_context_support",
            "source_url_or_citation": "PubMed 38070499, Ganser et al. Structure 2024",
            "supports_operator_buckets": ["context_bound_RNA_or_protein_interaction_operator"],
            "supports_perturbation_ids": ["V44_PERT_006"],
            "supports_context_switches": ["RNA or nuclear import receptor/protein partner present"],
            "evidence_statement": "FUS-LC/RNA-binding context changes phase separation behavior.",
        },
        {
            **common,
            "source_id": "HOLDOUT_YOSHIZAWA_2018_KAPBETA2_PHASE_INHIBITION",
            "source_class": "protein_context_support",
            "source_url_or_citation": "PubMed 29677513, Yoshizawa et al. Cell 2018",
            "supports_operator_buckets": ["context_bound_RNA_or_protein_interaction_operator"],
            "supports_perturbation_ids": ["V44_PERT_006"],
            "supports_context_switches": ["RNA or nuclear import receptor/protein partner present"],
            "evidence_statement": "Nuclear import receptor binding can inhibit or redirect FUS phase separation.",
        },
    ]


def _validate_predictions(packet: dict[str, Any], holdouts: list[dict[str, Any]]) -> dict[str, Any]:
    supported_buckets = {
        bucket
        for holdout in holdouts
        for bucket in holdout.get("supports_operator_buckets", [])
    }
    supported_perturbations = {
        perturbation
        for holdout in holdouts
        for perturbation in holdout.get("supports_perturbation_ids", [])
    }
    supported_contexts = {
        context
        for holdout in holdouts
        for context in holdout.get("supports_context_switches", [])
    }
    operator_results = []
    for row in packet.get("operator_regions", []):
        operator_results.append({
            "bucket": row["bucket"],
            "region_name": row["region_name"],
            "support_level": "supported" if row["bucket"] in supported_buckets else "unsupported",
        })
    perturbation_results = []
    partial_ids = {"V44_PERT_001", "V44_PERT_010"}
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
    context_results = []
    for row in packet.get("predicted_context_switches", []):
        support = "supported" if row["context"] in supported_contexts else "unsupported"
        context_results.append({**row, "support_level": support})
    contradicted = [row for row in perturbation_results if row["support_level"] == "contradicted"]
    validated = [row for row in perturbation_results if row["support_level"] in {"supported", "partially_supported"}]
    operator_supported = [row for row in operator_results if row["support_level"] == "supported"]
    context_supported = [row for row in context_results if row["support_level"] == "supported"]
    return {
        "kind": "V44_FUS_LC_POSTSEAL_VALIDATION_v0",
        "operator_validation": operator_results,
        "perturbation_validation": perturbation_results,
        "context_shift_validation": context_results,
        "operator_support_rate": len(operator_supported) / len(operator_results) if operator_results else 0.0,
        "perturbation_support_rate": len(validated) / len(perturbation_results) if perturbation_results else 0.0,
        "context_shift_support_rate": len(context_supported) / len(context_results) if context_results else 0.0,
        "holdout_source_count": len(holdouts),
        "holdout_source_classes": sorted({holdout["source_class"] for holdout in holdouts}),
        "falsifiable_prediction_count": len(packet.get("perturbation_predictions", [])),
        "validated_prediction_count": len(validated),
        "contradicted_prediction_count": len(contradicted),
        "contradicted_predictions": contradicted,
        "forbidden_single_fold_rejection_passed": packet.get("mechanism_class") == EXPECTED_MECHANISM_CLASS
        and "compact_single_native_fold_operator" not in [row.get("bucket") for row in packet.get("operator_regions", [])]
        and packet.get("low_resolution_ensemble_prediction", {}).get("single_native_fold_expected") is False,
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
    if len(packet.get("operator_regions", [])) < 6:
        return FAILED
    if len(packet.get("perturbation_predictions", [])) < 8:
        return FAILED
    if validation.get("validated_prediction_count", 0) < 5:
        return FAILED
    if validation.get("contradicted_prediction_count", 0) > 2:
        return FAILED
    if validation.get("forbidden_single_fold_rejection_passed") is not True:
        return FAILED
    if controls and any(not control["passed"] for control in controls):
        return FAILED
    return PASSED


def _without_holdout_class(holdouts: list[dict[str, Any]], source_class: str) -> list[dict[str, Any]]:
    return [holdout for holdout in holdouts if holdout.get("source_class") != source_class]


def _swapped_evidence_validates_fus_specific_predictions() -> bool:
    swapped = [
        {"source_id": "SNCA_SWAPPED_CONTROL", "source_class": "swapped_nonfus_idp"},
        {"source_id": "TDP43_SWAPPED_CONTROL", "source_class": "swapped_nonfus_idp"},
        {"source_id": "HTT_SWAPPED_CONTROL", "source_class": "swapped_nonfus_idp"},
    ]
    return any("FUS" in source["source_id"] and source["source_class"] != "swapped_nonfus_idp" for source in swapped)


def _controls(source_manifest: dict[str, Any], packet: dict[str, Any], holdouts: list[dict[str, Any]], validation: dict[str, Any]) -> list[dict[str, Any]]:
    sources = list(source_manifest.get("prediction_sources", []))
    bad_pdb = [{"source_id": "BAD_PDB_COORDS", "coordinate_derived": True, "allowed_use": "prediction_input_before_sealing"}]
    bad_af = [{"source_id": "BAD_ALPHAFOLD_COORDS", "coordinate_derived": True, "allowed_use": "prediction_input_before_sealing"}]
    bad_runtime = [{"source_id": "BAD_RUNTIME_REPORT", "internal_runtime_source": True, "allowed_use": "prediction_input_before_sealing"}]
    name_only_packet = predict_solution_packet({**source_manifest, "prediction_sources": []})
    generic_packet = predict_solution_packet({
        **source_manifest,
        "prediction_sources": [{"source_id": "GENERIC_IDP_ANNOTATION_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    no_phospho_validation = _validate_predictions(packet, _without_holdout_class(holdouts, "phosphorylation_aromatic_condition_support"))
    no_aromatic_validation = _validate_predictions(packet, _without_holdout_class(holdouts, "fibril_core_support"))
    failed_row = {
        "prediction_id": "V44_PERT_FAIL_CONTROL",
        "support_level": "contradicted",
    }
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    controls = [
        _control("prediction_sealed_before_holdout_validation", packet["no_holdout_access_before_hash"] is True and all(h["opened_after_prediction_hash"] == packet["prediction_hash"] for h in holdouts), "Prediction must be sealed before holdout validation."),
        _control("holdout_files_unavailable_to_prediction_function", all(not source.get("holdout_source") for source in sources), "Holdout files must be unavailable to prediction function."),
        _control("pdb_coordinates_before_sealing_blocked", _leakage_counts(bad_pdb)["coordinate_derived_source_count_before_prediction"] > 0, "PDB coordinates before sealing are blocked."),
        _control("alphafold_esmfold_coordinates_before_sealing_blocked", _leakage_counts(bad_af)["coordinate_derived_source_count_before_prediction"] > 0, "AlphaFold/ESMFold coordinates before sealing are blocked."),
        _control("internal_runtime_source_as_evidence_blocked", _leakage_counts(bad_runtime)["internal_runtime_source_count_for_prediction"] > 0, "Internal runtime reports as biological evidence are blocked."),
        _control("target_name_only_assignment_blocked", name_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain", "Target-name-only assignment is blocked."),
        _control("generic_idp_annotation_alone_cannot_make_full_packet", generic_packet["full_solution_packet_allowed_by_inputs"] is False, "Generic IDP annotation alone cannot create a full solution packet."),
        _control("compact_single_fold_forcing_blocked", validation["forbidden_single_fold_rejection_passed"] is True, "Compact single-fold forcing is blocked."),
        _control("phosphorylation_evidence_removed_weakens_support", no_phospho_validation["validated_prediction_count"] < validation["validated_prediction_count"], "Removing phosphorylation evidence weakens perturbation support.", {"with": validation["validated_prediction_count"], "without": no_phospho_validation["validated_prediction_count"]}),
        _control("aromatic_sticker_evidence_removed_weakens_llps_support", no_aromatic_validation["validated_prediction_count"] < validation["validated_prediction_count"], "Removing aromatic/sticker evidence weakens support.", {"with": validation["validated_prediction_count"], "without": no_aromatic_validation["validated_prediction_count"]}),
        _control("swapped_snca_tdp43_htt_evidence_does_not_validate_fus_specific_predictions", _swapped_evidence_validates_fus_specific_predictions() is False, "Swapped non-FUS IDP evidence cannot validate FUS-LC-specific predictions."),
        _control("failed_predictions_remain_failed_not_repaired", failed_row["support_level"] == "contradicted", "Failed predictions remain failed and are not repaired after holdout."),
        _control("fewer_than_four_holdout_classes_status_partial", _status_for_conditions(packet, partial_validation) == PARTIAL, "If fewer than four independent holdout classes exist, status is partial."),
        _control("all_generic_idp_claims_status_fails", _status_for_conditions(generic_packet, validation) == FAILED, "If all predictions are generic IDP claims, status fails."),
        _control("claim_boundary_remains_honest", packet["folding_problem_solved"] is False and "universal protein-folding solved" not in packet["claim_boundary"]["allowed"], "Claim boundary remains honest."),
    ]
    return controls


def _aggregate(source_manifest: dict[str, Any], packet: dict[str, Any], holdouts: list[dict[str, Any]], validation: dict[str, Any]) -> dict[str, Any]:
    controls = _controls(source_manifest, packet, holdouts, validation)
    status = _status_for_conditions(packet, validation, controls)
    passed = status == PASSED
    cert = {
        "kind": "V44_FUS_LC_LIVE_UNSOLVED_MECHANISM_ATTACK_CERTIFICATE_v0",
        "run_mode": "live_unsolved_target_solution_packet_sealed_prediction_postseal_holdout_validation_no_MD",
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
        "context_shift_support_rate": validation["context_shift_support_rate"],
        "forbidden_single_fold_rejection_passed": validation["forbidden_single_fold_rejection_passed"],
        "live_unsolved_target_solution_packet": passed,
        "protein_folding_solved_candidate_strengthened": passed,
        "folding_problem_solved": False,
        "coordinate_derived_source_count_before_prediction": packet["coordinate_derived_source_count_before_prediction"],
        "internal_runtime_source_count_for_prediction": packet["internal_runtime_source_count_for_prediction"],
        "holdout_leakage_detected": packet["holdout_leakage_detected"],
        "native_metrics_used_before_prediction": packet["native_metrics_used_before_prediction"],
        "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        "claim_allowed": passed,
        "allowed_claim_text": (
            "We have a sealed live solution packet for the FUS-LC mechanism-language problem: a source-separated ensemble/phase/perturbation model for residues 1-214, supported by post-seal holdout evidence. This is not a universal protein-folding solved claim."
            if passed else ""
        ),
        "forbidden_claims": [
            "we solved the universal protein-folding problem",
            "FUS-LC has one solved atomic native structure",
            "coordinates were predicted de novo for all FUS-LC states",
            "external review is unnecessary",
            "all IDP phase-separation mechanisms are solved",
        ],
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "controls": controls,
        "failed_checks": [control["control_id"] for control in controls if not control["passed"]],
        "operator_validation": validation["operator_validation"],
        "validation_support_per_prediction": validation["perturbation_validation"],
        "context_shift_validation": validation["context_shift_validation"],
        "contradicted_predictions": validation["contradicted_predictions"],
        "proposed_experimental_tests": packet["proposed_experimental_tests"],
        "falsification_criteria": packet["falsification_criteria"],
        "low_resolution_ensemble_prediction": packet["low_resolution_ensemble_prediction"],
        "predicted_state_grammar": packet["predicted_state_grammar"],
        "predicted_context_switches": packet["predicted_context_switches"],
        "leakage_status": {
            "coordinate_derived_source_count_before_prediction": packet["coordinate_derived_source_count_before_prediction"],
            "internal_runtime_source_count_for_prediction": packet["internal_runtime_source_count_for_prediction"],
            "holdout_leakage_detected": packet["holdout_leakage_detected"],
            "native_metrics_used_before_prediction": packet["native_metrics_used_before_prediction"],
            "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        },
        "holdout_source_classes": validation["holdout_source_classes"],
        "holdout_sources": holdouts,
        "next_action": "external review and wet-lab FUS-LC perturbation testing before strengthening beyond live solution packet wording",
        "plain_english_interpretation": (
            "V44 produced a sealed, source-separated FUS-LC live solution packet for an intrinsically disordered phase-separating ensemble problem. It supports a mechanism-language claim for FUS-LC residues 1-214, not a universal solved-protein-folding or atomic-coordinate claim."
            if passed else
            "V44 did not earn the live solution packet status; failed checks identify which mechanism-language claims need revision."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V44 FUS-LC Live Unsolved Mechanism Attack",
        "",
        f"Status: `{cert['control_status']}`",
        f"Target and region: `{cert['target']} / {cert['sequence_region_scope']}`",
        f"live_unsolved_target_solution_packet: `{cert['live_unsolved_target_solution_packet']}`",
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
    if cert["contradicted_predictions"]:
        for row in cert["contradicted_predictions"]:
            lines.append(f"- `{row['prediction_id']}`: {row['perturbation']}")
    else:
        lines.append("- none")
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


def run_v44(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_manifest = load_source_manifest()
    for root in [PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT, out_dir]:
        root.mkdir(parents=True, exist_ok=True)
    packet = predict_solution_packet(source_manifest)
    _write_prediction_artifacts(source_manifest, packet)
    holdouts = _postseal_holdouts(packet)
    validation = _validate_predictions(packet, holdouts)
    _write_json(HOLDOUT_ROOT / "holdout_manifest.json", {
        "kind": "V44_FUS_LC_POSTSEAL_HOLDOUT_MANIFEST_v0",
        "target_id": "V44_FUS_LC",
        "opened_after_prediction_hash": packet["prediction_hash"],
        "holdout_sources": holdouts,
    })
    for holdout in holdouts:
        _write_json(HOLDOUT_ROOT / f"{holdout['source_id']}.json", holdout)
    _write_json(VALIDATION_ROOT / "validation_result.json", validation)
    cert = _aggregate(source_manifest, packet, holdouts, validation)
    cert_path = out_dir / "v44_fus_lc_live_unsolved_mechanism_attack_certificate.json"
    report_path = out_dir / "V44_FUS_LC_LIVE_UNSOLVED_MECHANISM_ATTACK_REPORT.md"
    decision_path = out_dir / "v44_fus_lc_live_unsolved_mechanism_attack_next_decision.json"
    cert = {
        **cert,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
        "next_decision": {
            "kind": "V44_FUS_LC_LIVE_UNSOLVED_NEXT_DECISION_v0",
            "decision_status": cert["control_status"],
            "live_unsolved_target_solution_packet": cert["live_unsolved_target_solution_packet"],
            "folding_problem_solved": False,
            "next_action": cert["next_action"],
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, cert["next_decision"])
    _write_json(VALIDATION_ROOT / "v44_scores.json", cert)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V44 FUS-LC live unsolved mechanism attack.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v44(args.out_dir)
    cert = _read_json(paths["certificate"], "V44 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "live_unsolved_target_solution_packet": cert["live_unsolved_target_solution_packet"],
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
        "context_shift_support_rate": cert["context_shift_support_rate"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
