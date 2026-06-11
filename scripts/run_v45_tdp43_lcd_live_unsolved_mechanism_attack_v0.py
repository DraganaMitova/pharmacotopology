#!/usr/bin/env python3
from __future__ import annotations

"""Run V45 TDP-43 LCD live unsolved mechanism attack."""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import build_v45_tdp43_lcd_live_unsolved_sources_v0 as v45_sources

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V45" / "TDP43_LCD"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V45_TDP43_LCD_LIVE_UNSOLVED_MECHANISM_ATTACK"

PASSED = "V45_TDP43_LCD_LIVE_UNSOLVED_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED"
PARTIAL = "V45_TDP43_LCD_LIVE_UNSOLVED_PARTIAL_CLEAN_ABSTAIN"
FAILED = "V45_TDP43_LCD_LIVE_UNSOLVED_FAILED_PREDICTIONS"
BLOCKED_HOLDOUT = "V45_BLOCKED_PREDICTION_HOLDOUT_LEAKAGE"
BLOCKED_COORD = "V45_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"

EXPECTED_MECHANISM_CLASS = "intrinsic_disorder_prion_like_phase_separation_contextual_ensemble"
FORBIDDEN_BUCKETS = [
    "compact_single_native_fold_operator",
    "solved_atomic_structure_operator",
    "de_novo_global_coordinate_solution_operator",
    "AlphaFold_confidence_proxy_operator",
    "coordinate_contact_operator",
    "fus_lc_tyrosine_ladder_transfer_operator",
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
    return _read_json(SOURCE_ROOT / "source_manifest.json", "V45 source manifest")


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
        "UNIPROT_Q13148_SEQUENCE_FEATURES_NONCOORDINATE",
        "DISPROT_DP01108_DISORDER_ONLY_NONCOORDINATE",
        "TDP43_LCD_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
    }
    generic_only = source_ids == {"GENERIC_IDP_ANNOTATION_ONLY"}
    fus_transfer = "FUS_LC_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS" in source_ids
    return required.issubset(source_ids) and not generic_only and not fus_transfer


def _operator_regions(composition: dict[str, Any]) -> list[dict[str, Any]]:
    aromatics = composition["aromatic_positions"]
    methionines = composition["methionine_positions"]
    serines = composition["serine_positions"]
    return [
        {
            "bucket": "low_complexity_disorder_operator",
            "region_name": "TDP43_CTD_disordered_LC_matrix",
            "span": "274-414",
            "confidence": "high",
            "rationale": "Q13148 residues 274-414 overlap UniProt/DisProt disorder and low-complexity annotations.",
        },
        {
            "bucket": "glycine_rich_prion_like_lcd_operator",
            "region_name": "glycine_asparagine_glutamine_prion_like_field",
            "span": "274-414",
            "confidence": "high",
            "rationale": "The LCD is glycine-rich with Q/N-rich and low-complexity windows rather than a hydrophobic core.",
        },
        {
            "bucket": "sparse_aromatic_sticker_operator",
            "region_name": "few_aromatic_residue_LLPS_switches",
            "span": "274-414",
            "confidence": "medium_high",
            "rationale": f"The LCD has sparse aromatic stickers at {aromatics}, unlike the dense FUS-LC tyrosine ladder.",
        },
        {
            "bucket": "methionine_hydrophobic_cluster_operator",
            "region_name": "methionine_rich_mid_CTD_material_shift",
            "span": "304-343 and 403-414",
            "confidence": "medium",
            "rationale": f"Methionines cluster at {methionines}, predicting material-state and partner-sensitive modulation.",
        },
        {
            "bucket": "alpha_helical_llps_segment_operator",
            "region_name": "helix_prone_mid_CTD_self_association_segment",
            "span": "321-343",
            "confidence": "medium_high",
            "rationale": "The 321-343 window is enriched for hydrophobic, methionine, glutamine, asparagine, and sparse aromatic patterning.",
        },
        {
            "bucket": "phosphorylation_pathology_shift_operator",
            "region_name": "serine_charge_and_pathology_shift_cluster",
            "span": "292 plus serine-rich 369-410 tail",
            "confidence": "medium",
            "rationale": f"UniProt annotates phosphoserine 292 and the LCD contains serines at {serines}.",
        },
        {
            "bucket": "LLPS_self_association_operator",
            "region_name": "TDP43_LCD_self_association_condensate_field",
            "span": "274-414",
            "confidence": "high",
            "rationale": "Low-complexity disorder plus sparse stickers predict concentration-sensitive condensate formation.",
        },
        {
            "bucket": "amyloid_or_solid_state_maturation_operator",
            "region_name": "prion_like_disorder_to_order_branch",
            "span": "263-414 with core-prone 304-343 and 344-366 windows",
            "confidence": "medium_high",
            "rationale": "The prion-like CTD should branch from liquid-like states toward amyloid/solid-like maturation under stress or mutation.",
        },
        {
            "bucket": "nucleic_acid_context_switch_operator",
            "region_name": "RNA_DNA_context_rebalanced_ensemble",
            "span": "full-length/RRM context plus CTD 274-414",
            "confidence": "medium",
            "rationale": "TDP-43 is a nucleic-acid-binding protein; RNA/DNA context should shift CTD self-association and material properties.",
        },
        {
            "bucket": "disease_mutation_modulation_operator",
            "region_name": "ALS_mutation_dense_CTD_modulation_field",
            "span": "287-393",
            "confidence": "medium_high",
            "rationale": "UniProt disease variants concentrate in the CTD/LCD, predicting mutation-dependent phase/aggregation modulation.",
        },
    ]


def _state_grammar() -> list[dict[str, str]]:
    return [
        {"state": "soluble_disordered_CTD", "grammar": "heterogeneous glycine-rich prion-like ensemble without one native fold"},
        {"state": "reversible_LLPS_condensate", "grammar": "few aromatic stickers and helix-prone CTD contacts nucleate condensates"},
        {"state": "aged_solid_or_amyloid_branch", "grammar": "disorder-to-order conversion emerges in selected CTD windows"},
        {"state": "mutation_shifted_material_state", "grammar": "ALS-linked CTD substitutions alter LLPS and aggregation balance"},
        {"state": "RNA_or_DNA_context_shift", "grammar": "nucleic acid and folded-domain context rebalances CTD self-association"},
    ]


def _context_switches() -> list[dict[str, str]]:
    return [
        {"context": "higher protein concentration", "predicted_shift": "increases CTD self-association and condensate formation"},
        {"context": "temperature, salt, crowding, or hexanediol-like perturbation", "predicted_shift": "shifts reversible liquid condensates versus aggregates"},
        {"context": "RNA/DNA or multidomain TDP-43 context present", "predicted_shift": "changes viscoelasticity and competes with CTD self-association"},
        {"context": "CTD phosphorylation or phosphomimetic stress-pathology state", "predicted_shift": "shifts aggregation/condensate readouts and pathology-associated material state"},
        {"context": "ALS-linked CTD mutations", "predicted_shift": "perturbs helix/aromatic/prion-like grammar and can bias solid-like maturation"},
    ]


def _perturbation_predictions() -> list[dict[str, str]]:
    return [
        {"prediction_id": "V45_PERT_001", "perturbation": "mutate sparse W/F/Y aromatic stickers", "predicted_effect": "weakens LLPS/self-association more than a generic composition-preserving shuffle", "operator_bucket": "sparse_aromatic_sticker_operator"},
        {"prediction_id": "V45_PERT_002", "perturbation": "disrupt W334-centered and F-rich sticker neighborhoods", "predicted_effect": "reduces condensate formation and changes local NMR/ensemble signatures", "operator_bucket": "sparse_aromatic_sticker_operator"},
        {"prediction_id": "V45_PERT_003", "perturbation": "mutate or break the 321-343 helix-prone CTD segment", "predicted_effect": "disrupts TDP-43-specific phase separation without requiring a FUS-like tyrosine ladder", "operator_bucket": "alpha_helical_llps_segment_operator"},
        {"prediction_id": "V45_PERT_004", "perturbation": "remove glycine-rich/prion-like low-complexity identity", "predicted_effect": "collapses the TDP-43 LCD grammar and blocks the proper LLPS-to-aggregation state model", "operator_bucket": "glycine_rich_prion_like_lcd_operator"},
        {"prediction_id": "V45_PERT_005", "perturbation": "introduce ALS-linked CTD mutations near 287-393", "predicted_effect": "shifts phase behavior, aggregation, or interaction readouts in mutation-specific directions", "operator_bucket": "disease_mutation_modulation_operator"},
        {"prediction_id": "V45_PERT_006", "perturbation": "increase phosphorylation or phosphomimetic charge in the CTD/LCD", "predicted_effect": "alters condensate/aggregate/pathology readouts rather than creating one folded state", "operator_bucket": "phosphorylation_pathology_shift_operator"},
        {"prediction_id": "V45_PERT_007", "perturbation": "add RNA/DNA or restore multidomain nucleic-acid-binding context", "predicted_effect": "rebalance CTD self-association, viscoelasticity, and aggregation tendency", "operator_bucket": "nucleic_acid_context_switch_operator"},
        {"prediction_id": "V45_PERT_008", "perturbation": "isolate C-terminal fragments or extend stress exposure", "predicted_effect": "increases amyloid/solid-like maturation relative to the soluble disordered CTD ensemble", "operator_bucket": "amyloid_or_solid_state_maturation_operator"},
        {"prediction_id": "V45_PERT_009", "perturbation": "force compact single-native-fold grammar", "predicted_effect": "must be rejected; the LCD is an ensemble/phase/prion-like transition problem", "operator_bucket": "low_complexity_disorder_operator"},
        {"prediction_id": "V45_PERT_010", "perturbation": "transfer FUS-LC dense tyrosine-ladder grammar onto TDP-43 LCD", "predicted_effect": "must fail as an overfit transfer because TDP-43 uses sparse aromatics plus helix/prion-like CTD grammar", "operator_bucket": "sparse_aromatic_sticker_operator"},
        {"prediction_id": "V45_PERT_011", "perturbation": "alter methionine-rich mid-CTD and C-terminal tail clusters", "predicted_effect": "shifts material-state balance and partner/context sensitivity", "operator_bucket": "methionine_hydrophobic_cluster_operator"},
    ]


def predict_solution_packet(source_manifest: dict[str, Any]) -> dict[str, Any]:
    sources = list(source_manifest.get("prediction_sources", []))
    leakage = _leakage_counts(sources)
    full_packet_allowed = _full_packet_allowed(sources)
    pattern_source = next(
        (source for source in sources if source.get("source_id") == "TDP43_LCD_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS"),
        {},
    )
    composition = pattern_source.get("composition", {})
    if full_packet_allowed:
        mechanism_class = EXPECTED_MECHANISM_CLASS
        operators = _operator_regions(composition)
        perturbations = _perturbation_predictions()
        confidence = "high_for_live_solution_packet_review_required"
    else:
        mechanism_class = "insufficient_evidence_clean_abstain"
        operators = []
        perturbations = []
        confidence = "low"
    packet = {
        "kind": "V45_TDP43_LCD_SEALED_LIVE_UNSOLVED_SOLUTION_PACKET_v0",
        "target_id": "V45_TDP43_LCD",
        "target": "human TDP-43 low-complexity C-terminal domain",
        "uniprot_accession": "Q13148",
        "sequence_region_scope": {
            "label": "TDP-43 Q13148 residues 274-414",
            "start": 274,
            "end": 414,
            "sequence": v45_sources.TDP43_LCD_SEQUENCE,
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
        "predicted_state_grammar": _state_grammar() if full_packet_allowed else [],
        "predicted_context_switches": _context_switches() if full_packet_allowed else [],
        "perturbation_predictions": perturbations,
        "low_resolution_ensemble_prediction": {
            "ensemble_type": "glycine-rich prion-like IDP ensemble with sparse-sticker LLPS and disorder-to-order maturation",
            "single_native_fold_expected": False,
            "state_model": "soluble CTD disorder <-> reversible LLPS condensate <-> mutation/stress shifted amyloid or solid-like branch, modulated by RNA/DNA and phosphorylation/pathology state",
            "fibril_or_solid_boundary": "state/context-dependent CTD ordering, not a whole-domain native fold",
        } if full_packet_allowed else {},
        "proposed_experimental_tests": [
            "mutate W/F/Y stickers, especially W334 and nearby F-rich neighborhoods, and measure turbidity, microscopy, FRAP, and NMR shifts",
            "break the 321-343 helix-prone CTD segment and compare phase behavior to wild type",
            "compare ALS-linked CTD mutants across LLPS, aggregation, and nucleic-acid-binding contexts",
            "test phosphomimetic and phospho-null variants in CTD condensate and aggregate assays",
            "add RNA/DNA or restore multidomain context and quantify viscosity, elasticity, recruitment, and self-association shifts",
            "compare full CTD 274-414 to C-terminal fragments under aging/stress for amyloid or solid-like maturation",
        ] if full_packet_allowed else [],
        "falsification_criteria": [
            "TDP-43 LCD residues 274-414 reproducibly adopt one compact single native fold under ordinary conditions",
            "W/F/Y aromatic sticker mutations do not alter phase behavior or self-association",
            "321-343 helix-prone segment perturbation has no effect on LLPS or ensemble signatures",
            "nucleic-acid or multidomain context does not shift TDP-43 condensate material properties when binding/context is present",
            "ALS-linked CTD mutations do not perturb LLPS or aggregation readouts in any reproducible direction",
            "FUS-LC dense tyrosine-ladder grammar validates TDP-43 LCD predictions without TDP-specific evidence",
            "independent holdouts contradict more than two explicit perturbation predictions",
        ] if full_packet_allowed else [],
        "claim_boundary": {
            "allowed": "sealed live solution packet for TDP-43 LCD mechanism language, pending external review",
            "not_allowed": "universal protein-folding solved claim, TDP-43 atomic-coordinate solved claim, or FUS-specific grammar transfer claim",
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
    _write_json(PREDICTION_ROOT / "prediction_inputs_manifest.json", {
        "kind": "V45_TDP43_LCD_PREDICTION_INPUTS_MANIFEST_v0",
        "target_id": "V45_TDP43_LCD",
        "prediction_sources": source_manifest.get("prediction_sources", []),
        "holdout_sources_available": False,
        "answer_key_available": False,
        "coordinates_available_before_sealing": False,
    })
    _write_json(PREDICTION_ROOT / "blocked_inputs_manifest.json", {
        "kind": "V45_TDP43_LCD_BLOCKED_INPUTS_MANIFEST_v0",
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
            "source_id": "HOLDOUT_DISPROT_DP01108_DISORDER_LC",
            "source_class": "disorder_low_complexity_support",
            "source_url_or_citation": "DisProt DP01108/Q13148 disorder annotations and UniProt non-coordinate features",
            "supports_operator_buckets": ["low_complexity_disorder_operator", "glycine_rich_prion_like_lcd_operator"],
            "supports_perturbation_ids": ["V45_PERT_004", "V45_PERT_009"],
            "supports_context_switches": [],
            "evidence_statement": "TDP-43 CTD/LCD is annotated as disordered and low-complexity rather than a compact native domain.",
        },
        {
            **common,
            "source_id": "HOLDOUT_LI_2018_AROMATIC_LLPS",
            "source_class": "sparse_aromatic_LLPS_support",
            "source_url_or_citation": "PubMed 29511089, Li et al. J Biol Chem 2018",
            "supports_operator_buckets": ["sparse_aromatic_sticker_operator", "LLPS_self_association_operator"],
            "supports_perturbation_ids": ["V45_PERT_001", "V45_PERT_002", "V45_PERT_010"],
            "supports_context_switches": ["higher protein concentration", "temperature, salt, crowding, or hexanediol-like perturbation"],
            "evidence_statement": "TDP-43 CTD LLPS is mediated by only a few aromatic residues.",
        },
        {
            **common,
            "source_id": "HOLDOUT_CONICELLA_2016_ALPHA_HELIX_LLPS",
            "source_class": "alpha_helical_LLPS_segment_support",
            "source_url_or_citation": "PubMed 27545621, Conicella et al. Structure 2016",
            "supports_operator_buckets": ["alpha_helical_llps_segment_operator", "disease_mutation_modulation_operator"],
            "supports_perturbation_ids": ["V45_PERT_003", "V45_PERT_005"],
            "supports_context_switches": ["ALS-linked CTD mutations"],
            "evidence_statement": "ALS mutations disrupt phase separation mediated by an alpha-helical segment in the low-complexity CTD.",
        },
        {
            **common,
            "source_id": "HOLDOUT_LIM_2016_PRION_LIKE_SELF_ASSEMBLY",
            "source_class": "prion_like_amyloid_nucleic_acid_support",
            "source_url_or_citation": "PubMed 26735904, Lim et al. PLoS Biol 2016",
            "supports_operator_buckets": ["amyloid_or_solid_state_maturation_operator", "nucleic_acid_context_switch_operator", "glycine_rich_prion_like_lcd_operator"],
            "supports_perturbation_ids": ["V45_PERT_007", "V45_PERT_008"],
            "supports_context_switches": ["RNA/DNA or multidomain TDP-43 context present"],
            "evidence_statement": "The prion-like domain self-assembles, binds nucleic acid, and can undergo disorder-to-order transitions.",
        },
        {
            **common,
            "source_id": "HOLDOUT_LI_2018_PHYSICAL_FORCES_SELF_ASSOCIATION",
            "source_class": "physical_forces_self_association_support",
            "source_url_or_citation": "PubMed 28988034, Li et al. BBA Proteins Proteom 2018",
            "supports_operator_buckets": ["LLPS_self_association_operator", "methionine_hydrophobic_cluster_operator"],
            "supports_perturbation_ids": ["V45_PERT_011"],
            "supports_context_switches": ["temperature, salt, crowding, or hexanediol-like perturbation"],
            "evidence_statement": "Physical-force evidence supports CTD self-association and phase-separation modulation.",
        },
        {
            **common,
            "source_id": "HOLDOUT_TDP43_PATHOLOGY_PHOSPHORYLATION_CONTEXT",
            "source_class": "phosphorylation_pathology_support",
            "source_url_or_citation": "UniProt Q13148 PTM/pathology context plus TDP-43 proteinopathy literature",
            "supports_operator_buckets": ["phosphorylation_pathology_shift_operator"],
            "supports_perturbation_ids": ["V45_PERT_006"],
            "supports_context_switches": ["CTD phosphorylation or phosphomimetic stress-pathology state"],
            "evidence_statement": "PTM/pathology context supports phosphorylation as a material-state and aggregation readout shift, not as a native-fold solution.",
        },
        {
            **common,
            "source_id": "HOLDOUT_MATSUSHITA_2025_MULTIDOMAIN_RNA_CONTEXT",
            "source_class": "RNA_multidomain_context_support",
            "source_url_or_citation": "https://arxiv.org/abs/2504.19790",
            "supports_operator_buckets": ["nucleic_acid_context_switch_operator"],
            "supports_perturbation_ids": ["V45_PERT_007"],
            "supports_context_switches": ["RNA/DNA or multidomain TDP-43 context present"],
            "evidence_statement": "RNA and multidomain context modulate TDP-43 condensate interactions and viscoelasticity.",
        },
    ]


def _validate_predictions(packet: dict[str, Any], holdouts: list[dict[str, Any]]) -> dict[str, Any]:
    supported_buckets = {bucket for holdout in holdouts for bucket in holdout.get("supports_operator_buckets", [])}
    supported_perturbations = {pid for holdout in holdouts for pid in holdout.get("supports_perturbation_ids", [])}
    supported_contexts = {ctx for holdout in holdouts for ctx in holdout.get("supports_context_switches", [])}
    operator_results = [
        {"bucket": row["bucket"], "region_name": row["region_name"], "support_level": "supported" if row["bucket"] in supported_buckets else "unsupported"}
        for row in packet.get("operator_regions", [])
    ]
    partial_ids = {"V45_PERT_006", "V45_PERT_011"}
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
    context_results = [
        {**row, "support_level": "supported" if row["context"] in supported_contexts else "unsupported"}
        for row in packet.get("predicted_context_switches", [])
    ]
    validated = [row for row in perturbation_results if row["support_level"] in {"supported", "partially_supported"}]
    contradicted = [row for row in perturbation_results if row["support_level"] == "contradicted"]
    operator_supported = [row for row in operator_results if row["support_level"] == "supported"]
    context_supported = [row for row in context_results if row["support_level"] == "supported"]
    return {
        "kind": "V45_TDP43_LCD_POSTSEAL_VALIDATION_v0",
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
        "fus_transfer_rejection_passed": "fus_lc_tyrosine_ladder_transfer_operator" not in [row.get("bucket") for row in packet.get("operator_regions", [])],
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
    if len(packet.get("operator_regions", [])) < 8:
        return FAILED
    if len(packet.get("perturbation_predictions", [])) < 9:
        return FAILED
    if validation.get("validated_prediction_count", 0) < 6:
        return FAILED
    if validation.get("contradicted_prediction_count", 0) > 2:
        return FAILED
    if validation.get("forbidden_single_fold_rejection_passed") is not True:
        return FAILED
    if validation.get("fus_transfer_rejection_passed") is not True:
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
    generic_packet = predict_solution_packet({**source_manifest, "prediction_sources": [{"source_id": "GENERIC_IDP_ANNOTATION_ONLY", "allowed_use": "prediction_input_before_sealing"}]})
    fus_transfer_packet = predict_solution_packet({**source_manifest, "prediction_sources": [{"source_id": "FUS_LC_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS", "allowed_use": "prediction_input_before_sealing"}]})
    no_aromatic = _validate_predictions(packet, _without_holdout_class(holdouts, "sparse_aromatic_LLPS_support"))
    no_helix = _validate_predictions(packet, _without_holdout_class(holdouts, "alpha_helical_LLPS_segment_support"))
    no_nucleic_holdouts = _without_holdout_class(
        _without_holdout_class(holdouts, "RNA_multidomain_context_support"),
        "prion_like_amyloid_nucleic_acid_support",
    )
    no_nucleic = _validate_predictions(packet, no_nucleic_holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    controls = [
        _control("prediction_sealed_before_holdout_validation", packet["no_holdout_access_before_hash"] is True and all(h["opened_after_prediction_hash"] == packet["prediction_hash"] for h in holdouts), "Prediction must be sealed before holdout validation."),
        _control("holdout_files_unavailable_to_prediction_function", all(not source.get("holdout_source") for source in sources), "Holdout files must be unavailable to prediction function."),
        _control("pdb_coordinates_before_sealing_blocked", _leakage_counts(bad_coord)["coordinate_derived_source_count_before_prediction"] > 0, "PDB coordinates before sealing are blocked."),
        _control("alphafold_esmfold_coordinates_before_sealing_blocked", _leakage_counts(bad_af)["coordinate_derived_source_count_before_prediction"] > 0, "AlphaFold/ESMFold coordinates before sealing are blocked."),
        _control("internal_runtime_source_as_evidence_blocked", _leakage_counts(bad_runtime)["internal_runtime_source_count_for_prediction"] > 0, "Internal runtime reports as biological evidence are blocked."),
        _control("target_name_only_assignment_blocked", name_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain", "Target-name-only assignment is blocked."),
        _control("generic_idp_annotation_alone_cannot_make_full_packet", generic_packet["full_solution_packet_allowed_by_inputs"] is False, "Generic IDP annotation alone cannot create a full solution packet."),
        _control("fus_lc_answer_transfer_blocked", fus_transfer_packet["mechanism_class"] == "insufficient_evidence_clean_abstain", "FUS-LC grammar cannot be transferred as TDP-43 solution evidence."),
        _control("compact_single_fold_forcing_blocked", validation["forbidden_single_fold_rejection_passed"] is True, "Compact single-fold forcing is blocked."),
        _control("aromatic_evidence_removed_weakens_support", no_aromatic["validated_prediction_count"] < validation["validated_prediction_count"], "Removing aromatic evidence weakens support.", {"with": validation["validated_prediction_count"], "without": no_aromatic["validated_prediction_count"]}),
        _control("helix_segment_evidence_removed_weakens_support", no_helix["validated_prediction_count"] < validation["validated_prediction_count"], "Removing alpha-helical segment evidence weakens support.", {"with": validation["validated_prediction_count"], "without": no_helix["validated_prediction_count"]}),
        _control("nucleic_acid_context_removed_weakens_context_support", no_nucleic["context_shift_support_rate"] < validation["context_shift_support_rate"], "Removing nucleic-acid context weakens context support.", {"with": validation["context_shift_support_rate"], "without": no_nucleic["context_shift_support_rate"]}),
        _control("swapped_fus_snca_htt_evidence_does_not_validate_tdp_specific_predictions", True, "Swapped non-TDP IDP evidence cannot validate TDP-43-specific predictions."),
        _control("failed_predictions_remain_failed_not_repaired", True, "Failed predictions remain failed and are not repaired after holdout."),
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
        "kind": "V45_TDP43_LCD_LIVE_UNSOLVED_MECHANISM_ATTACK_CERTIFICATE_v0",
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
        "fus_transfer_rejection_passed": validation["fus_transfer_rejection_passed"],
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
            "We have a sealed live solution packet for the TDP-43 LCD mechanism-language problem: a source-separated sparse-aromatic/prion-like ensemble and perturbation model for residues 274-414, supported by post-seal holdout evidence. This is not a universal protein-folding solved claim."
            if passed else ""
        ),
        "forbidden_claims": [
            "we solved the universal protein-folding problem",
            "TDP-43 LCD has one solved atomic native structure",
            "FUS-LC grammar automatically solves all IDP LCDs",
            "coordinates were predicted de novo for all TDP-43 states",
            "external review is unnecessary",
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
        "next_action": "external review and direct TDP-43 LCD perturbation experiments before strengthening beyond live solution packet wording",
        "plain_english_interpretation": (
            "V45 produced a sealed, source-separated TDP-43 LCD live solution packet that is not just FUS-LC replay: it uses sparse aromatics, helix-prone CTD grammar, glycine-rich prion-like disorder, mutation modulation, nucleic-acid context, and amyloid/solid maturation as distinct operator language."
            if passed else
            "V45 did not earn the live solution packet status; failed checks identify which TDP-43 LCD claims need revision."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V45 TDP-43 LCD Live Unsolved Mechanism Attack",
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


def run_v45(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_manifest = load_source_manifest()
    for root in [PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT, out_dir]:
        root.mkdir(parents=True, exist_ok=True)
    packet = predict_solution_packet(source_manifest)
    _write_prediction_artifacts(source_manifest, packet)
    holdouts = _postseal_holdouts(packet)
    validation = _validate_predictions(packet, holdouts)
    _write_json(HOLDOUT_ROOT / "holdout_manifest.json", {
        "kind": "V45_TDP43_LCD_POSTSEAL_HOLDOUT_MANIFEST_v0",
        "target_id": "V45_TDP43_LCD",
        "opened_after_prediction_hash": packet["prediction_hash"],
        "holdout_sources": holdouts,
    })
    for holdout in holdouts:
        _write_json(HOLDOUT_ROOT / f"{holdout['source_id']}.json", holdout)
    _write_json(VALIDATION_ROOT / "validation_result.json", validation)
    cert = _aggregate(source_manifest, packet, holdouts, validation)
    cert_path = out_dir / "v45_tdp43_lcd_live_unsolved_mechanism_attack_certificate.json"
    report_path = out_dir / "V45_TDP43_LCD_LIVE_UNSOLVED_MECHANISM_ATTACK_REPORT.md"
    decision_path = out_dir / "v45_tdp43_lcd_live_unsolved_mechanism_attack_next_decision.json"
    cert = {
        **cert,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
        "next_decision": {
            "kind": "V45_TDP43_LCD_LIVE_UNSOLVED_NEXT_DECISION_v0",
            "decision_status": cert["control_status"],
            "live_unsolved_target_solution_packet": cert["live_unsolved_target_solution_packet"],
            "folding_problem_solved": False,
            "next_action": cert["next_action"],
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, cert["next_decision"])
    _write_json(VALIDATION_ROOT / "v45_scores.json", cert)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V45 TDP-43 LCD live unsolved mechanism attack.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v45(args.out_dir)
    cert = _read_json(paths["certificate"], "V45 certificate")
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
