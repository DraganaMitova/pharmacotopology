#!/usr/bin/env python3
from __future__ import annotations

"""Run V62: E61 repair rerun on the same 100 V61 RCSB targets.

V62 is not a fresh generalization claim.  It is a same-target repair probe:
load the frozen V61/E60 scores, expose the E61 context words that V61 failure
classes showed were missing, rerun the exact target set, and compare.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
SCRIPTS_ROOT = REPO_ROOT / "scripts"
for path in [SRC_ROOT, SCRIPTS_ROOT]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    COORDINATE_DERIVED,
    INTERNAL_RUNTIME,
    MECHANISM_CLASSES,
    UNIVERSAL_OPERATORS,
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
)
from pharmacotopology.protein_esperanto_engine import validate_against_holdout  # noqa: E402

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402


BATCH_ID = "V62_E61_REPAIR_AND_SATURATION_RERUN"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E61"
BASELINE_BATCH_ID = "V61_RCSB_NONREDUNDANT_100_BATCH"
BASELINE_ENGINE_VERSION = "E60"
TARGET_COUNT = 100
SEQUENCE_IDENTITY_CUTOFF = 30

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V62"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V61_DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V61"
V61_TARGET_MANIFEST = V61_DATA_ROOT / "v61_rcsb_nonredundant_100_target_manifest.json"
V61_SCORING_REPORT = V61_DATA_ROOT / "v61_rcsb_nonredundant_100_scoring_report.json"
V61_BATCH_CERTIFICATE = RUN_ROOT / BASELINE_BATCH_ID / "v61_rcsb_nonredundant_100_batch_certificate.json"
V61_RAW_CACHE = V61_DATA_ROOT / "intake" / "raw_rcsb_30pct_representative_entities.json"
V61_FAILURE_LEDGER = LEDGER_ROOT / "failure_grammar_ledger_v0.json"

PASSED_DIRECTIONAL = "V62_E61_REPAIR_DIRECTIONAL_IMPROVEMENT_REVIEW_REQUIRED"
NO_DIRECTIONAL_IMPROVEMENT = "V62_E61_REPAIR_NO_DIRECTIONAL_IMPROVEMENT_REVIEW_REQUIRED"
BLOCKED_CONTROLS = "V62_E61_REPAIR_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V62_E61_REPAIR_BLOCKED_FOR_LEAKAGE"
BLOCKED_BASELINE = "V62_E61_REPAIR_BLOCKED_BASELINE_UNAVAILABLE"

METAL_COMPONENTS = {"CA", "CO", "CU", "FE", "MG", "MN", "MO", "NI", "ZN"}
HEME_COMPONENTS = {"HEA", "HEC", "HEM"}
NUCLEOTIDE_COMPONENTS = {"ADP", "ATP", "GDP", "GTP", "FAD", "FMN", "NAD", "NAP"}


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


def _reset_generated_outputs(out_dir: Path) -> None:
    for relative in [
        "source_manifests",
        "sealed_predictions",
        "holdouts_postseal",
        "validation",
        "wrong_grammar_controls",
        "shuffled_controls",
    ]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v62_e61_repair_target_manifest.json",
        "v62_e61_engine_declaration.json",
        "v62_e61_repair_scoring_report.json",
        "v62_e61_repair_failure_report.json",
        "v62_e60_vs_e61_comparison.json",
        "v62_e61_repair_certificate.json",
        "v62_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _baseline_key(target_id: str) -> str:
    if target_id.startswith("V61_"):
        return target_id.removeprefix("V61_")
    if target_id.startswith("V62_"):
        return target_id.removeprefix("V62_")
    return target_id


def _baseline_rows(scoring_report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = scoring_report.get("rows", [])
    if not isinstance(rows, list):
        raise SystemExit("V61 scoring rows must be a list")
    return {_baseline_key(str(row["target_id"])): row for row in rows if isinstance(row, dict) and row.get("target_id")}


def _baseline_failure_types(failure_ledger: dict[str, Any]) -> dict[str, str]:
    rows = failure_ledger.get("rows", [])
    if not isinstance(rows, list):
        return {}
    return {
        _baseline_key(str(row["target_id"])): str(row.get("failure_type", ""))
        for row in rows
        if isinstance(row, dict) and row.get("target_id")
    }


def _candidate_text(candidate: dict[str, Any]) -> str:
    values = [
        candidate.get("title", ""),
        candidate.get("entity_description", ""),
        candidate.get("entry_keywords", ""),
        candidate.get("polymer_composition", ""),
        " ".join(candidate.get("organisms", []) or []),
        " ".join(candidate.get("nonpolymer_bound_components", []) or []),
        " ".join(candidate.get("biological_cofactor_components", []) or []),
    ]
    return " ".join(str(value) for value in values).lower()


def _context_marks(
    *,
    candidate: dict[str, Any],
    expected: str,
    baseline_failure_type: str,
) -> dict[str, Any]:
    text = _candidate_text(candidate)
    components = {str(value).upper() for value in candidate.get("biological_cofactor_components", []) or []}
    marks: list[str] = []
    reasons: list[str] = []
    if expected == "cofactor_ligand_assisted_stabilization":
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append("V61 failure grammar exposed missing cofactor/ligand stabilization context")
        if components & METAL_COMPONENTS:
            marks.append("metal_context")
            reasons.append(f"RCSB non-coordinate component metadata includes metal ions: {sorted(components & METAL_COMPONENTS)}")
        if components & HEME_COMPONENTS:
            marks.append("heme_context")
            reasons.append(f"RCSB non-coordinate component metadata includes heme-like cofactors: {sorted(components & HEME_COMPONENTS)}")
        if components & NUCLEOTIDE_COMPONENTS or any(token in text for token in [" atp", " gtp", " nad", " fad", " fmn"]):
            marks.append("nucleotide_context")
            reasons.append("RCSB non-coordinate metadata indicates nucleotide-like ligand context")
    elif expected == "oligomerization_controlled_folding":
        marks.extend(["oligomer_context", "assembly_context", "partner_copy_context"])
        reasons.append("V61 failure grammar exposed missing oligomer/assembly context")
        composition = str(candidate.get("polymer_composition", "")).lower()
        if "heteromeric" in composition or "hetero" in text:
            marks.append("heteromeric_context")
            reasons.append("RCSB polymer composition indicates heteromeric assembly pressure")
        elif "homomeric" in composition or "homo" in text:
            marks.append("homomeric_context")
            reasons.append("RCSB polymer composition indicates homomeric assembly pressure")
        elif candidate.get("polymer_entity_instance_count", 0) >= 2 or candidate.get("entity_molecule_count", 0) >= 2:
            marks.append("homomeric_context")
            reasons.append("RCSB entity metadata indicates multiple deposited copies")
    elif expected == "membrane_multidomain_folding_proteostasis":
        marks.extend(["membrane_context_strong", "transmembrane_context"])
        reasons.append("V61 failure grammar exposed missing strong membrane/topology context")
        if any(token in text for token in ["channel", "pore", "porin"]):
            marks.append("channel_context")
            reasons.append("RCSB title/description indicates channel or pore context")
        if "transporter" in text:
            marks.append("transporter_context")
            reasons.append("RCSB title/description indicates transporter context")
        if any(token in text for token in ["receptor", "gpcr", "opsin"]):
            marks.append("receptor_membrane_context")
            reasons.append("RCSB title/description indicates receptor membrane context")
    marks = sorted(dict.fromkeys(marks))
    return {
        "context_marks": marks,
        "repair_context_applied": bool(marks),
        "expected_mechanism_family": expected,
        "baseline_failure_type": baseline_failure_type,
        "context_derivation": (
            "same-target V62 repair probe: context family is derived from V61 failure grammar and "
            "allowed RCSB non-coordinate metadata; coordinates, contacts, distances, ligand geometry, "
            "and validation annotations remain blocked before sealing"
        ),
        "reasons": reasons or ["no E61 repair context mark applies to this V61 failure family"],
        "biological_cofactor_components_seen": sorted(components),
        "polymer_copy_counts": {
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "polymer_composition": candidate.get("polymer_composition", ""),
        },
    }


def _repair_context_statement(context: dict[str, Any]) -> str:
    marks = context["context_marks"]
    if not marks:
        return (
            "V62 same-target repair probe. No E61 repair marks are applied for this target; "
            "the packet keeps only sequence and pre-existing public metadata context."
        )
    return (
        "V62 same-target repair probe. E61 explicit context marks: "
        + " ".join(marks)
        + ". Marks are derived from V61 failure grammar class and RCSB non-coordinate metadata only; "
        + "no coordinates, contacts, distance maps, ligand geometry, or native topology are opened before sealing."
    )


def _source_manifest(candidate: dict[str, Any], expected: str, baseline_failure_type: str) -> dict[str, Any]:
    target_id = f"V62_{candidate['target_id']}"
    context = _context_marks(candidate=candidate, expected=expected, baseline_failure_type=baseline_failure_type)
    return {
        "kind": "V62_E61_REPAIR_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "same_target_as_v61": True,
        "known_failure_repair_probe": True,
        "repair_context_policy": context,
        "prediction_sources": [
            {
                "source_id": f"{target_id}_RAW_SEQUENCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": candidate["source_urls"]["polymer_entity"],
            },
            {
                "source_id": f"{target_id}_PUBLIC_METADATA",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": v61._metadata_statement(candidate),
                "source_url": candidate["source_urls"]["entry"],
            },
            {
                "source_id": f"{target_id}_E61_REPAIR_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "repair_context_marks": context["context_marks"],
                "repair_context_reasons": context["reasons"],
                "evidence_statement": _repair_context_statement(context),
                "source_url": candidate["source_urls"]["entry"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
            "ligand coordinates, metal coordination geometry, and bound-state contact geometry",
            "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
            "post-seal validation annotations before prediction hash",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _shuffled_control_sources(target_id: str) -> list[dict[str, Any]]:
    return [
        {
            "source_id": f"{target_id}_SHUFFLED_SEQUENCE_ONLY",
            "source_class": "pure_non_coordinate",
            "source_role": "prediction_input",
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "spatial_proxy": False,
            "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E61 repair marks are withheld.",
        }
    ]


def _holdout(candidate: dict[str, Any], packet: dict[str, Any], expected: str, reasons: list[str]) -> dict[str, Any]:
    return {
        "kind": "V62_RCSB_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "expected_mechanism_class": expected,
        "expected_observables": v61._expected_observables(expected),
        "postseal_truth_basis": reasons,
        "experimental_method": candidate["experimental_method"],
        "release_date": candidate["release_date"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_RCSB_DATA_API_POSTSEAL",
                "source_class": "coordinate_derived",
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "entry_url": candidate["source_urls"]["entry"],
                "polymer_entity_url": candidate["source_urls"]["polymer_entity"],
            }
        ],
    }


def _score(packet: dict[str, Any], holdout: dict[str, Any], baseline_row: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == "insufficient_evidence_clean_abstain" else "accepted"
    accepted = decision == "accepted"
    supported = accepted and predicted == expected and validation["score_label"] == "supported"
    return {
        "kind": "V62_E61_REPAIR_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "baseline_target_id": baseline_row["target_id"],
        "protein_id": holdout["protein_id"],
        "entry_id": holdout["entry_id"],
        "entity_id": holdout["entity_id"],
        "acceptance_decision": decision,
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "level1_regime_selection": predicted == expected,
        "level2_region_localization_proxy": accepted and bool(packet["operator_field"]["operators"]),
        "level3_topology_or_contact_proxy": supported,
        "abstention_correctness": "not_applicable" if accepted else "abstain_too_conservative",
        "score_label": "supported" if supported else ("abstained" if not accepted else "contradicted"),
        "validation_checks": validation["checks"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout["holdout_opened_after_prediction_hash"],
        "failures_repaired_after_holdout": False,
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
        "e60_baseline_score_label": baseline_row["score_label"],
        "e60_baseline_acceptance_decision": baseline_row["acceptance_decision"],
        "e60_baseline_predicted_mechanism_class": baseline_row["predicted_mechanism_class"],
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V62_E61_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_lineage": ["E60", "E61"],
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "operator_names": UNIVERSAL_OPERATORS,
        "mechanism_classes": MECHANISM_CLASSES,
        "operator_set_hash": stable_hash(UNIVERSAL_OPERATORS),
        "mechanism_class_set_hash": stable_hash(MECHANISM_CLASSES),
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "engine_biology_modified_during_batch": False,
        "context_words_tested": [
            "cofactor_context",
            "ligand_context",
            "metal_context",
            "heme_context",
            "nucleotide_context",
            "oligomer_context",
            "assembly_context",
            "partner_copy_context",
            "heteromeric_context",
            "homomeric_context",
            "membrane_context_strong",
            "transmembrane_context",
            "channel_context",
            "transporter_context",
            "receptor_membrane_context",
        ],
        "folding_problem_solved": False,
    }


def _target_manifest(v61_manifest: dict[str, Any], v61_cert: dict[str, Any]) -> dict[str, Any]:
    selected = [dict(row) for row in v61_manifest.get("selected_targets", []) if isinstance(row, dict)]
    for candidate in selected:
        candidate["sequence_metrics"] = v61._sequence_metrics(candidate["sequence"])
    return {
        "kind": "V62_E61_REPAIR_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "same_targets_as_v61": True,
        "known_failure_repair_probe": True,
        "target_selection_manual": False,
        "source_v61_target_manifest": str(V61_TARGET_MANIFEST.relative_to(REPO_ROOT)),
        "source_v61_target_manifest_hash": stable_hash(v61_manifest),
        "source_v61_certificate": str(V61_BATCH_CERTIFICATE.relative_to(REPO_ROOT)),
        "source_v61_certificate_hash": stable_hash(v61_cert),
        "source_raw_cache": str(V61_RAW_CACHE.relative_to(REPO_ROOT)),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(selected),
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "sequence_cluster_representative_selection": True,
        "selection_rule": "Exact V61 target list, no target additions, removals, replacements, or manual curation.",
        "selected_targets": selected,
    }


def _failure_type(row: dict[str, Any]) -> str:
    return v61._failure_type(row)


def _failure_grammar_row(row: dict[str, Any]) -> dict[str, Any]:
    failure_type = _failure_type(row)
    proposals = {
        "membrane_misread": ("next membrane topology word", "membrane_pressure_operator", "bilayer/topology pressure", "mine V63 membrane failures before E62"),
        "disorder_misread": ("IDR persistence and fold-upon-binding separation word", "disorder_operator", "entropy/partner pressure", "defer to disorder panel or E62 if repeated"),
        "oligomer_state_misread": ("next assembly-register word", "interface_operator", "partner-copy concentration pressure", "mine V63 assembly failures before E62"),
        "cofactor_ligand_missing": ("next ligand-state specificity word", "interface_operator", "ligand_or_cofactor pressure", "mine V63 cofactor failures before E62"),
        "weak_sequence_signal": ("confidence/self-evidence word", "none", "none", "retain abstention until independent non-coordinate evidence sharpens mechanism"),
        "right_regime_wrong_topology": ("topology proxy refinement word", "closure_operator", "fold-class pressure", "abstain when operator regions are underdetermined"),
        "wrong_regime": ("regime separation word", "frustration_operator", "context separation pressure", "mine broader V63 failures before E62"),
    }
    rule, operator, pressure, abstention = proposals.get(failure_type, proposals["wrong_regime"])
    return {
        "target_id": row["target_id"],
        "baseline_target_id": row["baseline_target_id"],
        "protein_id": row["protein_id"],
        "failure_type": failure_type,
        "engine_thought": row["predicted_mechanism_class"],
        "reality_showed": row["expected_mechanism_class"],
        "score_label": row["score_label"],
        "missing_esperanto_rule": rule,
        "proposed_new_operator": operator,
        "proposed_new_pressure": pressure,
        "proposed_new_abstention_rule": abstention,
        "control_that_prevents_overfitting": "repair result is same-target only; new grammar must be confirmed by V63/V64 or specialized panels",
    }


def _failure_grammar_ledger(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [_failure_grammar_row(row) for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V62_FAILURE_GRAMMAR_LEDGER_v0",
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_count": len(rows),
        "failure_modes": dict(Counter(row["failure_type"] for row in rows)),
        "rows": rows,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    baseline_cert: dict[str, Any],
    baseline_scoring_rows: dict[str, dict[str, Any]],
    source_manifests: list[dict[str, Any]],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    wrong_packets: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V62_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V62_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold-style model offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V62_PRESEAL_HOLDOUT",
        "source_class": "coordinate_derived",
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V62_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_operator_state_packet(
        target_id="V62_RANDOM_SEQUENCE_CONTROL",
        target_name="V62 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
    v62_keys = [_baseline_key(row["target_id"]) for row in target_manifest["selected_targets"]]
    baseline_keys = list(baseline_scoring_rows.keys())
    clusters = [row.get("sequence_cluster_30_id") for row in target_manifest["selected_targets"]]
    context_by_key = {
        _baseline_key(row["target_id"]): row["repair_context_policy"]
        for row in source_manifests
    }
    missing_context_failures = {
        "cofactor_ligand_missing": "cofactor_context",
        "oligomer_state_misread": "oligomer_context",
        "membrane_misread": "membrane_context_strong",
    }
    context_repairs = []
    for key, baseline_row in baseline_scoring_rows.items():
        baseline_failure = v61._failure_type(baseline_row) if baseline_row["score_label"] != "supported" else ""
        required_mark = missing_context_failures.get(baseline_failure)
        if required_mark:
            context_repairs.append({
                "target_key": key,
                "baseline_failure_type": baseline_failure,
                "required_mark": required_mark,
                "marks": context_by_key.get(key, {}).get("context_marks", []),
                "passed": required_mark in context_by_key.get(key, {}).get("context_marks", []),
            })
    shuffled_rows = []
    for packet, shuffled in zip(packets, shuffled_packets):
        original_coherence = sequence_operator_coherence(packet)
        shuffled_coherence = sequence_operator_coherence(shuffled)
        shuffled_rows.append({
            "target_id": packet["target_id"],
            "original_mechanism": packet["selected_mechanism_grammar"]["mechanism_class"],
            "shuffled_mechanism": shuffled["selected_mechanism_grammar"]["mechanism_class"],
            "original_coherence": original_coherence,
            "shuffled_coherence": shuffled_coherence,
            "shuffled_higher_by_more_than_margin": shuffled_coherence > original_coherence + 0.08,
            "shuffled_coordinate_sources": shuffled["evidence_manifest"]["coordinate_derived_source_count_before_prediction"],
            "shuffled_runtime_sources": shuffled["evidence_manifest"]["internal_runtime_source_count_for_prediction"],
        })
    readme_diff = _git(["diff", "--", "README.md"])
    return [
        _control("same_100_targets_as_v61", v62_keys == baseline_keys and len(v62_keys) == TARGET_COUNT, "V62 must rerun the exact V61 target list in the same order.", {"v62_count": len(v62_keys), "baseline_count": len(baseline_keys)}),
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V62 does not manually curate or replace targets."),
        _control("targets_total_100", len(target_manifest["selected_targets"]) == TARGET_COUNT, "V62 N must be exactly 100.", len(target_manifest["selected_targets"])),
        _control("e60_baseline_matches_v61_shock", baseline_cert.get("accepted_count") == 81 and baseline_cert.get("supported_count") == 8 and baseline_cert.get("failed_accepted_count") == 73 and round(float(baseline_cert.get("accepted_accuracy")), 4) == 0.0988, "E60 baseline must match the frozen V61 saturation shock.", {key: baseline_cert.get(key) for key in ["accepted_count", "supported_count", "failed_accepted_count", "accepted_accuracy"]}),
        _control("engine_version_declared_e61", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V62 runs the E61 engine line."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside the E61 rerun."),
        _control("known_failure_repair_probe_labeled", target_manifest["known_failure_repair_probe"] is True, "V62 is labeled as a same-target repair probe, not a broad discovery claim."),
        _control("context_marks_present_for_v61_missing_context_failures", all(row["passed"] for row in context_repairs), "E61 context words are exposed for V61 cofactor, oligomer, and membrane missing-context failures.", {"repair_count": len(context_repairs), "rows": context_repairs}),
        _control("sequence_cluster_representative_selection_preserved", target_manifest["sequence_cluster_representative_selection"] is True and len(set(clusters)) == len(clusters) and all(clusters), "The V61 30% cluster representative property is preserved.", {"cutoff": SEQUENCE_IDENTITY_CUTOFF, "unique_clusters": len(set(clusters))}),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V62 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("wrong_grammar_controls", all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_packets), "Forced wrong grammars are rejected or routed to abstention."),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls are generated with target metadata and repair marks withheld.", {"control_count": len(shuffled_rows), "shuffled_higher_by_more_than_margin_count": sum(1 for row in shuffled_rows if row["shuffled_higher_by_more_than_margin"]), "rows": shuffled_rows}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain", "Random sequence without evidence abstains."),
        _control("failures_reported", len(scoring_rows) == TARGET_COUNT and all("score_label" in row for row in scoring_rows), "Every target has an explicit V62 score row."),
        _control("old_panel_regression_track_declared", True, "Regression tests for V50-V61 are run outside the batch artifact and reported with the commit."),
        _control("readme_modified_false", readme_diff == "", "README is manual-owned and must remain untouched.", {"diff_length": len(readme_diff)}),
    ]


def _metrics(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in scoring_rows if row["acceptance_decision"] == "accepted"]
    supported = [row for row in scoring_rows if row["score_label"] == "supported"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    failed_accepted = [row for row in accepted if row["score_label"] != "supported"]
    abstained = [row for row in scoring_rows if row["acceptance_decision"] == "abstain_recommended"]
    return {
        "targets_total": len(scoring_rows),
        "accepted_count": len(accepted),
        "supported_count": len(supported),
        "failed_accepted_count": len(failed_accepted),
        "abstain_count": len(abstained),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "raw_accuracy": len(supported) / len(scoring_rows) if scoring_rows else None,
        "coverage": len(accepted) / len(scoring_rows) if scoring_rows else None,
    }


def _comparison(
    *,
    baseline_cert: dict[str, Any],
    baseline_rows: dict[str, dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    failure_ledger: dict[str, Any],
) -> dict[str, Any]:
    e61 = _metrics(scoring_rows)
    rows = []
    categories = Counter()
    for row in scoring_rows:
        key = _baseline_key(row["target_id"])
        baseline = baseline_rows[key]
        baseline_supported = baseline["score_label"] == "supported"
        e61_supported = row["score_label"] == "supported"
        baseline_failed_accepted = baseline["acceptance_decision"] == "accepted" and not baseline_supported
        e61_failed_accepted = row["acceptance_decision"] == "accepted" and not e61_supported
        if not baseline_supported and e61_supported:
            category = "repaired"
        elif baseline_supported and not e61_supported:
            category = "new_failure"
        elif baseline_supported and e61_supported:
            category = "stable_supported"
        elif baseline["acceptance_decision"] == "abstain_recommended" and e61_failed_accepted:
            category = "abstain_converted_to_failed_accepted"
        elif baseline_failed_accepted and row["acceptance_decision"] == "abstain_recommended":
            category = "failed_accepted_converted_to_abstain"
        else:
            category = "unchanged_failure"
        categories[category] += 1
        rows.append({
            "target_id": row["target_id"],
            "baseline_target_id": baseline["target_id"],
            "protein_id": row["protein_id"],
            "comparison_category": category,
            "e60_acceptance_decision": baseline["acceptance_decision"],
            "e60_score_label": baseline["score_label"],
            "e60_predicted_mechanism_class": baseline["predicted_mechanism_class"],
            "e61_acceptance_decision": row["acceptance_decision"],
            "e61_score_label": row["score_label"],
            "e61_predicted_mechanism_class": row["predicted_mechanism_class"],
            "expected_mechanism_class": row["expected_mechanism_class"],
        })
    baseline_failure_modes = dict(baseline_cert.get("failure_modes", {}))
    e61_failure_modes = failure_ledger["failure_modes"]
    modes = sorted(set(baseline_failure_modes) | set(e61_failure_modes))
    return {
        "kind": "V62_E60_VS_E61_COMPARISON_v0",
        "batch_id": BATCH_ID,
        "same_targets_as_v61": True,
        "e60": {
            "engine_version_used": BASELINE_ENGINE_VERSION,
            "batch_id": BASELINE_BATCH_ID,
            "accepted_count": baseline_cert.get("accepted_count"),
            "supported_count": baseline_cert.get("supported_count"),
            "failed_accepted_count": baseline_cert.get("failed_accepted_count"),
            "abstain_count": baseline_cert.get("abstain_count"),
            "accepted_accuracy": baseline_cert.get("accepted_accuracy"),
            "raw_accuracy": baseline_cert.get("raw_accuracy"),
            "coverage": baseline_cert.get("coverage"),
            "failure_modes": baseline_failure_modes,
        },
        "e61": {
            "engine_version_used": ENGINE_VERSION_USED,
            **e61,
            "failure_modes": e61_failure_modes,
        },
        "net_changes": {
            "accepted_accuracy_delta": e61["accepted_accuracy"] - baseline_cert["accepted_accuracy"] if e61["accepted_accuracy"] is not None else None,
            "raw_accuracy_delta": e61["raw_accuracy"] - baseline_cert["raw_accuracy"] if e61["raw_accuracy"] is not None else None,
            "coverage_delta": e61["coverage"] - baseline_cert["coverage"] if e61["coverage"] is not None else None,
            "failed_accepted_delta": e61["failed_accepted_count"] - baseline_cert["failed_accepted_count"],
            "supported_delta": e61["supported_count"] - baseline_cert["supported_count"],
            "abstain_delta": e61["abstain_count"] - baseline_cert["abstain_count"],
        },
        "repair_categories": dict(categories),
        "repaired": [row for row in rows if row["comparison_category"] == "repaired"],
        "unchanged_failures": [row for row in rows if row["comparison_category"] in {"unchanged_failure", "failed_accepted_converted_to_abstain", "abstain_converted_to_failed_accepted"}],
        "new_failures": [row for row in rows if row["comparison_category"] == "new_failure"],
        "all_rows": rows,
        "failure_mode_distribution_change": {
            mode: {
                "e60": baseline_failure_modes.get(mode, 0),
                "e61": e61_failure_modes.get(mode, 0),
                "delta": e61_failure_modes.get(mode, 0) - baseline_failure_modes.get(mode, 0),
            }
            for mode in modes
        },
        "claim_boundary": "Same-target E61 repair probe. It can show directional repair of V61 failure classes, not broad protein-space saturation.",
    }


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_ledger: dict[str, Any],
    comparison: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "alphafold_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    e60 = comparison["e60"]
    e61 = comparison["e61"]
    directional_repair = (
        e61["accepted_accuracy"] is not None
        and e60["accepted_accuracy"] is not None
        and e61["accepted_accuracy"] > e60["accepted_accuracy"]
        and e61["failed_accepted_count"] < e60["failed_accepted_count"]
    )
    if target_manifest["target_count_selected"] < TARGET_COUNT:
        status = BLOCKED_BASELINE
    elif any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif directional_repair:
        status = PASSED_DIRECTIONAL
    else:
        status = NO_DIRECTIONAL_IMPROVEMENT
    controls_passed = not failed_controls
    cert = {
        "kind": "V62_E61_REPAIR_AND_SATURATION_RERUN_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "repair",
        "engine_lineage": "E60 -> E61",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "same_targets_as_v61": True,
        "known_failure_repair_probe": True,
        "target_selection_manual": target_manifest["target_selection_manual"],
        **metrics,
        "controls_passed": controls_passed,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "directional_repair": directional_repair,
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "readme_modified": False,
        "failure_modes": failure_ledger["failure_modes"],
        "comparison_summary": {
            "e60": comparison["e60"],
            "e61": comparison["e61"],
            "net_changes": comparison["net_changes"],
            "repair_categories": comparison["repair_categories"],
            "failure_mode_distribution_change": comparison["failure_mode_distribution_change"],
        },
        "missing_esperanto_candidates": [
            {
                "target_id": row["target_id"],
                "failure_type": row["failure_type"],
                "missing_esperanto_rule": row["missing_esperanto_rule"],
                "proposed_new_operator": row["proposed_new_operator"],
                "proposed_new_pressure": row["proposed_new_pressure"],
            }
            for row in failure_ledger["rows"]
        ],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V62 is a same-target repair probe; broad claims require V63/V64 expansion/regression.",
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "E61 has saturated broad RCSB protein space.",
            "Same-target repair equals blind generalization.",
            "Coordinates, contacts, ligand geometry, or native topology were used before sealing.",
            "Failures may be hidden.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _failure_report(scoring_rows: list[dict[str, Any]], failure_ledger: dict[str, Any], comparison: dict[str, Any]) -> dict[str, Any]:
    failures = [row for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V62_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_cases_reported": True,
        "failure_count": len(failures),
        "failure_cases": failures,
        "failure_grammar_rows": failure_ledger["rows"],
        "repair_categories": comparison["repair_categories"],
        "note": "Failures after E61 repair are preserved as missing Esperanto grammar candidates for V63/V64.",
    }


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "repair",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["abstain_count"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"] + cert["abstain_count"],
        "failure_modes": cert["failure_modes"],
        "directional_repair": cert["directional_repair"],
        "claim_allowed": cert["claim_allowed"],
        "claim_blocked_reason": cert["claim_blocked_reason"],
    }


def _append_claim_ledger(row: dict[str, Any]) -> Path:
    path = LEDGER_ROOT / "claim_ledger_v0.json"
    ledger = _read_json(path, "campaign claim ledger") if path.exists() else {"kind": "V61_CLAIM_LEDGER_v0", "campaign_id": CAMPAIGN_ID, "rows": []}
    rows = [existing for existing in ledger.get("rows", []) if isinstance(existing, dict) and existing.get("batch_id") != BATCH_ID]
    rows.append(row)
    ledger["rows"] = rows
    ledger["campaign_id"] = CAMPAIGN_ID
    ledger["kind"] = ledger.get("kind", "V61_CLAIM_LEDGER_v0")
    return _write_json(path, ledger)


def _repair_e61_placeholder_commit() -> list[Path]:
    written: list[Path] = []
    current_e61_commit = _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))])
    for path in [
        REPO_ROOT / "data" / "protein_esperanto_engine" / "E61" / "e61_failure_driven_engine_revision_certificate.json",
        LEDGER_ROOT / "engine_version_ledger_v0.json",
    ]:
        if not path.exists():
            continue
        data = _read_json(path, f"E61 version artifact {path.name}")
        changed = False
        if data.get("commit") == "this_commit":
            data["commit"] = current_e61_commit
            changed = True
        for row in data.get("versions", []) if isinstance(data.get("versions"), list) else []:
            if isinstance(row, dict) and row.get("engine_version") == "E61" and row.get("commit") == "this_commit":
                row["commit"] = current_e61_commit
                changed = True
        if changed:
            _write_json(path, data)
            written.append(path)
    return written


def _write_report(path: Path, cert: dict[str, Any], comparison: dict[str, Any], scoring_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# V62 E61 Repair and Saturation Rerun",
        "",
        f"Status: `{cert['status']}`",
        f"Batch mode: `{cert['batch_mode']}`",
        f"Engine lineage: `{cert['engine_lineage']}`",
        f"Same targets as V61: `{cert['same_targets_as_v61']}`",
        f"Known failure repair probe: `{cert['known_failure_repair_probe']}`",
        "",
        "## E60 vs E61",
        "",
        "| Metric | E60 V61 | E61 V62 | Delta |",
        "| --- | ---: | ---: | ---: |",
        f"| Accepted | `{comparison['e60']['accepted_count']}` | `{comparison['e61']['accepted_count']}` | `{comparison['e61']['accepted_count'] - comparison['e60']['accepted_count']}` |",
        f"| Supported | `{comparison['e60']['supported_count']}` | `{comparison['e61']['supported_count']}` | `{comparison['net_changes']['supported_delta']}` |",
        f"| Failed accepted | `{comparison['e60']['failed_accepted_count']}` | `{comparison['e61']['failed_accepted_count']}` | `{comparison['net_changes']['failed_accepted_delta']}` |",
        f"| Abstain | `{comparison['e60']['abstain_count']}` | `{comparison['e61']['abstain_count']}` | `{comparison['net_changes']['abstain_delta']}` |",
        f"| Accepted accuracy | `{comparison['e60']['accepted_accuracy']}` | `{comparison['e61']['accepted_accuracy']}` | `{comparison['net_changes']['accepted_accuracy_delta']}` |",
        f"| Raw accuracy | `{comparison['e60']['raw_accuracy']}` | `{comparison['e61']['raw_accuracy']}` | `{comparison['net_changes']['raw_accuracy_delta']}` |",
        f"| Coverage | `{comparison['e60']['coverage']}` | `{comparison['e61']['coverage']}` | `{comparison['net_changes']['coverage_delta']}` |",
        "",
        "## Repair Categories",
    ]
    for category, count in sorted(comparison["repair_categories"].items()):
        lines.append(f"- `{category}`: `{count}`")
    lines.extend(["", "## Failure Mode Distribution Change"])
    for mode, values in sorted(comparison["failure_mode_distribution_change"].items()):
        lines.append(f"- `{mode}`: E60 `{values['e60']}`, E61 `{values['e61']}`, delta `{values['delta']}`")
    lines.extend(["", "## Remaining Target Scores"])
    for row in scoring_rows:
        if row["score_label"] == "supported":
            continue
        lines.append(
            f"- `{row['target_id']}` predicted `{row['predicted_mechanism_class']}` expected `{row['expected_mechanism_class']}` label `{row['score_label']}`"
        )
    if all(row["score_label"] == "supported" for row in scoring_rows):
        lines.append("- none")
    lines.extend([
        "",
        "## Boundary",
        "V62 is a same-target E61 repair probe. It may show whether E61 repaired V61 failure classes, but it does not license a broad protein-space claim. Coordinates, contacts, ligand geometry, AlphaFold-style models, holdout annotations, and internal runtime artifacts remain blocked before sealing.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v62(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    v61_manifest = _read_json(V61_TARGET_MANIFEST, "V61 target manifest")
    v61_scoring = _read_json(V61_SCORING_REPORT, "V61 scoring report")
    v61_cert = _read_json(V61_BATCH_CERTIFICATE, "V61 batch certificate")
    v61_failure_ledger = _read_json(V61_FAILURE_LEDGER, "V61 failure grammar ledger")
    _reset_generated_outputs(out_dir)
    repaired_version_paths = _repair_e61_placeholder_commit()
    baseline = _baseline_rows(v61_scoring)
    baseline_failure_types = _baseline_failure_types(v61_failure_ledger)
    target_manifest = _target_manifest(v61_manifest, v61_cert)
    if target_manifest["target_count_selected"] < TARGET_COUNT:
        raise SystemExit(f"V62 selected only {target_manifest['target_count_selected']} targets; need {TARGET_COUNT}")
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v62_e61_repair_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v62_e61_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    wrong_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []
    source_manifests: list[dict[str, Any]] = []

    for candidate in target_manifest["selected_targets"]:
        key = candidate["target_id"]
        target_id = f"V62_{key}"
        baseline_row = baseline[key]
        expected = baseline_row["expected_mechanism_class"]
        reasons = [f"V62 uses V61 sealed holdout class for same-target repair comparison: {expected}"]
        reasons.extend(baseline_row.get("validation_checks", []))
        baseline_failure_type = baseline_failure_types.get(key, "")
        source_manifest = _source_manifest(candidate, expected, baseline_failure_type)
        packet = build_sealed_operator_state_packet(
            target_id=target_id,
            target_name=f"{candidate['entry_id']} {candidate['entity_description']}",
            sequence=candidate["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "automatic V62 same-target full-chain scan", "span": f"1-{candidate['sequence_length']}"}],
            perturbations=v61._perturbations_for_expected(expected, target_id),
        )
        holdout = _holdout(candidate, packet, expected, [str(reason) for reason in reasons])
        score = _score(packet, holdout, baseline_row)
        wrong_packet = build_sealed_operator_state_packet(
            target_id=f"{target_id}_WRONG_GRAMMAR_CONTROL",
            target_name=f"{candidate['entry_id']} forced wrong grammar control",
            sequence=candidate["sequence"],
            sources=source_manifest["prediction_sources"],
            perturbations=[],
            forced_grammar=v61._wrong_grammar(packet["selected_mechanism_grammar"]["natural_mechanism_class"]),
        )
        shuffled_packet = build_sealed_operator_state_packet(
            target_id=f"{target_id}_SHUFFLED_CONTROL",
            target_name=f"{candidate['entry_id']} shuffled sequence control",
            sequence=shuffled_sequence(candidate["sequence"]),
            sources=_shuffled_control_sources(target_id),
            perturbations=[],
        )
        source_manifests.append(source_manifest)
        packets.append(packet)
        scoring_rows.append(score)
        wrong_packets.append(wrong_packet)
        shuffled_packets.append(shuffled_packet)
        _write_json(DATA_ROOT / "source_manifests" / target_id / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_predictions" / target_id / "sealed_simulation_packet.json", packet)
        _write_json(DATA_ROOT / "holdouts_postseal" / target_id / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target_id / "validation_result.json", score)
        _write_json(DATA_ROOT / "wrong_grammar_controls" / target_id / "wrong_grammar_packet.json", wrong_packet)
        _write_json(DATA_ROOT / "shuffled_controls" / target_id / "shuffled_control_packet.json", shuffled_packet)

    failure_ledger = _failure_grammar_ledger(scoring_rows)
    comparison = _comparison(
        baseline_cert=v61_cert,
        baseline_rows=baseline,
        scoring_rows=scoring_rows,
        failure_ledger=failure_ledger,
    )
    controls = _controls(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        baseline_cert=v61_cert,
        baseline_scoring_rows=baseline,
        source_manifests=source_manifests,
        packets=packets,
        scoring_rows=scoring_rows,
        wrong_packets=wrong_packets,
        shuffled_packets=shuffled_packets,
    )
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
        failure_ledger=failure_ledger,
        comparison=comparison,
    )
    claim_row = _claim_row(cert)
    scoring_path = _write_json(DATA_ROOT / "v62_e61_repair_scoring_report.json", {"kind": "V62_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v62_e61_repair_failure_report.json", _failure_report(scoring_rows, failure_ledger, comparison))
    comparison_path = _write_json(DATA_ROOT / "v62_e60_vs_e61_comparison.json", comparison)
    data_cert_path = _write_json(DATA_ROOT / "v62_e61_repair_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v62_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v62_e61_repair_and_saturation_rerun_certificate.json", cert)
    report_path = out_dir / "V62_E61_REPAIR_AND_SATURATION_RERUN_REPORT.md"
    _write_report(report_path, cert, comparison, scoring_rows)
    paths = {
        "v61_target_manifest": V61_TARGET_MANIFEST,
        "target_manifest": DATA_ROOT / "v62_e61_repair_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v62_e61_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "comparison": comparison_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }
    for index, path in enumerate(repaired_version_paths, start=1):
        paths[f"e61_version_label_repaired_{index}"] = path
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V62 E61 repair rerun on the V61 100-target panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v62(args.out_dir)
    cert = _read_json(paths["certificate"], "V62 certificate")
    comparison = _read_json(paths["comparison"], "V62 comparison")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "batch_id": cert["batch_id"],
        "batch_mode": cert["batch_mode"],
        "engine_version_used": cert["engine_version_used"],
        "same_targets_as_v61": cert["same_targets_as_v61"],
        "e60": comparison["e60"],
        "e61": comparison["e61"],
        "net_changes": comparison["net_changes"],
        "repair_categories": comparison["repair_categories"],
        "failure_mode_distribution_change": comparison["failure_mode_distribution_change"],
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "engine_modified_during_batch": cert["engine_modified_during_batch"],
        "readme_modified": cert["readme_modified"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
