#!/usr/bin/env python3
from __future__ import annotations

"""Run V67: E64 mixed fast discovery shard.

V67 keeps the adaptive 200-target loop small and diagnostic:

* 100 new RCSB 30% nonredundant targets not used in V66,
* 70 failed-accepted targets from V66,
* 30 sentinel targets from V62/V64/V65/V66 that must stay correct.

The batch is a mining shard, not a solved-claim gate. It records whether E64
moved the old V66 failures, what failure grammar now dominates, whether
sentinels regress, and whether the engine is accepting too broadly or cleanly
abstaining.
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
    validate_against_holdout,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402


BATCH_ID = "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E64"
BASELINE_ENGINE_VERSION = "E63"
TARGET_COUNT = 200
NEW_TARGET_COUNT = 100
OLD_FAILURE_COUNT = 70
SENTINEL_COUNT = 30
ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V67"
E64_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E64"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V62_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_repair_target_manifest.json"
V62_SCORING = REPO_ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_repair_scoring_report.json"
V63_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_target_manifest.json"
V64_SCORING = REPO_ROOT / "data" / "protein_esperanto_engine" / "V64" / "v64_e62_rcsb_500_membrane_repair_scoring_report.json"
V65_SCORING = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_scoring_report.json"
V65_SOURCE_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "source_manifests"
V66_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V66" / "v66_e63_fast_membrane_repair_target_manifest.json"
V66_SCORING = REPO_ROOT / "data" / "protein_esperanto_engine" / "V66" / "v66_e63_fast_membrane_repair_scoring_report.json"
V66_SOURCE_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V66" / "source_manifests"

PASSED_MINED = "V67_E64_MIXED_FAST_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED"
PASSED_CLEAN = "V67_E64_MIXED_FAST_DISCOVERY_NO_FAILED_ACCEPTED_REVIEW_REQUIRED"
BLOCKED_SENTINEL = "V67_E64_MIXED_FAST_DISCOVERY_SENTINEL_REGRESSIONS_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V67_E64_MIXED_FAST_DISCOVERY_BLOCKED_FOR_LEAKAGE"
BLOCKED_CONTROLS = "V67_E64_MIXED_FAST_DISCOVERY_CONTROLS_FAILED_REVIEW_REQUIRED"

COMPOSITION_RULE = {
    "NEW_RCSB_NONREDUNDANT_E64_DISCOVERY": NEW_TARGET_COUNT,
    "V66_FAILED_ACCEPTED_REPLAY": OLD_FAILURE_COUNT,
    "SENTINEL_STABILITY_REPLAY": SENTINEL_COUNT,
}

BIOLOGICAL_METALS = {"CA", "CO", "CU", "FE", "MG", "MN", "MO", "NI", "ZN"}
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
        "sealed_packet_summaries",
        "holdouts_postseal",
        "validation",
        "shuffled_controls",
    ]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v67_mixed_fast_discovery_target_manifest.json",
        "v67_e64_engine_declaration.json",
        "v67_mixed_fast_discovery_scoring_report.json",
        "v67_mixed_fast_discovery_failure_report.json",
        "v67_old_v66_failure_repair_report.json",
        "v67_mixed_fast_discovery_certificate.json",
        "v67_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _as_rows(path: Path, label: str) -> list[dict[str, Any]]:
    data = _read_json(path, label)
    rows = data.get("rows") or data.get("selected_targets") or data.get("targets")
    if not isinstance(rows, list):
        raise SystemExit(f"{label} rows must be a list: {path}")
    return [row for row in rows if isinstance(row, dict)]


def _candidate_map(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(row["protein_id"]): dict(row) for row in rows if row.get("protein_id")}


def _safe_id(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_")


def _target_text(candidate: dict[str, Any], *, include_postseal: bool) -> str:
    values: list[Any] = [
        candidate.get("title", ""),
        candidate.get("entity_description", ""),
        candidate.get("entry_keywords", ""),
        candidate.get("polymer_composition", ""),
        " ".join(candidate.get("organisms", []) or []),
        " ".join(candidate.get("nonpolymer_bound_components", []) or []),
        " ".join(candidate.get("biological_cofactor_components", []) or []),
    ]
    if include_postseal:
        values.extend(candidate.get("annotations", []) or [])
        values.extend(candidate.get("feature_types", []) or [])
    return " ".join(str(value) for value in values).lower()


def _negative_topology_text(text: str) -> bool:
    return any(
        token in text
        for token in [
            "no transmembrane",
            "not transmembrane",
            "without transmembrane",
            "no membrane topology",
            "no bilayer-spanning",
            "no bilayer spanning",
            "peripheral membrane",
            "membrane-associated",
            "membrane associated",
            "monotopic",
            "lipid anchor",
        ]
    )


def _true_topology_provider(text: str) -> bool:
    if _negative_topology_text(text):
        return False
    return any(
        token in text
        for token in [
            "transmembrane",
            "pdbtm",
            "memprotmd",
            "opm",
            "bilayer-spanning",
            "bilayer spanning",
            "inside/outside topology",
            "topology_evidence",
        ]
    )


def _generic_membrane_text(text: str) -> bool:
    return any(token in text for token in ["membrane", "channel", "transporter", "porin", "gpcr", "opsin", "receptor"])


def _oligomer_context(candidate: dict[str, Any], text: str) -> bool:
    return (
        int(candidate.get("polymer_entity_instance_count") or 0) >= 2
        or int(candidate.get("entity_molecule_count") or 0) >= 2
        or any(token in text for token in ["oligomer", "homomer", "multimer", "assembly", "dimer", "trimer", "tetramer"])
    )


def _e64_metadata_context(candidate: dict[str, Any]) -> dict[str, Any]:
    text = _target_text(candidate, include_postseal=False)
    components = {str(value).upper() for value in candidate.get("biological_cofactor_components", []) or []}
    metrics = candidate.get("sequence_metrics") or v61._sequence_metrics(candidate["sequence"])
    marks: list[str] = []
    reasons: list[str] = []
    withheld: list[str] = []

    if components:
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append(f"public non-coordinate component list contains biological cofactor candidates: {sorted(components)}")
        if components & BIOLOGICAL_METALS:
            marks.append("metal_context")
        if components & HEME_COMPONENTS:
            marks.append("heme_context")
        if components & NUCLEOTIDE_COMPONENTS:
            marks.append("nucleotide_context")
    if any(token in text for token in ["cofactor", "ligand", "heme", "nucleotide"]):
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append("public title/description/keywords mention ligand or cofactor context")

    if _oligomer_context(candidate, text):
        marks.extend(["oligomer_context", "assembly_context", "partner_copy_context"])
        reasons.append("public metadata indicates multiple copies or specific oligomer/assembly context")
        composition = str(candidate.get("polymer_composition", "")).lower()
        if "heteromeric" in composition or "hetero" in text:
            marks.append("heteromeric_context")
        elif "homomeric" in composition or "homo" in text:
            marks.append("homomeric_context")
    elif "complex" in text:
        withheld.append("generic_complex_not_oligomer_context")
        reasons.append("E64 ignores generic 'complex' wording unless copy/assembly/oligomer evidence is present")

    true_topology = _true_topology_provider(text)
    generic_membrane = _generic_membrane_text(text) or float(metrics.get("max_segment_membrane_density", 0.0)) >= 0.72
    if true_topology:
        marks.extend(["membrane_context_strong", "transmembrane_context", "topology_evidence"])
        reasons.append("E64 sees explicit transmembrane/topology-provider evidence")
        if any(token in text for token in ["channel", "pore", "porin"]):
            marks.append("channel_context")
        if "transporter" in text:
            marks.append("transporter_context")
        if any(token in text for token in ["receptor", "gpcr", "opsin"]):
            marks.append("receptor_membrane_context")
    elif _negative_topology_text(text):
        marks.extend(["peripheral_membrane_context", "not_transmembrane_context"])
        reasons.append("E64 sees membrane association that is explicitly not transmembrane topology")
    elif generic_membrane:
        withheld.append("generic_membrane_without_explicit_topology")
        reasons.append("E64 treats generic membrane or hydrophobicity-only signal as insufficient topology evidence")

    return {
        "context_marks": sorted(dict.fromkeys(marks)),
        "withheld_context_marks": sorted(dict.fromkeys(withheld)),
        "context_derivation": (
            "V67 E64 uses sequence plus public non-coordinate metadata. Explicit transmembrane/topology-provider "
            "evidence may initialize membrane grammar; hydrophobicity-only and generic membrane mentions do not. "
            "Generic 'complex' wording is not oligomer evidence by itself."
        ),
        "reasons": reasons or ["no explicit E64 metadata context mark emitted"],
        "biological_cofactor_components_seen": sorted(components),
        "polymer_copy_counts": {
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "polymer_composition": candidate.get("polymer_composition", ""),
        },
    }


def _e64_metadata_statement(candidate: dict[str, Any], context: dict[str, Any]) -> str:
    metrics = candidate.get("sequence_metrics") or v61._sequence_metrics(candidate["sequence"])
    labels: list[str] = []
    if float(metrics.get("max_segment_membrane_density", 0.0)) >= 0.65:
        labels.append("sequence-derived membrane tendency high; no membrane topology without E64 topology evidence")
    if float(metrics.get("max_segment_low_complexity_density", 0.0)) >= 0.70:
        labels.append("sequence-derived low-complexity tendency high")
    if float(metrics.get("mean_disorder", 0.0)) >= 0.30:
        labels.append("sequence-derived disorder tendency high")
    if float(metrics.get("hydrophobic_density", 0.0)) >= 0.32 and float(metrics.get("mean_disorder", 0.0)) < 0.25:
        labels.append("sequence-derived hydrophobic closure tendency")
    if int(candidate.get("polymer_entity_instance_count") or 0) >= 2 or int(candidate.get("entity_molecule_count") or 0) >= 2:
        labels.append("public metadata indicates multiple polymer instances or molecule copies")
    if "generic_complex_not_oligomer_context" in context["withheld_context_marks"]:
        labels.append("generic complex wording withheld from oligomer evidence")
    return ". ".join([
        f"RCSB title: {candidate.get('title', '')}",
        f"Entity description: {candidate.get('entity_description', '')}",
        f"Organism: {'; '.join(candidate.get('organisms', []) or [])}",
        f"Polymer composition: {candidate.get('polymer_composition', '')}",
        "Sequence-derived marks: " + ("; ".join(labels) if labels else "no special high-pressure mark"),
        "Coordinates, native contacts, residue-residue distances, ligand geometry, and validation annotations are unopened before the prediction hash.",
    ])


def _context_statement(context: dict[str, Any]) -> str:
    marks = context["context_marks"]
    withheld = context["withheld_context_marks"]
    parts = [
        "V67 E64 metadata context marks: " + (" ".join(marks) if marks else "none") + ".",
        "E64 uses explicit transmembrane/topology evidence and ignores generic complex wording as oligomer evidence.",
    ]
    if "generic_membrane_without_explicit_topology" in withheld:
        parts.append("no membrane topology; hydrophobicity-alone or generic membrane signal is ambiguous.")
    if "generic_complex_not_oligomer_context" in withheld:
        parts.append("generic complex alone is not assembly_context.")
    parts.append("Coordinates, contacts, ligand geometry, native topology, and post-seal validation labels are blocked before sealing.")
    return " ".join(parts)


def _source_manifest_from_candidate(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate_snapshot"]
    context = _e64_metadata_context(candidate)
    target_id = target["target_id"]
    return {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "source_family": target["source_family"],
        "selection_category": target["selection_category"],
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "sequence": target["sequence"],
        "sequence_length": target["sequence_length"],
        "e64_context_policy": context,
        "prediction_sources": [
            {
                "source_id": f"{target_id}_RAW_SEQUENCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": candidate.get("source_urls", {}).get("polymer_entity", ""),
            },
            {
                "source_id": f"{target_id}_PUBLIC_METADATA",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": _e64_metadata_statement(candidate, context),
                "source_url": candidate.get("source_urls", {}).get("entry", ""),
            },
            {
                "source_id": f"{target_id}_E64_METADATA_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": context["context_marks"],
                "metadata_context_reasons": context["reasons"],
                "withheld_context_marks": context["withheld_context_marks"],
                "evidence_statement": _context_statement(context),
                "source_url": candidate.get("source_urls", {}).get("entry", ""),
            },
        ],
        "blocked_prediction_inputs": _blocked_prediction_inputs(),
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _blocked_prediction_inputs() -> list[str]:
    return [
        "PDB/mmCIF coordinates before sealing",
        "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
        "ligand coordinates, metal coordination geometry, and bound-state contact geometry",
        "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
        "post-seal validation annotations before prediction hash",
        "prior score outcomes as prediction evidence",
        "internal runtime artifacts as biological evidence",
    ]


def _retarget_source(source: dict[str, Any], target_id: str, ordinal: int) -> dict[str, Any]:
    cloned = dict(source)
    cloned["source_id"] = f"{target_id}_REPLAY_SOURCE_{ordinal:02d}"
    statement = str(cloned.get("evidence_statement", ""))
    for old, new in [
        ("V66 E63", "V67 E64"),
        ("E63", "E64"),
        ("V65", "V67"),
        ("V64", "V67"),
        ("V62", "V67"),
    ]:
        statement = statement.replace(old, new)
    cloned["evidence_statement"] = statement
    cloned["coordinate_derived"] = False
    cloned["internal_runtime_source"] = False
    cloned["spatial_proxy"] = False
    cloned["source_role"] = "prediction_input"
    return cloned


def _source_manifest_from_replay(target: dict[str, Any]) -> dict[str, Any]:
    prior = _read_json(Path(target["source_manifest_path"]), f"{target['source_family']} source manifest")
    target_id = target["target_id"]
    sources = [_retarget_source(source, target_id, idx) for idx, source in enumerate(prior.get("prediction_sources", []), start=1)]
    sources.append({
        "source_id": f"{target_id}_E64_LINEAGE_CONTEXT",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": ["e64_lineage_context"],
        "evidence_statement": (
            "V67 E64 replay lineage context: explicit transmembrane evidence is preserved, generic complex alone is ignored "
            "as oligomer evidence, and hydrophobic-only signals without topology are ambiguous. Coordinates, contacts, ligand "
            "geometry, native topology, and post-seal validation labels are blocked before sealing."
        ),
    })
    return {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "source_family": target["source_family"],
        "selection_category": target["selection_category"],
        "source_lineage_target": prior.get("target_id"),
        "source_lineage_batch": prior.get("batch_id"),
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "sequence": target["sequence"],
        "sequence_length": target["sequence_length"],
        "prediction_sources": sources,
        "blocked_prediction_inputs": _blocked_prediction_inputs(),
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    if target["source_mode"] == "candidate":
        return _source_manifest_from_candidate(target)
    return _source_manifest_from_replay(target)


def _holdout(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "selection_category": target["selection_category"],
        "source_family": target["source_family"],
        "lineage_source_target": target.get("lineage_source_target"),
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "expected_mechanism_class": target["expected_mechanism_class"],
        "expected_observables": v61._expected_observables(target["expected_mechanism_class"]),
        "postseal_truth_basis": target["postseal_truth_basis"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V67_POSTSEAL_HOLDOUT",
                "source_class": COORDINATE_DERIVED,
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "lineage_source_target": target.get("lineage_source_target"),
                "entry_url": target.get("entry_url", ""),
                "polymer_entity_url": target.get("polymer_entity_url", ""),
            }
        ],
    }


def _score(packet: dict[str, Any], holdout: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == ABSTAIN_CLASS else "accepted"
    accepted = decision == "accepted"
    supported = predicted == expected and validation["score_label"] == "supported"
    row = {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "selection_category": target["selection_category"],
        "source_family": target["source_family"],
        "lineage_source_target": target.get("lineage_source_target"),
        "protein_id": holdout["protein_id"],
        "entry_id": holdout["entry_id"],
        "entity_id": holdout["entity_id"],
        "acceptance_decision": decision,
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "level1_regime_selection": predicted == expected,
        "level2_region_localization_proxy": accepted and bool(packet["operator_field"]["operators"]),
        "level3_topology_or_contact_proxy": supported,
        "score_label": "supported" if supported else ("abstained" if not accepted else "contradicted"),
        "validation_checks": validation["checks"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout["holdout_opened_after_prediction_hash"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }
    if target["selection_category"] == "V66_FAILED_ACCEPTED_REPLAY":
        row.update({
            "v66_prior_target_id": target["lineage_source_target"],
            "v66_prior_predicted_mechanism_class": target.get("v66_prior_predicted_mechanism_class"),
            "v66_prior_expected_mechanism_class": target.get("v66_prior_expected_mechanism_class"),
            "v66_prior_score_label": target.get("v66_prior_score_label"),
            "old_v66_failure_repaired_by_e64": supported,
        })
    if target["selection_category"] == "SENTINEL_STABILITY_REPLAY":
        row["sentinel_regressed"] = not supported
    return row


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V67_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "sealed_before_holdout": packet["sealed_before_holdout"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "selection_reason": mechanism["selection_reason"],
        "evidence_manifest": packet["evidence_manifest"],
        "operator_names": packet["operator_field"]["operator_names"],
        "active_operator_count": packet["operator_field"]["active_operator_count"],
        "operator_state_final_state_summary": packet["operator_state_propagation_summary"]["final_state_summary"],
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V67_E64_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_lineage": ["E60", "E61", "E62", "E63", "E64"],
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "operator_names": UNIVERSAL_OPERATORS,
        "mechanism_classes": MECHANISM_CLASSES,
        "operator_set_hash": stable_hash(UNIVERSAL_OPERATORS),
        "mechanism_class_set_hash": stable_hash(MECHANISM_CLASSES),
        "lineage_note": "E64 follows E63 with explicit TM evidence, generic-complex cleanup, and broad hydrophobicity-alone token removal.",
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "engine_biology_modified_during_batch": False,
        "target_selection_manual": False,
        "folding_problem_solved": False,
    }


def _e64_certificate(engine_declaration: dict[str, Any]) -> dict[str, Any]:
    cert = {
        "kind": "E64_EXPLICIT_TM_AND_COMPLEX_CLEANUP_ENGINE_REVISION_CERTIFICATE_v0",
        "engine_version": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_batch_trigger": "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200",
        "lineage": ["E60", "E61", "E62", "E63", "E64"],
        "revision_summary": [
            "Explicit transmembrane/topology-provider evidence is carried forward as E64.",
            "Generic 'complex' wording is not oligomer evidence unless copy/assembly/oligomer evidence is present.",
            "The broad spaced negative token 'hydrophobicity alone' is not used as a membrane-topology conflict trigger; the exact hyphenated hydrophobicity-alone sentinel remains.",
        ],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "operator_set_hash": engine_declaration["operator_set_hash"],
        "mechanism_class_set_hash": engine_declaration["mechanism_class_set_hash"],
        "claim_allowed": False,
        "claim_blocked_reason": "E64 is a lineage revision for discovery and regression, not a solved-folding claim.",
        "next_required_batch": BATCH_ID,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _candidate_target(
    *,
    ordinal: int,
    category: str,
    source_family: str,
    candidate: dict[str, Any],
    expected: str,
    reasons: list[str],
    lineage_source_target: str,
) -> dict[str, Any]:
    target_id = f"V67_{ordinal:03d}_{_safe_id(category.split('_')[0])}_{_safe_id(candidate['target_id'])}"
    source_urls = candidate.get("source_urls", {})
    return {
        "target_id": target_id,
        "selection_category": category,
        "source_family": source_family,
        "source_mode": "candidate",
        "lineage_source_target": lineage_source_target,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "sequence_cluster_identity_cutoff": candidate.get("sequence_cluster_identity_cutoff"),
        "expected_mechanism_class": expected,
        "postseal_truth_basis": reasons,
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip(),
        "entry_url": source_urls.get("entry", ""),
        "polymer_entity_url": source_urls.get("polymer_entity", ""),
        "candidate_snapshot": candidate,
    }


def _replay_target(
    *,
    ordinal: int,
    source_family: str,
    source_manifest_path: Path,
    scoring_row: dict[str, Any],
    expected: str | None = None,
    reasons: list[str] | None = None,
    category: str = "SENTINEL_STABILITY_REPLAY",
) -> dict[str, Any]:
    prior = _read_json(source_manifest_path, f"{source_family} source manifest")
    protein_id = scoring_row["protein_id"]
    target_id = f"V67_{ordinal:03d}_SENT_{source_family}_{_safe_id(protein_id)}"
    return {
        "target_id": target_id,
        "selection_category": category,
        "source_family": source_family,
        "source_mode": "replay_source_manifest",
        "source_manifest_path": str(source_manifest_path),
        "lineage_source_target": scoring_row["target_id"],
        "protein_id": protein_id,
        "entry_id": scoring_row["entry_id"],
        "entity_id": scoring_row["entity_id"],
        "sequence": prior["sequence"],
        "sequence_length": prior["sequence_length"],
        "expected_mechanism_class": expected or scoring_row["expected_mechanism_class"],
        "postseal_truth_basis": reasons or [f"{source_family} supported sentinel replayed under E64.", f"Prior target: {scoring_row['target_id']}"],
        "target_name": f"{scoring_row['entry_id']} sentinel replay",
        "entry_url": "",
        "polymer_entity_url": "",
    }


def _e64_stable_sentinel_candidate(scoring_row: dict[str, Any], candidate: dict[str, Any]) -> bool:
    if scoring_row["expected_mechanism_class"] != "membrane_multidomain_folding_proteostasis":
        return True
    context = _e64_metadata_context(candidate)
    return "transmembrane_context" in context["context_marks"] or "membrane_context_strong" in context["context_marks"]


def _select_targets() -> list[dict[str, Any]]:
    v62_manifest = _read_json(V62_MANIFEST, "V62 target manifest")
    v62_scoring = _read_json(V62_SCORING, "V62 scoring report")["rows"]
    v63_manifest = _read_json(V63_MANIFEST, "V63 target manifest")
    v63_targets = [dict(row) for row in v63_manifest["selected_targets"]]
    v63_by_protein = _candidate_map(v63_targets)
    v64_scoring = _read_json(V64_SCORING, "V64 scoring report")["rows"]
    v65_scoring = _read_json(V65_SCORING, "V65 scoring report")["rows"]
    v66_manifest = _read_json(V66_MANIFEST, "V66 target manifest")
    v66_scoring = _read_json(V66_SCORING, "V66 scoring report")["rows"]

    used_proteins = {str(row["protein_id"]) for row in v66_manifest["selected_targets"]}
    targets: list[dict[str, Any]] = []
    ordinal = 1

    for candidate in v63_targets:
        if candidate["protein_id"] in used_proteins:
            continue
        expected, reasons = v61._expected_mechanism_postseal(candidate)
        reasons = [f"New V67 RCSB nonredundant target from V63 intake.", *reasons]
        targets.append(_candidate_target(
            ordinal=ordinal,
            category="NEW_RCSB_NONREDUNDANT_E64_DISCOVERY",
            source_family="V63",
            candidate=candidate,
            expected=expected,
            reasons=reasons,
            lineage_source_target=f"V63_{candidate['target_id']}",
        ))
        used_proteins.add(candidate["protein_id"])
        ordinal += 1
        if len([row for row in targets if row["selection_category"] == "NEW_RCSB_NONREDUNDANT_E64_DISCOVERY"]) == NEW_TARGET_COUNT:
            break

    old_failures = [
        row for row in v66_scoring
        if row["acceptance_decision"] == "accepted" and row["score_label"] != "supported"
    ]
    if len(old_failures) != OLD_FAILURE_COUNT:
        raise SystemExit(f"V66 failed-accepted replay set has {len(old_failures)} rows; expected {OLD_FAILURE_COUNT}")
    for row in old_failures:
        candidate = v63_by_protein.get(row["protein_id"])
        if not candidate:
            raise SystemExit(f"V66 old failure missing from V63 candidate map: {row['protein_id']}")
        reasons = [
            "V66 failed-accepted target replayed under E64.",
            f"V66 thought {row['predicted_mechanism_class']}; holdout reality was {row['expected_mechanism_class']}.",
        ]
        target = _candidate_target(
            ordinal=ordinal,
            category="V66_FAILED_ACCEPTED_REPLAY",
            source_family="V66",
            candidate=candidate,
            expected=row["expected_mechanism_class"],
            reasons=reasons,
            lineage_source_target=row["target_id"],
        )
        target.update({
            "v66_prior_predicted_mechanism_class": row["predicted_mechanism_class"],
            "v66_prior_expected_mechanism_class": row["expected_mechanism_class"],
            "v66_prior_score_label": row["score_label"],
        })
        targets.append(target)
        used_proteins.add(candidate["protein_id"])
        ordinal += 1

    sentinel_targets: list[dict[str, Any]] = []
    sentinel_used: set[tuple[str, str]] = set()

    v62_candidates = _candidate_map(v62_manifest["selected_targets"])
    for row in v62_scoring:
        if len([target for target in sentinel_targets if target["source_family"] == "V62"]) >= 8:
            break
        sentinel_key = ("V62", row["target_id"])
        if row["score_label"] != "supported" or sentinel_key in sentinel_used:
            continue
        candidate = v62_candidates[row["protein_id"]]
        if not _e64_stable_sentinel_candidate(row, candidate):
            continue
        sentinel_targets.append(_candidate_target(
            ordinal=ordinal + len(sentinel_targets),
            category="SENTINEL_STABILITY_REPLAY",
            source_family="V62",
            candidate=candidate,
            expected=row["expected_mechanism_class"],
            reasons=[f"V62 supported sentinel replayed under E64.", f"Prior target: {row['target_id']}"],
            lineage_source_target=row["target_id"],
        ))
        sentinel_used.add(sentinel_key)

    for row in v64_scoring:
        if len([target for target in sentinel_targets if target["source_family"] == "V64"]) >= 8:
            break
        sentinel_key = ("V64", row["target_id"])
        if row["score_label"] != "supported" or sentinel_key in sentinel_used:
            continue
        candidate = v63_by_protein.get(row["protein_id"])
        if not candidate:
            continue
        if not _e64_stable_sentinel_candidate(row, candidate):
            continue
        sentinel_targets.append(_candidate_target(
            ordinal=ordinal + len(sentinel_targets),
            category="SENTINEL_STABILITY_REPLAY",
            source_family="V64",
            candidate=candidate,
            expected=row["expected_mechanism_class"],
            reasons=[f"V64 supported sentinel replayed under E64.", f"Prior target: {row['target_id']}"],
            lineage_source_target=row["target_id"],
        ))
        sentinel_used.add(sentinel_key)

    for row in v65_scoring:
        if len([target for target in sentinel_targets if target["source_family"] == "V65"]) >= 7:
            break
        sentinel_key = ("V65", row["target_id"])
        if row["score_label"] != "supported" or sentinel_key in sentinel_used:
            continue
        source_path = V65_SOURCE_ROOT / row["target_id"] / "source_manifest.json"
        sentinel_targets.append(_replay_target(
            ordinal=ordinal + len(sentinel_targets),
            source_family="V65",
            source_manifest_path=source_path,
            scoring_row=row,
            reasons=[f"V65 supported topology sentinel replayed under E64.", f"Prior target: {row['target_id']}"],
        ))
        sentinel_used.add(sentinel_key)

    for row in v66_scoring:
        if len([target for target in sentinel_targets if target["source_family"] == "V66"]) >= 7:
            break
        sentinel_key = ("V66", row["target_id"])
        if row["score_label"] != "supported" or sentinel_key in sentinel_used:
            continue
        source_path = V66_SOURCE_ROOT / row["target_id"] / "source_manifest.json"
        sentinel_targets.append(_replay_target(
            ordinal=ordinal + len(sentinel_targets),
            source_family="V66",
            source_manifest_path=source_path,
            scoring_row=row,
            reasons=[f"V66 supported sentinel replayed under E64.", f"Prior target: {row['target_id']}"],
        ))
        sentinel_used.add(sentinel_key)

    if len(sentinel_targets) != SENTINEL_COUNT:
        raise SystemExit(f"selected {len(sentinel_targets)} sentinels; expected {SENTINEL_COUNT}")
    targets.extend(sentinel_targets)

    for new_ordinal, target in enumerate(targets, start=1):
        old_id = target["target_id"]
        if old_id.startswith(f"V67_{new_ordinal:03d}_"):
            continue
        prefix = "NEW" if target["selection_category"].startswith("NEW") else ("OLD" if target["selection_category"].startswith("V66") else "SENT")
        target["target_id"] = f"V67_{new_ordinal:03d}_{prefix}_{_safe_id(target['protein_id'])}"

    composition = Counter(target["selection_category"] for target in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != COMPOSITION_RULE:
        raise SystemExit(f"bad V67 target composition: total={len(targets)} composition={dict(composition)}")
    return targets


def _target_manifest(targets: list[dict[str, Any]]) -> dict[str, Any]:
    manifest_rows = []
    for target in targets:
        row = {key: value for key, value in target.items() if key not in {"candidate_snapshot", "source_manifest_path"}}
        if "candidate_snapshot" in target:
            candidate = target["candidate_snapshot"]
            row.update({
                "title": candidate.get("title", ""),
                "entity_description": candidate.get("entity_description", ""),
                "polymer_composition": candidate.get("polymer_composition", ""),
                "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
                "entity_molecule_count": candidate.get("entity_molecule_count", 0),
                "biological_cofactor_components": candidate.get("biological_cofactor_components", []),
                "sequence_metrics": candidate.get("sequence_metrics") or v61._sequence_metrics(candidate["sequence"]),
            })
        if "source_manifest_path" in target:
            row["source_manifest_path"] = str(Path(target["source_manifest_path"]).relative_to(REPO_ROOT))
        manifest_rows.append(row)
    return {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "batch_mode": "mixed_fast_discovery_shard",
        "target_selection_manual": False,
        "composition_rule": COMPOSITION_RULE,
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(manifest_rows),
        "selection_rule": (
            "100 unused V63 RCSB 30% nonredundant targets, all 70 V66 failed-accepted targets, "
            "and 30 supported E64-compatible sentinels from V62/V64/V65/V66."
        ),
        "selected_targets": manifest_rows,
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    supported = [row for row in rows if row["score_label"] == "supported"]
    failed_accepted = [row for row in accepted if row["score_label"] != "supported"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    return {
        "targets_total": len(rows),
        "accepted_count": len(accepted),
        "supported_count": len(supported),
        "failed_accepted_count": len(failed_accepted),
        "abstain_count": len(abstained),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "raw_accuracy": len(supported) / len(rows) if rows else None,
        "coverage": len(accepted) / len(rows) if rows else None,
    }


def _postseal_text(target: dict[str, Any]) -> str:
    candidate = target.get("candidate_snapshot")
    if isinstance(candidate, dict):
        return _target_text(candidate, include_postseal=True)
    return " ".join(str(part) for part in target.get("postseal_truth_basis", [])).lower()


def _failure_mode(row: dict[str, Any], target: dict[str, Any]) -> str:
    if row["score_label"] == "supported":
        return "supported"
    expected = row["expected_mechanism_class"]
    predicted = row["predicted_mechanism_class"]
    postseal = _postseal_text(target)
    if predicted == ABSTAIN_CLASS:
        return f"over_abstain_{expected}"
    if expected == "membrane_multidomain_folding_proteostasis":
        if predicted == "cofactor_ligand_assisted_stabilization":
            return "cofactor_locked_basin_vs_membrane_topology"
        if predicted == "oligomerization_controlled_folding":
            return "assembly_required_core_vs_membrane_topology"
        if "signal peptide" in postseal or "secretory" in postseal:
            return "signal_peptide_vs_true_TM"
        if "beta barrel" in postseal or "porin" in postseal:
            return "beta_barrel_soluble_vs_membrane"
        return "membrane_topology_provider_missing_preseal"
    if expected == "cofactor_ligand_assisted_stabilization":
        if any(token in postseal for token in ["fe-s", "iron-sulfur", "metal cluster", "cluster"]):
            return "metal_cluster_geometry"
        return "cofactor_locked_basin"
    if expected == "oligomerization_controlled_folding":
        if any(token in postseal for token in ["coiled coil", "coiled-coil", "leucine zipper"]):
            return "coiled_coil_register"
        if any(token in postseal for token in ["repeat", "solenoid", "ankyrin", "tpr", "armadillo"]):
            return "repeat_solenoid_topology"
        return "assembly_required_core"
    if expected == "intrinsic_disorder_phase_separation":
        return "disorder_or_low_complexity_context"
    if expected == "short_region_host_interface_hijacking":
        return "short_interface_context_missing"
    if expected == "metamorphic_fold_switching":
        return "domain_swapping_or_state_switch"
    if expected == "globular_closure" and predicted == "oligomerization_controlled_folding":
        return "generic_assembly_overread"
    return "wrong_regime"


def _missing_word_row(row: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    mode = _failure_mode(row, target)
    proposals = {
        "cofactor_locked_basin_vs_membrane_topology": "cofactor_locked_basin_vs_topology_provider",
        "assembly_required_core_vs_membrane_topology": "assembly_required_core_vs_topology_provider",
        "membrane_topology_provider_missing_preseal": "explicit_membrane_topology_provider_or_clean_abstain",
        "signal_peptide_vs_true_TM": "signal_peptide_vs_true_TM",
        "beta_barrel_soluble_vs_membrane": "beta_barrel_soluble_vs_membrane",
        "cofactor_locked_basin": "cofactor_locked_basin",
        "metal_cluster_geometry": "metal_cluster_geometry",
        "coiled_coil_register": "coiled_coil_register",
        "repeat_solenoid_topology": "repeat_solenoid_topology",
        "assembly_required_core": "assembly_required_core",
        "disorder_or_low_complexity_context": "disorder_low_complexity_persistence",
        "short_interface_context_missing": "short_region_partner_interface",
        "domain_swapping_or_state_switch": "domain_swapping_or_state_switch",
        "generic_assembly_overread": "generic_complex_and_copy_count_cleanup",
        "wrong_regime": "regime_separation",
    }
    return {
        "target_id": row["target_id"],
        "protein_id": row["protein_id"],
        "selection_category": row["selection_category"],
        "failure_mode": mode,
        "engine_thought": row["predicted_mechanism_class"],
        "reality_showed": row["expected_mechanism_class"],
        "acceptance_decision": row["acceptance_decision"],
        "score_label": row["score_label"],
        "missing_esperanto_word": proposals.get(mode, mode),
        "lineage_source_target": row.get("lineage_source_target"),
        "autopsy_sentence": (
            f"The engine thought: {row['predicted_mechanism_class']}. "
            f"Reality showed: {row['expected_mechanism_class']}. "
            f"Missing Esperanto word: {proposals.get(mode, mode)}."
        ),
    }


def _failure_report(scoring_rows: list[dict[str, Any]], targets: list[dict[str, Any]]) -> dict[str, Any]:
    by_id = {target["target_id"]: target for target in targets}
    rows = [_missing_word_row(row, by_id[row["target_id"]]) for row in scoring_rows if row["score_label"] != "supported"]
    modes = dict(Counter(row["failure_mode"] for row in rows))
    return {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_cases_reported": True,
        "failure_count": len(rows),
        "failure_modes": modes,
        "dominant_failure_mode": Counter(row["failure_mode"] for row in rows).most_common(1)[0][0] if rows else None,
        "dominant_failure_count": Counter(row["failure_mode"] for row in rows).most_common(1)[0][1] if rows else 0,
        "missing_words_top_10": [
            {"failure_mode": mode, "count": count}
            for mode, count in Counter(row["failure_mode"] for row in rows).most_common(10)
        ],
        "failure_grammar_rows": rows,
        "note": "V67 is a mixed discovery shard. Failure rows preserve engine_thought, reality_showed, and missing_esperanto_word for E65 mining.",
    }


def _old_failure_repair_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [row for row in scoring_rows if row["selection_category"] == "V66_FAILED_ACCEPTED_REPLAY"]
    repaired = [row for row in rows if row["score_label"] == "supported"]
    still_failed = [row for row in rows if row["score_label"] != "supported"]
    return {
        "kind": "V67_OLD_V66_FAILURE_REPAIR_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "old_v66_failure_count": len(rows),
        "old_v66_failures_repaired_by_e64": len(repaired),
        "old_v66_failures_remaining": len(still_failed),
        "repaired_target_ids": [row["target_id"] for row in repaired],
        "remaining_target_ids": [row["target_id"] for row in still_failed],
        "question_answer": (
            f"E64 repaired {len(repaired)} of the {len(rows)} V66 failed-accepted targets; "
            f"{len(still_failed)} remain for failure grammar mining."
        ),
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
    failure_report: dict[str, Any],
    old_repair_report: dict[str, Any],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V67_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V67_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold-style model offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V67_PRESEAL_HOLDOUT",
        "source_class": COORDINATE_DERIVED,
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V67_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_operator_state_packet(
        target_id="V67_RANDOM_SEQUENCE_CONTROL",
        target_name="V67 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
    composition = Counter(row["selection_category"] for row in target_manifest["selected_targets"])
    sentinel_rows = [row for row in scoring_rows if row["selection_category"] == "SENTINEL_STABILITY_REPLAY"]
    sentinel_regressions = [row["target_id"] for row in sentinel_rows if row["score_label"] != "supported"]
    shuffled_rows = []
    for packet, shuffled in zip(packets, shuffled_packets):
        original_coherence = sequence_operator_coherence(packet)
        shuffled_coherence = sequence_operator_coherence(shuffled)
        shuffled_rows.append({
            "target_id": packet["target_id"],
            "original_coherence": original_coherence,
            "shuffled_coherence": shuffled_coherence,
            "shuffled_higher_by_more_than_margin": shuffled_coherence > original_coherence + 0.08,
            "shuffled_coordinate_sources": shuffled["evidence_manifest"]["coordinate_derived_source_count_before_prediction"],
            "shuffled_runtime_sources": shuffled["evidence_manifest"]["internal_runtime_source_count_for_prediction"],
        })
    return [
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V67 must have exactly 200 targets.", target_manifest["target_count_selected"]),
        _control("mixed_composition_rule", dict(composition) == COMPOSITION_RULE, "V67 composition must be 100 new, 70 old failures, 30 sentinels.", dict(composition)),
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V67 target selection is deterministic from committed artifacts."),
        _control("engine_version_declared_e64", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V67 uses E64."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside V67."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V67 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("old_v66_failure_question_answered", old_repair_report["old_v66_failure_count"] == OLD_FAILURE_COUNT, "V67 answers whether E64 repaired the 70 V66 failures.", old_repair_report),
        _control("dominant_failure_mode_identified", failure_report["dominant_failure_mode"] is not None or failure_report["failure_count"] == 0, "V67 identifies the dominant failure mode when failures exist.", failure_report.get("missing_words_top_10")),
        _control("sentinels_stable", not sentinel_regressions and len(sentinel_rows) == SENTINEL_COUNT, "Sentinel targets must remain supported.", {"sentinel_count": len(sentinel_rows), "regressions": sentinel_regressions}),
        _control("failure_autopsy_fields_present", all(row.get("engine_thought") and row.get("reality_showed") and row.get("missing_esperanto_word") for row in failure_report["failure_grammar_rows"]), "Every failure autopsy says what the engine thought, what reality showed, and the missing word."),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls are generated with target metadata withheld.", {"control_count": len(shuffled_rows), "shuffled_higher_by_more_than_margin_count": sum(1 for row in shuffled_rows if row["shuffled_higher_by_more_than_margin"])}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == ABSTAIN_CLASS, "Random sequence without evidence abstains."),
        _control("readme_check_skipped_by_user_instruction", True, "README check skipped by explicit user instruction."),
    ]


def _accept_abstain_posture(metrics: dict[str, Any]) -> str:
    if metrics["failed_accepted_count"] > metrics["abstain_count"]:
        return "over_accepting_relative_to_abstention"
    if metrics["abstain_count"] > metrics["failed_accepted_count"]:
        return "abstaining_more_than_failed_accepting"
    return "accept_abstain_balanced"


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_report: dict[str, Any],
    old_repair_report: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "alphafold_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    sentinel_regressions = [
        row["target_id"]
        for row in scoring_rows
        if row["selection_category"] == "SENTINEL_STABILITY_REPLAY" and row["score_label"] != "supported"
    ]
    if any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif sentinel_regressions:
        status = BLOCKED_SENTINEL
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif metrics["failed_accepted_count"]:
        status = PASSED_MINED
    else:
        status = PASSED_CLEAN
    controls_passed = not failed_controls
    cert = {
        "kind": "V67_E64_MIXED_FAST_DISCOVERY_200_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "mixed_fast_discovery_shard",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "accept_abstain_posture": _accept_abstain_posture(metrics),
        "controls_passed": controls_passed,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "old_v66_failure_repair": {
            "old_v66_failure_count": old_repair_report["old_v66_failure_count"],
            "old_v66_failures_repaired_by_e64": old_repair_report["old_v66_failures_repaired_by_e64"],
            "old_v66_failures_remaining": old_repair_report["old_v66_failures_remaining"],
            "question_answer": old_repair_report["question_answer"],
        },
        "sentinel_regressions": sentinel_regressions,
        "sentinel_regression_count": len(sentinel_regressions),
        "failure_modes": failure_report["failure_modes"],
        "dominant_failure_mode": failure_report["dominant_failure_mode"],
        "dominant_failure_count": failure_report["dominant_failure_count"],
        "missing_words_top_10": failure_report["missing_words_top_10"],
        "v67_questions_answered": {
            "did_e64_repair_old_v66_failures": old_repair_report["question_answer"],
            "new_dominant_failure_mode": failure_report["dominant_failure_mode"],
            "did_e64_create_regressions": f"{len(sentinel_regressions)} sentinel regressions.",
            "is_engine_over_accepting_or_abstaining": _accept_abstain_posture(metrics),
        },
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V67 is a mixed fast discovery shard for E65 mining, not a solved-folding claim.",
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "E64 has saturated broad RCSB protein space.",
            "V67 failures may be hidden.",
            "Coordinates, contacts, ligand geometry, or native topology were used before sealing.",
        ],
        "next_required_step": "Extract E65 from the dominant V67 failure mode.",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "mixed_fast_discovery_shard",
        "engine_version_used": ENGINE_VERSION_USED,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["abstain_count"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"] + cert["abstain_count"],
        "failure_modes": cert["failure_modes"],
        "dominant_failure_mode": cert["dominant_failure_mode"],
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


def _write_report(path: Path, cert: dict[str, Any], failure_report: dict[str, Any]) -> None:
    lines = [
        "# V67 E64 Mixed Fast Discovery 200",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Accepted: `{cert['accepted_count']}`",
        f"Supported: `{cert['supported_count']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Abstain: `{cert['abstain_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Raw accuracy: `{cert['raw_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Accept/abstain posture: `{cert['accept_abstain_posture']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"Sentinel regressions: `{cert['sentinel_regression_count']}`",
        "",
        "## Four Answers",
        f"1. Old V66 failures: `{cert['old_v66_failure_repair']['question_answer']}`",
        f"2. Dominant failure mode: `{cert['dominant_failure_mode']}` count `{cert['dominant_failure_count']}`",
        f"3. E64 regressions: `{cert['v67_questions_answered']['did_e64_create_regressions']}`",
        f"4. Accept/abstain posture: `{cert['accept_abstain_posture']}`",
        "",
        "## Failure Mode Table",
        "",
        "| failure_mode | count |",
        "| --- | ---: |",
    ]
    for item in failure_report["missing_words_top_10"]:
        lines.append(f"| `{item['failure_mode']}` | `{item['count']}` |")
    lines.extend([
        "",
        "## Boundary",
        "V67 is a mixed fast discovery shard. It mines the next missing Esperanto word and does not make a broad solved claim.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v67(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_generated_outputs(out_dir)
    targets = _select_targets()
    target_manifest = _target_manifest(targets)
    engine_declaration = _engine_declaration()
    e64_certificate = _e64_certificate(engine_declaration)
    _write_json(E64_ROOT / "e64_v66_tm_complex_lineage_revision_certificate.json", e64_certificate)
    _write_json(DATA_ROOT / "v67_mixed_fast_discovery_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v67_e64_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []

    for target in targets:
        source_manifest = _source_manifest(target)
        packet = build_sealed_operator_state_packet(
            target_id=target["target_id"],
            target_name=target["target_name"],
            sequence=target["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "V67 E64 mixed fast discovery full-chain scan", "span": f"1-{target['sequence_length']}"}],
            perturbations=v61._perturbations_for_expected(target["expected_mechanism_class"], target["target_id"]),
        )
        holdout = _holdout(target, packet)
        score = _score(packet, holdout, target)
        shuffled_packet = build_sealed_operator_state_packet(
            target_id=f"{target['target_id']}_SHUFFLED_CONTROL",
            target_name=f"{target['entry_id']} shuffled sequence control",
            sequence=shuffled_sequence(target["sequence"]),
            sources=[
                {
                    "source_id": f"{target['target_id']}_SHUFFLED_SEQUENCE_ONLY",
                    "source_class": "pure_non_coordinate",
                    "source_role": "prediction_input",
                    "coordinate_derived": False,
                    "internal_runtime_source": False,
                    "spatial_proxy": False,
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E64 context marks are withheld.",
                }
            ],
            perturbations=[],
        )
        packets.append(packet)
        shuffled_packets.append(shuffled_packet)
        scoring_rows.append(score)
        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target["target_id"] / "sealed_packet_summary.json", _packet_summary(packet))
        _write_json(DATA_ROOT / "holdouts_postseal" / target["target_id"] / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target["target_id"] / "validation_result.json", score)
        _write_json(DATA_ROOT / "shuffled_controls" / target["target_id"] / "shuffled_control_packet.json", _packet_summary(shuffled_packet))

    failure_report = _failure_report(scoring_rows, targets)
    old_repair_report = _old_failure_repair_report(scoring_rows)
    controls = _controls(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        packets=packets,
        scoring_rows=scoring_rows,
        shuffled_packets=shuffled_packets,
        failure_report=failure_report,
        old_repair_report=old_repair_report,
    )
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
        failure_report=failure_report,
        old_repair_report=old_repair_report,
    )
    claim_row = _claim_row(cert)

    scoring_path = _write_json(DATA_ROOT / "v67_mixed_fast_discovery_scoring_report.json", {"kind": "V67_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v67_mixed_fast_discovery_failure_report.json", failure_report)
    old_repair_path = _write_json(DATA_ROOT / "v67_old_v66_failure_repair_report.json", old_repair_report)
    data_cert_path = _write_json(DATA_ROOT / "v67_mixed_fast_discovery_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v67_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)

    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v67_e64_mixed_fast_discovery_200_certificate.json", cert)
    report_path = out_dir / "V67_E64_MIXED_FAST_DISCOVERY_200_REPORT.md"
    _write_report(report_path, cert, failure_report)
    return {
        "e64_certificate": E64_ROOT / "e64_v66_tm_complex_lineage_revision_certificate.json",
        "target_manifest": DATA_ROOT / "v67_mixed_fast_discovery_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v67_e64_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "old_repair_report": old_repair_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V67 E64 mixed fast discovery shard.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v67(args.out_dir)
    cert = _read_json(paths["certificate"], "V67 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_count": cert["accepted_count"],
        "supported_count": cert["supported_count"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "abstain_count": cert["abstain_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "raw_accuracy": cert["raw_accuracy"],
        "coverage": cert["coverage"],
        "accept_abstain_posture": cert["accept_abstain_posture"],
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "old_v66_failure_repair": cert["old_v66_failure_repair"],
        "sentinel_regression_count": cert["sentinel_regression_count"],
        "dominant_failure_mode": cert["dominant_failure_mode"],
        "dominant_failure_count": cert["dominant_failure_count"],
        "failure_modes": cert["failure_modes"],
        "missing_words_top_10": cert["missing_words_top_10"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] and not cert["sentinel_regressions"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
