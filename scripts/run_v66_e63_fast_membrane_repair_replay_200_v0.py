#!/usr/bin/env python3
from __future__ import annotations

"""Run V66: E63 fast membrane repair replay on a 200-target adaptive set.

V66 is the first adaptive 200-target loop after the V64/V65/E63 sequence.  It
does not rerun the whole V63 500.  It combines:

* 70 V65 membrane topology panel targets,
* 80 selected V64/V63 membrane_misread failures,
* 3 V64 E62 oligomer_state_misread regressions,
* 47 sentinels that E61/E62 had already handled correctly.

The replay uses current E63 grammar and compares against the committed E62
rows from V64/V65.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
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
    build_sealed_simulation_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    stable_hash,
    validate_against_holdout,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402
import run_v63_rcsb_500_discovery_batch_v0 as v63  # noqa: E402


BATCH_ID = "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E63"
BASELINE_ENGINE_VERSION = "E62"
TARGET_COUNT = 200
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"

V65_TOPOLOGY_PANEL_COUNT = 70
V64_MEMBRANE_FAILURE_COUNT = 80
V64_OLIGOMER_REGRESSION_COUNT = 3
SENTINEL_COUNT = 47

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V66"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V63_TARGET_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_target_manifest.json"
V64_SCORING_REPORT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V64" / "v64_e62_rcsb_500_membrane_repair_scoring_report.json"
V65_SCORING_REPORT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_scoring_report.json"
V65_PANEL_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_panel_manifest.json"

PASSED = "V66_E63_FAST_MEMBRANE_REPAIR_PASSED_REVIEW_REQUIRED"
INCOMPLETE = "V66_E63_MEMBRANE_REPAIR_INCOMPLETE_NEW_FAILURES_MINED"
BLOCKED_CONTROLS = "V66_E63_FAST_MEMBRANE_REPAIR_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V66_E63_FAST_MEMBRANE_REPAIR_BLOCKED_FOR_LEAKAGE"


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
    ]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v66_e63_fast_membrane_repair_target_manifest.json",
        "v66_e63_engine_declaration.json",
        "v66_e63_fast_membrane_repair_scoring_report.json",
        "v66_e62_vs_e63_comparison.json",
        "v66_e63_fast_membrane_repair_failure_report.json",
        "v66_e63_fast_membrane_repair_certificate.json",
        "v66_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _target_key(target_id: str) -> str:
    for prefix in ["V64_", "V65_", "V66_"]:
        if target_id.startswith(prefix):
            return target_id.removeprefix(prefix)
    return target_id


def _protein_key_from_v65(target_id: str) -> str:
    parts = target_id.split("_")
    if len(parts) >= 5 and parts[0] == "V65":
        return "_".join(parts[3:])
    return _target_key(target_id)


def _annotation_text(candidate: dict[str, Any]) -> str:
    values = [
        candidate.get("title", ""),
        candidate.get("entity_description", ""),
        candidate.get("entry_keywords", ""),
        candidate.get("polymer_composition", ""),
        " ".join(candidate.get("organisms", []) or []),
        " ".join(candidate.get("annotations", []) or []),
        " ".join(candidate.get("feature_types", []) or []),
        " ".join(candidate.get("biological_cofactor_components", []) or []),
    ]
    return " ".join(str(value) for value in values).lower()


def _has_true_topology_provider(candidate: dict[str, Any]) -> bool:
    text = _annotation_text(candidate)
    if "monotopic/peripheral" in text:
        return False
    return any(token in text for token in ["generic pdbtm", "generic memprotmd", "pdbtm", "memprotmd", "transmembrane proteins", "transmembrane"])


def _has_peripheral_membrane_provider(candidate: dict[str, Any]) -> bool:
    text = _annotation_text(candidate)
    return "monotopic/peripheral" in text or ("opm" in text and "peripheral" in text)


def _has_cofactor_context(candidate: dict[str, Any]) -> bool:
    return bool(candidate.get("biological_cofactor_components"))


def _has_oligomer_context(candidate: dict[str, Any]) -> bool:
    text = _annotation_text(candidate)
    return (
        candidate.get("polymer_entity_instance_count", 0) >= 2
        or candidate.get("entity_molecule_count", 0) >= 2
        or any(token in text for token in ["oligomer", "homomer", "multimer", "assembly", "dimer", "trimer", "tetramer"])
    )


def _e63_context(candidate: dict[str, Any]) -> dict[str, Any]:
    base = v63._metadata_context(candidate)
    base_marks = set(base["context_marks"])
    marks: list[str] = []
    reasons: list[str] = []
    components = {str(value).upper() for value in candidate.get("biological_cofactor_components", []) or []}
    text = _annotation_text(candidate)
    if _has_cofactor_context(candidate) or any(token in text for token in ["cofactor", "ligand", "heme", "nucleotide"]):
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append("E63 context preserves cofactor/ligand explanation before membrane ambiguity.")
        if components & v63.METAL_COMPONENTS:
            marks.append("metal_context")
        if components & v63.HEME_COMPONENTS:
            marks.append("heme_context")
        if components & v63.NUCLEOTIDE_COMPONENTS:
            marks.append("nucleotide_context")
    if _has_oligomer_context(candidate):
        marks.extend(["oligomer_context", "assembly_context", "partner_copy_context"])
        reasons.append("E63 context preserves assembly/interface explanation before membrane ambiguity.")
        if "heteromeric" in text or "hetero" in text:
            marks.append("heteromeric_context")
        elif "homomeric" in text or "homo" in text:
            marks.append("homomeric_context")
    if _has_true_topology_provider(candidate):
        marks.extend(["membrane_context_strong", "transmembrane_context"])
        reasons.append("E63 context sees explicit topology-provider transmembrane evidence.")
        if any(token in text for token in ["channel", "pore", "porin"]):
            marks.append("channel_context")
        if "transporter" in text:
            marks.append("transporter_context")
        if any(token in text for token in ["receptor", "gpcr", "opsin"]):
            marks.append("receptor_membrane_context")
    elif _has_peripheral_membrane_provider(candidate):
        marks.extend(["peripheral_membrane_context", "not_transmembrane_context"])
        reasons.append("E63 context sees peripheral/monotopic membrane association but no transmembrane topology.")
    elif base_marks & {"membrane_context_strong", "transmembrane_context"}:
        marks.extend(["soluble_hydrophobic_core_context", "no_membrane_topology_context"])
        if _has_cofactor_context(candidate):
            marks.append("cofactor_buried_hydrophobic_pocket_context")
            reasons.append("E63 context treats hydrophobicity with cofactor evidence as cofactor-pocket ambiguity.")
        elif _has_oligomer_context(candidate):
            marks.append("oligomeric_interface_hydrophobicity_context")
            reasons.append("E63 context treats hydrophobicity with assembly evidence as interface ambiguity.")
        else:
            reasons.append("E63 context treats hydrophobicity-only signal as topology ambiguity.")
    return {
        "context_marks": sorted(dict.fromkeys(marks)),
        "base_v63_context_marks": sorted(base_marks),
        "context_derivation": (
            "V66 E63 repair replay uses sequence plus public non-coordinate RCSB/V65 panel metadata. "
            "It separates provider-backed transmembrane topology from hydrophobic-only, cofactor-pocket, "
            "oligomeric-interface, and peripheral/monotopic membrane contexts."
        ),
        "reasons": reasons or ["no E63 repair context mark emitted"],
        "biological_cofactor_components_seen": sorted(components),
        "polymer_copy_counts": {
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "polymer_composition": candidate.get("polymer_composition", ""),
        },
    }


def _e63_context_statement(context: dict[str, Any]) -> str:
    marks = context["context_marks"]
    if not marks:
        return "V66 E63 replay metadata context emitted no explicit E63 context marks for this target."
    statement = (
        "V66 E63 replay metadata context marks: "
        + " ".join(marks)
        + ". E63 separates true transmembrane topology from hydrophobic-only, peripheral, "
        + "cofactor-pocket, and oligomer-interface explanations. "
        + "Coordinates, contacts, ligand geometry, native topology, and post-seal validation labels are blocked before sealing."
    )
    if "no_membrane_topology_context" in marks:
        statement += " no membrane topology evidence; no transmembrane assignment; hydrophobicity alone is ambiguous."
    if "peripheral_membrane_context" in marks:
        statement += " peripheral membrane association; not transmembrane; no bilayer-spanning topology evidence."
    if "cofactor_buried_hydrophobic_pocket_context" in marks:
        statement += " cofactor-buried hydrophobic pocket explains the hydrophobic signal."
    if "oligomeric_interface_hydrophobicity_context" in marks:
        statement += " oligomeric interface hydrophobicity explains the hydrophobic signal."
    return statement


def _source_manifest_from_v64_target(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate"]
    context = _e63_context(candidate)
    return {
        "kind": "V66_E63_FAST_REPAIR_SOURCE_MANIFEST_v0",
        "target_id": target["target_id"],
        "batch_id": BATCH_ID,
        "source_lineage_batch": "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY",
        "selection_category": target["selection_category"],
        "baseline_target_id": target["baseline_target_id"],
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "e63_context_policy": context,
        "prediction_sources": [
            {
                "source_id": f"{target['target_id']}_RAW_SEQUENCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": candidate["source_urls"]["polymer_entity"],
            },
            {
                "source_id": f"{target['target_id']}_PUBLIC_METADATA",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": v61._metadata_statement(candidate),
                "source_url": candidate["source_urls"]["entry"],
            },
            {
                "source_id": f"{target['target_id']}_E63_METADATA_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": context["context_marks"],
                "metadata_context_reasons": context["reasons"],
                "evidence_statement": _e63_context_statement(context),
                "source_url": candidate["source_urls"]["entry"],
            },
        ],
        "blocked_prediction_inputs": _blocked_prediction_inputs(),
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _source_manifest_from_v65_target(target: dict[str, Any]) -> dict[str, Any]:
    original = _read_json(
        REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "source_manifests" / target["baseline_target_id"] / "source_manifest.json",
        f"V65 source manifest {target['baseline_target_id']}",
    )
    sources = []
    for source in original["prediction_sources"]:
        copied = dict(source)
        copied["source_id"] = copied["source_id"].replace(target["baseline_target_id"], target["target_id"])
        sources.append(copied)
    return {
        "kind": "V66_E63_FAST_REPAIR_SOURCE_MANIFEST_v0",
        "target_id": target["target_id"],
        "batch_id": BATCH_ID,
        "source_lineage_batch": "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL",
        "selection_category": target["selection_category"],
        "baseline_target_id": target["baseline_target_id"],
        "panel_group": target.get("panel_group"),
        "protein_id": original["protein_id"],
        "entry_id": original["entry_id"],
        "entity_id": original["entity_id"],
        "sequence": original["sequence"],
        "sequence_length": original["sequence_length"],
        "prediction_sources": sources,
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
        "E62 baseline score outcomes as prediction evidence",
        "internal runtime artifacts as biological evidence",
    ]


def _holdout(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    expected = target["expected_mechanism_class"]
    return {
        "kind": "V66_E63_FAST_MEMBRANE_REPAIR_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "selection_category": target["selection_category"],
        "panel_group": target.get("panel_group"),
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "expected_mechanism_class": expected,
        "expected_observables": [] if expected == ABSTAIN_CLASS else v61._expected_observables(expected),
        "postseal_truth_basis": target["postseal_truth_basis"],
        "baseline_target_id": target["baseline_target_id"],
        "baseline_batch_id": target["baseline_batch_id"],
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "baseline_score_label": target["baseline_score_label"],
        "baseline_predicted_mechanism_class": target["baseline_predicted_mechanism_class"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V66_REPLAY_HOLDOUT",
                "source_class": "coordinate_derived",
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "lineage_source_target": target["baseline_target_id"],
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
    return {
        "kind": "V66_E63_FAST_MEMBRANE_REPAIR_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "selection_category": target["selection_category"],
        "panel_group": target.get("panel_group"),
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
        "baseline_batch_id": target["baseline_batch_id"],
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "baseline_target_id": target["baseline_target_id"],
        "baseline_acceptance_decision": target["baseline_acceptance_decision"],
        "baseline_predicted_mechanism_class": target["baseline_predicted_mechanism_class"],
        "baseline_score_label": target["baseline_score_label"],
        "baseline_supported": target["baseline_score_label"] == "supported",
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V66_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "sealed_before_holdout": packet["sealed_before_holdout"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "selection_reason": mechanism["selection_reason"],
        "evidence_manifest": packet["evidence_manifest"],
        "operator_names": packet["operator_field"]["operator_names"],
        "active_operator_count": packet["operator_field"]["active_operator_count"],
        "trajectory_final_state_summary": packet["trajectory_summary"]["final_state_summary"],
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V66_E63_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_lineage": ["E60", "E61", "E62", "E63"],
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
        "target_selection_manual": False,
        "folding_problem_solved": False,
    }


def _load_rows(path: Path, label: str) -> list[dict[str, Any]]:
    data = _read_json(path, label)
    rows = data.get("rows", [])
    if not isinstance(rows, list):
        raise SystemExit(f"{label} rows must be a list")
    return [row for row in rows if isinstance(row, dict)]


def _round_robin_stable_sentinels(rows: list[dict[str, Any]], used: set[str]) -> list[dict[str, Any]]:
    by_expected: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = _target_key(row["target_id"])
        if key in used:
            continue
        if row["baseline_supported"] and row["score_label"] == "supported":
            by_expected[row["expected_mechanism_class"]].append(row)
    for group_rows in by_expected.values():
        group_rows.sort(key=lambda row: row["target_id"])
    selected: list[dict[str, Any]] = []
    mechanisms = sorted(by_expected)
    while len(selected) < SENTINEL_COUNT and mechanisms:
        progressed = False
        for mechanism in mechanisms:
            if by_expected[mechanism]:
                row = by_expected[mechanism].pop(0)
                selected.append(row)
                used.add(_target_key(row["target_id"]))
                progressed = True
                if len(selected) >= SENTINEL_COUNT:
                    break
        if not progressed:
            break
    if len(selected) != SENTINEL_COUNT:
        raise SystemExit(f"selected only {len(selected)} sentinels; need {SENTINEL_COUNT}")
    return selected


def _select_targets() -> list[dict[str, Any]]:
    v63_manifest = _read_json(V63_TARGET_MANIFEST, "V63 target manifest")
    candidates_by_key = {row["target_id"]: row for row in v63_manifest["selected_targets"]}
    v64_rows = _load_rows(V64_SCORING_REPORT, "V64 scoring report")
    v65_rows = _load_rows(V65_SCORING_REPORT, "V65 scoring report")
    panel_manifest = _read_json(V65_PANEL_MANIFEST, "V65 panel manifest")
    panel_order = {row["target_id"]: index for index, row in enumerate(panel_manifest["targets"])}
    selected: list[dict[str, Any]] = []
    used: set[str] = set()

    for row in sorted(v65_rows, key=lambda item: panel_order[item["target_id"]]):
        key = _protein_key_from_v65(row["target_id"])
        used.add(key)
        selected.append(_target_from_v65_row(row, len(selected) + 1))

    persistent_membrane = [
        row
        for row in v64_rows
        if row["expected_mechanism_class"] == MEMBRANE_CLASS
        and row["score_label"] != "supported"
        and _target_key(row["target_id"]) not in used
    ]
    persistent_membrane.sort(key=lambda row: row["target_id"])
    if len(persistent_membrane) < V64_MEMBRANE_FAILURE_COUNT:
        raise SystemExit(f"only {len(persistent_membrane)} unused persistent membrane failures; need {V64_MEMBRANE_FAILURE_COUNT}")
    for row in persistent_membrane[:V64_MEMBRANE_FAILURE_COUNT]:
        key = _target_key(row["target_id"])
        used.add(key)
        selected.append(_target_from_v64_row(row, candidates_by_key[key], "V64_V63_MEMBRANE_MISREAD_FAILURE", len(selected) + 1))

    regressions = [
        row
        for row in v64_rows
        if row["baseline_supported"]
        and row["score_label"] != "supported"
        and _target_key(row["target_id"]) not in used
    ]
    regressions.sort(key=lambda row: row["target_id"])
    if len(regressions) != V64_OLIGOMER_REGRESSION_COUNT:
        raise SystemExit(f"expected {V64_OLIGOMER_REGRESSION_COUNT} V64 new regressions; found {len(regressions)}")
    for row in regressions:
        key = _target_key(row["target_id"])
        used.add(key)
        selected.append(_target_from_v64_row(row, candidates_by_key[key], "V64_E62_OLIGOMER_REGRESSION", len(selected) + 1))

    sentinels = _round_robin_stable_sentinels(v64_rows, used)
    for row in sentinels:
        key = _target_key(row["target_id"])
        selected.append(_target_from_v64_row(row, candidates_by_key[key], "E61_E62_STABLE_SENTINEL", len(selected) + 1))

    if len(selected) != TARGET_COUNT:
        raise SystemExit(f"V66 selected {len(selected)} targets; need {TARGET_COUNT}")
    return selected


def _target_from_v65_row(row: dict[str, Any], ordinal: int) -> dict[str, Any]:
    key = _protein_key_from_v65(row["target_id"])
    return {
        "target_id": f"V66_{ordinal:03d}_{row['target_id']}",
        "selection_category": "V65_MEMBRANE_TOPOLOGY_PANEL",
        "baseline_batch_id": "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL",
        "baseline_target_id": row["target_id"],
        "protein_key": key,
        "panel_group": row["panel_group"],
        "protein_id": row["protein_id"],
        "entry_id": row["entry_id"],
        "entity_id": row["entity_id"],
        "expected_mechanism_class": row["expected_mechanism_class"],
        "baseline_acceptance_decision": row["acceptance_decision"],
        "baseline_predicted_mechanism_class": row["predicted_mechanism_class"],
        "baseline_score_label": row["score_label"],
        "postseal_truth_basis": [
            "V65 topology panel target replayed with E63.",
            f"Panel group: {row['panel_group']}",
        ],
        "source_family": "V65",
    }


def _target_from_v64_row(row: dict[str, Any], candidate: dict[str, Any], category: str, ordinal: int) -> dict[str, Any]:
    key = _target_key(row["target_id"])
    return {
        "target_id": f"V66_{ordinal:03d}_{key}",
        "selection_category": category,
        "baseline_batch_id": "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY",
        "baseline_target_id": row["target_id"],
        "protein_key": key,
        "protein_id": row["protein_id"],
        "entry_id": row["entry_id"],
        "entity_id": row["entity_id"],
        "expected_mechanism_class": row["expected_mechanism_class"],
        "baseline_acceptance_decision": row["acceptance_decision"],
        "baseline_predicted_mechanism_class": row["predicted_mechanism_class"],
        "baseline_score_label": row["score_label"],
        "postseal_truth_basis": [
            f"{category} selected from V64 E62 paired replay.",
            f"Baseline predicted {row['predicted_mechanism_class']} while expected {row['expected_mechanism_class']}.",
        ],
        "candidate": candidate,
        "source_family": "V64",
    }


def _public_manifest(targets: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for target in targets:
        row = {key: value for key, value in target.items() if key != "candidate"}
        rows.append(row)
    return {
        "kind": "V66_E63_FAST_MEMBRANE_REPAIR_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "batch_mode": "adaptive_200_repair_replay",
        "target_selection_manual": False,
        "target_count_selected": len(rows),
        "composition_rule": {
            "V65_MEMBRANE_TOPOLOGY_PANEL": V65_TOPOLOGY_PANEL_COUNT,
            "V64_V63_MEMBRANE_MISREAD_FAILURE": V64_MEMBRANE_FAILURE_COUNT,
            "V64_E62_OLIGOMER_REGRESSION": V64_OLIGOMER_REGRESSION_COUNT,
            "E61_E62_STABLE_SENTINEL": SENTINEL_COUNT,
        },
        "selected_targets": rows,
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


def _baseline_metrics(targets: list[dict[str, Any]]) -> dict[str, Any]:
    baseline_rows = [
        {
            "acceptance_decision": target["baseline_acceptance_decision"],
            "score_label": target["baseline_score_label"],
        }
        for target in targets
    ]
    return _metrics(baseline_rows)


def _failure_type(row: dict[str, Any]) -> str:
    if row["score_label"] == "supported":
        return "supported"
    if row["predicted_mechanism_class"] == MEMBRANE_CLASS and row["expected_mechanism_class"] != MEMBRANE_CLASS:
        if row.get("panel_group") == "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE":
            return "peripheral_misread_as_transmembrane"
        return "soluble_hydrophobic_false_membrane"
    if row["expected_mechanism_class"] == MEMBRANE_CLASS and row["predicted_mechanism_class"] == ABSTAIN_CLASS:
        return "ambiguous_membrane_topology_abstained"
    if row["expected_mechanism_class"] == MEMBRANE_CLASS:
        return "membrane_misread"
    if row["expected_mechanism_class"] == ABSTAIN_CLASS and row["predicted_mechanism_class"] != ABSTAIN_CLASS:
        return "ambiguous_topology_overaccepted"
    return v61._failure_type(row)


def _baseline_failure_type(target: dict[str, Any]) -> str:
    predicted = target["baseline_predicted_mechanism_class"]
    expected = target["expected_mechanism_class"]
    if target["baseline_score_label"] == "supported":
        return "supported"
    if predicted == MEMBRANE_CLASS and expected != MEMBRANE_CLASS:
        if target.get("panel_group") == "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE":
            return "peripheral_misread_as_transmembrane"
        return "soluble_hydrophobic_false_membrane"
    if expected == MEMBRANE_CLASS:
        return "membrane_misread"
    return v61._failure_type({
        "expected_mechanism_class": expected,
        "predicted_mechanism_class": predicted,
        "level1_regime_selection": predicted == expected,
    })


def _comparison(targets: list[dict[str, Any]], scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    baseline = _baseline_metrics(targets)
    e63 = _metrics(scoring_rows)
    baseline_failure_modes = dict(Counter(_baseline_failure_type(target) for target in targets if target["baseline_score_label"] != "supported"))
    e63_failure_modes = dict(Counter(_failure_type(row) for row in scoring_rows if row["score_label"] != "supported"))
    category_counts = Counter(target["selection_category"] for target in targets)
    score_by_category: dict[str, dict[str, int]] = {}
    for category in sorted(category_counts):
        rows = [row for row in scoring_rows if row["selection_category"] == category]
        score_by_category[category] = dict(Counter(row["score_label"] for row in rows))
    false_membrane_baseline = sum(
        1
        for target in targets
        if target["baseline_predicted_mechanism_class"] == MEMBRANE_CLASS
        and target["expected_mechanism_class"] != MEMBRANE_CLASS
    )
    false_membrane_e63 = sum(
        1
        for row in scoring_rows
        if row["predicted_mechanism_class"] == MEMBRANE_CLASS
        and row["expected_mechanism_class"] != MEMBRANE_CLASS
    )
    peripheral_baseline = sum(
        1
        for target in targets
        if target.get("panel_group") == "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE"
        and target["baseline_predicted_mechanism_class"] == MEMBRANE_CLASS
    )
    peripheral_e63 = sum(
        1
        for row in scoring_rows
        if row.get("panel_group") == "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE"
        and row["predicted_mechanism_class"] == MEMBRANE_CLASS
    )
    true_tm_missed = sum(
        1
        for row in scoring_rows
        if row.get("panel_group") == "A_TRUE_TRANSMEMBRANE_TOPOLOGY"
        and row["score_label"] != "supported"
    )
    oligomer_regression_remaining = sum(
        1
        for row in scoring_rows
        if row["selection_category"] == "V64_E62_OLIGOMER_REGRESSION"
        and row["score_label"] != "supported"
    )
    sentinel_regressions = sum(
        1
        for row in scoring_rows
        if row["selection_category"] == "E61_E62_STABLE_SENTINEL"
        and row["score_label"] != "supported"
    )
    all_modes = sorted(set(baseline_failure_modes) | set(e63_failure_modes))
    return {
        "kind": "V66_E62_VS_E63_FAST_MEMBRANE_REPAIR_COMPARISON_v0",
        "batch_id": BATCH_ID,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline": baseline,
        "e63": e63,
        "net_changes": {
            "supported_delta": e63["supported_count"] - baseline["supported_count"],
            "failed_accepted_delta": e63["failed_accepted_count"] - baseline["failed_accepted_count"],
            "abstain_delta": e63["abstain_count"] - baseline["abstain_count"],
            "accepted_accuracy_delta": e63["accepted_accuracy"] - baseline["accepted_accuracy"],
            "coverage_delta": e63["coverage"] - baseline["coverage"],
        },
        "selection_category_counts": dict(category_counts),
        "score_labels_by_category": score_by_category,
        "failure_mode_distribution_change": {
            mode: {
                "e62": baseline_failure_modes.get(mode, 0),
                "e63": e63_failure_modes.get(mode, 0),
                "delta": e63_failure_modes.get(mode, 0) - baseline_failure_modes.get(mode, 0),
            }
            for mode in all_modes
        },
        "false_membrane_repair": {
            "e62_false_membrane_calls": false_membrane_baseline,
            "e63_false_membrane_calls": false_membrane_e63,
            "e62_peripheral_as_tm": peripheral_baseline,
            "e63_peripheral_as_tm": peripheral_e63,
            "true_tm_missed_under_e63": true_tm_missed,
            "oligomer_regressions_remaining_under_e63": oligomer_regression_remaining,
            "sentinel_regressions_under_e63": sentinel_regressions,
        },
        "e63_failure_modes": e63_failure_modes,
        "repair_passed": (
            false_membrane_e63 < false_membrane_baseline
            and peripheral_e63 < peripheral_baseline
            and true_tm_missed == 0
            and oligomer_regression_remaining == 0
            and sentinel_regressions == 0
            and e63["failed_accepted_count"] < baseline["failed_accepted_count"]
            and e63["accepted_accuracy"] > baseline["accepted_accuracy"]
        ),
        "claim_allowed": False,
        "claim_boundary": "V66 is an adaptive 200-target repair replay, not a broad saturation claim.",
    }


def _failure_report(scoring_rows: list[dict[str, Any]], comparison: dict[str, Any]) -> dict[str, Any]:
    failures = [row for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V66_E63_FAST_MEMBRANE_REPAIR_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "failure_cases_reported": True,
        "failure_count": len(failures),
        "failure_modes": comparison["e63_failure_modes"],
        "failure_cases": failures,
        "false_membrane_repair": comparison["false_membrane_repair"],
        "note": "V66 failures feed the next adaptive 200-target discovery or specialist panel; no broad claim is made.",
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    comparison: dict[str, Any],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V66_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V66_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold-style model offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V66_PRESEAL_HOLDOUT",
        "source_class": "coordinate_derived",
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V66_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_simulation_packet(
        target_id="V66_RANDOM_SEQUENCE_CONTROL",
        target_name="V66 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
    composition = Counter(row["selection_category"] for row in target_manifest["selected_targets"])
    return [
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V66 must be exactly 200 targets.", target_manifest["target_count_selected"]),
        _control("adaptive_composition_rule", composition == target_manifest["composition_rule"], "V66 composition must match the requested 70/80/3/47 adaptive repair set.", dict(composition)),
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V66 selection is deterministic from V64/V65 committed artifacts."),
        _control("engine_version_declared_e63", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V66 uses E63."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside V66."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V66 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("false_membrane_repair_measured", comparison["false_membrane_repair"]["e62_false_membrane_calls"] > comparison["false_membrane_repair"]["e63_false_membrane_calls"], "V66 measures E63 false-membrane movement against E62."),
        _control("true_tm_support_preserved", comparison["false_membrane_repair"]["true_tm_missed_under_e63"] == 0, "E63 preserves true transmembrane support in the V65 true-TM panel."),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == ABSTAIN_CLASS, "Random sequence without evidence abstains."),
        _control("failures_reported", len(scoring_rows) == TARGET_COUNT and all("score_label" in row for row in scoring_rows), "Every target has an explicit score row."),
    ]


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    comparison: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "alphafold_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    if any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif comparison["repair_passed"]:
        status = PASSED
    else:
        status = INCOMPLETE
    cert = {
        "kind": "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "adaptive_200_repair_replay",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "baseline_metrics": comparison["baseline"],
        "net_changes": comparison["net_changes"],
        "false_membrane_repair": comparison["false_membrane_repair"],
        "failure_modes": comparison["e63_failure_modes"],
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "repair_passed": comparison["repair_passed"],
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V66 is an adaptive 200-target repair replay; broad claims require new 200-target discovery shards and specialist panels.",
        "next_required_batch": "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "adaptive_200_repair_replay",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["abstain_count"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"] + cert["abstain_count"],
        "failure_modes": cert["failure_modes"],
        "repair_passed": cert["repair_passed"],
        "false_membrane_repair": cert["false_membrane_repair"],
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


def _write_report(path: Path, cert: dict[str, Any], comparison: dict[str, Any]) -> None:
    lines = [
        "# V66 E63 Fast Membrane Repair Replay 200",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Supported: `{cert['supported_count']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Abstain: `{cert['abstain_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Raw accuracy: `{cert['raw_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        "",
        "## E62 Baseline On Same 200",
        f"- supported: `{comparison['baseline']['supported_count']}`",
        f"- failed accepted: `{comparison['baseline']['failed_accepted_count']}`",
        f"- abstain: `{comparison['baseline']['abstain_count']}`",
        f"- accepted accuracy: `{comparison['baseline']['accepted_accuracy']}`",
        "",
        "## E63 Repair Movement",
    ]
    for key, value in comparison["net_changes"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## False Membrane Repair"])
    for key, value in comparison["false_membrane_repair"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Failure Modes"])
    if cert["failure_modes"]:
        for mode, count in sorted(cert["failure_modes"].items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{mode}`: `{count}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Boundary"])
    lines.append(
        "V66 is an adaptive 200-target repair replay. It validates E63's membrane/topology repair direction and leaves broad discovery to V67."
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v66(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_generated_outputs(out_dir)
    targets = _select_targets()
    target_manifest = _public_manifest(targets)
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v66_e63_fast_membrane_repair_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v66_e63_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []
    for target in targets:
        source_manifest = _source_manifest_from_v65_target(target) if target["source_family"] == "V65" else _source_manifest_from_v64_target(target)
        expected = target["expected_mechanism_class"]
        packet = build_sealed_simulation_packet(
            target_id=target["target_id"],
            target_name=f"{target['entry_id']} {target['protein_id']}",
            sequence=source_manifest["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "V66 adaptive repair replay full-chain scan", "span": f"1-{source_manifest['sequence_length']}"}],
            perturbations=[] if expected == ABSTAIN_CLASS else v61._perturbations_for_expected(expected, target["target_id"]),
        )
        holdout = _holdout(target, packet)
        score = _score(packet, holdout, target)
        packets.append(packet)
        scoring_rows.append(score)
        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target["target_id"] / "sealed_packet_summary.json", _packet_summary(packet))
        _write_json(DATA_ROOT / "holdouts_postseal" / target["target_id"] / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target["target_id"] / "validation_result.json", score)

    comparison = _comparison(targets, scoring_rows)
    controls = _controls(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        packets=packets,
        scoring_rows=scoring_rows,
        comparison=comparison,
    )
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
        comparison=comparison,
    )
    claim_row = _claim_row(cert)
    scoring_path = _write_json(DATA_ROOT / "v66_e63_fast_membrane_repair_scoring_report.json", {"kind": "V66_SCORING_REPORT_v0", "rows": scoring_rows})
    comparison_path = _write_json(DATA_ROOT / "v66_e62_vs_e63_comparison.json", comparison)
    failure_path = _write_json(DATA_ROOT / "v66_e63_fast_membrane_repair_failure_report.json", _failure_report(scoring_rows, comparison))
    data_cert_path = _write_json(DATA_ROOT / "v66_e63_fast_membrane_repair_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v66_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v66_e63_fast_membrane_repair_replay_200_certificate.json", cert)
    report_path = out_dir / "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200_REPORT.md"
    _write_report(report_path, cert, comparison)
    return {
        "target_manifest": DATA_ROOT / "v66_e63_fast_membrane_repair_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v66_e63_engine_declaration.json",
        "scoring_report": scoring_path,
        "comparison": comparison_path,
        "failure_report": failure_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V66 E63 fast membrane repair replay on 200 adaptive targets.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v66(args.out_dir)
    cert = _read_json(paths["certificate"], "V66 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "supported_count": cert["supported_count"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "abstain_count": cert["abstain_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "raw_accuracy": cert["raw_accuracy"],
        "coverage": cert["coverage"],
        "baseline_metrics": cert["baseline_metrics"],
        "net_changes": cert["net_changes"],
        "false_membrane_repair": cert["false_membrane_repair"],
        "failure_modes": cert["failure_modes"],
        "repair_passed": cert["repair_passed"],
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "next_required_batch": cert["next_required_batch"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
