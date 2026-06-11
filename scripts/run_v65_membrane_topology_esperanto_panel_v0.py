#!/usr/bin/env python3
from __future__ import annotations

"""Run V65: membrane topology Esperanto panel against the E62 grammar line.

V65 is a specialized topology panel.  It uses real RCSB/V63 sequences and
annotation-provider names already present in the V63 intake cache, then splits
the panel into true transmembrane proteins and the dangerous decoys:
hydrophobic soluble sequences, cofactor-buried hydrophobic pockets, oligomeric
interface hydrophobicity, and monotopic/peripheral membrane association.
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


BATCH_ID = "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E62"
TARGETS_PER_GROUP = 14
PANEL_GROUPS = [
    "A_TRUE_TRANSMEMBRANE_TOPOLOGY",
    "B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY",
    "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET",
    "D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY",
    "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE",
]

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"
V63_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_target_manifest.json"

MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"

PASSED = "V65_MEMBRANE_TOPOLOGY_PANEL_PASSED_REVIEW_REQUIRED"
FAILURES = "V65_MEMBRANE_TOPOLOGY_PANEL_FAILURES_REVIEW_REQUIRED"
BLOCKED_CONTROLS = "V65_MEMBRANE_TOPOLOGY_PANEL_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V65_MEMBRANE_TOPOLOGY_PANEL_BLOCKED_FOR_LEAKAGE"


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
        "v65_membrane_topology_panel_manifest.json",
        "v65_e62_engine_declaration.json",
        "v65_membrane_topology_scoring_report.json",
        "v65_membrane_topology_failure_report.json",
        "v65_membrane_topology_certificate.json",
        "v65_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _annotation_text(candidate: dict[str, Any]) -> str:
    values = [
        " ".join(candidate.get("annotations", []) or []),
        " ".join(candidate.get("feature_types", []) or []),
        candidate.get("title", ""),
        candidate.get("entity_description", ""),
        candidate.get("entry_keywords", ""),
        candidate.get("polymer_composition", ""),
    ]
    return " ".join(str(value) for value in values).lower()


def _has_true_tm_provider(candidate: dict[str, Any]) -> bool:
    text = _annotation_text(candidate)
    if "monotopic/peripheral" in text:
        return False
    return any(token in text for token in ["generic pdbtm", "generic memprotmd", "pdbtm", "memprotmd"])


def _has_peripheral_provider(candidate: dict[str, Any]) -> bool:
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


def _panel_source_statement(group: str, candidate: dict[str, Any]) -> str:
    if group == "A_TRUE_TRANSMEMBRANE_TOPOLOGY":
        text = _annotation_text(candidate)
        if "beta-barrel" in text:
            subtype = "transmembrane beta barrel"
        elif "alpha-helical" in text or "transmembrane" in text:
            subtype = "transmembrane helix"
        else:
            subtype = "provider-classified transmembrane segment"
        return (
            "V65 topology source from OPM/PDBTM/MemProtMD annotations. "
            f"membrane_context_strong transmembrane_context topology_evidence {subtype}; "
            "pore-facing versus lipid-facing residues and inside/outside topology must be treated as topology evidence, not hydrophobicity alone."
        )
    if group == "B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY":
        return (
            "V65 soluble hydrophobic decoy. soluble_hydrophobic_core_context hydrophobicity-alone signal; "
            "no membrane topology evidence, no OPM/PDBTM/MemProtMD transmembrane assignment, and no bilayer-spanning topology source."
        )
    if group == "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET":
        components = " ".join(str(value) for value in candidate.get("biological_cofactor_components", []) or [])
        return (
            "V65 soluble cofactor decoy. cofactor_context ligand_context cofactor-buried hydrophobic pocket; "
            f"biological cofactors/components: {components}. No membrane topology evidence is provided for prediction."
        )
    if group == "D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY":
        return (
            "V65 soluble oligomer decoy. oligomer_context assembly_context partner_copy_context oligomeric interface hydrophobicity; "
            "no membrane topology evidence is provided for prediction."
        )
    return (
        "V65 OPM monotopic/peripheral source. peripheral membrane association, amphipathic peripheral helix or lipid-facing surface; "
        "not transmembrane, no bilayer-spanning topology evidence, and no inside/outside transmembrane topology assignment."
    )


def _expected_for_group(group: str) -> str:
    if group == "A_TRUE_TRANSMEMBRANE_TOPOLOGY":
        return MEMBRANE_CLASS
    if group == "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET":
        return "cofactor_ligand_assisted_stabilization"
    if group == "D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY":
        return "oligomerization_controlled_folding"
    return ABSTAIN_CLASS


def _expected_observables(expected: str) -> list[dict[str, Any]]:
    if expected == ABSTAIN_CLASS:
        return []
    return v61._expected_observables(expected)


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate"]
    target_id = target["target_id"]
    statement = _panel_source_statement(target["panel_group"], candidate)
    return {
        "kind": "V65_MEMBRANE_TOPOLOGY_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "panel_group": target["panel_group"],
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "topology_source_policy": {
            "uses_membrane_specific_provider_names": target["panel_group"] in {
                "A_TRUE_TRANSMEMBRANE_TOPOLOGY",
                "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE",
            },
            "provider_names_seen_in_v63_annotations": [
                provider
                for provider in ["OPM", "PDBTM", "MemProtMD"]
                if provider.lower() in _annotation_text(candidate)
            ],
            "dangerous_confusion_tested": target["dangerous_confusion_tested"],
            "expected_panel_distinction": target["expected_panel_distinction"],
        },
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
                "source_id": f"{target_id}_V65_TOPOLOGY_PANEL_SOURCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": statement,
                "source_url": candidate["source_urls"]["entry"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
            "ligand coordinates, metal coordination geometry, and bound-state contact geometry",
            "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
            "post-seal validation annotations before prediction hash",
            "V64 score outcomes as prediction evidence",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _panel_target(group: str, candidate: dict[str, Any], ordinal: int) -> dict[str, Any]:
    expected = _expected_for_group(group)
    confusion = {
        "A_TRUE_TRANSMEMBRANE_TOPOLOGY": "true membrane topology should not be hidden by cofactor/oligomer context",
        "B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY": "hydrophobicity alone must not become membrane topology",
        "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET": "cofactor-buried hydrophobicity must not become membrane topology",
        "D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY": "oligomeric interface hydrophobicity must not become membrane topology",
        "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE": "peripheral membrane association must not become transmembrane grammar",
    }[group]
    distinction = {
        "A_TRUE_TRANSMEMBRANE_TOPOLOGY": "explicit topology evidence -> membrane grammar",
        "B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY": "hydrophobicity alone -> clean abstain",
        "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET": "cofactor pocket explains hydrophobicity -> cofactor grammar",
        "D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY": "interface explains hydrophobicity -> oligomer grammar",
        "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE": "peripheral/monotopic association -> clean abstain until a peripheral class exists",
    }[group]
    return {
        "target_id": f"V65_{group[0]}_{ordinal:02d}_{candidate['target_id']}",
        "panel_group": group,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "target_name": f"{candidate['entry_id']} {candidate['entity_description']}",
        "sequence_length": candidate["sequence_length"],
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "sequence_metrics": candidate["sequence_metrics"],
        "expected_mechanism_class": expected,
        "dangerous_confusion_tested": confusion,
        "expected_panel_distinction": distinction,
        "candidate": candidate,
    }


def _select_panel(v63_manifest: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = [dict(row) for row in v63_manifest.get("selected_targets", []) if isinstance(row, dict)]
    for candidate in candidates:
        candidate["sequence_metrics"] = v61._sequence_metrics(candidate["sequence"])
    selected: list[dict[str, Any]] = []
    used: set[str] = set()

    def add_group(group: str, pool: list[dict[str, Any]]) -> None:
        pool = [row for row in pool if row["target_id"] not in used]
        if len(pool) < TARGETS_PER_GROUP:
            raise SystemExit(f"V65 group {group} has only {len(pool)} candidates; need {TARGETS_PER_GROUP}")
        for ordinal, candidate in enumerate(pool[:TARGETS_PER_GROUP], start=1):
            used.add(candidate["target_id"])
            selected.append(_panel_target(group, candidate, ordinal))

    true_tm = sorted([row for row in candidates if _has_true_tm_provider(row)], key=lambda row: row["target_id"])
    peripheral = sorted([row for row in candidates if _has_peripheral_provider(row)], key=lambda row: row["target_id"])
    no_membrane_provider = [row for row in candidates if not _has_true_tm_provider(row) and not _has_peripheral_provider(row)]
    cofactor = sorted(
        [row for row in no_membrane_provider if _has_cofactor_context(row)],
        key=lambda row: (-row["sequence_metrics"]["max_segment_membrane_density"], row["target_id"]),
    )
    oligomer = sorted(
        [row for row in no_membrane_provider if _has_oligomer_context(row) and not _has_cofactor_context(row)],
        key=lambda row: (-row["sequence_metrics"]["max_segment_membrane_density"], row["target_id"]),
    )
    hydrophobic = sorted(
        [row for row in no_membrane_provider],
        key=lambda row: (-row["sequence_metrics"]["max_segment_membrane_density"], row["target_id"]),
    )

    add_group("A_TRUE_TRANSMEMBRANE_TOPOLOGY", true_tm)
    add_group("E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE", peripheral)
    add_group("C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET", cofactor)
    add_group("D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY", oligomer)
    add_group("B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY", hydrophobic)
    order = {group: index for index, group in enumerate(PANEL_GROUPS)}
    return sorted(selected, key=lambda row: (order[row["panel_group"]], row["target_id"]))


def _panel_manifest(v63_manifest: dict[str, Any]) -> dict[str, Any]:
    targets = _select_panel(v63_manifest)
    public_targets = []
    for target in targets:
        row = dict(target)
        row.pop("candidate", None)
        public_targets.append(row)
    groups = Counter(target["panel_group"] for target in targets)
    return {
        "kind": "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "source_manifest": str(V63_MANIFEST.relative_to(REPO_ROOT)),
        "target_selection_manual": False,
        "targets_per_group": TARGETS_PER_GROUP,
        "panel_target_count": len(targets),
        "panel_groups": dict(groups),
        "required_esperanto_distinctions": [
            "transmembrane helix",
            "transmembrane beta barrel",
            "re-entrant loop",
            "amphipathic peripheral helix",
            "signal peptide",
            "lipid anchor",
            "soluble hydrophobic core",
            "cofactor-buried hydrophobic pocket",
            "oligomeric interface hydrophobicity",
            "pore-facing versus lipid-facing residues",
            "inside/outside topology",
        ],
        "targets": public_targets,
        "_runtime_targets": targets,
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V65_E62_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_lineage": ["E60", "E61", "E62"],
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


def _holdout(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate"]
    expected = target["expected_mechanism_class"]
    return {
        "kind": "V65_MEMBRANE_TOPOLOGY_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "panel_group": target["panel_group"],
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "expected_mechanism_class": expected,
        "expected_observables": _expected_observables(expected),
        "postseal_truth_basis": [
            target["expected_panel_distinction"],
            target["dangerous_confusion_tested"],
        ],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V65_PANEL_HOLDOUT",
                "source_class": "coordinate_derived",
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "entry_url": candidate["source_urls"]["entry"],
                "polymer_entity_url": candidate["source_urls"]["polymer_entity"],
            }
        ],
    }


def _score(packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == ABSTAIN_CLASS else "accepted"
    supported = predicted == expected and validation["score_label"] == "supported"
    return {
        "kind": "V65_MEMBRANE_TOPOLOGY_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "panel_group": holdout["panel_group"],
        "protein_id": holdout["protein_id"],
        "entry_id": holdout["entry_id"],
        "entity_id": holdout["entity_id"],
        "acceptance_decision": decision,
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "level1_regime_selection": predicted == expected,
        "level2_region_localization_proxy": decision == "accepted" and bool(packet["operator_field"]["operators"]),
        "level3_topology_or_contact_proxy": supported,
        "score_label": "supported" if supported else ("abstained" if decision == "abstain_recommended" else "contradicted"),
        "validation_checks": validation["checks"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout["holdout_opened_after_prediction_hash"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V65_COMPACT_SEALED_PACKET_SUMMARY_v0",
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


def _metrics(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in scoring_rows if row["acceptance_decision"] == "accepted"]
    supported = [row for row in scoring_rows if row["score_label"] == "supported"]
    failed_accepted = [row for row in accepted if row["score_label"] != "supported"]
    abstained = [row for row in scoring_rows if row["acceptance_decision"] == "abstain_recommended"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    false_membrane = [
        row for row in scoring_rows
        if row["predicted_mechanism_class"] == MEMBRANE_CLASS and row["expected_mechanism_class"] != MEMBRANE_CLASS
    ]
    true_tm_missed = [
        row for row in scoring_rows
        if row["panel_group"] == "A_TRUE_TRANSMEMBRANE_TOPOLOGY" and row["score_label"] != "supported"
    ]
    peripheral_false_tm = [
        row for row in scoring_rows
        if row["panel_group"] == "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE"
        and row["predicted_mechanism_class"] == MEMBRANE_CLASS
    ]
    return {
        "targets_total": len(scoring_rows),
        "accepted_count": len(accepted),
        "supported_count": len(supported),
        "failed_accepted_count": len(failed_accepted),
        "abstain_count": len(abstained),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "raw_accuracy": len(supported) / len(scoring_rows) if scoring_rows else None,
        "coverage": len(accepted) / len(scoring_rows) if scoring_rows else None,
        "false_membrane_call_count": len(false_membrane),
        "true_transmembrane_missed_count": len(true_tm_missed),
        "peripheral_false_transmembrane_count": len(peripheral_false_tm),
    }


def _failure_type(row: dict[str, Any]) -> str:
    if row["score_label"] == "supported":
        return "supported"
    if row["predicted_mechanism_class"] == MEMBRANE_CLASS and row["expected_mechanism_class"] != MEMBRANE_CLASS:
        if row["panel_group"] == "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE":
            return "peripheral_misread_as_transmembrane"
        return "soluble_hydrophobic_false_membrane"
    if row["expected_mechanism_class"] == MEMBRANE_CLASS and row["predicted_mechanism_class"] != MEMBRANE_CLASS:
        return "true_transmembrane_missed"
    if row["expected_mechanism_class"] == ABSTAIN_CLASS and row["predicted_mechanism_class"] != ABSTAIN_CLASS:
        return "ambiguous_topology_overaccepted"
    return v61._failure_type(row)


def _failure_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    failures = [row for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V65_MEMBRANE_TOPOLOGY_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_cases_reported": True,
        "failure_count": len(failures),
        "failure_modes": dict(Counter(_failure_type(row) for row in failures)),
        "failure_cases": failures,
        "topology_confusions": [
            {
                "target_id": row["target_id"],
                "panel_group": row["panel_group"],
                "failure_type": _failure_type(row),
                "predicted": row["predicted_mechanism_class"],
                "expected": row["expected_mechanism_class"],
                "score_label": row["score_label"],
            }
            for row in failures
        ],
        "note": "V65 is a topology-specific failure panel; failures are E63 grammar-mining material.",
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    panel_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V65_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V65_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V65_PRESEAL_HOLDOUT",
        "source_class": "coordinate_derived",
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    random_packet = build_sealed_simulation_packet(
        target_id="V65_RANDOM_SEQUENCE_CONTROL",
        target_name="V65 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
    groups = panel_manifest["panel_groups"]
    distinctions = set(panel_manifest["required_esperanto_distinctions"])
    required = {
        "transmembrane helix",
        "transmembrane beta barrel",
        "amphipathic peripheral helix",
        "signal peptide",
        "lipid anchor",
        "soluble hydrophobic core",
        "cofactor-buried hydrophobic pocket",
        "oligomeric interface hydrophobicity",
        "inside/outside topology",
    }
    return [
        _control("panel_groups_present", all(groups.get(group) == TARGETS_PER_GROUP for group in PANEL_GROUPS), "V65 must split the panel into the five requested groups.", groups),
        _control("target_selection_manual_false", panel_manifest["target_selection_manual"] is False, "V65 target split is deterministic from V63 annotations and sequence-field metrics."),
        _control("provider_specific_topology_sources_present", groups.get("A_TRUE_TRANSMEMBRANE_TOPOLOGY") == TARGETS_PER_GROUP and groups.get("E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE") == TARGETS_PER_GROUP, "True and peripheral membrane groups use OPM/PDBTM/MemProtMD provider names from the V63 intake."),
        _control("esperanto_distinctions_enumerated", required <= distinctions, "V65 records the membrane grammar distinctions required before broad expansion."),
        _control("engine_version_declared_e62", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V65 uses E62."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside V65."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V65 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == ABSTAIN_CLASS, "Random sequence without evidence abstains."),
        _control("failures_reported", len(scoring_rows) == TARGETS_PER_GROUP * len(PANEL_GROUPS) and all("score_label" in row for row in scoring_rows), "Every V65 target has an explicit score row."),
    ]


def _aggregate_certificate(
    *,
    panel_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    if any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif metrics["failed_accepted_count"] or metrics["true_transmembrane_missed_count"]:
        status = FAILURES
    else:
        status = PASSED
    by_group: dict[str, dict[str, int]] = {}
    for group in PANEL_GROUPS:
        rows = [row for row in scoring_rows if row["panel_group"] == group]
        by_group[group] = dict(Counter(row["score_label"] for row in rows))
    confusion_matrix: dict[str, dict[str, int]] = defaultdict(dict)
    for row in scoring_rows:
        expected = row["expected_mechanism_class"]
        predicted = row["predicted_mechanism_class"]
        confusion_matrix[expected][predicted] = confusion_matrix[expected].get(predicted, 0) + 1
    failure_modes = dict(Counter(_failure_type(row) for row in scoring_rows if row["score_label"] != "supported"))
    cert = {
        "kind": "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "membrane_topology_panel",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": panel_manifest["target_selection_manual"],
        "panel_groups": panel_manifest["panel_groups"],
        **metrics,
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "score_labels_by_group": by_group,
        "failure_modes": failure_modes,
        "topology_confusion_matrix": {key: dict(value) for key, value in confusion_matrix.items()},
        "engine_revision_required": bool(failure_modes),
        "engine_revision_recommended": "E63_MEMBRANE_TOPOLOGY_GRAMMAR_REVISION" if failure_modes else None,
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V65 is a membrane topology panel and grammar-mining gate; broad claims require a repaired engine and later expansion.",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "membrane_topology_panel",
        "engine_version_used": ENGINE_VERSION_USED,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["abstain_count"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"] + cert["abstain_count"],
        "failure_modes": cert["failure_modes"],
        "engine_revision_required": cert["engine_revision_required"],
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


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V65 Membrane Topology Esperanto Panel",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Supported: `{cert['supported_count']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Abstain: `{cert['abstain_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Raw accuracy: `{cert['raw_accuracy']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"False membrane calls: `{cert['false_membrane_call_count']}`",
        f"Peripheral false transmembrane calls: `{cert['peripheral_false_transmembrane_count']}`",
        f"True transmembrane missed: `{cert['true_transmembrane_missed_count']}`",
        "",
        "## Groups",
    ]
    for group, count in cert["panel_groups"].items():
        labels = cert["score_labels_by_group"].get(group, {})
        lines.append(f"- `{group}`: `{count}` targets, labels `{labels}`")
    lines.extend(["", "## Failure Modes"])
    if cert["failure_modes"]:
        for mode, count in sorted(cert["failure_modes"].items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{mode}`: `{count}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Boundary"])
    lines.append(
        "V65 is a membrane-topology grammar panel. It distinguishes topology evidence from hydrophobicity, cofactor pockets, oligomeric interfaces, and peripheral membrane association; it does not make a broad folding claim."
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v65(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    v63_manifest = _read_json(V63_MANIFEST, "V63 target manifest")
    _reset_generated_outputs(out_dir)
    panel_manifest = _panel_manifest(v63_manifest)
    runtime_targets = panel_manifest.pop("_runtime_targets")
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v65_membrane_topology_panel_manifest.json", panel_manifest)
    _write_json(DATA_ROOT / "v65_e62_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []
    for target in runtime_targets:
        source_manifest = _source_manifest(target)
        packet = build_sealed_simulation_packet(
            target_id=target["target_id"],
            target_name=target["target_name"],
            sequence=target["candidate"]["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "V65 topology panel full-chain scan", "span": f"1-{target['sequence_length']}"}],
            perturbations=v61._perturbations_for_expected(target["expected_mechanism_class"], target["target_id"]),
        )
        holdout = _holdout(target, packet)
        score = _score(packet, holdout)
        packets.append(packet)
        scoring_rows.append(score)
        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target["target_id"] / "sealed_packet_summary.json", _packet_summary(packet))
        _write_json(DATA_ROOT / "holdouts_postseal" / target["target_id"] / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target["target_id"] / "validation_result.json", score)

    controls = _controls(
        panel_manifest=panel_manifest,
        engine_declaration=engine_declaration,
        packets=packets,
        scoring_rows=scoring_rows,
    )
    cert = _aggregate_certificate(
        panel_manifest=panel_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
    )
    claim_row = _claim_row(cert)
    scoring_path = _write_json(DATA_ROOT / "v65_membrane_topology_scoring_report.json", {"kind": "V65_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v65_membrane_topology_failure_report.json", _failure_report(scoring_rows))
    data_cert_path = _write_json(DATA_ROOT / "v65_membrane_topology_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v65_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v65_membrane_topology_esperanto_panel_certificate.json", cert)
    report_path = out_dir / "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL_REPORT.md"
    _write_report(report_path, cert)
    return {
        "panel_manifest": DATA_ROOT / "v65_membrane_topology_panel_manifest.json",
        "engine_declaration": DATA_ROOT / "v65_e62_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V65 membrane topology Esperanto panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v65(args.out_dir)
    cert = _read_json(paths["certificate"], "V65 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "supported_count": cert["supported_count"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "abstain_count": cert["abstain_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "raw_accuracy": cert["raw_accuracy"],
        "controls_passed": cert["controls_passed"],
        "false_membrane_call_count": cert["false_membrane_call_count"],
        "peripheral_false_transmembrane_count": cert["peripheral_false_transmembrane_count"],
        "true_transmembrane_missed_count": cert["true_transmembrane_missed_count"],
        "failure_modes": cert["failure_modes"],
        "engine_revision_required": cert["engine_revision_required"],
        "engine_revision_recommended": cert["engine_revision_recommended"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
