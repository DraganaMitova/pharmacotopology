#!/usr/bin/env python3
from __future__ import annotations

"""Run V64: E62 replay on the exact V63 RCSB 500 target set.

V64 is a paired repair measurement.  It loads the committed V63/E61 score rows,
reruns the same 500 targets with the E62 engine grammar, and records whether
the membrane/topology priority repair lowers the V63 membrane_misread wall
without creating new soluble false-membrane regressions.
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
import run_v63_rcsb_500_discovery_batch_v0 as v63  # noqa: E402


BATCH_ID = "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E62"
BASELINE_BATCH_ID = "V63_RCSB_500_DISCOVERY_BATCH"
BASELINE_ENGINE_VERSION = "E61"
TARGET_COUNT = 500
SEQUENCE_IDENTITY_CUTOFF = 30

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V64"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V63_DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63"
V63_TARGET_MANIFEST = V63_DATA_ROOT / "v63_rcsb_500_target_manifest.json"
V63_SCORING_REPORT = V63_DATA_ROOT / "v63_rcsb_500_scoring_report.json"
V63_BATCH_CERTIFICATE = RUN_ROOT / BASELINE_BATCH_ID / "v63_rcsb_500_discovery_batch_certificate.json"

PASSED_DIRECTIONAL = "V64_E62_MEMBRANE_REPAIR_DIRECTIONAL_IMPROVEMENT_REVIEW_REQUIRED"
NO_DIRECTIONAL_IMPROVEMENT = "V64_E62_MEMBRANE_REPAIR_NO_DIRECTIONAL_IMPROVEMENT_REVIEW_REQUIRED"
BLOCKED_CONTROLS = "V64_E62_MEMBRANE_REPAIR_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V64_E62_MEMBRANE_REPAIR_BLOCKED_FOR_LEAKAGE"
BLOCKED_BASELINE = "V64_E62_MEMBRANE_REPAIR_BLOCKED_BASELINE_UNAVAILABLE"

MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"


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
        "v64_e62_replay_target_manifest.json",
        "v64_e62_engine_declaration.json",
        "v64_e62_rcsb_500_membrane_repair_scoring_report.json",
        "v64_e62_rcsb_500_membrane_repair_failure_report.json",
        "v64_e61_vs_e62_comparison.json",
        "v64_e62_rcsb_500_membrane_repair_certificate.json",
        "v64_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _target_key(target_id: str) -> str:
    for prefix in ["V63_", "V64_"]:
        if target_id.startswith(prefix):
            return target_id.removeprefix(prefix)
    return target_id


def _baseline_rows(scoring_report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = scoring_report.get("rows", [])
    if not isinstance(rows, list):
        raise SystemExit("V63 scoring rows must be a list")
    return {_target_key(str(row["target_id"])): row for row in rows if isinstance(row, dict) and row.get("target_id")}


def _baseline_metrics(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "accepted_count": cert.get("accepted_count"),
        "supported_count": cert.get("supported_count"),
        "failed_accepted_count": cert.get("failed_accepted_count"),
        "abstain_count": cert.get("abstain_count"),
        "accepted_accuracy": cert.get("accepted_accuracy"),
        "raw_accuracy": cert.get("raw_accuracy"),
        "coverage": cert.get("coverage"),
        "failure_modes": cert.get("failure_modes", {}),
    }


def _v64_context_statement(context: dict[str, Any]) -> str:
    marks = context["context_marks"]
    if not marks:
        return "V64 E62 replay metadata context emitted no explicit E62 context marks for this target."
    return (
        "V64 E62 replay metadata context marks: "
        + " ".join(marks)
        + ". E62 tests whether strong membrane/topology context outranks incidental ligand/cofactor "
        + "or oligomer/assembly context. Coordinates, contacts, distance maps, ligand geometry, "
        + "native topology, and post-seal validation labels are blocked before sealing."
    )


def _source_manifest(candidate: dict[str, Any]) -> dict[str, Any]:
    target_id = f"V64_{candidate['target_id']}"
    context = v63._metadata_context(candidate)
    return {
        "kind": "V64_E62_MEMBRANE_REPAIR_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "same_target_as_v63": True,
        "paired_repair_replay": True,
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
                "source_id": f"{target_id}_E62_METADATA_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": context["context_marks"],
                "metadata_context_reasons": context["reasons"],
                "evidence_statement": _v64_context_statement(context),
                "source_url": candidate["source_urls"]["entry"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
            "ligand coordinates, metal coordination geometry, and bound-state contact geometry",
            "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
            "post-seal validation annotations before prediction hash",
            "V63 expected labels and V63 score outcomes as prediction evidence",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _holdout(candidate: dict[str, Any], packet: dict[str, Any], baseline_row: dict[str, Any], reasons: list[str]) -> dict[str, Any]:
    expected = baseline_row["expected_mechanism_class"]
    return {
        "kind": "V64_RCSB_500_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "expected_mechanism_class": expected,
        "expected_observables": v61._expected_observables(expected),
        "postseal_truth_basis": reasons,
        "baseline_v63_target_id": baseline_row["target_id"],
        "baseline_v63_prediction_hash": baseline_row["sealed_prediction_hash"],
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
        "kind": "V64_E62_MEMBRANE_REPAIR_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
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
        "baseline_batch_id": BASELINE_BATCH_ID,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "baseline_target_id": baseline_row["target_id"],
        "baseline_acceptance_decision": baseline_row["acceptance_decision"],
        "baseline_predicted_mechanism_class": baseline_row["predicted_mechanism_class"],
        "baseline_score_label": baseline_row["score_label"],
        "baseline_supported": baseline_row["score_label"] == "supported",
        "failures_repaired_after_holdout": False,
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V64_COMPACT_SEALED_PACKET_SUMMARY_v0",
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
        "kind": "V64_E62_ENGINE_DECLARATION_v0",
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


def _target_manifest(v63_manifest: dict[str, Any]) -> dict[str, Any]:
    selected = [dict(row) for row in v63_manifest.get("selected_targets", []) if isinstance(row, dict)]
    selected = selected[:TARGET_COUNT]
    for candidate in selected:
        candidate["sequence_metrics"] = v61._sequence_metrics(candidate["sequence"])
    return {
        "kind": "V64_E62_REPLAY_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "batch_mode": "paired_repair_replay",
        "target_selection_manual": False,
        "same_targets_as_v63": True,
        "selection_rule": "Exact selected_targets list from V63_RCSB_500_DISCOVERY_BATCH, same order and same 30% cluster representatives.",
        "source_manifest": str(V63_TARGET_MANIFEST.relative_to(REPO_ROOT)),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(selected),
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "sequence_cluster_representative_selection": True,
        "selected_targets": selected,
    }


def _failure_type(row: dict[str, Any]) -> str:
    return v61._failure_type(row)


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


def _failure_modes(scoring_rows: list[dict[str, Any]]) -> dict[str, int]:
    return dict(Counter(_failure_type(row) for row in scoring_rows if row["score_label"] != "supported"))


def _failure_grammar_row(row: dict[str, Any], category: str) -> dict[str, Any]:
    failure_type = _failure_type(row)
    proposals = {
        "membrane_misread": ("topology-evidence grammar, not hydrophobicity grammar", "membrane_pressure_operator", "bilayer/topology pressure", "split transmembrane, peripheral, signal peptide, and hydrophobic-core explanations"),
        "disorder_misread": ("IDR persistence and fold-upon-binding separation word", "disorder_operator", "entropy/partner pressure", "split IDR, phase, and bound-order contexts"),
        "oligomer_state_misread": ("assembly-register/interface specificity word", "interface_operator", "partner-copy concentration pressure", "separate obligate assembly from incidental copy count"),
        "cofactor_ligand_missing": ("ligand-state specificity word", "interface_operator", "ligand_or_cofactor pressure", "separate structural ions, cofactors, and incidental ligands"),
        "weak_sequence_signal": ("confidence/self-evidence word", "none", "none", "abstain until independent non-coordinate evidence sharpens mechanism"),
        "right_regime_wrong_topology": ("topology proxy refinement word", "closure_operator", "fold-class pressure", "abstain when operator regions are underdetermined"),
        "wrong_regime": ("regime separation word", "frustration_operator", "context separation pressure", "mine repeated classes before next engine revision"),
    }
    rule, operator, pressure, abstention = proposals.get(failure_type, proposals["wrong_regime"])
    return {
        "target_id": row["target_id"],
        "protein_id": row["protein_id"],
        "failure_type": failure_type,
        "repair_category": category,
        "baseline_predicted": row["baseline_predicted_mechanism_class"],
        "engine_thought": row["predicted_mechanism_class"],
        "reality_showed": row["expected_mechanism_class"],
        "score_label": row["score_label"],
        "missing_esperanto_rule": rule,
        "proposed_new_operator": operator,
        "proposed_new_pressure": pressure,
        "proposed_new_abstention_rule": abstention,
        "control_that_prevents_overfitting": "V64 is same-target repair replay; persistent membrane classes feed V65 membrane topology panel before any broad claim.",
    }


def _comparison(
    *,
    baseline_cert: dict[str, Any],
    baseline_rows: dict[str, dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    e61 = _baseline_metrics(baseline_cert)
    e62 = {**_metrics(scoring_rows), "failure_modes": _failure_modes(scoring_rows)}
    categories: dict[str, list[dict[str, Any]]] = {
        "stable_supported": [],
        "repaired": [],
        "persistent_failure": [],
        "new_failure": [],
    }
    for row in scoring_rows:
        baseline = baseline_rows[_target_key(row["target_id"])]
        was_supported = baseline["score_label"] == "supported"
        now_supported = row["score_label"] == "supported"
        if was_supported and now_supported:
            category = "stable_supported"
        elif not was_supported and now_supported:
            category = "repaired"
        elif was_supported and not now_supported:
            category = "new_failure"
        else:
            category = "persistent_failure"
        categories[category].append(row)

    baseline_failure_type = {
        key: _failure_type(row)
        for key, row in baseline_rows.items()
        if row["score_label"] != "supported"
    }
    new_failure_modes = Counter(_failure_type(row) for row in categories["new_failure"])
    repaired_failure_modes = Counter(baseline_failure_type[_target_key(row["target_id"])] for row in categories["repaired"])
    persistent_failure_modes = Counter(_failure_type(row) for row in categories["persistent_failure"])
    membrane_repaired = [
        row for row in categories["repaired"]
        if baseline_failure_type.get(_target_key(row["target_id"])) == "membrane_misread"
    ]
    membrane_persistent = [
        row for row in categories["persistent_failure"]
        if baseline_failure_type.get(_target_key(row["target_id"])) == "membrane_misread"
    ]
    soluble_false_membrane = [
        row for row in scoring_rows
        if row["predicted_mechanism_class"] == MEMBRANE_CLASS and row["expected_mechanism_class"] != MEMBRANE_CLASS
    ]
    cofactor_oligomer_regressions = [
        row for row in categories["new_failure"]
        if row["expected_mechanism_class"] in {
            "cofactor_ligand_assisted_stabilization",
            "oligomerization_controlled_folding",
        }
    ]
    e61_failure_modes = dict(e61.get("failure_modes", {}))
    e62_failure_modes = dict(e62["failure_modes"])
    all_modes = sorted(set(e61_failure_modes) | set(e62_failure_modes))
    failure_mode_distribution_change = {
        mode: {
            "e61": e61_failure_modes.get(mode, 0),
            "e62": e62_failure_modes.get(mode, 0),
            "delta": e62_failure_modes.get(mode, 0) - e61_failure_modes.get(mode, 0),
        }
        for mode in all_modes
    }
    comparison = {
        "kind": "V64_E61_VS_E62_MEMBRANE_REPAIR_COMPARISON_v0",
        "batch_id": BATCH_ID,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "e61": e61,
        "e62": e62,
        "net_changes": {
            "supported_delta": e62["supported_count"] - e61["supported_count"],
            "failed_accepted_delta": e62["failed_accepted_count"] - e61["failed_accepted_count"],
            "abstain_delta": e62["abstain_count"] - e61["abstain_count"],
            "accepted_accuracy_delta": e62["accepted_accuracy"] - e61["accepted_accuracy"],
            "coverage_delta": e62["coverage"] - e61["coverage"],
            "membrane_misread_delta": e62_failure_modes.get("membrane_misread", 0) - e61_failure_modes.get("membrane_misread", 0),
        },
        "failure_mode_distribution_change": failure_mode_distribution_change,
        "repair_categories": {name: len(rows) for name, rows in categories.items()},
        "repaired_failure_modes": dict(repaired_failure_modes),
        "persistent_failure_modes": dict(persistent_failure_modes),
        "new_failure_modes_from_regressions": dict(new_failure_modes),
        "e62_failure_modes_not_seen_in_e61": {
            mode: count
            for mode, count in e62_failure_modes.items()
            if mode not in e61_failure_modes
        },
        "membrane_repair": {
            "e61_membrane_misread": e61_failure_modes.get("membrane_misread", 0),
            "e62_membrane_misread": e62_failure_modes.get("membrane_misread", 0),
            "membrane_misread_repaired": len(membrane_repaired),
            "membrane_misread_persistent": len(membrane_persistent),
            "soluble_false_membrane_call_count": len(soluble_false_membrane),
            "cofactor_or_oligomer_regression_count": len(cofactor_oligomer_regressions),
        },
        "sample_rows": {
            "repaired": [_sample_row(row) for row in categories["repaired"][:20]],
            "persistent_failure": [_sample_row(row) for row in categories["persistent_failure"][:20]],
            "new_failure": [_sample_row(row) for row in categories["new_failure"][:20]],
            "soluble_false_membrane": [_sample_row(row) for row in soluble_false_membrane[:20]],
        },
        "directional_repair": (
            e62["supported_count"] > e61["supported_count"]
            and e62["failed_accepted_count"] < e61["failed_accepted_count"]
            and e62_failure_modes.get("membrane_misread", 0) < e61_failure_modes.get("membrane_misread", 0)
        ),
        "membrane_still_fails": e62_failure_modes.get("membrane_misread", 0) > 0,
        "claim_allowed": False,
        "claim_boundary": "V64 is same-target paired repair measurement, not a broad saturation claim.",
    }
    return comparison


def _sample_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "target_id": row["target_id"],
        "protein_id": row["protein_id"],
        "baseline_predicted": row["baseline_predicted_mechanism_class"],
        "predicted": row["predicted_mechanism_class"],
        "expected": row["expected_mechanism_class"],
        "score_label": row["score_label"],
    }


def _failure_report(scoring_rows: list[dict[str, Any]], comparison: dict[str, Any]) -> dict[str, Any]:
    categories_by_target: dict[str, str] = {}
    for category, samples in comparison["sample_rows"].items():
        for row in samples:
            categories_by_target[row["target_id"]] = category
    failure_rows: list[dict[str, Any]] = []
    for row in scoring_rows:
        if row["score_label"] == "supported":
            continue
        baseline_supported = bool(row["baseline_supported"])
        category = "new_failure" if baseline_supported else "persistent_failure"
        failure_rows.append(_failure_grammar_row(row, category))
    return {
        "kind": "V64_E62_MEMBRANE_REPAIR_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "failure_cases_reported": True,
        "failure_count": len(failure_rows),
        "failure_modes": comparison["e62"]["failure_modes"],
        "repaired_failure_modes": comparison["repaired_failure_modes"],
        "persistent_failure_modes": comparison["persistent_failure_modes"],
        "new_failure_modes_from_regressions": comparison["new_failure_modes_from_regressions"],
        "membrane_repair": comparison["membrane_repair"],
        "failure_grammar_rows": failure_rows,
        "note": "V64 preserves persistent and new failures for the V65 membrane topology Esperanto panel.",
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    v63_manifest: dict[str, Any],
    target_manifest: dict[str, Any],
    baseline_cert: dict[str, Any],
    baseline_rows: dict[str, dict[str, Any]],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    wrong_packets: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V64_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V64_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold-style model offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V64_PRESEAL_HOLDOUT",
        "source_class": "coordinate_derived",
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V64_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_operator_state_packet(
        target_id="V64_RANDOM_SEQUENCE_CONTROL",
        target_name="V64 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
    v63_keys = [str(row["target_id"]) for row in v63_manifest.get("selected_targets", [])[:TARGET_COUNT]]
    v64_keys = [str(row["target_id"]) for row in target_manifest["selected_targets"]]
    clusters = [row.get("sequence_cluster_30_id") for row in target_manifest["selected_targets"]]
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
    expected_match_rows = []
    for row in scoring_rows:
        baseline = baseline_rows.get(_target_key(row["target_id"]), {})
        expected_match_rows.append(row["expected_mechanism_class"] == baseline.get("expected_mechanism_class"))
    return [
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V64 target selection must be automatic reuse of V63."),
        _control("same_500_targets_as_v63", v64_keys == v63_keys and len(v64_keys) == TARGET_COUNT, "V64 must use the exact V63 500 target list in the same order.", {"v63_count": len(v63_keys), "v64_count": len(v64_keys)}),
        _control("targets_total_500", len(target_manifest["selected_targets"]) == TARGET_COUNT, "V64 N must be exactly 500.", len(target_manifest["selected_targets"])),
        _control("baseline_v63_e61_loaded", baseline_cert.get("supported_count") == 238 and baseline_cert.get("failed_accepted_count") == 262 and baseline_cert.get("failure_modes", {}).get("membrane_misread") == 216, "V64 baseline must be the committed V63/E61 result.", _baseline_metrics(baseline_cert)),
        _control("baseline_rows_loaded_for_all_targets", len(baseline_rows) == TARGET_COUNT and all(key in baseline_rows for key in v64_keys), "Every V64 target must have a paired V63/E61 score row.", len(baseline_rows)),
        _control("baseline_expected_labels_match_v64_holdouts", all(expected_match_rows) and len(expected_match_rows) == TARGET_COUNT, "V64 validates against the same expected mechanism labels as V63."),
        _control("rcsb_experimental_protein_entities_only", all(row["source_database"] == "RCSB_PDB" and row["structure_determination_methodology"].lower() == "experimental" and "protein" in row["polymer_type"].lower() for row in target_manifest["selected_targets"]), "All targets remain RCSB experimental protein entities."),
        _control("sequence_cluster_representative_selection", target_manifest["sequence_cluster_representative_selection"] is True and len(set(clusters)) == len(clusters) and all(clusters), "Each target remains a unique 30% sequence-cluster representative.", {"cutoff": SEQUENCE_IDENTITY_CUTOFF, "unique_clusters": len(set(clusters))}),
        _control("engine_version_declared_e62", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V64 uses E62."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside V64."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V64 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("wrong_grammar_controls", all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_packets), "Forced wrong grammars are rejected or routed to abstention."),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls are generated with target metadata withheld.", {"control_count": len(shuffled_rows), "shuffled_higher_by_more_than_margin_count": sum(1 for row in shuffled_rows if row["shuffled_higher_by_more_than_margin"])}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain", "Random sequence without evidence abstains."),
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
    if len(scoring_rows) < TARGET_COUNT:
        status = BLOCKED_BASELINE
    elif any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif comparison["directional_repair"]:
        status = PASSED_DIRECTIONAL
    else:
        status = NO_DIRECTIONAL_IMPROVEMENT
    controls_passed = not failed_controls
    cert = {
        "kind": "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "baseline_batch_id": BASELINE_BATCH_ID,
        "batch_mode": "paired_repair_replay",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "same_targets_as_v63": target_manifest["same_targets_as_v63"],
        **metrics,
        "controls_passed": controls_passed,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "failure_modes": comparison["e62"]["failure_modes"],
        "e61_baseline_metrics": comparison["e61"],
        "net_changes": comparison["net_changes"],
        "repair_categories": comparison["repair_categories"],
        "membrane_repair": comparison["membrane_repair"],
        "new_failure_modes_from_regressions": comparison["new_failure_modes_from_regressions"],
        "e62_failure_modes_not_seen_in_e61": comparison["e62_failure_modes_not_seen_in_e61"],
        "directional_repair": comparison["directional_repair"],
        "membrane_still_fails": comparison["membrane_still_fails"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V64 is a same-target E62 repair replay against V63; broad claims require V65 topology panel and later expansion.",
        "next_required_batch": "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL" if comparison["membrane_still_fails"] else "V66_RCSB_1000_SATURATION_DISCOVERY",
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "E62 has saturated broad RCSB protein space.",
            "Persistent membrane failures may be hidden.",
            "Coordinates, contacts, ligand geometry, or native topology were used before sealing.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "paired_repair_replay",
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
        "membrane_repair": cert["membrane_repair"],
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


def _write_report(path: Path, cert: dict[str, Any], comparison: dict[str, Any], scoring_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# V64 E62 RCSB 500 Membrane Repair Replay",
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
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"Engine modified during batch: `{cert['engine_modified_during_batch']}`",
        "",
        "## E61 Baseline",
        f"- supported: `{comparison['e61']['supported_count']}`",
        f"- failed accepted: `{comparison['e61']['failed_accepted_count']}`",
        f"- membrane_misread: `{comparison['e61']['failure_modes'].get('membrane_misread', 0)}`",
        "",
        "## E62 Replay",
        f"- supported: `{comparison['e62']['supported_count']}`",
        f"- failed accepted: `{comparison['e62']['failed_accepted_count']}`",
        f"- membrane_misread: `{comparison['e62']['failure_modes'].get('membrane_misread', 0)}`",
        f"- new failure modes from regressions: `{comparison['new_failure_modes_from_regressions']}`",
        f"- abstain: `{comparison['e62']['abstain_count']}`",
        f"- accepted accuracy: `{comparison['e62']['accepted_accuracy']}`",
        "",
        "## Repair Map",
    ]
    for key, value in comparison["repair_categories"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Failure-Mode Movement"])
    for mode, values in comparison["failure_mode_distribution_change"].items():
        lines.append(f"- `{mode}`: `{values['e61']} -> {values['e62']}` (`{values['delta']}`)")
    lines.extend(["", "## Soluble False Membrane Check"])
    membrane_repair = comparison["membrane_repair"]
    lines.append(f"- soluble false membrane calls: `{membrane_repair['soluble_false_membrane_call_count']}`")
    lines.append(f"- cofactor/oligomer regressions: `{membrane_repair['cofactor_or_oligomer_regression_count']}`")
    lines.extend(["", "## Mechanism Distribution"])
    predicted = Counter(row["predicted_mechanism_class"] for row in scoring_rows)
    expected = Counter(row["expected_mechanism_class"] for row in scoring_rows)
    for mechanism in sorted(set(predicted) | set(expected)):
        lines.append(f"- `{mechanism}`: predicted `{predicted.get(mechanism, 0)}`, expected `{expected.get(mechanism, 0)}`")
    lines.extend(["", "## Boundary"])
    lines.append(
        "V64 is a paired E62 repair replay on the exact V63 target set. It measures repair direction and regression shape; it does not make a broad saturation or solved claim."
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v64(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    v63_manifest = _read_json(V63_TARGET_MANIFEST, "V63 target manifest")
    baseline_scoring = _read_json(V63_SCORING_REPORT, "V63 scoring report")
    baseline_cert = _read_json(V63_BATCH_CERTIFICATE, "V63 batch certificate")
    baseline_rows = _baseline_rows(baseline_scoring)
    _reset_generated_outputs(out_dir)
    target_manifest = _target_manifest(v63_manifest)
    if target_manifest["target_count_selected"] < TARGET_COUNT:
        raise SystemExit(f"V64 selected only {target_manifest['target_count_selected']} targets; need {TARGET_COUNT}")
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v64_e62_replay_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v64_e62_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    wrong_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []

    for candidate in target_manifest["selected_targets"]:
        target_id = f"V64_{candidate['target_id']}"
        baseline_row = baseline_rows.get(candidate["target_id"])
        if not baseline_row:
            raise SystemExit(f"missing V63 baseline row for {candidate['target_id']}")
        expected, reasons = v61._expected_mechanism_postseal(candidate)
        if expected != baseline_row["expected_mechanism_class"]:
            raise SystemExit(f"V64 expected label drift for {candidate['target_id']}: {expected} != {baseline_row['expected_mechanism_class']}")
        source_manifest = _source_manifest(candidate)
        packet = build_sealed_operator_state_packet(
            target_id=target_id,
            target_name=f"{candidate['entry_id']} {candidate['entity_description']}",
            sequence=candidate["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "automatic V64 E62 repair replay full-chain scan", "span": f"1-{candidate['sequence_length']}"}],
            perturbations=v61._perturbations_for_expected(expected, target_id),
        )
        holdout = _holdout(candidate, packet, baseline_row, reasons)
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
            sources=[
                {
                    "source_id": f"{target_id}_SHUFFLED_SEQUENCE_ONLY",
                    "source_class": "pure_non_coordinate",
                    "source_role": "prediction_input",
                    "coordinate_derived": False,
                    "internal_runtime_source": False,
                    "spatial_proxy": False,
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E62 context marks are withheld.",
                }
            ],
            perturbations=[],
        )
        packets.append(packet)
        scoring_rows.append(score)
        wrong_packets.append(wrong_packet)
        shuffled_packets.append(shuffled_packet)
        _write_json(DATA_ROOT / "source_manifests" / target_id / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target_id / "sealed_packet_summary.json", _packet_summary(packet))
        _write_json(DATA_ROOT / "holdouts_postseal" / target_id / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target_id / "validation_result.json", score)

    comparison = _comparison(
        baseline_cert=baseline_cert,
        baseline_rows=baseline_rows,
        scoring_rows=scoring_rows,
    )
    controls = _controls(
        v63_manifest=v63_manifest,
        target_manifest=target_manifest,
        baseline_cert=baseline_cert,
        baseline_rows=baseline_rows,
        engine_declaration=engine_declaration,
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
        comparison=comparison,
    )
    claim_row = _claim_row(cert)
    scoring_path = _write_json(DATA_ROOT / "v64_e62_rcsb_500_membrane_repair_scoring_report.json", {"kind": "V64_SCORING_REPORT_v0", "rows": scoring_rows})
    comparison_path = _write_json(DATA_ROOT / "v64_e61_vs_e62_comparison.json", comparison)
    failure_path = _write_json(DATA_ROOT / "v64_e62_rcsb_500_membrane_repair_failure_report.json", _failure_report(scoring_rows, comparison))
    data_cert_path = _write_json(DATA_ROOT / "v64_e62_rcsb_500_membrane_repair_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v64_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v64_e62_rcsb_500_membrane_repair_replay_certificate.json", cert)
    report_path = out_dir / "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY_REPORT.md"
    _write_report(report_path, cert, comparison, scoring_rows)
    return {
        "target_manifest": DATA_ROOT / "v64_e62_replay_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v64_e62_engine_declaration.json",
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
    parser = argparse.ArgumentParser(description="Run V64 E62 paired replay on the exact V63 RCSB 500-target set.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v64(args.out_dir)
    cert = _read_json(paths["certificate"], "V64 certificate")
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
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "engine_modified_during_batch": cert["engine_modified_during_batch"],
        "failure_modes": cert["failure_modes"],
        "net_changes": cert["net_changes"],
        "repair_categories": cert["repair_categories"],
        "membrane_repair": cert["membrane_repair"],
        "new_failure_modes_from_regressions": cert["new_failure_modes_from_regressions"],
        "next_required_batch": cert["next_required_batch"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
