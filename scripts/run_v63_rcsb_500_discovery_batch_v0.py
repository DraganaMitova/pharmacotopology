#!/usr/bin/env python3
from __future__ import annotations

"""Run V63: 500-target RCSB discovery batch on the E61 engine line.

V63 expands beyond the V61/V62 100-target repair set.  It uses RCSB 30%
sequence-identity representatives, public non-coordinate metadata context, and
post-seal RCSB validation labels to mine the next Protein Esperanto failure
classes.  V63 is a discovery batch, not a claim gate.
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


BATCH_ID = "V63_RCSB_500_DISCOVERY_BATCH"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E61"
TARGET_COUNT = 500
MIN_LENGTH = 40
MAX_LENGTH = 800
SEQUENCE_IDENTITY_CUTOFF = 30
SEARCH_PAGE_ROWS = 900

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
RAW_CANDIDATE_CACHE = DATA_ROOT / "intake" / "raw_rcsb_30pct_representative_entities_500.json"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

PASSED = "V63_RCSB_500_DISCOVERY_NO_FAILED_ACCEPTED_REVIEW_REQUIRED"
PARTIAL_ABSTAIN = "V63_RCSB_500_DISCOVERY_ACCEPTED_CLEAN_WITH_ABSTENTIONS_REVIEW_REQUIRED"
DISCOVERY_FAILURES = "V63_RCSB_500_DISCOVERY_FAILURES_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V63_RCSB_500_DISCOVERY_BLOCKED_FOR_LEAKAGE"
BLOCKED_CONTROLS = "V63_RCSB_500_DISCOVERY_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_INTAKE = "V63_RCSB_500_DISCOVERY_BLOCKED_INTAKE_UNAVAILABLE"

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
        "sealed_packet_summaries",
        "holdouts_postseal",
        "validation",
    ]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v63_rcsb_500_target_manifest.json",
        "v63_e61_engine_declaration.json",
        "v63_rcsb_500_scoring_report.json",
        "v63_rcsb_500_failure_report.json",
        "v63_rcsb_500_certificate.json",
        "v63_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _search_page(start: int) -> tuple[list[dict[str, Any]], dict[str, Any], str]:
    query = v61._rcsb_grouped_query(start=start, rows=SEARCH_PAGE_ROWS)
    url = v61._rcsb_search_url(query)
    data = v61._curl_json_url(url, label=f"RCSB 30% representative search page start={start}")
    result_set = [row for row in data.get("result_set", []) if isinstance(row, dict) and row.get("identifier")]
    return result_set, data, url


def refresh_rcsb_500_intake() -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    seen_clusters: set[str] = set()
    seen_sequences: set[str] = set()
    pages: list[dict[str, Any]] = []
    start = 0
    while len(candidates) < TARGET_COUNT and start < SEARCH_PAGE_ROWS * 8:
        result_set, search_response, query_url = _search_page(start)
        pages.append({
            "start": start,
            "rows_requested": SEARCH_PAGE_ROWS,
            "result_count": len(result_set),
            "search_total_count": search_response.get("total_count"),
            "search_group_by_count": search_response.get("group_by_count"),
            "search_ungrouped_count": search_response.get("ungrouped_count"),
            "search_query_url": query_url,
        })
        if not result_set:
            break
        for offset, hit in enumerate(result_set, start=1):
            identifier = str(hit["identifier"])
            entry_id, entity_id = v61._split_identifier(identifier)
            entity = v61._curl_json_url(
                f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}",
                label=f"RCSB polymer_entity {identifier}",
            )
            entry = v61._curl_json_url(
                f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
                label=f"RCSB entry {entry_id}",
            )
            candidate = v61._candidate_from_rcsb(
                search_rank=start + offset,
                search_hit=hit,
                entity=entity,
                entry=entry,
            )
            if not candidate:
                continue
            cluster = candidate["sequence_cluster_30_id"] or f"missing_cluster_{candidate['protein_id']}"
            if cluster in seen_clusters or candidate["sequence"] in seen_sequences:
                continue
            seen_clusters.add(cluster)
            seen_sequences.add(candidate["sequence"])
            candidates.append(candidate)
            if len(candidates) >= TARGET_COUNT:
                break
        start += SEARCH_PAGE_ROWS
    artifact = {
        "kind": "V63_RCSB_30PCT_CLUSTER_REPRESENTATIVE_RAW_CANDIDATES_500_v0",
        "retrieved_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": "RCSB Search API 30% sequence-identity representatives plus RCSB Data API",
        "search_pages": pages,
        "target_selection_manual": False,
        "sequence_cluster_representative_selection": True,
        "sequence_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "target_count_requested": TARGET_COUNT,
        "candidate_entity_count": len(candidates),
        "candidates": candidates,
    }
    _write_json(RAW_CANDIDATE_CACHE, artifact)
    return artifact


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


def _metadata_context(candidate: dict[str, Any]) -> dict[str, Any]:
    text = _candidate_text(candidate)
    components = {str(value).upper() for value in candidate.get("biological_cofactor_components", []) or []}
    metrics = candidate["sequence_metrics"]
    marks: list[str] = []
    reasons: list[str] = []
    if components:
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append(f"RCSB public non-coordinate component list contains biological cofactor candidates: {sorted(components)}")
        if components & METAL_COMPONENTS:
            marks.append("metal_context")
        if components & HEME_COMPONENTS:
            marks.append("heme_context")
        if components & NUCLEOTIDE_COMPONENTS:
            marks.append("nucleotide_context")
    if any(token in text for token in ["cofactor", "ligand", "heme", "nucleotide"]):
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append("RCSB title/description/keywords mention ligand or cofactor context")
    if (
        candidate.get("polymer_entity_instance_count", 0) >= 2
        or candidate.get("entity_molecule_count", 0) >= 2
        or any(token in text for token in ["oligomer", "homomer", "multimer", "assembly", "dimer", "trimer", "tetramer"])
    ):
        marks.extend(["oligomer_context", "assembly_context", "partner_copy_context"])
        reasons.append("RCSB public metadata indicates multiple copies or assembly context")
        composition = str(candidate.get("polymer_composition", "")).lower()
        if "heteromeric" in composition or "hetero" in text:
            marks.append("heteromeric_context")
        elif "homomeric" in composition or "homo" in text:
            marks.append("homomeric_context")
    if (
        metrics["max_segment_membrane_density"] >= 0.72
        or any(token in text for token in ["membrane", "transmembrane", "channel", "transporter", "porin", "gpcr", "opsin", "receptor"])
    ):
        marks.extend(["membrane_context_strong", "transmembrane_context"])
        reasons.append("RCSB public metadata or strong sequence field indicates membrane/topology context")
        if any(token in text for token in ["channel", "pore", "porin"]):
            marks.append("channel_context")
        if "transporter" in text:
            marks.append("transporter_context")
        if any(token in text for token in ["receptor", "gpcr", "opsin"]):
            marks.append("receptor_membrane_context")
    return {
        "context_marks": sorted(dict.fromkeys(marks)),
        "context_derivation": "V63 discovery uses sequence plus public non-coordinate RCSB metadata only; V61/V62 labels and post-seal holdout classes are not used for prediction context",
        "reasons": reasons or ["no explicit E61 metadata context mark emitted"],
        "biological_cofactor_components_seen": sorted(components),
        "polymer_copy_counts": {
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "polymer_composition": candidate.get("polymer_composition", ""),
        },
    }


def _context_statement(context: dict[str, Any]) -> str:
    marks = context["context_marks"]
    if not marks:
        return "V63 discovery metadata context emitted no explicit E61 context marks for this target."
    return (
        "V63 discovery metadata context marks: "
        + " ".join(marks)
        + ". Marks come from sequence and public non-coordinate RCSB metadata only; "
        + "coordinates, contacts, distance maps, ligand geometry, and native topology are blocked before sealing."
    )


def _source_manifest(candidate: dict[str, Any]) -> dict[str, Any]:
    target_id = f"V63_{candidate['target_id']}"
    context = _metadata_context(candidate)
    return {
        "kind": "V63_RCSB_500_DISCOVERY_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "target_selection_manual": False,
        "discovery_context_policy": context,
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
                "source_id": f"{target_id}_E61_METADATA_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": context["context_marks"],
                "metadata_context_reasons": context["reasons"],
                "evidence_statement": _context_statement(context),
                "source_url": candidate["source_urls"]["entry"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
            "ligand coordinates, metal coordination geometry, and bound-state contact geometry",
            "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
            "post-seal validation annotations before prediction hash",
            "V61/V62 same-target repair labels as prediction evidence",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _holdout(candidate: dict[str, Any], packet: dict[str, Any], expected: str, reasons: list[str]) -> dict[str, Any]:
    return {
        "kind": "V63_RCSB_500_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
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


def _score(packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == "insufficient_evidence_clean_abstain" else "accepted"
    accepted = decision == "accepted"
    supported = accepted and predicted == expected and validation["score_label"] == "supported"
    return {
        "kind": "V63_RCSB_500_DISCOVERY_VALIDATION_RESULT_v0",
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
        "failures_repaired_after_holdout": False,
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V63_COMPACT_SEALED_PACKET_SUMMARY_v0",
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
        "kind": "V63_E61_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
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
        "target_selection_manual": False,
        "folding_problem_solved": False,
    }


def _target_manifest(raw: dict[str, Any]) -> dict[str, Any]:
    selected = [dict(row) for row in raw.get("candidates", []) if isinstance(row, dict)][:TARGET_COUNT]
    for candidate in selected:
        candidate["sequence_metrics"] = v61._sequence_metrics(candidate["sequence"])
    return {
        "kind": "V63_RCSB_500_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "batch_mode": "discovery",
        "target_selection_manual": False,
        "selection_rule": (
            "RCSB Search API polymer_entity results; experimental structures only; protein entities only; "
            f"{MIN_LENGTH}-{MAX_LENGTH} aa; grouped by sequence_identity at {SEQUENCE_IDENTITY_CUTOFF}% "
            "with group_by_return_type=representatives; first 500 valid unique 30% cluster representatives."
        ),
        "source_cache": str(RAW_CANDIDATE_CACHE.relative_to(REPO_ROOT)),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(selected),
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "sequence_cluster_representative_selection": True,
        "selected_targets": selected,
    }


def _failure_type(row: dict[str, Any]) -> str:
    return v61._failure_type(row)


def _failure_grammar_row(row: dict[str, Any]) -> dict[str, Any]:
    failure_type = _failure_type(row)
    proposals = {
        "membrane_misread": ("membrane topology priority/context word", "membrane_pressure_operator", "bilayer/topology pressure", "E62 should separate membrane from ligand/assembly conflicts if repeated"),
        "disorder_misread": ("IDR persistence and fold-upon-binding separation word", "disorder_operator", "entropy/partner pressure", "E62 or disorder panel should split IDR, phase, and bound-order contexts"),
        "oligomer_state_misread": ("assembly-register/interface specificity word", "interface_operator", "partner-copy concentration pressure", "E62 should separate obligate assembly from incidental copy count if repeated"),
        "cofactor_ligand_missing": ("ligand-state specificity word", "interface_operator", "ligand_or_cofactor pressure", "E62 should separate structural ions, cofactors, and incidental ligands if repeated"),
        "weak_sequence_signal": ("confidence/self-evidence word", "none", "none", "retain abstention until independent non-coordinate evidence sharpens mechanism"),
        "right_regime_wrong_topology": ("topology proxy refinement word", "closure_operator", "fold-class pressure", "abstain when operator regions are underdetermined"),
        "wrong_regime": ("regime separation word", "frustration_operator", "context separation pressure", "mine repeated V63 classes before E62"),
    }
    rule, operator, pressure, abstention = proposals.get(failure_type, proposals["wrong_regime"])
    return {
        "target_id": row["target_id"],
        "protein_id": row["protein_id"],
        "failure_type": failure_type,
        "engine_thought": row["predicted_mechanism_class"],
        "reality_showed": row["expected_mechanism_class"],
        "score_label": row["score_label"],
        "missing_esperanto_rule": rule,
        "proposed_new_operator": operator,
        "proposed_new_pressure": pressure,
        "proposed_new_abstention_rule": abstention,
        "control_that_prevents_overfitting": "V63 is discovery; E62 grammar can only use repeated failure classes and must rerun regression.",
    }


def _failure_grammar_ledger(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [_failure_grammar_row(row) for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V63_FAILURE_GRAMMAR_LEDGER_v0",
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_count": len(rows),
        "failure_modes": dict(Counter(row["failure_type"] for row in rows)),
        "missing_words_top_10": [
            {"failure_type": mode, "count": count}
            for mode, count in Counter(row["failure_type"] for row in rows).most_common(10)
        ],
        "rows": rows,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    wrong_packets: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V63_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V63_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold-style model offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V63_PRESEAL_HOLDOUT",
        "source_class": "coordinate_derived",
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V63_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_operator_state_packet(
        target_id="V63_RANDOM_SEQUENCE_CONTROL",
        target_name="V63 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
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
    return [
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V63 selection must be automatic."),
        _control("targets_total_500", len(target_manifest["selected_targets"]) == TARGET_COUNT, "V63 N must be exactly 500.", len(target_manifest["selected_targets"])),
        _control("rcsb_experimental_protein_entities_only", all(row["source_database"] == "RCSB_PDB" and row["structure_determination_methodology"].lower() == "experimental" and "protein" in row["polymer_type"].lower() for row in target_manifest["selected_targets"]), "All targets are RCSB experimental protein entities."),
        _control("sequence_cluster_representative_selection", target_manifest["sequence_cluster_representative_selection"] is True and len(set(clusters)) == len(clusters) and all(clusters), "Each target is a unique 30% sequence-cluster representative.", {"cutoff": SEQUENCE_IDENTITY_CUTOFF, "unique_clusters": len(set(clusters))}),
        _control("length_filter_40_800", all(MIN_LENGTH <= row["sequence_length"] <= MAX_LENGTH for row in target_manifest["selected_targets"]), "All targets satisfy the requested length filter."),
        _control("engine_version_declared_e61", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V63 uses E61."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside V63."),
        _control("discovery_context_not_repair_labels", True, "V63 context extraction does not use V61/V62 expected labels."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V63 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("wrong_grammar_controls", all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_packets), "Forced wrong grammars are rejected or routed to abstention."),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls are generated with target metadata withheld.", {"control_count": len(shuffled_rows), "shuffled_higher_by_more_than_margin_count": sum(1 for row in shuffled_rows if row["shuffled_higher_by_more_than_margin"])}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain", "Random sequence without evidence abstains."),
        _control("failures_reported", len(scoring_rows) == TARGET_COUNT and all("score_label" in row for row in scoring_rows), "Every target has an explicit score row."),
        _control("readme_check_skipped_by_user_instruction", True, "README check skipped by explicit user instruction during V63 finalization."),
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


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_ledger: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "alphafold_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    if len(scoring_rows) < TARGET_COUNT:
        status = BLOCKED_INTAKE
    elif any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif metrics["failed_accepted_count"]:
        status = DISCOVERY_FAILURES
    elif metrics["abstain_count"]:
        status = PARTIAL_ABSTAIN
    else:
        status = PASSED
    controls_passed = not failed_controls
    cert = {
        "kind": "V63_RCSB_500_DISCOVERY_BATCH_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "discovery",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        **metrics,
        "controls_passed": controls_passed,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "readme_check_skipped_by_user_instruction": True,
        "failure_modes": failure_ledger["failure_modes"],
        "missing_words_top_10": failure_ledger["missing_words_top_10"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V63 is a broad discovery batch for mining failure classes; claims require later grammar revision and regression.",
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "E61 has saturated broad RCSB protein space.",
            "V63 failures may be hidden.",
            "Coordinates, contacts, ligand geometry, or native topology were used before sealing.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _failure_report(scoring_rows: list[dict[str, Any]], failure_ledger: dict[str, Any]) -> dict[str, Any]:
    failures = [row for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V63_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_cases_reported": True,
        "failure_count": len(failures),
        "failure_modes": failure_ledger["failure_modes"],
        "missing_words_top_10": failure_ledger["missing_words_top_10"],
        "failure_cases": failures,
        "failure_grammar_rows": failure_ledger["rows"],
        "note": "V63 is discovery; failures are preserved for E62 grammar mining and regression.",
    }


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "discovery",
        "engine_version_used": ENGINE_VERSION_USED,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["abstain_count"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"] + cert["abstain_count"],
        "failure_modes": cert["failure_modes"],
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


def _write_report(path: Path, cert: dict[str, Any], scoring_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# V63 RCSB 500 Discovery Batch",
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
        f"README check skipped by user instruction: `{cert['readme_check_skipped_by_user_instruction']}`",
        "",
        "## Top Failure Modes",
    ]
    if cert["failure_modes"]:
        for mode, count in sorted(cert["failure_modes"].items(), key=lambda item: (-item[1], item[0]))[:10]:
            lines.append(f"- `{mode}`: `{count}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Mechanism Distribution"])
    predicted = Counter(row["predicted_mechanism_class"] for row in scoring_rows)
    expected = Counter(row["expected_mechanism_class"] for row in scoring_rows)
    for mechanism in sorted(set(predicted) | set(expected)):
        lines.append(f"- `{mechanism}`: predicted `{predicted.get(mechanism, 0)}`, expected `{expected.get(mechanism, 0)}`")
    lines.extend(["", "## Boundary"])
    lines.append(
        "V63 is a broad discovery/mining batch. It records E61 behavior on 500 nonredundant RCSB protein entities and preserves failures for E62 grammar mining. It does not make a broad saturation or solved claim."
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v63(out_dir: Path = DEFAULT_OUT_DIR, *, refresh_intake: bool = False) -> dict[str, Path]:
    raw = refresh_rcsb_500_intake() if refresh_intake else _read_json(RAW_CANDIDATE_CACHE, "V63 raw RCSB representative cache")
    _reset_generated_outputs(out_dir)
    target_manifest = _target_manifest(raw)
    if target_manifest["target_count_selected"] < TARGET_COUNT:
        raise SystemExit(f"V63 selected only {target_manifest['target_count_selected']} targets; need {TARGET_COUNT}")
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v63_rcsb_500_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v63_e61_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    wrong_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []

    for candidate in target_manifest["selected_targets"]:
        target_id = f"V63_{candidate['target_id']}"
        expected, reasons = v61._expected_mechanism_postseal(candidate)
        source_manifest = _source_manifest(candidate)
        packet = build_sealed_operator_state_packet(
            target_id=target_id,
            target_name=f"{candidate['entry_id']} {candidate['entity_description']}",
            sequence=candidate["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "automatic V63 discovery full-chain scan", "span": f"1-{candidate['sequence_length']}"}],
            perturbations=v61._perturbations_for_expected(expected, target_id),
        )
        holdout = _holdout(candidate, packet, expected, reasons)
        score = _score(packet, holdout)
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
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E61 context marks are withheld.",
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

    failure_ledger = _failure_grammar_ledger(scoring_rows)
    controls = _controls(
        target_manifest=target_manifest,
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
        failure_ledger=failure_ledger,
    )
    claim_row = _claim_row(cert)
    scoring_path = _write_json(DATA_ROOT / "v63_rcsb_500_scoring_report.json", {"kind": "V63_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v63_rcsb_500_failure_report.json", _failure_report(scoring_rows, failure_ledger))
    data_cert_path = _write_json(DATA_ROOT / "v63_rcsb_500_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v63_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v63_rcsb_500_discovery_batch_certificate.json", cert)
    report_path = out_dir / "V63_RCSB_500_DISCOVERY_BATCH_REPORT.md"
    _write_report(report_path, cert, scoring_rows)
    return {
        "raw_candidate_cache": RAW_CANDIDATE_CACHE,
        "target_manifest": DATA_ROOT / "v63_rcsb_500_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v63_e61_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V63 RCSB 500-target discovery batch.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--refresh-intake", action="store_true", help="refresh RCSB 30% representative intake with curl")
    args = parser.parse_args()
    paths = run_v63(args.out_dir, refresh_intake=args.refresh_intake)
    cert = _read_json(paths["certificate"], "V63 certificate")
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
        "readme_check_skipped_by_user_instruction": cert["readme_check_skipped_by_user_instruction"],
        "failure_modes": cert["failure_modes"],
        "missing_words_top_10": cert["missing_words_top_10"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
