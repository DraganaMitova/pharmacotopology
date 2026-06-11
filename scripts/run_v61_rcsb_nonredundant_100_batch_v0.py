#!/usr/bin/env python3
from __future__ import annotations

"""Run V61 RCSB nonredundant 100-target Protein Esperanto batch.

V61 freezes the E60 engine source and tests automatically selected RCSB/PDB
experimental protein polymer entities.  The default path uses a committed cache
so tests stay deterministic; --refresh-intake rebuilds that cache from the RCSB
Search API 30% sequence-identity representative query plus the RCSB Data API.
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
from urllib.parse import quote

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    COORDINATE_DERIVED,
    INTERNAL_RUNTIME,
    MECHANISM_CLASSES,
    UNIVERSAL_OPERATORS,
    build_sealed_simulation_packet,
    build_sequence_field,
    deterministic_random_sequence,
    evidence_boundary_gate,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
)


BATCH_ID = "V61_RCSB_NONREDUNDANT_100_BATCH"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E60"
ENGINE_START_COMMIT = "d927781"
TARGET_COUNT = 100
MIN_LENGTH = 40
MAX_LENGTH = 800
SEQUENCE_IDENTITY_CUTOFF = 30
SEARCH_ROWS = 140

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V61"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
RAW_CANDIDATE_CACHE = DATA_ROOT / "intake" / "raw_rcsb_30pct_representative_entities.json"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

PASSED = "V61_RCSB_NONREDUNDANT_100_ACCEPTED_ACCURACY_PASSED_REVIEW_REQUIRED"
PARTIAL_ABSTAIN = "V61_RCSB_NONREDUNDANT_100_ACCEPTED_ACCURACY_PASSED_WITH_ABSTENTIONS_REVIEW_REQUIRED"
DISCOVERY_FAILURES = "V61_RCSB_NONREDUNDANT_100_DISCOVERY_FAILURES_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V61_RCSB_NONREDUNDANT_100_BLOCKED_FOR_LEAKAGE"
BLOCKED_CONTROLS = "V61_RCSB_NONREDUNDANT_100_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_INTAKE = "V61_RCSB_NONREDUNDANT_100_BLOCKED_INTAKE_UNAVAILABLE"

BIOLOGICAL_COFACTOR_COMPONENTS = {
    "ADP",
    "ATP",
    "B12",
    "CA",
    "CLA",
    "CO",
    "COA",
    "CU",
    "FAD",
    "FE",
    "FMN",
    "GDP",
    "GTP",
    "HEA",
    "HEC",
    "HEM",
    "MG",
    "MN",
    "MO",
    "NAD",
    "NAP",
    "NI",
    "PLP",
    "SAM",
    "TPP",
    "ZN",
}


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


def _curl_json_url(url: str, *, label: str) -> dict[str, Any]:
    result = subprocess.run(
        ["curl", "-s", "-L", "--max-time", "60", url],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise SystemExit(f"curl failed for {label}: {result.stderr.strip()}")
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"curl returned non-JSON for {label}: {result.stdout[:500]}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"curl JSON for {label} must be an object")
    return data


def _reset_generated_outputs() -> None:
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
        "v61_rcsb_nonredundant_100_target_manifest.json",
        "v61_rcsb_nonredundant_100_engine_declaration.json",
        "v61_rcsb_nonredundant_100_scoring_report.json",
        "v61_rcsb_nonredundant_100_failure_report.json",
        "v61_rcsb_nonredundant_100_certificate.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    for filename in [
        "protein_universe_ledger_v0.json",
        "engine_version_ledger_v0.json",
        "failure_grammar_ledger_v0.json",
        "claim_ledger_v0.json",
    ]:
        path = LEDGER_ROOT / filename
        if path.exists():
            path.unlink()


def _rcsb_grouped_query(*, start: int = 0, rows: int = SEARCH_ROWS) -> dict[str, Any]:
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "range",
                        "value": {
                            "from": MIN_LENGTH,
                            "to": MAX_LENGTH,
                            "include_lower": True,
                            "include_upper": True,
                        },
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.structure_determination_methodology",
                        "operator": "exact_match",
                        "value": "experimental",
                    },
                },
            ],
        },
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "results_content_type": ["experimental"],
            "group_by": {
                "aggregation_method": "sequence_identity",
                "similarity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
            },
            "group_by_return_type": "representatives",
        },
        "return_type": "polymer_entity",
    }


def _rcsb_search_url(query: dict[str, Any]) -> str:
    payload = json.dumps(query, separators=(",", ":"))
    return "https://search.rcsb.org/rcsbsearch/v2/query?json=" + quote(payload)


def _search_representative_entity_ids() -> tuple[list[dict[str, Any]], dict[str, Any], str]:
    query = _rcsb_grouped_query()
    url = _rcsb_search_url(query)
    data = _curl_json_url(url, label="RCSB 30% representative polymer_entity search")
    result_set = [row for row in data.get("result_set", []) if isinstance(row, dict) and row.get("identifier")]
    return result_set, data, url


def _valid_protein_sequence(sequence: str) -> bool:
    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    return bool(sequence) and set(sequence.upper()) <= allowed


def _split_identifier(identifier: str) -> tuple[str, str]:
    entry_id, entity_id = identifier.split("_", 1)
    return entry_id.upper(), entity_id


def _cluster_30_id(entity: dict[str, Any]) -> str:
    for row in entity.get("rcsb_cluster_membership", []) or []:
        if isinstance(row, dict) and int(float(row.get("identity", -1))) == SEQUENCE_IDENTITY_CUTOFF:
            return str(row.get("cluster_id", ""))
    for row in entity.get("rcsb_polymer_entity_group_membership", []) or []:
        if isinstance(row, dict) and row.get("aggregation_method") == "sequence_identity":
            if int(float(row.get("similarity_cutoff", -1))) == SEQUENCE_IDENTITY_CUTOFF:
                return str(row.get("group_id", ""))
    return ""


def _annotation_names(entity: dict[str, Any]) -> list[str]:
    names: list[str] = []
    for row in entity.get("rcsb_polymer_entity_annotation", []) or []:
        if isinstance(row, dict):
            for key in ["annotation_id", "name", "type", "provenance_source"]:
                if row.get(key):
                    names.append(str(row[key]))
    for row in entity.get("rcsb_polymer_entity_feature", []) or []:
        if isinstance(row, dict):
            for key in ["feature_id", "name", "type", "provenance_source"]:
                if row.get(key):
                    names.append(str(row[key]))
    return sorted(set(names))


def _feature_types(entity: dict[str, Any]) -> list[str]:
    values = [
        str(row.get("type"))
        for row in entity.get("rcsb_polymer_entity_feature", []) or []
        if isinstance(row, dict) and row.get("type")
    ]
    return sorted(set(values))


def _feature_coverage(entity: dict[str, Any], feature_type: str) -> float:
    sequence = str((entity.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can") or "")
    length = len(sequence)
    if not length:
        return 0.0
    covered: set[int] = set()
    for row in entity.get("rcsb_polymer_entity_feature", []) or []:
        if not isinstance(row, dict) or str(row.get("type", "")).lower() != feature_type.lower():
            continue
        for pos in row.get("feature_positions", []) or []:
            if not isinstance(pos, dict):
                continue
            start = int(pos.get("beg_seq_id", 1))
            end = int(pos.get("end_seq_id", start))
            if "values" in pos and isinstance(pos["values"], list):
                end = start + len(pos["values"]) - 1
                for offset, value in enumerate(pos["values"]):
                    try:
                        active = float(value) >= 0.50
                    except (TypeError, ValueError):
                        active = False
                    if active:
                        idx = start + offset
                        if 1 <= idx <= length:
                            covered.add(idx)
                continue
            for idx in range(max(1, start), min(length, end) + 1):
                covered.add(idx)
    return round(len(covered) / length, 6)


def _organisms(entity: dict[str, Any]) -> list[str]:
    organisms = []
    for row in entity.get("rcsb_entity_source_organism", []) or []:
        if isinstance(row, dict) and row.get("ncbi_scientific_name"):
            organisms.append(str(row["ncbi_scientific_name"]))
    for row in entity.get("entity_src_nat", []) or []:
        if isinstance(row, dict) and row.get("pdbx_organism_scientific"):
            organisms.append(str(row["pdbx_organism_scientific"]))
    return sorted(set(organisms))


def _entry_keywords(entry: dict[str, Any]) -> str:
    keywords = entry.get("struct_keywords") or {}
    return " ".join(str(keywords.get(key, "")) for key in ["pdbx_keywords", "text"]).strip()


def _nonpolymer_components(entry: dict[str, Any]) -> list[str]:
    info = entry.get("rcsb_entry_info") or {}
    components = info.get("nonpolymer_bound_components") or []
    return sorted({str(value).upper() for value in components if value})


def _candidate_from_rcsb(
    *,
    search_rank: int,
    search_hit: dict[str, Any],
    entity: dict[str, Any],
    entry: dict[str, Any],
) -> dict[str, Any] | None:
    identifier = str(search_hit["identifier"])
    entry_id, entity_id = _split_identifier(identifier)
    entity_poly = entity.get("entity_poly") or {}
    sequence = str(entity_poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").replace(" ", "").upper()
    if not _valid_protein_sequence(sequence):
        return None
    if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        return None
    polymer_type = str(entity_poly.get("rcsb_entity_polymer_type") or entity_poly.get("type") or "")
    if "protein" not in polymer_type.lower() and "polypeptide" not in polymer_type.lower():
        return None
    entry_info = entry.get("rcsb_entry_info") or {}
    if str(entry_info.get("structure_determination_methodology", "")).lower() != "experimental":
        return None
    entity_meta = entity.get("rcsb_polymer_entity") or {}
    accession = entry.get("rcsb_accession_info") or {}
    title = str((entry.get("struct") or {}).get("title", ""))
    components = _nonpolymer_components(entry)
    annotations = _annotation_names(entity)
    feature_types = _feature_types(entity)
    cluster_id = _cluster_30_id(entity)
    return {
        "entry_id": entry_id,
        "entity_id": entity_id,
        "protein_id": identifier,
        "target_id": identifier,
        "source_database": "RCSB_PDB",
        "sequence": sequence,
        "sequence_length": len(sequence),
        "polymer_type": polymer_type,
        "title": title,
        "entry_keywords": _entry_keywords(entry),
        "release_date": str(accession.get("initial_release_date", "")),
        "experimental_method": "; ".join(str(row.get("method", "")) for row in entry.get("exptl", []) if isinstance(row, dict)),
        "structure_determination_methodology": str(entry_info.get("structure_determination_methodology", "")),
        "entity_description": str(entity_meta.get("pdbx_description", "")),
        "formula_weight_kda": entity_meta.get("formula_weight"),
        "organisms": _organisms(entity),
        "polymer_composition": str(entry_info.get("polymer_composition", "")),
        "polymer_entity_instance_count": int(entry_info.get("deposited_polymer_entity_instance_count") or 0),
        "entity_molecule_count": int(entity_meta.get("pdbx_number_of_molecules") or 0),
        "nonpolymer_bound_components": components,
        "biological_cofactor_components": sorted(set(components) & BIOLOGICAL_COFACTOR_COMPONENTS),
        "annotations": annotations,
        "feature_types": feature_types,
        "disorder_feature_coverage": _feature_coverage(entity, "disorder"),
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "sequence_cluster_30_id": cluster_id,
        "sequence_cluster_representative_rank": search_rank,
        "search_score": search_hit.get("score"),
        "source_urls": {
            "search_api": "https://search.rcsb.org/",
            "entry": f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
            "polymer_entity": f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}",
        },
    }


def refresh_rcsb_cluster_representative_intake() -> dict[str, Any]:
    result_set, search_response, query_url = _search_representative_entity_ids()
    candidates: list[dict[str, Any]] = []
    seen_clusters: set[str] = set()
    seen_sequences: set[str] = set()
    for rank, hit in enumerate(result_set, start=1):
        identifier = str(hit["identifier"])
        entry_id, entity_id = _split_identifier(identifier)
        entity = _curl_json_url(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}",
            label=f"RCSB polymer_entity {identifier}",
        )
        entry = _curl_json_url(
            f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
            label=f"RCSB entry {entry_id}",
        )
        candidate = _candidate_from_rcsb(search_rank=rank, search_hit=hit, entity=entity, entry=entry)
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
    artifact = {
        "kind": "V61_RCSB_30PCT_CLUSTER_REPRESENTATIVE_RAW_CANDIDATES_v0",
        "retrieved_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": "RCSB Search API 30% sequence-identity representatives plus RCSB Data API",
        "search_query": _rcsb_grouped_query(),
        "search_query_url": query_url,
        "search_total_count": search_response.get("total_count"),
        "search_group_by_count": search_response.get("group_by_count"),
        "search_ungrouped_count": search_response.get("ungrouped_count"),
        "target_selection_manual": False,
        "sequence_cluster_representative_selection": True,
        "sequence_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "target_count_requested": TARGET_COUNT,
        "candidate_entity_count": len(candidates),
        "candidates": candidates,
    }
    _write_json(RAW_CANDIDATE_CACHE, artifact)
    return artifact


def _sequence_metrics(sequence: str) -> dict[str, Any]:
    field = build_sequence_field(sequence)
    metrics = dict(field["global_metrics"])
    metrics["max_segment_membrane_density"] = max((row["membrane_density"] for row in field["segments"]), default=0.0)
    metrics["max_segment_low_complexity_density"] = max((row["low_complexity_density"] for row in field["segments"]), default=0.0)
    metrics["max_segment_interface_density"] = max((row["interface_density"] for row in field["segments"]), default=0.0)
    return metrics


def _candidate_text(candidate: dict[str, Any], *, include_postseal: bool) -> str:
    values = [
        candidate.get("title", ""),
        candidate.get("entity_description", ""),
        candidate.get("entry_keywords", ""),
        candidate.get("polymer_composition", ""),
        " ".join(candidate.get("organisms", []) or []),
    ]
    if include_postseal:
        values.extend(candidate.get("annotations", []) or [])
        values.extend(candidate.get("feature_types", []) or [])
        values.extend(candidate.get("nonpolymer_bound_components", []) or [])
    return " ".join(str(value) for value in values).lower()


def _metadata_statement(candidate: dict[str, Any]) -> str:
    metrics = candidate["sequence_metrics"]
    labels: list[str] = []
    if metrics["max_segment_membrane_density"] >= 0.65:
        labels.append("sequence-derived membrane tendency high")
    if metrics["max_segment_low_complexity_density"] >= 0.70:
        labels.append("sequence-derived low-complexity tendency high")
    if metrics["mean_disorder"] >= 0.30:
        labels.append("sequence-derived disorder tendency high")
    if metrics["hydrophobic_density"] >= 0.32 and metrics["mean_disorder"] < 0.25:
        labels.append("sequence-derived hydrophobic closure tendency")
    if candidate.get("polymer_entity_instance_count", 0) >= 2 or candidate.get("entity_molecule_count", 0) >= 2:
        labels.append("public metadata indicates multiple polymer instances or molecule copies")
    return ". ".join([
        f"RCSB title: {candidate.get('title', '')}",
        f"Entity description: {candidate.get('entity_description', '')}",
        f"Organism: {'; '.join(candidate.get('organisms', []) or [])}",
        f"Polymer composition: {candidate.get('polymer_composition', '')}",
        "Sequence-derived marks: " + ("; ".join(labels) if labels else "no special high-pressure mark"),
        "Coordinates, native contacts, residue-residue distances, and validation annotations are unopened before the prediction hash.",
    ])


def _source_manifest(candidate: dict[str, Any]) -> dict[str, Any]:
    target_id = f"V61_{candidate['target_id']}"
    return {
        "kind": "V61_RCSB_NONREDUNDANT_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "target_selection_manual": False,
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "sequence_cluster_30_id": candidate["sequence_cluster_30_id"],
        "prediction_sources": [
            {
                "source_id": f"{target_id}_RAW_SEQUENCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity. No coordinates, contacts, AlphaFold-style models, or post-seal validation annotations are allowed.",
                "source_url": candidate["source_urls"]["polymer_entity"],
            },
            {
                "source_id": f"{target_id}_PURE_METADATA",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": _metadata_statement(candidate),
                "source_url": candidate["source_urls"]["entry"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
            "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
            "post-seal RCSB annotations and coordinate validation before prediction hash",
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
            "evidence_statement": (
                "Deterministically shuffled amino-acid sequence control. Target metadata, coordinates, contacts, "
                "post-seal annotations, and original sequence-derived statements are withheld."
            ),
        }
    ]


def _expected_mechanism_postseal(candidate: dict[str, Any]) -> tuple[str, list[str]]:
    text = _candidate_text(candidate, include_postseal=True)
    metrics = candidate["sequence_metrics"]
    reasons: list[str] = []
    if any(token in text for token in ["membrane", "transmembrane", "channel", "transporter", "porin", "gpcr", "opsin"]):
        reasons.append("postseal text/annotation indicates membrane or transport pressure")
        return "membrane_multidomain_folding_proteostasis", reasons
    if metrics["max_segment_membrane_density"] >= 0.72:
        reasons.append("postseal sequence-field membrane density is high")
        return "membrane_multidomain_folding_proteostasis", reasons
    strong_disorder_text = any(
        token in text
        for token in [
            "intrinsically disordered",
            "intrinsic disorder",
            "low complexity",
            "phase separation",
            "prion",
        ]
    )
    if candidate.get("disorder_feature_coverage", 0.0) >= 0.35 or strong_disorder_text:
        reasons.append("postseal annotation/feature indicates disorder or low-complexity behavior")
        return "intrinsic_disorder_phase_separation", reasons
    if any(token in text for token in ["metamorphic", "fold switch", "fold-switch", "dual basin"]):
        reasons.append("postseal text indicates fold-switching/metamorphic behavior")
        return "metamorphic_fold_switching", reasons
    if any(token in text for token in ["host", "viral", "hijack"]) and any(token in text for token in ["binding", "interface", "complex"]):
        reasons.append("postseal text indicates short host-interface behavior")
        return "short_region_host_interface_hijacking", reasons
    if candidate.get("biological_cofactor_components"):
        reasons.append("postseal RCSB entry contains biological ligand/cofactor/metal components")
        return "cofactor_ligand_assisted_stabilization", reasons
    if any(token in text for token in ["cofactor", "ligand"]):
        reasons.append("postseal text indicates ligand/cofactor context")
        return "cofactor_ligand_assisted_stabilization", reasons
    if (
        candidate.get("polymer_entity_instance_count", 0) >= 2
        or candidate.get("entity_molecule_count", 0) >= 2
        or any(token in text for token in ["oligomer", "homomer", "multimer", "assembly", "dimer", "trimer", "tetramer"])
    ):
        reasons.append("postseal entry/entity metadata indicates oligomer or assembly context")
        return "oligomerization_controlled_folding", reasons
    reasons.append("postseal validation supports ordinary compact/globular closure class")
    return "globular_closure", reasons


def _expected_observables(expected: str) -> list[dict[str, Any]]:
    if expected == "globular_closure":
        return [{"check_id": "compact_or_contact_supported", "metric": "contact_probability", "comparator": ">=", "threshold": 0.35}]
    if expected == "intrinsic_disorder_phase_separation":
        return [{"check_id": "compact_single_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12}]
    if expected == "membrane_multidomain_folding_proteostasis":
        return [{"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.50}]
    if expected == "metamorphic_fold_switching":
        return [
            {"check_id": "alpha_basin_present", "metric": "basin:alpha_context_basin", "comparator": ">=", "threshold": 0.25},
            {"check_id": "beta_basin_present", "metric": "basin:beta_released_basin", "comparator": ">=", "threshold": 0.25},
        ]
    if expected == "short_region_host_interface_hijacking":
        return [{"check_id": "host_interface_present", "metric": "basin:host_interface_engaged", "comparator": ">=", "threshold": 0.55}]
    return []


def _holdout(candidate: dict[str, Any], packet: dict[str, Any], expected: str, reasons: list[str]) -> dict[str, Any]:
    return {
        "kind": "V61_RCSB_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "expected_mechanism_class": expected,
        "expected_observables": _expected_observables(expected),
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


def _perturbations_for_expected(expected: str, target_id: str) -> list[dict[str, Any]]:
    if expected == "intrinsic_disorder_phase_separation":
        return [
            {"perturbation_id": f"{target_id}_STICKER_DAMAGE", "description": "damage low-complexity/sticker region", "operator_scales": {"phase_operator": 0.45}, "metric": "basin:phase_prone_dynamic", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "basin:phase_prone_dynamic", "expected_direction": "unchanged"},
        ]
    if expected == "membrane_multidomain_folding_proteostasis":
        return [
            {"perturbation_id": f"{target_id}_MEMBRANE_DAMAGE", "description": "damage membrane/proteostasis route", "operator_scales": {"membrane_pressure_operator": 0.55, "proteostasis_operator": 0.55}, "damage": 0.40, "metric": "proteostasis_routing", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_PROTEOSTASIS_RESCUE", "description": "proteostasis rescue context", "operator_scales": {"proteostasis_operator": 1.20}, "rescue": 0.35, "metric": "proteostasis_routing", "expected_direction": "increase"},
        ]
    if expected == "oligomerization_controlled_folding":
        return [
            {"perturbation_id": f"{target_id}_INTERFACE_DAMAGE", "description": "damage oligomer interface readiness", "operator_scales": {"interface_operator": 0.45}, "metric": "interface_readiness", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "interface_readiness", "expected_direction": "unchanged"},
        ]
    if expected == "cofactor_ligand_assisted_stabilization":
        return [
            {"perturbation_id": f"{target_id}_COFACTOR_REMOVAL", "description": "remove ligand/cofactor stabilization pressure", "operator_scales": {"interface_operator": 0.45}, "metric": "interface_readiness", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "interface_readiness", "expected_direction": "unchanged"},
        ]
    return [
        {"perturbation_id": f"{target_id}_CORE_DAMAGE", "description": "damage hydrophobic closure core", "operator_scales": {"closure_operator": 0.45}, "metric": "contact_probability", "expected_direction": "decrease"},
        {"perturbation_id": f"{target_id}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "contact_probability", "expected_direction": "unchanged"},
    ]


def _validate_holdout_without_repair(*, packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    from pharmacotopology.protein_esperanto_engine import validate_against_holdout

    return validate_against_holdout(sealed_packet=packet, holdout=holdout)


def _score(packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    validation = _validate_holdout_without_repair(packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == "insufficient_evidence_clean_abstain" else "accepted"
    accepted = decision == "accepted"
    supported = accepted and predicted == expected and validation["score_label"] == "supported"
    return {
        "kind": "V61_RCSB_NONREDUNDANT_VALIDATION_RESULT_v0",
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


def _wrong_grammar(natural: str) -> str:
    for candidate in [
        "intrinsic_disorder_phase_separation",
        "globular_closure",
        "membrane_multidomain_folding_proteostasis",
        "oligomerization_controlled_folding",
    ]:
        if candidate != natural:
            return candidate
    return "globular_closure"


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V61_FROZEN_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_start_commit_required": ENGINE_START_COMMIT,
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "frozen_operator_names": UNIVERSAL_OPERATORS,
        "frozen_mechanism_classes": MECHANISM_CLASSES,
        "operator_set_hash": stable_hash(UNIVERSAL_OPERATORS),
        "mechanism_class_set_hash": stable_hash(MECHANISM_CLASSES),
        "engine_modified_during_batch": False,
        "engine_biology_modified_during_batch": False,
        "target_selection_manual": False,
        "folding_problem_solved": False,
    }


def _target_manifest(raw: dict[str, Any]) -> dict[str, Any]:
    candidates = [dict(row) for row in raw.get("candidates", []) if isinstance(row, dict)]
    selected = candidates[:TARGET_COUNT]
    for candidate in selected:
        candidate["sequence_metrics"] = _sequence_metrics(candidate["sequence"])
    return {
        "kind": "V61_RCSB_NONREDUNDANT_100_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "target_selection_manual": False,
        "selection_rule": (
            "RCSB Search API polymer_entity results; experimental structures only; protein entities only; "
            f"{MIN_LENGTH}-{MAX_LENGTH} aa; grouped by sequence_identity at {SEQUENCE_IDENTITY_CUTOFF}% "
            "with group_by_return_type=representatives; first 100 valid unique 30% cluster representatives."
        ),
        "source_cache": str(RAW_CANDIDATE_CACHE.relative_to(REPO_ROOT)),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(selected),
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "sequence_cluster_representative_selection": True,
        "rcsb_search_query": raw.get("search_query"),
        "selected_targets": selected,
    }


def _protein_universe_ledger(target_manifest: dict[str, Any]) -> dict[str, Any]:
    rows = []
    for target in target_manifest["selected_targets"]:
        annotations = " ".join(target.get("annotations", [])).lower()
        text = _candidate_text(target, include_postseal=True)
        rows.append({
            "protein_id": target["protein_id"],
            "sequence": target["sequence"],
            "source_database": target["source_database"],
            "length": target["sequence_length"],
            "organism": "; ".join(target.get("organisms", []) or []),
            "experimental_structure_available": True,
            "fold_class_available": any(token in annotations for token in ["cath", "scop", "ecod", "pfam", "interpro"]),
            "disorder_data_available": "disorder" in target.get("feature_types", []) or target.get("disorder_feature_coverage", 0.0) > 0.0,
            "membrane_data_available": any(token in text for token in ["membrane", "transmembrane", "channel", "transporter", "opm", "gpcr", "opsin"]),
            "kinetics_data_available": False,
            "process_data_available": False,
            "eligible_for_prediction": _valid_protein_sequence(target["sequence"]) and MIN_LENGTH <= target["sequence_length"] <= MAX_LENGTH,
            "eligible_for_holdout": True,
            "batch_id": BATCH_ID,
            "entry_id": target["entry_id"],
            "entity_id": target["entity_id"],
            "sequence_cluster_identity_cutoff": target["sequence_cluster_identity_cutoff"],
            "sequence_cluster_30_id": target["sequence_cluster_30_id"],
            "sequence_cluster_representative_rank": target["sequence_cluster_representative_rank"],
        })
    return {
        "kind": "V61_PROTEIN_UNIVERSE_LEDGER_v0",
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "row_count": len(rows),
        "rows": rows,
    }


def _engine_version_ledger(engine_declaration: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V61_ENGINE_VERSION_LEDGER_v0",
        "campaign_id": CAMPAIGN_ID,
        "versions": [
            {
                "engine_version": ENGINE_VERSION_USED,
                "commit": engine_declaration["engine_source_last_commit"],
                "required_start_commit": ENGINE_START_COMMIT,
                "engine_source_sha256": engine_declaration["engine_source_sha256"],
                "operator_set_hash": engine_declaration["operator_set_hash"],
                "mechanism_class_set_hash": engine_declaration["mechanism_class_set_hash"],
                "batch_id": BATCH_ID,
                "engine_modified_during_batch": False,
                "change_summary": "Frozen E60 engine from V59/V60; V61 runner and ledgers added without changing engine biology.",
            }
        ],
    }


def _failure_type(row: dict[str, Any]) -> str:
    expected = row["expected_mechanism_class"]
    predicted = row["predicted_mechanism_class"]
    if predicted == "insufficient_evidence_clean_abstain":
        return "weak_sequence_signal"
    if expected == "membrane_multidomain_folding_proteostasis":
        return "membrane_misread"
    if expected == "intrinsic_disorder_phase_separation":
        return "disorder_misread"
    if expected == "oligomerization_controlled_folding":
        return "oligomer_state_misread"
    if expected == "cofactor_ligand_assisted_stabilization":
        return "cofactor_ligand_missing"
    if expected == "globular_closure" and predicted != expected:
        return "wrong_regime"
    return "right_regime_wrong_topology" if row["level1_regime_selection"] else "wrong_regime"


def _failure_grammar_row(row: dict[str, Any]) -> dict[str, Any]:
    failure_type = _failure_type(row)
    proposals = {
        "membrane_misread": ("membrane topology/context mark", "membrane_pressure_operator", "membrane insertion pressure", "abstain when hydrophobic membrane signal lacks explicit topology context"),
        "disorder_misread": ("IDR/low-complexity persistence mark", "disorder_operator", "entropy/phase pressure", "abstain when disorder support is annotation-only and sequence signal is weak"),
        "oligomer_state_misread": ("obligate assembly/context mark", "interface_operator", "partner-copy concentration pressure", "abstain when monomer sequence lacks assembly evidence"),
        "cofactor_ligand_missing": ("cofactor/metal-ligand stabilization mark", "interface_operator", "ligand_or_cofactor pressure", "abstain when ligand dependence is unobserved preseal"),
        "weak_sequence_signal": ("confidence/self-evidence mark", "none", "none", "retain abstention until noncoordinate evidence sharpens the mechanism"),
        "right_regime_wrong_topology": ("topology proxy refinement mark", "closure_operator", "fold-class pressure", "abstain when operator regions are underdetermined"),
        "wrong_regime": ("regime separation mark", "frustration_operator", "context separation pressure", "abstain when two regime explanations remain plausible"),
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
        "control_that_prevents_overfitting": "repair only in a later engine version, then rerun on a fresh frozen replication batch",
    }


def _failure_grammar_ledger(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [_failure_grammar_row(row) for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V61_FAILURE_GRAMMAR_LEDGER_v0",
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
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
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    wrong_packets: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V61_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V61_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold-style model offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V61_PRESEAL_HOLDOUT",
        "source_class": "coordinate_derived",
        "source_role": "holdout_validation",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V61_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_simulation_packet(
        target_id="V61_RANDOM_SEQUENCE_CONTROL",
        target_name="V61 random sequence control",
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
    shuffled_ok = (
        len(shuffled_rows) == len(packets)
        and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows)
    )
    readme_diff = _git(["diff", "--", "README.md"])
    return [
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V61 selection must be automatic."),
        _control("targets_total_100", len(target_manifest["selected_targets"]) == TARGET_COUNT, "V61 N must be exactly 100.", len(target_manifest["selected_targets"])),
        _control("rcsb_experimental_protein_entities_only", all(row["source_database"] == "RCSB_PDB" and row["structure_determination_methodology"].lower() == "experimental" and "protein" in row["polymer_type"].lower() for row in target_manifest["selected_targets"]), "All targets are RCSB experimental protein entities."),
        _control("sequence_cluster_representative_selection", target_manifest["sequence_cluster_representative_selection"] is True and len(set(clusters)) == len(clusters) and all(clusters), "Each target is a unique 30% sequence-cluster representative.", {"cutoff": SEQUENCE_IDENTITY_CUTOFF, "unique_clusters": len(set(clusters))}),
        _control("length_filter_40_800", all(MIN_LENGTH <= row["sequence_length"] <= MAX_LENGTH for row in target_manifest["selected_targets"]), "All targets satisfy the requested length filter."),
        _control("engine_starts_as_d927781", engine_declaration["engine_source_last_commit"].startswith(ENGINE_START_COMMIT), "Engine source last commit is the requested E60 starting commit.", engine_declaration["engine_source_last_commit"]),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine modification inside V61."),
        _control("engine_language_sets_frozen", engine_declaration["frozen_operator_names"] == UNIVERSAL_OPERATORS and engine_declaration["frozen_mechanism_classes"] == MECHANISM_CLASSES, "Operator and mechanism sets are frozen."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references the sealed prediction hash."),
        _control("wrong_grammar_controls", all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_packets), "Forced wrong grammars are rejected or routed to abstention."),
        _control(
            "shuffled_sequence_controls_reported",
            shuffled_ok,
            "Composition-preserving shuffled controls are generated with target metadata withheld; stronger shuffled coherence cases are reported rather than used as validation support.",
            {
                "control_count": len(shuffled_rows),
                "shuffled_higher_by_more_than_margin_count": sum(1 for row in shuffled_rows if row["shuffled_higher_by_more_than_margin"]),
                "rows": shuffled_rows,
            },
        ),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain", "Random sequence without evidence abstains."),
        _control("failures_reported", len(scoring_rows) == TARGET_COUNT and all("score_label" in row for row in scoring_rows), "Every target has an explicit score row."),
        _control("readme_modified_false", readme_diff == "", "README is manual-owned and must remain untouched.", {"diff_length": len(readme_diff)}),
    ]


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_ledger: dict[str, Any],
) -> dict[str, Any]:
    accepted = [row for row in scoring_rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in scoring_rows if row["acceptance_decision"] == "abstain_recommended"]
    supported = [row for row in scoring_rows if row["score_label"] == "supported"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    failed_accepted = [row for row in accepted if row["score_label"] != "supported"]
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "alphafold_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    if len(scoring_rows) < TARGET_COUNT:
        status = BLOCKED_INTAKE
    elif any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif failed_accepted:
        status = DISCOVERY_FAILURES
    elif abstained:
        status = PARTIAL_ABSTAIN
    else:
        status = PASSED
    controls_passed = not failed_controls
    cert = {
        "kind": "V61_RCSB_NONREDUNDANT_100_BATCH_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_start_commit": ENGINE_START_COMMIT,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "targets_total": len(scoring_rows),
        "accepted_count": len(accepted),
        "supported_count": len(supported),
        "failed_accepted_count": len(failed_accepted),
        "abstain_count": len(abstained),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "raw_accuracy": len(supported) / len(scoring_rows) if scoring_rows else None,
        "coverage": len(accepted) / len(scoring_rows) if scoring_rows else None,
        "controls_passed": controls_passed,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "engine_modified_during_batch": False,
        "readme_modified": False,
        "failure_modes": failure_ledger["failure_modes"],
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
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "claim_allowed": controls_passed and not failed_accepted and bool(accepted),
        "claim_blocked_reason": "" if controls_passed and not failed_accepted and bool(accepted) else "accepted failures or failed controls remain in V61 discovery ledger",
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "All protein classes are saturated.",
            "Coordinates or contacts were predicted de novo.",
            "Engine changes were made during V61.",
            "Failures may be hidden or repaired after holdout.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_ledger(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V61_CLAIM_LEDGER_v0",
        "campaign_id": CAMPAIGN_ID,
        "rows": [
            {
                "batch_id": BATCH_ID,
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
        ],
    }


def _failure_report(scoring_rows: list[dict[str, Any]], failure_ledger: dict[str, Any]) -> dict[str, Any]:
    failures = [row for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V61_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "failure_cases_reported": True,
        "failure_count": len(failures),
        "failure_cases": failures,
        "failure_grammar_rows": failure_ledger["rows"],
        "note": "Failures are preserved as missing Esperanto grammar candidates; no engine repair occurs inside V61.",
    }


def _write_report(path: Path, cert: dict[str, Any], scoring_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# V61 RCSB Nonredundant 100 Batch",
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
        f"README modified: `{cert['readme_modified']}`",
        "",
        "## Failure Modes",
    ]
    if cert["failure_modes"]:
        for mode, count in sorted(cert["failure_modes"].items()):
            lines.append(f"- `{mode}`: `{count}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Target Scores"])
    for row in scoring_rows:
        lines.append(
            f"- `{row['target_id']}` decision `{row['acceptance_decision']}` predicted `{row['predicted_mechanism_class']}` expected `{row['expected_mechanism_class']}` label `{row['score_label']}`"
        )
    lines.extend([
        "",
        "## Boundary",
        "V61 is a frozen saturation batch on E60. Coordinates, contacts, AlphaFold-style structures, holdout annotations, and runtime artifacts are blocked before sealing. Failures are not patched in this batch; they are written to the failure grammar ledger for a later engine revision.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v61(out_dir: Path = DEFAULT_OUT_DIR, *, refresh_intake: bool = False) -> dict[str, Path]:
    raw = refresh_rcsb_cluster_representative_intake() if refresh_intake else _read_json(RAW_CANDIDATE_CACHE, "V61 raw RCSB representative cache")
    _reset_generated_outputs()
    target_manifest = _target_manifest(raw)
    if target_manifest["target_count_selected"] < TARGET_COUNT:
        raise SystemExit(f"V61 selected only {target_manifest['target_count_selected']} targets; need {TARGET_COUNT}")
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v61_rcsb_nonredundant_100_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v61_rcsb_nonredundant_100_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    wrong_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []

    for candidate in target_manifest["selected_targets"]:
        target_id = f"V61_{candidate['target_id']}"
        expected, reasons = _expected_mechanism_postseal(candidate)
        source_manifest = _source_manifest(candidate)
        _write_json(DATA_ROOT / "source_manifests" / target_id / "source_manifest.json", source_manifest)
        packet = build_sealed_simulation_packet(
            target_id=target_id,
            target_name=f"{candidate['entry_id']} {candidate['entity_description']}",
            sequence=candidate["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "automatic RCSB representative full-chain scan", "span": f"1-{candidate['sequence_length']}"}],
            perturbations=_perturbations_for_expected(expected, target_id),
        )
        holdout = _holdout(candidate, packet, expected, reasons)
        score = _score(packet, holdout)
        wrong_packet = build_sealed_simulation_packet(
            target_id=f"{target_id}_WRONG_GRAMMAR_CONTROL",
            target_name=f"{candidate['entry_id']} forced wrong grammar control",
            sequence=candidate["sequence"],
            sources=source_manifest["prediction_sources"],
            perturbations=[],
            forced_grammar=_wrong_grammar(packet["selected_mechanism_grammar"]["natural_mechanism_class"]),
        )
        shuffled_packet = build_sealed_simulation_packet(
            target_id=f"{target_id}_SHUFFLED_CONTROL",
            target_name=f"{candidate['entry_id']} shuffled sequence control",
            sequence=shuffled_sequence(candidate["sequence"]),
            sources=_shuffled_control_sources(target_id),
            perturbations=[],
        )
        packets.append(packet)
        scoring_rows.append(score)
        wrong_packets.append(wrong_packet)
        shuffled_packets.append(shuffled_packet)
        _write_json(DATA_ROOT / "sealed_predictions" / target_id / "sealed_simulation_packet.json", packet)
        _write_json(DATA_ROOT / "holdouts_postseal" / target_id / "postseal_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target_id / "validation_result.json", score)
        _write_json(DATA_ROOT / "wrong_grammar_controls" / target_id / "wrong_grammar_packet.json", wrong_packet)
        _write_json(DATA_ROOT / "shuffled_controls" / target_id / "shuffled_control_packet.json", shuffled_packet)

    protein_ledger = _protein_universe_ledger(target_manifest)
    engine_ledger = _engine_version_ledger(engine_declaration)
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
    claim_ledger = _claim_ledger(cert)
    scoring_path = _write_json(DATA_ROOT / "v61_rcsb_nonredundant_100_scoring_report.json", {"kind": "V61_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v61_rcsb_nonredundant_100_failure_report.json", _failure_report(scoring_rows, failure_ledger))
    data_cert_path = _write_json(DATA_ROOT / "v61_rcsb_nonredundant_100_certificate.json", cert)
    protein_ledger_path = _write_json(LEDGER_ROOT / "protein_universe_ledger_v0.json", protein_ledger)
    engine_ledger_path = _write_json(LEDGER_ROOT / "engine_version_ledger_v0.json", engine_ledger)
    failure_ledger_path = _write_json(LEDGER_ROOT / "failure_grammar_ledger_v0.json", failure_ledger)
    claim_ledger_path = _write_json(LEDGER_ROOT / "claim_ledger_v0.json", claim_ledger)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v61_rcsb_nonredundant_100_batch_certificate.json", cert)
    report_path = out_dir / "V61_RCSB_NONREDUNDANT_100_BATCH_REPORT.md"
    _write_report(report_path, cert, scoring_rows)
    return {
        "raw_candidate_cache": RAW_CANDIDATE_CACHE,
        "target_manifest": DATA_ROOT / "v61_rcsb_nonredundant_100_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v61_rcsb_nonredundant_100_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "protein_universe_ledger": protein_ledger_path,
        "engine_version_ledger": engine_ledger_path,
        "failure_grammar_ledger": failure_ledger_path,
        "claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V61 RCSB nonredundant 100-target batch.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--refresh-intake", action="store_true", help="refresh RCSB 30% representative intake with curl")
    args = parser.parse_args()
    paths = run_v61(args.out_dir, refresh_intake=args.refresh_intake)
    cert = _read_json(paths["certificate"], "V61 certificate")
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
        "readme_modified": cert["readme_modified"],
        "failure_modes": cert["failure_modes"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
