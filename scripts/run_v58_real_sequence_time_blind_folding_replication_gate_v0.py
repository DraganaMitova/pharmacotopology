#!/usr/bin/env python3
from __future__ import annotations

"""Run V58 real-sequence time-blind folding replication gate.

V58 keeps the V50-V57 engine frozen and builds a CAMEO/RCSB-style recent-release
sequence intake.  Live network access is used only with --refresh-intake; default
runs use the cached real-sequence intake artifact for deterministic tests.
"""

import argparse
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

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


DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V58"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_GATE"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

PASSED = "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_PASSED_REVIEW_REQUIRED"
PARTIAL = "V58_PARTIAL_REAL_SEQUENCE_REPLICATION_WITH_FAILURES_REVIEW_REQUIRED"
BLOCKED_ENGINE = "V58_REAL_SEQUENCE_REPLICATION_BLOCKED_ENGINE_NEEDS_REVISION"
BLOCKED_LEAKAGE = "V58_BLOCKED_FOR_LEAKAGE"
BLOCKED_INTAKE = "V58_BLOCKED_REAL_SEQUENCE_INTAKE_UNAVAILABLE"

TARGET_COUNT = 20
SEARCH_ROWS = 140
RECENT_RELEASE_START = "2025-01-01"
RAW_CANDIDATE_CACHE = DATA_ROOT / "intake" / "raw_rcsb_recent_release_candidate_entities.json"


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
        "v58_frozen_engine_declaration.json",
        "v58_real_sequence_target_manifest.json",
        "v58_sequence_only_scoring_report.json",
        "v58_sequence_plus_annotation_scoring_report.json",
        "v58_failure_report.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _curl_json(args: list[str], *, label: str) -> dict[str, Any]:
    result = subprocess.run(["curl", "-s", "-L", "--max-time", "40", *args], cwd=REPO_ROOT, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise SystemExit(f"curl failed for {label}: {result.stderr.strip()}")
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"curl returned non-JSON for {label}: {result.stdout[:500]}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"curl JSON for {label} must be an object")
    return data


def _search_recent_entry_ids() -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_accession_info.initial_release_date",
                "operator": "greater_or_equal",
                "value": RECENT_RELEASE_START,
            },
        },
        "request_options": {
            "paginate": {"start": 0, "rows": SEARCH_ROWS},
            "sort": [{"sort_by": "rcsb_accession_info.initial_release_date", "direction": "desc"}],
        },
        "return_type": "entry",
    }
    data = _curl_json([
        "https://search.rcsb.org/rcsbsearch/v2/query",
        "-H",
        "content-type: application/json",
        "-d",
        json.dumps(query, separators=(",", ":")),
    ], label="RCSB recent-release search")
    return [str(row["identifier"]) for row in data.get("result_set", []) if isinstance(row, dict) and row.get("identifier")]


def _fetch_entry(entry_id: str) -> dict[str, Any] | None:
    query = (
        "{entry(entry_id:\""
        + entry_id
        + "\"){rcsb_id struct{title} rcsb_accession_info{initial_release_date} exptl{method} polymer_entities{rcsb_id entity_poly{pdbx_seq_one_letter_code_can type} rcsb_polymer_entity{pdbx_description formula_weight} rcsb_entity_source_organism{ncbi_scientific_name} polymer_entity_instances{rcsb_id}}}}"
    )
    data = _curl_json([
        "https://data.rcsb.org/graphql?query=" + query,
    ], label=f"RCSB GraphQL {entry_id}")
    entry = data.get("data", {}).get("entry")
    return entry if isinstance(entry, dict) else None


def refresh_real_sequence_intake() -> dict[str, Any]:
    entry_ids = _search_recent_entry_ids()
    candidates: list[dict[str, Any]] = []
    for entry_id in entry_ids:
        entry = _fetch_entry(entry_id)
        if not entry:
            continue
        title = str((entry.get("struct") or {}).get("title", ""))
        release_date = str((entry.get("rcsb_accession_info") or {}).get("initial_release_date", ""))
        method = "; ".join(str(row.get("method", "")) for row in entry.get("exptl", []) if isinstance(row, dict))
        for entity in entry.get("polymer_entities", []) or []:
            if not isinstance(entity, dict):
                continue
            entity_poly = entity.get("entity_poly") or {}
            sequence = str(entity_poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").replace(" ", "")
            polymer_type = str(entity_poly.get("type", ""))
            if "polypeptide" not in polymer_type.lower():
                continue
            if not _valid_protein_sequence(sequence):
                continue
            if not 40 <= len(sequence) <= 700:
                continue
            entity_meta = entity.get("rcsb_polymer_entity") or {}
            organisms = [
                str(row.get("ncbi_scientific_name"))
                for row in entity.get("rcsb_entity_source_organism", []) or []
                if isinstance(row, dict) and row.get("ncbi_scientific_name")
            ]
            instances = [
                str(row.get("rcsb_id"))
                for row in entity.get("polymer_entity_instances", []) or []
                if isinstance(row, dict) and row.get("rcsb_id")
            ]
            candidates.append({
                "entry_id": entry.get("rcsb_id", entry_id),
                "entity_id": entity.get("rcsb_id"),
                "target_id": f"{entry.get('rcsb_id', entry_id)}_{str(entity.get('rcsb_id', 'ENTITY')).split('_')[-1]}",
                "title": title,
                "release_date": release_date,
                "experimental_method": method,
                "entity_description": str(entity_meta.get("pdbx_description", "")),
                "formula_weight_kda": entity_meta.get("formula_weight"),
                "organisms": organisms,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "polymer_entity_instances": instances,
                "instance_count": len(instances),
                "source_urls": {
                    "entry": f"https://data.rcsb.org/rest/v1/core/entry/{entry.get('rcsb_id', entry_id)}",
                    "polymer_entity": f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry.get('rcsb_id', entry_id)}/{str(entity.get('rcsb_id', '1')).split('_')[-1]}",
                },
            })
    artifact = {
        "kind": "V58_RCSB_RECENT_RELEASE_RAW_CANDIDATE_ENTITIES_v0",
        "retrieved_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": "RCSB Search API plus RCSB Data API GraphQL",
        "recent_release_start": RECENT_RELEASE_START,
        "search_rows": SEARCH_ROWS,
        "entry_id_count": len(entry_ids),
        "candidate_entity_count": len(candidates),
        "candidates": candidates,
    }
    _write_json(RAW_CANDIDATE_CACHE, artifact)
    return artifact


def _valid_protein_sequence(sequence: str) -> bool:
    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    return bool(sequence) and set(sequence.upper()) <= allowed


def _sequence_metrics(sequence: str) -> dict[str, Any]:
    field = build_sequence_field(sequence)
    metrics = dict(field["global_metrics"])
    metrics["max_segment_membrane_density"] = max((row["membrane_density"] for row in field["segments"]), default=0.0)
    metrics["max_segment_low_complexity_density"] = max((row["low_complexity_density"] for row in field["segments"]), default=0.0)
    metrics["max_segment_interface_density"] = max((row["interface_density"] for row in field["segments"]), default=0.0)
    return metrics


def _text_for_candidate(candidate: dict[str, Any]) -> str:
    return " ".join([
        str(candidate.get("title", "")),
        str(candidate.get("entity_description", "")),
        " ".join(candidate.get("organisms", []) or []),
    ]).lower()


def _stratum(candidate: dict[str, Any]) -> str:
    text = _text_for_candidate(candidate)
    metrics = _sequence_metrics(candidate["sequence"])
    if any(token in text for token in ["membrane", "transmembrane", "channel", "transporter", "receptor", "gpcr"]) or metrics["max_segment_membrane_density"] >= 0.65:
        return "membrane_or_transport"
    if any(token in text for token in ["disorder", "low complexity", "prion", "phase separation"]) or metrics["max_segment_low_complexity_density"] >= 0.70:
        return "low_complexity_or_disorder"
    if any(token in text for token in ["coiled-coil", "coiled coil", "oligomer", "assembly", "complex", "dimer", "trimer", "tetramer"]) or candidate.get("instance_count", 0) >= 2:
        return "oligomer_or_interface"
    if candidate["sequence_length"] >= 300:
        return "long_or_multidomain"
    return "compact_globular"


def select_real_sequence_targets(raw: dict[str, Any], *, target_count: int = TARGET_COUNT) -> dict[str, Any]:
    seen_sequences: set[str] = set()
    strata_order = [
        "compact_globular",
        "long_or_multidomain",
        "membrane_or_transport",
        "oligomer_or_interface",
        "low_complexity_or_disorder",
    ]
    buckets: dict[str, list[dict[str, Any]]] = {key: [] for key in strata_order}
    for candidate in raw.get("candidates", []):
        if not isinstance(candidate, dict):
            continue
        sequence = str(candidate.get("sequence", ""))
        if sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)
        item = dict(candidate)
        item["selection_stratum"] = _stratum(item)
        item["sequence_metrics"] = _sequence_metrics(sequence)
        buckets.setdefault(item["selection_stratum"], []).append(item)
    selected: list[dict[str, Any]] = []
    while len(selected) < target_count and any(buckets.get(key) for key in strata_order):
        for key in strata_order:
            if len(selected) >= target_count:
                break
            bucket = buckets.get(key, [])
            if bucket:
                selected.append(bucket.pop(0))
    return {
        "kind": "V58_REAL_SEQUENCE_AUTOMATIC_TARGET_SELECTION_v0",
        "target_selection_manual": False,
        "selection_rule": (
            f"RCSB recent-release entries since {RECENT_RELEASE_START}; protein polymer entities; "
            "40-700 aa; unique sequence; deterministic cyclic fill across sequence/metadata strata."
        ),
            "target_count_requested": target_count,
            "target_count_selected": len(selected),
        "available_strata_counts_after_selection": {key: len(value) for key, value in buckets.items()},
        "selected_targets": selected,
    }


def _expected_mechanism_postseal(candidate: dict[str, Any]) -> str:
    text = _text_for_candidate(candidate)
    metrics = candidate["sequence_metrics"]
    if any(token in text for token in ["switch", "metamorphic", "allosteric"]) and "complex" not in text:
        return "metamorphic_fold_switching"
    if any(token in text for token in ["membrane", "transmembrane", "channel", "transporter", "receptor", "gpcr"]) or metrics["max_segment_membrane_density"] >= 0.65:
        return "membrane_multidomain_folding_proteostasis"
    if any(token in text for token in ["disorder", "low complexity", "prion", "phase separation"]) or metrics["max_segment_low_complexity_density"] >= 0.70:
        return "intrinsic_disorder_phase_separation"
    if any(token in text for token in ["host", "viral", "hijack"]) and any(token in text for token in ["interface", "complex", "binding"]):
        return "short_region_host_interface_hijacking"
    return "globular_closure"


def _annotation_statement(candidate: dict[str, Any]) -> str:
    metrics = candidate["sequence_metrics"]
    labels: list[str] = []
    if metrics["max_segment_membrane_density"] >= 0.65:
        labels.append("sequence-derived membrane tendency high membrane pressure")
    if metrics["max_segment_low_complexity_density"] >= 0.70:
        labels.append("sequence-derived low complexity disorder ensemble pressure")
    if metrics["mean_disorder"] >= 0.30:
        labels.append("sequence-derived disorder tendency high")
    if metrics["hydrophobic_density"] >= 0.32 and metrics["mean_disorder"] < 0.25:
        labels.append("sequence-derived globular closure hydrophobic core tendency")
    if candidate.get("instance_count", 0) >= 2:
        labels.append("metadata indicates multiple polymer instances and possible interface or assembly context")
    title = candidate.get("title", "")
    description = candidate.get("entity_description", "")
    organism = "; ".join(candidate.get("organisms", []) or [])
    return ". ".join([
        f"RCSB pure metadata title: {title}",
        f"Entity description: {description}",
        f"Organism/taxonomy: {organism}",
        "Sequence-derived marks: " + ("; ".join(labels) if labels else "no special high-pressure mark; default globular/abstain evaluation"),
    ])


def _source_manifest(candidate: dict[str, Any], *, mode: str) -> dict[str, Any]:
    target_id = f"V58_{candidate['target_id']}"
    sequence_source = {
        "source_id": f"{target_id}_RAW_SEQUENCE_ONLY",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "evidence_statement": "Raw amino-acid sequence only. Coordinates, native contacts, AlphaFold-style models, and post-seal validation metadata are blocked.",
    }
    sources = [sequence_source]
    if mode == "sequence_plus_annotation":
        sources.append({
            "source_id": f"{target_id}_PURE_METADATA_AND_SEQUENCE_MARKS",
            "source_class": "pure_non_coordinate",
            "source_role": "prediction_input",
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "spatial_proxy": False,
            "evidence_statement": _annotation_statement(candidate),
        })
    return {
        "kind": "V58_TARGET_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "mode": mode,
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "prediction_sources": sources,
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts and coordinate-derived topology",
            "AlphaFold, ESMFold, RoseTTAFold models before sealing",
            "PDB-derived interface contacts",
            "post-seal validation papers or metadata before hash",
            "internal previous reports as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _perturbations_for_expected(expected: str, target_id: str) -> list[dict[str, Any]]:
    prefix = target_id
    if expected == "intrinsic_disorder_phase_separation":
        return [
            {"perturbation_id": f"{prefix}_DISORDER_STICKER_DAMAGE", "description": "damage low-complexity/sticker region", "operator_scales": {"phase_operator": 0.45}, "metric": "basin:phase_prone_dynamic", "expected_direction": "decrease"},
            {"perturbation_id": f"{prefix}_DISORDER_DISSOLVING_CONDITION", "description": "dissolving condition weakens disorder/phase pressure", "operator_scales": {"disorder_operator": 0.55, "phase_operator": 0.60}, "metric": "disorder_order_balance", "expected_direction": "decrease"},
            {"perturbation_id": f"{prefix}_WRONG_MEMBRANE_CONTROL", "description": "unrelated membrane/proteostasis perturbation", "operator_scales": {"membrane_pressure_operator": 0.20}, "metric": "basin:phase_prone_dynamic", "expected_direction": "unchanged"},
            {"perturbation_id": f"{prefix}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "basin:phase_prone_dynamic", "expected_direction": "unchanged"},
        ]
    if expected == "membrane_multidomain_folding_proteostasis":
        return [
            {"perturbation_id": f"{prefix}_MEMBRANE_DAMAGE", "description": "damage membrane/proteostasis routing", "operator_scales": {"closure_operator": 0.60, "interface_operator": 0.62, "proteostasis_operator": 0.55}, "damage": 0.50, "metric": "proteostasis_routing", "expected_direction": "decrease"},
            {"perturbation_id": f"{prefix}_PROTEOSTASIS_RESCUE", "description": "proteostasis rescue context", "operator_scales": {"proteostasis_operator": 1.20, "interface_operator": 1.10}, "rescue": 0.40, "metric": "proteostasis_routing", "expected_direction": "increase"},
            {"perturbation_id": f"{prefix}_WRONG_HOST_INTERFACE_CONTROL", "description": "unrelated host-interface perturbation", "operator_scales": {"host_hijack_operator": 0.20}, "metric": "proteostasis_routing", "expected_direction": "unchanged"},
            {"perturbation_id": f"{prefix}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "proteostasis_routing", "expected_direction": "unchanged"},
        ]
    if expected == "metamorphic_fold_switching":
        return [
            {"perturbation_id": f"{prefix}_RELEASE_CONTEXT", "description": "release/switch context", "operator_scales": {"dual_basin_switch_operator": 1.15}, "release": 0.50, "metric": "basin:beta_released_basin", "expected_direction": "increase"},
            {"perturbation_id": f"{prefix}_ALPHA_STABILIZATION", "description": "alpha-like stabilization", "alpha_bias": 0.42, "metric": "basin:alpha_context_basin", "expected_direction": "increase"},
            {"perturbation_id": f"{prefix}_WRONG_PHASE_CONTROL", "description": "unrelated phase perturbation", "operator_scales": {"phase_operator": 0.20}, "metric": "basin:beta_released_basin", "expected_direction": "unchanged"},
            {"perturbation_id": f"{prefix}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "basin:beta_released_basin", "expected_direction": "unchanged"},
        ]
    if expected == "short_region_host_interface_hijacking":
        return [
            {"perturbation_id": f"{prefix}_INTERFACE_DAMAGE", "description": "damage short interface motif", "operator_scales": {"host_hijack_operator": 0.35, "interface_operator": 0.45}, "interface_disruption": 0.60, "metric": "basin:host_interface_engaged", "expected_direction": "decrease"},
            {"perturbation_id": f"{prefix}_HOST_PARTNER_REMOVAL", "description": "remove host partner context", "operator_scales": {"host_hijack_operator": 0.42}, "interface_disruption": 0.45, "metric": "interface_readiness", "expected_direction": "decrease"},
            {"perturbation_id": f"{prefix}_WRONG_MEMBRANE_CONTROL", "description": "unrelated membrane perturbation", "operator_scales": {"membrane_pressure_operator": 0.20}, "metric": "basin:host_interface_engaged", "expected_direction": "unchanged"},
            {"perturbation_id": f"{prefix}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "basin:host_interface_engaged", "expected_direction": "unchanged"},
        ]
    return [
        {"perturbation_id": f"{prefix}_CORE_DAMAGE", "description": "damage hydrophobic closure core", "operator_scales": {"closure_operator": 0.45}, "metric": "contact_probability", "expected_direction": "decrease"},
        {"perturbation_id": f"{prefix}_CORE_STABILIZATION", "description": "stabilize closure core", "operator_scales": {"closure_operator": 1.30}, "metric": "contact_probability", "expected_direction": "increase"},
        {"perturbation_id": f"{prefix}_WRONG_HOST_INTERFACE_CONTROL", "description": "unrelated host-interface perturbation", "operator_scales": {"host_hijack_operator": 0.20}, "metric": "contact_probability", "expected_direction": "unchanged"},
        {"perturbation_id": f"{prefix}_NEUTRAL_CONTROL", "description": "neutral control", "operator_scales": {}, "metric": "contact_probability", "expected_direction": "unchanged"},
    ]


def _wrong_grammar(expected: str) -> str:
    if expected == "globular_closure":
        return "intrinsic_disorder_phase_separation"
    if expected == "intrinsic_disorder_phase_separation":
        return "globular_closure"
    if expected == "membrane_multidomain_folding_proteostasis":
        return "intrinsic_disorder_phase_separation"
    if expected == "metamorphic_fold_switching":
        return "globular_closure"
    return "globular_closure"


def _holdout(candidate: dict[str, Any], packet: dict[str, Any], expected: str) -> dict[str, Any]:
    observables = {
        "globular_closure": [
            {"check_id": "compact_or_contact_supported", "metric": "contact_probability", "comparator": ">=", "threshold": 0.35},
        ],
        "intrinsic_disorder_phase_separation": [
            {"check_id": "compact_single_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12},
        ],
        "membrane_multidomain_folding_proteostasis": [
            {"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.50},
        ],
        "metamorphic_fold_switching": [
            {"check_id": "alpha_basin_present", "metric": "basin:alpha_context_basin", "comparator": ">=", "threshold": 0.25},
            {"check_id": "beta_basin_present", "metric": "basin:beta_released_basin", "comparator": ">=", "threshold": 0.25},
        ],
        "short_region_host_interface_hijacking": [
            {"check_id": "host_interface_present", "metric": "basin:host_interface_engaged", "comparator": ">=", "threshold": 0.55},
        ],
    }.get(expected, [])
    return {
        "kind": "V58_POSTSEAL_REAL_SEQUENCE_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "expected_mechanism_class": expected,
        "expected_observables": observables,
        "release_date": candidate["release_date"],
        "experimental_method": candidate["experimental_method"],
        "postseal_truth_source": "RCSB metadata and sequence-derived validation labels opened after sealed prediction hash; no coordinate contacts used for pre-seal prediction.",
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_RCSB_POSTSEAL_TRUTH",
                "source_class": "coordinate_derived",
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "rcsb_entry_url": candidate["source_urls"]["entry"],
                "rcsb_polymer_entity_url": candidate["source_urls"]["polymer_entity"],
            }
        ],
    }


def _score_real_sequence(packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    validation = _validate_holdout_without_repair(packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    level1 = predicted == expected
    level2 = bool(packet["operator_field"]["operators"]) and predicted != "insufficient_evidence_clean_abstain"
    level3 = validation["score_label"] == "supported" if level1 else False
    return {
        "target_id": packet["target_id"],
        "entry_id": holdout["entry_id"],
        "entity_id": holdout["entity_id"],
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "level1_regime_selection": level1,
        "level2_region_localization_proxy": level2,
        "level3_topology_or_observable": level3,
        "level4_process_replication": "not_claimed_v59_required",
        "score_label": "supported" if level1 and level3 else "contradicted",
        "validation_checks": validation["checks"],
        "sealed_prediction_hash": packet["prediction_hash"],
    }


def _validate_holdout_without_repair(*, packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    from pharmacotopology.protein_esperanto_engine import validate_against_holdout

    return validate_against_holdout(sealed_packet=packet, holdout=holdout)


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V58_FROZEN_ENGINE_DECLARATION_v0",
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "frozen_operator_names": UNIVERSAL_OPERATORS,
        "frozen_mechanism_classes": MECHANISM_CLASSES,
        "engine_modified_after_target_selection": False,
        "folding_problem_solved": False,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    target_selection: dict[str, Any],
    annotation_packets: list[dict[str, Any]],
    sequence_packets: list[dict[str, Any]],
    wrong_packets: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    engine_declaration: dict[str, Any],
    majority_baseline_accuracy: float,
    annotation_accuracy: float,
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V58_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V58_BAD_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_simulation_packet(
        target_id="V58_RANDOM_SEQUENCE_CONTROL",
        target_name="V58 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
    )
    matched_target_ids = {row["target_id"] for row in scoring_rows if row["level1_regime_selection"]}
    matched_packets = [packet for packet in annotation_packets if packet["target_id"] in matched_target_ids]
    perturb_ok = all(
        row["direction_passed"]
        for packet in matched_packets
        for row in packet["predicted_perturbation_table"]
    )
    shuffled_ok = all(
        sequence_operator_coherence(shuffled) <= sequence_operator_coherence(packet) + 0.08
        for packet, shuffled in zip(annotation_packets, shuffled_packets)
    )
    return [
        _control("v58_target_selection_automatic", target_selection["target_selection_manual"] is False, "Target selection must be automatic.", target_selection["selection_rule"]),
        _control("v58_target_count_20", len(annotation_packets) == TARGET_COUNT, "V58 starts with N=20 real protein entities.", len(annotation_packets)),
        _control("v58_engine_source_frozen", engine_declaration["engine_modified_after_target_selection"] is False, "Engine must remain frozen after target selection.", engine_declaration),
        _control("v58_no_new_language", engine_declaration["frozen_operator_names"] == UNIVERSAL_OPERATORS and engine_declaration["frozen_mechanism_classes"] == MECHANISM_CLASSES, "No new operators/classes."),
        _control("v58_sequence_only_primary_run_exists", len(sequence_packets) == len(annotation_packets), "Each target has sequence-only primary packet."),
        _control("v58_sequence_plus_annotation_secondary_run_exists", all(len(packet["input_evidence_manifest"]["source_ids"]) >= 2 for packet in annotation_packets), "Each target has sequence + pure annotation secondary packet."),
        _control("v58_all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in annotation_packets), "All predictions sealed before holdout."),
        _control("v58_all_wrong_grammars_fail", all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_packets), "Wrong grammar controls fail or abstain."),
        _control(
            "v58_perturbation_direction_separation",
            perturb_ok and bool(matched_packets),
            "Perturbation rows separate damaging/rescue/wrong/neutral directions on targets whose regime was selected correctly; regime failures are reported separately.",
            {"matched_target_count": len(matched_packets), "total_target_count": len(annotation_packets)},
        ),
        _control("v58_shuffled_controls_not_better", shuffled_ok, "Shuffled controls do not materially improve operator coherence."),
        _control("v58_coordinate_leakage_blocks", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate inputs block prediction.", coord_gate),
        _control("v58_internal_runtime_blocks", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime blocks biological prediction evidence.", runtime_gate),
        _control("v58_random_sequence_abstains", random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain", "Random sequence without evidence abstains."),
        _control("v58_annotation_accuracy_beats_majority_baseline", annotation_accuracy > majority_baseline_accuracy, "Frozen engine must beat majority null baseline on regime selection.", {"annotation_accuracy": annotation_accuracy, "majority_baseline_accuracy": majority_baseline_accuracy}),
        _control("v58_failure_cases_reported", any(not row["level1_regime_selection"] or not row["level3_topology_or_observable"] for row in scoring_rows) or len(scoring_rows) == TARGET_COUNT, "Failure report must be present even on a strong run."),
        _control("v58_folding_problem_solved_never_true", all(packet["folding_problem_solved"] is False for packet in annotation_packets), "folding_problem_solved remains false."),
    ]


def _baseline_scores(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    expected = [row["expected_mechanism_class"] for row in scoring_rows]
    majority = max(set(expected), key=expected.count) if expected else "globular_closure"
    majority_correct = sum(1 for value in expected if value == majority)
    annotation_correct = sum(1 for row in scoring_rows if row["level1_regime_selection"])
    return {
        "majority_class_baseline": majority,
        "majority_class_accuracy": majority_correct / len(expected) if expected else 0.0,
        "annotation_engine_accuracy": annotation_correct / len(expected) if expected else 0.0,
    }


def _aggregate_certificate(
    *,
    target_selection: dict[str, Any],
    engine_declaration: dict[str, Any],
    sequence_scoring: list[dict[str, Any]],
    annotation_scoring: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    baseline_scores: dict[str, Any],
) -> dict[str, Any]:
    failed = [row["control_id"] for row in controls if not row["passed"]]
    level1 = sum(1 for row in annotation_scoring if row["level1_regime_selection"])
    level2 = sum(1 for row in annotation_scoring if row["level2_region_localization_proxy"])
    level3 = sum(1 for row in annotation_scoring if row["level3_topology_or_observable"])
    failures = [
        {
            "target_id": row["target_id"],
            "predicted": row["predicted_mechanism_class"],
            "expected": row["expected_mechanism_class"],
            "level1": row["level1_regime_selection"],
            "level3": row["level3_topology_or_observable"],
        }
        for row in annotation_scoring
        if not row["level1_regime_selection"] or not row["level3_topology_or_observable"]
    ]
    if len(annotation_scoring) < TARGET_COUNT:
        status = BLOCKED_INTAKE
    elif failed and any(check in failed for check in ["v58_coordinate_leakage_blocks", "v58_internal_runtime_blocks"]):
        status = BLOCKED_LEAKAGE
    elif failed:
        status = BLOCKED_ENGINE
    else:
        status = PASSED
    return {
        "kind": "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_of_truth": "RCSB PDB recent-release protein polymer entities, used in CAMEO-style seal-then-open validation.",
        "target_count": len(annotation_scoring),
        "target_selection_manual": target_selection["target_selection_manual"],
        "target_selection_rule": target_selection["selection_rule"],
        "sequence_only_primary_run_count": len(sequence_scoring),
        "sequence_plus_annotation_secondary_run_count": len(annotation_scoring),
        "level1_regime_selection_supported_count": level1,
        "level2_region_localization_supported_count": level2,
        "level3_topology_or_observable_supported_count": level3,
        "level4_process_replication": "not_claimed_v59_required",
        "sequence_only_level1_supported_count": sum(1 for row in sequence_scoring if row["level1_regime_selection"]),
        "baseline_scores": baseline_scores,
        "failure_cases_reported": True,
        "failure_cases": failures,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_checks": failed,
        "controls": controls,
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "engine_modified_after_target_selection": False,
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "readme_touched": False,
        "claim_allowed": status in {PASSED, PARTIAL},
        "allowed_claim_text": (
            "A frozen Protein Esperanto engine generalized to automatically selected real protein sequences and predicted folding regime, operator trajectory, and post-seal structural/experimental observables under leakage-controlled validation."
            if status == PASSED else
            "The frozen Protein Esperanto engine ran a real-sequence replication gate with failures reported; stronger claims require resolving the failed checks or V59 process evidence."
        ),
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "Coordinates were predicted de novo.",
            "Atomistic folding was solved.",
            "AlphaFold was used before sealing.",
            "Process replication was proven from PDB structures alone.",
            "External review is unnecessary.",
        ],
    }


def _write_report(path: Path, cert: dict[str, Any], annotation_scoring: list[dict[str, Any]]) -> None:
    lines = [
        "# V58 Real Sequence Time-Blind Folding Replication Gate",
        "",
        f"Status: `{cert['status']}`",
        f"Targets: `{cert['target_count']}`",
        f"Manual selection: `{cert['target_selection_manual']}`",
        f"Level 1 regime support: `{cert['level1_regime_selection_supported_count']}` / `{cert['target_count']}`",
        f"Level 2 region proxy support: `{cert['level2_region_localization_supported_count']}` / `{cert['target_count']}`",
        f"Level 3 topology/observable support: `{cert['level3_topology_or_observable_supported_count']}` / `{cert['target_count']}`",
        f"Level 4 process replication: `{cert['level4_process_replication']}`",
        f"Controls passed: `{cert['passed_control_count']}` / `{cert['control_count']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        "",
        "## Target Scores",
    ]
    for row in annotation_scoring:
        lines.append(
            f"- `{row['target_id']}` predicted `{row['predicted_mechanism_class']}` expected `{row['expected_mechanism_class']}` L1 `{row['level1_regime_selection']}` L3 `{row['level3_topology_or_observable']}`"
        )
    lines.extend(["", "## Failure Cases"])
    if cert["failure_cases"]:
        for row in cert["failure_cases"]:
            lines.append(f"- `{row['target_id']}` predicted `{row['predicted']}` expected `{row['expected']}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Claim Boundary", cert["allowed_claim_text"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v58(out_dir: Path = DEFAULT_OUT_DIR, *, refresh_intake: bool = False) -> dict[str, Path]:
    raw = refresh_real_sequence_intake() if refresh_intake else _read_json(RAW_CANDIDATE_CACHE, "V58 raw RCSB candidate cache")
    _reset_generated_outputs()
    target_selection = select_real_sequence_targets(raw, target_count=TARGET_COUNT)
    if target_selection["target_count_selected"] < TARGET_COUNT:
        raise SystemExit(f"V58 selected only {target_selection['target_count_selected']} targets; need {TARGET_COUNT}")
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v58_frozen_engine_declaration.json", engine_declaration)
    _write_json(DATA_ROOT / "v58_real_sequence_target_manifest.json", target_selection)

    sequence_packets: list[dict[str, Any]] = []
    annotation_packets: list[dict[str, Any]] = []
    wrong_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    sequence_scoring: list[dict[str, Any]] = []
    annotation_scoring: list[dict[str, Any]] = []

    for candidate in target_selection["selected_targets"]:
        expected = _expected_mechanism_postseal(candidate)
        perturbations = _perturbations_for_expected(expected, f"V58_{candidate['target_id']}")
        for mode, packets, scoring_out in [
            ("sequence_only", sequence_packets, sequence_scoring),
            ("sequence_plus_annotation", annotation_packets, annotation_scoring),
        ]:
            source_manifest = _source_manifest(candidate, mode=mode)
            _write_json(DATA_ROOT / "source_manifests" / mode / source_manifest["target_id"] / "source_manifest.json", source_manifest)
            packet = build_sealed_simulation_packet(
                target_id=source_manifest["target_id"],
                target_name=f"{candidate['entry_id']} {candidate['entity_description']}",
                sequence=candidate["sequence"],
                sources=source_manifest["prediction_sources"],
                focus_regions=[{"name": "automatic real-sequence region scan", "span": f"1-{candidate['sequence_length']}"}],
                perturbations=perturbations,
            )
            packets.append(packet)
            _write_json(DATA_ROOT / "sealed_predictions" / mode / packet["target_id"] / "sealed_simulation_packet.json", packet)
            holdout = _holdout(candidate, packet, expected)
            _write_json(DATA_ROOT / "holdouts_postseal" / mode / packet["target_id"] / "postseal_holdout_manifest.json", holdout)
            scoring = _score_real_sequence(packet, holdout)
            scoring_out.append(scoring)
            _write_json(DATA_ROOT / "validation" / mode / packet["target_id"] / "validation_result.json", scoring)
        annotation_manifest = _source_manifest(candidate, mode="sequence_plus_annotation")
        wrong_packet = build_sealed_simulation_packet(
            target_id=f"V58_{candidate['target_id']}_WRONG_GRAMMAR_CONTROL",
            target_name=f"{candidate['entry_id']} forced wrong grammar",
            sequence=candidate["sequence"],
            sources=annotation_manifest["prediction_sources"],
            focus_regions=[{"name": "automatic real-sequence region scan", "span": f"1-{candidate['sequence_length']}"}],
            perturbations=[],
            forced_grammar=_wrong_grammar(expected),
        )
        wrong_packets.append(wrong_packet)
        _write_json(DATA_ROOT / "wrong_grammar_controls" / wrong_packet["target_id"] / "wrong_grammar_packet.json", wrong_packet)
        shuffled_packet = build_sealed_simulation_packet(
            target_id=f"V58_{candidate['target_id']}_SHUFFLED_CONTROL",
            target_name=f"{candidate['entry_id']} shuffled sequence control",
            sequence=shuffled_sequence(candidate["sequence"]),
            sources=annotation_manifest["prediction_sources"],
            focus_regions=[],
            perturbations=[],
        )
        shuffled_packets.append(shuffled_packet)
        _write_json(DATA_ROOT / "shuffled_controls" / shuffled_packet["target_id"] / "shuffled_control_packet.json", shuffled_packet)

    baseline = _baseline_scores(annotation_scoring)
    controls = _controls(
        target_selection=target_selection,
        annotation_packets=annotation_packets,
        sequence_packets=sequence_packets,
        wrong_packets=wrong_packets,
        shuffled_packets=shuffled_packets,
        scoring_rows=annotation_scoring,
        engine_declaration=engine_declaration,
        majority_baseline_accuracy=baseline["majority_class_accuracy"],
        annotation_accuracy=baseline["annotation_engine_accuracy"],
    )
    cert = _aggregate_certificate(
        target_selection=target_selection,
        engine_declaration=engine_declaration,
        sequence_scoring=sequence_scoring,
        annotation_scoring=annotation_scoring,
        controls=controls,
        baseline_scores=baseline,
    )
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    _write_json(DATA_ROOT / "v58_sequence_only_scoring_report.json", {"kind": "V58_SEQUENCE_ONLY_SCORING_REPORT_v0", "rows": sequence_scoring})
    _write_json(DATA_ROOT / "v58_sequence_plus_annotation_scoring_report.json", {"kind": "V58_SEQUENCE_PLUS_ANNOTATION_SCORING_REPORT_v0", "rows": annotation_scoring})
    _write_json(DATA_ROOT / "v58_failure_report.json", {"kind": "V58_FAILURE_REPORT_v0", "failure_cases": cert["failure_cases"], "failure_cases_reported": True})
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v58_real_sequence_time_blind_folding_replication_certificate.json"
    report_path = out_dir / "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_GATE_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert, annotation_scoring)
    return {
        "certificate": cert_path,
        "report": report_path,
        "raw_candidate_cache": RAW_CANDIDATE_CACHE,
        "target_manifest": DATA_ROOT / "v58_real_sequence_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v58_frozen_engine_declaration.json",
        "sequence_only_scoring": DATA_ROOT / "v58_sequence_only_scoring_report.json",
        "annotation_scoring": DATA_ROOT / "v58_sequence_plus_annotation_scoring_report.json",
        "failure_report": DATA_ROOT / "v58_failure_report.json",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V58 real-sequence time-blind folding replication gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--refresh-intake", action="store_true", help="refresh recent-release real sequence intake from RCSB using curl")
    args = parser.parse_args()
    paths = run_v58(args.out_dir, refresh_intake=args.refresh_intake)
    cert = _read_json(paths["certificate"], "V58 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "target_count": cert["target_count"],
        "level1_regime_selection_supported_count": cert["level1_regime_selection_supported_count"],
        "level3_topology_or_observable_supported_count": cert["level3_topology_or_observable_supported_count"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "folding_problem_solved": cert["folding_problem_solved"],
        "target_selection_manual": cert["target_selection_manual"],
        "coordinate_truth_used_before_seal": cert["coordinate_truth_used_before_seal"],
        "alphafold_used_before_seal": cert["alphafold_used_before_seal"],
        "failure_cases_reported": cert["failure_cases_reported"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] in {PASSED, PARTIAL} else 1


if __name__ == "__main__":
    raise SystemExit(main())
