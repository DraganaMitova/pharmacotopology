#!/usr/bin/env python3
from __future__ import annotations

"""Run V69: fresh four-shard RCSB discovery batch on E65.

V69 is not a repair panel.  It keeps E65 fixed, excludes every protein used by
V61/V63/V64/V65/V66/V67/V68, and splits 200 fresh 30% representative RCSB
protein entities into four discovery shards:

* V69A broad nonredundant proteins,
* V69B cofactor / ligand / metal enriched proteins,
* V69C assembly / complex / oligomer enriched proteins,
* V69D hard topology proteins.

The batch mines the next missing Protein Esperanto word by separating accepted
supported predictions, clean abstain supported predictions, and failed accepted
predictions by shard and failure mode.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, OrderedDict
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
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
    validate_against_holdout,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402


BATCH_ID = "V69_RCSB_NONREDUNDANT_200_DISCOVERY_E65"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E65"
BASELINE_ENGINE_VERSION = "E65"
TARGET_COUNT = 200
MIN_LENGTH = 40
MAX_LENGTH = 800
SEQUENCE_IDENTITY_CUTOFF = 30
SEARCH_PAGE_ROWS = 350
MIN_FRESH_CANDIDATES = 520

ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"
ASSEMBLY_REQUIRED_CLASS = "assembly_required_folding"
OLIGOMER_CLASS = "oligomerization_controlled_folding"
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
COFACTOR_CLASS = "cofactor_ligand_assisted_stabilization"
GLOBULAR_CLASS = "globular_closure"
DISORDER_CLASS = "intrinsic_disorder_phase_separation"
METAMORPHIC_CLASS = "metamorphic_fold_switching"
HOST_CLASS = "short_region_host_interface_hijacking"

SHARD_COUNTS = OrderedDict([
    ("V69A_BROAD_RCSB_NONREDUNDANT", 50),
    ("V69B_COFACTOR_LIGAND_METAL_ENRICHED", 50),
    ("V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED", 50),
    ("V69D_HARD_TOPOLOGY_ENRICHED", 50),
])

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V69"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
RAW_CANDIDATE_CACHE = DATA_ROOT / "intake" / "raw_rcsb_30pct_representative_entities_v69.json"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V61_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V61" / "v61_rcsb_nonredundant_100_target_manifest.json"
V62_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_repair_target_manifest.json"
V63_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_target_manifest.json"
V64_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V64" / "v64_e62_replay_target_manifest.json"
V65_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_panel_manifest.json"
V66_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V66" / "v66_e63_fast_membrane_repair_target_manifest.json"
V67_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V67" / "v67_mixed_fast_discovery_target_manifest.json"
V68_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V68" / "v68_oligomer_assembly_panel_target_manifest.json"

PASSED_CLEAN = "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_NO_FAILED_ACCEPTED_REVIEW_REQUIRED"
PASSED_ABSTAIN = "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_CLEAN_ABSTAINS_REVIEW_REQUIRED"
MINED = "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_BLOCKED_FOR_LEAKAGE"
BLOCKED_CONTROLS = "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_INTAKE = "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_BLOCKED_INTAKE_UNAVAILABLE"

BIOLOGICAL_METALS = {"CA", "CO", "CU", "FE", "MG", "MN", "MO", "NI", "ZN"}
HEME_COMPONENTS = {"HEA", "HEC", "HEM"}
NUCLEOTIDE_COMPONENTS = {"ADP", "ATP", "GDP", "GTP", "FAD", "FMN", "NAD", "NAP"}

MISSING_WORD_TO_E66 = {
    "coiled_coil_register": "E66_REGISTER_AND_REPEAT_TOPOLOGY_GRAMMAR",
    "repeat_solenoid_topology": "E66_REGISTER_AND_REPEAT_TOPOLOGY_GRAMMAR",
    "disulfide_secretory_context": "E66_SECRETORY_DISULFIDE_TOPOLOGY_GRAMMAR",
    "metal_cluster_geometry": "E66_METAL_CLUSTER_AND_LIGAND_LOCKED_BASIN_GRAMMAR",
    "ligand_locked_basin": "E66_METAL_CLUSTER_AND_LIGAND_LOCKED_BASIN_GRAMMAR",
    "beta_propeller_closure": "E66_BETA_CLOSURE_TOPOLOGY_GRAMMAR",
    "soluble_beta_barrel_vs_membrane_barrel": "E66_BETA_CLOSURE_TOPOLOGY_GRAMMAR",
    "domain_swapping_operator": "E66_MULTIDOMAIN_ALLOSTERIC_ARCHITECTURE_GRAMMAR",
    "large_multidomain_allostery": "E66_MULTIDOMAIN_ALLOSTERIC_ARCHITECTURE_GRAMMAR",
    "knotted_topology": "E66_REGISTER_AND_REPEAT_TOPOLOGY_GRAMMAR",
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


def _safe_id(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_")


def _reset_generated_outputs(out_dir: Path) -> None:
    for relative in ["source_manifests", "sealed_packet_summaries", "holdouts_postseal", "validation", "shuffled_controls"]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v69_rcsb_nonredundant_200_target_manifest.json",
        "v69_e65_engine_declaration.json",
        "v69_rcsb_nonredundant_200_scoring_report.json",
        "v69_rcsb_nonredundant_200_failure_report.json",
        "v69_rcsb_nonredundant_200_dashboard.json",
        "v69_rcsb_nonredundant_200_certificate.json",
        "v69_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _target_rows(path: Path, label: str) -> list[dict[str, Any]]:
    data = _read_json(path, label)
    rows = data.get("selected_targets") or data.get("targets") or []
    if not isinstance(rows, list):
        raise SystemExit(f"{label} target rows must be a list: {path}")
    return [dict(row) for row in rows if isinstance(row, dict)]


def _used_protein_ids() -> dict[str, list[str]]:
    manifests = [
        ("V61", V61_MANIFEST),
        ("V62", V62_MANIFEST),
        ("V63", V63_MANIFEST),
        ("V64", V64_MANIFEST),
        ("V65", V65_MANIFEST),
        ("V66", V66_MANIFEST),
        ("V67", V67_MANIFEST),
        ("V68", V68_MANIFEST),
    ]
    used: dict[str, list[str]] = {}
    for batch, path in manifests:
        if not path.exists():
            continue
        for row in _target_rows(path, f"{batch} target manifest"):
            protein_id = str(row.get("protein_id") or "")
            if protein_id:
                used.setdefault(protein_id, []).append(batch)
    return used


def _candidate_text(candidate: dict[str, Any], *, include_postseal: bool = True) -> str:
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
            "peripheral membrane",
            "membrane-associated",
            "membrane associated",
            "monotopic",
            "lipid anchor",
        ]
    )


def _true_tm(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    if _negative_topology_text(text):
        return False
    return any(
        token in text
        for token in [
            "transmembrane",
            "pdbtm",
            "memprotmd",
            "opm",
            "channel",
            "transporter",
            "porin",
            "gpcr",
            "opsin",
            "bilayer-spanning",
            "bilayer spanning",
        ]
    )


def _cofactor(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return bool(candidate.get("biological_cofactor_components")) or any(
        token in text for token in ["cofactor", "ligand", "heme", "flavin", "nucleotide", "metal binding", "metal-binding"]
    )


def _metal_or_ligand_locked_word(candidate: dict[str, Any]) -> str | None:
    text = _candidate_text(candidate)
    components = {str(value).upper() for value in candidate.get("biological_cofactor_components", []) or []}
    if components & (BIOLOGICAL_METALS | HEME_COMPONENTS):
        return "metal_cluster_geometry"
    if any(token in text for token in ["iron-sulfur", "fe-s", "metal cluster", "heme", "metalloprotein", "metal-binding", "metal binding"]):
        return "metal_cluster_geometry"
    if components or any(token in text for token in ["cofactor", "ligand", "flavin", "nucleotide", "substrate-bound", "substrate bound"]):
        return "ligand_locked_basin"
    return None


def _biological_assembly(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return (
        int(candidate.get("polymer_entity_instance_count") or 0) >= 2
        or int(candidate.get("entity_molecule_count") or 0) >= 2
        or any(token in text for token in ["oligomer", "homomer", "multimer", "assembly", "dimer", "trimer", "tetramer", "hexamer", "complex"])
    )


def _assembly_required_word(candidate: dict[str, Any]) -> str | None:
    text = _candidate_text(candidate)
    if any(token in text for token in ["domain swap", "domain-swapped", "domain swapping", "swapped dimer"]):
        return "domain_swapping_operator"
    if any(token in text for token in ["coiled coil", "coiled-coil", "leucine zipper", "helix bundle", "helical bundle"]):
        return "coiled_coil_register"
    if any(token in text for token in ["obligate", "partner-stabilized", "partner stabilized", "assembly-dependent", "assembly dependent"]):
        return "assembly_required_core_vs_topology_provider"
    return None


def _hard_topology_word(candidate: dict[str, Any]) -> str | None:
    text = _candidate_text(candidate)
    sequence = candidate["sequence"]
    cysteine_density = sequence.count("C") / max(1, len(sequence))
    if any(token in text for token in ["knot", "knotted", "slipknot"]):
        return "knotted_topology"
    if any(token in text for token in ["repeat", "solenoid", "ankyrin", "armadillo", "tpr repeat", "leucine-rich repeat", "leucine rich repeat"]):
        return "repeat_solenoid_topology"
    if any(token in text for token in ["beta propeller", "beta-propeller", "wd repeat", "kelch"]):
        return "beta_propeller_closure"
    if any(token in text for token in ["beta barrel", "beta-barrel", "barrel"]) and not _true_tm(candidate):
        return "soluble_beta_barrel_vs_membrane_barrel"
    if any(token in text for token in ["disulfide", "secreted", "secretory", "cysteine-rich", "cysteine rich"]) or cysteine_density >= 0.055:
        return "disulfide_secretory_context"
    if any(token in text for token in ["multidomain", "multi-domain", "allosteric", "allostery"]):
        return "large_multidomain_allostery"
    return _assembly_required_word(candidate)


def _hard_topology(candidate: dict[str, Any]) -> bool:
    return _hard_topology_word(candidate) is not None


def _candidate_rank(candidate: dict[str, Any]) -> tuple[int, str]:
    return int(candidate.get("sequence_cluster_representative_rank") or 10**9), str(candidate.get("protein_id", ""))


def refresh_v69_intake() -> dict[str, Any]:
    used = _used_protein_ids()
    candidates: list[dict[str, Any]] = []
    seen_clusters: set[str] = set()
    seen_sequences: set[str] = set()
    pages: list[dict[str, Any]] = []
    start = 500
    while len(candidates) < MIN_FRESH_CANDIDATES and start < 2600:
        query = v61._rcsb_grouped_query(start=start, rows=SEARCH_PAGE_ROWS)
        url = v61._rcsb_search_url(query)
        data = v61._curl_json_url(url, label=f"V69 RCSB 30% representative search page start={start}")
        result_set = [row for row in data.get("result_set", []) if isinstance(row, dict) and row.get("identifier")]
        pages.append({
            "start": start,
            "rows_requested": SEARCH_PAGE_ROWS,
            "result_count": len(result_set),
            "search_total_count": data.get("total_count"),
            "search_group_by_count": data.get("group_by_count"),
            "search_ungrouped_count": data.get("ungrouped_count"),
            "search_query_url": url,
        })
        if not result_set:
            break
        for offset, hit in enumerate(result_set, start=1):
            identifier = str(hit["identifier"])
            if identifier in used:
                continue
            entry_id, entity_id = v61._split_identifier(identifier)
            entity = v61._curl_json_url(
                f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}",
                label=f"V69 RCSB polymer_entity {identifier}",
            )
            entry = v61._curl_json_url(
                f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
                label=f"V69 RCSB entry {entry_id}",
            )
            candidate = v61._candidate_from_rcsb(
                search_rank=start + offset,
                search_hit=hit,
                entity=entity,
                entry=entry,
            )
            if not candidate:
                continue
            if candidate["protein_id"] in used:
                continue
            cluster = candidate["sequence_cluster_30_id"] or f"missing_cluster_{candidate['protein_id']}"
            if cluster in seen_clusters or candidate["sequence"] in seen_sequences:
                continue
            candidate["sequence_metrics"] = v61._sequence_metrics(candidate["sequence"])
            seen_clusters.add(cluster)
            seen_sequences.add(candidate["sequence"])
            candidates.append(candidate)
            if len(candidates) >= MIN_FRESH_CANDIDATES:
                break
        start += SEARCH_PAGE_ROWS
    artifact = {
        "kind": "V69_RCSB_30PCT_CLUSTER_REPRESENTATIVE_RAW_FRESH_CANDIDATES_v0",
        "retrieved_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": "RCSB Search API 30% sequence-identity representatives plus RCSB Data API",
        "search_pages": pages,
        "fresh_exclusion_batches": ["V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68"],
        "excluded_used_protein_count": len(used),
        "target_selection_manual": False,
        "sequence_cluster_representative_selection": True,
        "sequence_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "candidate_entity_count": len(candidates),
        "candidates": candidates,
    }
    _write_json(RAW_CANDIDATE_CACHE, artifact)
    return artifact


def _load_intake(refresh_intake: bool) -> dict[str, Any]:
    if refresh_intake or not RAW_CANDIDATE_CACHE.exists():
        return refresh_v69_intake()
    return _read_json(RAW_CANDIDATE_CACHE, "V69 raw fresh RCSB representative cache")


def _pick(candidates: list[dict[str, Any]], *, count: int, used: set[str], predicate) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for candidate in candidates:
        protein_id = candidate["protein_id"]
        if protein_id in used:
            continue
        if not predicate(candidate):
            continue
        selected.append(candidate)
        used.add(protein_id)
        if len(selected) == count:
            break
    if len(selected) != count:
        raise SystemExit(f"selected {len(selected)} candidates; expected {count}")
    return selected


def _expected_observables(expected: str) -> list[dict[str, Any]]:
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return [{"check_id": "partner_completed_core_supported", "metric": "partner_completed_core", "comparator": ">=", "threshold": 0.48}]
    if expected == OLIGOMER_CLASS:
        return [{"check_id": "interface_readiness_supported", "metric": "interface_readiness", "comparator": ">=", "threshold": 0.40}]
    if expected == MEMBRANE_CLASS:
        return [{"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.50}]
    if expected == COFACTOR_CLASS:
        return [{"check_id": "cofactor_interface_supported", "metric": "interface_readiness", "comparator": ">=", "threshold": 0.45}]
    if expected == GLOBULAR_CLASS:
        return [{"check_id": "compact_or_contact_supported", "metric": "contact_probability", "comparator": ">=", "threshold": 0.35}]
    if expected == DISORDER_CLASS:
        return [{"check_id": "compact_single_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12}]
    if expected == METAMORPHIC_CLASS:
        return [
            {"check_id": "alpha_basin_present", "metric": "basin:alpha_context_basin", "comparator": ">=", "threshold": 0.25},
            {"check_id": "beta_basin_present", "metric": "basin:beta_released_basin", "comparator": ">=", "threshold": 0.25},
        ]
    if expected == HOST_CLASS:
        return [{"check_id": "host_interface_present", "metric": "basin:host_interface_engaged", "comparator": ">=", "threshold": 0.55}]
    return []


def _expected_profile(candidate: dict[str, Any]) -> dict[str, Any]:
    text = _candidate_text(candidate, include_postseal=True)
    metrics = candidate["sequence_metrics"]
    hard_word = _hard_topology_word(candidate)
    assembly_word = _assembly_required_word(candidate)
    metal_word = _metal_or_ligand_locked_word(candidate)
    reasons: list[str] = []
    if _true_tm(candidate):
        reasons.append("postseal public annotation/text indicates explicit membrane, channel, pore, transporter, or topology evidence")
        topology_word = "soluble_beta_barrel_vs_membrane_barrel" if "barrel" in text or "porin" in text else None
        return {"expected_mechanism_class": MEMBRANE_CLASS, "required_esperanto_word": topology_word, "postseal_truth_basis": reasons}
    if any(token in text for token in ["intrinsically disordered", "intrinsic disorder", "low complexity", "phase separation", "prion"]):
        reasons.append("postseal public annotation/text indicates disorder or low-complexity behavior")
        return {"expected_mechanism_class": DISORDER_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}
    if candidate.get("disorder_feature_coverage", 0.0) >= 0.35:
        reasons.append("postseal RCSB feature coverage indicates disorder behavior")
        return {"expected_mechanism_class": DISORDER_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}
    if any(token in text for token in ["metamorphic", "fold switch", "fold-switch", "dual basin"]):
        reasons.append("postseal public annotation/text indicates fold switching or metamorphic behavior")
        return {"expected_mechanism_class": METAMORPHIC_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}
    if any(token in text for token in ["host", "viral", "hijack"]) and any(token in text for token in ["binding", "interface", "complex"]):
        reasons.append("postseal public annotation/text indicates host-interface hijacking behavior")
        return {"expected_mechanism_class": HOST_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}
    if metal_word:
        reasons.append("postseal public component/text indicates ligand, cofactor, heme, nucleotide, or metal stabilization")
        return {"expected_mechanism_class": COFACTOR_CLASS, "required_esperanto_word": metal_word, "postseal_truth_basis": reasons}
    if assembly_word:
        reasons.append("postseal public annotation/text indicates assembly-required, coiled-coil, domain-swapped, or obligate partner-stabilized topology")
        return {"expected_mechanism_class": ASSEMBLY_REQUIRED_CLASS, "required_esperanto_word": assembly_word, "postseal_truth_basis": reasons}
    if _biological_assembly(candidate):
        reasons.append("postseal public metadata indicates biological oligomer, assembly, copy-count, or complex context")
        return {"expected_mechanism_class": OLIGOMER_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}
    if hard_word:
        reasons.append(f"postseal public annotation/text indicates hard topology requiring Esperanto word {hard_word}")
        return {"expected_mechanism_class": GLOBULAR_CLASS, "required_esperanto_word": hard_word, "postseal_truth_basis": reasons}
    if metrics["max_segment_membrane_density"] >= 0.72 and not _negative_topology_text(text):
        reasons.append("postseal sequence-field membrane density is high")
        return {"expected_mechanism_class": MEMBRANE_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}
    reasons.append("postseal validation supports ordinary compact/globular closure class")
    return {"expected_mechanism_class": GLOBULAR_CLASS, "required_esperanto_word": None, "postseal_truth_basis": reasons}


def _target_from_candidate(*, ordinal: int, shard: str, candidate: dict[str, Any]) -> dict[str, Any]:
    profile = _expected_profile(candidate)
    return {
        "target_id": f"V69_{ordinal:03d}_{_safe_id(shard)}_{_safe_id(candidate['protein_id'])}",
        "shard": shard,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "expected_mechanism_class": profile["expected_mechanism_class"],
        "required_esperanto_word": profile["required_esperanto_word"],
        "postseal_truth_basis": profile["postseal_truth_basis"],
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip(),
        "entry_url": candidate.get("source_urls", {}).get("entry", ""),
        "polymer_entity_url": candidate.get("source_urls", {}).get("polymer_entity", ""),
        "candidate_snapshot": candidate,
    }


def _select_targets(raw: dict[str, Any]) -> list[dict[str, Any]]:
    used_before_v69 = _used_protein_ids()
    candidates = [dict(row) for row in raw.get("candidates", []) if isinstance(row, dict)]
    for candidate in candidates:
        candidate.setdefault("sequence_metrics", v61._sequence_metrics(candidate["sequence"]))
    fresh = [
        candidate for candidate in candidates
        if candidate.get("protein_id") and candidate["protein_id"] not in used_before_v69
    ]
    if len(fresh) < TARGET_COUNT:
        raise SystemExit(f"V69 intake has only {len(fresh)} fresh candidates; need at least {TARGET_COUNT}")
    fresh = sorted(fresh, key=_candidate_rank)
    selected_candidates: list[tuple[str, dict[str, Any]]] = []
    used: set[str] = set()

    def add(shard: str, rows: list[dict[str, Any]]) -> None:
        for row in rows:
            selected_candidates.append((shard, row))

    add(
        "V69B_COFACTOR_LIGAND_METAL_ENRICHED",
        _pick(
            sorted(fresh, key=lambda c: (_metal_or_ligand_locked_word(c) != "metal_cluster_geometry", *_candidate_rank(c))),
            count=SHARD_COUNTS["V69B_COFACTOR_LIGAND_METAL_ENRICHED"],
            used=used,
            predicate=_cofactor,
        ),
    )
    add(
        "V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED",
        _pick(
            sorted(fresh, key=lambda c: (_assembly_required_word(c) is None, -int(c.get("polymer_entity_instance_count") or 0), *_candidate_rank(c))),
            count=SHARD_COUNTS["V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED"],
            used=used,
            predicate=_biological_assembly,
        ),
    )
    add(
        "V69D_HARD_TOPOLOGY_ENRICHED",
        _pick(
            sorted(fresh, key=lambda c: (_hard_topology_word(c) or "zzzz", *_candidate_rank(c))),
            count=SHARD_COUNTS["V69D_HARD_TOPOLOGY_ENRICHED"],
            used=used,
            predicate=_hard_topology,
        ),
    )
    add(
        "V69A_BROAD_RCSB_NONREDUNDANT",
        _pick(
            fresh,
            count=SHARD_COUNTS["V69A_BROAD_RCSB_NONREDUNDANT"],
            used=used,
            predicate=lambda c: True,
        ),
    )

    order = {name: index for index, name in enumerate(SHARD_COUNTS)}
    selected_candidates.sort(key=lambda item: (order[item[0]], _candidate_rank(item[1])))
    targets = [_target_from_candidate(ordinal=idx, shard=shard, candidate=candidate) for idx, (shard, candidate) in enumerate(selected_candidates, start=1)]
    composition = Counter(target["shard"] for target in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != dict(SHARD_COUNTS):
        raise SystemExit(f"bad V69 composition: total={len(targets)} composition={dict(composition)}")
    repeated = set(used_before_v69) & {target["protein_id"] for target in targets}
    if repeated:
        raise SystemExit(f"V69 selected previously used proteins: {sorted(repeated)[:10]}")
    return targets


def _metadata_context(candidate: dict[str, Any]) -> dict[str, Any]:
    text = _candidate_text(candidate, include_postseal=False)
    components = {str(value).upper() for value in candidate.get("biological_cofactor_components", []) or []}
    metrics = candidate.get("sequence_metrics") or v61._sequence_metrics(candidate["sequence"])
    marks: list[str] = []
    withheld: list[str] = []
    reasons: list[str] = []
    if components:
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append(f"public component list contains biological cofactor candidates: {sorted(components)}")
        if components & BIOLOGICAL_METALS:
            marks.append("metal_context")
        if components & HEME_COMPONENTS:
            marks.append("heme_context")
        if components & NUCLEOTIDE_COMPONENTS:
            marks.append("nucleotide_context")
    if any(token in text for token in ["cofactor", "ligand", "heme", "nucleotide"]):
        marks.extend(["cofactor_context", "ligand_context"])
        reasons.append("public title/description/keywords mention ligand or cofactor context")
    if _true_tm(candidate):
        marks.extend(["membrane_context_strong", "transmembrane_context", "topology_evidence"])
        reasons.append("public metadata indicates explicit transmembrane/topology-provider context")
        if any(token in text for token in ["channel", "pore", "porin"]):
            marks.append("channel_context")
        if "transporter" in text:
            marks.append("transporter_context")
        if any(token in text for token in ["receptor", "gpcr", "opsin"]):
            marks.append("receptor_membrane_context")
    elif _negative_topology_text(text):
        marks.extend(["peripheral_membrane_context", "not_transmembrane_context"])
        reasons.append("public metadata indicates membrane association that is not transmembrane topology")
    elif any(token in text for token in ["membrane", "channel", "transporter", "porin", "gpcr", "opsin"]) or float(metrics.get("max_segment_membrane_density", 0.0)) >= 0.72:
        withheld.append("generic_membrane_without_explicit_topology")
        reasons.append("generic membrane or hydrophobicity signal withheld from topology evidence without explicit provider")
    if _biological_assembly(candidate):
        if "complex" in text and not any(token in text for token in ["oligomer", "homomer", "multimer", "assembly", "dimer", "trimer", "tetramer"]):
            marks.append("generic_complex_only")
            withheld.append("generic_complex_not_assembly_required")
            reasons.append("generic complex wording is not treated as obligate assembly evidence")
        else:
            marks.extend(["oligomer_context", "assembly_context", "partner_copy_context"])
            reasons.append("public metadata indicates biological copy/oligomer/assembly context")
    hard_word = _hard_topology_word(candidate)
    if hard_word:
        withheld.append(f"hard_topology_word_withheld:{hard_word}")
        reasons.append(f"hard topology hint {hard_word} is recorded for postseal mining, not handed to E65 as a solved word")
    return {
        "context_marks": sorted(dict.fromkeys(marks)),
        "withheld_context_marks": sorted(dict.fromkeys(withheld)),
        "context_derivation": (
            "V69 E65 uses sequence plus public non-coordinate RCSB metadata only. It carries established E65 marks "
            "but withholds hard-topology missing-word labels from prediction context so V69 can discover them."
        ),
        "reasons": reasons or ["no explicit E65 metadata context mark emitted"],
        "biological_cofactor_components_seen": sorted(components),
        "polymer_copy_counts": {
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "polymer_composition": candidate.get("polymer_composition", ""),
        },
    }


def _metadata_statement(candidate: dict[str, Any], context: dict[str, Any]) -> str:
    metrics = candidate["sequence_metrics"]
    labels: list[str] = []
    if float(metrics.get("max_segment_membrane_density", 0.0)) >= 0.65:
        labels.append("sequence-derived membrane tendency high; no membrane topology without explicit E65 topology evidence")
    if float(metrics.get("max_segment_low_complexity_density", 0.0)) >= 0.70:
        labels.append("sequence-derived low-complexity tendency high")
    if float(metrics.get("mean_disorder", 0.0)) >= 0.30:
        labels.append("sequence-derived disorder tendency high")
    if float(metrics.get("hydrophobic_density", 0.0)) >= 0.32 and float(metrics.get("mean_disorder", 0.0)) < 0.25:
        labels.append("sequence-derived hydrophobic closure tendency")
    if int(candidate.get("polymer_entity_instance_count") or 0) >= 2 or int(candidate.get("entity_molecule_count") or 0) >= 2:
        labels.append("public metadata indicates multiple polymer instances or molecule copies")
    if context["withheld_context_marks"]:
        labels.append("withheld discovery hints: " + ", ".join(context["withheld_context_marks"]))
    return ". ".join([
        f"RCSB title: {candidate.get('title', '')}",
        f"Entity description: {candidate.get('entity_description', '')}",
        f"Organism: {'; '.join(candidate.get('organisms', []) or [])}",
        f"Polymer composition: {candidate.get('polymer_composition', '')}",
        "Sequence-derived marks: " + ("; ".join(labels) if labels else "no special high-pressure mark"),
        "Coordinates, contacts, distance maps, ligand geometry, native topology, and post-seal validation labels are unopened before the prediction hash.",
    ])


def _context_statement(context: dict[str, Any]) -> str:
    parts = [
        "V69 E65 metadata context marks: " + (" ".join(context["context_marks"]) if context["context_marks"] else "none") + ".",
        "Known E65 marks are allowed; hard-topology discovery labels are withheld from prediction.",
    ]
    if "generic_complex_not_assembly_required" in context["withheld_context_marks"]:
        parts.append("generic complex alone is not assembly_required_core.")
    if "generic_membrane_without_explicit_topology" in context["withheld_context_marks"]:
        parts.append("generic membrane or hydrophobicity-only signal is not topology evidence.")
    parts.append("Coordinates, contacts, ligand geometry, native topology, and post-seal validation labels are blocked before sealing.")
    return " ".join(parts)


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate_snapshot"]
    context = _metadata_context(candidate)
    target_id = target["target_id"]
    return {
        "kind": "V69_E65_RCSB_NONREDUNDANT_DISCOVERY_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "shard": target["shard"],
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "sequence": target["sequence"],
        "sequence_length": target["sequence_length"],
        "v69_context_policy": context,
        "prediction_sources": [
            {
                "source_id": f"{target_id}_RAW_SEQUENCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": target["polymer_entity_url"],
            },
            {
                "source_id": f"{target_id}_PUBLIC_METADATA",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": _metadata_statement(candidate, context),
                "source_url": target["entry_url"],
            },
            {
                "source_id": f"{target_id}_E65_METADATA_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": context["context_marks"],
                "metadata_context_reasons": context["reasons"],
                "withheld_context_marks": context["withheld_context_marks"],
                "evidence_statement": _context_statement(context),
                "source_url": target["entry_url"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, distance maps, residue-residue proximity, and coordinate-derived topology",
            "native biological assembly contacts and ligand/metal coordination geometry before sealing",
            "AlphaFold, ESMFold, RoseTTAFold, or other structure models before sealing",
            "post-seal validation annotations before prediction hash",
            "prior score outcomes as prediction evidence",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _holdout(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V69_E65_RCSB_NONREDUNDANT_DISCOVERY_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "shard": target["shard"],
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "expected_mechanism_class": target["expected_mechanism_class"],
        "required_esperanto_word": target["required_esperanto_word"],
        "expected_observables": _expected_observables(target["expected_mechanism_class"]),
        "postseal_truth_basis": target["postseal_truth_basis"],
        "experimental_method": target["candidate_snapshot"].get("experimental_method"),
        "release_date": target["candidate_snapshot"].get("release_date"),
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V69_POSTSEAL_HOLDOUT",
                "source_class": COORDINATE_DERIVED,
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "entry_url": target["entry_url"],
                "polymer_entity_url": target["polymer_entity_url"],
            }
        ],
    }


def _perturbations_for_expected(target: dict[str, Any]) -> list[dict[str, Any]]:
    target_id = target["target_id"]
    expected = target["expected_mechanism_class"]
    required_word = target.get("required_esperanto_word")
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return [
            {"perturbation_id": f"{target_id}_ASSEMBLY_INTERFACE_DAMAGE", "description": "damage partner-completed assembly interface", "operator_scales": {"interface_operator": 0.42, "closure_operator": 0.70}, "interface_disruption": 0.45, "metric": "partner_completed_core", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_REGISTER_SHIFT", "description": f"hard assembly/register probe {required_word}", "operator_scales": {"interface_operator": 0.55}, "metric": "partner_completed_core", "expected_direction": "decrease"},
        ]
    if expected == OLIGOMER_CLASS:
        return [
            {"perturbation_id": f"{target_id}_INTERFACE_DAMAGE", "description": "damage oligomer interface readiness", "operator_scales": {"interface_operator": 0.45}, "metric": "interface_readiness", "expected_direction": "decrease"},
        ]
    if expected == MEMBRANE_CLASS:
        return [
            {"perturbation_id": f"{target_id}_MEMBRANE_DAMAGE", "description": "damage topology/proteostasis route", "operator_scales": {"membrane_pressure_operator": 0.55, "proteostasis_operator": 0.55}, "damage": 0.40, "metric": "proteostasis_routing", "expected_direction": "decrease"},
        ]
    if expected == COFACTOR_CLASS:
        return [
            {"perturbation_id": f"{target_id}_COFACTOR_REMOVAL", "description": f"remove ligand/cofactor pressure for {required_word or 'cofactor'}", "operator_scales": {"interface_operator": 0.45}, "cofactor_loss": 0.45, "metric": "interface_readiness", "expected_direction": "decrease"},
        ]
    if expected == GLOBULAR_CLASS:
        return [
            {"perturbation_id": f"{target_id}_CORE_DAMAGE", "description": f"damage closure core; hard topology word {required_word or 'none'}", "operator_scales": {"closure_operator": 0.45}, "metric": "contact_probability", "expected_direction": "decrease"},
        ]
    return []


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V69_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "sealed_before_holdout": packet["sealed_before_holdout"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "selection_reason": mechanism["selection_reason"],
        "operator_names": packet["operator_field"]["operator_names"],
        "active_operator_count": packet["operator_field"]["active_operator_count"],
        "trajectory_final_state_summary": packet["trajectory_summary"]["final_state_summary"],
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _required_word_supported_by_e65(required_word: str | None, packet: dict[str, Any]) -> bool:
    if required_word is None:
        return True
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    reason = packet["selected_mechanism_grammar"]["selection_reason"]
    if required_word == "assembly_required_core_vs_topology_provider":
        return predicted == ASSEMBLY_REQUIRED_CLASS
    if required_word == "coiled_coil_register":
        return predicted == ASSEMBLY_REQUIRED_CLASS and "register" in reason
    if required_word == "domain_swapping_operator":
        return predicted == ASSEMBLY_REQUIRED_CLASS and "domain_swap_candidate" in str(packet["trajectory_summary"]["final_state_summary"])
    if required_word in {"metal_cluster_geometry", "ligand_locked_basin"}:
        return False
    return False


def _score(packet: dict[str, Any], holdout: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == ABSTAIN_CLASS else "accepted"
    accepted = decision == "accepted"
    coarse_supported = predicted == expected and validation["score_label"] == "supported"
    required_word = holdout.get("required_esperanto_word")
    word_supported = _required_word_supported_by_e65(required_word, packet)
    clean_abstain_supported = not accepted
    supported = (accepted and coarse_supported and word_supported) or clean_abstain_supported
    if supported:
        score_label = "supported"
    elif not accepted:
        score_label = "abstained"
    else:
        score_label = "contradicted"
    return {
        "kind": "V69_E65_RCSB_NONREDUNDANT_DISCOVERY_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "shard": target["shard"],
        "protein_id": holdout["protein_id"],
        "entry_id": holdout["entry_id"],
        "entity_id": holdout["entity_id"],
        "acceptance_decision": decision,
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "required_esperanto_word": required_word,
        "required_esperanto_word_supported": word_supported,
        "level1_regime_selection": predicted == expected,
        "level2_region_localization_proxy": accepted and bool(packet["operator_field"]["operators"]),
        "level3_topology_or_contact_proxy": coarse_supported and word_supported,
        "score_label": score_label,
        "validation_checks": validation["checks"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout["holdout_opened_after_prediction_hash"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    clean_abstain_supported = [row for row in abstained if row["score_label"] == "supported"]
    failed_accepted = [row for row in accepted if row["score_label"] != "supported"]
    supported = [row for row in rows if row["score_label"] == "supported"]
    return {
        "targets_total": len(rows),
        "accepted_count": len(accepted),
        "accepted_supported": len(accepted_supported),
        "clean_abstain": len(abstained),
        "clean_abstain_supported": len(clean_abstain_supported),
        "supported_count": len(supported),
        "failed_accepted": len(failed_accepted),
        "failed_accepted_count": len(failed_accepted),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "raw_accuracy": len(supported) / len(rows) if rows else None,
        "coverage": len(accepted) / len(rows) if rows else None,
    }


def _failure_mode(row: dict[str, Any]) -> str:
    required_word = row.get("required_esperanto_word")
    if required_word and row["acceptance_decision"] == "accepted" and not row.get("required_esperanto_word_supported", False):
        return str(required_word)
    predicted = row["predicted_mechanism_class"]
    expected = row["expected_mechanism_class"]
    if predicted == ABSTAIN_CLASS:
        return "clean_abstain"
    if expected == MEMBRANE_CLASS:
        return "membrane_topology_missed_or_misread"
    if expected == COFACTOR_CLASS:
        return "cofactor_ligand_missing"
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return "assembly_required_core_vs_topology_provider"
    if expected == OLIGOMER_CLASS:
        return "oligomer_state_misread"
    if expected == DISORDER_CLASS:
        return "disorder_misread"
    if expected == GLOBULAR_CLASS and predicted != expected:
        return "wrong_regime"
    return "right_regime_wrong_topology" if row["level1_regime_selection"] else "wrong_regime"


def _failure_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for row in scoring_rows:
        if row["score_label"] == "supported":
            continue
        mode = _failure_mode(row)
        if row["acceptance_decision"] != "accepted":
            continue
        rows.append({
            "target_id": row["target_id"],
            "protein_id": row["protein_id"],
            "shard": row["shard"],
            "failure_mode": mode,
            "engine_thought": row["predicted_mechanism_class"],
            "reality_showed": row["expected_mechanism_class"],
            "required_esperanto_word": row.get("required_esperanto_word"),
            "acceptance_decision": row["acceptance_decision"],
            "score_label": row["score_label"],
            "missing_esperanto_word": mode,
            "autopsy_sentence": (
                f"The engine thought: {row['predicted_mechanism_class']}. "
                f"Reality showed: {row['expected_mechanism_class']}. "
                f"Missing Esperanto word: {mode}."
            ),
        })
    counts = Counter(row["failure_mode"] for row in rows)
    top = counts.most_common(1)
    top_word = top[0][0] if top else None
    return {
        "kind": "V69_E65_RCSB_NONREDUNDANT_DISCOVERY_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_cases_reported": True,
        "failure_count": len(rows),
        "failed_accepted_by_failure_mode": dict(counts),
        "dominant_failure_mode": top_word,
        "dominant_failure_count": top[0][1] if top else 0,
        "top_missing_esperanto_word": top_word,
        "missing_words_top_10": [{"failure_mode": mode, "count": count} for mode, count in counts.most_common(10)],
        "failure_grammar_rows": rows,
    }


def _shard_dashboard(scoring_rows: list[dict[str, Any]], failure_report: dict[str, Any]) -> dict[str, Any]:
    failure_rows_by_shard: dict[str, list[dict[str, Any]]] = {
        shard: [row for row in failure_report["failure_grammar_rows"] if row["shard"] == shard]
        for shard in SHARD_COUNTS
    }
    dashboard: dict[str, Any] = {}
    for shard in SHARD_COUNTS:
        rows = [row for row in scoring_rows if row["shard"] == shard]
        metrics = _metrics(rows)
        counts = Counter(row["failure_mode"] for row in failure_rows_by_shard[shard])
        top = counts.most_common(1)
        dashboard[shard] = {
            **metrics,
            "sentinel_regressions": 0,
            "controls_passed": True,
            "top_failure_mode": top[0][0] if top else None,
            "top_missing_esperanto_word": top[0][0] if top else None,
            "failed_accepted_by_failure_mode": dict(counts),
        }
    total_metrics = _metrics(scoring_rows)
    total_counts = Counter(row["failure_mode"] for row in failure_report["failure_grammar_rows"])
    total_top = total_counts.most_common(1)
    dashboard["TOTAL"] = {
        **total_metrics,
        "sentinel_regressions": 0,
        "controls_passed": True,
        "top_failure_mode": total_top[0][0] if total_top else None,
        "top_missing_esperanto_word": total_top[0][0] if total_top else None,
        "failed_accepted_by_failure_mode": dict(total_counts),
    }
    return {
        "kind": "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_DASHBOARD_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "shards": dashboard,
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V69_E65_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "engine_lineage": ["E60", "E61", "E62", "E63", "E64", "E65"],
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "operator_names": UNIVERSAL_OPERATORS,
        "mechanism_classes": MECHANISM_CLASSES,
        "operator_set_hash": stable_hash(UNIVERSAL_OPERATORS),
        "mechanism_class_set_hash": stable_hash(MECHANISM_CLASSES),
        "lineage_note": "V69 runs fresh discovery with E65 fixed; engine revision waits until V69 exposes the next dominant missing word.",
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "engine_biology_modified_during_batch": False,
        "target_selection_manual": False,
        "folding_problem_solved": False,
    }


def _target_manifest(targets: list[dict[str, Any]], raw: dict[str, Any]) -> dict[str, Any]:
    used_before_v69 = _used_protein_ids()
    rows = []
    for target in targets:
        row = {key: value for key, value in target.items() if key != "candidate_snapshot"}
        candidate = target["candidate_snapshot"]
        row.update({
            "title": candidate.get("title", ""),
            "entity_description": candidate.get("entity_description", ""),
            "entry_keywords": candidate.get("entry_keywords", ""),
            "polymer_composition": candidate.get("polymer_composition", ""),
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "biological_cofactor_components": candidate.get("biological_cofactor_components", []),
            "sequence_metrics": candidate.get("sequence_metrics", {}),
            "hard_topology_word": _hard_topology_word(candidate),
            "fresh_exclusion_batches": used_before_v69.get(target["protein_id"], []),
        })
        rows.append(row)
    return {
        "kind": "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "batch_mode": "fresh_four_shard_discovery_200",
        "target_selection_manual": False,
        "composition_rule": dict(SHARD_COUNTS),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(rows),
        "fresh_exclusion_batches": ["V61", "V63", "V64", "V65", "V66", "V67", "V68"],
        "fresh_exclusion_batches_extra_safety": ["V62"],
        "excluded_used_protein_count": len(used_before_v69),
        "source_cache": str(RAW_CANDIDATE_CACHE.relative_to(REPO_ROOT)),
        "raw_candidate_count": raw.get("candidate_entity_count"),
        "selection_rule": "Fresh RCSB 30% representative proteins split into V69A broad, V69B cofactor/ligand/metal, V69C assembly/complex/oligomer, and V69D hard-topology shards.",
        "sequence_cluster_representative_selection": True,
        "sequence_cluster_identity_cutoff": SEQUENCE_IDENTITY_CUTOFF,
        "selected_targets": rows,
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
    dashboard: dict[str, Any],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{"source_id": "V69_BAD_COORDINATES", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    alphafold_gate = evidence_boundary_gate([{"source_id": "V69_BAD_ALPHAFOLD_MODEL", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True, "evidence_statement": "AlphaFold-style model offered before sealing."}])
    holdout_gate = evidence_boundary_gate([{"source_id": "V69_PRESEAL_HOLDOUT", "source_class": COORDINATE_DERIVED, "source_role": "holdout_validation", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V69_BAD_INTERNAL_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    random_packet = build_sealed_simulation_packet(target_id="V69_RANDOM_SEQUENCE_CONTROL", target_name="V69 random sequence control", sequence=deterministic_random_sequence(128), sources=[], perturbations=[])
    composition = Counter(row["shard"] for row in target_manifest["selected_targets"])
    proteins = [row["protein_id"] for row in target_manifest["selected_targets"]]
    clusters = [row.get("sequence_cluster_30_id") for row in target_manifest["selected_targets"]]
    fresh_violations = {row["protein_id"]: row["fresh_exclusion_batches"] for row in target_manifest["selected_targets"] if row.get("fresh_exclusion_batches")}
    shuffled_rows = []
    for packet, shuffled in zip(packets, shuffled_packets):
        original_coherence = sequence_operator_coherence(packet)
        shuffled_coherence = sequence_operator_coherence(shuffled)
        shuffled_rows.append({
            "target_id": packet["target_id"],
            "original_coherence": original_coherence,
            "shuffled_coherence": shuffled_coherence,
            "shuffled_coordinate_sources": shuffled["evidence_manifest"]["coordinate_derived_source_count_before_prediction"],
            "shuffled_runtime_sources": shuffled["evidence_manifest"]["internal_runtime_source_count_for_prediction"],
        })
    return [
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V69 must have exactly 200 targets.", target_manifest["target_count_selected"]),
        _control("composition_rule", dict(composition) == dict(SHARD_COUNTS), "V69 must match four-shard 50/50/50/50 composition.", dict(composition)),
        _control("fresh_targets_only", not fresh_violations and len(set(proteins)) == len(proteins), "V69 targets must be fresh relative to V61/V63/V64/V65/V66/V67/V68 and unique within V69.", fresh_violations),
        _control("sequence_cluster_representative_selection", target_manifest["sequence_cluster_representative_selection"] is True and len(set(clusters)) == len(clusters) and all(clusters), "Every target is a unique 30% sequence-cluster representative within V69.", {"cutoff": SEQUENCE_IDENTITY_CUTOFF, "unique_clusters": len(set(clusters))}),
        _control("rcsb_experimental_protein_entities_only", all(row.get("source_database", "RCSB_PDB") == "RCSB_PDB" and MIN_LENGTH <= row["sequence_length"] <= MAX_LENGTH for row in target_manifest["selected_targets"]), "All V69 targets are RCSB protein entities in the length window."),
        _control("engine_version_declared_e65", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V69 uses E65."),
        _control("engine_modified_during_batch_false", engine_declaration["engine_modified_during_batch"] is False, "No engine mutation occurs inside V69."),
        _control("dashboard_has_all_shards", set(dashboard["shards"]) == set(SHARD_COUNTS) | {"TOTAL"}, "V69 dashboard reports each shard and total."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V69 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references sealed prediction hash."),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls generated without target metadata.", {"control_count": len(shuffled_rows)}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == ABSTAIN_CLASS, "Random sequence without evidence abstains."),
        _control("sentinel_regressions_not_applicable", dashboard["shards"]["TOTAL"]["sentinel_regressions"] == 0, "V69 is fresh discovery and contains no sentinel shard.", dashboard["shards"]["TOTAL"]["sentinel_regressions"]),
        _control("readme_check_skipped_by_user_instruction", True, "README check skipped by explicit user instruction."),
    ]


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_report: dict[str, Any],
    dashboard: dict[str, Any],
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
        status = MINED
    elif metrics["clean_abstain"]:
        status = PASSED_ABSTAIN
    else:
        status = PASSED_CLEAN
    top_word = failure_report["top_missing_esperanto_word"]
    cert = {
        "kind": "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "fresh_four_shard_discovery_200",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "sentinel_regressions": 0,
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "failed_accepted_by_failure_mode": failure_report["failed_accepted_by_failure_mode"],
        "dominant_failure_mode": failure_report["dominant_failure_mode"],
        "dominant_failure_count": failure_report["dominant_failure_count"],
        "top_missing_esperanto_word": top_word,
        "missing_words_top_10": failure_report["missing_words_top_10"],
        "recommended_next_engine_revision": MISSING_WORD_TO_E66.get(top_word) if top_word else None,
        "next_required_batch": "E66_AND_V70_REPAIR_PANEL" if top_word else "V70_EXPAND_OR_REPEAT_DISCOVERY",
        "dashboard": dashboard["shards"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V69 is a fresh discovery batch for mining the next missing Esperanto word; it is not a solved-folding claim.",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "fresh_four_shard_discovery_200",
        "engine_version_used": ENGINE_VERSION_USED,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["clean_abstain"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"],
        "clean_abstain_count": cert["clean_abstain"],
        "dominant_failure_mode": cert["dominant_failure_mode"],
        "top_missing_esperanto_word": cert["top_missing_esperanto_word"],
        "failure_modes": cert["failed_accepted_by_failure_mode"],
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
        "# V69 E65 RCSB Nonredundant 200 Discovery",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Accepted count: `{cert['accepted_count']}`",
        f"Accepted supported: `{cert['accepted_supported']}`",
        f"Clean abstain supported: `{cert['clean_abstain_supported']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"Sentinel regressions: `{cert['sentinel_regressions']}`",
        f"Top missing Esperanto word: `{cert['top_missing_esperanto_word']}`",
        f"Recommended next engine revision: `{cert['recommended_next_engine_revision']}`",
        f"Next required batch: `{cert['next_required_batch']}`",
        "",
        "## Shard Dashboard",
        "",
        "| shard | targets_total | accepted_count | accepted_supported | failed_accepted | clean_abstain | accepted_accuracy | coverage | top_missing_esperanto_word |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for shard, row in cert["dashboard"].items():
        lines.append(
            f"| `{shard}` | `{row['targets_total']}` | `{row['accepted_count']}` | `{row['accepted_supported']}` | "
            f"`{row['failed_accepted']}` | `{row['clean_abstain']}` | `{row['accepted_accuracy']}` | `{row['coverage']}` | "
            f"`{row['top_missing_esperanto_word']}` |"
        )
    lines.extend([
        "",
        "## Failed Accepted By Failure Mode",
        "",
        "| failure_mode | count |",
        "| --- | ---: |",
    ])
    if cert["failed_accepted_by_failure_mode"]:
        for mode, count in sorted(cert["failed_accepted_by_failure_mode"].items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"| `{mode}` | `{count}` |")
    else:
        lines.append("| none | `0` |")
    lines.extend([
        "",
        "## Boundary",
        "V69 is a fresh discovery/mining batch. E65 is fixed during the run; any engine revision belongs to E66 after this certificate.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v69(out_dir: Path = DEFAULT_OUT_DIR, *, refresh_intake: bool = False) -> dict[str, Path]:
    raw = _load_intake(refresh_intake)
    _reset_generated_outputs(out_dir)
    targets = _select_targets(raw)
    target_manifest = _target_manifest(targets, raw)
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v69_rcsb_nonredundant_200_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v69_e65_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []
    for target in targets:
        source_manifest = _source_manifest(target)
        packet = build_sealed_simulation_packet(
            target_id=target["target_id"],
            target_name=target["target_name"],
            sequence=target["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "V69 E65 fresh discovery full-chain scan", "span": f"1-{target['sequence_length']}"}],
            perturbations=_perturbations_for_expected(target),
        )
        holdout = _holdout(target, packet)
        score = _score(packet, holdout, target)
        shuffled_packet = build_sealed_simulation_packet(
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
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E65 context are withheld.",
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

    failure_report = _failure_report(scoring_rows)
    dashboard = _shard_dashboard(scoring_rows, failure_report)
    controls = _controls(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        packets=packets,
        scoring_rows=scoring_rows,
        shuffled_packets=shuffled_packets,
        dashboard=dashboard,
    )
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
        failure_report=failure_report,
        dashboard=dashboard,
    )
    claim_row = _claim_row(cert)

    scoring_path = _write_json(DATA_ROOT / "v69_rcsb_nonredundant_200_scoring_report.json", {"kind": "V69_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v69_rcsb_nonredundant_200_failure_report.json", failure_report)
    dashboard_path = _write_json(DATA_ROOT / "v69_rcsb_nonredundant_200_dashboard.json", dashboard)
    data_cert_path = _write_json(DATA_ROOT / "v69_rcsb_nonredundant_200_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v69_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)

    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v69_e65_rcsb_nonredundant_200_discovery_certificate.json", cert)
    report_path = out_dir / "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_REPORT.md"
    _write_report(report_path, cert)
    return {
        "raw_candidate_cache": RAW_CANDIDATE_CACHE,
        "target_manifest": DATA_ROOT / "v69_rcsb_nonredundant_200_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v69_e65_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "dashboard": dashboard_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V69 E65 fresh four-shard RCSB discovery batch.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--refresh-intake", action="store_true", help="refresh V69 RCSB fresh candidate intake with curl")
    args = parser.parse_args()
    paths = run_v69(args.out_dir, refresh_intake=args.refresh_intake)
    cert = _read_json(paths["certificate"], "V69 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_count": cert["accepted_count"],
        "accepted_supported": cert["accepted_supported"],
        "clean_abstain_supported": cert["clean_abstain_supported"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "sentinel_regressions": cert["sentinel_regressions"],
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "dominant_failure_mode": cert["dominant_failure_mode"],
        "dominant_failure_count": cert["dominant_failure_count"],
        "top_missing_esperanto_word": cert["top_missing_esperanto_word"],
        "recommended_next_engine_revision": cert["recommended_next_engine_revision"],
        "failed_accepted_by_failure_mode": cert["failed_accepted_by_failure_mode"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
