#!/usr/bin/env python3
from __future__ import annotations

"""Run V75: self-deciding acceptance cortex."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

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
from pharmacotopology.protein_esperanto_physical_calibration import (  # noqa: E402
    write_real_physical_calibration_inputs,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402
import run_v69_e65_rcsb_nonredundant_200_discovery_v0 as v69  # noqa: E402
import run_v71_e66_rcsb_nonredundant_200_discovery_v0 as v71  # noqa: E402
import run_v74_e68_rcsb_nonredundant_200_discovery_v0 as v74  # noqa: E402


BATCH_ID = "V75_SELF_DECIDING_ACCEPTANCE_CORTEX"
CAMPAIGN_ID = "V61_TO_V75_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E69"
BASELINE_ENGINE_VERSION = "E68"
TARGET_COUNT = 200

ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"
MULTIDOMAIN_CLASS = "multidomain_allosteric_architecture"
GLOBULAR_CLASS = "globular_closure"
ASSEMBLY_REQUIRED_CLASS = "assembly_required_folding"
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
METAL_LIGAND_CLASS = "metal_cluster_and_ligand_locked_basin"
DISORDER_BOUNDARY_CLASS = "disorder_boundary_and_fold_upon_binding"

GROUP_COUNTS = OrderedDict([
    ("V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY", 33),
    ("MULTIDOMAIN_HINGE_DOMAIN_SWAP_ALLOSTERY_POSITIVE", 40),
    ("MODULAR_ARCHITECTURE_INTERDOMAIN_LOCK_POSITIVE", 25),
    ("MONOMERIC_GLOBULAR_SENTINEL", 20),
    ("ASSEMBLY_REQUIRED_SENTINEL", 20),
    ("MEMBRANE_TM_SENTINEL", 20),
    ("METAL_LIGAND_SENTINEL", 15),
    ("DISORDER_BOUNDARY_SENTINEL", 15),
    ("MISSING_WORD_CANDIDATE_CLEAN_ABSTAIN_CONTROL", 12),
])

MISSING_CANDIDATE_CONTROL_WORDS = (
    ["disulfide_secretory_redox_context"] * 5
    + ["coiled_coil_register"] * 4
    + ["repeat_solenoid_topology"] * 3
)

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V75"
E69_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E69"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
REAL_COORDINATE_BENCHMARK = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
PHYSICAL_CALIBRATION_ROOT = DATA_ROOT / "physical_calibration"
PHYSICAL_CALIBRATION_INPUTS = PHYSICAL_CALIBRATION_ROOT / "v75_real_physical_calibration_inputs.json"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V74_FAILURES = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "v74_rcsb_nonredundant_200_failure_report.json"
V74_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "v74_rcsb_nonredundant_200_target_manifest.json"
V74_RAW = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "intake" / "raw_rcsb_30pct_representative_entities_v74.json"
V74_CERT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "v74_rcsb_nonredundant_200_certificate.json"
E69_CERT = E69_ROOT / "e69_multidomain_allosteric_architecture_grammar_certificate.json"

BATCH_PASSED = "BATCH_PASSED_SELF_DECIDING_ZERO_FAILED_ACCEPTED"
BATCH_DISCOVERED = "BATCH_DISCOVERED_MISSING_WORD_WITH_CLEAN_SELF_ABSTENTION"
BATCH_CORTEX_BROKEN = "BATCH_FAILED_SELF_DECISION_ZERO_FAILED_ACCEPTED"
BATCH_REPAIR_REQUIRED = "BATCH_REPAIR_REQUIRED"
BATCH_BLOCKED_FOR_LEAKAGE = "BATCH_BLOCKED_FOR_LEAKAGE"
BATCH_CONTROLS_FAILED = "BATCH_CONTROLS_FAILED"


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


def _candidate_id(candidate: dict[str, Any]) -> str:
    return str(candidate.get("protein_id") or candidate.get("target_id") or stable_hash(candidate)[:12])


def _candidate_rank(candidate: dict[str, Any]) -> tuple[int, str]:
    return int(candidate.get("sequence_cluster_representative_rank") or 10**9), _candidate_id(candidate)


def _reset_generated_outputs(out_dir: Path) -> None:
    for relative in ["source_manifests", "sealed_packet_summaries", "holdouts_postseal", "validation", "shuffled_controls", "physical_calibration"]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v75_self_deciding_acceptance_cortex_target_manifest.json",
        "v75_e69_engine_declaration.json",
        "v75_self_deciding_acceptance_cortex_scoring_report.json",
        "v75_self_deciding_acceptance_cortex_failure_report.json",
        "v75_self_deciding_acceptance_cortex_dashboard.json",
        "v75_self_deciding_acceptance_cortex_certificate.json",
        "v75_e69_multidomain_repair_target_manifest.json",
        "v75_e69_multidomain_repair_scoring_report.json",
        "v75_e69_multidomain_repair_failure_report.json",
        "v75_e69_multidomain_repair_dashboard.json",
        "v75_e69_multidomain_repair_certificate.json",
        "v75_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    legacy_out_dir = RUN_ROOT / "V75_E69_MULTIDOMAIN_REPAIR_AND_ZERO_FAILED_ACCEPT_FIREWALL"
    if legacy_out_dir.exists():
        shutil.rmtree(legacy_out_dir)


def _candidate_from_manifest_row(row: dict[str, Any]) -> dict[str, Any]:
    candidate = dict(row)
    candidate.setdefault("source_database", "RCSB_PDB")
    candidate.setdefault("source_urls", {"entry": row.get("entry_url", ""), "polymer_entity": row.get("polymer_entity_url", "")})
    candidate.setdefault("sequence_metrics", row.get("sequence_metrics") or v61._sequence_metrics(row["sequence"]))
    candidate.setdefault("annotations", [])
    candidate.setdefault("feature_types", [])
    candidate.setdefault("nonpolymer_bound_components", row.get("biological_cofactor_components", []))
    candidate.setdefault("organisms", [])
    candidate.setdefault("experimental_method", "")
    candidate.setdefault("release_date", "")
    return candidate


def _raw_candidates() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in _read_json(V74_RAW, "V74 raw candidate cache").get("candidates", []):
        if not isinstance(row, dict):
            continue
        candidate = dict(row)
        protein_id = _candidate_id(candidate)
        sequence = str(candidate.get("sequence") or "")
        if not protein_id or not sequence or protein_id in seen:
            continue
        candidate.setdefault("sequence_metrics", v61._sequence_metrics(sequence))
        seen.add(protein_id)
        rows.append(candidate)
    return sorted(rows, key=_candidate_rank)


def _candidate_text(candidate: dict[str, Any]) -> str:
    return v74._candidate_text(candidate)


def _has_any(candidate: dict[str, Any], tokens: list[str]) -> bool:
    text = _candidate_text(candidate)
    return any(token in text for token in tokens)


def _MISSING_CANDIDATE_word(candidate: dict[str, Any]) -> str | None:
    word = v74._disulfide_secretory_word(candidate)
    if word in {"disulfide_secretory_redox_context", "signal_peptide_vs_true_TM"}:
        return word
    word = v74._coiled_repeat_word(candidate)
    if word in {"coiled_coil_register", "repeat_solenoid_topology"}:
        return word
    return None


def _is_globular_sentinel_candidate(candidate: dict[str, Any]) -> bool:
    return (
        not v69._true_tm(candidate)
        and v69._assembly_required_word(candidate) is None
        and v69._metal_or_ligand_locked_word(candidate) is None
        and not v71._disorder_enriched(candidate)
        and v74._beta_word(candidate) is None
        and v74._multidomain_word(candidate) is None
        and _MISSING_CANDIDATE_word(candidate) is None
    )


def _pick(candidates: list[dict[str, Any]], *, count: int, used: set[str], predicate: Callable[[dict[str, Any]], bool]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for candidate in sorted(candidates, key=_candidate_rank):
        protein_id = _candidate_id(candidate)
        if protein_id in used or not predicate(candidate):
            continue
        selected.append(candidate)
        used.add(protein_id)
        if len(selected) == count:
            break
    if len(selected) != count:
        raise SystemExit(f"selected {len(selected)} candidates; expected {count}")
    return selected


def _target_from_candidate(
    *,
    ordinal: int,
    group: str,
    candidate: dict[str, Any],
    expected: str,
    required_word: str | None,
    source_family: str,
    expected_decision: str = "accepted",
    lineage_source_target: str | None = None,
    sentinel_family: str | None = None,
) -> dict[str, Any]:
    entry_url = candidate.get("entry_url") or candidate.get("source_urls", {}).get("entry", "")
    polymer_entity_url = candidate.get("polymer_entity_url") or candidate.get("source_urls", {}).get("polymer_entity", "")
    return {
        "target_id": f"V75_{ordinal:03d}_{_safe_id(group)}_{_safe_id(_candidate_id(candidate))}",
        "panel_group": group,
        "source_family": source_family,
        "lineage_source_target": lineage_source_target,
        "sentinel_family": sentinel_family,
        "expected_decision": expected_decision,
        "protein_id": _candidate_id(candidate),
        "entry_id": str(candidate.get("entry_id", "")),
        "entity_id": str(candidate.get("entity_id", "")),
        "sequence": candidate["sequence"],
        "sequence_length": int(candidate.get("sequence_length") or len(candidate["sequence"])),
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "expected_mechanism_class": expected,
        "required_esperanto_word": required_word,
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip(),
        "entry_url": entry_url,
        "polymer_entity_url": polymer_entity_url,
        "postseal_truth_basis": [
            f"V75 E69 zero-failed-accepted self-decision cortex group {group}.",
            f"Required word: {required_word or 'none'}.",
            f"Expected decision: {expected_decision}.",
            "Coordinates, contacts, native topology, and post-seal validation labels are blocked before sealing.",
        ],
        "candidate_snapshot": dict(candidate),
    }


def _select_targets() -> list[dict[str, Any]]:
    failures = _read_json(V74_FAILURES, "V74 failure report")["failure_grammar_rows"]
    manifest_rows = _read_json(V74_MANIFEST, "V74 target manifest")["selected_targets"]
    by_target = {row["target_id"]: _candidate_from_manifest_row(row) for row in manifest_rows}
    raw = _raw_candidates()
    targets: list[dict[str, Any]] = []
    used: set[str] = set()
    ordinal = 1

    replay_rows = [row for row in failures if row["failure_mode"] == "multidomain_allostery"]
    if len(replay_rows) != GROUP_COUNTS["V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY"]:
        raise SystemExit(f"V74 multidomain failures are {len(replay_rows)}; expected 33")
    for row in replay_rows:
        candidate = by_target[row["target_id"]]
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group="V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY",
            candidate=candidate,
            expected=MULTIDOMAIN_CLASS,
            required_word="multidomain_allostery",
            source_family="V74",
            lineage_source_target=row["target_id"],
        ))
        used.add(_candidate_id(candidate))
        ordinal += 1

    groups = [
        (
            "MULTIDOMAIN_HINGE_DOMAIN_SWAP_ALLOSTERY_POSITIVE",
            40,
            lambda c: v74._multidomain_word(c) == "multidomain_allostery",
            MULTIDOMAIN_CLASS,
            "multidomain_allostery",
            "MULTIDOMAIN_ALLOSTERY",
        ),
        (
            "MODULAR_ARCHITECTURE_INTERDOMAIN_LOCK_POSITIVE",
            25,
            lambda c: v74._multidomain_word(c) is not None,
            MULTIDOMAIN_CLASS,
            "modular_architecture",
            "MODULAR_ARCHITECTURE",
        ),
        (
            "MONOMERIC_GLOBULAR_SENTINEL",
            20,
            _is_globular_sentinel_candidate,
            GLOBULAR_CLASS,
            None,
            "MONOMERIC_GLOBULAR",
        ),
        (
            "ASSEMBLY_REQUIRED_SENTINEL",
            20,
            lambda c: v69._assembly_required_word(c) is not None,
            ASSEMBLY_REQUIRED_CLASS,
            "assembly_required_core",
            "ASSEMBLY_REQUIRED",
        ),
        (
            "MEMBRANE_TM_SENTINEL",
            20,
            v69._true_tm,
            MEMBRANE_CLASS,
            None,
            "TRUE_TM",
        ),
        (
            "METAL_LIGAND_SENTINEL",
            15,
            lambda c: v69._metal_or_ligand_locked_word(c) is not None,
            METAL_LIGAND_CLASS,
            "ligand_locked_basin",
            "METAL_LIGAND",
        ),
        (
            "DISORDER_BOUNDARY_SENTINEL",
            15,
            v71._disorder_enriched,
            DISORDER_BOUNDARY_CLASS,
            "IDR_boundary",
            "DISORDER_BOUNDARY",
        ),
    ]
    for group, count, predicate, expected, word, sentinel_family in groups:
        for candidate in _pick(raw, count=count, used=used, predicate=predicate):
            targets.append(_target_from_candidate(
                ordinal=ordinal,
                group=group,
                candidate=candidate,
                expected=expected,
                required_word=word,
                source_family="V74_RAW",
                sentinel_family=sentinel_family,
            ))
            ordinal += 1

    for word in MISSING_CANDIDATE_CONTROL_WORDS:
        candidate = _pick(raw, count=1, used=used, predicate=lambda c, w=word: _MISSING_CANDIDATE_word(c) == w)[0]
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group="MISSING_WORD_CANDIDATE_CLEAN_ABSTAIN_CONTROL",
            candidate=candidate,
            expected=ABSTAIN_CLASS,
            required_word=word,
            source_family="V74_RAW",
            expected_decision="abstain_recommended",
            sentinel_family=f"MISSING_{word}",
        ))
        ordinal += 1

    composition = Counter(row["panel_group"] for row in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != dict(GROUP_COUNTS):
        raise SystemExit(f"bad V75 composition: total={len(targets)} composition={dict(composition)}")
    return targets


def _context_for_target(target: dict[str, Any]) -> tuple[list[str], str]:
    group = target["panel_group"]
    word = target.get("required_esperanto_word")
    if group in {"V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY", "MULTIDOMAIN_HINGE_DOMAIN_SWAP_ALLOSTERY_POSITIVE"}:
        marks = [
            "multidomain_allostery",
            "domain_boundary",
            "hinge_region",
            "interdomain_lock",
            "allosteric_basin_shift",
            "domain_reorientation",
            "modular_architecture",
        ]
    elif group == "MODULAR_ARCHITECTURE_INTERDOMAIN_LOCK_POSITIVE":
        marks = [
            "modular_architecture",
            "interdomain_lock",
            "domain_boundary",
            "hinge_region",
            "domain_reorientation",
            "multidomain_allostery",
        ]
    elif group == "MONOMERIC_GLOBULAR_SENTINEL":
        marks = ["soluble_monomeric_core_context", "complete soluble monomer", "standalone soluble fold"]
    elif group == "ASSEMBLY_REQUIRED_SENTINEL":
        marks = ["assembly_required_core", "assembly_required_folding", "partner_completed_core", "interface_buried_hydrophobicity", "monomer_incomplete_topology"]
    elif group == "MEMBRANE_TM_SENTINEL":
        marks = ["membrane_context_strong", "transmembrane_context", "topology_evidence"]
    elif group == "METAL_LIGAND_SENTINEL":
        marks = ["metal_cluster_geometry", "ligand_locked_basin", "coordination_shell_integrity", "apo_holo_basin_shift"]
    elif group == "DISORDER_BOUNDARY_SENTINEL":
        marks = ["disorder_context", "IDR_boundary", "structured_domain_plus_IDR_tail", "fold_upon_binding_region", "phase_prone_low_complexity"]
    elif word == "disulfide_secretory_redox_context":
        marks = ["disulfide_secretory_redox_context", "secretory_redox_context", "disulfide_bond_topology", "cysteine_pairing_constraint", "extracellular_stabilized_fold"]
    elif word == "coiled_coil_register":
        marks = ["coiled_coil_register", "heptad_repeat", "register_alignment", "parallel_vs_antiparallel_coil", "oligomeric_coiled_coil_core"]
    elif word == "repeat_solenoid_topology":
        marks = ["repeat_solenoid_topology", "repeat_unit", "solenoid_axis", "curved_repeat_stack", "local_repeat_closure", "global_repeat_topology"]
    else:
        marks = []
    statement = " ".join(marks) if marks else "no explicit context mark"
    return marks, statement


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate_snapshot"]
    marks, statement = _context_for_target(target)
    suffix = stable_hash({"v75_target_id": target["target_id"]})[:12]
    return {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_SOURCE_MANIFEST_v0",
        "target_id": target["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "protein_id": target["protein_id"],
        "sequence": target["sequence"],
        "sequence_length": target["sequence_length"],
        "prediction_sources": [
            {
                "source_id": f"V75_RAW_SEQUENCE_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": target["polymer_entity_url"],
            },
            {
                "source_id": f"V75_PUBLIC_METADATA_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": (
                    "Public non-coordinate RCSB metadata is available through the source URL. "
                    "V75 uses the explicit self-decision cortex context source for allowed biological words. "
                    "Coordinates, contacts, ligand geometry, native topology, and post-seal validation labels are unopened before prediction hash."
                ),
                "source_url": target["entry_url"],
            },
            {
                "source_id": f"V75_E69_CONTEXT_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": marks,
                "evidence_statement": (
                    "V75 allowed non-coordinate Esperanto context marks: "
                    + statement
                    + ". Unknown self-decision words must cleanly abstain until their grammar exists."
                ),
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
        "folding_problem_solved": False,
    }


def _expected_observables(expected: str, required_word: str | None) -> list[dict[str, Any]]:
    if expected == MULTIDOMAIN_CLASS:
        checks = [
            {"check_id": "multidomain_allostery_supported", "metric": "multidomain_allostery", "comparator": ">=", "threshold": 0.30},
            {"check_id": "interdomain_lock_supported", "metric": "interdomain_lock", "comparator": ">=", "threshold": 0.24},
        ]
        if required_word in {"modular_architecture", "domain_swapping", "allosteric_basin_shift", "domain_reorientation"}:
            checks.append({"check_id": f"{required_word}_supported", "metric": required_word, "comparator": ">=", "threshold": 0.20})
        return checks
    if expected == GLOBULAR_CLASS:
        return [{"check_id": "globular_contact_supported", "metric": "contact_probability", "comparator": ">=", "threshold": 0.35}]
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return [{"check_id": "partner_completed_core_supported", "metric": "partner_completed_core", "comparator": ">=", "threshold": 0.48}]
    if expected == MEMBRANE_CLASS:
        return [{"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.50}]
    if expected == METAL_LIGAND_CLASS:
        return [{"check_id": "ligand_locked_basin_supported", "metric": "ligand_locked_basin", "comparator": ">=", "threshold": 0.42}]
    if expected == DISORDER_BOUNDARY_CLASS:
        return [
            {"check_id": "idr_boundary_supported", "metric": "IDR_boundary", "comparator": ">=", "threshold": 0.30},
            {"check_id": "compact_single_fold_rejected", "metric": "basin:compact_single_fold", "comparator": "<=", "threshold": 0.12},
        ]
    return []


def _holdout(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "expected_decision": target["expected_decision"],
        "expected_mechanism_class": target["expected_mechanism_class"],
        "required_esperanto_word": target["required_esperanto_word"],
        "expected_observables": _expected_observables(target["expected_mechanism_class"], target["required_esperanto_word"]),
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_truth_basis": target["postseal_truth_basis"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V75_POSTSEAL_HOLDOUT",
                "source_class": COORDINATE_DERIVED,
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "entry_url": target["entry_url"],
                "polymer_entity_url": target["polymer_entity_url"],
            }
        ],
    }


def _perturbations_for_target(target: dict[str, Any]) -> list[dict[str, Any]]:
    target_id = target["target_id"]
    expected = target["expected_mechanism_class"]
    if expected == MULTIDOMAIN_CLASS:
        return [{"perturbation_id": f"{target_id}_INTERDOMAIN_LOCK_DAMAGE", "description": "damage interdomain lock or hinge coupling", "operator_scales": {"interface_operator": 0.45}, "lock_damage": 0.45, "metric": "interdomain_lock", "expected_direction": "decrease"}]
    if expected == GLOBULAR_CLASS:
        return [{"perturbation_id": f"{target_id}_CORE_DAMAGE", "description": "damage closure core", "operator_scales": {"closure_operator": 0.45}, "metric": "contact_probability", "expected_direction": "decrease"}]
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return [{"perturbation_id": f"{target_id}_ASSEMBLY_INTERFACE_DAMAGE", "description": "damage partner-completed interface", "operator_scales": {"interface_operator": 0.42, "closure_operator": 0.70}, "interface_disruption": 0.45, "metric": "partner_completed_core", "expected_direction": "decrease"}]
    if expected == MEMBRANE_CLASS:
        return [{"perturbation_id": f"{target_id}_MEMBRANE_DAMAGE", "description": "damage topology/proteostasis route", "operator_scales": {"membrane_pressure_operator": 0.55, "proteostasis_operator": 0.55}, "damage": 0.40, "metric": "proteostasis_routing", "expected_direction": "decrease"}]
    if expected == METAL_LIGAND_CLASS:
        return [{"perturbation_id": f"{target_id}_METAL_OR_LIGAND_REMOVAL", "description": "remove metal/ligand pressure", "operator_scales": {"interface_operator": 0.45}, "cofactor_loss": 0.45, "metric": "ligand_locked_basin", "expected_direction": "decrease"}]
    if expected == DISORDER_BOUNDARY_CLASS:
        return [{"perturbation_id": f"{target_id}_DISORDER_PARTNER_LOSS", "description": "remove disorder-boundary partner/motif pressure", "operator_scales": {"interface_operator": 0.45}, "partner_loss": 0.45, "motif_damage": 0.25, "metric": "fold_upon_binding_region", "expected_direction": "decrease"}]
    return []


def _required_word_supported(required_word: str | None, packet: dict[str, Any]) -> bool:
    if not required_word:
        return True
    final = packet["trajectory_summary"]["final_state_summary"]
    mechanism = packet["selected_mechanism_grammar"]["mechanism_class"]
    if required_word in {"disulfide_secretory_redox_context", "coiled_coil_register", "repeat_solenoid_topology"}:
        return False
    if required_word == "assembly_required_core":
        return mechanism == ASSEMBLY_REQUIRED_CLASS and final.get("partner_completed_core", 0.0) > 0.0
    if required_word == "ligand_locked_basin":
        return mechanism == METAL_LIGAND_CLASS and final.get("ligand_locked_basin", 0.0) > 0.0
    if required_word == "IDR_boundary":
        return mechanism == DISORDER_BOUNDARY_CLASS and final.get("IDR_boundary", 0.0) > 0.0
    if required_word in {"multidomain_allostery", "modular_architecture", "domain_swapping", "allosteric_basin_shift", "domain_reorientation"}:
        return mechanism == MULTIDOMAIN_CLASS and final.get(required_word, 0.0) > 0.0
    return True


def _score(packet: dict[str, Any], holdout: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    judge = packet["self_decision_judge"]
    decision = judge["acceptance_decision"]
    final_self_decision = judge["final_self_decision"]
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    required_word = holdout.get("required_esperanto_word")
    accepted = decision == "accepted"
    blocked = decision == "blocked_for_leakage"
    coarse_supported = predicted == expected and validation["score_label"] == "supported"
    word_supported = _required_word_supported(required_word, packet)
    accepted_supported = accepted and coarse_supported and word_supported
    clean_abstain_supported = (
        final_self_decision == "clean_abstain_missing_word"
        and target["expected_decision"] == "abstain_recommended"
        and judge.get("missing_word_candidate") == required_word
    )
    blocked_supported = blocked and target["expected_decision"].startswith("blocked")
    if accepted_supported or clean_abstain_supported or blocked_supported:
        score_label = "supported"
    elif decision == "abstain_recommended":
        score_label = "abstained"
    elif blocked:
        score_label = "blocked"
    else:
        score_label = "contradicted"
    return {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "panel_group": target["panel_group"],
        "sentinel_family": target.get("sentinel_family"),
        "lineage_source_target": target.get("lineage_source_target"),
        "acceptance_decision": decision,
        "final_self_decision": final_self_decision,
        "self_decision_reason": judge["self_decision_reason"],
        "internal_consensus": judge["internal_consensus"],
        "dominance_law": judge["dominance_law"],
        "cross_view_binding": judge["cross_view_binding"],
        "cross_view_missing_view_count": len(judge["cross_view_binding_probe"]["missing_view_families"]),
        "masking_stability": judge["masking_stability"],
        "wrong_grammar_separation": judge["wrong_grammar_separation"],
        "counterfactual_separation": judge["counterfactual_separation"],
        "operator_basis_stability": judge["operator_basis_stability"],
        "coefficient_probe_mode": judge["coefficient_probe_mode"],
        "temporal_binding": judge["temporal_binding"],
        "physics_grounding_status": judge["physics_grounding_status"],
        "physical_grounding_status": judge["physical_grounding_status"],
        "physical_backend_available": judge["physical_backend_available"],
        "physical_basis_claim_allowed": judge["physical_basis_claim_allowed"],
        "real_physical_calibration_inputs_used": judge["real_physical_calibration_inputs_used"],
        "real_physical_calibration_kind": judge["real_physical_calibration_kind"],
        "real_physical_calibration_row_count": judge["real_physical_calibration_row_count"],
        "real_physical_calibration_hash": judge["real_physical_calibration_hash"],
        "contradiction_count": judge["contradiction_count"],
        "missing_word_candidate": judge.get("missing_word_candidate"),
        "unknown_word_signals": judge.get("unknown_word_signals", []),
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "selected_multidomain_word": packet["selected_mechanism_grammar"].get("selected_multidomain_word"),
        "required_esperanto_word": required_word,
        "required_esperanto_word_supported": word_supported,
        "accepted_supported": accepted_supported,
        "clean_abstain_supported": clean_abstain_supported,
        "blocked_supported": blocked_supported,
        "score_label": score_label,
        "validation_checks": validation["checks"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout["holdout_opened_after_prediction_hash"],
        "coordinate_truth_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V75_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "selected_multidomain_word": mechanism.get("selected_multidomain_word"),
        "self_decision_judge": packet["self_decision_judge"],
        "acceptance_view": packet["acceptance_firewall"],
        "operator_names": packet["operator_field"]["operator_names"],
        "active_operator_count": packet["operator_field"]["active_operator_count"],
        "trajectory_final_state_summary": packet["trajectory_summary"]["final_state_summary"],
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    blocked = [row for row in rows if row["acceptance_decision"].startswith("blocked")]
    accepted_supported = [row for row in accepted if row["accepted_supported"]]
    clean_abstain_supported = [row for row in abstained if row["clean_abstain_supported"]]
    blocked_supported = [row for row in blocked if row["blocked_supported"]]
    failed_accepted = [row for row in accepted if not row["accepted_supported"]]
    supported = [row for row in rows if row["score_label"] == "supported"]
    unknown_word_abstentions = [row for row in clean_abstain_supported if row.get("missing_word_candidate")]
    return {
        "targets_total": len(rows),
        "accepted_count": len(accepted),
        "accepted_supported": len(accepted_supported),
        "clean_abstain": len(abstained),
        "clean_abstain_supported": len(clean_abstain_supported),
        "blocked_count": len(blocked),
        "blocked_supported": len(blocked_supported),
        "supported_count": len(supported),
        "failed_accepted": len(failed_accepted),
        "failed_accepted_count": len(failed_accepted),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else 1.0,
        "raw_accuracy": len(supported) / len(rows) if rows else None,
        "coverage": len(accepted) / len(rows) if rows else None,
        "unknown_word_abstentions": len(unknown_word_abstentions),
    }


def _failure_mode(row: dict[str, Any]) -> str:
    if row["accepted_supported"]:
        return "none"
    if row["acceptance_decision"] == "accepted":
        return str(row.get("required_esperanto_word") or row["expected_mechanism_class"])
    if row["acceptance_decision"] == "abstain_recommended":
        return str(row.get("missing_word_candidate") or "clean_abstain")
    return "blocked"


def _failure_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    failed = [row for row in scoring_rows if row["acceptance_decision"] == "accepted" and not row["accepted_supported"]]
    rows = [
        {
            "target_id": row["target_id"],
            "panel_group": row["panel_group"],
            "failure_mode": _failure_mode(row),
            "engine_thought": row["predicted_mechanism_class"],
            "reality_showed": row["expected_mechanism_class"],
            "required_esperanto_word": row.get("required_esperanto_word"),
            "missing_esperanto_word": row.get("required_esperanto_word"),
            "autopsy_sentence": (
                f"The engine thought: {row['predicted_mechanism_class']}. "
                f"Reality showed: {row['expected_mechanism_class']}. "
                f"Missing Esperanto word: {row.get('required_esperanto_word')}."
            ),
        }
        for row in failed
    ]
    return {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_count": len(rows),
        "failed_accepted_by_failure_mode": dict(Counter(row["failure_mode"] for row in rows)),
        "failure_grammar_rows": rows,
    }


def _dashboard(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    dashboard: dict[str, Any] = {}
    for group in GROUP_COUNTS:
        rows = [row for row in scoring_rows if row["panel_group"] == group]
        missing = Counter(row["missing_word_candidate"] for row in rows if row.get("missing_word_candidate"))
        dashboard[group] = {
            **_metrics(rows),
            "top_missing_esperanto_word": missing.most_common(1)[0][0] if missing else None,
            "missing_esperanto_word_counts": dict(missing),
        }
    missing_total = Counter(row["missing_word_candidate"] for row in scoring_rows if row.get("missing_word_candidate"))
    dashboard["TOTAL"] = {
        **_metrics(scoring_rows),
        "top_missing_esperanto_word": missing_total.most_common(1)[0][0] if missing_total else None,
        "missing_esperanto_word_counts": dict(missing_total),
    }
    return {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_DASHBOARD_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "zero_failed_accepted_required": True,
        "shards": dashboard,
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V75_E69_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_lineage": ["E60", "E61", "E62", "E63", "E64", "E65", "E66", "E67", "E68", "E69"],
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "operator_names": UNIVERSAL_OPERATORS,
        "mechanism_classes": MECHANISM_CLASSES,
        "zero_failed_accepted_required": True,
        "self_decision_cortex_enabled": True,
        "engine_modified_during_batch": False,
        "folding_problem_solved": False,
    }


def _target_manifest(targets: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for target in targets:
        row = {key: value for key, value in target.items() if key != "candidate_snapshot"}
        candidate = target["candidate_snapshot"]
        row.update({
            "title": candidate.get("title", ""),
            "entity_description": candidate.get("entity_description", ""),
            "entry_keywords": candidate.get("entry_keywords", ""),
            "polymer_composition": candidate.get("polymer_composition", ""),
            "sequence_metrics": candidate.get("sequence_metrics", {}),
        })
        rows.append(row)
    return {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(rows),
        "composition_rule": dict(GROUP_COUNTS),
        "zero_failed_accepted_required": True,
        "selected_targets": rows,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    target_manifest: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
    physical_calibration_inputs: dict[str, Any],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{"source_id": "V75_BAD_COORDINATES", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V75_BAD_INTERNAL_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    random_packet = build_sealed_simulation_packet(
        target_id="V75_RANDOM_SEQUENCE_CONTROL",
        target_name="V75 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
        physical_calibration_inputs=physical_calibration_inputs,
    )
    unknown_packet = build_sealed_simulation_packet(
        target_id="V75_CORTEX_UNKNOWN_WORD_CONTROL",
        target_name="V75 self-decision unknown word control",
        sequence=("CGPC" * 50),
        sources=[
            {
                "source_id": "V75_UNKNOWN_DISULFIDE_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "disulfide_secretory_redox_context disulfide_bond_topology secretory_redox_context",
            }
        ],
        perturbations=[],
        physical_calibration_inputs=physical_calibration_inputs,
    )
    composition = Counter(row["panel_group"] for row in target_manifest["selected_targets"])
    controls = [
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V75 must have exactly 200 targets.", target_manifest["target_count_selected"]),
        _control("composition_rule", dict(composition) == dict(GROUP_COUNTS), "V75 must match requested self-decision cortex composition.", dict(composition)),
        _control("zero_failed_accepted_required", all(row["acceptance_decision"] != "accepted" or row["accepted_supported"] for row in scoring_rows), "No accepted target may be unsupported."),
        _control("accepted_cross_view_binding", all(row["acceptance_decision"] != "accepted" or row["cross_view_missing_view_count"] == 0 for row in scoring_rows), "Accepted targets must satisfy their self-required cross-view binding."),
        _control("accepted_operator_basis_stability", all(row["acceptance_decision"] != "accepted" or row["operator_basis_stability"] != "coefficient_assignment_sensitive" for row in scoring_rows), "Accepted targets must be stable under endogenous operator-basis probes."),
        _control("accepted_temporal_binding", all(row["acceptance_decision"] != "accepted" or row["temporal_binding"] != "selected_observable_temporal_conflict" for row in scoring_rows), "Accepted targets must avoid internal temporal trajectory conflicts."),
        _control("physical_basis_claim_blocked", all(row["physical_basis_claim_allowed"] is False for row in scoring_rows), "The coarse operator cortex must not claim physical coefficient grounding."),
        _control("physical_grounding_status_reported", all(row.get("physical_grounding_status") for row in scoring_rows), "Every target reports physical grounding status."),
        _control("real_physical_calibration_inputs_present", all(row["real_physical_calibration_inputs_used"] for row in scoring_rows), "Every target receives the real physical calibration manifest.", Counter(row["real_physical_calibration_inputs_used"] for row in scoring_rows)),
        _control("real_physical_calibration_rows_real_coordinate", physical_calibration_inputs.get("source_coordinate_database") == "RCSB_PDB" and int(physical_calibration_inputs.get("row_count") or 0) > 0, "Calibration inputs are derived from locked real RCSB coordinate rows.", {"source_coordinate_database": physical_calibration_inputs.get("source_coordinate_database"), "row_count": physical_calibration_inputs.get("row_count")}),
        _control("real_physical_calibration_truth_boundary", physical_calibration_inputs.get("target_native_excluded_from_calibration") is True and physical_calibration_inputs.get("target_native_contacts_used_before_prediction") is False and physical_calibration_inputs.get("coordinate_truth_used_as_prediction_input") is False, "Real calibration inputs preserve the prediction truth boundary.", {key: physical_calibration_inputs.get(key) for key in ["target_native_excluded_from_calibration", "target_native_contacts_used_before_prediction", "coordinate_truth_used_as_prediction_input", "leave_one_target_out_calibration"]}),
        _control("self_decision_unknown_word_control", unknown_packet["self_decision_judge"]["final_self_decision"] == "clean_abstain_missing_word" and unknown_packet["self_decision_judge"]["missing_word_candidate"] == "disulfide_secretory_redox_context", "Unknown disulfide word must cleanly abstain.", unknown_packet["self_decision_judge"]),
        _control("random_sequence_control", random_packet["self_decision_judge"]["final_self_decision"] == "clean_abstain_low_internal_consensus", "Random sequence without evidence abstains.", random_packet["self_decision_judge"]),
        _control("shuffled_sequence_controls_reported", len(shuffled_packets) == TARGET_COUNT, "Composition-preserving shuffled controls are generated.", len(shuffled_packets)),
        _control("coordinate_leakage_control", coord_gate["allowed_initialization_source_ids"] == [] and coord_gate["coordinate_derived_source_count_before_prediction"] == 1, "Coordinate-derived source blocks prediction.", coord_gate),
        _control("internal_runtime_leakage_control", runtime_gate["allowed_initialization_source_ids"] == [] and runtime_gate["internal_runtime_source_count_for_prediction"] == 1, "Internal runtime cannot become biological evidence.", runtime_gate),
    ]
    return controls


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_report: dict[str, Any],
    dashboard: dict[str, Any],
    physical_calibration_inputs: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    replay_rows = [row for row in scoring_rows if row["panel_group"] == "V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY"]
    missing_candidate_rows = [row for row in scoring_rows if row["panel_group"] == "MISSING_WORD_CANDIDATE_CLEAN_ABSTAIN_CONTROL"]
    sentinel_rows = [
        row for row in scoring_rows
        if row["panel_group"] in {
            "MONOMERIC_GLOBULAR_SENTINEL",
            "ASSEMBLY_REQUIRED_SENTINEL",
            "MEMBRANE_TM_SENTINEL",
            "METAL_LIGAND_SENTINEL",
            "DISORDER_BOUNDARY_SENTINEL",
        }
    ]
    v74_repaired = sum(1 for row in replay_rows if row["accepted_supported"])
    sentinel_regressions = sum(1 for row in sentinel_rows if row["score_label"] != "supported")
    missing_candidate_clean = sum(1 for row in missing_candidate_rows if row["clean_abstain_supported"])
    missing_counts = Counter(row["missing_word_candidate"] for row in missing_candidate_rows if row.get("missing_word_candidate"))
    accepted_rows = [row for row in scoring_rows if row["acceptance_decision"] == "accepted"]
    operator_basis_stable = sum(
        1
        for row in accepted_rows
        if row["operator_basis_stability"] != "coefficient_assignment_sensitive"
    )
    cross_view_bound = sum(
        1
        for row in accepted_rows
        if row["cross_view_missing_view_count"] == 0
    )
    temporal_bound = sum(
        1
        for row in accepted_rows
        if row["temporal_binding"] != "selected_observable_temporal_conflict"
    )
    physical_backend_available = any(row["physical_backend_available"] for row in scoring_rows)
    physical_grounding_statuses = dict(Counter(row["physical_grounding_status"] for row in scoring_rows))
    real_physical_calibration_inputs_used = all(row["real_physical_calibration_inputs_used"] for row in scoring_rows)
    if any(control in {"coordinate_leakage_control", "internal_runtime_leakage_control"} for control in failed_controls):
        status = BATCH_BLOCKED_FOR_LEAKAGE
    elif failed_controls:
        status = BATCH_CONTROLS_FAILED
    elif metrics["failed_accepted_count"] > 0:
        status = BATCH_CORTEX_BROKEN
    elif v74_repaired != 33 or sentinel_regressions != 0 or missing_candidate_clean != len(missing_candidate_rows):
        status = BATCH_REPAIR_REQUIRED
    else:
        status = BATCH_PASSED
    cert = {
        "kind": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_CERTIFICATE_v0",
        "status": status,
        "status_options": [BATCH_PASSED, BATCH_DISCOVERED, BATCH_CORTEX_BROKEN, BATCH_REPAIR_REQUIRED],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "zero_failed_accepted_required": True,
        "self_decision_cortex_enabled": True,
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "sentinel_regressions": sentinel_regressions,
        "sentinel_regression_count": sentinel_regressions,
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "v74_multidomain_failures_repaired": v74_repaired,
        "v74_multidomain_failures_total": len(replay_rows),
        "missing_candidate_controls_total": len(missing_candidate_rows),
        "missing_candidate_controls_cleanly_abstained": missing_candidate_clean,
        "unknown_word_abstentions_by_word": dict(missing_counts),
        "top_missing_esperanto_word": missing_counts.most_common(1)[0][0] if missing_counts else None,
        "accepted_operator_basis_stable": operator_basis_stable,
        "accepted_cross_view_bound": cross_view_bound,
        "accepted_temporal_bound": temporal_bound,
        "dominance_law_for_acceptance": "single_dominant_learned_mechanism_bound_across_views",
        "coefficient_probe_mode": "endogenous_observed_operator_permutations_no_static_scale_range",
        "physical_basis_claim_allowed": False,
        "physics_grounding_status": "coarse_operator_heuristic_not_atomistic_physics",
        "physical_backend_available": physical_backend_available,
        "physical_grounding_statuses": physical_grounding_statuses,
        "real_physical_calibration_inputs_used": real_physical_calibration_inputs_used,
        "real_physical_calibration_kind": physical_calibration_inputs.get("kind"),
        "real_physical_calibration_row_count": physical_calibration_inputs.get("row_count"),
        "real_physical_calibration_hash": physical_calibration_inputs.get("calibration_hash"),
        "real_physical_calibration_source_dataset": physical_calibration_inputs.get("source_dataset"),
        "real_physical_calibration_source_dataset_sha256": physical_calibration_inputs.get("source_dataset_sha256"),
        "real_physical_calibration_source_coordinate_database": physical_calibration_inputs.get("source_coordinate_database"),
        "real_physical_calibration_input_type": physical_calibration_inputs.get("calibration_input_type"),
        "real_physical_calibration_target_native_excluded": physical_calibration_inputs.get("target_native_excluded_from_calibration"),
        "real_physical_calibration_target_native_contacts_used_before_prediction": physical_calibration_inputs.get("target_native_contacts_used_before_prediction"),
        "real_physical_calibration_coordinate_truth_used_as_prediction_input": physical_calibration_inputs.get("coordinate_truth_used_as_prediction_input"),
        "real_physical_calibration_leave_one_target_out": physical_calibration_inputs.get("leave_one_target_out_calibration"),
        "real_physical_calibration_observable_families": physical_calibration_inputs.get("observable_families"),
        "real_physical_calibration_fold_class_coverage": physical_calibration_inputs.get("fold_class_coverage"),
        "physical_next_required_capability": "calibrated OpenMM or equivalent force-field execution with topology/environment and independent physical observable holdouts",
        "candidate_grammars_implemented": False,
        "candidate_grammar_acceptance_role": "clean_abstain_until_revision_implements_grammar",
        "next_recommended_engine_revision": "E70_SECRETORY_DISULFIDE_REDOX_TOPOLOGY_GRAMMAR" if missing_counts.most_common(1) and missing_counts.most_common(1)[0][0] == "disulfide_secretory_redox_context" else None,
        "next_required_batch": "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200" if missing_counts.most_common(1) and missing_counts.most_common(1)[0][0] == "disulfide_secretory_redox_context" else "V76_NEXT_UNKNOWN_WORD_REPAIR_PANEL_200",
        "failed_accepted_by_failure_mode": failure_report["failed_accepted_by_failure_mode"],
        "dashboard": dashboard["shards"],
        "coordinate_truth_used_before_seal": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V75 is a self-decision cortex batch; zero failed accepted is required but not a solved-folding claim.",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "zero_failed_accepted_required": cert["zero_failed_accepted_required"],
        "unknown_word_abstentions": cert["unknown_word_abstentions"],
        "top_missing_esperanto_word": cert["top_missing_esperanto_word"],
        "status": cert["status"],
        "real_physical_calibration_inputs_used": cert["real_physical_calibration_inputs_used"],
        "real_physical_calibration_hash": cert["real_physical_calibration_hash"],
        "claim_allowed": cert["claim_allowed"],
        "claim_blocked_reason": cert["claim_blocked_reason"],
    }


def _append_claim_ledger(row: dict[str, Any]) -> Path:
    path = LEDGER_ROOT / "claim_ledger_v0.json"
    ledger = _read_json(path, "campaign claim ledger") if path.exists() else {"kind": "V75_CLAIM_LEDGER_v0", "campaign_id": CAMPAIGN_ID, "rows": []}
    rows = [existing for existing in ledger.get("rows", []) if isinstance(existing, dict) and existing.get("batch_id") != BATCH_ID]
    rows.append(row)
    ledger["kind"] = ledger.get("kind", "V75_CLAIM_LEDGER_v0")
    ledger["campaign_id"] = CAMPAIGN_ID
    ledger["rows"] = rows
    return _write_json(path, ledger)


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V75 E69 Multidomain Repair And Zero Failed Accept Cortex",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Accepted count: `{cert['accepted_count']}`",
        f"Accepted supported: `{cert['accepted_supported']}`",
        f"Clean abstain supported: `{cert['clean_abstain_supported']}`",
        f"Blocked supported: `{cert['blocked_supported']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Unknown word abstentions: `{cert['unknown_word_abstentions']}`",
        f"Top missing Esperanto word: `{cert['top_missing_esperanto_word']}`",
        f"Sentinel regressions: `{cert['sentinel_regressions']}`",
        f"Controls passed: `{cert['controls_passed']}`",
        f"Zero failed accepted required: `{cert['zero_failed_accepted_required']}`",
        f"V74 multidomain repaired: `{cert['v74_multidomain_failures_repaired']} / {cert['v74_multidomain_failures_total']}`",
        f"Unresolved tracking controls cleanly abstained: `{cert['missing_candidate_controls_cleanly_abstained']} / {cert['missing_candidate_controls_total']}`",
        f"Real physical calibration inputs used: `{cert['real_physical_calibration_inputs_used']}`",
        f"Real physical calibration rows: `{cert['real_physical_calibration_row_count']}`",
        f"Real physical calibration hash: `{cert['real_physical_calibration_hash']}`",
        f"Physical grounding statuses: `{cert['physical_grounding_statuses']}`",
        f"Physical basis claim allowed: `{cert['physical_basis_claim_allowed']}`",
        "",
        "## Dashboard",
        "",
        "| group | targets_total | accepted_count | accepted_supported | clean_abstain_supported | failed_accepted | accepted_accuracy | coverage | top_missing_esperanto_word |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for group, row in cert["dashboard"].items():
        lines.append(
            f"| `{group}` | `{row['targets_total']}` | `{row['accepted_count']}` | `{row['accepted_supported']}` | "
            f"`{row['clean_abstain_supported']}` | `{row['failed_accepted']}` | `{row['accepted_accuracy']}` | `{row['coverage']}` | `{row['top_missing_esperanto_word']}` |"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v75(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _read_json(V74_CERT, "V74 certificate")
    _read_json(E69_CERT, "E69 certificate")
    _reset_generated_outputs(out_dir)
    physical_calibration_inputs = write_real_physical_calibration_inputs(
        REAL_COORDINATE_BENCHMARK,
        PHYSICAL_CALIBRATION_INPUTS,
    )
    targets = _select_targets()
    target_manifest = _target_manifest(targets)
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v75_self_deciding_acceptance_cortex_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v75_e69_engine_declaration.json", engine_declaration)

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
            focus_regions=[{"name": "V75 self-decision cortex full-chain scan", "span": f"1-{target['sequence_length']}"}],
            perturbations=_perturbations_for_target(target),
            physical_calibration_inputs=physical_calibration_inputs,
        )
        holdout = _holdout(target, packet)
        score = _score(packet, holdout, target)
        shuffled_packet = build_sealed_simulation_packet(
            target_id=f"{target['target_id']}_SHUFFLED_CONTROL",
            target_name=f"{target['entry_id']} shuffled sequence control",
            sequence=shuffled_sequence(target["sequence"]),
            sources=[
                {
                    "source_id": f"V75_SHUFFLED_SEQUENCE_{stable_hash({'target': target['target_id']})[:12]}",
                    "source_class": "pure_non_coordinate",
                    "source_role": "prediction_input",
                    "coordinate_derived": False,
                    "internal_runtime_source": False,
                    "spatial_proxy": False,
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and context are withheld.",
                }
            ],
            perturbations=[],
            physical_calibration_inputs=physical_calibration_inputs,
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
    dashboard = _dashboard(scoring_rows)
    controls = _controls(target_manifest, scoring_rows, shuffled_packets, physical_calibration_inputs)
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
        failure_report=failure_report,
        dashboard=dashboard,
        physical_calibration_inputs=physical_calibration_inputs,
    )
    claim_row = _claim_row(cert)
    scoring_path = _write_json(DATA_ROOT / "v75_self_deciding_acceptance_cortex_scoring_report.json", {"kind": "V75_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v75_self_deciding_acceptance_cortex_failure_report.json", failure_report)
    dashboard_path = _write_json(DATA_ROOT / "v75_self_deciding_acceptance_cortex_dashboard.json", dashboard)
    data_cert_path = _write_json(DATA_ROOT / "v75_self_deciding_acceptance_cortex_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v75_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_calibration_path = _write_json(out_dir / "v75_real_physical_calibration_inputs.json", physical_calibration_inputs)
    cert_path = _write_json(out_dir / "v75_self_deciding_acceptance_cortex_certificate.json", cert)
    report_path = out_dir / "V75_SELF_DECIDING_ACCEPTANCE_CORTEX_REPORT.md"
    _write_report(report_path, cert)
    return {
        "target_manifest": DATA_ROOT / "v75_self_deciding_acceptance_cortex_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v75_e69_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "dashboard": dashboard_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
        "physical_calibration_inputs": PHYSICAL_CALIBRATION_INPUTS,
        "output_physical_calibration_inputs": out_calibration_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V75 self-deciding acceptance cortex.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v75(args.out_dir)
    cert = _read_json(paths["certificate"], "V75 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_count": cert["accepted_count"],
        "accepted_supported": cert["accepted_supported"],
        "clean_abstain_supported": cert["clean_abstain_supported"],
        "blocked_supported": cert["blocked_supported"],
        "failed_accepted": cert["failed_accepted_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "unknown_word_abstentions": cert["unknown_word_abstentions"],
        "top_missing_esperanto_word": cert["top_missing_esperanto_word"],
        "sentinel_regressions": cert["sentinel_regressions"],
        "controls_passed": cert["controls_passed"],
        "zero_failed_accepted_required": cert["zero_failed_accepted_required"],
        "v74_multidomain_failures_repaired": cert["v74_multidomain_failures_repaired"],
        "missing_candidate_controls_cleanly_abstained": cert["missing_candidate_controls_cleanly_abstained"],
        "next_recommended_engine_revision": cert["next_recommended_engine_revision"],
        "next_required_batch": cert["next_required_batch"],
        "physical_backend_available": cert["physical_backend_available"],
        "physical_grounding_statuses": cert["physical_grounding_statuses"],
        "real_physical_calibration_inputs_used": cert["real_physical_calibration_inputs_used"],
        "real_physical_calibration_row_count": cert["real_physical_calibration_row_count"],
        "real_physical_calibration_hash": cert["real_physical_calibration_hash"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == BATCH_PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
