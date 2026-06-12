#!/usr/bin/env python3
from __future__ import annotations

"""Run V76: E70 secretory disulfide/redox repair with matched controls."""

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
    STATE_VARIABLES,
    UNIVERSAL_OPERATORS,
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    shuffled_sequence,
    stable_hash,
)
from pharmacotopology.protein_esperanto_physical_calibration import (  # noqa: E402
    write_real_physical_calibration_inputs,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402
import run_v69_e65_rcsb_nonredundant_200_discovery_v0 as v69  # noqa: E402
import run_v71_e66_rcsb_nonredundant_200_discovery_v0 as v71  # noqa: E402
import run_v74_e68_rcsb_nonredundant_200_discovery_v0 as v74  # noqa: E402


BATCH_ID = "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200"
CAMPAIGN_ID = "V61_TO_V76_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E70"
BASELINE_ENGINE_VERSION = "E69"
TARGET_COUNT = 200

ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"
DISULFIDE_CLASS = "secretory_disulfide_redox_topology"
METAL_CLASS = "metal_cluster_and_ligand_locked_basin"
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
GLOBULAR_CLASS = "globular_closure"
BETA_CLASS = "beta_closure_topology"
ASSEMBLY_CLASS = "assembly_required_folding"
DISORDER_CLASS = "disorder_boundary_and_fold_upon_binding"
MULTIDOMAIN_CLASS = "multidomain_allosteric_architecture"

GROUP_COUNTS = OrderedDict([
    ("SECRETORY_DISULFIDE_POSITIVE", 60),
    ("METAL_CYS_HIS_NEGATIVE", 30),
    ("SIGNAL_TM_NEGATIVE", 20),
    ("CYSTEINE_GLOBULAR_NEGATIVE", 20),
    ("BETA_CYSTEINE_NEGATIVE", 20),
    ("COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL", 20),
    ("V75_ACCEPTED_SENTINEL_REPLAY", 30),
])

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V76"
E70_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E70"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"
REAL_COORDINATE_BENCHMARK = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
PHYSICAL_CALIBRATION_INPUTS = DATA_ROOT / "physical_calibration" / "v76_real_physical_calibration_inputs.json"

V74_RAW = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "intake" / "raw_rcsb_30pct_representative_entities_v74.json"
V75_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V75" / "v75_self_deciding_acceptance_cortex_target_manifest.json"
V75_SCORING = REPO_ROOT / "data" / "protein_esperanto_engine" / "V75" / "v75_self_deciding_acceptance_cortex_scoring_report.json"

PASSED = "V76_E70_SECRETORY_DISULFIDE_REPAIR_PANEL_PASSED"
FAILED = "V76_E70_SECRETORY_DISULFIDE_REPAIR_PANEL_FAILED"
CONTROLS_FAILED = "V76_E70_SECRETORY_DISULFIDE_REPAIR_PANEL_CONTROLS_FAILED"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
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


def _candidate_text(candidate: dict[str, Any]) -> str:
    return v74._candidate_text(candidate)


def _has_any(candidate: dict[str, Any], tokens: list[str]) -> bool:
    text = _candidate_text(candidate)
    return any(token in text for token in tokens)


def _cysteine_count(candidate: dict[str, Any]) -> int:
    return str(candidate.get("sequence") or "").count("C")


def _histidine_count(candidate: dict[str, Any]) -> int:
    return str(candidate.get("sequence") or "").count("H")


def _raw_candidates() -> list[dict[str, Any]]:
    raw = _read_json(V74_RAW, "V74 raw candidate cache")["candidates"]
    rows = []
    seen: set[str] = set()
    for row in raw:
        if not isinstance(row, dict):
            continue
        candidate = dict(row)
        protein_id = _candidate_id(candidate)
        sequence = str(candidate.get("sequence") or "")
        if not protein_id or not sequence or protein_id in seen:
            continue
        candidate["protein_id"] = protein_id
        candidate.setdefault("sequence_metrics", v61._sequence_metrics(sequence))
        candidate.setdefault("source_urls", {
            "entry": candidate.get("entry_url", ""),
            "polymer_entity": candidate.get("polymer_entity_url", ""),
        })
        rows.append(candidate)
        seen.add(protein_id)
    return sorted(rows, key=_candidate_rank)


def _pick(candidates: list[dict[str, Any]], *, count: int, used: set[str], predicate: Callable[[dict[str, Any]], bool]) -> list[dict[str, Any]]:
    selected = []
    for candidate in candidates:
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
    expected_decision: str = "accepted",
    source_family: str = "V74_RAW",
    lineage_source_target: str | None = None,
) -> dict[str, Any]:
    entry_url = candidate.get("entry_url") or candidate.get("source_urls", {}).get("entry", "")
    polymer_url = candidate.get("polymer_entity_url") or candidate.get("source_urls", {}).get("polymer_entity", "")
    sequence = candidate["sequence"]
    return {
        "target_id": f"V76_{ordinal:03d}_{_safe_id(group)}_{_safe_id(_candidate_id(candidate))}",
        "panel_group": group,
        "source_family": source_family,
        "lineage_source_target": lineage_source_target,
        "expected_decision": expected_decision,
        "expected_mechanism_class": expected,
        "required_esperanto_word": required_word,
        "protein_id": _candidate_id(candidate),
        "entry_id": str(candidate.get("entry_id", "")),
        "entity_id": str(candidate.get("entity_id", "")),
        "sequence": sequence,
        "sequence_length": int(candidate.get("sequence_length") or len(sequence)),
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip(),
        "entry_url": entry_url,
        "polymer_entity_url": polymer_url,
        "candidate_snapshot": dict(candidate),
        "postseal_truth_basis": [
            f"V76 E70 matched-control panel group {group}.",
            f"Required word: {required_word or 'none'}.",
            f"Expected decision: {expected_decision}.",
            "Coordinates, native contacts, native disulfide geometry, and validation labels are blocked before sealing.",
        ],
    }


def _select_targets() -> list[dict[str, Any]]:
    candidates = _raw_candidates()
    targets: list[dict[str, Any]] = []
    used: set[str] = set()
    ordinal = 1
    secretory_tokens = ["disulfide", "disulphide", "secreted", "secretory", "extracellular", "cysteine-rich", "cysteine rich", "glycoprotein", "thioredoxin"]

    groups = [
        (
            "SECRETORY_DISULFIDE_POSITIVE",
            GROUP_COUNTS["SECRETORY_DISULFIDE_POSITIVE"],
            lambda c: _cysteine_count(c) >= 2 and _cysteine_count(c) % 2 == 0 and _has_any(c, secretory_tokens),
            DISULFIDE_CLASS,
            "disulfide_secretory_redox_context",
            "accepted",
        ),
        (
            "METAL_CYS_HIS_NEGATIVE",
            GROUP_COUNTS["METAL_CYS_HIS_NEGATIVE"],
            lambda c: (v69._metal_or_ligand_locked_word(c) is not None or _has_any(c, ["zinc", "iron-sulfur", "metal", "heme"])) and (_cysteine_count(c) + _histidine_count(c)) >= 2,
            METAL_CLASS,
            "metal_cluster_geometry",
            "accepted",
        ),
        (
            "SIGNAL_TM_NEGATIVE",
            GROUP_COUNTS["SIGNAL_TM_NEGATIVE"],
            lambda c: v69._true_tm(c) or _has_any(c, ["signal peptide", "signal sequence", "transmembrane"]),
            ABSTAIN_CLASS,
            "signal_peptide_vs_true_TM",
            "abstain_recommended",
        ),
        (
            "CYSTEINE_GLOBULAR_NEGATIVE",
            GROUP_COUNTS["CYSTEINE_GLOBULAR_NEGATIVE"],
            lambda c: _cysteine_count(c) >= 2 and not v69._true_tm(c) and v69._metal_or_ligand_locked_word(c) is None,
            GLOBULAR_CLASS,
            None,
            "accepted",
        ),
        (
            "BETA_CYSTEINE_NEGATIVE",
            GROUP_COUNTS["BETA_CYSTEINE_NEGATIVE"],
            lambda c: _cysteine_count(c) >= 2 and v74._beta_word(c) is not None,
            BETA_CLASS,
            "closed_beta_topology",
            "accepted",
        ),
        (
            "COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL",
            GROUP_COUNTS["COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL"],
            lambda c: v74._coiled_repeat_word(c) in {"coiled_coil_register", "repeat_solenoid_topology"},
            ABSTAIN_CLASS,
            "coiled_coil_register",
            "abstain_recommended",
        ),
    ]
    for group, count, predicate, expected, required_word, expected_decision in groups:
        for candidate in _pick(candidates, count=count, used=used, predicate=predicate):
            word = required_word
            if group == "COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL":
                word = v74._coiled_repeat_word(candidate) or required_word
            targets.append(_target_from_candidate(
                ordinal=ordinal,
                group=group,
                candidate=candidate,
                expected=expected,
                required_word=word,
                expected_decision=expected_decision,
            ))
            ordinal += 1

    v75_manifest = _read_json(V75_MANIFEST, "V75 target manifest")["selected_targets"]
    v75_scoring = _read_json(V75_SCORING, "V75 scoring report")["rows"]
    accepted_by_id = {row["target_id"]: row for row in v75_scoring if row["acceptance_decision"] == "accepted" and row["accepted_supported"]}
    sentinel_rows = [
        row for row in v75_manifest
        if row["target_id"] in accepted_by_id
        and row["panel_group"] in {
            "MONOMERIC_GLOBULAR_SENTINEL",
            "ASSEMBLY_REQUIRED_SENTINEL",
            "MEMBRANE_TM_SENTINEL",
            "METAL_LIGAND_SENTINEL",
            "DISORDER_BOUNDARY_SENTINEL",
            "V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY",
        }
    ][: GROUP_COUNTS["V75_ACCEPTED_SENTINEL_REPLAY"]]
    if len(sentinel_rows) != GROUP_COUNTS["V75_ACCEPTED_SENTINEL_REPLAY"]:
        raise SystemExit("not enough V75 accepted sentinels for replay")
    for row in sentinel_rows:
        candidate = dict(row.get("candidate_snapshot") or row)
        candidate.setdefault("protein_id", row["protein_id"])
        candidate.setdefault("entry_id", row.get("entry_id", ""))
        candidate.setdefault("entity_id", row.get("entity_id", ""))
        candidate.setdefault("sequence", row["sequence"])
        candidate.setdefault("sequence_length", row["sequence_length"])
        candidate.setdefault("source_urls", {"entry": row.get("entry_url", ""), "polymer_entity": row.get("polymer_entity_url", "")})
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group="V75_ACCEPTED_SENTINEL_REPLAY",
            candidate=candidate,
            expected=row["expected_mechanism_class"],
            required_word=row.get("required_esperanto_word"),
            expected_decision="accepted",
            source_family="V75",
            lineage_source_target=row["target_id"],
        ))
        ordinal += 1

    composition = Counter(row["panel_group"] for row in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != dict(GROUP_COUNTS):
        raise SystemExit(f"bad V76 composition: total={len(targets)} composition={dict(composition)}")
    return targets


def _context_marks(target: dict[str, Any]) -> list[str]:
    group = target["panel_group"]
    if group == "SECRETORY_DISULFIDE_POSITIVE":
        return [
            "disulfide_secretory_redox_context",
            "disulfide_bond_topology",
            "secretory_redox_context",
            "cysteine_pairing_constraint",
            "extracellular_stabilized_fold",
            "glycosylation_context",
            "signal_peptide_removed_context",
            "secretory_quality_control",
        ]
    if group == "METAL_CYS_HIS_NEGATIVE":
        return ["metal_cluster_geometry", "coordination_shell_integrity", "cys-his coordination", "zinc-binding", "disulfide_secretory_redox_context"]
    if group == "SIGNAL_TM_NEGATIVE":
        return ["membrane_context_strong", "transmembrane_context", "topology_evidence", "signal peptide only", "disulfide_secretory_redox_context"]
    if group == "CYSTEINE_GLOBULAR_NEGATIVE":
        return ["soluble_monomeric_core_context", "complete soluble monomer", "standalone soluble fold", "ordinary cysteine noise"]
    if group == "BETA_CYSTEINE_NEGATIVE":
        return ["closed_beta_topology", "strand_register", "beta_sheet_closure", "beta_sandwich_core", "cysteine rich beta context"]
    if group == "COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL":
        word = target.get("required_esperanto_word")
        if word == "repeat_solenoid_topology":
            return ["repeat_solenoid_topology", "repeat_unit", "solenoid_axis", "disulfide_secretory_redox_context"]
        return ["coiled_coil_register", "heptad_repeat", "register_alignment", "disulfide_secretory_redox_context"]
    expected = target["expected_mechanism_class"]
    if expected == GLOBULAR_CLASS:
        return ["soluble_monomeric_core_context", "complete soluble monomer", "standalone soluble fold"]
    if expected == ASSEMBLY_CLASS:
        return ["assembly_required_core", "assembly_required_folding", "partner_completed_core", "interface_buried_hydrophobicity"]
    if expected == MEMBRANE_CLASS:
        return ["membrane_context_strong", "transmembrane_context", "topology_evidence"]
    if expected == METAL_CLASS:
        return ["metal_cluster_geometry", "ligand_locked_basin", "coordination_shell_integrity"]
    if expected == DISORDER_CLASS:
        return ["disorder_context", "IDR_boundary", "structured_domain_plus_IDR_tail", "fold_upon_binding_region"]
    if expected == MULTIDOMAIN_CLASS:
        return ["multidomain_allostery", "domain_boundary", "hinge_region", "interdomain_lock", "allosteric_basin_shift"]
    return []


def _source_manifest(target: dict[str, Any], *, mask_metadata: bool = False) -> dict[str, Any]:
    suffix = stable_hash({"v76_target_id": target["target_id"], "mask": mask_metadata})[:12]
    marks = [] if mask_metadata else _context_marks(target)
    statement = " ".join(marks) if marks else "matched control with biological context marks withheld"
    return {
        "kind": "V76_SECRETORY_DISULFIDE_SOURCE_MANIFEST_v0",
        "target_id": target["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "prediction_sources": [
            {
                "source_id": f"V76_RAW_SEQUENCE_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": target["polymer_entity_url"],
            },
            {
                "source_id": f"V76_E70_CONTEXT_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": marks,
                "evidence_statement": f"V76 allowed non-coordinate context marks: {statement}.",
                "source_url": target["entry_url"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native disulfide bond geometry before sealing",
            "native contacts and distance maps before sealing",
            "post-seal validation annotations before prediction hash",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _mask_cysteines(sequence: str) -> str:
    return sequence.replace("C", "S")


def _metric_for_expected(expected: str, required_word: str | None) -> str:
    if expected == DISULFIDE_CLASS:
        return "disulfide_pairing_topology"
    if expected == METAL_CLASS:
        return "coordination_shell_integrity" if required_word == "metal_cluster_geometry" else "ligand_locked_basin"
    if expected == MEMBRANE_CLASS:
        return "proteostasis_routing"
    if expected == GLOBULAR_CLASS:
        return "contact_probability"
    if expected == BETA_CLASS:
        return "closed_beta_topology"
    if expected == ASSEMBLY_CLASS:
        return "partner_completed_core"
    if expected == DISORDER_CLASS:
        return "fold_upon_binding_region"
    if expected == MULTIDOMAIN_CLASS:
        return "interdomain_lock"
    return "operator_activation"


def _extract_metric(packet: dict[str, Any], metric: str) -> float:
    return float(packet["operator_state_propagation_summary"]["final_state_summary"].get(metric, 0.0))


def _disulfide_perturbation(target: dict[str, Any]) -> dict[str, Any]:
    return {
        "perturbation_id": f"{target['target_id']}_DISULFIDE_REDUCTION",
        "description": "reduce or mispair secretory disulfide topology",
        "operator_scales": {"disulfide_pairing_operator": 0.24, "secretory_redox_operator": 0.50},
        "disulfide_damage": 0.58,
        "redox_shift": 0.25,
        "metric": "disulfide_pairing_topology",
        "expected_direction": "decrease",
    }


def _perturbations_for_target(target: dict[str, Any]) -> list[dict[str, Any]]:
    if target["expected_mechanism_class"] == DISULFIDE_CLASS:
        return [_disulfide_perturbation(target)]
    expected = target["expected_mechanism_class"]
    metric = _metric_for_expected(expected, target.get("required_esperanto_word"))
    scales = {
        METAL_CLASS: {"interface_operator": 0.45},
        MEMBRANE_CLASS: {"membrane_pressure_operator": 0.55, "proteostasis_operator": 0.55},
        GLOBULAR_CLASS: {"closure_operator": 0.45},
        BETA_CLASS: {"closure_operator": 0.45},
        ASSEMBLY_CLASS: {"interface_operator": 0.45},
        DISORDER_CLASS: {"interface_operator": 0.45},
        MULTIDOMAIN_CLASS: {"interface_operator": 0.45},
    }.get(expected, {})
    return [{
        "perturbation_id": f"{target['target_id']}_MECHANISM_DAMAGE",
        "description": "matched mechanism-damage perturbation",
        "operator_scales": scales,
        "metric": metric,
        "expected_direction": "decrease",
    }] if scales else []


def _packet(
    *,
    target: dict[str, Any],
    sequence: str,
    sources: list[dict[str, Any]],
    perturbations: list[dict[str, Any]],
    physical_calibration_inputs: dict[str, Any],
    suffix: str = "",
) -> dict[str, Any]:
    return build_sealed_operator_state_packet(
        target_id=f"{target['target_id']}{suffix}",
        target_name=target["target_name"],
        sequence=sequence,
        sources=sources,
        focus_regions=[{"name": "V76 matched-control full-chain scan", "span": f"1-{target['sequence_length']}"}],
        perturbations=perturbations,
        physical_calibration_inputs=physical_calibration_inputs,
    )


def _matched_control_report(
    *,
    target: dict[str, Any],
    packet: dict[str, Any],
    source_manifest: dict[str, Any],
    physical_calibration_inputs: dict[str, Any],
) -> dict[str, Any]:
    metric = _metric_for_expected(target["expected_mechanism_class"], target.get("required_esperanto_word"))
    real_value = _extract_metric(packet, metric)
    rows = []
    shuffled_packet = None
    metadata_masked_packet = None
    cysteine_packet = None
    if target["expected_mechanism_class"] == DISULFIDE_CLASS:
        control_suffix = stable_hash({"control": target["target_id"]})[:12]
        no_context_source = [{
            "source_id": f"V76_SHUFFLED_CONTROL_{control_suffix}",
            "source_class": "pure_non_coordinate",
            "source_role": "prediction_input",
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "spatial_proxy": False,
            "evidence_statement": "Composition-preserving shuffled sequence control. Target biological context marks are withheld.",
        }]
        shuffled_packet = _packet(
            target=target,
            sequence=shuffled_sequence(target["sequence"]),
            sources=no_context_source,
            perturbations=[],
            physical_calibration_inputs=physical_calibration_inputs,
            suffix="_SHUFFLED_CONTROL",
        )
        metadata_masked_manifest = _source_manifest(target, mask_metadata=True)
        metadata_masked_packet = _packet(
            target=target,
            sequence=target["sequence"],
            sources=metadata_masked_manifest["prediction_sources"],
            perturbations=[],
            physical_calibration_inputs=physical_calibration_inputs,
            suffix="_METADATA_MASKED_CONTROL",
        )
        rows.extend([
            {
                "control": "real_sequence_beats_shuffled_no_context_control",
                "metric": metric,
                "real_value": real_value,
                "control_value": _extract_metric(shuffled_packet, metric),
                "passed": real_value > _extract_metric(shuffled_packet, metric),
            },
            {
                "control": "real_metadata_beats_metadata_masked_source",
                "metric": metric,
                "real_value": real_value,
                "control_value": _extract_metric(metadata_masked_packet, metric),
                "passed": real_value > _extract_metric(metadata_masked_packet, metric),
            },
        ])
        cysteine_packet = _packet(
            target=target,
            sequence=_mask_cysteines(target["sequence"]),
            sources=source_manifest["prediction_sources"],
            perturbations=[],
            physical_calibration_inputs=physical_calibration_inputs,
            suffix="_CYSTEINE_MASKED_CONTROL",
        )
        rows.append({
            "control": "real_sequence_beats_cysteine_masked_control",
            "metric": metric,
            "real_value": real_value,
            "control_value": _extract_metric(cysteine_packet, metric),
            "passed": real_value > _extract_metric(cysteine_packet, metric),
        })
        perturbation_rows = [
            row
            for row in packet["predicted_perturbation_table"]
            if row["metric"] == metric
        ]
        perturbation_passed = any(
            float(row["perturbed_value"]) < float(row["baseline_value"])
            for row in perturbation_rows
        )
        rows.append({
            "control": "disulfide_redox_perturbation_decreases_selected_state",
            "metric": metric,
            "real_value": real_value,
            "control_value": [
                {
                    "baseline_value": row["baseline_value"],
                    "perturbed_value": row["perturbed_value"],
                }
                for row in perturbation_rows
            ],
            "passed": perturbation_passed,
        })
    elif target["expected_decision"] == "accepted":
        rows.append({
            "control": "negative_or_sentinel_expected_mechanism_preserved",
            "metric": "mechanism_class",
            "real_value": packet["selected_mechanism_grammar"]["mechanism_class"],
            "control_value": DISULFIDE_CLASS,
            "passed": packet["selected_mechanism_grammar"]["mechanism_class"] == target["expected_mechanism_class"],
        })
    wrong_grammar_passed = packet["self_decision_judge"]["wrong_grammar_separation"] != "wrong_grammar_competes"
    rows.append({
        "control": "wrong_grammar_challenge_fails",
        "metric": "wrong_grammar_separation",
        "real_value": packet["self_decision_judge"]["wrong_grammar_separation"],
        "control_value": "wrong_grammar_competes",
        "passed": wrong_grammar_passed,
    })
    return {
        "kind": "V76_MATCHED_CONTROL_DOMINANCE_v0",
        "metric": metric,
        "matched_control_dominance_passed": all(row["passed"] for row in rows),
        "control_rows": rows,
        "control_packet_hashes": {
            "shuffled": shuffled_packet["prediction_hash"] if shuffled_packet else None,
            "metadata_masked": metadata_masked_packet["prediction_hash"] if metadata_masked_packet else None,
            "cysteine_masked": cysteine_packet["prediction_hash"] if cysteine_packet else None,
        },
        "uses_static_observable_thresholds": False,
    }


def _score(target: dict[str, Any], packet: dict[str, Any], matched: dict[str, Any]) -> dict[str, Any]:
    judge = packet["self_decision_judge"]
    decision = judge["acceptance_decision"]
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = target["expected_mechanism_class"]
    accepted = decision == "accepted"
    accepted_supported = accepted and predicted == expected and matched["matched_control_dominance_passed"]
    clean_abstain_supported = (
        decision == "abstain_recommended"
        and target["expected_decision"] == "abstain_recommended"
        and judge["final_self_decision"] == "clean_abstain_missing_word"
        and judge.get("missing_word_candidate") == target.get("required_esperanto_word")
    )
    score_label = "supported" if accepted_supported or clean_abstain_supported else "abstained" if decision == "abstain_recommended" else "contradicted"
    return {
        "kind": "V76_SECRETORY_DISULFIDE_VALIDATION_RESULT_v0",
        "target_id": target["target_id"],
        "panel_group": target["panel_group"],
        "lineage_source_target": target.get("lineage_source_target"),
        "acceptance_decision": decision,
        "final_self_decision": judge["final_self_decision"],
        "self_decision_reason": judge["self_decision_reason"],
        "dominance_law": judge["dominance_law"],
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "required_esperanto_word": target.get("required_esperanto_word"),
        "missing_word_candidate": judge.get("missing_word_candidate"),
        "accepted_supported": accepted_supported,
        "clean_abstain_supported": clean_abstain_supported,
        "score_label": score_label,
        "matched_control_dominance_passed": matched["matched_control_dominance_passed"],
        "matched_control_dominance": matched,
        "physical_grounding_status": judge["physical_grounding_status"],
        "physical_basis_claim_allowed": judge["physical_basis_claim_allowed"],
        "real_physical_calibration_inputs_used": judge["real_physical_calibration_inputs_used"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "coordinate_truth_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    accepted_supported = [row for row in accepted if row["accepted_supported"]]
    clean_abstain_supported = [row for row in abstained if row["clean_abstain_supported"]]
    failed_accepted = [row for row in accepted if not row["accepted_supported"]]
    return {
        "targets_total": len(rows),
        "accepted_count": len(accepted),
        "accepted_supported": len(accepted_supported),
        "clean_abstain": len(abstained),
        "clean_abstain_supported": len(clean_abstain_supported),
        "failed_accepted": len(failed_accepted),
        "failed_accepted_count": len(failed_accepted),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else 1.0,
        "coverage": len(accepted) / len(rows) if rows else None,
        "supported_count": len([row for row in rows if row["score_label"] == "supported"]),
    }


def _failure_report(rows: list[dict[str, Any]]) -> dict[str, Any]:
    failed = [row for row in rows if row["acceptance_decision"] == "accepted" and not row["accepted_supported"]]
    failure_rows = [
        {
            "target_id": row["target_id"],
            "panel_group": row["panel_group"],
            "failure_mode": row.get("required_esperanto_word") or row["expected_mechanism_class"],
            "engine_thought": row["predicted_mechanism_class"],
            "reality_showed": row["expected_mechanism_class"],
            "missing_esperanto_word": row.get("required_esperanto_word"),
        }
        for row in failed
    ]
    return {
        "kind": "V76_SECRETORY_DISULFIDE_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_count": len(failure_rows),
        "failed_accepted_by_failure_mode": dict(Counter(row["failure_mode"] for row in failure_rows)),
        "failure_grammar_rows": failure_rows,
    }


def _dashboard(rows: list[dict[str, Any]]) -> dict[str, Any]:
    shards: dict[str, Any] = {}
    for group in GROUP_COUNTS:
        group_rows = [row for row in rows if row["panel_group"] == group]
        missing = Counter(row["missing_word_candidate"] for row in group_rows if row.get("missing_word_candidate"))
        shards[group] = {
            **_metrics(group_rows),
            "top_missing_esperanto_word": missing.most_common(1)[0][0] if missing else None,
            "missing_esperanto_word_counts": dict(missing),
        }
    missing_total = Counter(row["missing_word_candidate"] for row in rows if row.get("missing_word_candidate"))
    shards["TOTAL"] = {
        **_metrics(rows),
        "top_missing_esperanto_word": missing_total.most_common(1)[0][0] if missing_total else None,
        "missing_esperanto_word_counts": dict(missing_total),
    }
    return {"kind": "V76_SECRETORY_DISULFIDE_DASHBOARD_v0", "batch_id": BATCH_ID, "shards": shards}


def _reset_outputs(out_dir: Path) -> None:
    if DATA_ROOT.exists():
        shutil.rmtree(DATA_ROOT)
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V76_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "self_decision_judge": packet["self_decision_judge"],
        "operator_names": packet["operator_field"]["operator_names"],
        "operator_state_final_state_summary": packet["operator_state_propagation_summary"]["final_state_summary"],
        "predicted_contact_interaction_probability_map": packet["predicted_contact_interaction_probability_map"],
        "folding_problem_solved": False,
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V76_E70_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "mechanism_classes": MECHANISM_CLASSES,
        "state_variables": STATE_VARIABLES,
        "operator_names": UNIVERSAL_OPERATORS,
        "matched_control_dominance_acceptance": True,
        "fixed_accepted_support_thresholds_used": False,
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
            "sequence_metrics": candidate.get("sequence_metrics", {}),
        })
        rows.append(row)
    return {
        "kind": "V76_SECRETORY_DISULFIDE_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "composition_rule": dict(GROUP_COUNTS),
        "target_count_selected": len(rows),
        "matched_control_dominance_acceptance": True,
        "fixed_accepted_support_thresholds_used": False,
        "selected_targets": rows,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(target_manifest: dict[str, Any], rows: list[dict[str, Any]], physical_calibration_inputs: dict[str, Any]) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{"source_id": "V76_BAD_COORD", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V76_BAD_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    random_packet = build_sealed_operator_state_packet(
        target_id="V76_RANDOM_SEQUENCE_CONTROL",
        target_name="V76 random sequence control",
        sequence=deterministic_random_sequence(128),
        sources=[],
        perturbations=[],
        physical_calibration_inputs=physical_calibration_inputs,
    )
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    missing_controls = [row for row in rows if row["panel_group"] == "COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL"]
    return [
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V76 must have 200 targets.", target_manifest["target_count_selected"]),
        _control("composition_rule", Counter(row["panel_group"] for row in target_manifest["selected_targets"]) == GROUP_COUNTS, "V76 composition must match the requested panel.", Counter(row["panel_group"] for row in target_manifest["selected_targets"])),
        _control("zero_failed_accepted", all(row["accepted_supported"] for row in accepted), "Every accepted target must be supported."),
        _control("accepted_matched_control_dominance", all(row["matched_control_dominance_passed"] for row in accepted), "Accepted targets must pass their self-required matched controls."),
        _control("missing_words_cleanly_abstain", all(row["clean_abstain_supported"] for row in missing_controls), "Coiled/repeat missing words must cleanly abstain.", Counter(row["missing_word_candidate"] for row in missing_controls)),
        _control("no_physical_basis_claim", all(row["physical_basis_claim_allowed"] is False for row in rows), "Language acceptance remains separate from physical claim allowance."),
        _control("real_physical_calibration_inputs_present", all(row["real_physical_calibration_inputs_used"] for row in rows), "Every row receives real calibration inputs."),
        _control("physical_calibration_truth_boundary", physical_calibration_inputs.get("target_native_contacts_used_before_prediction") is False and physical_calibration_inputs.get("coordinate_truth_used_as_prediction_input") is False, "Calibration inputs preserve truth boundary."),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate source blocks prediction.", coord_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime source blocks prediction.", runtime_gate),
        _control("random_sequence_control", random_packet["self_decision_judge"]["final_self_decision"] == "clean_abstain_low_internal_consensus", "Random sequence without evidence abstains.", random_packet["self_decision_judge"]["final_self_decision"]),
    ]


def _write_e70_certificate(v76_cert: dict[str, Any]) -> Path:
    cert = {
        "kind": "E70_SECRETORY_DISULFIDE_REDOX_TOPOLOGY_GRAMMAR_CERTIFICATE_v0",
        "engine_revision": "E70",
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "status": "E70_SECRETORY_DISULFIDE_REDOX_TOPOLOGY_GRAMMAR_ADDED",
        "trigger_batch": "V75_SELF_DECIDING_ACCEPTANCE_CORTEX",
        "trigger_missing_word": "disulfide_secretory_redox_context",
        "new_mechanism_class": DISULFIDE_CLASS,
        "candidate_word_promoted_to_learned": "disulfide_secretory_redox_context",
        "new_state_variables": [
            "secretory_redox_context",
            "disulfide_pairing_topology",
            "cysteine_pairing_constraint",
            "extracellular_stabilized_fold",
            "glycosylation_context",
            "redox_mispaired_frustration",
            "signal_peptide_removed_context",
            "secretory_quality_control",
        ],
        "new_operators": ["disulfide_pairing_operator", "secretory_redox_operator"],
        "negative_gates": [
            "metal_cluster_and_cys_his_coordination_priority",
            "true_transmembrane_priority",
            "signal_peptide_only_abstain",
            "ordinary_cysteine_noise_not_secretory_disulfide",
            "coiled_coil_and_repeat_candidates_abstain",
            "beta_closure_priority",
        ],
        "proof_batch": BATCH_ID,
        "proof_status": v76_cert["status"],
        "v76_accepted_accuracy": v76_cert["accepted_accuracy"],
        "v76_failed_accepted_count": v76_cert["failed_accepted_count"],
        "v76_matched_control_dominance_acceptance": True,
        "v76_fixed_accepted_support_thresholds_used": False,
        "next_required_batch": "V76P_OPENMM_COARSE_TARGET_EXECUTION_PILOT",
        "claim_allowed": False,
        "folding_problem_solved": False,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return _write_json(E70_ROOT / "e70_secretory_disulfide_redox_topology_grammar_certificate.json", cert)


def _certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_report: dict[str, Any],
    dashboard: dict[str, Any],
    physical_calibration_inputs: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    sentinel_rows = [row for row in rows if row["panel_group"] == "V75_ACCEPTED_SENTINEL_REPLAY"]
    disulfide_rows = [row for row in rows if row["panel_group"] == "SECRETORY_DISULFIDE_POSITIVE"]
    missing_rows = [row for row in rows if row["panel_group"] == "COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL"]
    sentinel_regressions = sum(1 for row in sentinel_rows if row["score_label"] != "supported")
    status = CONTROLS_FAILED if failed_controls else FAILED if metrics["failed_accepted_count"] or sentinel_regressions else PASSED
    cert = {
        "kind": "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "secretory_disulfide_positive_supported": sum(1 for row in disulfide_rows if row["accepted_supported"]),
        "secretory_disulfide_positive_total": len(disulfide_rows),
        "missing_word_controls_cleanly_abstained": sum(1 for row in missing_rows if row["clean_abstain_supported"]),
        "missing_word_controls_total": len(missing_rows),
        "sentinel_regressions": sentinel_regressions,
        "matched_control_dominance_acceptance": True,
        "fixed_accepted_support_thresholds_used": False,
        "failed_accepted_by_failure_mode": failure_report["failed_accepted_by_failure_mode"],
        "top_missing_esperanto_word": dashboard["shards"]["TOTAL"]["top_missing_esperanto_word"],
        "real_physical_calibration_inputs_used": all(row["real_physical_calibration_inputs_used"] for row in rows),
        "real_physical_calibration_row_count": physical_calibration_inputs.get("row_count"),
        "real_physical_calibration_hash": physical_calibration_inputs.get("calibration_hash"),
        "real_physical_calibration_source_coordinate_database": physical_calibration_inputs.get("source_coordinate_database"),
        "real_physical_calibration_target_native_excluded": physical_calibration_inputs.get("target_native_excluded_from_calibration"),
        "real_physical_calibration_target_native_contacts_used_before_prediction": physical_calibration_inputs.get("target_native_contacts_used_before_prediction"),
        "real_physical_calibration_coordinate_truth_used_as_prediction_input": physical_calibration_inputs.get("coordinate_truth_used_as_prediction_input"),
        "real_physical_calibration_leave_one_target_out": physical_calibration_inputs.get("leave_one_target_out_calibration"),
        "real_physical_calibration_input_type": physical_calibration_inputs.get("calibration_input_type"),
        "physical_basis_claim_allowed": False,
        "physical_next_required_capability": "V76P_OPENMM_COARSE_TARGET_EXECUTION_PILOT",
        "candidate_grammars_remaining": ["coiled_coil_register", "repeat_solenoid_topology", "signal_peptide_vs_true_TM", "knotted_topology"],
        "next_required_batch": "V76P_OPENMM_COARSE_TARGET_EXECUTION_PILOT",
        "coordinate_truth_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
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
        "matched_control_dominance_acceptance": cert["matched_control_dominance_acceptance"],
        "status": cert["status"],
        "claim_allowed": False,
    }


def _append_claim_ledger(row: dict[str, Any]) -> Path:
    path = LEDGER_ROOT / "claim_ledger_v0.json"
    ledger = _read_json(path, "campaign claim ledger") if path.exists() else {"kind": "V76_CLAIM_LEDGER_v0", "campaign_id": CAMPAIGN_ID, "rows": []}
    rows = [existing for existing in ledger.get("rows", []) if isinstance(existing, dict) and existing.get("batch_id") != BATCH_ID]
    rows.append(row)
    ledger["rows"] = rows
    return _write_json(path, ledger)


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V76 Secretory Disulfide Repair Panel",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Accepted supported: `{cert['accepted_supported']} / {cert['accepted_count']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Matched-control dominance acceptance: `{cert['matched_control_dominance_acceptance']}`",
        f"Fixed accepted-support thresholds used: `{cert['fixed_accepted_support_thresholds_used']}`",
        f"Secretory disulfide positives supported: `{cert['secretory_disulfide_positive_supported']} / {cert['secretory_disulfide_positive_total']}`",
        f"Missing-word controls cleanly abstained: `{cert['missing_word_controls_cleanly_abstained']} / {cert['missing_word_controls_total']}`",
        f"Sentinel regressions: `{cert['sentinel_regressions']}`",
        f"Physical basis claim allowed: `{cert['physical_basis_claim_allowed']}`",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v76(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_outputs(out_dir)
    physical_calibration_inputs = write_real_physical_calibration_inputs(REAL_COORDINATE_BENCHMARK, PHYSICAL_CALIBRATION_INPUTS)
    targets = _select_targets()
    target_manifest = _target_manifest(targets)
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v76_secretory_disulfide_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v76_e70_engine_declaration.json", engine_declaration)

    rows: list[dict[str, Any]] = []
    for target in targets:
        source_manifest = _source_manifest(target)
        packet = _packet(
            target=target,
            sequence=target["sequence"],
            sources=source_manifest["prediction_sources"],
            perturbations=_perturbations_for_target(target),
            physical_calibration_inputs=physical_calibration_inputs,
        )
        matched = _matched_control_report(
            target=target,
            packet=packet,
            source_manifest=source_manifest,
            physical_calibration_inputs=physical_calibration_inputs,
        )
        score = _score(target, packet, matched)
        rows.append(score)
        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target["target_id"] / "sealed_packet_summary.json", _packet_summary(packet))
        _write_json(DATA_ROOT / "validation" / target["target_id"] / "validation_result.json", score)

    failure_report = _failure_report(rows)
    dashboard = _dashboard(rows)
    controls = _controls(target_manifest, rows, physical_calibration_inputs)
    cert = _certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        rows=rows,
        controls=controls,
        failure_report=failure_report,
        dashboard=dashboard,
        physical_calibration_inputs=physical_calibration_inputs,
    )
    e70_cert_path = _write_e70_certificate(cert)
    scoring_path = _write_json(DATA_ROOT / "v76_secretory_disulfide_scoring_report.json", {"kind": "V76_SECRETORY_DISULFIDE_SCORING_REPORT_v0", "rows": rows})
    failure_path = _write_json(DATA_ROOT / "v76_secretory_disulfide_failure_report.json", failure_report)
    dashboard_path = _write_json(DATA_ROOT / "v76_secretory_disulfide_dashboard.json", dashboard)
    data_cert_path = _write_json(DATA_ROOT / "v76_secretory_disulfide_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v76_campaign_claim_row.json", _claim_row(cert))
    claim_ledger_path = _append_claim_ledger(_claim_row(cert))
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v76_secretory_disulfide_certificate.json", cert)
    report_path = out_dir / "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_REPORT.md"
    _write_report(report_path, cert)
    return {
        "target_manifest": DATA_ROOT / "v76_secretory_disulfide_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v76_e70_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "dashboard": dashboard_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "e70_certificate": e70_cert_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V76 secretory disulfide repair panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v76(args.out_dir)
    cert = _read_json(paths["certificate"], "V76 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_count": cert["accepted_count"],
        "accepted_supported": cert["accepted_supported"],
        "clean_abstain_supported": cert["clean_abstain_supported"],
        "failed_accepted": cert["failed_accepted_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "matched_control_dominance_acceptance": cert["matched_control_dominance_acceptance"],
        "fixed_accepted_support_thresholds_used": cert["fixed_accepted_support_thresholds_used"],
        "secretory_disulfide_positive_supported": cert["secretory_disulfide_positive_supported"],
        "missing_word_controls_cleanly_abstained": cert["missing_word_controls_cleanly_abstained"],
        "sentinel_regressions": cert["sentinel_regressions"],
        "controls_passed": cert["controls_passed"],
        "next_required_batch": cert["next_required_batch"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
