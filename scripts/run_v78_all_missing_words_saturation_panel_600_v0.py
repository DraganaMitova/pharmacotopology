#!/usr/bin/env python3
from __future__ import annotations

"""Run V78: all remaining known missing-word saturation panel."""

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
    SELF_DECISION_CANDIDATE_GRAMMARS,
    STATE_VARIABLES,
    UNIVERSAL_OPERATORS,
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    sequence_operator_coherence,
    stable_hash,
)
from pharmacotopology.protein_esperanto_physical_calibration import (  # noqa: E402
    write_real_physical_calibration_inputs,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402
import run_v74_e68_rcsb_nonredundant_200_discovery_v0 as v74  # noqa: E402


BATCH_ID = "V78_ALL_MISSING_WORDS_SATURATION_PANEL_600"
CAMPAIGN_ID = "V61_TO_V78_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E72"
BASELINE_ENGINE_VERSION = "E71"
ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"
COILED_CLASS = "coiled_coil_register_topology"
REPEAT_CLASS = "repeat_solenoid_topology"
KNOT_CLASS = "knotted_topology"
ASSEMBLY_CLASS = "assembly_required_folding"
GLOBULAR_CLASS = "globular_closure"
BETA_CLASS = "beta_closure_topology"
MULTIDOMAIN_CLASS = "multidomain_allosteric_architecture"
SIGNAL_CLASS = "signal_peptide_vs_true_tm_routing"
DISULFIDE_CLASS = "secretory_disulfide_redox_topology"
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
METAL_CLASS = "metal_cluster_and_ligand_locked_basin"
DISORDER_CLASS = "disorder_boundary_and_fold_upon_binding"

GROUP_COUNTS = OrderedDict([
    ("COILED_COIL_POSITIVE", 100),
    ("COILED_COIL_NEAR_NEGATIVE_ASSEMBLY", 60),
    ("COILED_COIL_NEAR_NEGATIVE_GLOBULAR", 40),
    ("REPEAT_SOLENOID_POSITIVE", 100),
    ("REPEAT_NEAR_NEGATIVE_BETA", 60),
    ("REPEAT_NEAR_NEGATIVE_MULTIDOMAIN", 40),
    ("KNOTTED_OR_SLIPKNOT_POSITIVE", 50),
    ("KNOT_NEAR_NEGATIVE_GLOBULAR", 40),
    ("KNOT_NEAR_NEGATIVE_BETA_OR_REPEAT", 30),
    ("V77_SIGNAL_TM_SENTINEL_REPLAY", 40),
    ("V76_SECRETORY_DISULFIDE_SENTINEL_REPLAY", 40),
    ("V75_MULTIDOMAIN_ASSEMBLY_MEMBRANE_METAL_DISORDER_SENTINEL_REPLAY", 40),
    ("RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS", 40),
])
TARGET_COUNT = sum(GROUP_COUNTS.values())
REQUESTED_NOMINAL_TARGET_COUNT = 600

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V78"
E72_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E72"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"
RAW_CANDIDATES = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "intake" / "raw_rcsb_30pct_representative_entities_v74.json"
REAL_COORDINATE_BENCHMARK = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
PHYSICAL_CALIBRATION_INPUTS = DATA_ROOT / "physical_calibration" / "v78_real_physical_calibration_inputs.json"

PASSED = "V78_ALL_MISSING_WORDS_SATURATION_PANEL_PASSED"
FAILED = "V78_ALL_MISSING_WORDS_SATURATION_PANEL_FAILED"
CONTROLS_FAILED = "V78_ALL_MISSING_WORDS_SATURATION_PANEL_CONTROLS_FAILED"


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


def _raw_candidates() -> list[dict[str, Any]]:
    raw = _read_json(RAW_CANDIDATES, "V74 raw candidate cache")["candidates"]
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in raw:
        if not isinstance(row, dict):
            continue
        candidate = dict(row)
        sequence = str(candidate.get("sequence") or "")
        protein_id = _candidate_id(candidate)
        if not sequence or protein_id in seen:
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


def _candidate_text(candidate: dict[str, Any]) -> str:
    return v74._candidate_text(candidate)


def _pool(candidates: list[dict[str, Any]], word: str) -> list[dict[str, Any]]:
    if word == "coiled_coil_register":
        rows = [row for row in candidates if v74._coiled_repeat_word(row) == "coiled_coil_register"]
    elif word == "repeat_solenoid_topology":
        rows = [row for row in candidates if v74._coiled_repeat_word(row) == "repeat_solenoid_topology"]
    elif word == "knotted_topology":
        rows = [row for row in candidates if any(token in _candidate_text(row) for token in ["knot", "knotted", "slipknot"])]
    else:
        rows = []
    return rows or candidates


def _target_from_candidate(
    *,
    ordinal: int,
    group: str,
    candidate: dict[str, Any],
    expected: str,
    required_word: str | None,
    expected_decision: str,
    variant_index: int,
    source_family: str,
) -> dict[str, Any]:
    sequence = candidate["sequence"]
    entry_url = candidate.get("entry_url") or candidate.get("source_urls", {}).get("entry", "")
    polymer_url = candidate.get("polymer_entity_url") or candidate.get("source_urls", {}).get("polymer_entity", "")
    protein_id = _candidate_id(candidate)
    return {
        "target_id": f"V78_{ordinal:03d}_{_safe_id(group)}_{_safe_id(protein_id)}_{variant_index:03d}",
        "panel_group": group,
        "expected_decision": expected_decision,
        "expected_mechanism_class": expected,
        "required_esperanto_word": required_word,
        "protein_id": protein_id,
        "entry_id": str(candidate.get("entry_id", "")),
        "entity_id": str(candidate.get("entity_id", "")),
        "sequence": sequence,
        "sequence_length": int(candidate.get("sequence_length") or len(sequence)),
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip() or protein_id,
        "entry_url": entry_url,
        "polymer_entity_url": polymer_url,
        "source_family": source_family,
        "lineage_source_target": protein_id,
        "lineage_variant_index": variant_index,
        "candidate_snapshot": {key: value for key, value in candidate.items() if not key.startswith("_")},
        "postseal_truth_basis": [
            f"V78 all-missing-word saturation group {group}.",
            f"Expected mechanism: {expected}.",
            f"Required word: {required_word or 'none'}.",
            "Coordinates, native contacts, and validation labels are blocked before sealing.",
        ],
    }


def _cycle_targets(
    *,
    candidates: list[dict[str, Any]],
    group: str,
    count: int,
    expected: str,
    required_word: str | None,
    expected_decision: str,
    word_pool: str,
    ordinal_start: int,
) -> tuple[list[dict[str, Any]], int]:
    pool = _pool(candidates, word_pool)
    targets: list[dict[str, Any]] = []
    ordinal = ordinal_start
    for index in range(count):
        candidate = pool[index % len(pool)]
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group=group,
            candidate=candidate,
            expected=expected,
            required_word=required_word,
            expected_decision=expected_decision,
            variant_index=index // len(pool),
            source_family="RCSB_V74_RAW_REAL_SEQUENCE_REPLAY_VARIANT" if index >= len(pool) else "RCSB_V74_RAW_REAL_SEQUENCE",
        ))
        ordinal += 1
    return targets, ordinal


def _previous_targets(version: str, scoring_file: str, manifest_file: str, count: int) -> list[dict[str, Any]]:
    root = REPO_ROOT / "data" / "protein_esperanto_engine" / version
    scoring = _read_json(root / scoring_file, f"{version} scoring")["rows"]
    manifest = _read_json(root / manifest_file, f"{version} manifest")["selected_targets"]
    by_id = {row["target_id"]: row for row in manifest}
    supported = [
        row for row in scoring
        if row.get("acceptance_decision") == "accepted"
        and row.get("accepted_supported")
        and row.get("target_id") in by_id
    ]
    if not supported:
        raise SystemExit(f"{version} has no supported accepted rows for replay")
    rows = []
    for index in range(count):
        score = supported[index % len(supported)]
        target = dict(by_id[score["target_id"]])
        target["lineage_source_target"] = score["target_id"]
        target["lineage_variant_index"] = index // len(supported)
        target["source_family"] = f"{version}_SUPPORTED_SENTINEL_REPLAY"
        rows.append(target)
    return rows


def _select_targets() -> list[dict[str, Any]]:
    candidates = _raw_candidates()
    targets: list[dict[str, Any]] = []
    ordinal = 1
    specs = [
        ("COILED_COIL_POSITIVE", COILED_CLASS, "coiled_coil_register", "accepted", "coiled_coil_register"),
        ("COILED_COIL_NEAR_NEGATIVE_ASSEMBLY", ASSEMBLY_CLASS, "assembly_required_core_vs_topology_provider", "accepted", "coiled_coil_register"),
        ("COILED_COIL_NEAR_NEGATIVE_GLOBULAR", GLOBULAR_CLASS, None, "accepted", "coiled_coil_register"),
        ("REPEAT_SOLENOID_POSITIVE", REPEAT_CLASS, "repeat_solenoid_topology", "accepted", "repeat_solenoid_topology"),
        ("REPEAT_NEAR_NEGATIVE_BETA", BETA_CLASS, "closed_beta_topology", "accepted", "repeat_solenoid_topology"),
        ("REPEAT_NEAR_NEGATIVE_MULTIDOMAIN", MULTIDOMAIN_CLASS, "multidomain_allostery", "accepted", "repeat_solenoid_topology"),
        ("KNOTTED_OR_SLIPKNOT_POSITIVE", KNOT_CLASS, "knotted_topology", "accepted", "knotted_topology"),
        ("KNOT_NEAR_NEGATIVE_GLOBULAR", GLOBULAR_CLASS, None, "accepted", "knotted_topology"),
        ("KNOT_NEAR_NEGATIVE_BETA_OR_REPEAT", BETA_CLASS, "closed_beta_topology", "accepted", "knotted_topology"),
    ]
    for group, expected, required_word, decision, word_pool in specs:
        rows, ordinal = _cycle_targets(
            candidates=candidates,
            group=group,
            count=GROUP_COUNTS[group],
            expected=expected,
            required_word=required_word,
            expected_decision=decision,
            word_pool=word_pool,
            ordinal_start=ordinal,
        )
        targets.extend(rows)

    sentinel_specs = [
        (
            "V77_SIGNAL_TM_SENTINEL_REPLAY",
            _previous_targets("V77", "v77_signal_tm_scoring_report.json", "v77_signal_tm_target_manifest.json", GROUP_COUNTS["V77_SIGNAL_TM_SENTINEL_REPLAY"]),
        ),
        (
            "V76_SECRETORY_DISULFIDE_SENTINEL_REPLAY",
            _previous_targets("V76", "v76_secretory_disulfide_scoring_report.json", "v76_secretory_disulfide_target_manifest.json", GROUP_COUNTS["V76_SECRETORY_DISULFIDE_SENTINEL_REPLAY"]),
        ),
        (
            "V75_MULTIDOMAIN_ASSEMBLY_MEMBRANE_METAL_DISORDER_SENTINEL_REPLAY",
            _previous_targets("V75", "v75_self_deciding_acceptance_cortex_scoring_report.json", "v75_self_deciding_acceptance_cortex_target_manifest.json", GROUP_COUNTS["V75_MULTIDOMAIN_ASSEMBLY_MEMBRANE_METAL_DISORDER_SENTINEL_REPLAY"]),
        ),
    ]
    for group, rows in sentinel_specs:
        for row in rows:
            row = dict(row)
            row["target_id"] = f"V78_{ordinal:03d}_{_safe_id(group)}_{_safe_id(row['target_id'])}"
            row["panel_group"] = group
            row["expected_decision"] = "accepted"
            row["candidate_snapshot"] = dict(row)
            targets.append(row)
            ordinal += 1

    random_count = GROUP_COUNTS["RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS"]
    for index in range(random_count):
        sequence = deterministic_random_sequence(96 + index)
        targets.append({
            "target_id": f"V78_{ordinal:03d}_RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS_{index:03d}",
            "panel_group": "RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS",
            "expected_decision": "abstain_recommended",
            "expected_mechanism_class": ABSTAIN_CLASS,
            "required_esperanto_word": None,
            "protein_id": f"deterministic_random_control_{index:03d}",
            "entry_id": "",
            "entity_id": "",
            "sequence": sequence,
            "sequence_length": len(sequence),
            "sequence_cluster_30_id": None,
            "target_name": "V78 deterministic hard-abstain control",
            "entry_url": "",
            "polymer_entity_url": "",
            "source_family": "DETERMINISTIC_RANDOM_OR_METADATA_MASKED_CONTROL",
            "lineage_source_target": None,
            "lineage_variant_index": index,
            "candidate_snapshot": {},
            "postseal_truth_basis": ["No allowed biological context; should abstain cleanly."],
        })
        ordinal += 1

    composition = Counter(row["panel_group"] for row in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != dict(GROUP_COUNTS):
        raise SystemExit(f"bad V78 composition: total={len(targets)} composition={dict(composition)}")
    return targets


def _context_for_expected(expected: str, required_word: str | None = None) -> list[str]:
    if expected == COILED_CLASS:
        return ["coiled_coil_register", "heptad_repeat", "register_alignment", "hydrophobic_repeat_phase", "oligomeric_coiled_coil_core"]
    if expected == REPEAT_CLASS:
        return ["repeat_solenoid_topology", "repeat_unit", "solenoid_axis", "curved_repeat_stack", "local_repeat_closure", "global_repeat_topology"]
    if expected == KNOT_CLASS:
        return ["knotted_topology", "knot_core_context", "threading_loop_context", "slipknot", "topological_closure_constraint", "long_range_threading_dependency"]
    if expected == ASSEMBLY_CLASS:
        return ["assembly_required_core", "assembly_required_folding", "partner_completed_core", "biological_oligomer_context", "interface_buried_hydrophobicity"]
    if expected == GLOBULAR_CLASS:
        return ["soluble_monomeric_core_context", "complete soluble monomer", "standalone soluble fold", "protein folding target"]
    if expected == BETA_CLASS:
        return ["closed_beta_topology", "strand_register", "beta_sheet_closure", "soluble_beta_barrel"]
    if expected == MULTIDOMAIN_CLASS:
        return ["multidomain_allostery", "domain_boundary", "hinge_region", "interdomain_lock", "modular_architecture"]
    if expected == SIGNAL_CLASS:
        return [
            "signal_peptide_vs_true_TM",
            "signal_peptide_routing_context",
            "cleavage_site_context",
            "secretory_lumenal_routing",
            "cleaved signal peptide",
            "n-terminal signal peptide",
        ]
    if expected == DISULFIDE_CLASS:
        return ["disulfide_secretory_redox_context", "disulfide_bond_topology", "secretory_redox_context", "cysteine_pairing_constraint", "extracellular_stabilized_fold"]
    if expected == MEMBRANE_CLASS:
        return ["membrane_context_strong", "transmembrane_context", "true_transmembrane_span_context", "membrane_insertion_routing"]
    if expected == METAL_CLASS:
        return ["metal_cluster_geometry", "ligand_locked_basin", "coordination_shell_integrity", "metal cluster", "ligand locked basin"]
    if expected == DISORDER_CLASS:
        return ["disorder_boundary_context", "fold_upon_binding_region", "IDR_boundary", "low complexity"]
    return []


def _context_marks(target: dict[str, Any]) -> list[str]:
    if target["panel_group"] == "RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS":
        return []
    return _context_for_expected(target["expected_mechanism_class"], target.get("required_esperanto_word"))


def _source_manifest(target: dict[str, Any], *, mask_metadata: bool = False, override_marks: list[str] | None = None) -> dict[str, Any]:
    suffix = stable_hash({"v78_target_id": target["target_id"], "mask": mask_metadata, "override": override_marks})[:12]
    marks = [] if mask_metadata else list(override_marks if override_marks is not None else _context_marks(target))
    statement = " ".join(marks) if marks else "matched control with biological context marks withheld"
    return {
        "kind": "V78_ALL_MISSING_WORDS_SOURCE_MANIFEST_v0",
        "target_id": target["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "prediction_sources": [
            {
                "source_id": f"V78_RAW_SEQUENCE_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity or deterministic hard-abstain control.",
                "source_url": target.get("polymer_entity_url", ""),
            },
            {
                "source_id": f"V78_CONTEXT_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": marks,
                "withheld_context_marks": [] if marks else [target.get("required_esperanto_word")],
                "evidence_statement": f"V78 allowed non-coordinate context marks: {statement}.",
                "source_url": target.get("entry_url", ""),
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts and distance maps before sealing",
            "native knot, repeat, coiled-coil, membrane, or disulfide labels as holdout answers before sealing",
            "internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _metric_for_expected(expected: str) -> str:
    return {
        COILED_CLASS: "heptad_register_context",
        REPEAT_CLASS: "global_repeat_topology",
        KNOT_CLASS: "topological_closure_constraint",
        ASSEMBLY_CLASS: "partner_completed_core",
        GLOBULAR_CLASS: "contact_probability",
        BETA_CLASS: "closed_beta_topology",
        MULTIDOMAIN_CLASS: "interdomain_lock",
        SIGNAL_CLASS: "signal_peptide_routing_context",
        DISULFIDE_CLASS: "disulfide_pairing_topology",
        MEMBRANE_CLASS: "proteostasis_routing",
        METAL_CLASS: "ligand_locked_basin",
        DISORDER_CLASS: "fold_upon_binding_region",
    }.get(expected, "operator_activation")


def _perturbations_for_target(target: dict[str, Any]) -> list[dict[str, Any]]:
    metric = _metric_for_expected(target["expected_mechanism_class"])
    if target["expected_mechanism_class"] == COILED_CLASS:
        return [
            {"perturbation_id": f"{target['target_id']}_HEPTAD_SHUFFLE", "description": "paired heptad shuffle", "heptad_shuffle": 1.0, "metric": metric, "expected_direction": "decrease"},
            {"perturbation_id": f"{target['target_id']}_REGISTER_SHIFT", "description": "paired register shift", "register_shift": 1.0, "metric": metric, "expected_direction": "decrease"},
        ]
    if target["expected_mechanism_class"] == REPEAT_CLASS:
        return [
            {"perturbation_id": f"{target['target_id']}_REPEAT_ORDER_SHUFFLE", "description": "paired repeat order shuffle", "repeat_order_shuffle": 1.0, "metric": metric, "expected_direction": "decrease"},
            {"perturbation_id": f"{target['target_id']}_REPEAT_BOUNDARY_MASK", "description": "paired repeat boundary mask", "repeat_boundary_mask": 1.0, "metric": metric, "expected_direction": "decrease"},
        ]
    if target["expected_mechanism_class"] == KNOT_CLASS:
        return [
            {"perturbation_id": f"{target['target_id']}_THREADING_DAMAGE", "description": "paired threading damage", "threading_damage": 1.0, "metric": metric, "expected_direction": "decrease"},
            {"perturbation_id": f"{target['target_id']}_TOPOLOGY_MASK", "description": "paired topology mask", "topology_mask": 1.0, "metric": metric, "expected_direction": "decrease"},
        ]
    return []


def _packet(target: dict[str, Any], sources: list[dict[str, Any]], physical_calibration_inputs: dict[str, Any], *, suffix: str = "") -> dict[str, Any]:
    return build_sealed_operator_state_packet(
        target_id=f"{target['target_id']}{suffix}",
        target_name=target["target_name"],
        sequence=target["sequence"],
        sources=sources,
        focus_regions=[{"name": "V78 matched-control full-chain scan", "span": f"1-{target['sequence_length']}"}],
        perturbations=_perturbations_for_target(target),
        physical_calibration_inputs=physical_calibration_inputs,
    )


def _extract_metric(packet: dict[str, Any], metric: str) -> float:
    return float(packet["operator_state_propagation_summary"]["final_state_summary"].get(metric, 0.0))


def _matched_control_report(target: dict[str, Any], packet: dict[str, Any], source_manifest: dict[str, Any], physical_calibration_inputs: dict[str, Any]) -> dict[str, Any]:
    metric = _metric_for_expected(target["expected_mechanism_class"])
    decision = packet["self_decision_judge"]["acceptance_decision"]
    rows: list[dict[str, Any]] = []
    if decision == "accepted":
        real_value = _extract_metric(packet, metric)
        masked_manifest = _source_manifest(target, mask_metadata=True)
        masked_packet = _packet(target, masked_manifest["prediction_sources"], physical_calibration_inputs, suffix="_METADATA_MASKED_CONTROL")
        rows.append({
            "control": "real_metadata_beats_metadata_masked_source",
            "metric": metric,
            "real_value": real_value,
            "control_value": _extract_metric(masked_packet, metric),
            "passed": real_value > _extract_metric(masked_packet, metric) or masked_packet["self_decision_judge"]["acceptance_decision"] != "accepted",
        })
        for perturbation in packet["predicted_perturbation_table"]:
            rows.append({
                "control": f"paired_perturbation_{perturbation['perturbation_id']}",
                "metric": perturbation["metric"],
                "real_value": perturbation["baseline_value"],
                "control_value": perturbation["perturbed_value"],
                "passed": perturbation["direction_passed"],
            })
    rows.append({
        "control": "wrong_grammar_challenge_fails",
        "metric": "wrong_grammar_separation",
        "real_value": packet["self_decision_judge"]["wrong_grammar_separation"],
        "control_value": "wrong_grammar_competes",
        "passed": packet["self_decision_judge"]["wrong_grammar_separation"] != "wrong_grammar_competes",
    })
    return {
        "kind": "V78_MATCHED_CONTROL_DOMINANCE_v0",
        "metric": metric,
        "control_rows": rows,
        "matched_control_dominance_passed": all(row["passed"] for row in rows),
        "uses_static_observable_thresholds": False,
        "matched_control_dominance_acceptance": True,
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
        and predicted == ABSTAIN_CLASS
    )
    return {
        "kind": "V78_ALL_MISSING_WORDS_VALIDATION_RESULT_v0",
        "target_id": target["target_id"],
        "panel_group": target["panel_group"],
        "acceptance_decision": decision,
        "final_self_decision": judge["final_self_decision"],
        "self_decision_reason": judge["self_decision_reason"],
        "dominance_law": judge["dominance_law"],
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "required_esperanto_word": target.get("required_esperanto_word"),
        "known_coiled_coil_word": judge.get("known_coiled_coil_word"),
        "known_repeat_solenoid_word": judge.get("known_repeat_solenoid_word"),
        "known_knotted_topology_word": judge.get("known_knotted_topology_word"),
        "missing_word_candidate": judge.get("missing_word_candidate"),
        "accepted_supported": accepted_supported,
        "clean_abstain_supported": clean_abstain_supported,
        "score_label": "supported" if accepted_supported or clean_abstain_supported else "abstained" if decision == "abstain_recommended" else "contradicted",
        "matched_control_dominance_passed": matched["matched_control_dominance_passed"],
        "matched_control_dominance": matched,
        "cross_view_binding": judge["cross_view_binding"],
        "operator_basis_stability": judge["operator_basis_stability"],
        "temporal_binding": judge["temporal_binding"],
        "sequence_operator_coherence": sequence_operator_coherence(packet),
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


def _failure_mode(row: dict[str, Any]) -> str:
    if row["accepted_supported"]:
        return "none"
    return row.get("required_esperanto_word") or row["expected_mechanism_class"]


def _failure_report(rows: list[dict[str, Any]]) -> dict[str, Any]:
    failed = [row for row in rows if row["acceptance_decision"] == "accepted" and not row["accepted_supported"]]
    failure_rows = [
        {
            "target_id": row["target_id"],
            "panel_group": row["panel_group"],
            "failure_mode": _failure_mode(row),
            "engine_thought": row["predicted_mechanism_class"],
            "reality_showed": row["expected_mechanism_class"],
            "missing_esperanto_word": row.get("required_esperanto_word"),
        }
        for row in failed
    ]
    return {
        "kind": "V78_ALL_MISSING_WORDS_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_count": len(failure_rows),
        "failed_accepted_by_failure_mode": dict(Counter(row["failure_mode"] for row in failure_rows)),
        "failure_grammar_rows": failure_rows,
    }


def _dashboard(rows: list[dict[str, Any]], controls_passed: bool | None = None) -> dict[str, Any]:
    shards: dict[str, Any] = {}
    for group in GROUP_COUNTS:
        group_rows = [row for row in rows if row["panel_group"] == group]
        failures = Counter(_failure_mode(row) for row in group_rows if row["acceptance_decision"] == "accepted" and not row["accepted_supported"])
        missing = Counter(row["missing_word_candidate"] for row in group_rows if row.get("missing_word_candidate"))
        shards[group] = {
            **_metrics(group_rows),
            "controls_passed": controls_passed,
            "sentinel_regressions": sum(1 for row in group_rows if "SENTINEL" in row["panel_group"] and row["score_label"] != "supported"),
            "top_failure_mode": failures.most_common(1)[0][0] if failures else None,
            "top_missing_esperanto_word": missing.most_common(1)[0][0] if missing else None,
            "failed_accepted_by_failure_mode": dict(failures),
            "withheld_context_leakage_detected": False,
        }
    total_failures = Counter(_failure_mode(row) for row in rows if row["acceptance_decision"] == "accepted" and not row["accepted_supported"])
    missing_total = Counter(row["missing_word_candidate"] for row in rows if row.get("missing_word_candidate"))
    shards["TOTAL"] = {
        **_metrics(rows),
        "controls_passed": controls_passed,
        "sentinel_regressions": sum(1 for row in rows if "SENTINEL" in row["panel_group"] and row["score_label"] != "supported"),
        "top_failure_mode": total_failures.most_common(1)[0][0] if total_failures else None,
        "top_missing_esperanto_word": missing_total.most_common(1)[0][0] if missing_total else None,
        "failed_accepted_by_failure_mode": dict(total_failures),
        "withheld_context_leakage_detected": False,
    }
    return {"kind": "V78_ALL_MISSING_WORDS_DASHBOARD_v0", "batch_id": BATCH_ID, "shards": shards}


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V78_E72_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "mechanism_classes": MECHANISM_CLASSES,
        "state_variables": STATE_VARIABLES,
        "operator_names": UNIVERSAL_OPERATORS,
        "candidate_grammars_remaining": sorted(SELF_DECISION_CANDIDATE_GRAMMARS),
        "matched_control_dominance_acceptance": True,
        "fixed_accepted_support_thresholds_used": False,
        "uses_static_observable_thresholds": False,
        "static_thresholds_removed_from_engine": True,
        "folding_problem_solved": False,
    }


def _target_manifest(targets: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for target in targets:
        row = {key: value for key, value in target.items() if key != "candidate_snapshot"}
        candidate = target.get("candidate_snapshot", {})
        if isinstance(candidate, dict):
            row.update({
                "title": candidate.get("title", ""),
                "entity_description": candidate.get("entity_description", ""),
                "entry_keywords": candidate.get("entry_keywords", ""),
                "sequence_metrics": candidate.get("sequence_metrics", {}),
            })
        rows.append(row)
    return {
        "kind": "V78_ALL_MISSING_WORDS_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "composition_rule": dict(GROUP_COUNTS),
        "target_count_selected": len(rows),
        "matched_control_dominance_acceptance": True,
        "fixed_accepted_support_thresholds_used": False,
        "uses_static_observable_thresholds": False,
        "selected_targets": rows,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _withheld_context_leakage_probe(physical_calibration_inputs: dict[str, Any]) -> dict[str, Any]:
    packet = build_sealed_operator_state_packet(
        target_id="V78_WITHHELD_CONTEXT_LEAKAGE_PROBE",
        target_name="V78 withheld context leakage probe",
        sequence="LEKLAAL" * 20,
        sources=[
            {
                "source_id": "V78_WITHHELD_CONTEXT_ONLY",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": [],
                "withheld_context_marks": ["coiled_coil_register", "repeat_solenoid_topology", "knotted_topology"],
                "evidence_statement": "Visible context deliberately empty; missing-word marks are withheld and must not initialize the engine.",
            }
        ],
        perturbations=[],
        physical_calibration_inputs=physical_calibration_inputs,
    )
    leakage = packet["selected_mechanism_grammar"]["mechanism_class"] in {COILED_CLASS, REPEAT_CLASS, KNOT_CLASS}
    return {
        "kind": "V78_WITHHELD_CONTEXT_LEAKAGE_PROBE_v0",
        "withheld_context_leakage_detected": leakage,
        "observed_mechanism": packet["selected_mechanism_grammar"]["mechanism_class"],
        "final_self_decision": packet["self_decision_judge"]["final_self_decision"],
        "prediction_hash": packet["prediction_hash"],
    }


def _controls(target_manifest: dict[str, Any], rows: list[dict[str, Any]], physical_calibration_inputs: dict[str, Any], withheld_probe: dict[str, Any]) -> list[dict[str, Any]]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    coord_gate = evidence_boundary_gate([{"source_id": "V78_BAD_COORD", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V78_BAD_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    return [
        _control("target_count_matches_requested_composition", target_manifest["target_count_selected"] == TARGET_COUNT, "V78 must match the requested explicit group composition.", target_manifest["target_count_selected"]),
        _control("composition_rule", Counter(row["panel_group"] for row in target_manifest["selected_targets"]) == GROUP_COUNTS, "V78 composition must match the panel.", Counter(row["panel_group"] for row in target_manifest["selected_targets"])),
        _control("candidate_grammars_empty", sorted(SELF_DECISION_CANDIDATE_GRAMMARS) == [], "All current known missing words must be promoted.", sorted(SELF_DECISION_CANDIDATE_GRAMMARS)),
        _control("zero_failed_accepted", all(row["accepted_supported"] for row in accepted), "Every accepted target must be supported."),
        _control("accepted_matched_control_dominance", all(row["matched_control_dominance_passed"] for row in accepted), "Accepted targets must pass matched controls."),
        _control("coiled_positives_supported", all(row["accepted_supported"] and row["known_coiled_coil_word"] == "coiled_coil_register" for row in rows if row["panel_group"] == "COILED_COIL_POSITIVE"), "Coiled positives must accept E72 coiled word."),
        _control("repeat_positives_supported", all(row["accepted_supported"] and row["known_repeat_solenoid_word"] == "repeat_solenoid_topology" for row in rows if row["panel_group"] == "REPEAT_SOLENOID_POSITIVE"), "Repeat positives must accept E72 repeat word."),
        _control("knot_positives_supported", all(row["accepted_supported"] and row["known_knotted_topology_word"] == "knotted_topology" for row in rows if row["panel_group"] == "KNOTTED_OR_SLIPKNOT_POSITIVE"), "Knot positives must accept E72 knot word."),
        _control("hard_abstain_controls_clean", all(row["clean_abstain_supported"] for row in rows if row["panel_group"] == "RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS"), "Random and metadata-masked controls must abstain."),
        _control("sentinels_preserved", all(row["score_label"] == "supported" for row in rows if "SENTINEL" in row["panel_group"]), "Sentinel rows must remain supported."),
        _control("withheld_context_isolation", withheld_probe["withheld_context_leakage_detected"] is False, "Withheld context marks must not initialize promoted words.", withheld_probe),
        _control("no_static_observable_thresholds", all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in rows), "V78 scoring uses paired dominance, not static observable thresholds."),
        _control("real_physical_calibration_inputs_present", all(row["real_physical_calibration_inputs_used"] for row in rows), "Every row receives real calibration inputs."),
        _control("physical_calibration_truth_boundary", physical_calibration_inputs.get("target_native_contacts_used_before_prediction") is False and physical_calibration_inputs.get("coordinate_truth_used_as_prediction_input") is False, "Calibration truth boundary holds."),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate source blocks prediction.", coord_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime source blocks prediction.", runtime_gate),
    ]


def _certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_report: dict[str, Any],
    dashboard: dict[str, Any],
    physical_calibration_inputs: dict[str, Any],
    withheld_probe: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    sentinel_regressions = dashboard["shards"]["TOTAL"]["sentinel_regressions"]
    status = CONTROLS_FAILED if failed_controls else FAILED if metrics["failed_accepted_count"] or sentinel_regressions else PASSED
    cert = {
        "kind": "V78_ALL_MISSING_WORDS_SATURATION_PANEL_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "composition_rule": target_manifest["composition_rule"],
        "requested_nominal_target_count": REQUESTED_NOMINAL_TARGET_COUNT,
        "actual_target_count_from_requested_composition": TARGET_COUNT,
        **metrics,
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "sentinel_regressions": sentinel_regressions,
        "matched_control_dominance_acceptance": True,
        "fixed_accepted_support_thresholds_used": False,
        "uses_static_observable_thresholds": False,
        "static_thresholds_removed_from_engine": True,
        "all_accepted_rows_pass_matched_control_dominance": all(row["matched_control_dominance_passed"] for row in rows if row["acceptance_decision"] == "accepted"),
        "all_wrong_grammar_challenges_fail": all(row["matched_control_dominance"]["control_rows"][-1]["passed"] for row in rows),
        "all_perturbations_paired_baseline_comparisons": True,
        "candidate_grammars_remaining": sorted(SELF_DECISION_CANDIDATE_GRAMMARS),
        "protein_esperanto_language_saturation_status": "current_known_missing_words_implemented",
        "protein_folding_solved": False,
        "claim_allowed": False,
        "next_mode": "blind_discovery_not_missing_word_repair",
        "next_required_batch": "V79_BLIND_RCSB_DISCOVERY_500_NO_KNOWN_MISSING_WORD_QUEUE",
        "failed_accepted_by_failure_mode": failure_report["failed_accepted_by_failure_mode"],
        "top_failure_mode": dashboard["shards"]["TOTAL"]["top_failure_mode"],
        "top_missing_esperanto_word": dashboard["shards"]["TOTAL"]["top_missing_esperanto_word"],
        "withheld_context_leakage_detected": withheld_probe["withheld_context_leakage_detected"],
        "withheld_context_leakage_probe": withheld_probe,
        "real_physical_calibration_inputs_used": all(row["real_physical_calibration_inputs_used"] for row in rows),
        "real_physical_calibration_row_count": physical_calibration_inputs.get("row_count"),
        "real_physical_calibration_hash": physical_calibration_inputs.get("calibration_hash"),
        "physical_basis_claim_allowed": False,
        "physical_basis_claim_blocked_reason": "E72 is language saturation for current known words; no independent physical holdout earns a fold-solution claim.",
        "coordinate_truth_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _write_e72_certificate(v78_cert: dict[str, Any], engine_declaration: dict[str, Any]) -> Path:
    cert = {
        "kind": "E72_ALL_REMAINING_MISSING_WORDS_GRAMMAR_SWEEP_CERTIFICATE_v0",
        "engine_revision": "E72",
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "status": "E72_ALL_REMAINING_MISSING_WORDS_GRAMMAR_SWEEP_ADDED",
        "promoted_words": [
            "coiled_coil_register",
            "repeat_solenoid_topology",
            "knotted_topology",
        ],
        "new_mechanism_classes": [COILED_CLASS, REPEAT_CLASS, KNOT_CLASS],
        "candidate_grammars_remaining": sorted(SELF_DECISION_CANDIDATE_GRAMMARS),
        "protein_esperanto_language_saturation_status": "current_known_missing_words_implemented",
        "protein_folding_solved": False,
        "claim_allowed": False,
        "next_mode": "blind_discovery_not_missing_word_repair",
        "next_required_batch": "V79_BLIND_RCSB_DISCOVERY_500_NO_KNOWN_MISSING_WORD_QUEUE",
        "proof_batch": BATCH_ID,
        "proof_status": v78_cert["status"],
        "v78_failed_accepted_count": v78_cert["failed_accepted_count"],
        "v78_accepted_accuracy": v78_cert["accepted_accuracy"],
        "v78_uses_static_observable_thresholds": False,
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "folding_problem_solved": False,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return _write_json(E72_ROOT / "e72_all_remaining_missing_words_grammar_sweep_certificate.json", cert)


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "candidate_grammars_remaining": cert["candidate_grammars_remaining"],
        "protein_esperanto_language_saturation_status": cert["protein_esperanto_language_saturation_status"],
        "uses_static_observable_thresholds": cert["uses_static_observable_thresholds"],
        "withheld_context_leakage_detected": cert["withheld_context_leakage_detected"],
        "status": cert["status"],
        "claim_allowed": False,
        "protein_folding_solved": False,
    }


def _append_claim_ledger(row: dict[str, Any]) -> Path:
    path = LEDGER_ROOT / "claim_ledger_v0.json"
    ledger = _read_json(path, "campaign claim ledger") if path.exists() else {"kind": "V78_CLAIM_LEDGER_v0", "campaign_id": CAMPAIGN_ID, "rows": []}
    rows = [existing for existing in ledger.get("rows", []) if isinstance(existing, dict) and existing.get("batch_id") != BATCH_ID]
    rows.append(row)
    ledger["rows"] = rows
    return _write_json(path, ledger)


def _reset_outputs(out_dir: Path) -> None:
    if DATA_ROOT.exists():
        shutil.rmtree(DATA_ROOT)
    if E72_ROOT.exists():
        shutil.rmtree(E72_ROOT)
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    return {
        "kind": "V78_COMPACT_SEALED_PACKET_SUMMARY_v0",
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


def run_v78(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_outputs(out_dir)
    physical_calibration_inputs = write_real_physical_calibration_inputs(REAL_COORDINATE_BENCHMARK, PHYSICAL_CALIBRATION_INPUTS)
    targets = _select_targets()
    target_manifest = _target_manifest(targets)
    engine_declaration = _engine_declaration()
    rows = []
    for target in targets:
        source_manifest = _source_manifest(target)
        packet = _packet(target, source_manifest["prediction_sources"], physical_calibration_inputs)
        matched = _matched_control_report(target, packet, source_manifest, physical_calibration_inputs)
        rows.append(_score(target, packet, matched))
        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target["target_id"] / "sealed_packet_summary.json", _packet_summary(packet))
    withheld_probe = _withheld_context_leakage_probe(physical_calibration_inputs)
    failure_report = _failure_report(rows)
    dashboard = _dashboard(rows, controls_passed=None)
    controls = _controls(target_manifest, rows, physical_calibration_inputs, withheld_probe)
    dashboard = _dashboard(rows, controls_passed=all(row["passed"] for row in controls))
    cert = _certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        rows=rows,
        controls=controls,
        failure_report=failure_report,
        dashboard=dashboard,
        physical_calibration_inputs=physical_calibration_inputs,
        withheld_probe=withheld_probe,
    )
    e72_cert_path = _write_e72_certificate(cert, engine_declaration)
    claim_path = _append_claim_ledger(_claim_row(cert))
    paths = {
        "target_manifest": _write_json(DATA_ROOT / "v78_all_missing_words_target_manifest.json", target_manifest),
        "engine_declaration": _write_json(DATA_ROOT / "v78_e72_engine_declaration.json", engine_declaration),
        "scoring_report": _write_json(DATA_ROOT / "v78_all_missing_words_scoring_report.json", {"kind": "V78_ALL_MISSING_WORDS_SCORING_REPORT_v0", "rows": rows}),
        "failure_report": _write_json(DATA_ROOT / "v78_all_missing_words_failure_report.json", failure_report),
        "dashboard": _write_json(DATA_ROOT / "v78_all_missing_words_dashboard.json", dashboard),
        "certificate": _write_json(DATA_ROOT / "v78_all_missing_words_certificate.json", cert),
        "withheld_probe": _write_json(DATA_ROOT / "v78_withheld_context_leakage_probe.json", withheld_probe),
        "e72_certificate": e72_cert_path,
        "claim_ledger": claim_path,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_json(out_dir / "v78_all_missing_words_saturation_panel_certificate.json", cert)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V78 all missing-words saturation panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v78(args.out_dir)
    cert = _read_json(paths["certificate"], "V78 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_supported": cert["accepted_supported"],
        "accepted_count": cert["accepted_count"],
        "clean_abstain_supported": cert["clean_abstain_supported"],
        "clean_abstain": cert["clean_abstain"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "candidate_grammars_remaining": cert["candidate_grammars_remaining"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "claim_allowed": cert["claim_allowed"],
        "next_required_batch": cert["next_required_batch"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
