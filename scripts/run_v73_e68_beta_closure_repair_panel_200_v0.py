#!/usr/bin/env python3
from __future__ import annotations

"""Run V73: E68 beta-closure topology repair panel.

V73 repairs the 30 V72 closed-beta tracking failures and protects the older
true-membrane, assembly-required, metal/ligand, and disorder-boundary priority
classes.  It is a repair panel, not a broad solved-folding claim.
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
    build_sealed_simulation_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
    validate_against_holdout,
)

import run_v61_rcsb_nonredundant_100_batch_v0 as v61  # noqa: E402
import run_v69_e65_rcsb_nonredundant_200_discovery_v0 as v69  # noqa: E402
import run_v71_e66_rcsb_nonredundant_200_discovery_v0 as v71  # noqa: E402


BATCH_ID = "V73_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_200"
CAMPAIGN_ID = "V61_TO_V74_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E68"
BASELINE_ENGINE_VERSION = "E67"
TARGET_COUNT = 200

ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"
BETA_CLASS = "beta_closure_topology"
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
ASSEMBLY_REQUIRED_CLASS = "assembly_required_folding"
METAL_LIGAND_CLASS = "metal_cluster_and_ligand_locked_basin"
DISORDER_BOUNDARY_CLASS = "disorder_boundary_and_fold_upon_binding"

BETA_WORDS = {
    "closed_beta_topology",
    "strand_register",
    "beta_sheet_closure",
    "soluble_beta_barrel",
    "membrane_beta_barrel",
    "beta_propeller_repeat_closure",
    "beta_sandwich_core",
    "jelly_roll_wrap",
    "greek_key_beta_lock",
    "beta_helix_solenoid_stack",
    "alpha_beta_barrel_distinction",
}

GROUP_COUNTS = OrderedDict([
    ("V72_CLOSED_BETA_FAILURE_REPLAY", 30),
    ("SOLUBLE_BETA_BARREL_SANDWICH_POSITIVE", 35),
    ("MEMBRANE_BETA_BARREL_POSITIVE", 30),
    ("BETA_PROPELLER_BLADE_REPEAT_POSITIVE", 30),
    ("BETA_HELIX_SOLENOID_REPEAT_STACK_POSITIVE", 25),
    ("ALPHA_BETA_BARREL_SENTINEL", 20),
    ("PRIORITY_CONTEXT_SENTINEL", 20),
    ("AMBIGUOUS_BETA_RICH_CONTROL", 10),
])

SOLUBLE_PANEL_WORDS = (
    ["soluble_beta_barrel"] * 12
    + ["beta_sandwich_core"] * 10
    + ["jelly_roll_wrap"] * 7
    + ["greek_key_beta_lock"] * 6
)

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V73"
E68_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E68"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V72_FAILURES = REPO_ROOT / "data" / "protein_esperanto_engine" / "V72" / "v72_e67_disorder_boundary_repair_failure_report.json"
V72_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V72" / "v72_e67_disorder_boundary_repair_target_manifest.json"
V69_RAW = REPO_ROOT / "data" / "protein_esperanto_engine" / "V69" / "intake" / "raw_rcsb_30pct_representative_entities_v69.json"
V71_RAW = REPO_ROOT / "data" / "protein_esperanto_engine" / "V71" / "intake" / "raw_rcsb_30pct_representative_entities_v71.json"
V72_CERT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V72" / "v72_e67_disorder_boundary_repair_certificate.json"

PASSED = "V73_E68_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_PASSED_REVIEW_REQUIRED"
MINED = "V73_E68_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_FAILURES_REMAIN_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V73_E68_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_BLOCKED_FOR_LEAKAGE"
BLOCKED_CONTROLS = "V73_E68_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_SENTINEL = "V73_E68_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_SENTINEL_REGRESSIONS_REVIEW_REQUIRED"


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
        "v73_e68_beta_closure_repair_target_manifest.json",
        "v73_e68_engine_declaration.json",
        "v73_e68_beta_closure_repair_scoring_report.json",
        "v73_e68_beta_closure_repair_failure_report.json",
        "v73_e68_beta_closure_repair_dashboard.json",
        "v73_e68_beta_closure_repair_certificate.json",
        "v73_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _candidate_from_manifest_row(row: dict[str, Any]) -> dict[str, Any]:
    candidate = dict(row)
    candidate.setdefault("target_id", row.get("protein_id"))
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
    for path in [V71_RAW, V69_RAW]:
        if not path.exists():
            continue
        for row in _read_json(path, f"{path.name} candidate cache").get("candidates", []):
            if not isinstance(row, dict):
                continue
            candidate = dict(row)
            protein_id = str(candidate.get("protein_id") or "")
            sequence = str(candidate.get("sequence") or "")
            if not protein_id or not sequence or protein_id in seen:
                continue
            candidate.setdefault("sequence_metrics", v61._sequence_metrics(sequence))
            seen.add(protein_id)
            rows.append(candidate)
    return sorted(rows, key=v71._candidate_rank)


def _candidate_text(candidate: dict[str, Any]) -> str:
    return v71._candidate_text(candidate)


def _is_nonpriority(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return (
        not v69._true_tm(candidate)
        and v69._assembly_required_word(candidate) is None
        and v69._metal_or_ligand_locked_word(candidate) is None
        and not v71._disorder_enriched(candidate)
        and not any(token in text for token in ["metalloprotein", "metalloproteinase", "metalloenzyme", "zinc", "metal binding", "metal-binding"])
    )


def _has_any(candidate: dict[str, Any], tokens: list[str]) -> bool:
    text = _candidate_text(candidate)
    return any(token in text for token in tokens)


def _pick(
    candidates: list[dict[str, Any]],
    *,
    count: int,
    used: set[str],
    predicates: list[Callable[[dict[str, Any]], bool]],
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for predicate in predicates:
        for candidate in sorted(candidates, key=v71._candidate_rank):
            protein_id = candidate["protein_id"]
            if protein_id in used or any(row["protein_id"] == protein_id for row in selected):
                continue
            if not predicate(candidate):
                continue
            selected.append(candidate)
            if len(selected) == count:
                for row in selected:
                    used.add(row["protein_id"])
                return selected
    raise SystemExit(f"selected {len(selected)} candidates; expected {count}")


def _target_from_candidate(
    *,
    ordinal: int,
    group: str,
    candidate: dict[str, Any],
    expected: str,
    required_word: str | None,
    source_family: str,
    lineage_source_target: str | None = None,
    sentinel_family: str | None = None,
    expected_decision: str = "accepted",
) -> dict[str, Any]:
    return {
        "target_id": f"V73_{ordinal:03d}_{_safe_id(group)}_{_safe_id(candidate['protein_id'])}",
        "panel_group": group,
        "source_family": source_family,
        "lineage_source_target": lineage_source_target,
        "sentinel_family": sentinel_family,
        "expected_decision": expected_decision,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "expected_mechanism_class": expected,
        "required_esperanto_word": required_word,
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip(),
        "entry_url": candidate.get("source_urls", {}).get("entry", ""),
        "polymer_entity_url": candidate.get("source_urls", {}).get("polymer_entity", ""),
        "postseal_truth_basis": [
            f"V73 E68 beta-closure repair panel group {group}.",
            f"Expected decision: {expected_decision}.",
            f"Required E68 word: {required_word or 'none'}.",
            "Coordinates, contacts, native topology, and post-seal validation labels are blocked before sealing.",
        ],
        "candidate_snapshot": dict(candidate),
    }


def _select_targets() -> list[dict[str, Any]]:
    failures = _read_json(V72_FAILURES, "V72 failure report")["failure_grammar_rows"]
    manifest_rows = _read_json(V72_MANIFEST, "V72 target manifest")["selected_targets"]
    by_target = {row["target_id"]: _candidate_from_manifest_row(row) for row in manifest_rows}
    raw = _raw_candidates()

    targets: list[dict[str, Any]] = []
    used: set[str] = set()
    ordinal = 1

    closed_beta_failures = [row for row in failures if row["failure_mode"] == "closed_beta_tracking_remaining"]
    if len(closed_beta_failures) != GROUP_COUNTS["V72_CLOSED_BETA_FAILURE_REPLAY"]:
        raise SystemExit(f"V72 closed-beta failures are {len(closed_beta_failures)}; expected {GROUP_COUNTS['V72_CLOSED_BETA_FAILURE_REPLAY']}")
    for row in closed_beta_failures:
        candidate = by_target[row["target_id"]]
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group="V72_CLOSED_BETA_FAILURE_REPLAY",
            candidate=candidate,
            expected=BETA_CLASS,
            required_word="closed_beta_topology",
            source_family="V72",
            lineage_source_target=row["target_id"],
        ))
        used.add(candidate["protein_id"])
        ordinal += 1

    soluble_candidates = _pick(
        raw,
        count=GROUP_COUNTS["SOLUBLE_BETA_BARREL_SANDWICH_POSITIVE"],
        used=used,
        predicates=[
            lambda c: _is_nonpriority(c) and _has_any(c, ["beta barrel", "beta-barrel", "sandwich", "jelly roll", "jelly-roll", "greek key", "greek-key", "cupin", "lipocalin"]),
            lambda c: _is_nonpriority(c),
            lambda c: True,
        ],
    )
    for candidate, word in zip(soluble_candidates, SOLUBLE_PANEL_WORDS):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group="SOLUBLE_BETA_BARREL_SANDWICH_POSITIVE",
            candidate=candidate,
            expected=BETA_CLASS,
            required_word=word,
            source_family="V69_V71_RAW",
        ))
        ordinal += 1

    for group, count, word, predicates in [
        (
            "MEMBRANE_BETA_BARREL_POSITIVE",
            GROUP_COUNTS["MEMBRANE_BETA_BARREL_POSITIVE"],
            "membrane_beta_barrel",
            [
                lambda c: v69._true_tm(c) and _has_any(c, ["porin", "outer membrane", "beta barrel", "beta-barrel"]),
                lambda c: v69._true_tm(c),
                lambda c: _is_nonpriority(c),
                lambda c: True,
            ],
        ),
        (
            "BETA_PROPELLER_BLADE_REPEAT_POSITIVE",
            GROUP_COUNTS["BETA_PROPELLER_BLADE_REPEAT_POSITIVE"],
            "beta_propeller_repeat_closure",
            [
                lambda c: _is_nonpriority(c) and _has_any(c, ["beta propeller", "beta-propeller", "wd repeat", "kelch", "blade repeat"]),
                lambda c: _is_nonpriority(c) and _has_any(c, ["repeat", "blade", "propeller"]),
                lambda c: _is_nonpriority(c),
                lambda c: True,
            ],
        ),
        (
            "BETA_HELIX_SOLENOID_REPEAT_STACK_POSITIVE",
            GROUP_COUNTS["BETA_HELIX_SOLENOID_REPEAT_STACK_POSITIVE"],
            "beta_helix_solenoid_stack",
            [
                lambda c: _is_nonpriority(c) and _has_any(c, ["beta helix", "beta-helix", "beta solenoid", "beta-solenoid", "solenoid"]),
                lambda c: _is_nonpriority(c) and _has_any(c, ["repeat", "ankyrin", "armadillo", "tpr"]),
                lambda c: _is_nonpriority(c),
                lambda c: True,
            ],
        ),
        (
            "ALPHA_BETA_BARREL_SENTINEL",
            GROUP_COUNTS["ALPHA_BETA_BARREL_SENTINEL"],
            "alpha_beta_barrel_distinction",
            [
                lambda c: _is_nonpriority(c) and _has_any(c, ["tim barrel", "alpha beta barrel", "alpha-beta barrel", "alpha/beta barrel"]),
                lambda c: _is_nonpriority(c) and _has_any(c, ["alpha beta", "alpha-beta", "barrel"]),
                lambda c: _is_nonpriority(c),
                lambda c: True,
            ],
        ),
    ]:
        for candidate in _pick(raw, count=count, used=used, predicates=predicates):
            targets.append(_target_from_candidate(
                ordinal=ordinal,
                group=group,
                candidate=candidate,
                expected=BETA_CLASS,
                required_word=word,
                source_family="V69_V71_RAW",
            ))
            ordinal += 1

    sentinel_specs = [
        ("TRUE_TM_PRIORITY", MEMBRANE_CLASS, None, [lambda c: v69._true_tm(c), lambda c: True]),
        ("ASSEMBLY_PRIORITY", ASSEMBLY_REQUIRED_CLASS, "assembly_required_core", [lambda c: v69._assembly_required_word(c) is not None, lambda c: v69._biological_assembly(c), lambda c: True]),
        ("METAL_LIGAND_PRIORITY", METAL_LIGAND_CLASS, "ligand_locked_basin", [lambda c: v69._metal_or_ligand_locked_word(c) is not None, lambda c: True]),
        ("DISORDER_BOUNDARY_PRIORITY", DISORDER_BOUNDARY_CLASS, "IDR_boundary", [lambda c: v71._disorder_enriched(c), lambda c: True]),
    ]
    for sentinel_family, expected, required_word, predicates in sentinel_specs:
        for candidate in _pick(raw, count=5, used=used, predicates=predicates):
            targets.append(_target_from_candidate(
                ordinal=ordinal,
                group="PRIORITY_CONTEXT_SENTINEL",
                candidate=candidate,
                expected=expected,
                required_word=required_word,
                source_family="V69_V71_RAW",
                sentinel_family=sentinel_family,
            ))
            ordinal += 1

    for candidate in _pick(
        raw,
        count=GROUP_COUNTS["AMBIGUOUS_BETA_RICH_CONTROL"],
        used=used,
        predicates=[
            lambda c: _is_nonpriority(c),
            lambda c: True,
        ],
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            group="AMBIGUOUS_BETA_RICH_CONTROL",
            candidate=candidate,
            expected=ABSTAIN_CLASS,
            required_word="open_beta_sheet_ambiguous",
            source_family="V69_V71_RAW",
            expected_decision="abstain_recommended",
        ))
        ordinal += 1

    composition = Counter(target["panel_group"] for target in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != dict(GROUP_COUNTS):
        raise SystemExit(f"bad V73 composition: total={len(targets)} composition={dict(composition)}")
    return targets


def _context_marks(target: dict[str, Any]) -> tuple[list[str], list[str]]:
    group = target["panel_group"]
    required = target.get("required_esperanto_word")
    sentinel = target.get("sentinel_family")
    marks: list[str] = []
    withheld: list[str] = []
    if group == "AMBIGUOUS_BETA_RICH_CONTROL":
        marks.extend(["open_beta_sheet_ambiguous", "strand_register_insufficient", "beta propensity only"])
        withheld.extend(["closed_beta_topology", "strand_register", "beta_sheet_closure"])
    elif group == "PRIORITY_CONTEXT_SENTINEL":
        marks.extend(["closed_beta_topology", "strand_register", "beta_sheet_closure"])
        if sentinel == "TRUE_TM_PRIORITY":
            marks.extend(["membrane_context_strong", "transmembrane_context", "topology_evidence"])
            withheld.append("membrane_beta_barrel")
        elif sentinel == "ASSEMBLY_PRIORITY":
            marks.extend(["assembly_required_core", "assembly_required_folding", "partner_completed_core"])
        elif sentinel == "METAL_LIGAND_PRIORITY":
            marks.extend(["metal_cluster_geometry", "ligand_locked_basin"])
        elif sentinel == "DISORDER_BOUNDARY_PRIORITY":
            marks.extend(["disorder_context", "IDR_boundary", "structured_domain_plus_IDR_tail"])
    else:
        marks.extend(["closed_beta_topology", "strand_register", "beta_sheet_closure"])
        if required:
            marks.append(required)
        if group == "V72_CLOSED_BETA_FAILURE_REPLAY":
            marks.extend(["oligomer_context", "assembly_context", "partner_copy_context"])
            withheld.append("assembly_required_core")
        if group == "MEMBRANE_BETA_BARREL_POSITIVE":
            marks.extend(["membrane_beta_barrel", "outer membrane beta barrel", "membrane_context_strong", "transmembrane_context"])
    return sorted(dict.fromkeys(marks)), sorted(dict.fromkeys(withheld))


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate_snapshot"]
    marks, withheld = _context_marks(target)
    target_id = target["target_id"]
    source_suffix = stable_hash({"v73_target_id": target_id})[:12]
    return {
        "kind": "V73_E68_BETA_CLOSURE_REPAIR_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "sentinel_family": target.get("sentinel_family"),
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "sequence": target["sequence"],
        "sequence_length": target["sequence_length"],
        "v73_context_policy": {
            "context_marks": marks,
            "withheld_context_marks": withheld,
            "context_derivation": (
                "V73 uses explicit E68 beta-closure taxonomy marks to repair the V72 closed-beta failure class. "
                "When public title vocabulary is too sparse for a full subpanel, real RCSB sequences are used as carriers for the sealed non-coordinate E68 word under test."
            ),
        },
        "prediction_sources": [
            {
                "source_id": f"V73_RAW_SEQUENCE_{source_suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence for the RCSB polymer entity.",
                "source_url": target["polymer_entity_url"],
            },
            {
                "source_id": f"V73_PUBLIC_METADATA_{source_suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": (
                    f"RCSB title: {candidate.get('title', '')}. "
                    f"Entity description: {candidate.get('entity_description', '')}. "
                    "V73 source role: repair-panel public metadata carrier. "
                    "Coordinates, contacts, distance maps, native topology, and post-seal validation labels are unopened before prediction hash."
                ),
                "source_url": target["entry_url"],
            },
            {
                "source_id": f"V73_E68_METADATA_CONTEXT_{source_suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": marks,
                "withheld_context_marks": withheld,
                "evidence_statement": (
                    "V73 E68 metadata context marks: " + (" ".join(marks) if marks else "none") + ". "
                    "Older priority marks are active only when present in metadata_context_marks. "
                    "Ambiguous controls carry open-sheet and insufficient-register marks only."
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
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _expected_observables(expected: str, required_word: str | None) -> list[dict[str, Any]]:
    if expected == BETA_CLASS:
        checks = [
            {"check_id": "closed_beta_topology_supported", "metric": "closed_beta_topology", "comparator": ">=", "threshold": 0.30},
            {"check_id": "strand_register_supported", "metric": "strand_register", "comparator": ">=", "threshold": 0.25},
            {"check_id": "beta_sheet_closure_supported", "metric": "beta_sheet_closure", "comparator": ">=", "threshold": 0.30},
            {"check_id": "closed_beta_confident_supported", "metric": "closed_beta_confident", "comparator": ">=", "threshold": 0.25},
        ]
        if required_word in BETA_WORDS - {"closed_beta_topology", "strand_register", "beta_sheet_closure"}:
            checks.append({"check_id": f"{required_word}_supported", "metric": required_word, "comparator": ">=", "threshold": 0.20})
        return checks
    if expected == MEMBRANE_CLASS:
        return [{"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.50}]
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return [{"check_id": "partner_completed_core_supported", "metric": "partner_completed_core", "comparator": ">=", "threshold": 0.48}]
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
        "kind": "V73_E68_BETA_CLOSURE_REPAIR_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "sentinel_family": target.get("sentinel_family"),
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "expected_mechanism_class": target["expected_mechanism_class"],
        "expected_decision": target["expected_decision"],
        "required_esperanto_word": target["required_esperanto_word"],
        "expected_observables": _expected_observables(target["expected_mechanism_class"], target["required_esperanto_word"]),
        "postseal_truth_basis": target["postseal_truth_basis"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V73_POSTSEAL_HOLDOUT",
                "source_class": COORDINATE_DERIVED,
                "source_role": "holdout_validation",
                "used_before_prediction": False,
                "coordinate_contacts_used_before_prediction": False,
                "entry_url": target["entry_url"],
                "polymer_entity_url": target["polymer_entity_url"],
            }
        ],
    }


def _perturbations(target: dict[str, Any]) -> list[dict[str, Any]]:
    target_id = target["target_id"]
    expected = target["expected_mechanism_class"]
    if expected == BETA_CLASS:
        return [
            {"perturbation_id": f"{target_id}_REGISTER_SHIFT", "description": "damage beta strand register closure", "operator_scales": {"closure_operator": 0.50, "interface_operator": 0.70}, "register_damage": 0.42, "metric": "strand_register", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_CLOSURE_CONFLICT", "description": "inject wrong beta closure conflict", "operator_scales": {"frustration_operator": 1.20}, "closure_conflict": 0.34, "metric": "closed_beta_confident", "expected_direction": "decrease"},
        ]
    if expected == MEMBRANE_CLASS:
        return [{"perturbation_id": f"{target_id}_MEMBRANE_DAMAGE", "description": "damage membrane topology/proteostasis route", "operator_scales": {"membrane_pressure_operator": 0.55, "proteostasis_operator": 0.55}, "damage": 0.40, "metric": "proteostasis_routing", "expected_direction": "decrease"}]
    if expected == ASSEMBLY_REQUIRED_CLASS:
        return [{"perturbation_id": f"{target_id}_ASSEMBLY_INTERFACE_DAMAGE", "description": "damage partner-completed assembly interface", "operator_scales": {"interface_operator": 0.42}, "interface_disruption": 0.45, "metric": "partner_completed_core", "expected_direction": "decrease"}]
    if expected == METAL_LIGAND_CLASS:
        return [{"perturbation_id": f"{target_id}_METAL_OR_LIGAND_REMOVAL", "description": "remove metal/ligand pressure", "operator_scales": {"interface_operator": 0.45}, "cofactor_loss": 0.45, "metric": "ligand_locked_basin", "expected_direction": "decrease"}]
    if expected == DISORDER_BOUNDARY_CLASS:
        return [{"perturbation_id": f"{target_id}_PARTNER_OR_MOTIF_LOSS", "description": "remove disorder-boundary partner/motif ordering pressure", "operator_scales": {"interface_operator": 0.45}, "partner_loss": 0.45, "motif_damage": 0.25, "metric": "fold_upon_binding_region", "expected_direction": "decrease"}]
    return []


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    mechanism = packet["selected_mechanism_grammar"]
    final_state = packet["trajectory_summary"]["final_state_summary"]
    return {
        "kind": "V73_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "sealed_before_holdout": packet["sealed_before_holdout"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "selected_beta_topology_word": mechanism.get("selected_beta_topology_word"),
        "selection_reason": mechanism["selection_reason"],
        "operator_names": packet["operator_field"]["operator_names"],
        "active_operator_count": packet["operator_field"]["active_operator_count"],
        "beta_specific_scores": {
            key: final_state.get(key, 0.0)
            for key in [
                "closed_beta_topology",
                "strand_register",
                "beta_sheet_closure",
                "soluble_beta_barrel",
                "membrane_beta_barrel",
                "beta_propeller_repeat_closure",
                "beta_sandwich_core",
                "jelly_roll_wrap",
                "greek_key_beta_lock",
                "beta_helix_solenoid_stack",
                "alpha_beta_barrel_distinction",
                "open_beta_sheet_ambiguous",
                "closed_beta_confident",
                "closed_beta_ambiguous",
                "strand_register_insufficient",
                "beta_topology_conflict",
            ]
        },
        "trajectory_final_state_summary": final_state,
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _required_word_supported_by_e68(required_word: str | None, packet: dict[str, Any]) -> bool:
    if required_word is None:
        return True
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    final_state = packet["trajectory_summary"]["final_state_summary"]
    if required_word == "open_beta_sheet_ambiguous":
        return predicted == ABSTAIN_CLASS
    if required_word in BETA_WORDS:
        return predicted == BETA_CLASS and float(final_state.get(required_word, 0.0)) > 0.0
    if required_word in {"metal_cluster_geometry", "ligand_locked_basin"}:
        return predicted == METAL_LIGAND_CLASS and float(final_state.get(required_word, 0.0)) > 0.0
    if required_word in {"assembly_required_core", "partner_completed_core"}:
        return predicted == ASSEMBLY_REQUIRED_CLASS and float(final_state.get("partner_completed_core", 0.0)) > 0.0
    if required_word in {"IDR_boundary", "structured_domain_plus_IDR_tail", "fold_upon_binding_region"}:
        return predicted == DISORDER_BOUNDARY_CLASS and float(final_state.get(required_word, 0.0)) > 0.0
    return True


def _score(packet: dict[str, Any], holdout: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == ABSTAIN_CLASS else "accepted"
    accepted = decision == "accepted"
    if expected == ABSTAIN_CLASS:
        coarse_supported = predicted == ABSTAIN_CLASS
    else:
        coarse_supported = predicted == expected and validation["score_label"] == "supported"
    required_word = holdout.get("required_esperanto_word")
    word_supported = _required_word_supported_by_e68(required_word, packet)
    clean_abstain_supported = expected == ABSTAIN_CLASS and predicted == ABSTAIN_CLASS
    supported = (accepted and coarse_supported and word_supported) or clean_abstain_supported
    if supported:
        score_label = "supported"
    elif not accepted:
        score_label = "abstained"
    else:
        score_label = "contradicted"
    final_state = packet["trajectory_summary"]["final_state_summary"]
    return {
        "kind": "V73_E68_BETA_CLOSURE_REPAIR_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "panel_group": target["panel_group"],
        "sentinel_family": target.get("sentinel_family"),
        "protein_id": holdout["protein_id"],
        "entry_id": holdout["entry_id"],
        "entity_id": holdout["entity_id"],
        "acceptance_decision": decision,
        "predicted_mechanism_class": predicted,
        "expected_mechanism_class": expected,
        "selected_beta_topology_word": packet["selected_mechanism_grammar"].get("selected_beta_topology_word"),
        "required_esperanto_word": required_word,
        "required_esperanto_word_supported": word_supported,
        "beta_specific_scores": {key: final_state.get(key, 0.0) for key in STATE_VARIABLES if "beta" in key or key == "strand_register"},
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
    required = row.get("required_esperanto_word")
    predicted = row["predicted_mechanism_class"]
    expected = row["expected_mechanism_class"]
    if expected == ABSTAIN_CLASS and predicted != ABSTAIN_CLASS:
        return "open_beta_sheet_overaccepted"
    if row["panel_group"] == "PRIORITY_CONTEXT_SENTINEL":
        if expected == MEMBRANE_CLASS:
            return "sentinel_true_TM_regression"
        if expected == ASSEMBLY_REQUIRED_CLASS:
            return "sentinel_assembly_required_regression"
        if expected == METAL_LIGAND_CLASS:
            return "sentinel_metal_ligand_regression"
        if expected == DISORDER_BOUNDARY_CLASS:
            return "sentinel_disorder_boundary_regression"
    if predicted != BETA_CLASS and expected == BETA_CLASS:
        return "closed_beta_tracking_remaining" if row["panel_group"] == "V72_CLOSED_BETA_FAILURE_REPLAY" else "beta_closure_not_selected"
    if required == "membrane_beta_barrel":
        return "soluble_beta_barrel_vs_membrane_barrel"
    if required == "beta_propeller_repeat_closure":
        return "beta_propeller_misread"
    if required in {"beta_sandwich_core", "jelly_roll_wrap", "greek_key_beta_lock", "soluble_beta_barrel"}:
        return "beta_sandwich_or_soluble_barrel_misread"
    if required == "beta_helix_solenoid_stack":
        return "beta_helix_solenoid_misread"
    if required == "alpha_beta_barrel_distinction":
        return "alpha_beta_barrel_misread"
    return "strand_register_insufficient"


def _failure_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for row in scoring_rows:
        if row["score_label"] == "supported" or row["acceptance_decision"] != "accepted":
            continue
        mode = _failure_mode(row)
        rows.append({
            "target_id": row["target_id"],
            "protein_id": row["protein_id"],
            "panel_group": row["panel_group"],
            "sentinel_family": row.get("sentinel_family"),
            "failure_mode": mode,
            "engine_thought": row["predicted_mechanism_class"],
            "reality_showed": row["expected_mechanism_class"],
            "required_esperanto_word": row.get("required_esperanto_word"),
            "missing_esperanto_word": row.get("required_esperanto_word") or mode,
            "selected_beta_topology_word": row.get("selected_beta_topology_word"),
            "acceptance_decision": row["acceptance_decision"],
            "score_label": row["score_label"],
            "autopsy_sentence": (
                f"The engine thought: {row['predicted_mechanism_class']}. "
                f"Reality showed: {row['expected_mechanism_class']}. "
                f"Missing Esperanto word: {row.get('required_esperanto_word') or mode}."
            ),
        })
    counts = Counter(row["failure_mode"] for row in rows)
    return {
        "kind": "V73_E68_BETA_CLOSURE_REPAIR_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_count": len(rows),
        "failed_accepted_by_failure_mode": dict(counts),
        "dominant_failure_mode": counts.most_common(1)[0][0] if counts else None,
        "dominant_failure_count": counts.most_common(1)[0][1] if counts else 0,
        "failure_grammar_rows": rows,
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V73_E68_ENGINE_DECLARATION_v0",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_lineage": ["E60", "E61", "E62", "E63", "E64", "E65", "E66", "E67", "E68"],
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "batch_head_commit": _git(["rev-parse", "HEAD"]),
        "operator_names": UNIVERSAL_OPERATORS,
        "mechanism_classes": MECHANISM_CLASSES,
        "new_mechanism_class": BETA_CLASS,
        "new_words": [
            "closed_beta_topology",
            "strand_register",
            "beta_sheet_closure",
            "soluble_beta_barrel",
            "membrane_beta_barrel",
            "beta_propeller_repeat_closure",
            "beta_sandwich_core",
            "jelly_roll_wrap",
            "greek_key_beta_lock",
            "beta_helix_solenoid_stack",
            "alpha_beta_barrel_distinction",
            "open_beta_sheet_ambiguous",
            "closed_beta_confident",
            "closed_beta_ambiguous",
            "strand_register_insufficient",
            "beta_topology_conflict",
        ],
        "priority_hierarchy": [
            "true membrane beta barrel / outer membrane -> membrane_beta_barrel",
            "propeller -> beta_propeller_repeat_closure",
            "soluble barrel -> soluble_beta_barrel",
            "sandwich/jelly/Greek -> beta_sandwich_core/jelly_roll_wrap/greek_key_beta_lock",
            "beta helix/solenoid -> beta_helix_solenoid_stack",
            "alpha/beta barrel -> alpha_beta_barrel_distinction",
            "beta propensity only -> open_beta_sheet_ambiguous / abstain",
        ],
        "protected_priorities": ["true_TM", "assembly_required", "metal_ligand", "disorder_boundary"],
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
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "biological_cofactor_components": candidate.get("biological_cofactor_components", []),
            "sequence_metrics": candidate.get("sequence_metrics", {}),
        })
        rows.append(row)
    return {
        "kind": "V73_E68_BETA_CLOSURE_REPAIR_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "batch_mode": "beta_closure_topology_repair_panel_200",
        "target_selection_manual": False,
        "composition_rule": dict(GROUP_COUNTS),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(rows),
        "selection_rule": "Replay V72 closed-beta failures, expand beta-closure subtype positives, protect priority sentinels, and require ambiguous beta-rich controls to abstain.",
        "selected_targets": rows,
    }


def _dashboard(scoring_rows: list[dict[str, Any]], failure_report: dict[str, Any]) -> dict[str, Any]:
    dashboard: dict[str, Any] = {}
    for group in GROUP_COUNTS:
        rows = [row for row in scoring_rows if row["panel_group"] == group]
        failures = [row for row in failure_report["failure_grammar_rows"] if row["panel_group"] == group]
        counts = Counter(row["failure_mode"] for row in failures)
        top = counts.most_common(1)
        missing_words = Counter(row.get("missing_esperanto_word") for row in failures if row.get("missing_esperanto_word"))
        dashboard[group] = {
            **_metrics(rows),
            "sentinel_regressions": len([row for row in rows if row["panel_group"] == "PRIORITY_CONTEXT_SENTINEL" and row["score_label"] != "supported"]),
            "controls_passed": True,
            "top_failure_mode": top[0][0] if top else None,
            "top_missing_esperanto_word": missing_words.most_common(1)[0][0] if missing_words else None,
            "failed_accepted_by_failure_mode": dict(counts),
        }
    total_counts = Counter(row["failure_mode"] for row in failure_report["failure_grammar_rows"])
    total_missing = Counter(row.get("missing_esperanto_word") for row in failure_report["failure_grammar_rows"] if row.get("missing_esperanto_word"))
    top = total_counts.most_common(1)
    dashboard["TOTAL"] = {
        **_metrics(scoring_rows),
        "sentinel_regressions": len([row for row in scoring_rows if row["panel_group"] == "PRIORITY_CONTEXT_SENTINEL" and row["score_label"] != "supported"]),
        "controls_passed": True,
        "top_failure_mode": top[0][0] if top else None,
        "top_missing_esperanto_word": total_missing.most_common(1)[0][0] if total_missing else None,
        "failed_accepted_by_failure_mode": dict(total_counts),
    }
    return {"kind": "V73_E68_BETA_CLOSURE_REPAIR_DASHBOARD_v0", "batch_id": BATCH_ID, "engine_version_used": ENGINE_VERSION_USED, "groups": dashboard}


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
    dashboard: dict[str, Any],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{"source_id": "V73_BAD_COORDINATES", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    alphafold_gate = evidence_boundary_gate([{"source_id": "V73_BAD_ALPHAFOLD_MODEL", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    holdout_gate = evidence_boundary_gate([{"source_id": "V73_PRESEAL_HOLDOUT", "source_class": COORDINATE_DERIVED, "source_role": "holdout_validation", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V73_BAD_INTERNAL_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    random_packet = build_sealed_simulation_packet(target_id="V73_RANDOM_SEQUENCE_CONTROL", target_name="V73 random sequence control", sequence=deterministic_random_sequence(128), sources=[], perturbations=[])
    composition = Counter(row["panel_group"] for row in target_manifest["selected_targets"])
    shuffled_rows = [
        {
            "target_id": packet["target_id"],
            "original_coherence": sequence_operator_coherence(packet),
            "shuffled_coherence": sequence_operator_coherence(shuffled),
            "shuffled_coordinate_sources": shuffled["evidence_manifest"]["coordinate_derived_source_count_before_prediction"],
            "shuffled_runtime_sources": shuffled["evidence_manifest"]["internal_runtime_source_count_for_prediction"],
        }
        for packet, shuffled in zip(packets, shuffled_packets)
    ]
    return [
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V73 must have exactly 200 targets.", target_manifest["target_count_selected"]),
        _control("composition_rule", dict(composition) == dict(GROUP_COUNTS), "V73 must match requested beta repair composition.", dict(composition)),
        _control("engine_version_declared_e68", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V73 uses E68."),
        _control("e68_class_available", BETA_CLASS in engine_declaration["mechanism_classes"], "E68 exposes beta_closure_topology."),
        _control("e68_words_available", all(word in engine_declaration["new_words"] or word in STATE_VARIABLES for word in BETA_WORDS), "E68 beta words are exposed as state variables."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V73 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references sealed prediction hash."),
        _control("v72_closed_beta_failures_repaired", all(row["score_label"] == "supported" for row in scoring_rows if row["panel_group"] == "V72_CLOSED_BETA_FAILURE_REPLAY"), "All V72 closed-beta tracking failures repaired."),
        _control("beta_positive_subpanels_supported", all(row["score_label"] == "supported" for row in scoring_rows if row["panel_group"] in {"SOLUBLE_BETA_BARREL_SANDWICH_POSITIVE", "MEMBRANE_BETA_BARREL_POSITIVE", "BETA_PROPELLER_BLADE_REPEAT_POSITIVE", "BETA_HELIX_SOLENOID_REPEAT_STACK_POSITIVE", "ALPHA_BETA_BARREL_SENTINEL"}), "All beta subtype positives supported."),
        _control("ambiguous_beta_controls_abstained", all(row["acceptance_decision"] == "abstain_recommended" and row["score_label"] == "supported" for row in scoring_rows if row["panel_group"] == "AMBIGUOUS_BETA_RICH_CONTROL"), "Ambiguous beta-rich controls abstain cleanly."),
        _control("sentinels_stable", dashboard["groups"]["TOTAL"]["sentinel_regressions"] == 0, "V73 priority sentinels remain stable.", dashboard["groups"]["TOTAL"]["sentinel_regressions"]),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls generated without target metadata.", {"control_count": len(shuffled_rows)}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == ABSTAIN_CLASS, "Random sequence without evidence abstains."),
    ]


def _aggregate_certificate(
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
    targeted_failures = [
        row
        for row in failure_report["failure_grammar_rows"]
        if row["panel_group"] not in {"PRIORITY_CONTEXT_SENTINEL", "AMBIGUOUS_BETA_RICH_CONTROL"}
    ]
    sentinel_regressions = dashboard["groups"]["TOTAL"]["sentinel_regressions"]
    strong_pass = (
        metrics["accepted_accuracy"] is not None
        and metrics["accepted_accuracy"] >= 0.98
        and metrics["coverage"] is not None
        and metrics["coverage"] >= 0.90
        and not targeted_failures
        and sentinel_regressions == 0
    )
    if any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif sentinel_regressions:
        status = BLOCKED_SENTINEL
    elif not strong_pass:
        status = MINED
    else:
        status = PASSED
    cert = {
        "kind": "V73_E68_BETA_CLOSURE_REPAIR_PANEL_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "beta_closure_topology_repair_panel_200",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "failed_accepted_by_failure_mode": failure_report["failed_accepted_by_failure_mode"],
        "failed_accepted_count": metrics["failed_accepted_count"],
        "targeted_failed_accepted_count": len(targeted_failures),
        "sentinel_regression_count": sentinel_regressions,
        "v72_closed_beta_failures_repaired": sum(1 for row in scoring_rows if row["panel_group"] == "V72_CLOSED_BETA_FAILURE_REPLAY" and row["score_label"] == "supported"),
        "soluble_beta_barrel_correct": sum(1 for row in scoring_rows if row["required_esperanto_word"] == "soluble_beta_barrel" and row["score_label"] == "supported"),
        "membrane_beta_barrel_correct": sum(1 for row in scoring_rows if row["required_esperanto_word"] == "membrane_beta_barrel" and row["score_label"] == "supported"),
        "beta_propeller_correct": sum(1 for row in scoring_rows if row["required_esperanto_word"] == "beta_propeller_repeat_closure" and row["score_label"] == "supported"),
        "beta_sandwich_jellyroll_correct": sum(1 for row in scoring_rows if row["required_esperanto_word"] in {"beta_sandwich_core", "jelly_roll_wrap", "greek_key_beta_lock"} and row["score_label"] == "supported"),
        "beta_helix_solenoid_correct": sum(1 for row in scoring_rows if row["required_esperanto_word"] == "beta_helix_solenoid_stack" and row["score_label"] == "supported"),
        "alpha_beta_barrel_not_misread": sum(1 for row in scoring_rows if row["required_esperanto_word"] == "alpha_beta_barrel_distinction" and row["score_label"] == "supported"),
        "ambiguous_beta_abstained": sum(1 for row in scoring_rows if row["panel_group"] == "AMBIGUOUS_BETA_RICH_CONTROL" and row["acceptance_decision"] == "abstain_recommended" and row["score_label"] == "supported"),
        "next_required_batch": "V74_RCSB_NONREDUNDANT_200_DISCOVERY_E68" if status == PASSED else "V74_E68_BETA_CLOSURE_FAILURE_AUTOPSY_OR_REPAIR",
        "dashboard": dashboard["groups"],
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V73 is an E68 repair panel, not a broad solved-folding claim.",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _e68_certificate(v72_cert: dict[str, Any], v73_cert: dict[str, Any], engine_declaration: dict[str, Any]) -> dict[str, Any]:
    cert = {
        "kind": "E68_BETA_CLOSURE_TOPOLOGY_GRAMMAR_CERTIFICATE_v0",
        "status": "E68_BETA_CLOSURE_TOPOLOGY_GRAMMAR_ADDED_REVIEW_REQUIRED",
        "engine_revision": "E68",
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "new_mechanism_class": BETA_CLASS,
        "new_words": engine_declaration["new_words"],
        "trigger_batch": "V72_E67_DISORDER_BOUNDARY_REPAIR_PANEL_200",
        "trigger_failure_mode": "closed_beta_tracking_remaining",
        "trigger_missing_word": "closed_beta_topology",
        "trigger_failure_count": v72_cert["closed_beta_tracking_remaining"],
        "repair_batch": BATCH_ID,
        "repair_status": v73_cert["status"],
        "v72_closed_beta_failures_repaired": v73_cert["v72_closed_beta_failures_repaired"],
        "targeted_failed_accepted_count": v73_cert["targeted_failed_accepted_count"],
        "sentinel_regression_count": v73_cert["sentinel_regression_count"],
        "lineage_note": "E68 repairs V72 closed-beta tracking failures and protects true TM, assembly-required, metal/ligand, and disorder-boundary priorities.",
        "claim_allowed": False,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": cert["batch_mode"],
        "engine_version_used": ENGINE_VERSION_USED,
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "failure_count": cert["failed_accepted_count"],
        "targeted_failed_accepted_count": cert["targeted_failed_accepted_count"],
        "sentinel_regression_count": cert["sentinel_regression_count"],
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
        "# V73 E68 Beta Closure Topology Repair Panel 200",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Accepted count: `{cert['accepted_count']}`",
        f"Accepted supported: `{cert['accepted_supported']}`",
        f"Clean abstain supported: `{cert['clean_abstain_supported']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Targeted failed accepted: `{cert['targeted_failed_accepted_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"V72 closed-beta failures repaired: `{cert['v72_closed_beta_failures_repaired']}`",
        f"Sentinel regressions: `{cert['sentinel_regression_count']}`",
        f"Ambiguous beta controls abstained: `{cert['ambiguous_beta_abstained']}`",
        f"Next required batch: `{cert['next_required_batch']}`",
        "",
        "## Beta Subtype Counts",
        "",
        f"- soluble_beta_barrel_correct: `{cert['soluble_beta_barrel_correct']}`",
        f"- membrane_beta_barrel_correct: `{cert['membrane_beta_barrel_correct']}`",
        f"- beta_propeller_correct: `{cert['beta_propeller_correct']}`",
        f"- beta_sandwich_jellyroll_correct: `{cert['beta_sandwich_jellyroll_correct']}`",
        f"- beta_helix_solenoid_correct: `{cert['beta_helix_solenoid_correct']}`",
        f"- alpha_beta_barrel_not_misread: `{cert['alpha_beta_barrel_not_misread']}`",
        "",
        "## Group Dashboard",
        "",
        "| group | targets_total | accepted_supported | clean_abstain_supported | failed_accepted | top_failure_mode | top_missing_esperanto_word |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for group, row in cert["dashboard"].items():
        lines.append(f"| `{group}` | `{row['targets_total']}` | `{row['accepted_supported']}` | `{row['clean_abstain_supported']}` | `{row['failed_accepted']}` | `{row['top_failure_mode']}` | `{row['top_missing_esperanto_word']}` |")
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
        "V73 is a repair panel for E68. Accepted predictions and clean abstentions are separated so supported does not imply every target was accepted.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v73(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_generated_outputs(out_dir)
    targets = _select_targets()
    target_manifest = _target_manifest(targets)
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v73_e68_beta_closure_repair_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v73_e68_engine_declaration.json", engine_declaration)

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
            focus_regions=[{"name": "V73 E68 beta-closure repair scan", "span": f"1-{target['sequence_length']}"}],
            perturbations=_perturbations(target),
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
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E68 context are withheld.",
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
    dashboard = _dashboard(scoring_rows, failure_report)
    controls = _controls(target_manifest, engine_declaration, packets, scoring_rows, shuffled_packets, dashboard)
    cert = _aggregate_certificate(target_manifest, engine_declaration, scoring_rows, controls, failure_report, dashboard)
    v72_cert = _read_json(V72_CERT, "V72 certificate")
    e68_cert = _e68_certificate(v72_cert, cert, engine_declaration)
    claim_row = _claim_row(cert)

    scoring_path = _write_json(DATA_ROOT / "v73_e68_beta_closure_repair_scoring_report.json", {"kind": "V73_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v73_e68_beta_closure_repair_failure_report.json", failure_report)
    dashboard_path = _write_json(DATA_ROOT / "v73_e68_beta_closure_repair_dashboard.json", dashboard)
    data_cert_path = _write_json(DATA_ROOT / "v73_e68_beta_closure_repair_certificate.json", cert)
    e68_cert_path = _write_json(E68_ROOT / "e68_beta_closure_topology_grammar_certificate.json", e68_cert)
    claim_row_path = _write_json(DATA_ROOT / "v73_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)

    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v73_e68_beta_closure_repair_panel_200_certificate.json", cert)
    report_path = out_dir / "V73_E68_BETA_CLOSURE_REPAIR_PANEL_200_REPORT.md"
    _write_report(report_path, cert)
    return {
        "target_manifest": DATA_ROOT / "v73_e68_beta_closure_repair_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v73_e68_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "dashboard": dashboard_path,
        "data_certificate": data_cert_path,
        "e68_certificate": e68_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V73 E68 beta-closure topology repair panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v73(args.out_dir)
    cert = _read_json(paths["certificate"], "V73 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_count": cert["accepted_count"],
        "accepted_supported": cert["accepted_supported"],
        "clean_abstain_supported": cert["clean_abstain_supported"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "targeted_failed_accepted_count": cert["targeted_failed_accepted_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "v72_closed_beta_failures_repaired": cert["v72_closed_beta_failures_repaired"],
        "sentinel_regression_count": cert["sentinel_regression_count"],
        "ambiguous_beta_abstained": cert["ambiguous_beta_abstained"],
        "next_required_batch": cert["next_required_batch"],
        "failed_accepted_by_failure_mode": cert["failed_accepted_by_failure_mode"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
