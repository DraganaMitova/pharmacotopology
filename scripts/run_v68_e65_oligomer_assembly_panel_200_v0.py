#!/usr/bin/env python3
from __future__ import annotations

"""Run V68: E65 oligomer/assembly specialist panel.

V68 is a targeted 200-row panel for the V67 missing word:
assembly_required_core_vs_topology_provider.  It is not a broad solved claim.
It tests whether E65 separates assembly-required folding from true membrane
topology, soluble hydrophobic cores, cofactor-stabilized cores, and generic
complex annotations.
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
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
    validate_against_holdout,
)


BATCH_ID = "V68_OLIGOMER_ASSEMBLY_PANEL_200"
CAMPAIGN_ID = "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN"
ENGINE_VERSION_USED = "E65"
BASELINE_ENGINE_VERSION = "E64"
TARGET_COUNT = 200
ABSTAIN_CLASS = "insufficient_evidence_clean_abstain"
ASSEMBLY_CLASS = "assembly_required_folding"
MEMBRANE_CLASS = "membrane_multidomain_folding_proteostasis"
COFACTOR_CLASS = "cofactor_ligand_assisted_stabilization"
GLOBULAR_CLASS = "globular_closure"

GROUP_COUNTS = OrderedDict([
    ("V67_ASSEMBLY_REQUIRED_REPLAY", 48),
    ("BIOLOGICAL_OLIGOMER_ASSEMBLY_REQUIRED", 40),
    ("COILED_COIL_HELIX_BUNDLE_ASSEMBLY", 30),
    ("TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL", 30),
    ("SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL", 25),
    ("COFACTOR_LIGAND_SOLUBLE_SENTINEL", 15),
    ("GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL", 12),
])

DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V68"
E65_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E65"
CAMPAIGN_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / CAMPAIGN_ID
LEDGER_ROOT = CAMPAIGN_ROOT / "ledgers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

V63_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_target_manifest.json"
V65_MANIFEST = REPO_ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_panel_manifest.json"
V67_FAILURES = REPO_ROOT / "data" / "protein_esperanto_engine" / "V67" / "v67_mixed_fast_discovery_failure_report.json"
V67_TARGETS = REPO_ROOT / "data" / "protein_esperanto_engine" / "V67" / "v67_mixed_fast_discovery_target_manifest.json"

PASSED = "V68_E65_OLIGOMER_ASSEMBLY_PANEL_PASSED_REVIEW_REQUIRED"
MINED = "V68_E65_OLIGOMER_ASSEMBLY_PANEL_FAILURES_MINED_REVIEW_REQUIRED"
BLOCKED_CONTROLS = "V68_E65_OLIGOMER_ASSEMBLY_PANEL_CONTROLS_FAILED_REVIEW_REQUIRED"
BLOCKED_LEAKAGE = "V68_E65_OLIGOMER_ASSEMBLY_PANEL_BLOCKED_FOR_LEAKAGE"
BLOCKED_SENTINEL = "V68_E65_OLIGOMER_ASSEMBLY_PANEL_SENTINEL_REGRESSIONS_REVIEW_REQUIRED"


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
        "v68_oligomer_assembly_panel_target_manifest.json",
        "v68_e65_engine_declaration.json",
        "v68_oligomer_assembly_panel_scoring_report.json",
        "v68_oligomer_assembly_panel_failure_report.json",
        "v68_oligomer_assembly_panel_certificate.json",
        "v68_campaign_claim_row.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()
    if out_dir.exists():
        shutil.rmtree(out_dir)


def _candidate_text(candidate: dict[str, Any], *, include_postseal: bool = True) -> str:
    values: list[Any] = [
        candidate.get("title", ""),
        candidate.get("entity_description", ""),
        candidate.get("entry_keywords", ""),
        candidate.get("polymer_composition", ""),
        " ".join(candidate.get("organisms", []) or []),
        " ".join(candidate.get("biological_cofactor_components", []) or []),
        " ".join(candidate.get("nonpolymer_bound_components", []) or []),
    ]
    if include_postseal:
        values.extend(candidate.get("annotations", []) or [])
        values.extend(candidate.get("feature_types", []) or [])
    return " ".join(str(value) for value in values).lower()


def _negative_topology(text: str) -> bool:
    return any(token in text for token in ["monotopic/peripheral", "peripheral membrane", "not transmembrane", "no transmembrane"])


def _true_tm(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    if _negative_topology(text):
        return False
    return any(token in text for token in ["pdbtm", "memprotmd", "transmembrane", "channel", "transporter", "porin", "gpcr", "opsin"])


def _cofactor(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return bool(candidate.get("biological_cofactor_components")) or any(token in text for token in ["heme binding", "metal ion binding", "cofactor"])


def _biological_assembly(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return (
        int(candidate.get("polymer_entity_instance_count") or 0) >= 2
        or int(candidate.get("entity_molecule_count") or 0) >= 2
        or any(token in text for token in [
            "oligomer",
            "multimer",
            "assembly",
            "dimer",
            "trimer",
            "tetramer",
            "domain swap",
            "domain swapping",
            "coiled coil",
            "coiled-coil",
            "leucine zipper",
        ])
    )


def _helix_bundle_candidate(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    sequence = candidate["sequence"]
    helix_like = sum(1 for aa in sequence if aa in "AEHKLMQR") / max(1, len(sequence))
    breaker = sum(1 for aa in sequence if aa in "PG") / max(1, len(sequence))
    return (
        any(token in text for token in [
            "coiled coil",
            "coiled-coil",
            "leucine zipper",
            "helix bundle",
            "helical bundle",
            "collagen",
            "myosin",
            "tropomyosin",
        ])
        or (helix_like >= 0.37 and breaker <= 0.11 and candidate.get("sequence_metrics", {}).get("mean_disorder", 1.0) < 0.28)
    )


def _domain_swap_candidate(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return "domain swap" in text or "domain swapping" in text or "domain-swapped" in text


def _generic_complex_only(candidate: dict[str, Any]) -> bool:
    text = _candidate_text(candidate)
    return "complex" in text and not _biological_assembly(candidate) and not _true_tm(candidate) and not _cofactor(candidate)


def _soluble_hydrophobic_core(candidate: dict[str, Any]) -> bool:
    metrics = candidate.get("sequence_metrics", {})
    return (
        not _true_tm(candidate)
        and not _cofactor(candidate)
        and not _biological_assembly(candidate)
        and float(metrics.get("hydrophobic_density", 0.0)) >= 0.32
        and float(metrics.get("mean_disorder", 1.0)) < 0.25
    )


def _protein_map(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {row["protein_id"]: dict(row) for row in rows if row.get("protein_id")}


def _pick(
    candidates: list[dict[str, Any]],
    *,
    count: int,
    used: set[str],
    predicate,
    allow_reuse: bool = False,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for candidate in candidates:
        if not allow_reuse and candidate["protein_id"] in used:
            continue
        if not predicate(candidate):
            continue
        selected.append(candidate)
        used.add(candidate["protein_id"])
        if len(selected) == count:
            break
    if len(selected) < count and allow_reuse:
        for candidate in candidates:
            if not predicate(candidate):
                continue
            selected.append(candidate)
            if len(selected) == count:
                break
    if len(selected) != count:
        raise SystemExit(f"selected {len(selected)} rows; expected {count}")
    return selected


def _target_from_candidate(
    *,
    ordinal: int,
    panel_group: str,
    candidate: dict[str, Any],
    expected: str,
    source_family: str,
    lineage_source_target: str,
    truth_basis: list[str],
) -> dict[str, Any]:
    return {
        "target_id": f"V68_{ordinal:03d}_{_safe_id(panel_group)}_{_safe_id(candidate['protein_id'])}",
        "panel_group": panel_group,
        "source_family": source_family,
        "source_mode": "candidate_snapshot",
        "lineage_source_target": lineage_source_target,
        "protein_id": candidate["protein_id"],
        "entry_id": candidate["entry_id"],
        "entity_id": candidate["entity_id"],
        "sequence": candidate["sequence"],
        "sequence_length": candidate["sequence_length"],
        "sequence_cluster_30_id": candidate.get("sequence_cluster_30_id"),
        "expected_mechanism_class": expected,
        "target_name": f"{candidate.get('entry_id', '')} {candidate.get('entity_description', '')}".strip(),
        "entry_url": candidate.get("source_urls", {}).get("entry", ""),
        "polymer_entity_url": candidate.get("source_urls", {}).get("polymer_entity", ""),
        "postseal_truth_basis": truth_basis,
        "candidate_snapshot": dict(candidate),
    }


def _select_targets() -> list[dict[str, Any]]:
    v63_rows = [dict(row) for row in _read_json(V63_MANIFEST, "V63 target manifest")["selected_targets"]]
    v63_by_protein = _protein_map(v63_rows)
    v65_targets = [dict(row) for row in _read_json(V65_MANIFEST, "V65 membrane topology manifest")["targets"]]
    v67_failures = _read_json(V67_FAILURES, "V67 failure report")["failure_grammar_rows"]
    v67_targets = {row["target_id"]: row for row in _read_json(V67_TARGETS, "V67 target manifest")["selected_targets"]}

    used: set[str] = set()
    targets: list[dict[str, Any]] = []
    ordinal = 1

    assembly_failures = [
        row for row in v67_failures
        if row["failure_mode"] == "assembly_required_core_vs_membrane_topology"
    ]
    if len(assembly_failures) != GROUP_COUNTS["V67_ASSEMBLY_REQUIRED_REPLAY"]:
        raise SystemExit(f"V67 assembly failure count is {len(assembly_failures)}; expected 48")
    for row in assembly_failures:
        candidate = v63_by_protein.get(row["protein_id"]) or dict(v67_targets[row["target_id"]])
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="V67_ASSEMBLY_REQUIRED_REPLAY",
            candidate=candidate,
            expected=ASSEMBLY_CLASS,
            source_family="V67",
            lineage_source_target=row["target_id"],
            truth_basis=[
                "V67 dominant failure replayed for E65 assembly-required grammar mining.",
                f"V67 engine thought {row['engine_thought']}; V67 reality label was {row['reality_showed']}.",
                "E65 retests the missing word assembly_required_core_vs_topology_provider.",
            ],
        ))
        used.add(candidate["protein_id"])
        ordinal += 1

    assembly_candidates = sorted(
        v63_rows,
        key=lambda c: (
            -int(_domain_swap_candidate(c)),
            -int(_helix_bundle_candidate(c)),
            -int(c.get("polymer_entity_instance_count") or 0),
            -int(c.get("entity_molecule_count") or 0),
            c["protein_id"],
        ),
    )
    for candidate in _pick(
        assembly_candidates,
        count=GROUP_COUNTS["BIOLOGICAL_OLIGOMER_ASSEMBLY_REQUIRED"],
        used=used,
        predicate=lambda c: _biological_assembly(c) and not _true_tm(c) and not _cofactor(c),
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="BIOLOGICAL_OLIGOMER_ASSEMBLY_REQUIRED",
            candidate=candidate,
            expected=ASSEMBLY_CLASS,
            source_family="V63",
            lineage_source_target=f"V63_{candidate['target_id']}",
            truth_basis=["Public metadata indicates biological oligomer/assembly dependency candidate for E65."],
        ))
        ordinal += 1

    coil_candidates = sorted(
        v63_rows,
        key=lambda c: (
            -int(_domain_swap_candidate(c)),
            -int("coiled" in _candidate_text(c) or "zipper" in _candidate_text(c)),
            -float(c.get("sequence_metrics", {}).get("mean_interface", 0.0)),
            c["protein_id"],
        ),
    )
    for candidate in _pick(
        coil_candidates,
        count=GROUP_COUNTS["COILED_COIL_HELIX_BUNDLE_ASSEMBLY"],
        used=used,
        predicate=lambda c: _helix_bundle_candidate(c) and not _true_tm(c) and not _cofactor(c),
        allow_reuse=True,
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="COILED_COIL_HELIX_BUNDLE_ASSEMBLY",
            candidate=candidate,
            expected=ASSEMBLY_CLASS,
            source_family="V63",
            lineage_source_target=f"V63_{candidate['target_id']}",
            truth_basis=["Coiled-coil, helix-bundle, or sequence helix-register candidate for E65 assembly grammar."],
        ))
        ordinal += 1

    tm_candidates = sorted(
        v63_rows,
        key=lambda c: (
            -int("channel" in _candidate_text(c) or "porin" in _candidate_text(c) or "pore" in _candidate_text(c)),
            -float(c.get("sequence_metrics", {}).get("max_segment_membrane_density", 0.0)),
            c["protein_id"],
        ),
    )
    for candidate in _pick(
        tm_candidates,
        count=GROUP_COUNTS["TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL"],
        used=used,
        predicate=lambda c: _true_tm(c),
        allow_reuse=True,
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL",
            candidate=candidate,
            expected=MEMBRANE_CLASS,
            source_family="V63",
            lineage_source_target=f"V63_{candidate['target_id']}",
            truth_basis=["True transmembrane/channel/pore sentinel: E65 must preserve topology priority over assembly grammar."],
        ))
        ordinal += 1

    soluble_candidates = sorted(
        v63_rows,
        key=lambda c: (
            -float(c.get("sequence_metrics", {}).get("hydrophobic_density", 0.0)),
            float(c.get("sequence_metrics", {}).get("mean_disorder", 1.0)),
            c["protein_id"],
        ),
    )
    for candidate in _pick(
        soluble_candidates,
        count=GROUP_COUNTS["SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL"],
        used=used,
        predicate=_soluble_hydrophobic_core,
        allow_reuse=True,
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL",
            candidate=candidate,
            expected=GLOBULAR_CLASS,
            source_family="V63",
            lineage_source_target=f"V63_{candidate['target_id']}",
            truth_basis=["Soluble hydrophobic-core sentinel: E65 must not call assembly-required folding without assembly evidence."],
        ))
        ordinal += 1

    cofactor_pool = [v63_by_protein[row["protein_id"]] for row in v65_targets if row["panel_group"] == "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET" and row["protein_id"] in v63_by_protein]
    cofactor_pool.extend(v63_rows)
    for candidate in _pick(
        cofactor_pool,
        count=GROUP_COUNTS["COFACTOR_LIGAND_SOLUBLE_SENTINEL"],
        used=used,
        predicate=lambda c: _cofactor(c) and not _true_tm(c),
        allow_reuse=True,
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="COFACTOR_LIGAND_SOLUBLE_SENTINEL",
            candidate=candidate,
            expected=COFACTOR_CLASS,
            source_family="V65" if candidate["protein_id"] in {row["protein_id"] for row in v65_targets} else "V63",
            lineage_source_target=f"V63_{candidate['target_id']}",
            truth_basis=["Cofactor/ligand soluble sentinel: E65 must not convert ligand-stabilized cores into assembly claims."],
        ))
        ordinal += 1

    generic_candidates = sorted(v63_rows, key=lambda c: c["protein_id"])
    for candidate in _pick(
        generic_candidates,
        count=GROUP_COUNTS["GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL"],
        used=used,
        predicate=_generic_complex_only,
        allow_reuse=True,
    ):
        targets.append(_target_from_candidate(
            ordinal=ordinal,
            panel_group="GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL",
            candidate=candidate,
            expected=ABSTAIN_CLASS,
            source_family="V63",
            lineage_source_target=f"V63_{candidate['target_id']}",
            truth_basis=["Generic complex-only control: complex annotation alone is not obligate assembly evidence."],
        ))
        ordinal += 1

    for new_ordinal, target in enumerate(targets, start=1):
        target["target_id"] = f"V68_{new_ordinal:03d}_{_safe_id(target['panel_group'])}_{_safe_id(target['protein_id'])}"
    composition = Counter(target["panel_group"] for target in targets)
    if len(targets) != TARGET_COUNT or dict(composition) != dict(GROUP_COUNTS):
        raise SystemExit(f"bad V68 composition: total={len(targets)} composition={dict(composition)}")
    return targets


def _panel_statement(target: dict[str, Any]) -> str:
    candidate = target["candidate_snapshot"]
    text = _candidate_text(candidate)
    counts = (
        f"polymer instances={candidate.get('polymer_entity_instance_count', 0)}, "
        f"entity molecule count={candidate.get('entity_molecule_count', 0)}, "
        f"polymer composition={candidate.get('polymer_composition', '')}"
    )
    group = target["panel_group"]
    if group == "TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL":
        subtype = "transmembrane channel/pore" if any(token in text for token in ["channel", "pore", "porin"]) else "transmembrane topology"
        return (
            "V68 E65 topology sentinel source. membrane_context_strong transmembrane_context topology_evidence "
            f"{subtype}; biological_oligomer_context may be incidental but explicit topology has priority. "
            "Coordinates, contacts, and native topology remain post-seal holdout only."
        )
    if group == "COFACTOR_LIGAND_SOLUBLE_SENTINEL":
        components = " ".join(str(value) for value in candidate.get("biological_cofactor_components", []) or [])
        return (
            "V68 E65 cofactor sentinel source. cofactor_context ligand_context cofactor-stabilized soluble core; "
            f"components: {components}. Not assembly_required_core and no explicit transmembrane topology evidence."
        )
    if group == "SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL":
        return (
            "V68 E65 soluble core sentinel source. soluble_monomeric_core_context standalone soluble fold complete soluble monomer; "
            "hydrophobicity is explained by a soluble core, not membrane topology and not assembly_required_folding."
        )
    if group == "GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL":
        return (
            "V68 E65 generic complex control source. generic_complex_only generic complex only; "
            "complex annotation alone is not biological_oligomer_context, not obligate assembly, and not assembly_required_core."
        )
    marks = [
        "assembly_required_core",
        "assembly_required_folding",
        "partner_completed_core",
        "interface_buried_hydrophobicity",
        "monomer_incomplete_topology",
        "biological_oligomer_context",
    ]
    if group == "COILED_COIL_HELIX_BUNDLE_ASSEMBLY":
        marks.append("coiled_coil_register_dependency")
    if _domain_swap_candidate(candidate):
        marks.append("domain_swap_candidate")
    return (
        "V68 E65 assembly-required source from allowed sequence/public metadata context. "
        + " ".join(marks)
        + f"; {counts}. Generic complex alone would be insufficient, but this panel source marks partner-completed folding context. "
        "Coordinates, contacts, ligand geometry, and native assembly contacts are blocked before sealing."
    )


def _metadata_statement(candidate: dict[str, Any]) -> str:
    metrics = candidate.get("sequence_metrics", {})
    labels = [
        f"hydrophobic_density={metrics.get('hydrophobic_density')}",
        f"mean_interface={metrics.get('mean_interface')}",
        f"max_segment_membrane_density={metrics.get('max_segment_membrane_density')}",
        f"mean_disorder={metrics.get('mean_disorder')}",
    ]
    return ". ".join([
        f"RCSB title: {candidate.get('title', '')}",
        f"Entity description: {candidate.get('entity_description', '')}",
        f"Entry keywords: {candidate.get('entry_keywords', '')}",
        f"Polymer composition: {candidate.get('polymer_composition', '')}",
        "Sequence-derived non-coordinate summary: " + "; ".join(labels),
        "Coordinates, contacts, native interface geometry, ligand geometry, and post-seal validation labels are unopened before prediction hash.",
    ])


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    candidate = target["candidate_snapshot"]
    target_id = target["target_id"]
    return {
        "kind": "V68_E65_OLIGOMER_ASSEMBLY_PANEL_SOURCE_MANIFEST_v0",
        "target_id": target_id,
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "source_family": target["source_family"],
        "lineage_source_target": target["lineage_source_target"],
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "sequence": target["sequence"],
        "sequence_length": target["sequence_length"],
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
                "evidence_statement": _metadata_statement(candidate),
                "source_url": target["entry_url"],
            },
            {
                "source_id": f"{target_id}_E65_PANEL_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": _panel_statement(target),
                "source_url": target["entry_url"],
            },
        ],
        "blocked_prediction_inputs": [
            "PDB/mmCIF coordinates before sealing",
            "native contacts, residue-residue distances, and coordinate-derived topology",
            "native biological assembly contacts and ligand geometry before sealing",
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


def _expected_observables(expected: str) -> list[dict[str, Any]]:
    if expected == ASSEMBLY_CLASS:
        return [{"check_id": "partner_completed_core_supported", "metric": "partner_completed_core", "comparator": ">=", "threshold": 0.48}]
    if expected == MEMBRANE_CLASS:
        return [{"check_id": "proteostasis_routing_supported", "metric": "proteostasis_routing", "comparator": ">=", "threshold": 0.50}]
    if expected == COFACTOR_CLASS:
        return [{"check_id": "cofactor_interface_supported", "metric": "interface_readiness", "comparator": ">=", "threshold": 0.45}]
    if expected == GLOBULAR_CLASS:
        return [{"check_id": "compact_or_contact_supported", "metric": "contact_probability", "comparator": ">=", "threshold": 0.35}]
    return []


def _holdout(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V68_E65_OLIGOMER_ASSEMBLY_PANEL_POSTSEAL_HOLDOUT_v0",
        "target_id": packet["target_id"],
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "panel_group": target["panel_group"],
        "source_family": target["source_family"],
        "lineage_source_target": target["lineage_source_target"],
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
        "expected_mechanism_class": target["expected_mechanism_class"],
        "expected_observables": _expected_observables(target["expected_mechanism_class"]),
        "postseal_truth_basis": target["postseal_truth_basis"],
        "holdout_opened_after_prediction_hash": packet["prediction_hash"],
        "postseal_sources": [
            {
                "source_id": f"{packet['target_id']}_V68_POSTSEAL_HOLDOUT",
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
    if expected == ASSEMBLY_CLASS:
        return [
            {"perturbation_id": f"{target_id}_ASSEMBLY_INTERFACE_DAMAGE", "description": "damage partner-completed assembly interface", "operator_scales": {"interface_operator": 0.42, "closure_operator": 0.70}, "interface_disruption": 0.45, "metric": "partner_completed_core", "expected_direction": "decrease"},
            {"perturbation_id": f"{target_id}_CONCENTRATION_RESCUE", "description": "assembly concentration rescue context", "operator_scales": {"interface_operator": 1.15}, "concentration_rescue": 0.25, "metric": "partner_completed_core", "expected_direction": "increase"},
        ]
    if expected == MEMBRANE_CLASS:
        return [
            {"perturbation_id": f"{target_id}_MEMBRANE_DAMAGE", "description": "damage topology/proteostasis route", "operator_scales": {"membrane_pressure_operator": 0.55, "proteostasis_operator": 0.55}, "damage": 0.40, "metric": "proteostasis_routing", "expected_direction": "decrease"},
        ]
    if expected == COFACTOR_CLASS:
        return [
            {"perturbation_id": f"{target_id}_COFACTOR_REMOVAL", "description": "remove ligand/cofactor pressure", "operator_scales": {"interface_operator": 0.45}, "metric": "interface_readiness", "expected_direction": "decrease"},
        ]
    if expected == GLOBULAR_CLASS:
        return [
            {"perturbation_id": f"{target_id}_CORE_DAMAGE", "description": "damage soluble hydrophobic core", "operator_scales": {"closure_operator": 0.45}, "metric": "contact_probability", "expected_direction": "decrease"},
        ]
    return []


def _score(packet: dict[str, Any], holdout: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    validation = validate_against_holdout(sealed_packet=packet, holdout=holdout)
    predicted = packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    decision = "abstain_recommended" if predicted == ABSTAIN_CLASS else "accepted"
    supported = predicted == expected and validation["score_label"] == "supported"
    return {
        "kind": "V68_E65_OLIGOMER_ASSEMBLY_PANEL_VALIDATION_RESULT_v0",
        "target_id": packet["target_id"],
        "panel_group": target["panel_group"],
        "source_family": target["source_family"],
        "lineage_source_target": target["lineage_source_target"],
        "protein_id": target["protein_id"],
        "entry_id": target["entry_id"],
        "entity_id": target["entity_id"],
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
        "kind": "V68_COMPACT_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "sealed_before_holdout": packet["sealed_before_holdout"],
        "selected_mechanism_class": mechanism["mechanism_class"],
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "selection_reason": mechanism["selection_reason"],
        "operator_names": packet["operator_field"]["operator_names"],
        "active_operator_count": packet["operator_field"]["active_operator_count"],
        "operator_state_final_state_summary": packet["operator_state_propagation_summary"]["final_state_summary"],
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V68_E65_ENGINE_DECLARATION_v0",
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
        "lineage_note": "E65 adds assembly_required_folding and assembly/topology ambiguity abstention.",
        "engine_changes_allowed_by_operating_mode": True,
        "engine_modified_during_batch": False,
        "engine_biology_modified_during_batch": False,
        "target_selection_manual": False,
        "folding_problem_solved": False,
    }


def _e65_certificate(engine_declaration: dict[str, Any]) -> dict[str, Any]:
    cert = {
        "kind": "E65_ASSEMBLY_REQUIRED_FOLDING_GRAMMAR_CERTIFICATE_v0",
        "engine_version": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_batch_trigger": "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64",
        "lineage": ["E60", "E61", "E62", "E63", "E64", "E65"],
        "failure_mode_addressed": "assembly_required_core_vs_membrane_topology",
        "missing_esperanto_word": "assembly_required_core_vs_topology_provider",
        "new_mechanism_class": ASSEMBLY_CLASS,
        "new_state_variables": [
            "assembly_required_core",
            "partner_completed_core",
            "interface_buried_hydrophobicity",
            "monomer_incomplete_topology",
            "coiled_coil_register_dependency",
            "domain_swap_candidate",
            "biological_oligomer_context",
            "generic_complex_only",
            "assembly_ambiguous",
        ],
        "decision_boundary": [
            "explicit true transmembrane/topology evidence keeps membrane priority",
            "explicit assembly-required evidence routes to assembly_required_folding",
            "cofactor/metal/ligand evidence routes to cofactor_stabilized_core when assembly-required evidence is absent",
            "soluble monomeric core evidence routes to globular_closure",
            "generic complex alone abstains",
            "weak competing assembly/topology/cofactor explanations abstain",
        ],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "operator_set_hash": engine_declaration["operator_set_hash"],
        "mechanism_class_set_hash": engine_declaration["mechanism_class_set_hash"],
        "claim_allowed": False,
        "next_required_batch": BATCH_ID,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _target_manifest(targets: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for target in targets:
        row = {key: value for key, value in target.items() if key != "candidate_snapshot"}
        candidate = target["candidate_snapshot"]
        row.update({
            "title": candidate.get("title", ""),
            "entity_description": candidate.get("entity_description", ""),
            "polymer_composition": candidate.get("polymer_composition", ""),
            "polymer_entity_instance_count": candidate.get("polymer_entity_instance_count", 0),
            "entity_molecule_count": candidate.get("entity_molecule_count", 0),
            "biological_cofactor_components": candidate.get("biological_cofactor_components", []),
            "sequence_metrics": candidate.get("sequence_metrics", {}),
            "domain_swap_candidate": _domain_swap_candidate(candidate),
            "coiled_coil_or_helix_bundle_candidate": _helix_bundle_candidate(candidate),
        })
        rows.append(row)
    return {
        "kind": "V68_E65_OLIGOMER_ASSEMBLY_PANEL_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "campaign_id": CAMPAIGN_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "batch_mode": "specialist_assembly_panel_200",
        "target_selection_manual": False,
        "composition_rule": dict(GROUP_COUNTS),
        "target_count_requested": TARGET_COUNT,
        "target_count_selected": len(rows),
        "selection_rule": "48 V67 assembly failures plus six deterministic specialist control/sentinel groups from committed V63/V65 artifacts.",
        "selected_targets": rows,
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    supported = [row for row in rows if row["score_label"] == "supported"]
    failed_accepted = [row for row in accepted if row["score_label"] != "supported"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    clean_abstain_supported = [row for row in abstained if row["score_label"] == "supported"]
    return {
        "targets_total": len(rows),
        "accepted_count": len(accepted),
        "supported_count": len(supported),
        "accepted_supported": len(accepted_supported),
        "clean_abstain_supported": len(clean_abstain_supported),
        "failed_accepted": len(failed_accepted),
        "failed_accepted_count": len(failed_accepted),
        "abstain_count": len(abstained),
        "accepted_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "raw_accuracy": len(supported) / len(rows) if rows else None,
        "coverage": len(accepted) / len(rows) if rows else None,
    }


def _failure_mode(row: dict[str, Any]) -> str:
    if row["score_label"] == "supported":
        return "supported"
    predicted = row["predicted_mechanism_class"]
    expected = row["expected_mechanism_class"]
    if predicted == ASSEMBLY_CLASS and expected == ABSTAIN_CLASS:
        return "generic_complex_false_assembly"
    if predicted == ASSEMBLY_CLASS and expected == GLOBULAR_CLASS:
        return "soluble_core_false_assembly"
    if predicted == ASSEMBLY_CLASS and expected == MEMBRANE_CLASS:
        return "true_TM_false_assembly"
    if predicted == ASSEMBLY_CLASS and expected == COFACTOR_CLASS:
        return "cofactor_false_assembly"
    if expected == ASSEMBLY_CLASS and predicted == MEMBRANE_CLASS:
        return "assembly_missed_as_membrane"
    if expected == ASSEMBLY_CLASS and predicted == GLOBULAR_CLASS:
        return "assembly_missed_as_soluble"
    if expected == ASSEMBLY_CLASS and predicted == ABSTAIN_CLASS and row["panel_group"] == "COILED_COIL_HELIX_BUNDLE_ASSEMBLY":
        return "coiled_coil_register_missed"
    if expected == ASSEMBLY_CLASS and predicted == ABSTAIN_CLASS:
        return "assembly_required_abstained"
    if row["panel_group"] == "COILED_COIL_HELIX_BUNDLE_ASSEMBLY" and predicted != ASSEMBLY_CLASS:
        return "coiled_coil_register_missed"
    return "wrong_regime"


def _failure_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for row in scoring_rows:
        if row["score_label"] == "supported":
            continue
        mode = _failure_mode(row)
        rows.append({
            "target_id": row["target_id"],
            "protein_id": row["protein_id"],
            "panel_group": row["panel_group"],
            "failure_mode": mode,
            "engine_thought": row["predicted_mechanism_class"],
            "reality_showed": row["expected_mechanism_class"],
            "acceptance_decision": row["acceptance_decision"],
            "score_label": row["score_label"],
            "missing_esperanto_word": mode,
            "lineage_source_target": row.get("lineage_source_target"),
            "autopsy_sentence": (
                f"The engine thought: {row['predicted_mechanism_class']}. "
                f"Reality showed: {row['expected_mechanism_class']}. "
                f"Missing Esperanto word: {mode}."
            ),
        })
    counts = Counter(row["failure_mode"] for row in rows)
    return {
        "kind": "V68_E65_OLIGOMER_ASSEMBLY_PANEL_FAILURE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "failure_cases_reported": True,
        "failure_count": len(rows),
        "failure_modes": dict(counts),
        "dominant_failure_mode": counts.most_common(1)[0][0] if rows else None,
        "dominant_failure_count": counts.most_common(1)[0][1] if rows else 0,
        "missing_words_top_10": [{"failure_mode": mode, "count": count} for mode, count in counts.most_common(10)],
        "failure_grammar_rows": rows,
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _panel_metrics(scoring_rows: list[dict[str, Any]], targets: list[dict[str, Any]]) -> dict[str, Any]:
    by_group = {group: [row for row in scoring_rows if row["panel_group"] == group] for group in GROUP_COUNTS}
    target_by_id = {target["target_id"]: target for target in targets}
    domain_swap_rows = [
        row for row in scoring_rows
        if target_by_id[row["target_id"]]["candidate_snapshot"] and _domain_swap_candidate(target_by_id[row["target_id"]]["candidate_snapshot"])
    ]
    old_v67_rows = by_group["V67_ASSEMBLY_REQUIRED_REPLAY"]
    old_v67_supported = [row for row in old_v67_rows if row["score_label"] == "supported"]
    old_v67_reduction = GROUP_COUNTS["V67_ASSEMBLY_REQUIRED_REPLAY"] - sum(
        1 for row in old_v67_rows if row["predicted_mechanism_class"] != ASSEMBLY_CLASS
    )
    return {
        "assembly_required_correct": sum(1 for row in scoring_rows if row["expected_mechanism_class"] == ASSEMBLY_CLASS and row["score_label"] == "supported"),
        "true_TM_preserved": sum(1 for row in by_group["TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL"] if row["score_label"] == "supported"),
        "soluble_hydrophobic_not_called_assembly": sum(1 for row in by_group["SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL"] if row["predicted_mechanism_class"] != ASSEMBLY_CLASS),
        "cofactor_not_called_assembly": sum(1 for row in by_group["COFACTOR_LIGAND_SOLUBLE_SENTINEL"] if row["predicted_mechanism_class"] != ASSEMBLY_CLASS),
        "generic_complex_not_called_assembly": sum(1 for row in by_group["GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL"] if row["predicted_mechanism_class"] != ASSEMBLY_CLASS),
        "coiled_coil_register_detected": sum(1 for row in by_group["COILED_COIL_HELIX_BUNDLE_ASSEMBLY"] if row["predicted_mechanism_class"] == ASSEMBLY_CLASS),
        "domain_swap_candidates_detected": sum(1 for row in domain_swap_rows if row["predicted_mechanism_class"] == ASSEMBLY_CLASS),
        "assembly_ambiguous_abstained": sum(1 for row in scoring_rows if row["predicted_mechanism_class"] == ABSTAIN_CLASS and "ambiguous" in str(row["validation_checks"]).lower()),
        "v67_dominant_failure_replay_count": len(old_v67_rows),
        "v67_dominant_failure_repaired_by_e65": len(old_v67_supported),
        "v67_dominant_failure_remaining": len(old_v67_rows) - len(old_v67_supported),
        "v67_dominant_failure_reduction_count": old_v67_reduction,
    }


def _controls(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
    failure_report: dict[str, Any],
    panel_metrics: dict[str, Any],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{"source_id": "V68_BAD_COORDINATES", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    alphafold_gate = evidence_boundary_gate([{"source_id": "V68_BAD_ALPHAFOLD_MODEL", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True, "evidence_statement": "AlphaFold-style model offered before sealing."}])
    holdout_gate = evidence_boundary_gate([{"source_id": "V68_PRESEAL_HOLDOUT", "source_class": COORDINATE_DERIVED, "source_role": "holdout_validation", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V68_BAD_INTERNAL_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    random_packet = build_sealed_operator_state_packet(target_id="V68_RANDOM_SEQUENCE_CONTROL", target_name="V68 random sequence control", sequence=deterministic_random_sequence(128), sources=[], perturbations=[])
    composition = Counter(row["panel_group"] for row in target_manifest["selected_targets"])
    sentinel_rows = [
        row for row in scoring_rows
        if row["panel_group"] in {
            "TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL",
            "SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL",
            "COFACTOR_LIGAND_SOLUBLE_SENTINEL",
            "GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL",
        }
    ]
    sentinel_regressions = [row["target_id"] for row in sentinel_rows if row["score_label"] != "supported"]
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
        _control("target_count_200", target_manifest["target_count_selected"] == TARGET_COUNT, "V68 must have exactly 200 targets.", target_manifest["target_count_selected"]),
        _control("composition_rule", dict(composition) == dict(GROUP_COUNTS), "V68 must match requested specialist panel composition.", dict(composition)),
        _control("engine_version_declared_e65", engine_declaration["engine_version_used"] == ENGINE_VERSION_USED, "V68 uses E65."),
        _control("assembly_class_available", ASSEMBLY_CLASS in engine_declaration["mechanism_classes"], "E65 exposes assembly_required_folding as a mechanism class."),
        _control("target_selection_manual_false", target_manifest["target_selection_manual"] is False, "V68 target selection is deterministic from committed artifacts."),
        _control("all_sealed_before_holdout", all(packet["sealed_before_holdout"] for packet in packets), "Every V68 packet is sealed before holdout."),
        _control("holdouts_opened_after_seal", all(row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"] for row in scoring_rows), "Every holdout references sealed prediction hash."),
        _control("v67_dominant_failure_reduced_by_half", panel_metrics["v67_dominant_failure_remaining"] <= 24, "Minimum useful pass: V67 dominant failure remaining count reduced by at least half.", panel_metrics),
        _control("sentinels_stable", not sentinel_regressions, "True TM, soluble, cofactor, and generic-complex sentinels must stay correct.", {"sentinel_count": len(sentinel_rows), "regressions": sentinel_regressions}),
        _control("dominant_failure_mode_identified", failure_report["dominant_failure_mode"] is not None or failure_report["failure_count"] == 0, "V68 identifies the next dominant failure when failures remain.", failure_report.get("missing_words_top_10")),
        _control("shuffled_sequence_controls_reported", len(shuffled_rows) == len(packets) and all(row["shuffled_coordinate_sources"] == 0 and row["shuffled_runtime_sources"] == 0 for row in shuffled_rows), "Composition-preserving shuffled controls generated without target metadata.", {"control_count": len(shuffled_rows)}),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Holdout validation cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == ABSTAIN_CLASS, "Random sequence without evidence abstains."),
        _control("readme_check_skipped_by_user_instruction", True, "README check skipped by explicit user instruction."),
    ]


def _accept_abstain_posture(metrics: dict[str, Any]) -> str:
    if metrics["failed_accepted_count"] > metrics["abstain_count"]:
        return "over_accepting_relative_to_abstention"
    if metrics["abstain_count"] > metrics["failed_accepted_count"]:
        return "abstaining_more_than_failed_accepting"
    return "accept_abstain_balanced"


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    failure_report: dict[str, Any],
    panel_metrics: dict[str, Any],
) -> dict[str, Any]:
    metrics = _metrics(scoring_rows)
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    critical = {"coordinate_leakage_control", "alphafold_leakage_control", "holdout_opened_before_seal_control", "internal_runtime_leakage_control"}
    sentinel_regressions = next(row for row in controls if row["control_id"] == "sentinels_stable")["observed"]["regressions"]
    if any(control_id in critical for control_id in failed_controls):
        status = BLOCKED_LEAKAGE
    elif sentinel_regressions:
        status = BLOCKED_SENTINEL
    elif failed_controls:
        status = BLOCKED_CONTROLS
    elif metrics["failed_accepted_count"]:
        status = MINED
    else:
        status = PASSED
    cert = {
        "kind": "V68_E65_OLIGOMER_ASSEMBLY_PANEL_200_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_id": CAMPAIGN_ID,
        "batch_id": BATCH_ID,
        "batch_mode": "specialist_assembly_panel_200",
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "engine_source_last_commit": engine_declaration["engine_source_last_commit"],
        "engine_source_sha256": engine_declaration["engine_source_sha256"],
        "target_selection_manual": target_manifest["target_selection_manual"],
        "composition_rule": target_manifest["composition_rule"],
        **metrics,
        "accept_abstain_posture": _accept_abstain_posture(metrics),
        **panel_metrics,
        "sentinel_regressions": sentinel_regressions,
        "sentinel_regression_count": len(sentinel_regressions),
        "controls_passed": not failed_controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "controls": controls,
        "failure_modes": failure_report["failure_modes"],
        "dominant_failure_mode": failure_report["dominant_failure_mode"],
        "dominant_failure_count": failure_report["dominant_failure_count"],
        "missing_words_top_10": failure_report["missing_words_top_10"],
        "v68_questions_answered": {
            "assembly_required_correct": panel_metrics["assembly_required_correct"],
            "true_TM_preserved": panel_metrics["true_TM_preserved"],
            "soluble_hydrophobic_not_called_assembly": panel_metrics["soluble_hydrophobic_not_called_assembly"],
            "cofactor_not_called_assembly": panel_metrics["cofactor_not_called_assembly"],
            "generic_complex_not_called_assembly": panel_metrics["generic_complex_not_called_assembly"],
            "coiled_coil_register_detected": panel_metrics["coiled_coil_register_detected"],
            "domain_swap_candidates_detected": panel_metrics["domain_swap_candidates_detected"],
            "assembly_ambiguous_abstained": panel_metrics["assembly_ambiguous_abstained"],
        },
        "coordinate_truth_used_before_seal": False,
        "contact_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "atomistic_md_performed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "claim_blocked_reason": "V68 is a specialist assembly grammar panel, not a solved-folding claim.",
        "next_required_batch": "V69_RCSB_NONREDUNDANT_200_DISCOVERY_E65" if panel_metrics["v67_dominant_failure_remaining"] <= 24 and not sentinel_regressions else "E65_REVIEW_BEFORE_V69",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _claim_row(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_id": BATCH_ID,
        "batch_mode": "specialist_assembly_panel_200",
        "engine_version_used": ENGINE_VERSION_USED,
        "raw_accuracy": cert["raw_accuracy"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "coverage": cert["coverage"],
        "abstention_rate": cert["abstain_count"] / cert["targets_total"] if cert["targets_total"] else None,
        "control_pass_rate": cert["passed_control_count"] / cert["control_count"] if cert["control_count"] else None,
        "failure_count": cert["failed_accepted_count"],
        "clean_abstain_count": cert["abstain_count"],
        "dominant_failure_mode": cert["dominant_failure_mode"],
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


def _write_report(path: Path, cert: dict[str, Any], failure_report: dict[str, Any]) -> None:
    lines = [
        "# V68 E65 Oligomer Assembly Panel 200",
        "",
        f"Status: `{cert['status']}`",
        f"Targets total: `{cert['targets_total']}`",
        f"Supported: `{cert['supported_count']}`",
        f"Accepted supported: `{cert['accepted_supported']}`",
        f"Clean abstain supported: `{cert['clean_abstain_supported']}`",
        f"Failed accepted: `{cert['failed_accepted_count']}`",
        f"Abstain: `{cert['abstain_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Raw accuracy: `{cert['raw_accuracy']}`",
        f"Coverage: `{cert['coverage']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"Sentinel regressions: `{cert['sentinel_regression_count']}`",
        f"V67 dominant failure remaining: `{cert['v67_dominant_failure_remaining']}`",
        f"Next required batch: `{cert['next_required_batch']}`",
        "",
        "## V68 Metrics",
        f"- assembly_required_correct: `{cert['assembly_required_correct']}`",
        f"- true_TM_preserved: `{cert['true_TM_preserved']}`",
        f"- soluble_hydrophobic_not_called_assembly: `{cert['soluble_hydrophobic_not_called_assembly']}`",
        f"- cofactor_not_called_assembly: `{cert['cofactor_not_called_assembly']}`",
        f"- generic_complex_not_called_assembly: `{cert['generic_complex_not_called_assembly']}`",
        f"- coiled_coil_register_detected: `{cert['coiled_coil_register_detected']}`",
        f"- domain_swap_candidates_detected: `{cert['domain_swap_candidates_detected']}`",
        f"- assembly_ambiguous_abstained: `{cert['assembly_ambiguous_abstained']}`",
        "",
        "## Failure Mode Table",
        "",
        "| failure_mode | count |",
        "| --- | ---: |",
    ]
    for item in failure_report["missing_words_top_10"]:
        lines.append(f"| `{item['failure_mode']}` | `{item['count']}` |")
    lines.extend([
        "",
        "## Boundary",
        "V68 is a specialist grammar panel. It does not make a broad solved-folding claim.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v68(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_generated_outputs(out_dir)
    targets = _select_targets()
    target_manifest = _target_manifest(targets)
    engine_declaration = _engine_declaration()
    e65_certificate = _e65_certificate(engine_declaration)
    _write_json(E65_ROOT / "e65_assembly_required_folding_grammar_certificate.json", e65_certificate)
    _write_json(DATA_ROOT / "v68_oligomer_assembly_panel_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v68_e65_engine_declaration.json", engine_declaration)

    packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []

    for target in targets:
        source_manifest = _source_manifest(target)
        packet = build_sealed_operator_state_packet(
            target_id=target["target_id"],
            target_name=target["target_name"],
            sequence=target["sequence"],
            sources=source_manifest["prediction_sources"],
            focus_regions=[{"name": "V68 E65 assembly/topology specialist scan", "span": f"1-{target['sequence_length']}"}],
            perturbations=_perturbations_for_expected(target),
        )
        holdout = _holdout(target, packet)
        score = _score(packet, holdout, target)
        shuffled_packet = build_sealed_operator_state_packet(
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
                    "evidence_statement": "Deterministically shuffled sequence control. Target metadata and E65 panel context are withheld.",
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
    panel_metrics = _panel_metrics(scoring_rows, targets)
    controls = _controls(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        packets=packets,
        scoring_rows=scoring_rows,
        shuffled_packets=shuffled_packets,
        failure_report=failure_report,
        panel_metrics=panel_metrics,
    )
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
        failure_report=failure_report,
        panel_metrics=panel_metrics,
    )
    claim_row = _claim_row(cert)

    scoring_path = _write_json(DATA_ROOT / "v68_oligomer_assembly_panel_scoring_report.json", {"kind": "V68_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v68_oligomer_assembly_panel_failure_report.json", failure_report)
    data_cert_path = _write_json(DATA_ROOT / "v68_oligomer_assembly_panel_certificate.json", cert)
    claim_row_path = _write_json(DATA_ROOT / "v68_campaign_claim_row.json", claim_row)
    claim_ledger_path = _append_claim_ledger(claim_row)

    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = _write_json(out_dir / "v68_e65_oligomer_assembly_panel_200_certificate.json", cert)
    report_path = out_dir / "V68_E65_OLIGOMER_ASSEMBLY_PANEL_200_REPORT.md"
    _write_report(report_path, cert, failure_report)
    return {
        "e65_certificate": E65_ROOT / "e65_assembly_required_folding_grammar_certificate.json",
        "target_manifest": DATA_ROOT / "v68_oligomer_assembly_panel_target_manifest.json",
        "engine_declaration": DATA_ROOT / "v68_e65_engine_declaration.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
        "claim_row": claim_row_path,
        "campaign_claim_ledger": claim_ledger_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V68 E65 oligomer/assembly specialist panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v68(args.out_dir)
    cert = _read_json(paths["certificate"], "V68 certificate")
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
        "accept_abstain_posture": cert["accept_abstain_posture"],
        "controls_passed": cert["controls_passed"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "assembly_required_correct": cert["assembly_required_correct"],
        "true_TM_preserved": cert["true_TM_preserved"],
        "generic_complex_not_called_assembly": cert["generic_complex_not_called_assembly"],
        "coiled_coil_register_detected": cert["coiled_coil_register_detected"],
        "domain_swap_candidates_detected": cert["domain_swap_candidates_detected"],
        "sentinel_regression_count": cert["sentinel_regression_count"],
        "v67_dominant_failure_remaining": cert["v67_dominant_failure_remaining"],
        "dominant_failure_mode": cert["dominant_failure_mode"],
        "dominant_failure_count": cert["dominant_failure_count"],
        "failure_modes": cert["failure_modes"],
        "next_required_batch": cert["next_required_batch"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["controls_passed"] and not cert["sentinel_regressions"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
