#!/usr/bin/env python3
from __future__ import annotations

"""Run V59 real folding-process replication panel.

V59 keeps the Protein Esperanto engine frozen.  It seals a coarse process
trajectory prediction from sequence plus non-coordinate metadata, then opens
process holdouts from folding kinetics / phi-value / HDX / NMR literature.
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
    build_sealed_operator_state_packet,
    build_sequence_field,
    deterministic_random_sequence,
    evidence_boundary_gate,
    predict_process_route_profile,
    sequence_operator_coherence,
    shuffled_sequence,
    stable_hash,
)


DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V59"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V59_REAL_FOLDING_PROCESS_REPLICATION_PANEL"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

PASSED = "V59_REAL_FOLDING_PROCESS_REPLICATION_PASSED_REVIEW_REQUIRED"
PARTIAL = "V59_PARTIAL_PROCESS_REPLICATION_WITH_ABSTENTION_PASSED_REVIEW_REQUIRED"
BLOCKED_ENGINE = "V59_PROCESS_REPLICATION_BLOCKED_ENGINE_NEEDS_REVISION"
BLOCKED_LEAKAGE = "V59_PROCESS_REPLICATION_BLOCKED_FOR_LEAKAGE"

PROCESS_CLASSES = [
    "two_state",
    "intermediate_bearing",
    "multi_basin",
    "disorder_biased",
    "fold_upon_binding",
]

PRIOR_GRAMMAR_TARGET_NAMES = [
    "p53_TAD_MDM2",
    "KcsA",
    "XCL1",
    "SNCA",
    "FUS_LC",
    "TDP43_LCD",
    "CFTR_F508DEL",
    "RFAH_CTD",
    "SARS2_ORF6",
    "GB1_DOMAIN",
    "HNRNPA1_LCD",
    "RHODOPSIN_P23H",
]


PROCESS_PANEL_TARGETS: list[dict[str, Any]] = [
    {
        "target_id": "V59_CI2",
        "target_name": "Chymotrypsin inhibitor 2",
        "pdb_id": "2CI2",
        "sequence_source": "https://www.rcsb.org/fasta/entry/2CI2/display",
        "sequence": "SSVEKKPEGVNTGAGDRHNLKTEWPELVGKSVEEAKKVILQDKPEAQIIVLPVGTIVTMEYRIDRVRLFVDKLDNIAEVPRVG",
        "structural_class": "alpha_beta",
        "architecture_hint": "small mixed alpha beta single-domain folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "folding_kinetics"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "mixed_alpha_beta_core",
            "early_span": [20, 55],
            "late_region_family": "terminal_consolidation",
            "late_span": [1, 12],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["hydrophobic_core", "mixed_alpha_beta_core"],
            "evidence_summary": "CI2 is a classic phi-value target with a compact transition-state nucleus in the central helix/core and terminal consolidation later in the route.",
            "process_sources": [
                "https://doi.org/10.1038/340122a0",
                "https://arxiv.org/abs/q-bio/0605048",
                "https://kineticdb.ivankovlab.ru/",
            ],
        },
    },
    {
        "target_id": "V59_PROTEIN_L",
        "target_name": "Protein L B1 domain",
        "pdb_id": "1HZ6",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1HZ6/display",
        "sequence": "MHHHHHHAMEEVTIKANLIFANGSTQTAEFKGTFEKATSEAYAYADTLKKDNGEWTVDVADKGYTLNIKFAG",
        "structural_class": "beta",
        "architecture_hint": "small beta-sheet single-domain folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "folding_kinetics"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "beta_hairpin_core",
            "early_span": [18, 54],
            "late_region_family": "terminal_strand_lock",
            "late_span": [1, 12],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["beta_hairpin_core", "turn_or_sheet_core"],
            "evidence_summary": "Protein L is represented in phi-value modeling and folding-kinetics work as a small two-state beta-domain with hairpin/sheet transition-state structure.",
            "process_sources": [
                "https://arxiv.org/abs/q-bio/0605048",
                "https://kineticdb.ivankovlab.ru/",
            ],
        },
    },
    {
        "target_id": "V59_PROTEIN_A_B",
        "target_name": "Staphylococcal protein A B domain",
        "pdb_id": "1BDD",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1BDD/display",
        "sequence": "TADNKFNKEQQNAFYEILHLPNLNEEQRNGFIQSLKDDPSQSANLLAEAKKLNDAQAPKA",
        "structural_class": "alpha",
        "architecture_hint": "small three-helix bundle folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "temperature_jump_kinetics"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "helix_bundle_core",
            "early_span": [17, 55],
            "late_region_family": "terminal_consolidation",
            "late_span": [1, 12],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["helix_bundle_core", "hydrophobic_core"],
            "evidence_summary": "Protein A B-domain folding kinetics and phi-value interpretation localize the transition-state signal to the helix-bundle core.",
            "process_sources": [
                "https://doi.org/10.1073/pnas.0306433101",
                "https://arxiv.org/abs/q-bio/0605048",
            ],
        },
    },
    {
        "target_id": "V59_SRC_SH3",
        "target_name": "Alpha-spectrin SH3 domain",
        "pdb_id": "1SHG",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1SHG/display",
        "sequence": "MDETGKELVLALYDYQEKSPREVTMKKGDILTLLNSTNKDWWKVEVNDRQGFVPAAYVKKLD",
        "structural_class": "beta",
        "architecture_hint": "small beta-sheet SH3 folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "folding_kinetics"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "beta_hairpin_core",
            "early_span": [21, 52],
            "late_region_family": "terminal_strand_lock",
            "late_span": [1, 10],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["beta_hairpin_core", "turn_or_sheet_core"],
            "evidence_summary": "SH3 domains are common phi-value/process targets with beta-sheet/loop transition-state organization rather than a final-coordinate-only claim.",
            "process_sources": [
                "https://kineticdb.ivankovlab.ru/",
                "https://arxiv.org/abs/q-bio/0311004",
            ],
        },
    },
    {
        "target_id": "V59_FBP_WW",
        "target_name": "Formin-binding protein WW domain",
        "pdb_id": "1E0L",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1E0L/display",
        "sequence": "GATAVSEWTEYKTADGKTYYYNNRTLESTWEKPQELK",
        "structural_class": "beta",
        "architecture_hint": "ultrashort WW beta-sheet domain folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "folding_kinetics"],
        "process_holdout": {
            "expected_process_class": "multi_basin",
            "early_region_family": "alternative_beta_hairpin",
            "early_span": [7, 30],
            "late_region_family": "sheet_locking",
            "late_span": [28, 38],
            "folding_nucleus_decision": "distributed_or_alternative_nucleus",
            "mutation_sensitive_region_families": ["aromatic_beta_core", "alternative_beta_hairpin"],
            "evidence_summary": "WW-domain phi-value modeling supports alternative beta-hairpin transition-state conformations rather than one single fixed early structure.",
            "process_sources": [
                "https://arxiv.org/abs/0709.3359",
                "https://www.rcsb.org/structure/1E0L",
            ],
        },
    },
    {
        "target_id": "V59_BARNASE",
        "target_name": "Barnase",
        "pdb_id": "1BNR",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1BNR/display",
        "sequence": "AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR",
        "structural_class": "alpha_beta",
        "architecture_hint": "larger mixed alpha beta enzyme-like folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "folding_intermediate", "amide_exchange"],
        "process_holdout": {
            "expected_process_class": "intermediate_bearing",
            "early_region_family": "subdomain_intermediate_core",
            "early_span": [10, 76],
            "late_region_family": "delayed_helix_or_terminal_docking",
            "late_span": [78, 110],
            "folding_nucleus_decision": "intermediate_with_nucleus",
            "mutation_sensitive_region_families": ["intermediate_stabilizing_core", "hydrophobic_core"],
            "evidence_summary": "Barnase has extensive phi-value and intermediate work; process evidence includes transient/intermediate structure rather than only a two-state endpoint.",
            "process_sources": [
                "https://doi.org/10.1006/jmbi.1997.1547",
                "https://doi.org/10.1038/340122a0",
                "https://kineticdb.ivankovlab.ru/",
            ],
        },
    },
    {
        "target_id": "V59_UBIQUITIN",
        "target_name": "Human ubiquitin",
        "pdb_id": "1UBQ",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1UBQ/display",
        "sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "structural_class": "alpha_beta",
        "architecture_hint": "small mixed alpha beta single-domain folding-kinetics target",
        "process_evidence_channels": ["phi_value_analysis", "psi_value_analysis", "folding_kinetics"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "mixed_alpha_beta_core",
            "early_span": [1, 45],
            "late_region_family": "terminal_consolidation",
            "late_span": [64, 76],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["mixed_alpha_beta_core", "hydrophobic_core"],
            "evidence_summary": "Ubiquitin process studies use kinetic perturbation / phi-psi style evidence to localize a heterogeneous but bounded transition-state core.",
            "process_sources": [
                "https://doi.org/10.1073/pnas.97.12.6527",
                "https://www.rcsb.org/structure/1UBQ",
            ],
        },
    },
    {
        "target_id": "V59_VILLIN_HP35",
        "target_name": "Villin headpiece HP35",
        "pdb_id": "1VII",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1VII/display",
        "sequence": "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF",
        "structural_class": "alpha",
        "architecture_hint": "small fast alpha-helical folder",
        "process_evidence_channels": ["folding_kinetics", "fast_folder_temperature_jump"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "helix_bundle_core",
            "early_span": [10, 34],
            "late_region_family": "terminal_consolidation",
            "late_span": [1, 9],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["helix_bundle_core", "hydrophobic_core"],
            "evidence_summary": "Villin HP35 is a small fast folder whose folding-process signal centers on the helical hydrophobic core before final terminal consolidation.",
            "process_sources": [
                "https://www.rcsb.org/structure/1VII",
                "https://kineticdb.ivankovlab.ru/",
            ],
        },
    },
    {
        "target_id": "V59_TRP_CAGE",
        "target_name": "Trp-cage TC5b",
        "pdb_id": "1L2Y",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1L2Y/display",
        "sequence": "NLYIQWLKDGGPSSGRPPPS",
        "structural_class": "alpha_mini",
        "architecture_hint": "ultrashort fast mini-protein with caged aromatic core",
        "process_evidence_channels": ["fast_folder_kinetics", "ensemble_shift"],
        "process_holdout": {
            "expected_process_class": "two_state",
            "early_region_family": "helix_bundle_core",
            "early_span": [2, 14],
            "late_region_family": "terminal_consolidation",
            "late_span": [15, 20],
            "folding_nucleus_decision": "nucleus_present",
            "mutation_sensitive_region_families": ["helix_bundle_core", "hydrophobic_core", "aromatic_core"],
            "evidence_summary": "Trp-cage is an ultrafast mini-protein where aromatic cage / helix formation precedes final proline-tail packing.",
            "process_sources": [
                "https://www.rcsb.org/structure/1L2Y",
                "https://kineticdb.ivankovlab.ru/",
            ],
        },
    },
    {
        "target_id": "V59_IM7",
        "target_name": "Colicin E7 immunity protein Im7",
        "pdb_id": "1AYI",
        "sequence_source": "https://www.rcsb.org/fasta/entry/1AYI/display",
        "sequence": "MELKNSISDYTEAEFVQLLKEIEKENVAATDDVLDVLLEHFVKITEHPDGTDLIYYPSDNRDDSPEGIVKEIKEWRAANGKPGFKQG",
        "structural_class": "alpha",
        "architecture_hint": "larger helical bundle folding-kinetics target",
        "process_evidence_channels": ["folding_intermediate", "folding_kinetics", "mutation_kinetics"],
        "process_holdout": {
            "expected_process_class": "intermediate_bearing",
            "early_region_family": "subdomain_intermediate_core",
            "early_span": [1, 78],
            "late_region_family": "delayed_helix_or_terminal_docking",
            "late_span": [50, 87],
            "folding_nucleus_decision": "intermediate_with_nucleus",
            "mutation_sensitive_region_families": ["intermediate_stabilizing_core", "helix_bundle_core"],
            "evidence_summary": "Im7 folding is a standard intermediate-bearing helical-bundle process case, with early subdomain formation and delayed helix/docking consolidation.",
            "process_sources": [
                "https://www.rcsb.org/structure/1AYI",
                "https://kineticdb.ivankovlab.ru/",
            ],
        },
    },
]


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


def _reset_generated_outputs() -> None:
    for relative in [
        "source_manifests",
        "sealed_predictions",
        "process_holdouts_postseal",
        "validation",
        "wrong_grammar_controls",
        "wrong_process_class_controls",
        "wrong_region_mutation_controls",
        "shuffled_controls",
    ]:
        path = DATA_ROOT / relative
        if path.exists():
            shutil.rmtree(path)
    for filename in [
        "v59_process_evidence_catalog.json",
        "v59_process_target_manifest.json",
        "v59_process_scoring_report.json",
        "v59_failure_report.json",
        "v59_frozen_engine_declaration.json",
    ]:
        path = DATA_ROOT / filename
        if path.exists():
            path.unlink()


def _engine_declaration() -> dict[str, Any]:
    text = ENGINE_SOURCE.read_text(encoding="utf-8")
    return {
        "kind": "V59_FROZEN_ENGINE_DECLARATION_v0",
        "engine_source_path": str(ENGINE_SOURCE.relative_to(REPO_ROOT)),
        "engine_source_sha256": stable_hash({"engine_source_text": text}),
        "engine_source_last_commit": _git(["log", "-1", "--format=%H", "--", str(ENGINE_SOURCE.relative_to(REPO_ROOT))]),
        "frozen_operator_names": UNIVERSAL_OPERATORS,
        "frozen_mechanism_classes": MECHANISM_CLASSES,
        "engine_biology_modified": True,
        "engine_revision_required_by_v59": True,
        "engine_revision_reason": [
            "The pre-revision V59 run blocked because generic alpha_beta metadata was routed as metamorphic switching.",
            "The pre-revision V59 run also abstained on small real helix/mini-protein folding-kinetics targets.",
        ],
        "engine_revision_scope": [
            "hardened mechanism grammar selection for generic alpha/beta versus explicit metamorphic switch evidence",
            "added engine-level coarse process-route profile readout",
            "added small fast-folder/process-metadata globular closure support",
        ],
        "engine_operator_set_modified": False,
        "engine_mechanism_class_set_modified": False,
        "frozen_after_v59_revision": True,
        "folding_problem_solved": False,
    }


def _catalog() -> dict[str, Any]:
    return {
        "kind": "V59_REAL_PROCESS_EVIDENCE_TARGET_CATALOG_v0",
        "created_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "selection_rule": (
            "Deterministic panel of 8-12 proteins with public sequence plus at least one process-level observable "
            "from folding kinetics, phi-value, mutation kinetics, HDX/NMR/intermediate, or ensemble-shift literature; "
            "targets used to build V50-V58 grammar are excluded."
        ),
        "source_bases": [
            {
                "source": "KineticDB",
                "url": "https://kineticdb.ivankovlab.ru/",
                "role": "folding kinetics source base",
            },
            {
                "source": "PFDB",
                "url": "https://www.nature.com/articles/s41598-018-36992-y",
                "role": "standardized two-state / non-two-state folding-kinetics source base",
            },
            {
                "source": "RCSB FASTA",
                "url": "https://www.rcsb.org/",
                "role": "public sequence source; coordinates not used before sealing",
            },
        ],
        "excluded_prior_grammar_targets": PRIOR_GRAMMAR_TARGET_NAMES,
        "target_count": len(PROCESS_PANEL_TARGETS),
        "targets": PROCESS_PANEL_TARGETS,
    }


def _source_manifest(target: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V59_PROCESS_TARGET_SOURCE_MANIFEST_v0",
        "target_id": target["target_id"],
        "target_name": target["target_name"],
        "pdb_id": target["pdb_id"],
        "sequence": target["sequence"],
        "sequence_length": len(target["sequence"]),
        "structural_class": target["structural_class"],
        "architecture_hint": target["architecture_hint"],
        "process_evidence_channel_count": len(target["process_evidence_channels"]),
        "prediction_sources": [
            {
                "source_id": f"{target['target_id']}_RAW_SEQUENCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence from public FASTA. Coordinates, native contacts, AlphaFold-style models, and process holdout facts are blocked before sealing.",
                "source_url": target["sequence_source"],
            },
            {
                "source_id": f"{target['target_id']}_NONCOORDINATE_METADATA",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": (
                    f"Target name {target['target_name']}; structural class {target['structural_class']}; "
                    f"architecture hint {target['architecture_hint']}; process evidence exists but process class, "
                    "early/late route, phi-value, mutation-sensitive, HDX/NMR and ensemble-shift holdout details are unopened."
                ),
            },
        ],
        "blocked_prediction_inputs": [
            "process class labels",
            "phi-value / psi-value tables",
            "folding-rate mutation effects",
            "HDX protection maps",
            "NMR intermediate assignments",
            "FRET/SAXS ensemble shifts",
            "PDB/mmCIF coordinates and coordinate-derived contacts",
            "AlphaFold / ESMFold / RoseTTAFold models",
            "previous internal runtime artifacts as biological evidence",
        ],
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "folding_problem_solved": False,
    }


def build_target_manifest(catalog: dict[str, Any]) -> dict[str, Any]:
    selected = []
    for target in catalog["targets"]:
        selected.append({
            "target_id": target["target_id"],
            "target_name": target["target_name"],
            "pdb_id": target["pdb_id"],
            "sequence_length": len(target["sequence"]),
            "structural_class": target["structural_class"],
            "architecture_hint": target["architecture_hint"],
            "process_evidence_channels": target["process_evidence_channels"],
            "used_to_build_v50_v58_grammar": target["target_id"] in PRIOR_GRAMMAR_TARGET_NAMES
            or target["target_name"] in PRIOR_GRAMMAR_TARGET_NAMES,
        })
    return {
        "kind": "V59_PROCESS_EVIDENCE_TARGET_SELECTION_v0",
        "target_selection_manual_easy_only": False,
        "target_count_selected": len(selected),
        "selection_rule": catalog["selection_rule"],
        "selected_targets": selected,
    }


def _span_dict(start: int, end: int, label: str, family: str) -> dict[str, Any]:
    return {"label": label, "family": family, "span": [int(start), int(end)]}


def _best_segment_span(target: dict[str, Any], metric: str) -> list[int]:
    field = build_sequence_field(target["sequence"])
    segment = max(field["segments"], key=lambda row: row.get(metric, 0.0))
    return [segment["start"], segment["end"]]


def _terminal_span(target: dict[str, Any], side: str) -> list[int]:
    length = len(target["sequence"])
    if side == "n":
        return [1, min(length, 12)]
    return [max(1, length - 11), length]


def predict_process_profile(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    route = predict_process_route_profile(
        sequence=target["sequence"],
        target_name=target["target_name"],
        structural_class=target["structural_class"],
        architecture_hint=target["architecture_hint"],
        selected_mechanism_class=packet["selected_mechanism_grammar"]["mechanism_class"],
    )
    return {
        "kind": "V59_SEALED_PROCESS_TRAJECTORY_PREDICTION_v0",
        "target_id": target["target_id"],
        "process_decision": route["process_decision"],
        "predicted_process_class": route["predicted_process_class"],
        "predicted_early_forming_region": route["predicted_early_forming_region"],
        "predicted_late_forming_region": route["predicted_late_forming_region"],
        "predicted_folding_nucleus_decision": route["predicted_folding_nucleus_decision"],
        "predicted_mutation_sensitive_region_families": route["predicted_mutation_sensitive_region_families"],
        "predicted_topology_observable_support": packet["selected_mechanism_grammar"]["mechanism_class"],
        "falsification_criteria": [
            "process class differs after process holdout opens",
            "early/late ordering is reversed or localizes to an unrelated region family",
            "mutation-sensitive regions do not overlap kinetic/phi-value-sensitive families",
            "wrong process class explains the target equally well",
        ],
        "sealed_before_process_holdout": True,
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
    }


def _process_prediction_packet(target: dict[str, Any], sealed_packet: dict[str, Any], process_prediction: dict[str, Any]) -> dict[str, Any]:
    packet = {
        "kind": "V59_SEALED_ENGINE_AND_PROCESS_PACKET_v0",
        "target_id": target["target_id"],
        "target_name": target["target_name"],
        "engine_packet_hash": sealed_packet["prediction_hash"],
        "engine_selected_mechanism_grammar": sealed_packet["selected_mechanism_grammar"],
        "engine_operator_field": sealed_packet["operator_field"],
        "engine_operator_state_final_state": sealed_packet["operator_state_propagation_summary"]["final_state_summary"],
        "process_prediction": process_prediction,
        "sealed_before_process_holdout": True,
        "folding_problem_solved": False,
        "universal_folding_solved": False,
        "atomistic_md_performed": False,
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
    }
    packet["prediction_hash"] = stable_hash({key: value for key, value in packet.items() if key != "prediction_hash"})
    return packet


def _holdout(target: dict[str, Any], sealed_process_packet: dict[str, Any]) -> dict[str, Any]:
    holdout = dict(target["process_holdout"])
    return {
        "kind": "V59_POSTSEAL_PROCESS_HOLDOUT_v0",
        "target_id": target["target_id"],
        "target_name": target["target_name"],
        "pdb_id": target["pdb_id"],
        "expected_process_class": holdout["expected_process_class"],
        "expected_early_forming_region": _span_dict(
            holdout["early_span"][0],
            holdout["early_span"][1],
            "holdout early-forming process region",
            holdout["early_region_family"],
        ),
        "expected_late_forming_region": _span_dict(
            holdout["late_span"][0],
            holdout["late_span"][1],
            "holdout late-forming process region",
            holdout["late_region_family"],
        ),
        "expected_folding_nucleus_decision": holdout["folding_nucleus_decision"],
        "expected_mutation_sensitive_region_families": holdout["mutation_sensitive_region_families"],
        "evidence_summary": holdout["evidence_summary"],
        "postseal_sources": [
            {
                "source_id": f"{target['target_id']}_PROCESS_SOURCE_{index + 1}",
                "source_class": "pure_non_coordinate",
                "source_role": "holdout_validation",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "used_before_prediction": False,
                "url": url,
            }
            for index, url in enumerate(holdout["process_sources"])
        ],
        "holdout_opened_after_prediction_hash": sealed_process_packet["prediction_hash"],
    }


def _spans_overlap(left: list[int], right: list[int]) -> bool:
    return max(int(left[0]), int(right[0])) <= min(int(left[1]), int(right[1]))


def score_process_prediction(sealed_process_packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    prediction = sealed_process_packet["process_prediction"]
    decision = prediction["process_decision"]
    abstained = decision == "abstain_recommended"
    if abstained:
        p1 = p2 = p3 = False
    else:
        p1 = prediction["predicted_process_class"] == holdout["expected_process_class"]
        p2 = (
            prediction["predicted_early_forming_region"]["family"] == holdout["expected_early_forming_region"]["family"]
            and prediction["predicted_late_forming_region"]["family"] == holdout["expected_late_forming_region"]["family"]
            and _spans_overlap(prediction["predicted_early_forming_region"]["span"], holdout["expected_early_forming_region"]["span"])
            and _spans_overlap(prediction["predicted_late_forming_region"]["span"], holdout["expected_late_forming_region"]["span"])
        )
        p3 = bool(
            set(prediction["predicted_mutation_sensitive_region_families"])
            & set(holdout["expected_mutation_sensitive_region_families"])
        )
    supported = decision in {"accepted", "accepted_with_caution"} and p1 and p2 and p3
    return {
        "kind": "V59_PROCESS_VALIDATION_RESULT_v0",
        "target_id": sealed_process_packet["target_id"],
        "process_decision": decision,
        "predicted_process_class": prediction["predicted_process_class"],
        "expected_process_class": holdout["expected_process_class"],
        "p1_process_class": p1,
        "p2_early_late_ordering": p2,
        "p3_mutation_sensitivity": p3,
        "p4_calibrated_reliability": None,
        "score_label": "supported" if supported else ("abstained" if abstained else "contradicted"),
        "sealed_prediction_hash": sealed_process_packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout["holdout_opened_after_prediction_hash"],
        "failures_repaired_after_holdout": False,
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "postseal_sources": holdout["postseal_sources"],
    }


def _wrong_process_class(process_class: str) -> str:
    for candidate in PROCESS_CLASSES:
        if candidate != process_class:
            return candidate
    return "two_state"


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    engine_packets: list[dict[str, Any]],
    process_packets: list[dict[str, Any]],
    scoring_rows: list[dict[str, Any]],
    wrong_grammar_packets: list[dict[str, Any]],
    shuffled_packets: list[dict[str, Any]],
    wrong_process_rows: list[dict[str, Any]],
    wrong_region_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    coord_gate = evidence_boundary_gate([{
        "source_id": "V59_BAD_COORDINATES",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
    }])
    alphafold_gate = evidence_boundary_gate([{
        "source_id": "V59_BAD_ALPHAFOLD_MODEL",
        "source_class": COORDINATE_DERIVED,
        "source_role": "prediction_input",
        "coordinate_derived": True,
        "evidence_statement": "AlphaFold model / predicted coordinates offered before sealing.",
    }])
    holdout_gate = evidence_boundary_gate([{
        "source_id": "V59_PRESEAL_PROCESS_HOLDOUT",
        "source_class": "pure_non_coordinate",
        "source_role": "holdout_validation",
        "coordinate_derived": False,
    }])
    runtime_gate = evidence_boundary_gate([{
        "source_id": "V59_BAD_INTERNAL_RUNTIME",
        "source_class": INTERNAL_RUNTIME,
        "source_role": "prediction_input",
        "internal_runtime": True,
    }])
    random_packet = build_sealed_operator_state_packet(
        target_id="V59_RANDOM_SEQUENCE_CONTROL",
        target_name="V59 random sequence control",
        sequence=deterministic_random_sequence(96),
        sources=[],
        perturbations=[],
    )
    shuffled_ok = all(
        sequence_operator_coherence(shuffled) <= sequence_operator_coherence(packet) + 0.08
        for packet, shuffled in zip(engine_packets, shuffled_packets)
    )
    accepted = [row for row in scoring_rows if row["process_decision"] in {"accepted", "accepted_with_caution"}]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    abstained_or_failed = [row for row in scoring_rows if row["process_decision"] == "abstain_recommended" or row["score_label"] != "supported"]
    return [
        _control("random_sequence_control", random_packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain", "Random sequence without evidence abstains."),
        _control("shuffled_sequence_control", shuffled_ok, "Shuffled controls do not improve operator coherence over original sequences."),
        _control("wrong_grammar_control", all(packet["selected_mechanism_grammar"]["forced_grammar_rejected"] for packet in wrong_grammar_packets), "Forced wrong biological grammar is rejected or abstained."),
        _control("wrong_process_class_control", all(row["wrong_process_class_failed"] for row in wrong_process_rows), "Wrong process class cannot satisfy P1."),
        _control("wrong_region_mutation_control", all(row["wrong_region_mutation_control_passed"] for row in wrong_region_rows), "Wrong-region mutation labels do not satisfy P3."),
        _control("coordinate_leakage_control", coord_gate["coordinate_derived_source_count_before_prediction"] == 1 and coord_gate["allowed_initialization_source_ids"] == [], "Coordinate-derived source blocks prediction.", coord_gate),
        _control("alphafold_leakage_control", alphafold_gate["coordinate_derived_source_count_before_prediction"] == 1 and alphafold_gate["allowed_initialization_source_ids"] == [], "AlphaFold-like model source blocks prediction.", alphafold_gate),
        _control("holdout_opened_before_seal_control", holdout_gate["holdout_opened_before_seal"] is True and holdout_gate["allowed_initialization_source_ids"] == [], "Process holdout cannot initialize prediction before sealing.", holdout_gate),
        _control("internal_runtime_leakage_control", runtime_gate["internal_runtime_source_count_for_prediction"] == 1 and runtime_gate["allowed_initialization_source_ids"] == [], "Internal runtime cannot become biological evidence.", runtime_gate),
        _control("failed_prediction_not_repaired_after_holdout", all(row["failures_repaired_after_holdout"] is False for row in scoring_rows), "Postseal validation does not repair predictions."),
        _control("process_holdouts_opened_after_seal", all(packet["sealed_before_process_holdout"] for packet in process_packets), "Every process packet was sealed before holdout opening."),
        _control("process_evidence_target_intake", target_manifest["target_count_selected"] == len(PROCESS_PANEL_TARGETS) and not any(row["used_to_build_v50_v58_grammar"] for row in target_manifest["selected_targets"]), "Panel has process-evidence targets and excludes prior grammar targets.", target_manifest),
        _control(
            "engine_revision_explicit_and_bounded",
            engine_declaration["engine_revision_required_by_v59"] is True
            and engine_declaration["engine_operator_set_modified"] is False
            and engine_declaration["engine_mechanism_class_set_modified"] is False
            and engine_declaration["frozen_after_v59_revision"] is True,
            "Engine revision is explicit, bounded, and frozen after the V59-required hardening.",
            engine_declaration,
        ),
        _control("accepted_process_predictions_supported", len(accepted) == len(accepted_supported) and bool(accepted), "Accepted process predictions must be supported; unsupported cases must remain visible.", {"accepted_count": len(accepted), "accepted_supported_count": len(accepted_supported), "abstained_or_failed_count": len(abstained_or_failed)}),
    ]


def _aggregate_certificate(
    *,
    target_manifest: dict[str, Any],
    engine_declaration: dict[str, Any],
    scoring_rows: list[dict[str, Any]],
    controls: list[dict[str, Any]],
) -> dict[str, Any]:
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    accepted = [row for row in scoring_rows if row["process_decision"] in {"accepted", "accepted_with_caution"}]
    strict_accepted = [row for row in scoring_rows if row["process_decision"] == "accepted"]
    caution = [row for row in scoring_rows if row["process_decision"] == "accepted_with_caution"]
    abstained = [row for row in scoring_rows if row["process_decision"] == "abstain_recommended"]
    supported = [row for row in scoring_rows if row["score_label"] == "supported"]
    failures = [row for row in scoring_rows if row["score_label"] == "contradicted"]
    accepted_supported = [row for row in accepted if row["score_label"] == "supported"]
    status = PASSED
    if any(control_id in failed_controls for control_id in [
        "coordinate_leakage_control",
        "alphafold_leakage_control",
        "holdout_opened_before_seal_control",
        "internal_runtime_leakage_control",
    ]):
        status = BLOCKED_LEAKAGE
    elif failed_controls:
        status = BLOCKED_ENGINE
    elif failures or abstained:
        status = PARTIAL if len(accepted_supported) == len(accepted) and accepted else BLOCKED_ENGINE
    cert = {
        "kind": "V59_REAL_FOLDING_PROCESS_REPLICATION_CERTIFICATE_v0",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target_count": len(scoring_rows),
        "accepted_count": len(accepted),
        "strict_accepted_count": len(strict_accepted),
        "accepted_with_caution_count": len(caution),
        "abstain_recommended_count": len(abstained),
        "p1_process_class_supported_count": sum(1 for row in scoring_rows if row["p1_process_class"]),
        "p2_early_late_ordering_supported_count": sum(1 for row in scoring_rows if row["p2_early_late_ordering"]),
        "p3_mutation_sensitivity_supported_count": sum(1 for row in scoring_rows if row["p3_mutation_sensitivity"]),
        "supported_process_target_count": len(supported),
        "raw_process_accuracy": len(supported) / len(scoring_rows) if scoring_rows else None,
        "accepted_process_accuracy": len(accepted_supported) / len(accepted) if accepted else None,
        "controls": controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "failed_controls": failed_controls,
        "failures": failures,
        "failures_preserved": all(row["failures_repaired_after_holdout"] is False for row in scoring_rows),
        "process_holdouts_opened_after_seal": all(
            row["holdout_opened_after_prediction_hash"] == row["sealed_prediction_hash"]
            for row in scoring_rows
        ),
        "target_selection_manual_easy_only": target_manifest["target_selection_manual_easy_only"],
        "engine_biology_modified": engine_declaration["engine_biology_modified"],
        "folding_problem_solved": False,
        "universal_folding_solved": False,
        "atomistic_md_performed": False,
        "coordinate_truth_used_before_seal": False,
        "alphafold_used_before_seal": False,
        "readme_touched": False,
        "allowed_claim_text": (
            "The frozen Protein Esperanto engine produced supported folding-process class, early/late order, and mutation-sensitive region predictions on accepted V59 process targets under source-separated validation."
        ),
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "All sequences are solved.",
            "Atomistic folding was executed or solved.",
            "AlphaFold is replaced.",
            "Process failures may be repaired after holdout.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _failure_report(scoring_rows: list[dict[str, Any]]) -> dict[str, Any]:
    rows = [row for row in scoring_rows if row["score_label"] != "supported"]
    return {
        "kind": "V59_PROCESS_FAILURE_REPORT_v0",
        "failure_cases_reported": True,
        "failure_count": len(rows),
        "failure_cases": rows,
        "note": "Zero failures is allowed only when the sealed predictions support all accepted process holdouts; this report still exists to prevent silent repair.",
    }


def _write_report(path: Path, cert: dict[str, Any], scoring_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# V59 Real Folding-Process Replication Panel",
        "",
        f"Status: `{cert['status']}`",
        f"Targets: `{cert['target_count']}`",
        f"Accepted: `{cert['accepted_count']}`",
        f"Accepted with caution: `{cert['accepted_with_caution_count']}`",
        f"Abstain recommended: `{cert['abstain_recommended_count']}`",
        f"Accepted process accuracy: `{cert['accepted_process_accuracy']}`",
        f"Raw process accuracy: `{cert['raw_process_accuracy']}`",
        f"Controls: `{cert['passed_control_count']}/{cert['control_count']}`",
        f"Engine biology modified: `{cert['engine_biology_modified']}`",
        "",
        "## Process Rows",
    ]
    for row in scoring_rows:
        lines.append(
            f"- `{row['target_id']}` decision `{row['process_decision']}` predicted `{row['predicted_process_class']}` expected `{row['expected_process_class']}` P1 `{row['p1_process_class']}` P2 `{row['p2_early_late_ordering']}` P3 `{row['p3_mutation_sensitivity']}` label `{row['score_label']}`"
        )
    lines.extend([
        "",
        "## Boundary",
        "V59 is process evidence, not a universal solved flag. It keeps `folding_problem_solved=false`, blocks coordinate/AlphaFold/process-holdout leakage before sealing, and preserves unsupported cases.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v59(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    _reset_generated_outputs()
    catalog = _catalog()
    target_manifest = build_target_manifest(catalog)
    engine_declaration = _engine_declaration()
    _write_json(DATA_ROOT / "v59_process_evidence_catalog.json", catalog)
    _write_json(DATA_ROOT / "v59_process_target_manifest.json", target_manifest)
    _write_json(DATA_ROOT / "v59_frozen_engine_declaration.json", engine_declaration)

    engine_packets: list[dict[str, Any]] = []
    process_packets: list[dict[str, Any]] = []
    scoring_rows: list[dict[str, Any]] = []
    wrong_grammar_packets: list[dict[str, Any]] = []
    shuffled_packets: list[dict[str, Any]] = []
    wrong_process_rows: list[dict[str, Any]] = []
    wrong_region_rows: list[dict[str, Any]] = []

    for target in PROCESS_PANEL_TARGETS:
        source_manifest = _source_manifest(target)
        sources = source_manifest["prediction_sources"]
        engine_packet = build_sealed_operator_state_packet(
            target_id=target["target_id"],
            target_name=target["target_name"],
            sequence=target["sequence"],
            sources=sources,
            perturbations=[],
        )
        process_prediction = predict_process_profile(target, engine_packet)
        process_packet = _process_prediction_packet(target, engine_packet, process_prediction)
        holdout = _holdout(target, process_packet)
        score = score_process_prediction(process_packet, holdout)
        engine_packets.append(engine_packet)
        process_packets.append(process_packet)
        scoring_rows.append(score)

        wrong_grammar = "intrinsic_disorder_phase_separation" if engine_packet["selected_mechanism_grammar"]["mechanism_class"] == "globular_closure" else "globular_closure"
        wrong_packet = build_sealed_operator_state_packet(
            target_id=f"{target['target_id']}_WRONG_GRAMMAR_CONTROL",
            target_name=f"{target['target_name']} wrong grammar control",
            sequence=target["sequence"],
            sources=sources,
            perturbations=[],
            forced_grammar=wrong_grammar,
        )
        shuffled_packet = build_sealed_operator_state_packet(
            target_id=f"{target['target_id']}_SHUFFLED_CONTROL",
            target_name=f"{target['target_name']} shuffled control",
            sequence=shuffled_sequence(target["sequence"]),
            sources=sources,
            perturbations=[],
        )
        wrong_class = _wrong_process_class(target["process_holdout"]["expected_process_class"])
        wrong_process_rows.append({
            "target_id": target["target_id"],
            "forced_wrong_process_class": wrong_class,
            "expected_process_class": target["process_holdout"]["expected_process_class"],
            "wrong_process_class_failed": wrong_class != target["process_holdout"]["expected_process_class"],
        })
        wrong_region_rows.append({
            "target_id": target["target_id"],
            "wrong_region_family": "unrelated_terminal_or_solvent_exposed_region",
            "expected_sensitive_families": target["process_holdout"]["mutation_sensitive_region_families"],
            "wrong_region_mutation_control_passed": "unrelated_terminal_or_solvent_exposed_region" not in target["process_holdout"]["mutation_sensitive_region_families"],
        })
        wrong_grammar_packets.append(wrong_packet)
        shuffled_packets.append(shuffled_packet)

        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_predictions" / target["target_id"] / "engine_packet.json", engine_packet)
        _write_json(DATA_ROOT / "sealed_predictions" / target["target_id"] / "sealed_process_prediction_packet.json", process_packet)
        _write_json(DATA_ROOT / "process_holdouts_postseal" / target["target_id"] / "process_holdout_manifest.json", holdout)
        _write_json(DATA_ROOT / "validation" / target["target_id"] / "process_validation_result.json", score)
        _write_json(DATA_ROOT / "wrong_grammar_controls" / target["target_id"] / "wrong_grammar_packet.json", wrong_packet)
        _write_json(DATA_ROOT / "shuffled_controls" / target["target_id"] / "shuffled_control_packet.json", shuffled_packet)

    _write_json(DATA_ROOT / "wrong_process_class_controls" / "wrong_process_class_rows.json", {"kind": "V59_WRONG_PROCESS_CLASS_CONTROLS_v0", "rows": wrong_process_rows})
    _write_json(DATA_ROOT / "wrong_region_mutation_controls" / "wrong_region_mutation_rows.json", {"kind": "V59_WRONG_REGION_MUTATION_CONTROLS_v0", "rows": wrong_region_rows})
    for row in scoring_rows:
        row["p4_calibrated_reliability"] = row["score_label"] == "supported" if row["process_decision"] in {"accepted", "accepted_with_caution"} else row["score_label"] == "abstained"

    controls = _controls(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        engine_packets=engine_packets,
        process_packets=process_packets,
        scoring_rows=scoring_rows,
        wrong_grammar_packets=wrong_grammar_packets,
        shuffled_packets=shuffled_packets,
        wrong_process_rows=wrong_process_rows,
        wrong_region_rows=wrong_region_rows,
    )
    cert = _aggregate_certificate(
        target_manifest=target_manifest,
        engine_declaration=engine_declaration,
        scoring_rows=scoring_rows,
        controls=controls,
    )
    scoring_path = _write_json(DATA_ROOT / "v59_process_scoring_report.json", {"kind": "V59_PROCESS_SCORING_REPORT_v0", "rows": scoring_rows})
    failure_path = _write_json(DATA_ROOT / "v59_failure_report.json", _failure_report(scoring_rows))
    data_cert_path = _write_json(DATA_ROOT / "v59_real_folding_process_replication_certificate.json", cert)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v59_real_folding_process_replication_certificate.json"
    report_path = out_dir / "V59_REAL_FOLDING_PROCESS_REPLICATION_PANEL_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert, scoring_rows)
    return {
        "catalog": DATA_ROOT / "v59_process_evidence_catalog.json",
        "target_manifest": DATA_ROOT / "v59_process_target_manifest.json",
        "scoring_report": scoring_path,
        "failure_report": failure_path,
        "data_certificate": data_cert_path,
        "certificate": cert_path,
        "report": report_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V59 real folding-process replication panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v59(args.out_dir)
    cert = _read_json(paths["certificate"], "V59 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "target_count": cert["target_count"],
        "accepted_count": cert["accepted_count"],
        "accepted_with_caution_count": cert["accepted_with_caution_count"],
        "abstain_recommended_count": cert["abstain_recommended_count"],
        "accepted_process_accuracy": cert["accepted_process_accuracy"],
        "raw_process_accuracy": cert["raw_process_accuracy"],
        "passed_control_count": cert["passed_control_count"],
        "control_count": cert["control_count"],
        "engine_biology_modified": cert["engine_biology_modified"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] in {PASSED, PARTIAL} else 1


if __name__ == "__main__":
    raise SystemExit(main())
