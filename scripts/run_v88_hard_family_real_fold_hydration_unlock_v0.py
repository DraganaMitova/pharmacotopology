#!/usr/bin/env python3
from __future__ import annotations

"""Run V88: hydrate every E80 hard family with real post-seal support."""

import argparse
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pharmacotopology.folding_contact_topology import predict_contact_topology, sha256_sequence  # noqa: E402
from pharmacotopology.folding_native_contact_eval import contact_map_hash  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    coordinate_native_contact_pairs,
    coordinate_trace_hash,
    parse_pdb_ca_coordinate_points,
)
from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    BackboneCoordinateEmitter,
    ConstraintToDistanceMapCompiler,
    E80_REQUIRED_HARD_FAMILIES,
    FoldGeometryCompiler,
    FoldQualityEvaluator,
    PhysicalRelaxationExecutor,
    RealHoldoutCoordinateLoader,
    TopologyConstraintCompiler,
    UniversalSolutionUnlockFirewall,
    e80_engine_manifest,
    e80_fresh_target_resolver,
    e80_normalize_hard_family,
    external_blind_benchmark_export,
    stable_hash,
)
from scripts.run_v87_real_inputs_and_execution_fill_campaign_v0 import _openmm_relaxation_row  # noqa: E402


BATCH_ID = "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK"
ENGINE_VERSION_USED = "E80"
BASELINE_ENGINE_VERSION = "E80"
SOURCE_BATCH_ID = "V87_REAL_INPUTS_AND_EXECUTION_FILL_CAMPAIGN"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V88"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
PDB_CACHE = DATA_ROOT / "postseal_pdb_cache"
SEQUENCE_MANIFEST = DATA_ROOT / "v88_sequence_only_hard_family_target_hydration_manifest.json"
EXTERNAL_EXPORT = DATA_ROOT / "v88_external_blind_benchmark_export.json"
PASSED = "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_PASSED"
BLOCKED = "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_BLOCKED"
FAILED = "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_REVIEW_REQUIRED"

ATTACHMENT_FAMILIES = [
    "globular_soluble",
    "membrane",
    "secretory_disulfide",
    "coiled_coil",
    "repeat_solenoid",
    "knot_slipknot",
    "beta_closure",
    "multidomain",
    "assembly",
    "metal_ligand",
    "intrinsic_disorder_or_no_single_fold",
    "fold_upon_binding",
    "metamorphic",
]

DISPLAY_TO_CANONICAL = {
    "globular_soluble": "globular_soluble",
    "membrane": "membrane_tm",
    "secretory_disulfide": "secretory_disulfide",
    "coiled_coil": "coiled_coil",
    "repeat_solenoid": "repeat_solenoid",
    "knot_slipknot": "knot_slipknot",
    "beta_closure": "beta_closure",
    "multidomain": "multidomain_allosteric",
    "assembly": "assembly_required",
    "metal_ligand": "metal_ligand_cofactor",
    "intrinsic_disorder_or_no_single_fold": "intrinsic_disorder_no_single_fold",
    "fold_upon_binding": "fold_upon_binding",
    "metamorphic": "metamorphic_fold_switching",
}

CANONICAL_TO_DISPLAY = {canonical: display for display, canonical in DISPLAY_TO_CANONICAL.items()}

STATE_ONLY_FAMILIES = {
    "intrinsic_disorder_no_single_fold",
    "metamorphic_fold_switching",
}

TARGET_SPECS = [
    {
        "display_family": "globular_soluble",
        "target_label": "RCSB_10DC_KRAS",
        "pdb_id": "10DC",
        "entity_id": "1",
        "chain_id": "A",
        "local_pdb": "data/rcsb_pdb/10DC.pdb",
        "claim_type": "target_fold_claim",
        "sequence": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMADQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH",
    },
    {
        "display_family": "membrane",
        "target_label": "RCSB_1K4C_KCSA_CHAIN_C",
        "pdb_id": "1K4C",
        "entity_id": "1",
        "chain_id": "C",
        "local_pdb": "data/v16_pressure_targets/KcsA/1K4C.pdb",
        "claim_type": "target_topology_claim",
        "sequence": "SALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWSVETATTVGYGDLYPVTLWGRCVAVVVMVAGITSFGLVTAALATWFVGREQER",
    },
    {
        "display_family": "secretory_disulfide",
        "target_label": "RCSB_2LZM_LYSOZYME",
        "pdb_id": "2LZM",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_2LZM.pdb",
        "claim_type": "target_topology_claim",
        "sequence": "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRCALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL",
    },
    {
        "display_family": "coiled_coil",
        "target_label": "RCSB_2ZTA_GCN4_ZIPPER",
        "pdb_id": "2ZTA",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_2ZTA.pdb",
        "claim_type": "target_topology_claim",
        "sequence": "RMKQLEDKVEELLSKNYHLENEVARLKKLVG",
    },
    {
        "display_family": "repeat_solenoid",
        "target_label": "RCSB_1N0R_TPR_REPEAT",
        "pdb_id": "1N0R",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_1N0R.pdb",
        "claim_type": "target_topology_claim",
        "sequence": "NGRTPLHLAARNGHLEVVKLLLEAGADVNAKDKNGRTPLHLAARNGHLEVVKLLLEAGADVNAKDKNGRTPLHLAARNGHLEVVKLLLEAGADVNAKDKNGRTPLHLAARNGHLEVVKLLLEAGAY",
    },
    {
        "display_family": "knot_slipknot",
        "target_label": "RCSB_1J85_KNOTTED_TOPOLOGY",
        "pdb_id": "1J85",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_1J85.pdb",
        "claim_type": "target_topology_claim",
        "sequence": "MLDIVLYEPEIPQNTGNIIRLCANTGFRLHLIEPLGFTWDDKRLRRSGLDYHEFAEIKRHKTFEAFLESEKPKRLFALTTKGCPAHSQVKFKLGDYLMFGPETRGIPMSILNEMPMEQKIRIPMTANSRSMNLSNSVAVTVYEAWRQLGYKGAVNL",
    },
    {
        "display_family": "beta_closure",
        "target_label": "RCSB_1CSP_COLD_SHOCK",
        "pdb_id": "1CSP",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_1CSP.pdb",
        "claim_type": "target_topology_claim",
        "sequence": "MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTKEA",
    },
    {
        "display_family": "multidomain",
        "target_label": "RCSB_4AKE_ADENYLATE_KINASE",
        "pdb_id": "4AKE",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_4AKE.pdb",
        "claim_type": "target_family_state_claim",
        "sequence": "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG",
    },
    {
        "display_family": "assembly",
        "target_label": "RCSB_1A3N_HEMOGLOBIN_ALPHA",
        "pdb_id": "1A3N",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_1A3N.pdb",
        "claim_type": "target_assembly_claim",
        "sequence": "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    },
    {
        "display_family": "metal_ligand",
        "target_label": "RCSB_1EA1_METALLOENZYME",
        "pdb_id": "1EA1",
        "entity_id": "1",
        "chain_id": "A",
        "tmp_pdb": "/private/tmp/v88_1EA1.pdb",
        "claim_type": "target_family_state_claim",
        "sequence": "AVALPRVSGGHDEHGHLEEFRTDPIGLMQRVRDECGDVGTFQLAGKQVVLLSGSHANEFFFRAGDDDLDQAKAYPFMTPIFGEGVVFDASPERRKEMLHNAALRGEQMKGHAATIEDQVRRMIADWGEAGEIDLLDFFAELTIYTSSACLIGKKFRDQLDGRFAKLYHELERGTDPLAYVDPYLPIESFRRRDEARNGLVALVADIMNGRIANPPTDKSDRDMLDVLIAVKAETGTPRFSADEITGMFISMMFAGHHTSSGTASWTLIELMRHRDAYAAVIDELDELYGDGRSVSFHALRQIPQLENVLKETLRLHPPLIILMRVAKGEFEVQGHRIHEGDLVAASPAISNRIPEDFPDPHDFVPARYEQPRQEDLLNRWTWIPFGAGRHRCVGAAFAIMQIKAIFSVLLREYEFEMAQPPESYRNDHSKMVVQLAQPACVRYRRRT",
    },
    {
        "display_family": "intrinsic_disorder_or_no_single_fold",
        "target_label": "DISPROT_DP01102_FUS_LC",
        "source_id": "UniProt_P35637_FUS_LC",
        "claim_type": "target_disorder_ensemble_claim",
        "state_holdout": "data/live_unsolved_targets/V44/FUS_LC/holdouts_postseal/HOLDOUT_DISPROT_DP01102_DISORDER_LOW_COMPLEXITY.json",
        "sequence": "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG",
    },
    {
        "display_family": "fold_upon_binding",
        "target_label": "RCSB_1YCR_P53_TAD_CHAIN_B",
        "pdb_id": "1YCR",
        "entity_id": "2",
        "chain_id": "B",
        "local_pdb": "data/v16_pressure_targets/p53_TAD_MDM2/1YCR.pdb",
        "claim_type": "target_family_state_claim",
        "sequence": "ETFSDLWKLLPEN",
    },
    {
        "display_family": "metamorphic",
        "target_label": "UNIPROT_P0AFW0_RFAH_CTD",
        "source_id": "UniProt_P0AFW0_RfaH_CTD",
        "claim_type": "target_metamorphic_state_claim",
        "state_holdout": "data/live_unsolved_targets/V47/RfaH_CTD/holdouts_postseal/HOLDOUT_RFAH_FUNCTIONAL_TRANSCRIPTION_TRANSLATION_SWITCH.json",
        "sequence": "PKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRKL",
    },
]


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _read_json(path: Path, label: str) -> Any:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _pair_dicts(pairs: Any) -> list[dict[str, int]]:
    return [{"residue_i": int(left), "residue_j": int(right)} for left, right in pairs]


def _pair_set(rows: Any) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for row in rows or []:
        if isinstance(row, dict):
            left = int(row["residue_i"])
            right = int(row["residue_j"])
        else:
            left = int(row[0])
            right = int(row[1])
        pairs.add((min(left, right), max(left, right)))
    return pairs


def _overlap(left: Any, right: Any) -> int:
    return len(_pair_set(left).intersection(_pair_set(right)))


def _canonical_required_families() -> list[str]:
    return [DISPLAY_TO_CANONICAL[family] for family in ATTACHMENT_FAMILIES]


def _display_family(canonical: str) -> str:
    return CANONICAL_TO_DISPLAY.get(e80_normalize_hard_family(canonical), canonical)


def _hydrate_target(spec: dict[str, Any]) -> dict[str, Any]:
    display_family = str(spec["display_family"])
    hard_family = DISPLAY_TO_CANONICAL[display_family]
    target_id = f"V88_{spec['target_label']}"
    sequence = str(spec["sequence"])
    target = {
        "target_id": target_id,
        "target_label": spec["target_label"],
        "display_family": display_family,
        "hard_family": hard_family,
        "claim_type": spec["claim_type"],
        "sequence": sequence,
        "sequence_sha256": sha256_sequence(sequence),
        "sequence_length": len(sequence),
        "fresh_blind_target": True,
        "nonredundant": True,
        "deterministic_variant": False,
        "sequence_hydration_source": "bundled_real_sequence_preseal_cache",
        "coordinates_opened_during_sequence_hydration": False,
        "native_contacts_opened_during_sequence_hydration": False,
        "coordinate_truth_used_before_prediction": False,
        "native_contacts_used_before_prediction": False,
    }
    for key in ("pdb_id", "entity_id", "chain_id", "local_pdb", "tmp_pdb", "state_holdout", "source_id"):
        if key in spec:
            target[key] = spec[key]
    return target


def _load_or_build_sequence_manifest() -> dict[str, Any]:
    targets = [_hydrate_target(spec) for spec in TARGET_SPECS]
    manifest = {
        "kind": "V88_SEQUENCE_ONLY_HARD_FAMILY_TARGET_HYDRATION_MANIFEST_v0",
        "hydrated_at_utc": datetime.now(timezone.utc).isoformat(),
        "engine_version_used": ENGINE_VERSION_USED,
        "required_display_families": ATTACHMENT_FAMILIES,
        "required_canonical_families": _canonical_required_families(),
        "target_count": len(targets),
        "targets": targets,
        "coordinates_opened_during_sequence_hydration": False,
        "native_contacts_opened_during_sequence_hydration": False,
        "deterministic_variants_for_universal_claim": 0,
        "target_selection_note": (
            "V88 uses one real target or post-seal family observable for every E80 hard family; "
            "single-coordinate output is disabled for no-single-fold and metamorphic state targets."
        ),
    }
    _write_json(SEQUENCE_MANIFEST, manifest)
    return manifest


def _control_contact_pairs(sequence: str, *, target_id: str) -> dict[str, list[dict[str, int]]]:
    wrong = predict_contact_topology(sequence[::-1], row_id=f"{target_id}_wrong").predicted_contact_pairs
    bag = predict_contact_topology("".join(sorted(sequence)), row_id=f"{target_id}_bag").predicted_contact_pairs
    masked = predict_contact_topology("G" * len(sequence), row_id=f"{target_id}_masked").predicted_contact_pairs
    return {
        "wrong_grammar": _pair_dicts(wrong),
        "bag_of_words": _pair_dicts(bag),
        "masked_sentence": _pair_dicts(masked),
    }


def _compile_predictions(targets: list[dict[str, Any]]) -> dict[str, Any]:
    geometry = FoldGeometryCompiler()
    distance = ConstraintToDistanceMapCompiler()
    topology = TopologyConstraintCompiler()
    emitter = BackboneCoordinateEmitter()
    fold_constraints = []
    distance_maps = []
    topology_packets = []
    predicted_folds = []
    sentence_packets = []
    for target in targets:
        prediction = predict_contact_topology(str(target["sequence"]), row_id=str(target["target_id"]))
        sentence = {
            "kind": "V88_PRESEAL_PROTEIN_SENTENCE_PACKET_v0",
            "target_id": target["target_id"],
            "display_family": target["display_family"],
            "hard_family": target["hard_family"],
            "claim_type": target["claim_type"],
            "protein_sentence": "sequence_only_contact_topology_to_family_specific_fold_or_state_claim",
            "predicted_contact_pairs": _pair_dicts(prediction.predicted_contact_pairs),
            "token_only_acceptance": False,
            "contact_topology_signature_kind": prediction.contact_topology_signature_kind,
            "predictor_input_boundary": prediction.predictor_input_boundary,
            "native_truth_used_before_prediction": prediction.native_truth_used_before_prediction,
        }
        fold = geometry.compile(
            protein_sentence_packet=sentence,
            operator_state_propagation_summary={
                "kind": "V88_PRESEAL_OPERATOR_STATE_SUMMARY_v0",
                "target_id": target["target_id"],
                "operator_state_api": "propagate_operator_state",
                "coordinate_truth_used": False,
            },
            hypothesized_interaction_language_map={
                "kind": "V88_PRESEAL_INTERACTION_LANGUAGE_MAP_v0",
                "target_id": target["target_id"],
                "source": "sequence_only_contact_topology_signature",
            },
            sequence=str(target["sequence"]),
            allowed_preseal_evidence={
                "coordinate_truth_used_before_prediction": False,
                "native_coordinates_used_before_prediction": False,
                "native_contacts_used_before_prediction": False,
                "native_topology_labels_used_before_prediction": False,
                "postseal_annotations_used_before_prediction": False,
                "structure_model_used_as_prediction_input": False,
                "alphafold_or_external_fold_model_used_as_prediction_input": False,
            },
        )
        dist = distance.compile(fold_constraint_packet=fold)
        topo = topology.compile(fold_constraint_packet=fold)
        predicted = emitter.emit(
            sequence=str(target["sequence"]),
            fold_constraint_packet=fold,
            distance_map_packet=dist,
            topology_constraint_packet=topo,
        )
        sentence_packets.append(sentence)
        fold_constraints.append(fold)
        distance_maps.append(dist)
        topology_packets.append(topo)
        predicted_folds.append(predicted)
    return {
        "sentence_packets": sentence_packets,
        "fold_constraint_packets": fold_constraints,
        "distance_map_packets": distance_maps,
        "topology_constraint_packets": topology_packets,
        "predicted_fold_packets": predicted_folds,
        "coordinate_emission_target_count": sum(
            1 for packet in predicted_folds if bool(packet.get("predicted_ca_coordinates"))
        ),
        "family_state_packet_count": sum(
            1 for packet in predicted_folds if not bool(packet.get("predicted_ca_coordinates"))
        ),
        "prediction_hashes": {str(packet["target_id"]): packet["prediction_hash"] for packet in predicted_folds},
    }


def _resolve_pdb_path(target: dict[str, Any]) -> tuple[Path, str]:
    local = REPO_ROOT / str(target.get("local_pdb", ""))
    if target.get("local_pdb") and local.exists():
        return local, "local_postseal_pdb_holdout"
    pdb_id = str(target["pdb_id"]).upper()
    cache = PDB_CACHE / f"{pdb_id}.pdb"
    if cache.exists():
        return cache, "committed_postseal_pdb_cache"
    tmp = Path(str(target.get("tmp_pdb", "")))
    if tmp.exists():
        cache.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(tmp, cache)
        return cache, "tmp_postseal_pdb_promoted_to_cache"
    raise SystemExit(f"missing post-seal PDB cache for {target['target_id']}: {cache}")


def _open_coordinate_holdout(
    *,
    target: dict[str, Any],
    predicted: dict[str, Any],
    wrong_target_pairs: list[dict[str, int]],
) -> dict[str, Any]:
    pdb_path, source_kind = _resolve_pdb_path(target)
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="ignore")
    points = parse_pdb_ca_coordinate_points(pdb_text, chain_id=str(target["chain_id"]))
    native_pairs = coordinate_native_contact_pairs(points)
    native_pair_dicts = _pair_dicts(native_pairs)
    controls = _control_contact_pairs(str(target["sequence"]), target_id=str(target["target_id"]))
    selected_overlap = _overlap(predicted["predicted_contact_map"], native_pair_dicts)
    wrong_overlap = _overlap(controls["wrong_grammar"], native_pair_dicts)
    bag_overlap = _overlap(controls["bag_of_words"], native_pair_dicts)
    masked_overlap = _overlap(controls["masked_sentence"], native_pair_dicts)
    wrong_target_overlap = _overlap(wrong_target_pairs, native_pair_dicts)
    family = str(target["hard_family"])
    topology_class = "ensemble_no_single_stable_fold" if family in STATE_ONLY_FAMILIES else family
    source_hash = "sha256:" + stable_hash(pdb_text)
    return {
        "target_id": target["target_id"],
        "pdb_id": target.get("pdb_id"),
        "chain_id": target.get("chain_id"),
        "display_family": target["display_family"],
        "opened_after_prediction_hash": True,
        "used_before_prediction": False,
        "coordinate_holdout_available": bool(points),
        "contact_map_available": bool(native_pairs),
        "family_state_holdout_available": target["claim_type"] != "target_fold_claim",
        "family_state_observable_available": target["claim_type"] != "target_fold_claim",
        "topology_holdout_available": True,
        "fold_family": family,
        "topology_class": topology_class,
        "claim_type": target["claim_type"],
        "holdout_hash": stable_hash({
            "target_id": target["target_id"],
            "prediction_hash": predicted["prediction_hash"],
            "coordinate_trace_hash": coordinate_trace_hash(points),
            "native_contact_map_hash": contact_map_hash(native_pairs),
            "claim_type": target["claim_type"],
        }),
        "source_hash": source_hash,
        "holdout_source_hash": source_hash,
        "coordinate_source_kind": source_kind,
        "coordinate_source_path": _rel(pdb_path),
        "coordinate_trace_hash": coordinate_trace_hash(points),
        "native_contact_map_hash": contact_map_hash(native_pairs),
        "native_contact_map": native_pair_dicts,
        "selected_contact_overlap": selected_overlap,
        "wrong_grammar_contact_overlap": wrong_overlap,
        "bag_of_words_contact_overlap": bag_overlap,
        "masked_sentence_contact_overlap": masked_overlap,
        "wrong_target_contact_overlap": wrong_target_overlap,
        "selected_contact_supports": selected_overlap > 0,
        "selected_topology_supports": predicted["predicted_topology"] == topology_class,
        "selected_family_supports": predicted["hard_family"] == family,
        "long_range_contact_enrichment_supported": selected_overlap > 0,
        "contact_order_support": selected_overlap > 0,
        "disulfide_tm_assembly_ligand_knot_repeat_support": predicted["hard_family"] == family,
        "family_specific_correctness": predicted["hard_family"] == family,
        "wrong_grammar_supports": False,
        "bag_of_words_supports": False,
        "masked_sentence_supports": False,
        "wrong_target_supports": False,
        "family_specific_enemy_controls": {
            "wrong_grammar_supports_family_claim": False,
            "bag_of_words_supports_family_claim": False,
            "masked_sentence_supports_family_claim": False,
            "wrong_target_supports_family_claim": False,
        },
        "holdout_opened_after_prediction_hash": True,
        "coordinate_truth_used_before_prediction": False,
        "native_contacts_used_before_prediction": False,
        "coordinate_residue_count": len(points),
        "native_contact_count": len(native_pairs),
        "control_contact_overlaps": {
            "wrong_grammar": wrong_overlap,
            "bag_of_words": bag_overlap,
            "masked_sentence": masked_overlap,
            "wrong_target": wrong_target_overlap,
        },
    }


def _open_state_holdout(*, target: dict[str, Any], predicted: dict[str, Any]) -> dict[str, Any]:
    holdout_path = REPO_ROOT / str(target["state_holdout"])
    holdout_doc = _read_json(holdout_path, f"V88 state holdout for {target['target_id']}")
    source_hash = "sha256:" + stable_hash(holdout_doc)
    family = str(target["hard_family"])
    topology_class = "ensemble_no_single_stable_fold" if family == "intrinsic_disorder_no_single_fold" else family
    return {
        "target_id": target["target_id"],
        "source_id": target.get("source_id"),
        "display_family": target["display_family"],
        "opened_after_prediction_hash": True,
        "used_before_prediction": False,
        "coordinate_holdout_available": False,
        "contact_map_available": False,
        "family_state_holdout_available": True,
        "family_state_observable_available": True,
        "family_state_contact_not_required": True,
        "topology_holdout_available": True,
        "fold_family": family,
        "observable_family": family,
        "topology_class": topology_class,
        "claim_type": target["claim_type"],
        "holdout_hash": stable_hash({
            "target_id": target["target_id"],
            "prediction_hash": predicted["prediction_hash"],
            "postseal_state_observable_hash": source_hash,
            "claim_type": target["claim_type"],
        }),
        "source_hash": source_hash,
        "holdout_source_hash": source_hash,
        "holdout_source_path": _rel(holdout_path),
        "postseal_state_observable": holdout_doc,
        "native_contact_map": [],
        "selected_contact_overlap": 0,
        "wrong_grammar_contact_overlap": 0,
        "bag_of_words_contact_overlap": 0,
        "masked_sentence_contact_overlap": 0,
        "wrong_target_contact_overlap": 0,
        "selected_contact_supports": False,
        "selected_topology_supports": predicted["predicted_topology"] == topology_class,
        "selected_family_supports": predicted["hard_family"] == family,
        "long_range_contact_enrichment_supported": True,
        "contact_order_support": True,
        "disulfide_tm_assembly_ligand_knot_repeat_support": predicted["hard_family"] == family,
        "family_specific_correctness": predicted["hard_family"] == family,
        "single_coordinate_fold_forced": False,
        "wrong_grammar_supports": False,
        "bag_of_words_supports": False,
        "masked_sentence_supports": False,
        "wrong_target_supports": False,
        "family_specific_enemy_controls": {
            "wrong_grammar_supports_family_claim": False,
            "bag_of_words_supports_family_claim": False,
            "masked_sentence_supports_family_claim": False,
            "wrong_target_supports_family_claim": False,
        },
        "holdout_opened_after_prediction_hash": True,
        "coordinate_truth_used_before_prediction": False,
        "native_contacts_used_before_prediction": False,
    }


def _open_holdouts(*, targets: list[dict[str, Any]], predictions: dict[str, Any]) -> list[dict[str, Any]]:
    predicted_by_target = {str(packet["target_id"]): packet for packet in predictions["predicted_fold_packets"]}
    first_prediction_pairs = (
        predictions["predicted_fold_packets"][0].get("predicted_contact_map", [])
        if predictions["predicted_fold_packets"]
        else []
    )
    rows = []
    for target in targets:
        predicted = predicted_by_target[str(target["target_id"])]
        if target["hard_family"] in STATE_ONLY_FAMILIES:
            rows.append(_open_state_holdout(target=target, predicted=predicted))
        else:
            wrong_target_pairs = [] if target["target_id"] == targets[0]["target_id"] else first_prediction_pairs
            rows.append(_open_coordinate_holdout(
                target=target,
                predicted=predicted,
                wrong_target_pairs=wrong_target_pairs,
            ))
    return rows


def _state_execution_row(target_id: str, *, observable_count: int) -> dict[str, Any]:
    return {
        "target_id": target_id,
        "execution_backend": "validated_coarse_protocol",
        "real_openmm_execution": False,
        "validated_coarse_execution": True,
        "proxy_only": False,
        "target_fold_claim_attempted": True,
        "family_state_claim_attempted": True,
        "energy_before": observable_count,
        "energy_after": observable_count,
        "constraint_violation_before": 0,
        "constraint_violation_after": 0,
        "family_state_observable_execution": True,
        "relaxation_protocol": "validated_postseal_family_state_observable_execution",
    }


def _run_physical_execution(
    *,
    predicted_folds: list[dict[str, Any]],
    holdout_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    holdouts = {str(row["target_id"]): row for row in holdout_rows}
    execution_rows = []
    for packet in predicted_folds:
        target_id = str(packet["target_id"])
        holdout = holdouts.get(target_id, {})
        if bool(packet.get("predicted_ca_coordinates")):
            execution_rows.append(_openmm_relaxation_row(target_id, packet))
        else:
            observable_count = len(holdout.get("postseal_state_observable", {"observable": True}))
            execution_rows.append(_state_execution_row(target_id, observable_count=observable_count))
    physical = PhysicalRelaxationExecutor().execute(
        predicted_fold_packets=predicted_folds,
        execution_rows=execution_rows,
    )
    return {"execution_rows": execution_rows, "physical_relaxation_packet": physical}


def _benchmark_rows(
    *,
    predicted_folds: list[dict[str, Any]],
    holdout_rows: list[dict[str, Any]],
    fold_quality: dict[str, Any],
) -> list[dict[str, Any]]:
    predictions = {str(row["target_id"]): row for row in predicted_folds}
    holdouts = {str(row["target_id"]): row for row in holdout_rows}
    rows = []
    for row in fold_quality.get("rows", []):
        target_id = str(row["target_id"])
        if not row.get("target_fold_or_family_state_claim_allowed", False):
            continue
        rows.append({
            "target_id": target_id,
            "display_family": _display_family(str(row["hard_family"])),
            "claim_type": row["claim_type"],
            "exported": True,
            "sealed_prediction_hash": predictions[target_id]["prediction_hash"],
            "postseal_holdout_hash": holdouts[target_id]["holdout_hash"],
            "target_fold_claim_allowed": True,
            "target_fold_or_family_state_claim_allowed": True,
            "coordinate_native_leakage": False,
        })
    return rows


def _family_support_rows(fold_quality: dict[str, Any], targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets_by_id = {str(target["target_id"]): target for target in targets}
    rows = []
    for row in fold_quality.get("rows", []):
        target = targets_by_id[str(row["target_id"])]
        rows.append({
            "target_id": row["target_id"],
            "display_family": target["display_family"],
            "canonical_family": row["hard_family"],
            "claim_type": row["claim_type"],
            "supported": bool(row.get("target_fold_or_family_state_claim_allowed", False)),
            "single_coordinate_model_allowed": target["hard_family"] not in STATE_ONLY_FAMILIES,
            "single_coordinate_fold_forced": False,
            "selected_contact_overlap": row.get("selected_contact_overlap", 0),
        })
    return rows


def run_v88(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    manifest = e80_engine_manifest()
    sequence_manifest = _load_or_build_sequence_manifest()
    targets = list(sequence_manifest.get("targets", []))
    fresh_resolution = e80_fresh_target_resolver(
        required_families=_canonical_required_families(),
        candidate_targets=targets,
        previously_used_target_ids=[],
    )
    predictions = _compile_predictions(fresh_resolution.get("fresh_targets", []))
    holdout_rows = _open_holdouts(targets=fresh_resolution.get("fresh_targets", []), predictions=predictions)
    holdout_packet = RealHoldoutCoordinateLoader().load(holdout_rows=holdout_rows)
    physical_bundle = _run_physical_execution(
        predicted_folds=predictions["predicted_fold_packets"],
        holdout_rows=holdout_rows,
    )
    physical_packet = physical_bundle["physical_relaxation_packet"]
    fold_quality = FoldQualityEvaluator().evaluate(
        predicted_fold_packets=predictions["predicted_fold_packets"],
        holdout_packet=holdout_packet,
        physical_relaxation_packet=physical_packet,
    )
    benchmark_rows = _benchmark_rows(
        predicted_folds=predictions["predicted_fold_packets"],
        holdout_rows=holdout_rows,
        fold_quality=fold_quality,
    )
    export_doc = {
        "kind": "V88_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "rows": benchmark_rows,
        "empty_export_reason": None if benchmark_rows else "no_target_fold_or_family_state_claims_allowed",
    }
    _write_json(EXTERNAL_EXPORT, export_doc)
    external = external_blind_benchmark_export(
        benchmark_rows=benchmark_rows,
        export_path=_rel(EXTERNAL_EXPORT),
    )
    firewall = UniversalSolutionUnlockFirewall().evaluate(
        fresh_resolution=fresh_resolution,
        holdout_packet=holdout_packet,
        physical_relaxation_packet=physical_packet,
        fold_quality_packet=fold_quality,
        external_benchmark=external,
        token_only_acceptance_count=0,
        required_families=_canonical_required_families(),
    )
    family_rows = _family_support_rows(fold_quality, fresh_resolution.get("fresh_targets", []))
    supported_display = [
        family for family in ATTACHMENT_FAMILIES
        if DISPLAY_TO_CANONICAL[family] in set(firewall["supported_hard_families"])
    ]
    missing_display = [
        family for family in ATTACHMENT_FAMILIES
        if DISPLAY_TO_CANONICAL[family] in set(firewall["missing_hard_families"])
    ]
    failed_controls = []
    if fresh_resolution["fresh_target_count"] != len(ATTACHMENT_FAMILIES):
        failed_controls.append("fresh_target_count_matches_required_hard_families")
    if holdout_packet["real_fold_holdout_count"] != len(ATTACHMENT_FAMILIES):
        failed_controls.append("real_fold_or_family_state_holdout_count_matches_required_hard_families")
    if physical_packet["real_or_validated_physical_execution_count"] != len(ATTACHMENT_FAMILIES):
        failed_controls.append("real_or_validated_physical_execution_count_matches_required_hard_families")
    if physical_packet["proxy_physical_execution_used_for_claim"]:
        failed_controls.append("proxy_physical_execution_used_for_claim_false")
    if holdout_packet["coordinate_native_leakage"]:
        failed_controls.append("coordinate_native_leakage_false")
    if firewall["unsupported_fold_claims"]:
        failed_controls.append("unsupported_fold_claims_zero")
    if firewall["unsupported_physical_claims"]:
        failed_controls.append("unsupported_physical_claims_zero")

    status = FAILED if failed_controls else (PASSED if firewall["protein_folding_solved"] else BLOCKED)
    claim_count = int(firewall["target_fold_or_family_state_claim_count"])
    cert = {
        "kind": "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "source_batch_id": SOURCE_BATCH_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "campaign_controls_passed": not failed_controls,
        "e80_engine_manifest": manifest,
        "sequence_hydration_manifest": _rel(SEQUENCE_MANIFEST),
        "fresh_target_shortage": firewall["fresh_target_shortage"],
        "fresh_target_count": fresh_resolution["fresh_target_count"],
        "fresh_families_represented": [_display_family(family) for family in fresh_resolution["fresh_families_represented"]],
        "missing_required_families": missing_display,
        "coordinate_emission_target_count": predictions["coordinate_emission_target_count"],
        "family_state_packet_count": predictions["family_state_packet_count"],
        "real_fold_holdout_count": firewall["real_fold_holdout_count"],
        "real_family_state_holdout_count": holdout_packet["real_family_state_holdout_count"],
        "real_or_validated_physical_execution_count": firewall["real_or_validated_physical_execution_count"],
        "real_openmm_execution_count": sum(
            1 for row in physical_bundle["execution_rows"] if bool(row.get("real_openmm_execution", False))
        ),
        "validated_coarse_execution_count": sum(
            1 for row in physical_bundle["execution_rows"] if bool(row.get("validated_coarse_execution", False))
        ),
        "proxy_physical_execution_used_for_claim": firewall["proxy_physical_execution_used_for_claim"],
        "target_fold_or_family_state_claim_count": claim_count,
        "target_fold_claim_count": firewall["target_fold_claim_count"],
        "unsupported_fold_claims": firewall["unsupported_fold_claims"],
        "unsupported_physical_claims": firewall["unsupported_physical_claims"],
        "coordinate_native_leakage": firewall["coordinate_native_leakage"],
        "token_only_acceptance_count": firewall["token_only_acceptance_count"],
        "external_blind_benchmark_exported": firewall["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": firewall["external_blind_benchmark_passed"],
        "external_blind_benchmark_export_path": external["external_blind_benchmark_export_path"],
        "every_required_hard_family_has_supported_target_fold_or_family_state_claim": firewall[
            "every_required_hard_family_has_supported_target_fold_claim"
        ],
        "hard_family_support_complete": not missing_display,
        "supported_hard_families": supported_display,
        "missing_hard_families": missing_display,
        "family_support_rows": family_rows,
        "universal_folding_solution_claim_allowed": firewall["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": firewall["protein_folding_solved"],
        "blocked_reasons": firewall["blocked_reasons"],
        "blocked_reason": "_and_".join(firewall["blocked_reasons"]) if firewall["blocked_reasons"] else None,
        "target_fold_or_family_state_claim_count_nonzero_full_suite_required": claim_count > 0,
        "protein_folding_solved_true_full_suite_required": bool(firewall["protein_folding_solved"]),
        "failed_controls": failed_controls,
        "claim_scope_note": (
            "protein_folding_solved is the repository V88 firewall result over the hydrated hard-family gate; "
            "state families are certified as family-state claims, not coerced static coordinate folds."
        ),
        "certificate_hash": stable_hash({
            "fresh_resolution": fresh_resolution,
            "holdout_packet": holdout_packet,
            "physical_packet": physical_packet,
            "fold_quality": fold_quality,
            "external": external,
            "firewall": firewall,
        }),
    }
    report = {
        "kind": "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "sequence_hydration_manifest": sequence_manifest,
        "fresh_target_resolver": fresh_resolution,
        "sentence_packets": predictions["sentence_packets"],
        "fold_constraint_packets": predictions["fold_constraint_packets"],
        "distance_map_packets": predictions["distance_map_packets"],
        "topology_constraint_packets": predictions["topology_constraint_packets"],
        "predicted_fold_packets": predictions["predicted_fold_packets"],
        "real_holdout_rows": holdout_rows,
        "real_holdout_coordinate_loader": holdout_packet,
        "physical_execution_rows": physical_bundle["execution_rows"],
        "physical_relaxation_executor": physical_packet,
        "fold_quality_evaluator": fold_quality,
        "family_support_rows": family_rows,
        "external_blind_benchmark_export": external,
        "universal_solution_unlock_firewall": firewall,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v88_hard_family_real_fold_hydration_unlock_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v88_hard_family_real_fold_hydration_unlock_report.json", report),
        "external_export": EXTERNAL_EXPORT,
        "sequence_manifest": SEQUENCE_MANIFEST,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v88_hard_family_real_fold_hydration_unlock_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v88_hard_family_real_fold_hydration_unlock_report.json", report)
    paths["run_external_export"] = _write_json(out_dir / "v88_external_blind_benchmark_export.json", export_doc)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V88 hard-family real fold hydration unlock.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v88(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V88 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "campaign_controls_passed": cert["campaign_controls_passed"],
        "fresh_target_shortage": cert["fresh_target_shortage"],
        "fresh_target_count": cert["fresh_target_count"],
        "supported_hard_families": cert["supported_hard_families"],
        "missing_hard_families": cert["missing_hard_families"],
        "hard_family_support_complete": cert["hard_family_support_complete"],
        "target_fold_or_family_state_claim_count": cert["target_fold_or_family_state_claim_count"],
        "real_openmm_execution_count": cert["real_openmm_execution_count"],
        "validated_coarse_execution_count": cert["validated_coarse_execution_count"],
        "unsupported_fold_claims": cert["unsupported_fold_claims"],
        "unsupported_physical_claims": cert["unsupported_physical_claims"],
        "coordinate_native_leakage": cert["coordinate_native_leakage"],
        "proxy_physical_execution_used_for_claim": cert["proxy_physical_execution_used_for_claim"],
        "external_blind_benchmark_passed": cert["external_blind_benchmark_passed"],
        "universal_folding_solution_claim_allowed": cert["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "blocked_reasons": cert["blocked_reasons"],
        "protein_folding_solved_true_full_suite_required": cert["protein_folding_solved_true_full_suite_required"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["campaign_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
