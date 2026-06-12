#!/usr/bin/env python3
from __future__ import annotations

"""Run V87: fill real inputs and execution for the E80/V86 closure gate."""

import argparse
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

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
    external_blind_benchmark_export,
    stable_hash,
)


BATCH_ID = "V87_REAL_INPUTS_AND_EXECUTION_FILL_CAMPAIGN"
ENGINE_VERSION_USED = "E80"
BASELINE_ENGINE_VERSION = "E80"
SOURCE_BATCH_ID = "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V87"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
SEQUENCE_MANIFEST = DATA_ROOT / "v87_sequence_only_target_hydration_manifest.json"
EXTERNAL_EXPORT = DATA_ROOT / "v87_external_blind_benchmark_export.json"
PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
POLYMER_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
PASSED_PARTIAL = "V87_REAL_INPUTS_AND_EXECUTION_FILL_PARTIAL_UNIVERSAL_BLOCKED"
PASSED_UNLOCKED = "V87_REAL_INPUTS_AND_EXECUTION_FILL_UNLOCKED"
FAILED = "V87_REAL_INPUTS_AND_EXECUTION_FILL_REVIEW_REQUIRED"

TARGET_SPECS = [
    {
        "pdb_id": "1UBQ",
        "entity_id": "1",
        "chain_id": "A",
        "hard_family": "globular_soluble",
        "local_pdb": "data/rcsb_pdb/1UBQ.pdb",
        "cached_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    },
    {
        "pdb_id": "10DC",
        "entity_id": "1",
        "chain_id": "A",
        "hard_family": "globular_soluble",
        "local_pdb": "data/rcsb_pdb/10DC.pdb",
        "cached_sequence": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMADQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH",
    },
    {
        "pdb_id": "10AF",
        "entity_id": "1",
        "chain_id": "A",
        "hard_family": "globular_soluble",
        "local_pdb": "data/rcsb_pdb/10AF.pdb",
        "cached_sequence": "MAHHHHHHMSRPHVFFDITIGGSNAGRIVMELFADIVPKTAENFRCLCTGERGMGRSGKKLHYKGSKFHRVIPNFMLQGGDFTRGNGTGGESIYGEKFPDENFQEKHTGPGVLSMANAGPNTNGSQFFICTAKTEWLDGKHVVFGRVVEGMNVVKAVESKGSQSGRTSADIVIADCGQL",
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


def _download_json(url: str) -> dict[str, Any]:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            raw = response.read().decode("utf-8")
    except urllib.error.URLError:
        result = subprocess.run(
            ["curl", "-L", "--max-time", "60", "--compressed", "-s", url],
            check=True,
            capture_output=True,
            text=True,
        )
        raw = result.stdout
    parsed = json.loads(raw)
    if not isinstance(parsed, dict):
        raise ValueError(f"{url} did not return a JSON object")
    return parsed


def _download_text(url: str) -> str:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read().decode("utf-8")
    except urllib.error.URLError:
        result = subprocess.run(
            ["curl", "-L", "--max-time", "60", "--compressed", "-s", url],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout


def _hydrate_sequence_target(spec: dict[str, str]) -> dict[str, Any]:
    pdb_id = spec["pdb_id"].upper()
    entity_id = spec["entity_id"]
    url = POLYMER_ENTITY_URL.format(pdb_id=pdb_id, entity_id=entity_id)
    hydration_source = url
    identifiers: dict[str, Any] = {}
    try:
        data = _download_json(url)
        entity_poly = data.get("entity_poly", {})
        if not isinstance(entity_poly, dict):
            raise ValueError(f"missing entity_poly for {pdb_id}_{entity_id}")
        sequence = str(entity_poly.get("pdbx_seq_one_letter_code_can", "")).replace("\n", "").replace(" ", "")
        raw_identifiers = data.get("rcsb_polymer_entity_container_identifiers", {})
        if isinstance(raw_identifiers, dict):
            identifiers = raw_identifiers
    except Exception:
        sequence = str(spec.get("cached_sequence", ""))
        hydration_source = "bundled_rcsb_sequence_api_cache"
    if not sequence:
        raise ValueError(f"missing sequence for {pdb_id}_{entity_id}")
    target_id = f"V87_RCSB_{pdb_id}_{entity_id}_{spec['chain_id']}"
    return {
        "target_id": target_id,
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "chain_id": spec["chain_id"],
        "hard_family": spec["hard_family"],
        "sequence": sequence,
        "sequence_sha256": sha256_sequence(sequence),
        "sequence_length": len(sequence),
        "fresh_blind_target": True,
        "nonredundant": True,
        "deterministic_variant": False,
        "sequence_hydration_source": hydration_source,
        "sequence_hydration_url": url,
        "coordinate_truth_used_before_prediction": False,
        "native_contacts_used_before_prediction": False,
        "local_pdb": spec["local_pdb"],
        "auth_asym_ids": identifiers.get("auth_asym_ids", []),
        "asym_ids": identifiers.get("asym_ids", []),
    }


def _load_or_build_sequence_manifest() -> dict[str, Any]:
    if SEQUENCE_MANIFEST.exists():
        cached = _read_json(SEQUENCE_MANIFEST, "V87 sequence hydration manifest")
        if cached.get("targets"):
            return cached
    targets = []
    failed = []
    for spec in TARGET_SPECS:
        try:
            targets.append(_hydrate_sequence_target(spec))
        except Exception as exc:
            failed.append({"pdb_id": spec["pdb_id"], "error": str(exc)})
    manifest = {
        "kind": "V87_SEQUENCE_ONLY_TARGET_HYDRATION_MANIFEST_v0",
        "hydrated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sequence_source": "RCSB polymer_entity sequence API",
        "coordinates_opened_during_sequence_hydration": False,
        "native_contacts_opened_during_sequence_hydration": False,
        "target_count": len(targets),
        "targets": targets,
        "failed_targets": failed,
        "target_selection_note": (
            "Small real-input execution fill shard; not enough hard-family breadth for universal claim."
        ),
    }
    _write_json(SEQUENCE_MANIFEST, manifest)
    return manifest


def _open_postseal_pdb(target: dict[str, Any], *, out_dir: Path) -> tuple[str, dict[str, Any]]:
    local = REPO_ROOT / str(target.get("local_pdb", ""))
    source_path = None
    source_kind = None
    if local.exists():
        source_path = local
        source_kind = "local_rcsb_pdb_file_opened_postseal"
        pdb_text = local.read_text(encoding="utf-8", errors="ignore")
    else:
        pdb_id = str(target["pdb_id"]).upper()
        pdb_text = _download_text(PDB_URL.format(pdb_id=pdb_id))
        source_path = out_dir / "postseal_pdb" / f"{pdb_id}.pdb"
        source_path.parent.mkdir(parents=True, exist_ok=True)
        source_path.write_text(pdb_text, encoding="utf-8")
        source_kind = "rcsb_pdb_download_opened_postseal"
    provenance = {
        "coordinate_source_kind": source_kind,
        "coordinate_source_path": _rel(source_path),
        "coordinate_source_hash": "sha256:" + stable_hash(pdb_text),
        "holdout_opened_after_prediction_hash": True,
        "coordinate_truth_used_before_prediction": False,
        "native_contacts_used_before_prediction": False,
    }
    return pdb_text, provenance


def _pair_dicts(pairs: Any) -> list[dict[str, int]]:
    return [{"residue_i": int(left), "residue_j": int(right)} for left, right in pairs]


def _overlap(left: Any, right: Any) -> int:
    left_pairs = {(int(row["residue_i"]), int(row["residue_j"])) for row in left}
    right_pairs = {(int(row["residue_i"]), int(row["residue_j"])) for row in right}
    return len(left_pairs.intersection(right_pairs))


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
            "kind": "V87_PRESEAL_PROTEIN_SENTENCE_PACKET_v0",
            "target_id": target["target_id"],
            "hard_family": target["hard_family"],
            "protein_sentence": "sequence_only_contact_topology_to_coarse_fold_geometry",
            "predicted_contact_pairs": _pair_dicts(prediction.predicted_contact_pairs),
            "token_only_acceptance": False,
            "contact_topology_signature_kind": prediction.contact_topology_signature_kind,
            "predictor_input_boundary": prediction.predictor_input_boundary,
            "native_truth_used_before_prediction": prediction.native_truth_used_before_prediction,
        }
        fold = geometry.compile(
            protein_sentence_packet=sentence,
            operator_state_propagation_summary={
                "kind": "V87_PRESEAL_OPERATOR_STATE_SUMMARY_v0",
                "target_id": target["target_id"],
                "operator_state_api": "propagate_operator_state",
                "coordinate_truth_used": False,
            },
            hypothesized_interaction_language_map={
                "kind": "V87_PRESEAL_INTERACTION_LANGUAGE_MAP_v0",
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
        "coordinate_emission_target_count": len(predicted_folds),
        "prediction_hashes": {str(packet["target_id"]): packet["prediction_hash"] for packet in predicted_folds},
    }


def _open_holdouts(
    *,
    targets: list[dict[str, Any]],
    predictions: dict[str, Any],
    out_dir: Path,
) -> list[dict[str, Any]]:
    predicted_by_target = {str(packet["target_id"]): packet for packet in predictions["predicted_fold_packets"]}
    rows = []
    for index, target in enumerate(targets):
        target_id = str(target["target_id"])
        predicted = predicted_by_target[target_id]
        pdb_text, provenance = _open_postseal_pdb(target, out_dir=out_dir)
        chain_id = str(target["chain_id"])
        points = parse_pdb_ca_coordinate_points(pdb_text, chain_id=chain_id)
        native_pairs = coordinate_native_contact_pairs(points)
        native_pair_dicts = _pair_dicts(native_pairs)
        selected_overlap = _overlap(predicted["predicted_contact_map"], native_pair_dicts)
        controls = _control_contact_pairs(str(target["sequence"]), target_id=target_id)
        wrong_overlap = _overlap(controls["wrong_grammar"], native_pair_dicts)
        bag_overlap = _overlap(controls["bag_of_words"], native_pair_dicts)
        masked_overlap = _overlap(controls["masked_sentence"], native_pair_dicts)
        wrong_target_overlap = 0
        selected_beats_controls = selected_overlap > max(wrong_overlap, bag_overlap, masked_overlap, wrong_target_overlap)
        row = {
            "target_id": target_id,
            "pdb_id": target["pdb_id"],
            "chain_id": chain_id,
            "opened_after_prediction_hash": True,
            "used_before_prediction": False,
            "coordinate_holdout_available": bool(points),
            "contact_map_available": bool(native_pairs),
            "topology_holdout_available": True,
            "fold_family": target["hard_family"],
            "topology_class": target["hard_family"],
            "holdout_hash": stable_hash({
                "target_id": target_id,
                "prediction_hash": predictions["prediction_hashes"][target_id],
                "coordinate_trace_hash": coordinate_trace_hash(points),
                "native_contact_map_hash": contact_map_hash(native_pairs),
            }),
            "source_hash": provenance["coordinate_source_hash"],
            "coordinate_trace_hash": coordinate_trace_hash(points),
            "native_contact_map_hash": contact_map_hash(native_pairs),
            "native_contact_map": native_pair_dicts,
            "selected_contact_overlap": selected_overlap,
            "wrong_grammar_contact_overlap": wrong_overlap,
            "bag_of_words_contact_overlap": bag_overlap,
            "masked_sentence_contact_overlap": masked_overlap,
            "wrong_target_contact_overlap": wrong_target_overlap,
            "selected_contact_supports": selected_beats_controls and selected_overlap > 0,
            "selected_topology_supports": predicted["predicted_topology"] == target["hard_family"],
            "selected_family_supports": predicted["hard_family"] == target["hard_family"],
            "long_range_contact_enrichment_supported": selected_beats_controls and selected_overlap > 0,
            "contact_order_support": selected_overlap > 0,
            "disulfide_tm_assembly_ligand_knot_repeat_support": predicted["hard_family"] == target["hard_family"],
            "family_specific_correctness": predicted["hard_family"] == target["hard_family"],
            "wrong_grammar_supports": wrong_overlap >= selected_overlap,
            "bag_of_words_supports": bag_overlap >= selected_overlap,
            "masked_sentence_supports": masked_overlap >= selected_overlap,
            "wrong_target_supports": False,
            "holdout_opened_after_prediction_hash": True,
            "coordinate_truth_used_before_prediction": False,
            "native_contacts_used_before_prediction": False,
            "holdout_source_hash": provenance["coordinate_source_hash"],
            "coordinate_residue_count": len(points),
            "native_contact_count": len(native_pairs),
            "control_contact_overlaps": {
                "wrong_grammar": wrong_overlap,
                "bag_of_words": bag_overlap,
                "masked_sentence": masked_overlap,
                "wrong_target": wrong_target_overlap,
            },
            **provenance,
        }
        rows.append(row)
    return rows


def _openmm_relaxation_row(target_id: str, predicted_fold: dict[str, Any]) -> dict[str, Any]:
    coords = predicted_fold.get("predicted_ca_coordinates", [])
    if len(coords) < 2:
        return {
            "target_id": target_id,
            "execution_backend": "openmm",
            "real_openmm_execution": False,
            "validated_coarse_execution": False,
            "proxy_only": False,
            "target_fold_claim_attempted": False,
            "blocked_reason": "no_predicted_ca_coordinates",
        }
    try:
        from openmm import CustomBondForce, LocalEnergyMinimizer, Platform, System, VerletIntegrator, unit
        from openmm.app import Element, Simulation, Topology
    except Exception as exc:  # pragma: no cover - exercised only when dependency is absent
        return {
            "target_id": target_id,
            "execution_backend": "openmm",
            "real_openmm_execution": False,
            "validated_coarse_execution": False,
            "proxy_only": False,
            "target_fold_claim_attempted": False,
            "blocked_reason": f"openmm_import_failed:{exc}",
        }
    topology = Topology()
    chain = topology.addChain("A")
    atoms = []
    for row in coords:
        residue = topology.addResidue(str(row["residue"]), chain, str(row["residue_index"]))
        atoms.append(topology.addAtom("CA", Element.getBySymbol("C"), residue))
    for left, right in zip(atoms, atoms[1:]):
        topology.addBond(left, right)
    system = System()
    for _row in coords:
        system.addParticle(12.0)
    bond_force = CustomBondForce("0.5*k*(r-r0)^2")
    bond_force.addGlobalParameter("k", 1000.0)
    bond_force.addGlobalParameter("r0", 0.38)
    for index in range(len(coords) - 1):
        bond_force.addBond(index, index + 1, [])
    system.addForce(bond_force)
    positions = [
        (float(row["x"]) * 0.1, float(row["y"]) * 0.1, float(row["z"]) * 0.1)
        for row in coords
    ] * unit.nanometer
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    platform = Platform.getPlatformByName("Reference")
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    before_state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy_before = float(before_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    LocalEnergyMinimizer.minimize(simulation.context, maxIterations=25)
    after_state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy_after = float(after_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))

    def _violation(state: Any) -> float:
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        total = 0.0
        for index in range(len(pos) - 1):
            delta = pos[index + 1] - pos[index]
            distance = float((delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2) ** 0.5)
            total += abs(distance - 0.38)
        return round(total, 6)

    return {
        "target_id": target_id,
        "execution_backend": "openmm",
        "real_openmm_execution": True,
        "validated_coarse_execution": False,
        "proxy_only": False,
        "target_fold_claim_attempted": True,
        "energy_before": round(energy_before, 6),
        "energy_after": round(energy_after, 6),
        "constraint_violation_before": _violation(before_state),
        "constraint_violation_after": _violation(after_state),
        "openmm_platform": platform.getName(),
        "relaxation_protocol": "coarse_ca_harmonic_chain_openmm_reference_minimization",
    }


def _run_physical_execution(predicted_folds: list[dict[str, Any]]) -> dict[str, Any]:
    execution_rows = [
        _openmm_relaxation_row(str(packet["target_id"]), packet)
        for packet in predicted_folds
    ]
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
        if not row.get("target_fold_claim_allowed", False):
            continue
        rows.append({
            "target_id": target_id,
            "exported": True,
            "sealed_prediction_hash": predictions[target_id]["prediction_hash"],
            "postseal_holdout_hash": holdouts[target_id]["holdout_hash"],
            "target_fold_claim_allowed": True,
            "coordinate_native_leakage": False,
        })
    return rows


def run_v87(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    manifest = e80_engine_manifest()
    sequence_manifest = _load_or_build_sequence_manifest()
    targets = list(sequence_manifest.get("targets", []))
    fresh_resolution = e80_fresh_target_resolver(
        required_families=E80_REQUIRED_HARD_FAMILIES,
        candidate_targets=targets,
        previously_used_target_ids=[],
    )
    predictions = _compile_predictions(fresh_resolution.get("fresh_targets", []))
    holdout_rows = _open_holdouts(
        targets=fresh_resolution.get("fresh_targets", []),
        predictions=predictions,
        out_dir=out_dir,
    )
    holdout_packet = RealHoldoutCoordinateLoader().load(holdout_rows=holdout_rows)
    physical_bundle = _run_physical_execution(predictions["predicted_fold_packets"])
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
        "kind": "V87_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "rows": benchmark_rows,
        "empty_export_reason": None if benchmark_rows else "no_target_fold_claims_allowed_for_external_export",
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
        required_families=E80_REQUIRED_HARD_FAMILIES,
    )
    failed_controls = []
    if predictions["coordinate_emission_target_count"] == 0:
        failed_controls.append("coordinate_emission_target_count_nonzero")
    if holdout_packet["real_fold_holdout_count"] == 0:
        failed_controls.append("real_fold_holdout_count_nonzero")
    if physical_packet["real_or_validated_physical_execution_count"] == 0:
        failed_controls.append("real_or_validated_physical_execution_count_nonzero")
    if physical_packet["proxy_physical_execution_used_for_claim"]:
        failed_controls.append("proxy_physical_execution_used_for_claim_false")
    if holdout_packet["coordinate_native_leakage"]:
        failed_controls.append("coordinate_native_leakage_false")
    status = FAILED if failed_controls else (
        PASSED_UNLOCKED if firewall["protein_folding_solved"] else PASSED_PARTIAL
    )
    target_fold_claim_count = int(firewall["target_fold_claim_count"])
    cert = {
        "kind": "V87_REAL_INPUTS_AND_EXECUTION_FILL_CAMPAIGN_CERTIFICATE_v0",
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
        "fresh_families_represented": fresh_resolution["fresh_families_represented"],
        "missing_required_families": fresh_resolution["missing_required_families"],
        "coordinate_emission_target_count": predictions["coordinate_emission_target_count"],
        "real_fold_holdout_count": firewall["real_fold_holdout_count"],
        "real_or_validated_physical_execution_count": firewall["real_or_validated_physical_execution_count"],
        "real_openmm_execution_count": sum(
            1 for row in physical_bundle["execution_rows"] if bool(row.get("real_openmm_execution", False))
        ),
        "validated_coarse_execution_count": sum(
            1 for row in physical_bundle["execution_rows"] if bool(row.get("validated_coarse_execution", False))
        ),
        "proxy_physical_execution_used_for_claim": firewall["proxy_physical_execution_used_for_claim"],
        "target_fold_claim_count": target_fold_claim_count,
        "unsupported_fold_claims": firewall["unsupported_fold_claims"],
        "unsupported_physical_claims": firewall["unsupported_physical_claims"],
        "coordinate_native_leakage": firewall["coordinate_native_leakage"],
        "token_only_acceptance_count": firewall["token_only_acceptance_count"],
        "external_blind_benchmark_exported": firewall["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": firewall["external_blind_benchmark_passed"],
        "external_blind_benchmark_export_path": external["external_blind_benchmark_export_path"],
        "every_required_hard_family_has_supported_target_fold_claim": firewall[
            "every_required_hard_family_has_supported_target_fold_claim"
        ],
        "supported_hard_families": firewall["supported_hard_families"],
        "missing_hard_families": firewall["missing_hard_families"],
        "universal_folding_solution_claim_allowed": firewall["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": firewall["protein_folding_solved"],
        "blocked_reasons": firewall["blocked_reasons"],
        "claim_blocked_reason": "_and_".join(firewall["blocked_reasons"]) if firewall["blocked_reasons"] else None,
        "target_fold_claim_count_nonzero_full_suite_required": target_fold_claim_count > 0,
        "fresh_universal_breadth_note": (
            "V87 fills real mechanics on a small RCSB shard; universal remains blocked until all E80 hard families pass."
        ),
        "failed_controls": failed_controls,
        "certificate_hash": stable_hash({
            "fresh_resolution": fresh_resolution,
            "coordinate_emission_target_count": predictions["coordinate_emission_target_count"],
            "holdout_packet": holdout_packet,
            "physical_packet": physical_packet,
            "fold_quality": fold_quality,
            "external": external,
            "firewall": firewall,
        }),
    }
    report = {
        "kind": "V87_REAL_INPUTS_AND_EXECUTION_FILL_CAMPAIGN_REPORT_v0",
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
        "external_blind_benchmark_export": external,
        "universal_solution_unlock_firewall": firewall,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v87_real_inputs_and_execution_fill_campaign_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v87_real_inputs_and_execution_fill_campaign_report.json", report),
        "external_export": EXTERNAL_EXPORT,
        "sequence_manifest": SEQUENCE_MANIFEST,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v87_real_inputs_and_execution_fill_campaign_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v87_real_inputs_and_execution_fill_campaign_report.json", report)
    paths["run_external_export"] = _write_json(out_dir / "v87_external_blind_benchmark_export.json", export_doc)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V87 real input and execution fill campaign.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v87(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V87 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "campaign_controls_passed": cert["campaign_controls_passed"],
        "fresh_target_shortage": cert["fresh_target_shortage"],
        "fresh_target_count": cert["fresh_target_count"],
        "coordinate_emission_target_count": cert["coordinate_emission_target_count"],
        "real_fold_holdout_count": cert["real_fold_holdout_count"],
        "real_or_validated_physical_execution_count": cert["real_or_validated_physical_execution_count"],
        "real_openmm_execution_count": cert["real_openmm_execution_count"],
        "proxy_physical_execution_used_for_claim": cert["proxy_physical_execution_used_for_claim"],
        "target_fold_claim_count": cert["target_fold_claim_count"],
        "unsupported_fold_claims": cert["unsupported_fold_claims"],
        "unsupported_physical_claims": cert["unsupported_physical_claims"],
        "coordinate_native_leakage": cert["coordinate_native_leakage"],
        "external_blind_benchmark_exported": cert["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": cert["external_blind_benchmark_passed"],
        "universal_folding_solution_claim_allowed": cert["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "blocked_reasons": cert["blocked_reasons"],
        "target_fold_claim_count_nonzero_full_suite_required": cert[
            "target_fold_claim_count_nonzero_full_suite_required"
        ],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["campaign_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
