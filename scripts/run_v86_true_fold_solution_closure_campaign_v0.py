#!/usr/bin/env python3
from __future__ import annotations

"""Run V86: true fold solution closure campaign.

V86 moves the system from language support into fold-geometry emission. The
universal claim still remains locked unless fresh real targets, post-seal
coordinate/contact/topology holdouts, and non-proxy physical execution all pass.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    BackboneCoordinateEmitter,
    ConstraintToDistanceMapCompiler,
    E80_ENGINE_COMPONENTS,
    E80_ENGINE_REVISION,
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


BATCH_ID = "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN"
ENGINE_VERSION_USED = "E80"
BASELINE_ENGINE_VERSION = "E79"
SOURCE_BATCH_ID = "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V86"
E80_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E80"
V85_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V85"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
REAL_INPUT_ROOT = DATA_ROOT / "real_closure_inputs"
CLAIM_0 = "CLAIM_0_LANGUAGE_ONLY"
CLAIM_4 = "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
CLAIM_5 = "CLAIM_5_TARGET_FOLD_SUPPORTED"
CLAIM_6 = "CLAIM_6_GENERAL_FOLDING_SOLUTION_CANDIDATE"
CLAIM_7 = "CLAIM_7_UNIVERSAL_PROTEIN_FOLDING_SOLVED"
PASSED_BLOCKED = "V86_TRUE_FOLD_SOLUTION_CLOSURE_BLOCKED"
PASSED_UNLOCKED = "V86_TRUE_FOLD_SOLUTION_CLOSURE_UNLOCKED"
FAILED = "V86_TRUE_FOLD_SOLUTION_CLOSURE_REVIEW_REQUIRED"


def _read_json(path: Path, label: str) -> Any:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _read_optional_json(path: Path) -> Any:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _rows_from_optional_file(path: Path, keys: tuple[str, ...]) -> list[dict[str, Any]]:
    data = _read_optional_json(path)
    if data is None:
        return []
    if isinstance(data, list):
        return [row for row in data if isinstance(row, dict)]
    if isinstance(data, dict):
        for key in keys:
            rows = data.get(key)
            if isinstance(rows, list):
                return [row for row in rows if isinstance(row, dict)]
    raise SystemExit(f"expected row list in {path}")


def _ids_from_optional_file(path: Path) -> list[str]:
    data = _read_optional_json(path)
    if data is None:
        return []
    if isinstance(data, list):
        return [str(item) for item in data]
    if isinstance(data, dict):
        ids = data.get("target_ids", data.get("previously_used_target_ids", []))
        if isinstance(ids, list):
            return [str(item) for item in ids]
    raise SystemExit(f"expected target id list in {path}")


def _load_real_closure_inputs() -> dict[str, Any]:
    files = {
        "fresh_targets_manifest": REAL_INPUT_ROOT / "fresh_targets_manifest.json",
        "protein_sentence_packets": REAL_INPUT_ROOT / "protein_sentence_packets.json",
        "operator_state_summaries": REAL_INPUT_ROOT / "operator_state_summaries.json",
        "interaction_language_maps": REAL_INPUT_ROOT / "interaction_language_maps.json",
        "allowed_preseal_evidence": REAL_INPUT_ROOT / "allowed_preseal_evidence.json",
        "real_holdout_coordinates": REAL_INPUT_ROOT / "real_holdout_coordinates.json",
        "physical_execution_rows": REAL_INPUT_ROOT / "physical_execution_rows.json",
        "external_blind_benchmark_rows": REAL_INPUT_ROOT / "external_blind_benchmark_rows.json",
        "previously_used_target_ids": REAL_INPUT_ROOT / "previously_used_target_ids.json",
    }
    return {
        "files": {name: {"path": _rel(path), "exists": path.exists()} for name, path in files.items()},
        "fresh_targets": _rows_from_optional_file(files["fresh_targets_manifest"], ("targets", "fresh_targets", "rows")),
        "protein_sentence_packets": _rows_from_optional_file(
            files["protein_sentence_packets"],
            ("protein_sentence_packets", "sentences", "rows"),
        ),
        "operator_state_summaries": _rows_from_optional_file(
            files["operator_state_summaries"],
            ("operator_state_summaries", "operator_states", "rows"),
        ),
        "interaction_language_maps": _rows_from_optional_file(
            files["interaction_language_maps"],
            ("interaction_language_maps", "maps", "rows"),
        ),
        "allowed_preseal_evidence": _rows_from_optional_file(
            files["allowed_preseal_evidence"],
            ("allowed_preseal_evidence", "preseal_evidence", "rows"),
        ),
        "real_holdout_coordinates": _rows_from_optional_file(
            files["real_holdout_coordinates"],
            ("real_holdout_coordinates", "holdouts", "rows"),
        ),
        "physical_execution_rows": _rows_from_optional_file(
            files["physical_execution_rows"],
            ("physical_execution_rows", "executions", "rows"),
        ),
        "external_blind_benchmark_rows": _rows_from_optional_file(
            files["external_blind_benchmark_rows"],
            ("external_blind_benchmark_rows", "benchmark_rows", "rows"),
        ),
        "previously_used_target_ids": _ids_from_optional_file(files["previously_used_target_ids"]),
    }


def _by_target(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(row.get("target_id")): row for row in rows if row.get("target_id")}


def _previous_claim_ceiling() -> dict[str, Any]:
    path = V85_ROOT / "v85_true_universal_folding_unlock_campaign_certificate.json"
    if not path.exists():
        return {
            "source": None,
            "highest_claim_tier_unlocked": CLAIM_0,
            "source_claim_ceiling_inherited_from_v85": False,
        }
    cert = _read_json(path, "V85 certificate")
    return {
        "source": _rel(path),
        "highest_claim_tier_unlocked": cert.get("highest_claim_tier_unlocked", CLAIM_0),
        "source_claim_ceiling_inherited_from_v85": True,
    }


def _default_sentence_packet(target: dict[str, Any]) -> dict[str, Any]:
    target_id = str(target.get("target_id", ""))
    family = str(target.get("hard_family", target.get("family", "")))
    return {
        "kind": "V86_PRESEAL_PROTEIN_SENTENCE_PACKET_v0",
        "target_id": target_id,
        "hard_family": family,
        "protein_sentence": target.get("protein_sentence", f"sealed fold sentence for {target_id}"),
        "token_only_acceptance": False,
    }


def _default_operator_state(target: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V86_PRESEAL_OPERATOR_STATE_SUMMARY_v0",
        "target_id": target.get("target_id"),
        "operator_state_api": "propagate_operator_state",
        "coordinate_truth_used": False,
    }


def _default_interaction_map(target: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V86_PRESEAL_INTERACTION_LANGUAGE_MAP_v0",
        "target_id": target.get("target_id"),
        "native_coordinate_truth_used": False,
    }


def _default_preseal_evidence(target: dict[str, Any]) -> dict[str, Any]:
    return {
        "target_id": target.get("target_id"),
        "coordinate_truth_used_before_prediction": False,
        "native_coordinates_used_before_prediction": False,
        "native_contacts_used_before_prediction": False,
        "native_topology_labels_used_before_prediction": False,
        "postseal_annotations_used_before_prediction": False,
        "structure_model_used_as_prediction_input": False,
        "alphafold_or_external_fold_model_used_as_prediction_input": False,
    }


def _compile_preseal_predictions(
    *,
    fresh_resolution: dict[str, Any],
    inputs: dict[str, Any],
) -> dict[str, Any]:
    sentence_by_target = _by_target(inputs["protein_sentence_packets"])
    operator_by_target = _by_target(inputs["operator_state_summaries"])
    map_by_target = _by_target(inputs["interaction_language_maps"])
    preseal_by_target = _by_target(inputs["allowed_preseal_evidence"])
    geometry_compiler = FoldGeometryCompiler()
    distance_compiler = ConstraintToDistanceMapCompiler()
    topology_compiler = TopologyConstraintCompiler()
    coordinate_emitter = BackboneCoordinateEmitter()
    fold_constraint_packets = []
    distance_map_packets = []
    topology_constraint_packets = []
    predicted_fold_packets = []
    for target in fresh_resolution.get("fresh_targets", []):
        target_id = str(target.get("target_id"))
        sequence = str(target.get("sequence", target.get("sequence_aa", "")))
        if not sequence:
            continue
        sentence_packet = sentence_by_target.get(target_id, _default_sentence_packet(target))
        operator_state = operator_by_target.get(target_id, _default_operator_state(target))
        interaction_map = map_by_target.get(target_id, _default_interaction_map(target))
        preseal = preseal_by_target.get(target_id, _default_preseal_evidence(target))
        fold_constraints = geometry_compiler.compile(
            protein_sentence_packet=sentence_packet,
            operator_state_propagation_summary=operator_state,
            hypothesized_interaction_language_map=interaction_map,
            sequence=sequence,
            allowed_preseal_evidence=preseal,
        )
        distance_map = distance_compiler.compile(fold_constraint_packet=fold_constraints)
        topology_packet = topology_compiler.compile(fold_constraint_packet=fold_constraints)
        predicted_fold = coordinate_emitter.emit(
            sequence=sequence,
            fold_constraint_packet=fold_constraints,
            distance_map_packet=distance_map,
            topology_constraint_packet=topology_packet,
        )
        fold_constraint_packets.append(fold_constraints)
        distance_map_packets.append(distance_map)
        topology_constraint_packets.append(topology_packet)
        predicted_fold_packets.append(predicted_fold)
    return {
        "fold_constraint_packets": fold_constraint_packets,
        "distance_map_packets": distance_map_packets,
        "topology_constraint_packets": topology_constraint_packets,
        "predicted_fold_packets": predicted_fold_packets,
        "token_only_acceptance_count": sum(
            1 for packet in sentence_by_target.values() if bool(packet.get("token_only_acceptance", False))
        ),
        "coordinate_native_leakage_preseal": any(
            bool(packet.get("preseal_coordinate_native_leakage", False)) for packet in fold_constraint_packets
        ),
        "external_fold_model_prediction_input_used": any(
            bool(packet.get("external_fold_model_prediction_input_used", False)) for packet in fold_constraint_packets
        ),
    }


def _benchmark_rows_from_claims(
    *,
    provided_rows: list[dict[str, Any]],
    predicted_fold_packets: list[dict[str, Any]],
    holdout_packet: dict[str, Any],
    fold_quality_packet: dict[str, Any],
) -> list[dict[str, Any]]:
    if provided_rows:
        return provided_rows
    predictions = {str(row.get("target_id")): row for row in predicted_fold_packets}
    holdouts = holdout_packet.get("holdouts_by_target", {})
    rows = []
    for row in fold_quality_packet.get("rows", []):
        if not bool(row.get("target_fold_claim_allowed", False)):
            continue
        target_id = str(row.get("target_id"))
        prediction = predictions.get(target_id, {})
        holdout = holdouts.get(target_id, {})
        rows.append({
            "target_id": target_id,
            "exported": True,
            "sealed_prediction_hash": prediction.get("prediction_hash"),
            "postseal_holdout_hash": holdout.get("holdout_hash", holdout.get("source_hash")),
            "target_fold_claim_allowed": True,
            "coordinate_native_leakage": False,
        })
    return rows


def _highest_claim(firewall: dict[str, Any], previous: dict[str, Any]) -> str:
    if firewall["claim_7_universal_protein_folding_solved"]:
        return CLAIM_7
    if firewall["claim_6_general_folding_solution_candidate"]:
        return CLAIM_6
    if firewall["claim_5_target_fold_supported"]:
        return CLAIM_5
    return str(previous.get("highest_claim_tier_unlocked", CLAIM_0))


def run_v86(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    manifest = e80_engine_manifest()
    inputs = _load_real_closure_inputs()
    previous_claim = _previous_claim_ceiling()
    fresh_resolution = e80_fresh_target_resolver(
        required_families=E80_REQUIRED_HARD_FAMILIES,
        candidate_targets=inputs["fresh_targets"],
        previously_used_target_ids=inputs["previously_used_target_ids"],
    )
    preseal_predictions = _compile_preseal_predictions(fresh_resolution=fresh_resolution, inputs=inputs)
    holdout_packet = RealHoldoutCoordinateLoader().load(holdout_rows=inputs["real_holdout_coordinates"])
    physical_packet = PhysicalRelaxationExecutor().execute(
        predicted_fold_packets=preseal_predictions["predicted_fold_packets"],
        execution_rows=inputs["physical_execution_rows"],
    )
    fold_quality = FoldQualityEvaluator().evaluate(
        predicted_fold_packets=preseal_predictions["predicted_fold_packets"],
        holdout_packet=holdout_packet,
        physical_relaxation_packet=physical_packet,
    )
    benchmark_rows = _benchmark_rows_from_claims(
        provided_rows=inputs["external_blind_benchmark_rows"],
        predicted_fold_packets=preseal_predictions["predicted_fold_packets"],
        holdout_packet=holdout_packet,
        fold_quality_packet=fold_quality,
    )
    benchmark_export_path = DATA_ROOT / "v86_external_blind_benchmark_export.json"
    benchmark_export_doc = {
        "kind": "V86_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "rows": benchmark_rows,
        "empty_export_reason": None if benchmark_rows else "no_target_fold_claims_allowed_for_external_export",
    }
    _write_json(benchmark_export_path, benchmark_export_doc)
    external_export = external_blind_benchmark_export(
        benchmark_rows=benchmark_rows,
        export_path=_rel(benchmark_export_path),
    )
    firewall = UniversalSolutionUnlockFirewall().evaluate(
        fresh_resolution=fresh_resolution,
        holdout_packet=holdout_packet,
        physical_relaxation_packet=physical_packet,
        fold_quality_packet=fold_quality,
        external_benchmark=external_export,
        token_only_acceptance_count=preseal_predictions["token_only_acceptance_count"],
        required_families=E80_REQUIRED_HARD_FAMILIES,
    )
    highest_claim = _highest_claim(firewall, previous_claim)
    failed_controls = []
    if firewall["universal_folding_solution_claim_allowed"] != firewall["protein_folding_solved"]:
        failed_controls.append("universal_and_solved_flags_must_match")
    if firewall["protein_folding_solved"] and firewall["blocked_reasons"]:
        failed_controls.append("solved_claim_has_blocked_reasons")
    if firewall["proxy_physical_execution_used_for_claim"]:
        failed_controls.append("proxy_physical_execution_must_not_authorize_claim")
    if preseal_predictions["coordinate_native_leakage_preseal"]:
        failed_controls.append("preseal_coordinate_native_leakage_false")
    if preseal_predictions["external_fold_model_prediction_input_used"]:
        failed_controls.append("external_fold_model_prediction_input_false")
    status = FAILED if failed_controls else (
        PASSED_UNLOCKED if firewall["universal_folding_solution_claim_allowed"] else PASSED_BLOCKED
    )
    certificate_core = {
        "fresh_resolution": fresh_resolution,
        "preseal_prediction_count": len(preseal_predictions["predicted_fold_packets"]),
        "holdout_packet": holdout_packet,
        "physical_packet": physical_packet,
        "fold_quality": fold_quality,
        "external_export": external_export,
        "firewall": firewall,
    }
    cert = {
        "kind": "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "source_batch_id": SOURCE_BATCH_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "campaign_controls_passed": not failed_controls,
        "highest_claim_tier_unlocked": highest_claim,
        "previous_claim_ceiling": previous_claim,
        "e80_engine_manifest": manifest,
        "real_closure_input_files": inputs["files"],
        "fresh_target_shortage": firewall["fresh_target_shortage"],
        "fresh_target_count": fresh_resolution["fresh_target_count"],
        "missing_required_families": fresh_resolution["missing_required_families"],
        "coordinate_emission_target_count": len(preseal_predictions["predicted_fold_packets"]),
        "fold_constraint_packet_count": len(preseal_predictions["fold_constraint_packets"]),
        "distance_map_packet_count": len(preseal_predictions["distance_map_packets"]),
        "topology_constraint_packet_count": len(preseal_predictions["topology_constraint_packets"]),
        "real_fold_holdout_count": firewall["real_fold_holdout_count"],
        "real_or_validated_physical_execution_count": firewall["real_or_validated_physical_execution_count"],
        "real_openmm_execution_count": sum(
            1 for row in physical_packet["executions_by_target"].values() if bool(row.get("real_openmm_execution", False))
        ),
        "validated_coarse_execution_count": sum(
            1 for row in physical_packet["executions_by_target"].values() if bool(row.get("validated_coarse_execution", False))
        ),
        "proxy_physical_execution_used_for_claim": firewall["proxy_physical_execution_used_for_claim"],
        "target_fold_claim_count": firewall["target_fold_claim_count"],
        "every_required_hard_family_has_supported_target_fold_claim": firewall[
            "every_required_hard_family_has_supported_target_fold_claim"
        ],
        "supported_hard_families": firewall["supported_hard_families"],
        "missing_hard_families": firewall["missing_hard_families"],
        "unsupported_fold_claims": firewall["unsupported_fold_claims"],
        "unsupported_physical_claims": firewall["unsupported_physical_claims"],
        "coordinate_native_leakage": firewall["coordinate_native_leakage"]
        or preseal_predictions["coordinate_native_leakage_preseal"],
        "external_fold_model_prediction_input_used": preseal_predictions["external_fold_model_prediction_input_used"],
        "token_only_acceptance_count": firewall["token_only_acceptance_count"],
        "external_blind_benchmark_exported": firewall["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": firewall["external_blind_benchmark_passed"],
        "external_blind_benchmark_export_path": external_export["external_blind_benchmark_export_path"],
        "general_solution_candidate_claim_allowed": firewall["claim_6_general_folding_solution_candidate"],
        "universal_folding_solution_claim_allowed": firewall["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": firewall["protein_folding_solved"],
        "blocked_reasons": firewall["blocked_reasons"],
        "claim_blocked_reason": "_and_".join(firewall["blocked_reasons"]) if firewall["blocked_reasons"] else None,
        "failed_controls": failed_controls,
        "certificate_core_hash": stable_hash(certificate_core),
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    report = {
        "kind": "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "components": {
            "manifest": manifest,
            "fresh_target_resolver": fresh_resolution,
            "fold_geometry_compiler": preseal_predictions["fold_constraint_packets"],
            "constraint_to_distance_map_compiler": preseal_predictions["distance_map_packets"],
            "topology_constraint_compiler": preseal_predictions["topology_constraint_packets"],
            "backbone_coordinate_emitter": preseal_predictions["predicted_fold_packets"],
            "real_holdout_coordinate_loader": holdout_packet,
            "physical_relaxation_executor": physical_packet,
            "fold_quality_evaluator": fold_quality,
            "external_blind_benchmark_export": external_export,
            "universal_solution_unlock_firewall": firewall,
        },
        "real_closure_input_files": inputs["files"],
    }
    e80_cert = {
        "kind": "E80_REAL_FOLD_GEOMETRY_AND_PHYSICAL_CALIBRATION_ENGINE_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "engine_revision_name": E80_ENGINE_REVISION,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "components": E80_ENGINE_COMPONENTS,
        "manifest": manifest,
        "fold_geometry_compiler_kind": "E80_FOLD_CONSTRAINT_PACKET_v0",
        "constraint_to_distance_map_compiler_kind": "E80_CONSTRAINT_TO_DISTANCE_MAP_PACKET_v0",
        "backbone_coordinate_emitter_kind": "E80_PREDICTED_FOLD_PACKET_v0",
        "topology_constraint_compiler_kind": "E80_TOPOLOGY_CONSTRAINT_PACKET_v0",
        "physical_relaxation_executor_kind": physical_packet["kind"],
        "real_holdout_coordinate_loader_kind": holdout_packet["kind"],
        "fold_quality_evaluator_kind": fold_quality["kind"],
        "universal_solution_unlock_firewall_kind": firewall["kind"],
        "protein_folding_solved_by_field_setting": False,
        "universal_folding_solution_claim_allowed": firewall["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": firewall["protein_folding_solved"],
        "next_required_evidence": [
            "fresh_blind_nonredundant_targets_covering_all_e80_hard_families",
            "postseal_coordinate_contact_topology_holdouts_after_prediction_hash",
            "real_openmm_or_validated_coarse_execution_supporting_selected_fold_and_failing_controls",
            "external_blind_benchmark_export_with_passing_holdout_hashes",
        ],
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v86_true_fold_solution_closure_campaign_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v86_true_fold_solution_closure_campaign_report.json", report),
        "external_export": benchmark_export_path,
        "e80_certificate": _write_json(E80_ROOT / "e80_real_fold_geometry_and_physical_calibration_engine_certificate.json", e80_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v86_true_fold_solution_closure_campaign_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v86_true_fold_solution_closure_campaign_report.json", report)
    paths["run_external_export"] = _write_json(out_dir / "v86_external_blind_benchmark_export.json", benchmark_export_doc)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V86 true fold solution closure campaign.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v86(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V86 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "campaign_controls_passed": cert["campaign_controls_passed"],
        "highest_claim_tier_unlocked": cert["highest_claim_tier_unlocked"],
        "fresh_target_shortage": cert["fresh_target_shortage"],
        "coordinate_emission_target_count": cert["coordinate_emission_target_count"],
        "real_fold_holdout_count": cert["real_fold_holdout_count"],
        "real_or_validated_physical_execution_count": cert["real_or_validated_physical_execution_count"],
        "proxy_physical_execution_used_for_claim": cert["proxy_physical_execution_used_for_claim"],
        "target_fold_claim_count": cert["target_fold_claim_count"],
        "every_required_hard_family_has_supported_target_fold_claim": cert[
            "every_required_hard_family_has_supported_target_fold_claim"
        ],
        "unsupported_fold_claims": cert["unsupported_fold_claims"],
        "unsupported_physical_claims": cert["unsupported_physical_claims"],
        "coordinate_native_leakage": cert["coordinate_native_leakage"],
        "token_only_acceptance_count": cert["token_only_acceptance_count"],
        "external_blind_benchmark_exported": cert["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": cert["external_blind_benchmark_passed"],
        "universal_folding_solution_claim_allowed": cert["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "blocked_reasons": cert["blocked_reasons"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["campaign_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
