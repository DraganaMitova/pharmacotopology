#!/usr/bin/env python3
from __future__ import annotations

"""Run V85: true universal folding unlock campaign.

The campaign is intentionally evidence-gated. It can unlock CLAIM_5/6/7 only
from fresh blind targets, post-seal fold holdouts, and real or validated
physical executions. Missing evidence is a successful block, not a failure.
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
    E79_ENGINE_REVISION,
    E79_HARD_REGIME_FAMILIES,
    E79_UNLOCK_COMPONENTS,
    atomistic_or_validated_physical_executor,
    e79_unlock_engine_manifest,
    e79_universal_claim_firewall,
    external_blind_benchmark_export,
    family_generalization_evaluator,
    fresh_target_resolver,
    real_fold_holdout_loader,
    stable_hash,
    target_fold_evaluator,
)


BATCH_ID = "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN"
ENGINE_VERSION_USED = "E79"
BASELINE_ENGINE_VERSION = "E78"
SOURCE_BATCH_ID = "V84_FINAL_CLOSED_LOOP_PROTEIN_FOLDING_CAMPAIGN"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V85"
E79_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E79"
V84_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V84"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
REAL_INPUT_ROOT = DATA_ROOT / "real_unlock_inputs"
CALIBRATION_ROOTS = [
    REPO_ROOT / "data" / "protein_esperanto_engine" / "V75" / "physical_calibration",
    REPO_ROOT / "data" / "protein_esperanto_engine" / "V76" / "physical_calibration",
    REPO_ROOT / "data" / "protein_esperanto_engine" / "V77" / "physical_calibration",
    REPO_ROOT / "data" / "protein_esperanto_engine" / "V78" / "physical_calibration",
]
CLAIM_4 = "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
CLAIM_5 = "CLAIM_5_TARGET_FOLD_SUPPORTED"
CLAIM_6 = "CLAIM_6_GENERAL_SOLUTION_CANDIDATE"
CLAIM_7 = "CLAIM_7_UNIVERSAL_PROTEIN_FOLDING_SOLVED"
CLAIM_0 = "CLAIM_0_LANGUAGE_ONLY"
PASSED_BLOCKED = "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_BLOCKED"
PASSED_UNLOCKED = "V85_TRUE_UNIVERSAL_FOLDING_UNLOCKED"
FAILED = "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_REVIEW_REQUIRED"


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
    raise SystemExit(f"expected rows list in {path}")


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


def _load_real_unlock_inputs() -> dict[str, Any]:
    files = {
        "fresh_targets_manifest": REAL_INPUT_ROOT / "fresh_targets_manifest.json",
        "sealed_predictions": REAL_INPUT_ROOT / "sealed_predictions.json",
        "real_fold_holdouts": REAL_INPUT_ROOT / "real_fold_holdouts.json",
        "real_physical_executions": REAL_INPUT_ROOT / "real_physical_executions.json",
        "sentinels": REAL_INPUT_ROOT / "sentinels.json",
        "external_blind_benchmark_rows": REAL_INPUT_ROOT / "external_blind_benchmark_rows.json",
        "previously_used_target_ids": REAL_INPUT_ROOT / "previously_used_target_ids.json",
        "unresolved_classes": REAL_INPUT_ROOT / "unresolved_classes.json",
    }
    return {
        "files": {name: {"path": _rel(path), "exists": path.exists()} for name, path in files.items()},
        "fresh_targets": _rows_from_optional_file(files["fresh_targets_manifest"], ("targets", "fresh_targets", "rows")),
        "sealed_predictions": _rows_from_optional_file(files["sealed_predictions"], ("predictions", "sealed_predictions", "rows")),
        "real_fold_holdouts": _rows_from_optional_file(files["real_fold_holdouts"], ("holdouts", "real_fold_holdouts", "rows")),
        "real_physical_executions": _rows_from_optional_file(
            files["real_physical_executions"],
            ("executions", "real_physical_executions", "rows"),
        ),
        "sentinels": _rows_from_optional_file(files["sentinels"], ("sentinels", "rows")),
        "external_blind_benchmark_rows": _rows_from_optional_file(
            files["external_blind_benchmark_rows"],
            ("benchmark_rows", "rows"),
        ),
        "previously_used_target_ids": _ids_from_optional_file(files["previously_used_target_ids"]),
        "unresolved_classes": _ids_from_optional_file(files["unresolved_classes"]),
    }


def _load_prior_calibration_inputs() -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for root in CALIBRATION_ROOTS:
        for path in sorted(root.glob("*.json")):
            data = _read_optional_json(path)
            if not isinstance(data, dict):
                continue
            summaries.append({
                "path": _rel(path),
                "kind": data.get("kind"),
                "calibration_input_type": data.get("calibration_input_type"),
                "row_count": data.get("row_count", len(data.get("calibration_rows", []))),
                "fold_class_coverage": data.get("fold_class_coverage", []),
                "coordinate_truth_used_as_prediction_input": bool(
                    data.get("coordinate_truth_used_as_prediction_input", True)
                ),
                "target_native_contacts_used_before_prediction": bool(
                    data.get("target_native_contacts_used_before_prediction", True)
                ),
                "universal_physical_law_claim_allowed": bool(data.get("universal_physical_law_claim_allowed", True)),
                "calibration_hash": data.get("calibration_hash"),
            })
    return summaries


def _previous_claim_ceiling() -> dict[str, Any]:
    path = V84_ROOT / "v84_final_closed_loop_protein_folding_campaign_certificate.json"
    if not path.exists():
        return {
            "source": None,
            "highest_claim_tier_unlocked": CLAIM_0,
            "source_claim_ceiling_inherited_from_v84": False,
        }
    cert = _read_json(path, "V84 certificate")
    ledger = cert.get("final_claim_ledger", {})
    return {
        "source": _rel(path),
        "highest_claim_tier_unlocked": ledger.get("highest_claim_tier_unlocked", cert.get("highest_claim_tier_unlocked", CLAIM_0)),
        "source_claim_ceiling_inherited_from_v84": True,
    }


def _sentinels_preserved(rows: list[dict[str, Any]]) -> bool:
    return bool(rows) and all(bool(row.get("preserved", row.get("sentinel_preserved", False))) for row in rows)


def _benchmark_rows_from_target_claims(
    *,
    provided_rows: list[dict[str, Any]],
    sealed_predictions: list[dict[str, Any]],
    fold_holdouts: dict[str, Any],
    target_fold_evaluation: dict[str, Any],
) -> list[dict[str, Any]]:
    if provided_rows:
        return provided_rows
    predictions = {str(row.get("target_id")): row for row in sealed_predictions}
    holdouts = fold_holdouts.get("holdouts_by_target", {})
    rows = []
    for row in target_fold_evaluation.get("rows", []):
        if not bool(row.get("target_fold_claim_allowed", False)):
            continue
        target_id = str(row.get("target_id"))
        prediction = predictions.get(target_id, {})
        holdout = holdouts.get(target_id, {})
        rows.append({
            "target_id": target_id,
            "exported": True,
            "sealed_prediction_hash": prediction.get("prediction_hash"),
            "postseal_holdout_hash": holdout.get("source_hash"),
            "target_fold_claim_allowed": True,
            "coordinate_native_leakage": False,
        })
    return rows


def _highest_claim(firewall: dict[str, Any], previous: dict[str, Any]) -> str:
    if firewall["claim_7_universal_protein_folding_solved"]:
        return CLAIM_7
    if firewall["claim_6_general_solution_candidate"]:
        return CLAIM_6
    if firewall["claim_5_target_fold_supported"]:
        return CLAIM_5
    return str(previous.get("highest_claim_tier_unlocked", CLAIM_0))


def run_v85(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    manifest = e79_unlock_engine_manifest()
    inputs = _load_real_unlock_inputs()
    calibration_inputs = _load_prior_calibration_inputs()
    previous_claim = _previous_claim_ceiling()
    fresh_resolution = fresh_target_resolver(
        required_families=E79_HARD_REGIME_FAMILIES,
        candidate_targets=inputs["fresh_targets"],
        previously_used_target_ids=inputs["previously_used_target_ids"],
    )
    fold_holdouts = real_fold_holdout_loader(holdout_rows=inputs["real_fold_holdouts"])
    physical_execution = atomistic_or_validated_physical_executor(execution_rows=inputs["real_physical_executions"])
    target_fold = target_fold_evaluator(
        fresh_resolution=fresh_resolution,
        sealed_predictions=inputs["sealed_predictions"],
        fold_holdouts=fold_holdouts,
        physical_execution=physical_execution,
    )
    unsupported_claims = target_fold["unsupported_fold_claims"] + physical_execution["unsupported_physical_claims"]
    family_generalization = family_generalization_evaluator(
        required_families=E79_HARD_REGIME_FAMILIES,
        target_fold_evaluation=target_fold,
        sentinels_preserved=_sentinels_preserved(inputs["sentinels"]),
        failed_accepted_count=0,
        unsupported_claims=unsupported_claims,
    )
    unresolved_classes = sorted(set(inputs["unresolved_classes"] + family_generalization["missing_families"]))
    benchmark_rows = _benchmark_rows_from_target_claims(
        provided_rows=inputs["external_blind_benchmark_rows"],
        sealed_predictions=inputs["sealed_predictions"],
        fold_holdouts=fold_holdouts,
        target_fold_evaluation=target_fold,
    )
    benchmark_export_path = DATA_ROOT / "v85_external_blind_benchmark_export.json"
    benchmark_export_doc = {
        "kind": "V85_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0",
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
    firewall = e79_universal_claim_firewall(
        fresh_resolution=fresh_resolution,
        physical_execution=physical_execution,
        target_fold_evaluation=target_fold,
        family_generalization=family_generalization,
        external_benchmark=external_export,
        unresolved_classes=unresolved_classes,
    )
    highest_claim = _highest_claim(firewall, previous_claim)
    unlocked = bool(firewall["claim_7_universal_protein_folding_solved"])
    failed_controls = []
    if firewall["universal_folding_solution_claim_allowed"] != firewall["protein_folding_solved"]:
        failed_controls.append("universal_and_solved_flags_must_match")
    if firewall["protein_folding_solved"] and firewall["blocked_reasons"]:
        failed_controls.append("solved_claim_has_blocked_reasons")
    if firewall["general_solution_candidate_claim_allowed"] and not firewall["claim_5_target_fold_supported"]:
        failed_controls.append("general_claim_requires_target_fold_support")
    if firewall["universal_folding_solution_claim_allowed"] and not external_export["external_blind_benchmark_passed"]:
        failed_controls.append("universal_claim_requires_external_blind_benchmark")
    if firewall["proxy_physical_execution_used_for_claim"]:
        failed_controls.append("proxy_physical_execution_must_not_authorize_claim")
    if target_fold["coordinate_native_leakage"]:
        failed_controls.append("coordinate_native_leakage_false")
    status = FAILED if failed_controls else (PASSED_UNLOCKED if unlocked else PASSED_BLOCKED)
    certificate_core = {
        "fresh_resolution": fresh_resolution,
        "physical_execution_summary": {
            "execution_count": physical_execution["execution_count"],
            "real_or_validated_execution_count": physical_execution["real_or_validated_execution_count"],
            "proxy_physical_execution_used_for_claim": physical_execution["proxy_physical_execution_used_for_claim"],
        },
        "target_fold_claim_count": target_fold["target_fold_claim_count"],
        "family_generalization": family_generalization,
        "external_export": external_export,
        "firewall": firewall,
    }
    cert = {
        "kind": "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "source_batch_id": SOURCE_BATCH_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "campaign_controls_passed": not failed_controls,
        "highest_claim_tier_unlocked": highest_claim,
        "previous_claim_ceiling": previous_claim,
        "e79_engine_manifest": manifest,
        "real_unlock_input_files": inputs["files"],
        "prior_calibration_inputs": calibration_inputs,
        "prior_calibration_inputs_count": len(calibration_inputs),
        "prior_calibration_inputs_used_for_universal_unlock": False,
        "fresh_target_shortage": firewall["fresh_target_shortage"],
        "fresh_target_count": fresh_resolution["fresh_target_count"],
        "missing_required_families": fresh_resolution["missing_required_families"],
        "deterministic_variant_count": fresh_resolution["deterministic_variant_count"],
        "deterministic_variants_for_universal_claim_allowed": False,
        "real_fold_holdout_count": fold_holdouts["loaded_holdout_count"],
        "real_or_validated_physical_execution_count": physical_execution["real_or_validated_execution_count"],
        "proxy_physical_execution_used_for_claim": firewall["proxy_physical_execution_used_for_claim"],
        "target_fold_claim_count": firewall["target_fold_claim_count"],
        "general_solution_candidate_claim_allowed": firewall["general_solution_candidate_claim_allowed"],
        "universal_folding_solution_claim_allowed": firewall["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": firewall["protein_folding_solved"],
        "unsupported_fold_claims": firewall["unsupported_fold_claims"],
        "unsupported_physical_claims": firewall["unsupported_physical_claims"],
        "unsupported_claims": firewall["unsupported_fold_claims"] + firewall["unsupported_physical_claims"],
        "coordinate_native_leakage": firewall["coordinate_native_leakage"],
        "external_blind_benchmark_exported": firewall["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": external_export["external_blind_benchmark_passed"],
        "external_blind_benchmark_export_path": external_export["external_blind_benchmark_export_path"],
        "sentinels_preserved": family_generalization["sentinels_preserved"],
        "failed_accepted_count": family_generalization["failed_accepted_count"],
        "unresolved_hard_classes": unresolved_classes,
        "blocked_reasons": firewall["blocked_reasons"],
        "claim_blocked_reason": "_and_".join(firewall["blocked_reasons"]) if firewall["blocked_reasons"] else None,
        "failed_controls": failed_controls,
        "certificate_core_hash": stable_hash(certificate_core),
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    report = {
        "kind": "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "components": {
            "manifest": manifest,
            "fresh_target_resolver": fresh_resolution,
            "real_fold_holdout_loader": fold_holdouts,
            "atomistic_or_validated_physical_executor": physical_execution,
            "target_fold_evaluator": target_fold,
            "family_generalization_evaluator": family_generalization,
            "external_blind_benchmark_export": external_export,
            "universal_claim_firewall": firewall,
        },
        "real_unlock_input_files": inputs["files"],
        "prior_calibration_inputs": calibration_inputs,
    }
    e79_cert = {
        "kind": "E79_UNIVERSAL_FOLDING_CLAIM_UNLOCK_ENGINE_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "engine_revision_name": E79_ENGINE_REVISION,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "components": E79_UNLOCK_COMPONENTS,
        "manifest": manifest,
        "universal_claim_firewall_kind": firewall["kind"],
        "real_fold_holdout_loader_kind": fold_holdouts["kind"],
        "fresh_target_resolver_kind": fresh_resolution["kind"],
        "atomistic_or_validated_physical_executor_kind": physical_execution["kind"],
        "target_fold_evaluator_kind": target_fold["kind"],
        "family_generalization_evaluator_kind": family_generalization["kind"],
        "external_blind_benchmark_export_kind": external_export["kind"],
        "universal_folding_solution_claim_allowed_by_field_setting": False,
        "protein_folding_solved_by_field_setting": False,
        "universal_folding_solution_claim_allowed": firewall["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": firewall["protein_folding_solved"],
        "next_required_evidence": [
            "fresh_nonredundant_blind_targets_covering_all_hard_regime_families",
            "postseal_coordinate_contact_topology_holdouts_for_selected_targets",
            "real_or_validated_physical_execution_rows_that_support_selected_and_fail_controls",
            "external_blind_benchmark_rows_exported_after_sealed_prediction_hashes",
        ],
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v85_true_universal_folding_unlock_campaign_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v85_true_universal_folding_unlock_campaign_report.json", report),
        "external_export": benchmark_export_path,
        "e79_certificate": _write_json(E79_ROOT / "e79_universal_folding_claim_unlock_engine_certificate.json", e79_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v85_true_universal_folding_unlock_campaign_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v85_true_universal_folding_unlock_campaign_report.json", report)
    paths["run_external_export"] = _write_json(out_dir / "v85_external_blind_benchmark_export.json", benchmark_export_doc)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V85 true universal folding unlock campaign.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v85(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V85 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "campaign_controls_passed": cert["campaign_controls_passed"],
        "highest_claim_tier_unlocked": cert["highest_claim_tier_unlocked"],
        "fresh_target_shortage": cert["fresh_target_shortage"],
        "real_fold_holdout_count": cert["real_fold_holdout_count"],
        "real_or_validated_physical_execution_count": cert["real_or_validated_physical_execution_count"],
        "proxy_physical_execution_used_for_claim": cert["proxy_physical_execution_used_for_claim"],
        "target_fold_claim_count": cert["target_fold_claim_count"],
        "general_solution_candidate_claim_allowed": cert["general_solution_candidate_claim_allowed"],
        "universal_folding_solution_claim_allowed": cert["universal_folding_solution_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "unsupported_fold_claims": cert["unsupported_fold_claims"],
        "unsupported_physical_claims": cert["unsupported_physical_claims"],
        "coordinate_native_leakage": cert["coordinate_native_leakage"],
        "external_blind_benchmark_exported": cert["external_blind_benchmark_exported"],
        "external_blind_benchmark_passed": cert["external_blind_benchmark_passed"],
        "blocked_reasons": cert["blocked_reasons"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["campaign_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
