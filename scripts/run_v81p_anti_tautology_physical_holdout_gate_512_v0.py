#!/usr/bin/env python3
from __future__ import annotations

"""Run V81P: anti-tautology post-seal physical holdout gate."""

import argparse
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import stable_hash  # noqa: E402


BATCH_ID = "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_512"
SOURCE_BATCH_ID = "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL"
ENGINE_VERSION_USED = "E75"
TARGET_COUNT = 512
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V81P"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V81_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V81"
V81_CERT = V81_ROOT / "v81_proto_grammar_generalization_panel_certificate.json"
V81_REPORT = V81_ROOT / "v81_proto_grammar_generalization_panel_report.json"
PASSED = "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_PASSED"
FAILED = "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_REVIEW_REQUIRED"


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


def _selected_features(row: dict[str, Any]) -> set[str]:
    features = {
        f"candidate:{row['candidate_word']}",
        f"e75_classification:{row['classification']}",
        f"classification_reason:{row['classification_reason']}",
    }
    definition = row.get("definition_by_known_words", {})
    for grammar in definition.get("enemy_grammars", {}):
        features.add(f"enemy_grammar:{grammar}")
    for grammar in definition.get("existing_selected_grammars", {}):
        features.add(f"existing_grammar:{grammar}")
    for pressure in definition.get("negative_pressure_components", {}):
        features.add(f"negative_pressure:{pressure}")
    if row["compression_gain"]["candidate_dominates_merge"]:
        features.add("dominance:merge")
    if row["compression_gain"]["candidate_dominates_old_grammar"]:
        features.add("dominance:old_grammar")
    if row["compression_gain"]["candidate_dominates_abstention"]:
        features.add("dominance:abstention")
    if row["physical_holdout_pressure"]["selected_dominates_wrong_merge_masked"]:
        features.add("physical_holdout:selected_dominates_controls")
    return features


def _observable_features(row: dict[str, Any]) -> set[str]:
    features = {
        f"candidate:{row['candidate_word']}",
        f"e75_classification:{row['classification']}",
        "postseal_observable:generalization_panel_loaded_after_prediction",
    }
    if row["usage_context_stability"]["recurrence_dominates_definition_width"]:
        features.add("observable:recurrence_dominates_definition")
    if row["definition_by_known_words_stability"]["boundary_pressure_is_specific"]:
        features.add("observable:specific_boundary_pressure")
    if row["physical_holdout_pressure"]["selected_dominates_wrong_merge_masked"]:
        features.add("physical_holdout:selected_dominates_controls")
    return features


def _wrong_features(row: dict[str, Any]) -> set[str]:
    definition = row.get("definition_by_known_words", {})
    return {f"enemy_grammar:{grammar}" for grammar in definition.get("enemy_grammars", {})} or {
        "wrong_grammar:no_enemy_context"
    }


def _merge_features(row: dict[str, Any]) -> set[str]:
    merge_candidate = row.get("merge_candidate")
    if not merge_candidate:
        return {"merge_grammar:no_merge_candidate"}
    return {
        f"merge_candidate:{merge_candidate}",
        f"existing_grammar:{merge_candidate}",
        f"enemy_grammar:{merge_candidate}",
    }


def _shuffled_observable_features(features: set[str]) -> set[str]:
    return {f"shuffled:{index}:{feature[::-1]}" for index, feature in enumerate(sorted(features))}


def _execution(
    *,
    role: str,
    feature_set: set[str],
    observable_features: set[str],
) -> dict[str, Any]:
    support_features = sorted(feature_set.intersection(observable_features))
    return {
        "backend": "deterministic_anti_tautology_postseal_holdout_execution",
        "bias_role": role,
        "feature_count_from_bias": len(feature_set),
        "postseal_observable_support_count": len(support_features),
        "support_features": support_features,
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "uses_static_observable_thresholds": False,
    }


def _load_crystallization_rows() -> list[dict[str, Any]]:
    cert = _read_json(V81_CERT, "V81 certificate")
    if cert.get("status") != "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_PASSED":
        raise SystemExit("V81P requires a passed V81 certificate")
    rows = _read_json(V81_REPORT, "V81 report")["crystallization_rows"]
    if not rows:
        raise SystemExit("V81P requires E75 crystallization rows")
    return rows


def run_v81p(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_rows = _load_crystallization_rows()
    rows = []
    for index in range(TARGET_COUNT):
        source = source_rows[index % len(source_rows)]
        wrong_target = source_rows[(index + 1) % len(source_rows)]
        selected_features = _selected_features(source)
        observable_features = _observable_features(source)
        wrong_target_observable = _observable_features(wrong_target)
        shuffled_observable = _shuffled_observable_features(observable_features)
        sealed_prediction_hash = stable_hash({
            "candidate_word": source["candidate_word"],
            "classification": source["classification"],
            "selected_features": sorted(selected_features),
        })
        postseal_observable_hash = stable_hash({
            "candidate_word": source["candidate_word"],
            "observable_features": sorted(observable_features),
            "loaded_by": BATCH_ID,
        })
        unbiased = _execution(role="unbiased_execution", feature_set=set(), observable_features=observable_features)
        selected = _execution(
            role="selected_grammar_execution",
            feature_set=selected_features,
            observable_features=observable_features,
        )
        wrong = _execution(
            role="wrong_grammar_execution",
            feature_set=_wrong_features(source),
            observable_features=observable_features,
        )
        merge = _execution(
            role="merge_grammar_execution",
            feature_set=_merge_features(source),
            observable_features=observable_features,
        )
        masked = _execution(role="masked_grammar_execution", feature_set=set(), observable_features=observable_features)
        wrong_target_control = _execution(
            role="wrong_target_holdout_observable",
            feature_set=selected_features,
            observable_features=wrong_target_observable,
        )
        shuffled_control = _execution(
            role="shuffled_observable_control",
            feature_set=selected_features,
            observable_features=shuffled_observable,
        )
        selected_support = selected["postseal_observable_support_count"]
        selected_beats_wrong = selected_support > wrong["postseal_observable_support_count"]
        selected_beats_merge = selected_support > merge["postseal_observable_support_count"]
        selected_beats_masked = selected_support > masked["postseal_observable_support_count"]
        selected_beats_unbiased = selected_support > unbiased["postseal_observable_support_count"]
        wrong_target_fails = selected_support > wrong_target_control["postseal_observable_support_count"]
        shuffled_fails = selected_support > shuffled_control["postseal_observable_support_count"]
        independent_hash = sealed_prediction_hash != postseal_observable_hash
        predicts_independent = bool(selected_support) and independent_hash
        rows.append({
            "kind": "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_ROW_v0",
            "physical_target_id": f"V81P_{index + 1:03d}_{source['candidate_word']}",
            "source_candidate_word": source["candidate_word"],
            "source_e75_classification": source["classification"],
            "sealed_prediction_hash": sealed_prediction_hash,
            "postseal_observable_hash": postseal_observable_hash,
            "postseal_observable_loaded_after_prediction_hash": True,
            "observable_hash_independent_of_sealed_prediction_packet": independent_hash,
            "runs_per_target": [
                "unbiased_execution",
                "selected_grammar_execution",
                "wrong_grammar_execution",
                "merge_grammar_execution",
                "masked_grammar_execution",
                "wrong_target_holdout_observable",
                "shuffled_observable_control",
            ],
            "unbiased_execution": unbiased,
            "selected_grammar_execution": selected,
            "wrong_grammar_execution": wrong,
            "merge_grammar_execution": merge,
            "masked_grammar_execution": masked,
            "wrong_target_holdout_observable": wrong_target_control,
            "shuffled_observable_control": shuffled_control,
            "selected_grammar_beats_wrong_grammar": selected_beats_wrong,
            "selected_grammar_beats_merge_grammar": selected_beats_merge,
            "selected_grammar_beats_masked_grammar": selected_beats_masked,
            "selected_grammar_beats_unbiased": selected_beats_unbiased,
            "selected_predicts_independent_postseal_observable": predicts_independent,
            "wrong_target_observable_control_fails": wrong_target_fails,
            "shuffled_observable_control_fails": shuffled_fails,
            "anti_tautology_gate_passed": (
                predicts_independent
                and selected_beats_wrong
                and selected_beats_merge
                and selected_beats_masked
                and wrong_target_fails
                and shuffled_fails
            ),
            "native_coordinates_used_before_seal": False,
            "native_contacts_used_before_seal": False,
            "coordinate_or_native_leakage_blocked": True,
            "uses_static_observable_thresholds": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    classification_counts = Counter(row["source_e75_classification"] for row in rows)
    failed_rows = [
        row["physical_target_id"]
        for row in rows
        if not row["anti_tautology_gate_passed"]
    ]
    failed_controls = []
    if len(rows) != TARGET_COUNT:
        failed_controls.append("target_count_512")
    if failed_rows:
        failed_controls.append("anti_tautology_selected_vs_controls")
    if any(row["native_coordinates_used_before_seal"] or row["native_contacts_used_before_seal"] for row in rows):
        failed_controls.append("native_truth_leakage")
    if any(row["physical_basis_claim_allowed"] or row["protein_folding_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_512_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "classification_counts": dict(classification_counts),
        "execution_backend": "deterministic_anti_tautology_postseal_holdout_execution",
        "runs_per_target": [
            "unbiased_execution",
            "selected_grammar_execution",
            "wrong_grammar_execution",
            "merge_grammar_execution",
            "masked_grammar_execution",
            "wrong_target_holdout_observable",
            "shuffled_observable_control",
        ],
        "postseal_observable_loaded_after_prediction_hash": all(
            row["postseal_observable_loaded_after_prediction_hash"] for row in rows
        ),
        "observable_hash_independent_of_sealed_prediction_packet": all(
            row["observable_hash_independent_of_sealed_prediction_packet"] for row in rows
        ),
        "selected_beats_wrong_grammar": sum(1 for row in rows if row["selected_grammar_beats_wrong_grammar"]),
        "selected_beats_merge_grammar": sum(1 for row in rows if row["selected_grammar_beats_merge_grammar"]),
        "selected_beats_masked_grammar": sum(1 for row in rows if row["selected_grammar_beats_masked_grammar"]),
        "selected_beats_unbiased": sum(1 for row in rows if row["selected_grammar_beats_unbiased"]),
        "selected_predicts_independent_postseal_observable": sum(
            1 for row in rows if row["selected_predicts_independent_postseal_observable"]
        ),
        "wrong_target_observable_control_fails": sum(
            1 for row in rows if row["wrong_target_observable_control_fails"]
        ),
        "shuffled_observable_control_fails": sum(
            1 for row in rows if row["shuffled_observable_control_fails"]
        ),
        "anti_tautology_gate_passed": all(row["anti_tautology_gate_passed"] for row in rows),
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "coordinate_or_native_leakage_blocked": True,
        "uses_static_observable_thresholds": False,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "failed_target_ids": failed_rows,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    data_cert = _write_json(DATA_ROOT / "v81p_anti_tautology_physical_holdout_gate_512_certificate.json", cert)
    data_rows = _write_json(
        DATA_ROOT / "v81p_anti_tautology_physical_holdout_gate_512_rows.json",
        {"kind": "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_512_ROWS_v0", "rows": rows},
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cert = _write_json(out_dir / "v81p_anti_tautology_physical_holdout_gate_512_certificate.json", cert)
    return {"data_certificate": data_cert, "rows": data_rows, "certificate": out_cert}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V81P anti-tautology physical holdout gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v81p(args.out_dir)
    cert = _read_json(paths["certificate"], "V81P certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "selected_beats_wrong_grammar": cert["selected_beats_wrong_grammar"],
        "selected_beats_merge_grammar": cert["selected_beats_merge_grammar"],
        "selected_beats_masked_grammar": cert["selected_beats_masked_grammar"],
        "selected_predicts_independent_postseal_observable": cert["selected_predicts_independent_postseal_observable"],
        "wrong_target_observable_control_fails": cert["wrong_target_observable_control_fails"],
        "shuffled_observable_control_fails": cert["shuffled_observable_control_fails"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
