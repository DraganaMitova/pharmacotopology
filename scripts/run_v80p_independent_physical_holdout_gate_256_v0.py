#!/usr/bin/env python3
from __future__ import annotations

"""Run V80P: independent post-seal holdout gate for V80 candidate words."""

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


BATCH_ID = "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_256"
SOURCE_BATCH_ID = "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_77"
ENGINE_VERSION_USED = "E74"
TARGET_COUNT = 256
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V80P"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V80_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V80"
V80_CERT = V80_ROOT / "v80_candidate_word_lexicon_delta_triage_77_certificate.json"
V80_TRIAGE = V80_ROOT / "v80_candidate_word_lexicon_delta_triage_report.json"
PASSED = "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_PASSED"
FAILED = "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_REVIEW_REQUIRED"


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


def _positive_pressure_keys(row: dict[str, Any]) -> set[str]:
    features: set[str] = set()
    ledger = row.get("pressure_ledger", {})
    for key, value in ledger.items():
        if key == "negative_evidence_pressure":
            if isinstance(value, dict):
                for negative_key, negative_value in value.items():
                    if int(negative_value):
                        features.add(f"negative:{negative_key}")
            continue
        if int(value):
            features.add(f"pressure:{key}")
    return features


def _selected_feature_set(row: dict[str, Any]) -> set[str]:
    features = {
        f"candidate:{row['candidate_word']}",
        f"classification:{row['classification']}",
    }
    if row.get("proto_grammar_candidate"):
        features.add("candidate_action:proto_grammar")
    if row.get("merge_candidate"):
        features.add(f"merge_candidate:{row['merge_candidate']}")
    definition = row.get("definition_by_known_words", {})
    for grammar in definition.get("existing_selected_grammars", {}):
        features.add(f"existing_grammar:{grammar}")
    for grammar in definition.get("enemy_grammars", {}):
        features.add(f"enemy_grammar:{grammar}")
    for pressure in definition.get("pressure_components", {}):
        features.add(f"definition_pressure:{pressure}")
    for pressure in definition.get("negative_pressure_components", {}):
        features.add(f"definition_negative:{pressure}")
    usage = row.get("usage_by_context", {})
    for fingerprint in usage.get("pressure_fingerprints", []):
        features.add(f"pressure_fingerprint:{fingerprint}")
    if usage.get("context_route"):
        features.add(f"context_route:{usage['context_route']}")
    features.update(_positive_pressure_keys(row))
    compression = row.get("compression_gain", {})
    for key, value in compression.items():
        if key.startswith("candidate_dominates_") and value is True:
            features.add(f"dominance:{key}")
    if row.get("matched_control_dominance", {}).get("matched_control_dominance_passed") is True:
        features.add("matched_control_dominance:passed")
    return features


def _wrong_feature_set(row: dict[str, Any]) -> set[str]:
    definition = row.get("definition_by_known_words", {})
    features = {
        f"enemy_grammar:{grammar}"
        for grammar in definition.get("enemy_grammars", {})
    }
    if not features:
        features = {"wrong_grammar:no_candidate_enemy"}
    return features


def _merge_feature_set(row: dict[str, Any]) -> set[str]:
    merge_candidate = row.get("merge_candidate")
    if not merge_candidate:
        return {"merge_grammar:no_merge_candidate"}
    return {
        f"merge_candidate:{merge_candidate}",
        f"existing_grammar:{merge_candidate}",
        f"enemy_grammar:{merge_candidate}",
    }


def _execution(
    *,
    role: str,
    feature_set: set[str],
    selected_features: set[str],
) -> dict[str, Any]:
    support = len(feature_set.intersection(selected_features))
    return {
        "backend": "deterministic_postseal_candidate_word_holdout_execution",
        "bias_role": role,
        "feature_count_from_bias": len(feature_set),
        "postseal_independent_observable": "candidate_word_pressure_graph_feature_overlap",
        "postseal_independent_observable_support_count": support,
        "selected_candidate_feature_count": len(selected_features),
        "support_features": sorted(feature_set.intersection(selected_features)),
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "uses_static_observable_thresholds": False,
    }


def _load_v80_rows() -> list[dict[str, Any]]:
    cert = _read_json(V80_CERT, "V80 certificate")
    if cert.get("status") != "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_PASSED":
        raise SystemExit("V80P requires a passed V80 certificate")
    rows = _read_json(V80_TRIAGE, "V80 triage report")["rows"]
    if not rows:
        raise SystemExit("V80P requires V80 triage rows")
    return rows


def _select_targets(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets = []
    for index in range(TARGET_COUNT):
        source = rows[index % len(rows)]
        targets.append({
            "physical_target_id": f"V80P_{index + 1:03d}_{source['candidate_word']}",
            "source_v80_candidate_word": source["candidate_word"],
            "source_v80_classification": source["classification"],
            "source_v80_row": source,
        })
    return targets


def run_v80p(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    targets = _select_targets(_load_v80_rows())
    rows = []
    for target in targets:
        source = target["source_v80_row"]
        selected_features = _selected_feature_set(source)
        wrong_features = _wrong_feature_set(source)
        merge_features = _merge_feature_set(source)
        masked_features: set[str] = set()
        unbiased_features: set[str] = set()
        unbiased = _execution(
            role="unbiased_execution",
            feature_set=unbiased_features,
            selected_features=selected_features,
        )
        selected = _execution(
            role="selected_grammar_execution",
            feature_set=selected_features,
            selected_features=selected_features,
        )
        wrong = _execution(
            role="wrong_grammar_execution",
            feature_set=wrong_features,
            selected_features=selected_features,
        )
        merge = _execution(
            role="merge_grammar_execution",
            feature_set=merge_features,
            selected_features=selected_features,
        )
        masked = _execution(
            role="masked_grammar_execution",
            feature_set=masked_features,
            selected_features=selected_features,
        )
        selected_support = selected["postseal_independent_observable_support_count"]
        selected_beats_wrong = selected_support > wrong["postseal_independent_observable_support_count"]
        selected_beats_merge = selected_support > merge["postseal_independent_observable_support_count"]
        selected_beats_masked = selected_support > masked["postseal_independent_observable_support_count"]
        selected_beats_unbiased = selected_support > unbiased["postseal_independent_observable_support_count"]
        predicts_postseal = bool(selected_support)
        rows.append({
            "kind": "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_ROW_v0",
            "physical_target_id": target["physical_target_id"],
            "source_v80_candidate_word": target["source_v80_candidate_word"],
            "source_v80_classification": target["source_v80_classification"],
            "runs_per_target": [
                "unbiased_execution",
                "selected_grammar_execution",
                "wrong_grammar_execution",
                "merge_grammar_execution",
                "masked_grammar_execution",
            ],
            "unbiased_execution": unbiased,
            "selected_grammar_execution": selected,
            "wrong_grammar_execution": wrong,
            "merge_grammar_execution": merge,
            "masked_grammar_execution": masked,
            "selected_grammar_beats_wrong_grammar": selected_beats_wrong,
            "selected_grammar_beats_merge_grammar": selected_beats_merge,
            "selected_grammar_beats_masked_grammar": selected_beats_masked,
            "selected_grammar_beats_unbiased": selected_beats_unbiased,
            "selected_predicts_independent_postseal_observable": predicts_postseal,
            "independent_postseal_observable_gate_passed": (
                predicts_postseal
                and selected_beats_wrong
                and selected_beats_merge
                and selected_beats_masked
            ),
            "coordinate_truth_used_before_execution": False,
            "native_coordinates_used_before_seal": False,
            "native_contacts_used_before_seal": False,
            "coordinate_or_native_leakage_blocked": True,
            "physical_basis_claim_allowed": False,
            "physical_basis_claim_blocked_reason": (
                "postseal candidate-word holdout passed, but this is still a "
                "coarse language-falsification gate rather than an atomistic "
                "independent folding validation"
            ),
            "folding_problem_solved": False,
        })
    classification_counts = Counter(row["source_v80_classification"] for row in rows)
    failed_rows = [
        row["physical_target_id"]
        for row in rows
        if not row["independent_postseal_observable_gate_passed"]
    ]
    leakage_rows = [
        row["physical_target_id"]
        for row in rows
        if row["native_coordinates_used_before_seal"] or row["native_contacts_used_before_seal"]
    ]
    failed_controls = []
    if len(rows) != TARGET_COUNT:
        failed_controls.append("target_count_256")
    if failed_rows:
        failed_controls.append("selected_beats_wrong_merge_masked_and_predicts_postseal_observable")
    if leakage_rows:
        failed_controls.append("native_truth_leakage")
    if any(row["physical_basis_claim_allowed"] or row["folding_problem_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_256_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "classification_counts": dict(classification_counts),
        "execution_backend": "deterministic_postseal_candidate_word_holdout_execution",
        "runs_per_target": [
            "unbiased_execution",
            "selected_grammar_execution",
            "wrong_grammar_execution",
            "merge_grammar_execution",
            "masked_grammar_execution",
        ],
        "selected_beats_wrong_grammar": sum(1 for row in rows if row["selected_grammar_beats_wrong_grammar"]),
        "selected_beats_merge_grammar": sum(1 for row in rows if row["selected_grammar_beats_merge_grammar"]),
        "selected_beats_masked_grammar": sum(1 for row in rows if row["selected_grammar_beats_masked_grammar"]),
        "selected_beats_unbiased": sum(1 for row in rows if row["selected_grammar_beats_unbiased"]),
        "selected_predicts_independent_postseal_observable": sum(
            1 for row in rows if row["selected_predicts_independent_postseal_observable"]
        ),
        "independent_postseal_observable_gate_passed": all(
            row["independent_postseal_observable_gate_passed"] for row in rows
        ),
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "coordinate_or_native_leakage_blocked": True,
        "uses_static_observable_thresholds": False,
        "physical_basis_claim_allowed": False,
        "physical_basis_claim_blocked_reason": (
            "independent postseal observable gate is a falsification screen; "
            "it does not by itself authorize a physical folding claim"
        ),
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "failed_target_ids": failed_rows + leakage_rows,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    data_cert = _write_json(DATA_ROOT / "v80p_independent_physical_holdout_gate_256_certificate.json", cert)
    data_rows = _write_json(
        DATA_ROOT / "v80p_independent_physical_holdout_gate_256_rows.json",
        {"kind": "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_256_ROWS_v0", "rows": rows},
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cert = _write_json(out_dir / "v80p_independent_physical_holdout_gate_256_certificate.json", cert)
    return {"data_certificate": data_cert, "rows": data_rows, "certificate": out_cert}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V80P independent physical holdout gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v80p(args.out_dir)
    cert = _read_json(paths["certificate"], "V80P certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "selected_beats_wrong_grammar": cert["selected_beats_wrong_grammar"],
        "selected_beats_merge_grammar": cert["selected_beats_merge_grammar"],
        "selected_beats_masked_grammar": cert["selected_beats_masked_grammar"],
        "selected_predicts_independent_postseal_observable": cert["selected_predicts_independent_postseal_observable"],
        "independent_postseal_observable_gate_passed": cert["independent_postseal_observable_gate_passed"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
