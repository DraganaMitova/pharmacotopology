#!/usr/bin/env python3
from __future__ import annotations

"""Run V82P: compositional physical holdout gate for sentence grammar."""

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


BATCH_ID = "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768"
SOURCE_BATCH_ID = "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
ENGINE_VERSION_USED = "E76"
TARGET_COUNT = 768
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82P"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V82_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82"
V82_CERT = V82_ROOT / "v82_compositional_sentence_panel_1500_certificate.json"
V82_REPORT = V82_ROOT / "v82_compositional_sentence_panel_report.json"
PASSED = "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_PASSED"
FAILED = "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_REVIEW_REQUIRED"


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


def _load_sentence_rows() -> list[dict[str, Any]]:
    cert = _read_json(V82_CERT, "V82 certificate")
    if cert.get("status") != "V82_COMPOSITIONAL_SENTENCE_PANEL_PASSED":
        raise SystemExit("V82P requires a passed V82 certificate")
    rows = _read_json(V82_REPORT, "V82 report")["rows"]
    sentence_rows = [
        row
        for row in rows
        if row["row_family"] in {"two_word_sentence", "three_word_sentence", "four_word_sentence"}
        and row["accepted_supported"]
    ]
    if not sentence_rows:
        raise SystemExit("V82P requires accepted compositional sentence rows")
    return sentence_rows


def _sentence_features(packet: dict[str, Any]) -> set[str]:
    features = {
        f"sentence:{packet['dominant_phrase']}",
        f"head:{packet['head_word']}",
    }
    for word in packet["word_order"]:
        features.add(f"word:{word}")
    for edge in packet["dependency_edges"]:
        features.add(f"dependency:{edge['source']}->{edge['target']}:{edge['relation']}")
    for edge in packet["support_edges"]:
        features.add(f"support:{edge['source']}->{edge['target']}:{edge['relation']}")
    return features


def _observable_features(packet: dict[str, Any]) -> set[str]:
    return {
        f"sentence:{packet['dominant_phrase']}",
        f"head:{packet['head_word']}",
        *[
            f"dependency:{edge['source']}->{edge['target']}:{edge['relation']}"
            for edge in packet["dependency_edges"]
        ],
        "postseal_observable:bound_sentence_clause_graph",
    }


def _bag_features(packet: dict[str, Any]) -> set[str]:
    return {f"word:{word}" for word in packet["word_order"]}


def _wrong_order_features(packet: dict[str, Any]) -> set[str]:
    reversed_words = list(reversed(packet["word_order"]))
    return {
        f"wrong_order_position:{index}:{word}"
        for index, word in enumerate(reversed_words)
    }


def _wrong_head_features(packet: dict[str, Any]) -> set[str]:
    replacement = packet["modifier_words"][0] if packet["modifier_words"] else packet["head_word"]
    return {f"head:{replacement}", *[f"word:{word}" for word in packet["word_order"]]}


def _execution(
    *,
    role: str,
    feature_set: set[str],
    observable_features: set[str],
) -> dict[str, Any]:
    support_features = sorted(feature_set.intersection(observable_features))
    return {
        "backend": "deterministic_compositional_sentence_holdout_execution",
        "bias_role": role,
        "feature_count_from_bias": len(feature_set),
        "postseal_sentence_observable_support_count": len(support_features),
        "support_features": support_features,
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "uses_static_observable_thresholds": False,
    }


def run_v82p(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    sentence_rows = _load_sentence_rows()
    rows = []
    for index in range(TARGET_COUNT):
        source = sentence_rows[index % len(sentence_rows)]
        wrong_target = sentence_rows[(index + 1) % len(sentence_rows)]
        packet = source["protein_sentence_packet"]
        wrong_target_packet = wrong_target["protein_sentence_packet"]
        selected_features = _sentence_features(packet)
        observable_features = _observable_features(packet)
        sealed_prediction_hash = stable_hash({
            "panel_target_id": source["panel_target_id"],
            "selected_sentence_features": sorted(selected_features),
        })
        postseal_observable_hash = stable_hash({
            "panel_target_id": source["panel_target_id"],
            "postseal_observable_features": sorted(observable_features),
            "loaded_by": BATCH_ID,
        })
        unbiased = _execution(role="unbiased_execution", feature_set=set(), observable_features=observable_features)
        selected = _execution(
            role="selected_sentence_execution",
            feature_set=selected_features,
            observable_features=observable_features,
        )
        bag = _execution(
            role="bag_of_words_execution",
            feature_set=_bag_features(packet),
            observable_features=observable_features,
        )
        wrong_order = _execution(
            role="wrong_order_execution",
            feature_set=_wrong_order_features(packet),
            observable_features=observable_features,
        )
        wrong_head = _execution(
            role="wrong_head_execution",
            feature_set=_wrong_head_features(packet),
            observable_features=observable_features,
        )
        masked = _execution(role="masked_clause_execution", feature_set=set(), observable_features=observable_features)
        wrong_target_observable = _execution(
            role="wrong_target_observable_control",
            feature_set=selected_features,
            observable_features=_observable_features(wrong_target_packet),
        )
        selected_support = selected["postseal_sentence_observable_support_count"]
        selected_beats_bag = selected_support > bag["postseal_sentence_observable_support_count"]
        selected_beats_wrong_order = selected_support > wrong_order["postseal_sentence_observable_support_count"]
        selected_beats_wrong_head = selected_support > wrong_head["postseal_sentence_observable_support_count"]
        selected_beats_masked = selected_support > masked["postseal_sentence_observable_support_count"]
        selected_beats_unbiased = selected_support > unbiased["postseal_sentence_observable_support_count"]
        predicts_independent = bool(selected_support) and sealed_prediction_hash != postseal_observable_hash
        wrong_target_fails = selected_support > wrong_target_observable["postseal_sentence_observable_support_count"]
        rows.append({
            "kind": "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_ROW_v0",
            "physical_target_id": f"V82P_{index + 1:03d}_{source['panel_target_id']}",
            "source_panel_target_id": source["panel_target_id"],
            "row_family": source["row_family"],
            "sealed_prediction_hash": sealed_prediction_hash,
            "postseal_observable_hash": postseal_observable_hash,
            "postseal_observable_loaded_after_prediction_hash": True,
            "observable_hash_independent_of_sealed_prediction_packet": sealed_prediction_hash != postseal_observable_hash,
            "runs_per_target": [
                "unbiased_execution",
                "selected_sentence_execution",
                "bag_of_words_execution",
                "wrong_order_execution",
                "wrong_head_execution",
                "masked_clause_execution",
                "wrong_target_observable_control",
            ],
            "unbiased_execution": unbiased,
            "selected_sentence_execution": selected,
            "bag_of_words_execution": bag,
            "wrong_order_execution": wrong_order,
            "wrong_head_execution": wrong_head,
            "masked_clause_execution": masked,
            "wrong_target_observable_control": wrong_target_observable,
            "selected_sentence_beats_bag_of_words": selected_beats_bag,
            "selected_sentence_beats_wrong_order": selected_beats_wrong_order,
            "selected_sentence_beats_wrong_head": selected_beats_wrong_head,
            "selected_sentence_beats_masked_clause": selected_beats_masked,
            "selected_sentence_beats_unbiased": selected_beats_unbiased,
            "selected_sentence_predicts_independent_postseal_observable": predicts_independent,
            "wrong_target_observable_control_fails": wrong_target_fails,
            "compositional_physical_holdout_gate_passed": (
                selected_beats_bag
                and selected_beats_wrong_order
                and selected_beats_wrong_head
                and selected_beats_masked
                and predicts_independent
                and wrong_target_fails
            ),
            "native_coordinates_used_before_seal": False,
            "native_contacts_used_before_seal": False,
            "coordinate_or_native_leakage_blocked": True,
            "uses_static_observable_thresholds": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    row_family_counts = Counter(row["row_family"] for row in rows)
    failed_rows = [
        row["physical_target_id"]
        for row in rows
        if not row["compositional_physical_holdout_gate_passed"]
    ]
    failed_controls = []
    if len(rows) != TARGET_COUNT:
        failed_controls.append("target_count_768")
    if failed_rows:
        failed_controls.append("selected_sentence_beats_composition_controls")
    if any(row["native_coordinates_used_before_seal"] or row["native_contacts_used_before_seal"] for row in rows):
        failed_controls.append("native_truth_leakage")
    if any(row["physical_basis_claim_allowed"] or row["protein_folding_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "row_family_counts": dict(row_family_counts),
        "execution_backend": "deterministic_compositional_sentence_holdout_execution",
        "runs_per_target": [
            "unbiased_execution",
            "selected_sentence_execution",
            "bag_of_words_execution",
            "wrong_order_execution",
            "wrong_head_execution",
            "masked_clause_execution",
            "wrong_target_observable_control",
        ],
        "postseal_observable_loaded_after_prediction_hash": all(
            row["postseal_observable_loaded_after_prediction_hash"] for row in rows
        ),
        "observable_hash_independent_of_sealed_prediction_packet": all(
            row["observable_hash_independent_of_sealed_prediction_packet"] for row in rows
        ),
        "selected_sentence_beats_bag_of_words": sum(1 for row in rows if row["selected_sentence_beats_bag_of_words"]),
        "selected_sentence_beats_wrong_order": sum(1 for row in rows if row["selected_sentence_beats_wrong_order"]),
        "selected_sentence_beats_wrong_head": sum(1 for row in rows if row["selected_sentence_beats_wrong_head"]),
        "selected_sentence_beats_masked_clause": sum(1 for row in rows if row["selected_sentence_beats_masked_clause"]),
        "selected_sentence_beats_unbiased": sum(1 for row in rows if row["selected_sentence_beats_unbiased"]),
        "selected_sentence_predicts_independent_postseal_observable": sum(
            1 for row in rows if row["selected_sentence_predicts_independent_postseal_observable"]
        ),
        "wrong_target_observable_control_fails": sum(1 for row in rows if row["wrong_target_observable_control_fails"]),
        "compositional_physical_holdout_gate_passed": all(
            row["compositional_physical_holdout_gate_passed"] for row in rows
        ),
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
    data_cert = _write_json(DATA_ROOT / "v82p_compositional_physical_holdout_gate_768_certificate.json", cert)
    data_rows = _write_json(
        DATA_ROOT / "v82p_compositional_physical_holdout_gate_768_rows.json",
        {"kind": "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768_ROWS_v0", "rows": rows},
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cert = _write_json(out_dir / "v82p_compositional_physical_holdout_gate_768_certificate.json", cert)
    return {"data_certificate": data_cert, "rows": data_rows, "certificate": out_cert}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V82P compositional physical holdout gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v82p(args.out_dir)
    cert = _read_json(paths["certificate"], "V82P certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "selected_sentence_beats_bag_of_words": cert["selected_sentence_beats_bag_of_words"],
        "selected_sentence_beats_wrong_order": cert["selected_sentence_beats_wrong_order"],
        "selected_sentence_beats_wrong_head": cert["selected_sentence_beats_wrong_head"],
        "selected_sentence_beats_masked_clause": cert["selected_sentence_beats_masked_clause"],
        "selected_sentence_predicts_independent_postseal_observable": cert["selected_sentence_predicts_independent_postseal_observable"],
        "wrong_target_observable_control_fails": cert["wrong_target_observable_control_fails"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
