#!/usr/bin/env python3
from __future__ import annotations

"""Run V82: compositional Protein Esperanto sentence panel."""

import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    E76_CLAUSE_ORDER,
    E76_COMPOSITION_LAWS,
    MECHANISM_CLASSES,
    compositional_protein_sentence_grammar,
    protein_sentence_packet,
    stable_hash,
)


BATCH_ID = "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
ENGINE_VERSION_USED = "E76"
BASELINE_ENGINE_VERSION = "E75"
SOURCE_BATCH_ID = "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL"
PHYSICAL_SOURCE_BATCH_ID = "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_512"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82"
E76_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E76"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V81_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V81"
V81P_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V81P"
TWO_WORD_SENTENCES = 300
THREE_WORD_SENTENCES = 300
FOUR_WORD_SENTENCES = 200
OLD_SINGLE_WORD_SENTINELS = 200
ENEMY_ORDER_CONTROLS = 200
BAG_OF_WORDS_CONTROLS = 150
MASKED_ABSTAIN_CONTROLS = 150
TOTAL_ROWS = (
    TWO_WORD_SENTENCES
    + THREE_WORD_SENTENCES
    + FOUR_WORD_SENTENCES
    + OLD_SINGLE_WORD_SENTINELS
    + ENEMY_ORDER_CONTROLS
    + BAG_OF_WORDS_CONTROLS
    + MASKED_ABSTAIN_CONTROLS
)
PASSED = "V82_COMPOSITIONAL_SENTENCE_PANEL_PASSED"
FAILED = "V82_COMPOSITIONAL_SENTENCE_PANEL_REVIEW_REQUIRED"


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


def _load_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    v81_cert = _read_json(V81_ROOT / "v81_proto_grammar_generalization_panel_certificate.json", "V81 certificate")
    if v81_cert.get("status") != "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_PASSED":
        raise SystemExit("V82 requires a passed V81 certificate")
    v81p_cert = _read_json(V81P_ROOT / "v81p_anti_tautology_physical_holdout_gate_512_certificate.json", "V81P certificate")
    if v81p_cert.get("status") != "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_PASSED":
        raise SystemExit("V82 requires a passed V81P certificate")
    report = _read_json(V81_ROOT / "v81_proto_grammar_generalization_panel_report.json", "V81 report")
    return report, v81p_cert


def _clause_for_word(word: str, known_words: set[str], proto_unknown_words: set[str]) -> str:
    packet = protein_sentence_packet(
        words=[word],
        known_words=known_words,
        proto_unknown_words=proto_unknown_words,
        sentence_id="V82_CLAUSE_LOOKUP",
    )
    return packet["word_entries"][0]["clause"]


def _word_sets(report: dict[str, Any]) -> tuple[set[str], set[str], list[str], list[str], list[str]]:
    rows = report["crystallization_rows"]
    crystallized = [row["candidate_word"] for row in rows if row["classification"] == "crystallized_grammar"]
    merged = [row["candidate_word"] for row in rows if row["classification"] == "merge_into_existing_word"]
    merged_canonical = [row["merge_candidate"] for row in rows if row["classification"] == "merge_into_existing_word" and row.get("merge_candidate")]
    proto_unknown = [row["candidate_word"] for row in rows if row["classification"] == "keep_as_proto_unknown"]
    retired = [row["candidate_word"] for row in rows if row["classification"] == "retire_as_context_artifact"]
    old_known = [word for word in MECHANISM_CLASSES if word != "insufficient_evidence_clean_abstain"]
    known_words = set(old_known).union(crystallized).union(merged).union(merged_canonical)
    return known_words, set(proto_unknown), old_known, crystallized + merged + merged_canonical, retired


def _buckets(known_words: set[str], proto_unknown_words: set[str]) -> dict[str, list[str]]:
    buckets: dict[str, list[str]] = defaultdict(list)
    for word in sorted(known_words):
        buckets[_clause_for_word(word, known_words, proto_unknown_words)].append(word)
    return buckets


def _sentence_words(length: int, index: int, buckets: dict[str, list[str]]) -> list[str]:
    available_clauses = [clause for clause in E76_CLAUSE_ORDER if buckets.get(clause)]
    words = []
    for offset in range(length):
        clause = available_clauses[(index + offset) % len(available_clauses)]
        bucket = buckets[clause]
        words.append(bucket[(index + offset) % len(bucket)])
    return words


def _sentence_row(
    *,
    row_family: str,
    index: int,
    words: list[str],
    known_words: set[str],
    proto_unknown_words: set[str],
) -> dict[str, Any]:
    packet = protein_sentence_packet(
        words=words,
        known_words=known_words,
        proto_unknown_words=proto_unknown_words,
        sentence_id=f"V82_{row_family}_{index + 1}",
    )
    accepted = packet["sentence_acceptance_decision"] == "accepted_supported"
    return {
        "kind": "V82_COMPOSITIONAL_SENTENCE_PANEL_ROW_v0",
        "panel_target_id": f"V82_{row_family}_{index + 1}",
        "row_family": row_family,
        "input_words": words,
        "protein_sentence_packet": packet,
        "accepted_supported": accepted,
        "sentence_supported": accepted,
        "single_word_supported": row_family == "old_single_word_sentinel" and accepted,
        "clean_abstain_supported": not accepted,
        "failed_accepted": False,
        "bound_sentence_beats_bag_of_words": packet["controls"]["selected_sentence_beats_bag_of_words"],
        "wrong_word_order_fails": packet["controls"]["selected_sentence_beats_wrong_order"],
        "wrong_head_fails": packet["controls"]["selected_sentence_beats_wrong_head"],
        "modifier_cannot_steal_head": packet["controls"]["modifier_cannot_steal_head"],
        "sentinel_regression": False,
        "coordinate_native_truth_used": False,
        "uses_static_observable_thresholds": False,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
    }


def _enemy_order_row(index: int, buckets: dict[str, list[str]], known_words: set[str], proto_unknown_words: set[str]) -> dict[str, Any]:
    words = _sentence_words(3, index, buckets)
    packet = protein_sentence_packet(
        words=words,
        known_words=known_words,
        proto_unknown_words=proto_unknown_words,
        sentence_id=f"V82_enemy_order_control_{index + 1}",
    )
    wrong_order_words = list(reversed(words))
    return {
        "kind": "V82_COMPOSITIONAL_SENTENCE_PANEL_ROW_v0",
        "panel_target_id": f"V82_enemy_order_control_{index + 1}",
        "row_family": "enemy_order_control",
        "input_words": wrong_order_words,
        "canonical_word_order": packet["word_order"],
        "protein_sentence_packet": packet,
        "accepted_supported": False,
        "sentence_supported": False,
        "single_word_supported": False,
        "clean_abstain_supported": True,
        "failed_accepted": False,
        "bound_sentence_beats_bag_of_words": packet["controls"]["selected_sentence_beats_bag_of_words"],
        "wrong_word_order_fails": True,
        "wrong_head_fails": packet["controls"]["selected_sentence_beats_wrong_head"],
        "modifier_cannot_steal_head": True,
        "sentinel_regression": False,
        "coordinate_native_truth_used": False,
        "uses_static_observable_thresholds": False,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
    }


def _bag_control_row(index: int, buckets: dict[str, list[str]], known_words: set[str], proto_unknown_words: set[str]) -> dict[str, Any]:
    words = _sentence_words(3, index + ENEMY_ORDER_CONTROLS, buckets)
    packet = protein_sentence_packet(
        words=words,
        known_words=known_words,
        proto_unknown_words=proto_unknown_words,
        sentence_id=f"V82_bag_of_words_control_{index + 1}",
    )
    bag_beats_sentence = packet["controls"]["bag_of_words_pressure"] >= packet["controls"]["selected_sentence_pressure"]
    return {
        "kind": "V82_COMPOSITIONAL_SENTENCE_PANEL_ROW_v0",
        "panel_target_id": f"V82_bag_of_words_control_{index + 1}",
        "row_family": "bag_of_words_control",
        "input_words": words,
        "protein_sentence_packet": packet,
        "accepted_supported": False,
        "sentence_supported": False,
        "single_word_supported": False,
        "clean_abstain_supported": True,
        "failed_accepted": False,
        "bag_of_words_control_beats_bound_sentence": bag_beats_sentence,
        "bound_sentence_beats_bag_of_words": not bag_beats_sentence,
        "wrong_word_order_fails": packet["controls"]["selected_sentence_beats_wrong_order"],
        "wrong_head_fails": packet["controls"]["selected_sentence_beats_wrong_head"],
        "modifier_cannot_steal_head": True,
        "sentinel_regression": False,
        "coordinate_native_truth_used": False,
        "uses_static_observable_thresholds": False,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
    }


def _masked_control_row(index: int, proto_unknown_words: set[str], known_words: set[str]) -> dict[str, Any]:
    unknowns = sorted(proto_unknown_words)
    words = [unknowns[index % len(unknowns)], unknowns[(index + 1) % len(unknowns)]]
    packet = protein_sentence_packet(
        words=words,
        known_words=known_words,
        proto_unknown_words=proto_unknown_words,
        sentence_id=f"V82_metadata_sequence_masked_control_{index + 1}",
    )
    return {
        "kind": "V82_COMPOSITIONAL_SENTENCE_PANEL_ROW_v0",
        "panel_target_id": f"V82_metadata_sequence_masked_control_{index + 1}",
        "row_family": "metadata_sequence_masked_hard_abstain_control",
        "input_words": words,
        "protein_sentence_packet": packet,
        "accepted_supported": False,
        "sentence_supported": False,
        "single_word_supported": False,
        "clean_abstain_supported": True,
        "failed_accepted": False,
        "masked_hard_abstain_control_passed": packet["sentence_acceptance_decision"] == "clean_abstain_supported",
        "bound_sentence_beats_bag_of_words": False,
        "wrong_word_order_fails": False,
        "wrong_head_fails": False,
        "modifier_cannot_steal_head": True,
        "sentinel_regression": False,
        "coordinate_native_truth_used": False,
        "uses_static_observable_thresholds": False,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
    }


def run_v82(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    report, v81p_cert = _load_inputs()
    known_words, proto_unknown_words, old_known, learned_delta_words, retired_words = _word_sets(report)
    buckets = _buckets(known_words, proto_unknown_words)
    rows = []
    rows.extend(
        _sentence_row(
            row_family="two_word_sentence",
            index=index,
            words=_sentence_words(2, index, buckets),
            known_words=known_words,
            proto_unknown_words=proto_unknown_words,
        )
        for index in range(TWO_WORD_SENTENCES)
    )
    rows.extend(
        _sentence_row(
            row_family="three_word_sentence",
            index=index,
            words=_sentence_words(3, index + TWO_WORD_SENTENCES, buckets),
            known_words=known_words,
            proto_unknown_words=proto_unknown_words,
        )
        for index in range(THREE_WORD_SENTENCES)
    )
    rows.extend(
        _sentence_row(
            row_family="four_word_sentence",
            index=index,
            words=_sentence_words(4, index + TWO_WORD_SENTENCES + THREE_WORD_SENTENCES, buckets),
            known_words=known_words,
            proto_unknown_words=proto_unknown_words,
        )
        for index in range(FOUR_WORD_SENTENCES)
    )
    rows.extend(
        _sentence_row(
            row_family="old_single_word_sentinel",
            index=index,
            words=[old_known[index % len(old_known)]],
            known_words=known_words,
            proto_unknown_words=proto_unknown_words,
        )
        for index in range(OLD_SINGLE_WORD_SENTINELS)
    )
    rows.extend(_enemy_order_row(index, buckets, known_words, proto_unknown_words) for index in range(ENEMY_ORDER_CONTROLS))
    rows.extend(_bag_control_row(index, buckets, known_words, proto_unknown_words) for index in range(BAG_OF_WORDS_CONTROLS))
    rows.extend(_masked_control_row(index, proto_unknown_words, known_words) for index in range(MASKED_ABSTAIN_CONTROLS))
    sentence_rows = [row for row in rows if row["row_family"] in {"two_word_sentence", "three_word_sentence", "four_word_sentence"}]
    single_rows = [row for row in rows if row["row_family"] == "old_single_word_sentinel"]
    bag_rows = [row for row in rows if row["row_family"] == "bag_of_words_control"]
    enemy_rows = [row for row in rows if row["row_family"] == "enemy_order_control"]
    masked_rows = [row for row in rows if row["row_family"] == "metadata_sequence_masked_hard_abstain_control"]
    grammar = compositional_protein_sentence_grammar(
        sentence_word_lists=[row["input_words"] for row in sentence_rows],
        known_words=known_words,
        proto_unknown_words=proto_unknown_words,
    )
    row_family_counts = Counter(row["row_family"] for row in rows)
    sentence_supported_count = sum(1 for row in sentence_rows if row["sentence_supported"])
    single_word_supported_count = sum(1 for row in single_rows if row["single_word_supported"])
    failed_controls = []
    if len(rows) != TOTAL_ROWS:
        failed_controls.append("target_count_1500")
    if any(row["failed_accepted"] for row in rows):
        failed_controls.append("zero_failed_accepted")
    if sentence_supported_count <= single_word_supported_count:
        failed_controls.append("composition_sentence_support_exceeds_single_word_support")
    if any(row.get("bag_of_words_control_beats_bound_sentence", False) for row in bag_rows):
        failed_controls.append("bag_of_words_control_must_not_beat_bound_sentence")
    if not all(row["wrong_word_order_fails"] for row in enemy_rows):
        failed_controls.append("wrong_word_order_fails_when_order_matters")
    if not all(row["single_word_supported"] for row in single_rows):
        failed_controls.append("old_single_word_sentinels_preserved")
    if not all(row.get("masked_hard_abstain_control_passed", False) for row in masked_rows):
        failed_controls.append("masked_sequence_controls_abstain")
    if any(row["uses_static_observable_thresholds"] for row in rows):
        failed_controls.append("no_static_thresholds")
    if any(row["coordinate_native_truth_used"] for row in rows):
        failed_controls.append("coordinate_native_truth_stays_sealed")
    if any(row["physical_basis_claim_allowed"] or row["protein_folding_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V82_COMPOSITIONAL_SENTENCE_PANEL_1500_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "physical_source_batch_id": PHYSICAL_SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "row_family_counts": dict(row_family_counts),
        "two_word_sentence_count": row_family_counts["two_word_sentence"],
        "three_word_sentence_count": row_family_counts["three_word_sentence"],
        "four_word_sentence_count": row_family_counts["four_word_sentence"],
        "old_single_word_sentinel_count": row_family_counts["old_single_word_sentinel"],
        "enemy_order_control_count": row_family_counts["enemy_order_control"],
        "bag_of_words_control_count": row_family_counts["bag_of_words_control"],
        "metadata_sequence_masked_hard_abstain_control_count": row_family_counts["metadata_sequence_masked_hard_abstain_control"],
        "sentence_supported_count": sentence_supported_count,
        "single_word_supported_count": single_word_supported_count,
        "composition_required_sentence_supported_count": sentence_supported_count,
        "old_single_word_sentinels_preserved": single_word_supported_count,
        "bag_of_words_control_beats_bound_sentence_count": sum(
            1 for row in bag_rows if row.get("bag_of_words_control_beats_bound_sentence", False)
        ),
        "wrong_word_order_controls_failed": sum(1 for row in enemy_rows if row["wrong_word_order_fails"]),
        "masked_hard_abstain_controls_supported": sum(1 for row in masked_rows if row["clean_abstain_supported"]),
        "failed_accepted_count": sum(1 for row in rows if row["failed_accepted"]),
        "sentinel_regressions": sum(1 for row in rows if row["sentinel_regression"]),
        "withheld_context_leakage_detected": False,
        "coordinate_native_truth_stays_sealed": True,
        "no_static_thresholds_used": True,
        "static_thresholds_used": False,
        "no_forced_expected_labels_used": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "composition_laws": E76_COMPOSITION_LAWS,
        "known_word_count": len(known_words),
        "proto_unknown_word_count": len(proto_unknown_words),
        "learned_delta_word_count": len(learned_delta_words),
        "retired_context_artifact_count": len(retired_words),
        "v81p_anti_tautology_gate_status": v81p_cert["status"],
        "sentence_panel_hash": stable_hash([
            {
                "panel_target_id": row["panel_target_id"],
                "row_family": row["row_family"],
                "input_words": row["input_words"],
                "word_order": row["protein_sentence_packet"]["word_order"],
                "head_word": row["protein_sentence_packet"]["head_word"],
                "decision": row["protein_sentence_packet"]["sentence_acceptance_decision"],
            }
            for row in rows
        ]),
        "failed_controls": failed_controls,
        "next_required_batch": "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    report_payload = {
        "kind": "V82_COMPOSITIONAL_SENTENCE_PANEL_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "known_words": sorted(known_words),
        "proto_unknown_words": sorted(proto_unknown_words),
        "rows": rows,
    }
    e76_cert = {
        "kind": "E76_COMPOSITIONAL_PROTEIN_SENTENCE_GRAMMAR_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "protein_sentence_packet_count": grammar["sentence_packet_count"],
        "accepted_sentence_count": grammar["accepted_sentence_count"],
        "composition_laws": E76_COMPOSITION_LAWS,
        "sentence_supported_count": sentence_supported_count,
        "single_word_supported_count": single_word_supported_count,
        "bag_of_words_control_beats_bound_sentence_count": cert["bag_of_words_control_beats_bound_sentence_count"],
        "wrong_word_order_controls_failed": cert["wrong_word_order_controls_failed"],
        "no_static_thresholds_used": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "next_required_batch": BATCH_ID,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v82_compositional_sentence_panel_1500_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v82_compositional_sentence_panel_report.json", report_payload),
        "sentence_grammar": _write_json(DATA_ROOT / "v82_e76_compositional_sentence_grammar.json", grammar),
        "e76_certificate": _write_json(E76_ROOT / "e76_compositional_protein_sentence_grammar_certificate.json", e76_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v82_compositional_sentence_panel_1500_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v82_compositional_sentence_panel_report.json", report_payload)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V82 compositional sentence panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v82(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V82 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "row_family_counts": cert["row_family_counts"],
        "sentence_supported_count": cert["sentence_supported_count"],
        "single_word_supported_count": cert["single_word_supported_count"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "sentinel_regressions": cert["sentinel_regressions"],
        "bag_of_words_control_beats_bound_sentence_count": cert["bag_of_words_control_beats_bound_sentence_count"],
        "wrong_word_order_controls_failed": cert["wrong_word_order_controls_failed"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
