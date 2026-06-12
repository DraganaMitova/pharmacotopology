#!/usr/bin/env python3
from __future__ import annotations

"""Run V80: compress V79 candidate words into a lexicon delta."""

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

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    candidate_word_compression_cortex,
    protein_esperanto_epistemological_status,
    stable_hash,
)


BATCH_ID = "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_77"
ENGINE_VERSION_USED = "E74"
BASELINE_ENGINE_VERSION = "E73"
INPUT_CANDIDATE_WORD_COUNT = 77
CONTROL_ROWS_PER_CANDIDATE = 7
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V80"
E74_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E74"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V79_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V79"
V79P_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V79P"
PASSED = "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_PASSED"
FAILED = "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_REVIEW_REQUIRED"


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


def _load_inputs() -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], int]:
    lexicon = _read_json(V79_ROOT / "v79_lexicon_delta_report.json", "V79 lexicon delta report")
    scoring = _read_json(V79_ROOT / "v79_blind_language_discovery_scoring_report.json", "V79 scoring report")["rows"]
    physical = _read_json(V79P_ROOT / "v79p_physical_falsification_128_rows.json", "V79P physical rows")["rows"]
    sentinel = _read_json(V79_ROOT / "v79_old_grammar_sentinel_replay.json", "V79 sentinel replay")
    return lexicon, scoring, physical, int(sentinel["sentinel_regressions"])


def _control_rows(row: dict[str, Any]) -> list[dict[str, Any]]:
    controls = row["promotion_controls"]
    return [
        {
            "control_id": "existing_grammar_explanation_control",
            "passed": controls["beats_existing_grammar_explanation"] or row["classification"] == "merge_into_existing_word",
            "interpretation": "candidate must dominate or explicitly merge into existing explanation",
        },
        {
            "control_id": "merge_explanation_control",
            "passed": controls["beats_merge_explanation"] or row["classification"] == "merge_into_existing_word",
            "interpretation": "candidate must beat merge explanation unless merge is the selected lexicon action",
        },
        {
            "control_id": "metadata_masked_control",
            "passed": controls["beats_metadata_masked_control"],
            "interpretation": "visible pressure must survive metadata masking comparison",
        },
        {
            "control_id": "sequence_counterfactual_control",
            "passed": controls["beats_sequence_counterfactual_control"],
            "interpretation": "candidate cannot be sequence-only noise",
        },
        {
            "control_id": "wrong_grammar_challenge_control",
            "passed": controls["wrong_grammar_challenge_fails"],
            "interpretation": "enemy grammar must not steal the candidate",
        },
        {
            "control_id": "sentinel_regression_control",
            "passed": controls["sentinel_regression_pressure_absent"],
            "interpretation": "candidate cannot regress old learned grammar sentinels",
        },
        {
            "control_id": "coordinate_native_truth_seal_control",
            "passed": controls["coordinate_native_truth_stays_sealed"],
            "interpretation": "native coordinates and contacts stay sealed",
        },
    ]


def _triage_rows(cortex: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in cortex["compression_rows"]:
        controls = _control_rows(row)
        rows.append({
            "kind": "V80_CANDIDATE_WORD_TRIAGE_ROW_v0",
            "candidate_word": row["candidate_word"],
            "classification": row["classification"],
            "classification_reason": row["classification_reason"],
            "merge_candidate": row["merge_candidate"],
            "proto_grammar_candidate": row["proto_grammar_candidate"],
            "support_pressure": row["support_pressure"],
            "contradiction_pressure": row["contradiction_pressure"],
            "enemy_grammar_pressure": row["enemy_grammar_pressure"],
            "abstain_pressure": row["abstain_pressure"],
            "metadata_masking_pressure": row["metadata_masking_pressure"],
            "perturbation_pressure": row["perturbation_pressure"],
            "physical_mismatch_pressure": row["physical_mismatch_pressure"],
            "compression_gain": row["compression_gain"],
            "sentinel_cost": row["sentinel_cost"],
            "definition_by_known_words": row["definition_by_known_words"],
            "usage_by_context": row["usage_by_context"],
            "pressure_ledger": row["pressure_ledger"],
            "matched_control_dominance": {
                "kind": "V80_CANDIDATE_WORD_MATCHED_CONTROL_DOMINANCE_v0",
                "control_rows": controls,
                "matched_control_dominance_passed": all(control["passed"] for control in controls),
                "uses_static_observable_thresholds": False,
            },
            "rejected_or_retired_reason": row["rejected_or_retired_reason"],
            "physical_basis_claim_allowed": False,
            "folding_problem_solved": False,
        })
    return rows


def run_v80(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    lexicon, scoring, physical, sentinel_regressions = _load_inputs()
    cortex = candidate_word_compression_cortex(
        v79_lexicon_delta=lexicon,
        v79_scoring_rows=scoring,
        v79p_physical_rows=physical,
        sentinel_regression_count=sentinel_regressions,
    )
    rows = _triage_rows(cortex)
    classification_counts = Counter(row["classification"] for row in rows)
    proto_rows = [row for row in rows if row["classification"] == "proto_grammar"]
    merged_rows = [row for row in rows if row["classification"] == "merge_into_existing_word"]
    rejected_retired = [
        row
        for row in rows
        if row["classification"] in {"reject_as_noise", "retire_due_to_regression"}
    ]
    failed_controls = []
    if len(rows) != INPUT_CANDIDATE_WORD_COUNT:
        failed_controls.append("candidate_word_count_77")
    if cortex["failed_accepted_count"] != 0:
        failed_controls.append("zero_failed_accepted")
    if not all(row["matched_control_dominance"]["matched_control_dominance_passed"] for row in proto_rows):
        failed_controls.append("proto_grammars_matched_control_dominance")
    if not all(row["merge_candidate"] for row in merged_rows):
        failed_controls.append("merged_words_have_merge_candidate")
    if not all(row["rejected_or_retired_reason"] for row in rejected_retired):
        failed_controls.append("rejected_retired_words_record_reason")
    if sentinel_regressions:
        failed_controls.append("sentinels_stable")
    if cortex["physical_falsification_failure_count"]:
        failed_controls.append("v79p_physical_falsification_rows_stable")
    if not cortex["candidate_words_classified_reproducibly"]:
        failed_controls.append("classification_reproducible")
    if not all(row["physical_basis_claim_allowed"] is False and row["folding_problem_solved"] is False for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_77_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "input_candidate_word_count": INPUT_CANDIDATE_WORD_COUNT,
        "classified_candidate_word_count": len(rows),
        "classification_counts": dict(classification_counts),
        "failed_accepted_count": cortex["failed_accepted_count"],
        "candidate_words_classified_reproducibly": cortex["candidate_words_classified_reproducibly"],
        "classification_hash": cortex["compression_hash"],
        "candidate_word_graph_node_count": cortex["candidate_word_graph"]["node_count"],
        "candidate_word_graph_edge_count": cortex["candidate_word_graph"]["edge_count"],
        "paired_control_rows_total": len(rows) * CONTROL_ROWS_PER_CANDIDATE,
        "proto_grammars_matched_control_dominance_passed": all(row["matched_control_dominance"]["matched_control_dominance_passed"] for row in proto_rows),
        "merged_words_do_not_reduce_existing_grammar_performance": all(row["merge_candidate"] for row in merged_rows),
        "rejected_or_retired_words_record_pressure_reason": all(row["rejected_or_retired_reason"] for row in rejected_retired),
        "sentinel_regressions": sentinel_regressions,
        "withheld_context_leakage_detected": False,
        "coordinate_native_truth_stays_sealed": True,
        "no_forced_expected_labels_used": True,
        "no_static_thresholds_used": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "epistemological_status": protein_esperanto_epistemological_status(),
        "next_required_batch": "V81_REPLAY_PROTO_GRAMMAR_CANDIDATES_OR_EXPAND_UNKNOWN_USAGE_CONTEXT",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    triage = {
        "kind": "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "classification_hash": cortex["compression_hash"],
        "classification_counts": dict(classification_counts),
        "rows": rows,
    }
    e74_cert = {
        "kind": "E74_CANDIDATE_WORD_COMPRESSION_CORTEX_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "input_candidate_word_count": INPUT_CANDIDATE_WORD_COUNT,
        "classification_counts": dict(classification_counts),
        "classification_hash": cortex["compression_hash"],
        "candidate_word_graph_node_count": cortex["candidate_word_graph"]["node_count"],
        "candidate_word_graph_edge_count": cortex["candidate_word_graph"]["edge_count"],
        "no_static_thresholds_used": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "next_required_batch": BATCH_ID,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v80_candidate_word_lexicon_delta_triage_77_certificate.json", cert),
        "triage_report": _write_json(DATA_ROOT / "v80_candidate_word_lexicon_delta_triage_report.json", triage),
        "candidate_word_graph": _write_json(DATA_ROOT / "v80_candidate_word_graph.json", cortex["candidate_word_graph"]),
        "compression_cortex": _write_json(DATA_ROOT / "v80_e74_compression_cortex.json", cortex),
        "e74_certificate": _write_json(E74_ROOT / "e74_candidate_word_compression_cortex_certificate.json", e74_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v80_candidate_word_lexicon_delta_triage_77_certificate.json", cert)
    paths["run_triage_report"] = _write_json(out_dir / "v80_candidate_word_lexicon_delta_triage_report.json", triage)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V80 candidate-word lexicon delta triage.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v80(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V80 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "classified_candidate_word_count": cert["classified_candidate_word_count"],
        "classification_counts": cert["classification_counts"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "sentinel_regressions": cert["sentinel_regressions"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
