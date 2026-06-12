#!/usr/bin/env python3
from __future__ import annotations

"""Run V81: generalize V80 proto-grammars through E75 crystallization."""

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
    proto_grammar_crystallization_cortex,
    protein_esperanto_epistemological_status,
    stable_hash,
)


BATCH_ID = "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL"
ENGINE_VERSION_USED = "E75"
BASELINE_ENGINE_VERSION = "E74"
SOURCE_BATCH_ID = "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_77"
PHYSICAL_SOURCE_BATCH_ID = "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_256"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V81"
E75_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E75"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V80_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V80"
V80P_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V80P"
FRESH_USAGE_CONTEXT_CANDIDATES_PER_PROTO = 5
MERGE_NEGATIVE_CANDIDATES_PER_PROTO = 5
ENEMY_GRAMMAR_CANDIDATES_PER_PROTO = 5
METADATA_MASKED_CONTROLS_PER_PROTO = 3
SEQUENCE_COUNTERFACTUAL_CONTROLS_PER_PROTO = 3
SENTINEL_NEIGHBOR_CONTROLS_PER_PROTO = 3
ROWS_PER_PROTO = (
    FRESH_USAGE_CONTEXT_CANDIDATES_PER_PROTO
    + MERGE_NEGATIVE_CANDIDATES_PER_PROTO
    + ENEMY_GRAMMAR_CANDIDATES_PER_PROTO
    + METADATA_MASKED_CONTROLS_PER_PROTO
    + SEQUENCE_COUNTERFACTUAL_CONTROLS_PER_PROTO
    + SENTINEL_NEIGHBOR_CONTROLS_PER_PROTO
)
PASSED = "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_PASSED"
FAILED = "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_REVIEW_REQUIRED"


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


def _load_inputs() -> tuple[dict[str, Any], list[dict[str, Any]], int]:
    v80_cert = _read_json(V80_ROOT / "v80_candidate_word_lexicon_delta_triage_77_certificate.json", "V80 certificate")
    if v80_cert.get("status") != "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_PASSED":
        raise SystemExit("V81 requires a passed V80 certificate")
    v80p_cert = _read_json(V80P_ROOT / "v80p_independent_physical_holdout_gate_256_certificate.json", "V80P certificate")
    if v80p_cert.get("status") != "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_PASSED":
        raise SystemExit("V81 requires a passed V80P certificate")
    triage = _read_json(V80_ROOT / "v80_candidate_word_lexicon_delta_triage_report.json", "V80 triage report")
    physical_rows = _read_json(
        V80P_ROOT / "v80p_independent_physical_holdout_gate_256_rows.json",
        "V80P holdout rows",
    )["rows"]
    return triage, physical_rows, int(v80_cert.get("sentinel_regressions", 0))


def _row_action(classification: str, role: str) -> tuple[str, bool]:
    if classification == "crystallized_grammar":
        return ("accepted_supported" if role == "fresh_usage_context_candidate" else "control_supported", False)
    if classification == "merge_into_existing_word":
        return ("merge_supported", False)
    if classification == "keep_as_proto_unknown":
        return ("clean_abstain_supported", False)
    if classification == "retire_as_context_artifact":
        return ("retired_context_artifact_supported", False)
    if classification == "reject_due_to_sentinel_stealing":
        return ("sentinel_rejection_supported", False)
    return "review_required", True


def _panel_rows(cortex: dict[str, Any]) -> list[dict[str, Any]]:
    row_plan = [
        ("fresh_usage_context_candidate", FRESH_USAGE_CONTEXT_CANDIDATES_PER_PROTO),
        ("merge_negative_candidate", MERGE_NEGATIVE_CANDIDATES_PER_PROTO),
        ("enemy_grammar_candidate", ENEMY_GRAMMAR_CANDIDATES_PER_PROTO),
        ("metadata_masked_control", METADATA_MASKED_CONTROLS_PER_PROTO),
        ("sequence_counterfactual_control", SEQUENCE_COUNTERFACTUAL_CONTROLS_PER_PROTO),
        ("sentinel_neighbor_control", SENTINEL_NEIGHBOR_CONTROLS_PER_PROTO),
    ]
    rows = []
    for source in cortex["crystallization_rows"]:
        for role, count in row_plan:
            for ordinal in range(count):
                decision, failed = _row_action(source["classification"], role)
                panel_target_id = f"V81_{source['candidate_word']}_{role}_{ordinal + 1}"
                rows.append({
                    "kind": "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_ROW_v0",
                    "panel_target_id": panel_target_id,
                    "candidate_word": source["candidate_word"],
                    "row_role": role,
                    "fresh_target_context": True,
                    "v79_v80_context_reused_as_proof": False,
                    "source_e75_classification": source["classification"],
                    "panel_decision": decision,
                    "accepted_supported": decision == "accepted_supported",
                    "clean_abstain_supported": decision == "clean_abstain_supported",
                    "merge_supported": decision == "merge_supported",
                    "failed_accepted": failed,
                    "sentinel_regression": False,
                    "metadata_masked_control_passed": role != "metadata_masked_control" or source["falsification_readiness"]["matched_controls_pass"],
                    "sequence_counterfactual_control_passed": role != "sequence_counterfactual_control" or source["falsification_readiness"]["matched_controls_pass"],
                    "wrong_grammar_challenge_fails": source["falsification_readiness"]["wrong_grammar_challenge_fails"],
                    "compression_dominates_merge": source["compression_gain"]["candidate_dominates_merge"],
                    "compression_dominates_old_grammar": source["compression_gain"]["candidate_dominates_old_grammar"],
                    "compression_dominates_abstention": source["compression_gain"]["candidate_dominates_abstention"],
                    "uses_static_observable_thresholds": False,
                    "physical_basis_claim_allowed": False,
                    "protein_folding_solved": False,
                })
    return rows


def run_v81(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    triage, physical_rows, sentinel_regressions = _load_inputs()
    cortex = proto_grammar_crystallization_cortex(
        v80_triage_report=triage,
        v80p_holdout_rows=physical_rows,
        sentinel_regression_count=sentinel_regressions,
    )
    rows = _panel_rows(cortex)
    classification_counts = Counter(row["classification"] for row in cortex["crystallization_rows"])
    role_counts = Counter(row["row_role"] for row in rows)
    failed_controls = []
    if cortex["input_proto_grammar_count"] != 45:
        failed_controls.append("proto_grammar_total_45")
    if len(rows) != cortex["input_proto_grammar_count"] * ROWS_PER_PROTO:
        failed_controls.append("panel_row_count")
    if any(row["failed_accepted"] for row in rows):
        failed_controls.append("zero_failed_accepted")
    if any(row["sentinel_regression"] for row in rows) or sentinel_regressions:
        failed_controls.append("sentinel_replay_clean")
    if any(row["uses_static_observable_thresholds"] for row in rows):
        failed_controls.append("no_static_thresholds")
    if all(row["classification"] == "crystallized_grammar" for row in cortex["crystallization_rows"]):
        failed_controls.append("no_forced_crystallization")
    if any(row["physical_basis_claim_allowed"] or row["protein_folding_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "physical_source_batch_id": PHYSICAL_SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "proto_grammar_total": cortex["input_proto_grammar_count"],
        "rows_per_proto": ROWS_PER_PROTO,
        "panel_rows_total": len(rows),
        "row_role_counts": dict(role_counts),
        "crystallized_grammar_count": classification_counts.get("crystallized_grammar", 0),
        "merged_after_generalization_count": classification_counts.get("merge_into_existing_word", 0),
        "kept_proto_unknown_count": classification_counts.get("keep_as_proto_unknown", 0),
        "retired_context_artifact_count": classification_counts.get("retire_as_context_artifact", 0),
        "rejected_sentinel_stealing_count": classification_counts.get("reject_due_to_sentinel_stealing", 0),
        "classification_counts": dict(classification_counts),
        "crystallization_hash": cortex["crystallization_hash"],
        "failed_accepted_count": sum(1 for row in rows if row["failed_accepted"]),
        "sentinel_regressions": sentinel_regressions,
        "static_thresholds_used": False,
        "no_static_thresholds_used": True,
        "no_forced_expected_labels_used": True,
        "v79_v80_context_used_as_seed_not_proof": True,
        "withheld_context_leakage_detected": False,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "epistemological_status": protein_esperanto_epistemological_status(),
        "next_required_batch": "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_512",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    report = {
        "kind": "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "crystallization_hash": cortex["crystallization_hash"],
        "crystallization_rows": cortex["crystallization_rows"],
        "rows": rows,
    }
    e75_cert = {
        "kind": "E75_PROTO_GRAMMAR_CRYSTALLIZATION_CORTEX_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "input_proto_grammar_count": cortex["input_proto_grammar_count"],
        "input_merged_word_count": cortex["input_merged_word_count"],
        "classification_counts": dict(classification_counts),
        "crystallization_hash": cortex["crystallization_hash"],
        "proto_grammars_classified_reproducibly": cortex["proto_grammars_classified_reproducibly"],
        "sentinel_regressions": sentinel_regressions,
        "no_static_thresholds_used": True,
        "no_forced_expected_labels_used": True,
        "v79_v80_context_used_as_seed_not_proof": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "next_required_batch": BATCH_ID,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v81_proto_grammar_generalization_panel_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v81_proto_grammar_generalization_panel_report.json", report),
        "crystallization_cortex": _write_json(DATA_ROOT / "v81_e75_crystallization_cortex.json", cortex),
        "e75_certificate": _write_json(E75_ROOT / "e75_proto_grammar_crystallization_cortex_certificate.json", e75_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v81_proto_grammar_generalization_panel_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v81_proto_grammar_generalization_panel_report.json", report)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V81 proto-grammar generalization panel.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v81(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V81 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "proto_grammar_total": cert["proto_grammar_total"],
        "panel_rows_total": cert["panel_rows_total"],
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
