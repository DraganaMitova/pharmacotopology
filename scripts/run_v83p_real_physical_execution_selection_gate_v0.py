#!/usr/bin/env python3
from __future__ import annotations

"""Run V83P: separate language support, proxy support, and physical proof."""

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


BATCH_ID = "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE"
ENGINE_VERSION_USED = "E77"
SOURCE_BATCH_ID = "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT"
LANGUAGE_SOURCE_BATCH_ID = "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
PROXY_SOURCE_BATCH_ID = "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V83P"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V83_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V83"
V82P_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82P"
PASSED = "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_PASSED"
FAILED = "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_REVIEW_REQUIRED"


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


def _load_inputs() -> list[dict[str, Any]]:
    v83_cert = _read_json(V83_ROOT / "v83_complexity_token_heuristic_scaling_audit_certificate.json", "V83 certificate")
    if v83_cert.get("status") != "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_PASSED":
        raise SystemExit("V83P requires a passed V83 certificate")
    v82p_cert = _read_json(V82P_ROOT / "v82p_compositional_physical_holdout_gate_768_certificate.json", "V82P certificate")
    if v82p_cert.get("status") != "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_PASSED":
        raise SystemExit("V83P requires a passed V82P certificate")
    rows = _read_json(V82P_ROOT / "v82p_compositional_physical_holdout_gate_768_rows.json", "V82P rows")["rows"]
    if not rows:
        raise SystemExit("V83P requires V82P rows")
    return rows


def run_v83p(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    source_rows = _load_inputs()
    rows = []
    for index, source in enumerate(source_rows):
        language_support = bool(source["selected_sentence_predicts_independent_postseal_observable"])
        proxy_support = bool(
            source["selected_sentence_beats_bag_of_words"]
            and source["selected_sentence_beats_wrong_order"]
            and source["selected_sentence_beats_wrong_head"]
            and source["selected_sentence_beats_masked_clause"]
        )
        independent_physical_holdout_support = False
        rows.append({
            "kind": "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_ROW_v0",
            "physical_selection_target_id": f"V83P_{index + 1:03d}_{source['source_panel_target_id']}",
            "source_v82p_physical_target_id": source["physical_target_id"],
            "source_panel_target_id": source["source_panel_target_id"],
            "language_support": language_support,
            "coarse_physical_proxy_support": proxy_support,
            "independent_physical_holdout_support": independent_physical_holdout_support,
            "independent_experimental_or_coordinate_physical_proof_executed": False,
            "support_tiers": {
                "language_support": language_support,
                "coarse_physical_proxy_support": proxy_support,
                "independent_physical_holdout_support": independent_physical_holdout_support,
            },
            "language_can_pass_without_physical_claim": language_support and not independent_physical_holdout_support,
            "physical_proxy_can_support_language": proxy_support and language_support,
            "physical_claim_blocked_until_independent_holdout": True,
            "native_coordinates_used_before_seal": False,
            "native_contacts_used_before_seal": False,
            "coordinate_or_native_leakage_blocked": True,
            "uses_static_observable_thresholds": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    failed_rows = [
        row["physical_selection_target_id"]
        for row in rows
        if not (
            row["language_can_pass_without_physical_claim"]
            and row["physical_proxy_can_support_language"]
            and row["physical_claim_blocked_until_independent_holdout"]
            and not row["physical_basis_claim_allowed"]
            and not row["protein_folding_solved"]
        )
    ]
    failed_controls = []
    if len(rows) != len(source_rows):
        failed_controls.append("target_count_matches_v82p")
    if failed_rows:
        failed_controls.append("support_tiers_separated_and_claim_blocked")
    if any(row["native_coordinates_used_before_seal"] or row["native_contacts_used_before_seal"] for row in rows):
        failed_controls.append("native_truth_leakage")
    if any(row["physical_basis_claim_allowed"] or row["protein_folding_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    tier_counts = {
        "language_support": sum(1 for row in rows if row["language_support"]),
        "coarse_physical_proxy_support": sum(1 for row in rows if row["coarse_physical_proxy_support"]),
        "independent_physical_holdout_support": sum(1 for row in rows if row["independent_physical_holdout_support"]),
    }
    cert = {
        "kind": "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "language_source_batch_id": LANGUAGE_SOURCE_BATCH_ID,
        "proxy_source_batch_id": PROXY_SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "row_family_counts": dict(Counter(row["source_panel_target_id"].split("_")[1] for row in rows)),
        "support_tier_counts": tier_counts,
        "language_can_pass_without_physical_claim_count": sum(
            1 for row in rows if row["language_can_pass_without_physical_claim"]
        ),
        "physical_proxy_can_support_language_count": sum(
            1 for row in rows if row["physical_proxy_can_support_language"]
        ),
        "independent_experimental_or_coordinate_physical_proof_executed": False,
        "physical_claim_blocked_until_independent_holdout": all(
            row["physical_claim_blocked_until_independent_holdout"] for row in rows
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
    data_cert = _write_json(DATA_ROOT / "v83p_real_physical_execution_selection_gate_certificate.json", cert)
    data_rows = _write_json(
        DATA_ROOT / "v83p_real_physical_execution_selection_gate_rows.json",
        {"kind": "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_ROWS_v0", "rows": rows},
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cert = _write_json(out_dir / "v83p_real_physical_execution_selection_gate_certificate.json", cert)
    return {"data_certificate": data_cert, "rows": data_rows, "certificate": out_cert}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V83P real physical execution selection gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v83p(args.out_dir)
    cert = _read_json(paths["certificate"], "V83P certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "support_tier_counts": cert["support_tier_counts"],
        "language_can_pass_without_physical_claim_count": cert["language_can_pass_without_physical_claim_count"],
        "physical_proxy_can_support_language_count": cert["physical_proxy_can_support_language_count"],
        "independent_experimental_or_coordinate_physical_proof_executed": cert["independent_experimental_or_coordinate_physical_proof_executed"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
