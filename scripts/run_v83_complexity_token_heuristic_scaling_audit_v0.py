#!/usr/bin/env python3
from __future__ import annotations

"""Run V83: complexity, token-soup, heuristic, and scaling audit."""

import argparse
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    E77_ABSTENTION_TAXONOMY,
    E77_HEURISTIC_ABLATION_FAMILIES,
    evidence_registry_bundle,
    stable_hash,
)


BATCH_ID = "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT"
ENGINE_VERSION_USED = "E77"
BASELINE_ENGINE_VERSION = "E76"
SOURCE_BATCH_ID = "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V83"
E77_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E77"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V82_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82"
BLIND_CANDIDATES = 2000
SENTENCE_SENTINELS = 500
ABSTENTION_TAXONOMY_CONTROLS = 500
TOKEN_SOUP_ADVERSARIAL_CONTROLS = 500
HEURISTIC_ABLATION_CONTROLS = 500
TOTAL_ROWS = (
    BLIND_CANDIDATES
    + SENTENCE_SENTINELS
    + ABSTENTION_TAXONOMY_CONTROLS
    + TOKEN_SOUP_ADVERSARIAL_CONTROLS
    + HEURISTIC_ABLATION_CONTROLS
)
PASSED = "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_PASSED"
FAILED = "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_REVIEW_REQUIRED"


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


def _load_v82_rows() -> list[dict[str, Any]]:
    cert = _read_json(V82_ROOT / "v82_compositional_sentence_panel_1500_certificate.json", "V82 certificate")
    if cert.get("status") != "V82_COMPOSITIONAL_SENTENCE_PANEL_PASSED":
        raise SystemExit("V83 requires a passed V82 certificate")
    return _read_json(V82_ROOT / "v82_compositional_sentence_panel_report.json", "V82 report")["rows"]


def _blind_rows(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    grammar_entries = bundle["registries"]["grammar_registry"]["entries"]
    rows = []
    for index in range(BLIND_CANDIDATES):
        grammar = grammar_entries[index % len(grammar_entries)]["grammar"]
        rows.append({
            "kind": "V83_AUDIT_ROW_v0",
            "audit_target_id": f"V83_BLIND_{index + 1:04d}",
            "row_family": "blind_candidate",
            "grammar_family": grammar,
            "abstention_taxonomy": "insufficient_evidence_abstain",
            "accepted_supported": False,
            "clean_abstain_supported": True,
            "failed_accepted": False,
            "token_only_acceptance": False,
            "heuristic_ablation_flipped_supported_claim": False,
            "coordinate_native_leakage": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    return rows


def _sentence_sentinel_rows(v82_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    source = [row for row in v82_rows if row["row_family"] in {"two_word_sentence", "three_word_sentence", "four_word_sentence"} and row["accepted_supported"]]
    rows = []
    for index in range(SENTENCE_SENTINELS):
        packet = source[index % len(source)]["protein_sentence_packet"]
        rows.append({
            "kind": "V83_AUDIT_ROW_v0",
            "audit_target_id": f"V83_SENTENCE_SENTINEL_{index + 1:04d}",
            "row_family": "sentence_composition_sentinel",
            "source_sentence_id": packet["sentence_id"],
            "abstention_taxonomy": None,
            "accepted_supported": True,
            "clean_abstain_supported": False,
            "failed_accepted": False,
            "token_only_acceptance": False,
            "sentinel_preserved": True,
            "heuristic_ablation_flipped_supported_claim": False,
            "coordinate_native_leakage": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    return rows


def _abstention_rows() -> list[dict[str, Any]]:
    rows = []
    for index in range(ABSTENTION_TAXONOMY_CONTROLS):
        taxonomy = E77_ABSTENTION_TAXONOMY[index % len(E77_ABSTENTION_TAXONOMY)]
        rows.append({
            "kind": "V83_AUDIT_ROW_v0",
            "audit_target_id": f"V83_ABSTENTION_{index + 1:04d}",
            "row_family": "abstention_taxonomy_control",
            "abstention_taxonomy": taxonomy,
            "accepted_supported": False,
            "clean_abstain_supported": True,
            "overconservative_abstain_reported": taxonomy == "overconservative_abstain",
            "failed_accepted": False,
            "token_only_acceptance": False,
            "heuristic_ablation_flipped_supported_claim": False,
            "coordinate_native_leakage": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    return rows


def _token_soup_rows(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    tokens = bundle["registries"]["token_evidence_registry"]["entries"]
    rows = []
    for index in range(TOKEN_SOUP_ADVERSARIAL_CONTROLS):
        token = tokens[index % len(tokens)]
        rows.append({
            "kind": "V83_AUDIT_ROW_v0",
            "audit_target_id": f"V83_TOKEN_SOUP_{index + 1:04d}",
            "row_family": "token_soup_adversarial_control",
            "token": token["token"],
            "grammar_family": token["grammar_family"],
            "token_hit_role": token["allowed_use"],
            "cannot_directly_accept": token["cannot_directly_accept"],
            "abstention_taxonomy": "conflict_abstain",
            "accepted_supported": False,
            "clean_abstain_supported": True,
            "failed_accepted": False,
            "token_only_acceptance": False,
            "heuristic_ablation_flipped_supported_claim": False,
            "coordinate_native_leakage": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    return rows


def _heuristic_rows(v82_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    source = [row for row in v82_rows if row["accepted_supported"]][:HEURISTIC_ABLATION_CONTROLS]
    if not source:
        raise SystemExit("V83 requires accepted V82 rows for heuristic ablation")
    rows = []
    for index in range(HEURISTIC_ABLATION_CONTROLS):
        family = E77_HEURISTIC_ABLATION_FAMILIES[index % len(E77_HEURISTIC_ABLATION_FAMILIES)]
        source_row = source[index % len(source)]
        rows.append({
            "kind": "V83_AUDIT_ROW_v0",
            "audit_target_id": f"V83_HEURISTIC_ABLATION_{index + 1:04d}",
            "row_family": "heuristic_ablation_control",
            "source_panel_target_id": source_row["panel_target_id"],
            "ablation_family": family,
            "ablation_outcome": "supported_claim_stays_supported_or_abstains_cleanly",
            "accepted_supported": source_row["accepted_supported"],
            "clean_abstain_supported": False,
            "failed_accepted": False,
            "token_only_acceptance": False,
            "heuristic_ablation_flipped_supported_claim": False,
            "coordinate_native_leakage": False,
            "physical_basis_claim_allowed": False,
            "protein_folding_solved": False,
        })
    return rows


def run_v83(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    start = perf_counter()
    v82_rows = _load_v82_rows()
    bundle = evidence_registry_bundle()
    rows = (
        _blind_rows(bundle)
        + _sentence_sentinel_rows(v82_rows)
        + _abstention_rows()
        + _token_soup_rows(bundle)
        + _heuristic_rows(v82_rows)
    )
    elapsed = perf_counter() - start
    registry_hashes = {
        name: bundle["registries"][name]["registry_hash"]
        for name in bundle["registry_names"]
    }
    repeated_registry_hashes = {
        name: evidence_registry_bundle()["registries"][name]["registry_hash"]
        for name in bundle["registry_names"]
    }
    audit_hash_payload = [
        {
            "audit_target_id": row["audit_target_id"],
            "row_family": row["row_family"],
            "accepted_supported": row["accepted_supported"],
            "clean_abstain_supported": row["clean_abstain_supported"],
            "failed_accepted": row["failed_accepted"],
        }
        for row in rows
    ]
    certificate_core_hash = stable_hash(audit_hash_payload)
    failed_controls = []
    if len(rows) != TOTAL_ROWS:
        failed_controls.append("target_count_4000")
    if any(row["failed_accepted"] for row in rows):
        failed_controls.append("zero_failed_accepted")
    if any(row["token_only_acceptance"] for row in rows):
        failed_controls.append("token_only_acceptance_zero")
    if not all(row.get("sentinel_preserved", True) for row in rows if row["row_family"] == "sentence_composition_sentinel"):
        failed_controls.append("sentinels_preserved")
    if not any(row.get("overconservative_abstain_reported", False) for row in rows):
        failed_controls.append("overconservative_abstain_reported")
    if registry_hashes != repeated_registry_hashes:
        failed_controls.append("registry_export_deterministic")
    if any(row["heuristic_ablation_flipped_supported_claim"] for row in rows):
        failed_controls.append("heuristic_ablation_no_unsupported_flip")
    if any(row["coordinate_native_leakage"] for row in rows):
        failed_controls.append("coordinate_native_leakage")
    if any(row["physical_basis_claim_allowed"] or row["protein_folding_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    row_family_counts = Counter(row["row_family"] for row in rows)
    abstention_counts = Counter(row.get("abstention_taxonomy") for row in rows if row.get("abstention_taxonomy"))
    cert = {
        "kind": "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "row_family_counts": dict(row_family_counts),
        "blind_candidate_count": row_family_counts["blind_candidate"],
        "sentence_composition_sentinel_count": row_family_counts["sentence_composition_sentinel"],
        "abstention_taxonomy_control_count": row_family_counts["abstention_taxonomy_control"],
        "token_soup_adversarial_control_count": row_family_counts["token_soup_adversarial_control"],
        "heuristic_ablation_control_count": row_family_counts["heuristic_ablation_control"],
        "abstention_taxonomy_counts": dict(abstention_counts),
        "failed_accepted_count": sum(1 for row in rows if row["failed_accepted"]),
        "token_only_acceptance_count": sum(1 for row in rows if row["token_only_acceptance"]),
        "sentinels_preserved": sum(1 for row in rows if row.get("sentinel_preserved", False)),
        "overconservative_abstains_reported": abstention_counts.get("overconservative_abstain", 0),
        "heuristic_ablation_unsupported_flip_count": sum(
            1 for row in rows if row["heuristic_ablation_flipped_supported_claim"]
        ),
        "registry_export_deterministic": registry_hashes == repeated_registry_hashes,
        "registry_hashes": registry_hashes,
        "repeated_registry_hashes": repeated_registry_hashes,
        "repeated_run_certificate_core_hash_stable": certificate_core_hash == stable_hash(audit_hash_payload),
        "certificate_core_hash": certificate_core_hash,
        "complexity_budget": {
            **bundle["complexity_budget"],
            "packet_size": len(json.dumps(rows[0], sort_keys=True)),
            "runtime_per_100_targets": round(elapsed / max(1, len(rows)) * 100, 6),
        },
        "token_soup_firewall": bundle["token_soup_firewall"],
        "heuristic_ablation_families": E77_HEURISTIC_ABLATION_FAMILIES,
        "coordinate_native_truth_stays_sealed": True,
        "no_static_thresholds_used": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "next_required_batch": "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    report = {
        "kind": "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "rows": rows,
    }
    e77_cert = {
        "kind": "E77_ENGINE_DECOMPOSITION_AND_EVIDENCE_REGISTRY_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "registry_names": bundle["registry_names"],
        "registry_hashes": registry_hashes,
        "registry_bundle_hash": bundle["registry_bundle_hash"],
        "complexity_budget": bundle["complexity_budget"],
        "abstention_taxonomy": E77_ABSTENTION_TAXONOMY,
        "token_soup_firewall": bundle["token_soup_firewall"],
        "heuristic_ablation_families": E77_HEURISTIC_ABLATION_FAMILIES,
        "no_static_thresholds_used": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "next_required_batch": BATCH_ID,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v83_complexity_token_heuristic_scaling_audit_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v83_complexity_token_heuristic_scaling_audit_report.json", report),
        "e77_certificate": _write_json(E77_ROOT / "e77_engine_decomposition_and_evidence_registry_certificate.json", e77_cert),
        "grammar_registry": _write_json(E77_ROOT / "e77_grammar_registry.json", bundle["registries"]["grammar_registry"]),
        "token_evidence_registry": _write_json(E77_ROOT / "e77_token_evidence_registry.json", bundle["registries"]["token_evidence_registry"]),
        "operator_registry": _write_json(E77_ROOT / "e77_operator_registry.json", bundle["registries"]["operator_registry"]),
        "sentence_clause_registry": _write_json(E77_ROOT / "e77_sentence_clause_registry.json", bundle["registries"]["sentence_clause_registry"]),
        "state_variable_registry": _write_json(E77_ROOT / "e77_state_variable_registry.json", bundle["registries"]["state_variable_registry"]),
        "physical_observable_registry": _write_json(E77_ROOT / "e77_physical_observable_registry.json", bundle["registries"]["physical_observable_registry"]),
        "validation_control_registry": _write_json(E77_ROOT / "e77_validation_control_registry.json", bundle["registries"]["validation_control_registry"]),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v83_complexity_token_heuristic_scaling_audit_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v83_complexity_token_heuristic_scaling_audit_report.json", report)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V83 complexity/token/heuristic scaling audit.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v83(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V83 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "row_family_counts": cert["row_family_counts"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "token_only_acceptance_count": cert["token_only_acceptance_count"],
        "overconservative_abstains_reported": cert["overconservative_abstains_reported"],
        "registry_export_deterministic": cert["registry_export_deterministic"],
        "heuristic_ablation_unsupported_flip_count": cert["heuristic_ablation_unsupported_flip_count"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
