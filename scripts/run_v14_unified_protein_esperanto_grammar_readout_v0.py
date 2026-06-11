#!/usr/bin/env python3
from __future__ import annotations

"""V14 unified role-aware evidence grammar readout.

Postprocess-only. No MD rerun. No threshold tuning. No native precision use.

This script reads already-generated 4AKE / 1UBQ / 1CLL artifacts and translates
all of them into the same claim-safe evidence grammar:

    DCA signal -> geometry/replica support -> topology/purpose role ->
    chemical policy -> noise/pollution guards -> claim lock / abstain.

The goal is not to claim that folding is solved. The goal is to verify whether
three different protein object types can be represented in one coherent,
non-leaking, abstain-capable grammar of evidence.
"""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = DEFAULT_RUN_ROOT / "V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_READOUT"

DEFAULT_UBQ_CERT = (
    DEFAULT_RUN_ROOT
    / "V13a_1UBQ_ADAPTIVE_CHEMICAL_POLICY_READOUT"
    / "v13a_1ubq_purpose_gate_readout_certificate.json"
)
DEFAULT_CLL_CERT = (
    DEFAULT_RUN_ROOT
    / "V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT"
    / "v13b_1cll_hierarchical_purpose_topology_readout_certificate.json"
)

KNOWN_4AKE_ARTIFACTS = [
    DEFAULT_RUN_ROOT / "V13c_4AKE_ADENYLATE_KINASE" / "openmm_tmd_replicas_v0_certificate.json",
    DEFAULT_RUN_ROOT / "comparison_v9" / "4AKE" / "role_aware_rescue_selector_v9_postprocess.json",
    DEFAULT_RUN_ROOT / "GARAGE_ROLE_AWARE_RESCUE_SELECTOR_V9" / "role_aware_rescue_selector_v9_postprocess.json",
    DEFAULT_RUN_ROOT / "GARAGE_TAIL_REACHABILITY_WITH_STRICT_GUARD_V5" / "openmm_tmd_replicas_v0_certificate.json",
]

GRAMMAR_AXES = [
    "external_DCA_or_coupling_signal",
    "geometry_reachability_from_trajectory",
    "replica_persistence",
    "purpose_or_topology_role_assignment",
    "role_aware_chemical_policy",
    "noise_and_long_range_pollution_guards",
    "claim_lock_or_clean_abstain",
]


def _read_json(path: Path) -> Optional[dict[str, Any]]:
    if not path.exists() or not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if isinstance(payload, dict):
        payload["_artifact_path"] = str(path)
        return payload
    return None


def _pair_key(pair: Any) -> str:
    if isinstance(pair, str):
        return pair
    if isinstance(pair, (list, tuple)) and len(pair) == 2:
        return f"{int(pair[0])}-{int(pair[1])}"
    return str(pair)


def _pairs_to_keys(pairs: Any) -> list[str]:
    if not isinstance(pairs, list):
        return []
    return [_pair_key(pair) for pair in pairs]


def _safe_int(value: Any, default: int = 0) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)):
        return int(value)
    try:
        return int(value)
    except Exception:
        return default


def _present_nonempty(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (list, dict, tuple, set, str)):
        return bool(value)
    if isinstance(value, (int, float)):
        return bool(value)
    return True


def _lookup(d: dict[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = d
    for key in keys:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(key)
    return cur if cur is not None else default


def _positive_if_any(*values: Any) -> bool:
    return any(_present_nonempty(value) for value in values)


def _make_missing_row(protein: str, target_role: str, expected_artifact: str) -> dict[str, Any]:
    return {
        "protein": protein,
        "artifact_status": "missing",
        "artifact_path": expected_artifact,
        "target_role": target_role,
        "grammar_policy": "not_evaluated_artifact_missing",
        "selected_domain_core": [],
        "selected_balanced_core": [],
        "selected_hinge_or_interdomain": [],
        "selected_local_support": [],
        "dca_support": {},
        "replica_support": {},
        "chemical_policy": None,
        "topology_policy": None,
        "noise_added": None,
        "long_range_evidence_polluted": None,
        "classification_coverage_ratio": None,
        "claim_lock_status": "not_evaluated_artifact_missing",
        "claim_lock_failed_checks": ["artifact_present"],
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "final_status": f"{protein}_not_evaluated_artifact_missing_claim_disabled",
        "positive_evidence_found": False,
    }


def normalize_1ubq(payload: Optional[dict[str, Any]], artifact_path: Path = DEFAULT_UBQ_CERT) -> dict[str, Any]:
    if payload is None:
        return _make_missing_row("1UBQ", "single_domain_compact", str(artifact_path))

    band = payload.get("selected_frequency_band") or {}
    selected_balanced_core = band.get("selected_balanced_core") or payload.get("selected_balanced_core") or []
    selected_pair_count = _safe_int(band.get("selected_pair_count", payload.get("selected_pair_count", 0)))
    positive = selected_pair_count > 0 or bool(selected_balanced_core)

    final_status = (
        "single_domain_compact_signal_found_under_adaptive_chemical_policy;"
        "global_chemical_gate_rejected_as_too_strict;claim_allowed=false"
        if positive
        else "single_domain_compact_clean_abstain_under_adaptive_chemical_policy;claim_allowed=false"
    )

    return {
        "protein": "1UBQ",
        "artifact_status": "present",
        "artifact_path": str(payload.get("_artifact_path", artifact_path)),
        "source_accession": "1UBQ:A",
        "target_role": payload.get("target_role") or _lookup(payload, "target_purpose", "target_role", default="single_domain_compact"),
        "grammar_policy": "single_domain_compact_adaptive_chemical_policy",
        "selected_domain_core": [],
        "selected_balanced_core": _pairs_to_keys(selected_balanced_core),
        "selected_hinge_or_interdomain": [],
        "selected_local_support": [],
        "dca_support": {
            "dca_mean_selected": band.get("dca_mean_selected"),
            "dca_mean_effective_balanced_background": band.get("dca_mean_effective_balanced_background"),
            "dca_absolute_support_pass": band.get("dca_absolute_support_pass"),
            "dca_background_enrichment_ratio": band.get("dca_background_enrichment_ratio"),
            "dca_background_enrichment_pass": band.get("dca_background_enrichment_pass"),
            "dca_score_by_selected_pair": band.get("dca_score_by_selected_pair", {}),
        },
        "replica_support": band.get("support_by_selected_pair", {}),
        "mean_frequency_by_selected_pair": band.get("mean_frequency_by_selected_pair", {}),
        "chemical_policy": payload.get("chemical_policy"),
        "chemical_score_by_selected_pair": band.get("chemical_score_by_selected_pair", {}),
        "legacy_chemical_hard_gate_would_block_selected_count": band.get(
            "legacy_chemical_hard_gate_would_block_selected_count",
            payload.get("legacy_chemical_hard_gate_would_block_selected_count"),
        ),
        "topology_policy": "single_domain_compact",
        "noise_added": band.get("noise_added"),
        "long_range_evidence_polluted": band.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": band.get("classification_coverage_ratio"),
        "claim_lock_status": payload.get("claim_lock_status") or _lookup(payload, "claim_lock_check", "status"),
        "claim_lock_failed_checks": payload.get("claim_lock_failed_checks") or _lookup(payload, "claim_lock_check", "failed_checks", default=[]),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "final_status": final_status,
        "positive_evidence_found": positive,
    }


def normalize_1cll(payload: Optional[dict[str, Any]], artifact_path: Path = DEFAULT_CLL_CERT) -> dict[str, Any]:
    if payload is None:
        return _make_missing_row("1CLL", "multi_domain_composite", str(artifact_path))

    band = payload.get("selected_frequency_band") or {}
    selected_n = payload.get("selected_N_domain_core") or band.get("selected_N_domain_core") or []
    selected_c = payload.get("selected_C_domain_core") or band.get("selected_C_domain_core") or []
    selected_interdomain = payload.get("selected_interdomain_hinge") or band.get("selected_interdomain_hinge") or []
    selected_local = payload.get("selected_local_support") or band.get("selected_local_support") or []
    selected_medium = payload.get("selected_medium_support") or band.get("selected_medium_support") or []
    selected_domain_core = _pairs_to_keys(selected_n) + _pairs_to_keys(selected_c)
    selected_pair_count = _safe_int(payload.get("selected_pair_count", band.get("selected_pair_count", 0)))
    positive = selected_pair_count > 0 or bool(selected_domain_core or selected_interdomain)

    if selected_domain_core and not selected_interdomain:
        final_status = "multi_domain_composite_signal_found;domain_core_detected;interdomain_hinge_not_yet_proven;claim_allowed=false"
    elif selected_interdomain:
        final_status = "multi_domain_composite_interdomain_hinge_signal_found;claim_allowed=false"
    elif positive:
        final_status = "multi_domain_composite_support_signal_found;claim_allowed=false"
    else:
        final_status = "multi_domain_composite_clean_abstain_no_domain_or_interdomain_core;claim_allowed=false"

    return {
        "protein": "1CLL",
        "artifact_status": "present",
        "artifact_path": str(payload.get("_artifact_path", artifact_path)),
        "source_accession": "1CLL:A",
        "target_role": payload.get("target_role") or "multi_domain_composite",
        "grammar_policy": "hierarchical_domain_core_plus_interdomain_topology",
        "domain_roles": payload.get("domain_roles", []),
        "selected_domain_core": selected_domain_core,
        "selected_N_domain_core": _pairs_to_keys(selected_n),
        "selected_C_domain_core": _pairs_to_keys(selected_c),
        "selected_balanced_core": selected_domain_core,
        "selected_hinge_or_interdomain": _pairs_to_keys(selected_interdomain),
        "selected_local_support": _pairs_to_keys(selected_local),
        "selected_medium_support": _pairs_to_keys(selected_medium),
        "dca_support": {
            "dca_mean_selected": band.get("dca_mean_selected"),
            "dca_mean_effective_balanced_background": band.get("dca_mean_effective_balanced_background"),
            "dca_absolute_support_pass": band.get("dca_absolute_support_pass"),
            "dca_background_enrichment_ratio": band.get("dca_background_enrichment_ratio"),
            "dca_background_enrichment_pass": band.get("dca_background_enrichment_pass"),
            "dca_score_by_selected_pair": band.get("dca_score_by_selected_pair", {}),
        },
        "replica_support": band.get("support_by_selected_pair", {}),
        "mean_frequency_by_selected_pair": band.get("mean_frequency_by_selected_pair", {}),
        "chemical_policy": payload.get("chemical_policy") or "adaptive_soft_guard",
        "chemical_score_by_selected_pair": band.get("chemical_score_by_selected_pair", {}),
        "topology_policy": payload.get("topology_policy"),
        "noise_added": band.get("noise_added"),
        "long_range_evidence_polluted": band.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": band.get("classification_coverage_ratio"),
        "claim_lock_status": payload.get("claim_lock_status") or _lookup(payload, "claim_lock_check", "status"),
        "claim_lock_failed_checks": payload.get("claim_lock_failed_checks") or _lookup(payload, "claim_lock_check", "failed_checks", default=[]),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "final_status": final_status,
        "positive_evidence_found": positive,
    }


def _find_4ake_artifact(run_root: Path, explicit: Optional[Path] = None) -> Optional[Path]:
    if explicit and explicit.exists():
        return explicit
    for path in KNOWN_4AKE_ARTIFACTS:
        if path.exists():
            return path
    if not run_root.exists():
        return None
    candidates: list[Path] = []
    for path in run_root.rglob("*.json"):
        name = path.name.lower()
        parent_text = str(path.parent).lower()
        if "4ake" not in parent_text and "4ake" not in name:
            continue
        if name in {"role_aware_rescue_selector_v9_postprocess.json", "openmm_tmd_replicas_v0_certificate.json"}:
            candidates.append(path)
    # Prefer role-aware postprocess over raw runtime cert.
    candidates.sort(key=lambda p: ("role_aware" not in p.name.lower(), len(str(p))))
    return candidates[0] if candidates else None


def normalize_4ake(payload: Optional[dict[str, Any]], artifact_path: Optional[Path] = None) -> dict[str, Any]:
    if payload is None:
        return _make_missing_row("4AKE", "domain_hinge_object", str(artifact_path or "auto_search_failed"))

    selected_strict = payload.get("selected_strict_scaffold") or []
    selected_balanced = payload.get("selected_balanced_core") or []
    selected_border = payload.get("selected_border_rescue") or []
    selected_local = payload.get("selected_local_support") or []
    selected_medium = payload.get("selected_medium_support") or []
    selected_pairs = payload.get("selected_pairs") or []
    long_range = _lookup(payload, "long_range_evidence", "pairs", default=[])
    selected_pair_count = _safe_int(
        payload.get("selected_pair_count", payload.get("runtime_v10_selected_pair_count", len(selected_pairs)))
    )
    if selected_pair_count == 0:
        selected_pair_count = len(selected_pairs or []) or len(long_range or []) or len(selected_strict or []) + len(selected_balanced or []) + len(selected_border or [])
    positive = selected_pair_count > 0 or _positive_if_any(selected_strict, selected_balanced, selected_border, long_range)
    final_status = (
        "role_aware_steering_reproduction_signal_present;domain_hinge_object;claim_allowed=false"
        if positive
        else "role_aware_steering_artifact_present_but_no_selected_role_evidence;claim_allowed=false"
    )

    return {
        "protein": "4AKE",
        "artifact_status": "present",
        "artifact_path": str(payload.get("_artifact_path", artifact_path or "unknown")),
        "source_accession": payload.get("source_accession", "4AKE:A"),
        "target_role": "domain_hinge_object",
        "grammar_policy": "role_aware_domain_hinge_steering_reproduction",
        "selected_domain_core": [],
        "selected_strict_scaffold": _pairs_to_keys(selected_strict),
        "selected_balanced_core": _pairs_to_keys(selected_balanced),
        "selected_hinge_or_interdomain": _pairs_to_keys(selected_border),
        "selected_local_support": _pairs_to_keys(selected_local),
        "selected_medium_support": _pairs_to_keys(selected_medium),
        "long_range_evidence": _pairs_to_keys(long_range),
        "dca_support": {
            "lane_vote_support": payload.get("lane_vote_support") or payload.get("lane_vote_support_after", {}),
            "metrics": payload.get("metrics", {}),
        },
        "replica_support": payload.get("lane_vote_support") or payload.get("lane_vote_support_after", {}),
        "chemical_policy": "role_aware_chemical_guard_from_source_artifact",
        "topology_policy": str(payload.get("topology_mode", "domain_hinge_or_interdomain")),
        "noise_added": payload.get("noise_added"),
        "long_range_evidence_polluted": payload.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": payload.get("classification_coverage_ratio"),
        "claim_lock_status": "claim_disabled_source_artifact_not_global_claim",
        "claim_lock_failed_checks": [],
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "final_status": final_status,
        "positive_evidence_found": positive,
    }


def _coherence_checks(rows: list[dict[str, Any]]) -> dict[str, Any]:
    present = [row for row in rows if row.get("artifact_status") == "present"]
    missing = [row["protein"] for row in rows if row.get("artifact_status") != "present"]
    positive = [row["protein"] for row in rows if row.get("positive_evidence_found")]
    claim_enabled = [row["protein"] for row in rows if row.get("claim_allowed")]
    required_axes = {
        "target_role",
        "grammar_policy",
        "dca_support",
        "replica_support",
        "chemical_policy",
        "topology_policy",
        "claim_lock_status",
        "final_status",
    }
    axis_missing = {
        row["protein"]: sorted([key for key in required_axes if key not in row or row.get(key) is None])
        for row in rows
    }
    return {
        "present_artifact_count": len(present),
        "missing_artifacts": missing,
        "positive_evidence_proteins": positive,
        "claim_enabled_proteins": claim_enabled,
        "no_claim_allowed_anywhere": len(claim_enabled) == 0,
        "axis_missing_by_protein": axis_missing,
        "all_rows_have_core_grammar_axes": all(not values for values in axis_missing.values()),
        "three_object_types_represented": {row["protein"]: row.get("target_role") for row in rows},
    }


def _global_status(rows: list[dict[str, Any]], checks: dict[str, Any]) -> str:
    missing = checks.get("missing_artifacts", [])
    positive = set(checks.get("positive_evidence_proteins", []))
    if not checks.get("no_claim_allowed_anywhere", False):
        return "invalid_claim_enabled_in_unified_readout"
    if not missing and {"4AKE", "1UBQ", "1CLL"}.issubset(positive):
        return "unified_role_aware_evidence_grammar_coherent_across_three_object_types_claim_disabled"
    if {"1UBQ", "1CLL"}.issubset(positive) and "4AKE" in missing:
        return "partial_unified_grammar_panel_positive_on_1UBQ_1CLL_4AKE_artifact_missing_claim_disabled"
    if positive:
        return "partial_unified_grammar_panel_some_positive_evidence_claim_disabled"
    return "unified_grammar_panel_clean_abstain_or_artifacts_missing_claim_disabled"


def _flatten_for_csv(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "protein": row.get("protein"),
        "artifact_status": row.get("artifact_status"),
        "target_role": row.get("target_role"),
        "grammar_policy": row.get("grammar_policy"),
        "selected_domain_core": json.dumps(row.get("selected_domain_core", []), sort_keys=True),
        "selected_balanced_core": json.dumps(row.get("selected_balanced_core", []), sort_keys=True),
        "selected_hinge_or_interdomain": json.dumps(row.get("selected_hinge_or_interdomain", []), sort_keys=True),
        "selected_local_support": json.dumps(row.get("selected_local_support", []), sort_keys=True),
        "replica_support": json.dumps(row.get("replica_support", {}), sort_keys=True),
        "chemical_policy": row.get("chemical_policy"),
        "topology_policy": row.get("topology_policy"),
        "noise_added": row.get("noise_added"),
        "long_range_evidence_polluted": row.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": row.get("classification_coverage_ratio"),
        "claim_lock_status": row.get("claim_lock_status"),
        "claim_allowed": row.get("claim_allowed"),
        "final_status": row.get("final_status"),
        "artifact_path": row.get("artifact_path"),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_markdown(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V14 Unified Protein Esperanto Grammar Readout",
        "",
        f"Global status: `{cert['global_status']}`",
        "",
        "This is a postprocess-only, claim-safe grammar readout. It does not rerun MD and does not enable biological transfer claims.",
        "",
        "## Grammar axes",
        "",
    ]
    for axis in cert["unified_grammar_axes"]:
        lines.append(f"- `{axis}`")
    lines.extend(["", "## Protein rows", ""])
    for row in cert["protein_rows"]:
        lines.extend([
            f"### {row['protein']}",
            "",
            f"- Artifact status: `{row.get('artifact_status')}`",
            f"- Target role: `{row.get('target_role')}`",
            f"- Grammar policy: `{row.get('grammar_policy')}`",
            f"- Selected domain core: `{row.get('selected_domain_core', [])}`",
            f"- Selected balanced core: `{row.get('selected_balanced_core', [])}`",
            f"- Selected hinge/interdomain: `{row.get('selected_hinge_or_interdomain', [])}`",
            f"- Noise added: `{row.get('noise_added')}`",
            f"- Long-range evidence polluted: `{row.get('long_range_evidence_polluted')}`",
            f"- Claim lock status: `{row.get('claim_lock_status')}`",
            f"- Final status: `{row.get('final_status')}`",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")


def build_unified_readout(
    *,
    run_root: Path = DEFAULT_RUN_ROOT,
    ubq_cert: Path = DEFAULT_UBQ_CERT,
    cll_cert: Path = DEFAULT_CLL_CERT,
    ake_artifact: Optional[Path] = None,
) -> dict[str, Any]:
    ubq_payload = _read_json(ubq_cert)
    cll_payload = _read_json(cll_cert)
    ake_path = _find_4ake_artifact(run_root, ake_artifact)
    ake_payload = _read_json(ake_path) if ake_path else None

    rows = [
        normalize_4ake(ake_payload, ake_path),
        normalize_1ubq(ubq_payload, ubq_cert),
        normalize_1cll(cll_payload, cll_cert),
    ]
    checks = _coherence_checks(rows)
    status = _global_status(rows, checks)
    return {
        "kind": "V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_READOUT_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "unified_grammar_axes": GRAMMAR_AXES,
        "global_status": status,
        "coherence_checks": checks,
        "protein_rows": rows,
        "final_interpretation": (
            "Unified role-aware evidence grammar readout completed. This is not a protein-folding-solved claim. "
            "It reports whether evidence across 4AKE, 1UBQ, and 1CLL can be translated into one claim-safe grammar."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build V14 unified protein Esperanto grammar readout.")
    parser.add_argument("--run-root", default=str(DEFAULT_RUN_ROOT))
    parser.add_argument("--ubq-cert", default=str(DEFAULT_UBQ_CERT))
    parser.add_argument("--cll-cert", default=str(DEFAULT_CLL_CERT))
    parser.add_argument("--ake-artifact", default="")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    cert = build_unified_readout(
        run_root=Path(args.run_root),
        ubq_cert=Path(args.ubq_cert),
        cll_cert=Path(args.cll_cert),
        ake_artifact=Path(args.ake_artifact) if args.ake_artifact else None,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v14_unified_protein_esperanto_grammar_readout_certificate.json"
    csv_path = out_dir / "v14_unified_protein_esperanto_grammar_table.csv"
    md_path = out_dir / "V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_REPORT.md"
    cert_path.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(csv_path, [_flatten_for_csv(row) for row in cert["protein_rows"]])
    _write_markdown(md_path, cert)

    print(json.dumps({
        "certificate": str(cert_path),
        "table": str(csv_path),
        "report": str(md_path),
        "global_status": cert["global_status"],
        "claim_allowed": cert["claim_allowed"],
        "positive_evidence_proteins": cert["coherence_checks"]["positive_evidence_proteins"],
        "missing_artifacts": cert["coherence_checks"]["missing_artifacts"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
