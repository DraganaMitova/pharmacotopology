#!/usr/bin/env python3
from __future__ import annotations

"""V15 dynamic separation grammar readout.

Postprocess-only.  No MD rerun.  No native precision use.  No fixed residue
count decides whether a pair is valid evidence.

Earlier layers used sequence distance as a flat proxy for "non-local" signal.
V15 retires that as a claim/control rule.  Sequence separation remains reported
as descriptive metadata, but role assignment is driven by contextual evidence:
protein purpose, domain relation, DCA support, geometry/trajectory support,
replica persistence, adaptive chemical policy, and noise/pollution guards.
"""

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = DEFAULT_RUN_ROOT / "V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT"
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
DEFAULT_V14_CERT = (
    DEFAULT_RUN_ROOT
    / "V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_READOUT"
    / "v14_unified_protein_esperanto_grammar_readout_certificate.json"
)
DEFAULT_4AKE_BRIDGE_CERT = (
    DEFAULT_RUN_ROOT
    / "V15_4AKE_DYNAMIC_GRAMMAR_BRIDGE"
    / "v15_4ake_dynamic_grammar_bridge_certificate.json"
)

GRAMMAR_AXES = [
    "external_DCA_or_coupling_signal",
    "geometry_reachability_from_trajectory",
    "replica_persistence",
    "purpose_or_topology_role_assignment",
    "role_aware_chemical_policy",
    "dynamic_separation_context_no_fixed_residue_cutoff",
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


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _pair_key(pair: Any) -> str:
    if isinstance(pair, str):
        return pair
    if isinstance(pair, (list, tuple)) and len(pair) == 2:
        return f"{int(pair[0])}-{int(pair[1])}"
    return str(pair)


def _parse_pair_key(pair: Any) -> Optional[tuple[int, int]]:
    if isinstance(pair, (list, tuple)) and len(pair) == 2:
        left, right = int(pair[0]), int(pair[1])
        return (left, right) if left < right else (right, left)
    if isinstance(pair, str) and "-" in pair:
        left_text, right_text = pair.split("-", maxsplit=1)
        try:
            left, right = int(left_text), int(right_text)
        except ValueError:
            return None
        return (left, right) if left < right else (right, left)
    return None


def _pairs_to_keys(pairs: Any) -> list[str]:
    if not isinstance(pairs, list):
        return []
    return [_pair_key(pair) for pair in pairs]


def _lookup(d: dict[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = d
    for key in keys:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(key)
    return cur if cur is not None else default


def _safe_int(value: Any, default: int = 0) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)):
        return int(value)
    try:
        return int(value)
    except Exception:
        return default


def _sequence_context(
    pair: tuple[int, int],
    *,
    sequence_length: Optional[int],
    domain_relation: str,
) -> dict[str, Any]:
    left, right = pair
    sep = right - left
    denom = max(1, (sequence_length or right) - 1)
    return {
        "sequence_separation": sep,
        "normalized_sequence_separation": round(sep / float(denom), 6),
        "domain_relation": domain_relation,
        "separation_filter_applied": False,
        "fixed_residue_cutoff_used": False,
        "interpretation": "descriptive_context_only_not_a_gate",
    }


def normalize_1ubq_dynamic(payload: Optional[dict[str, Any]], artifact_path: Path = DEFAULT_UBQ_CERT) -> dict[str, Any]:
    if payload is None:
        return {
            "protein": "1UBQ",
            "artifact_status": "missing",
            "target_role": "single_domain_compact",
            "grammar_policy": "dynamic_separation_not_evaluated_artifact_missing",
            "selected_pairs": [],
            "dynamic_pair_roles": {},
            "claim_allowed": False,
            "positive_evidence_found": False,
            "final_status": "1UBQ_dynamic_separation_not_evaluated_artifact_missing_claim_disabled",
        }

    band = payload.get("selected_frequency_band") or {}
    selected_pairs = _pairs_to_keys(band.get("selected_balanced_core") or payload.get("selected_balanced_core") or [])
    support = band.get("support_by_selected_pair", {}) if isinstance(band.get("support_by_selected_pair"), dict) else {}
    mean_freq = band.get("mean_frequency_by_selected_pair", {}) if isinstance(band.get("mean_frequency_by_selected_pair"), dict) else {}
    chem = band.get("chemical_score_by_selected_pair", {}) if isinstance(band.get("chemical_score_by_selected_pair"), dict) else {}
    dca = band.get("dca_score_by_selected_pair", {}) if isinstance(band.get("dca_score_by_selected_pair"), dict) else {}
    sequence_length = _safe_int(payload.get("sequence_length"), 76)

    pair_roles: dict[str, dict[str, Any]] = {}
    for key in selected_pairs:
        pair = _parse_pair_key(key)
        if pair is None:
            continue
        pair_roles[key] = {
            **_sequence_context(pair, sequence_length=sequence_length, domain_relation="single_domain"),
            "role_decision": "single_domain_compact_balanced_core",
            "evidence_class": "single_domain_compact_core_evidence",
            "selected": True,
            "support": support.get(key),
            "mean_frequency": mean_freq.get(key),
            "chemical_score": chem.get(key),
            "dca_score": dca.get(key) or band.get("dca_mean_selected"),
            "role_assignment_basis": [
                "single_domain_compact_target_role",
                "selected_balanced_core_evidence",
                "replica_support",
                "trajectory_frequency",
                "DCA_support",
                "adaptive_chemical_policy",
            ],
        }

    positive = bool(selected_pairs)
    return {
        "protein": "1UBQ",
        "artifact_status": "present",
        "artifact_path": str(payload.get("_artifact_path", artifact_path)),
        "target_role": payload.get("target_role") or _lookup(payload, "target_purpose", "target_role", default="single_domain_compact"),
        "grammar_policy": "single_domain_compact_dynamic_separation_adaptive_chemical_policy",
        "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
        "fixed_residue_cutoff_used": False,
        "selected_pairs": selected_pairs,
        "selected_balanced_core": selected_pairs,
        "dynamic_pair_roles": pair_roles,
        "replica_support": support,
        "mean_frequency_by_selected_pair": mean_freq,
        "chemical_policy": payload.get("chemical_policy"),
        "chemical_score_by_selected_pair": chem,
        "dca_support": {
            "dca_score_by_selected_pair": dca,
            "dca_mean_selected": band.get("dca_mean_selected"),
            "dca_background_enrichment_ratio": band.get("dca_background_enrichment_ratio"),
            "dca_background_enrichment_pass": band.get("dca_background_enrichment_pass"),
            "dca_absolute_support_pass": band.get("dca_absolute_support_pass"),
        },
        "noise_added": band.get("noise_added"),
        "long_range_evidence_polluted": band.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": band.get("classification_coverage_ratio"),
        "claim_lock_status": payload.get("claim_lock_status") or _lookup(payload, "claim_lock_check", "status"),
        "claim_lock_failed_checks": payload.get("claim_lock_failed_checks") or _lookup(payload, "claim_lock_check", "failed_checks", default=[]),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "positive_evidence_found": positive,
        "final_status": (
            "single_domain_compact_signal_found_with_dynamic_separation_context;claim_allowed=false"
            if positive
            else "single_domain_compact_clean_abstain_with_dynamic_separation_context;claim_allowed=false"
        ),
    }


def normalize_1cll_dynamic(payload: Optional[dict[str, Any]], artifact_path: Path = DEFAULT_CLL_CERT) -> dict[str, Any]:
    if payload is None:
        return {
            "protein": "1CLL",
            "artifact_status": "missing",
            "target_role": "multi_domain_composite",
            "grammar_policy": "dynamic_separation_not_evaluated_artifact_missing",
            "selected_pairs": [],
            "dynamic_pair_roles": {},
            "claim_allowed": False,
            "positive_evidence_found": False,
            "final_status": "1CLL_dynamic_separation_not_evaluated_artifact_missing_claim_disabled",
        }

    band = payload.get("selected_frequency_band") or {}
    selected_n = payload.get("selected_N_domain_core") or band.get("selected_N_domain_core") or []
    selected_c = payload.get("selected_C_domain_core") or band.get("selected_C_domain_core") or []
    selected_interdomain = payload.get("selected_interdomain_hinge") or band.get("selected_interdomain_hinge") or []
    selected_local = payload.get("selected_local_support") or band.get("selected_local_support") or []
    selected_medium = payload.get("selected_medium_support") or band.get("selected_medium_support") or []
    selected_pairs = _pairs_to_keys(selected_n) + _pairs_to_keys(selected_c) + _pairs_to_keys(selected_interdomain) + _pairs_to_keys(selected_local) + _pairs_to_keys(selected_medium)

    role_by_candidate_pair = payload.get("role_by_candidate_pair") if isinstance(payload.get("role_by_candidate_pair"), dict) else {}
    support = band.get("support_by_selected_pair", {}) if isinstance(band.get("support_by_selected_pair"), dict) else {}
    mean_freq = band.get("mean_frequency_by_selected_pair", {}) if isinstance(band.get("mean_frequency_by_selected_pair"), dict) else {}
    chem = band.get("chemical_score_by_selected_pair", {}) if isinstance(band.get("chemical_score_by_selected_pair"), dict) else {}
    dca = band.get("dca_score_by_selected_pair", {}) if isinstance(band.get("dca_score_by_selected_pair"), dict) else {}
    sequence_length = _safe_int(payload.get("sequence_length"), 148)

    pair_roles: dict[str, dict[str, Any]] = {}
    all_keys = sorted(set(role_by_candidate_pair) | set(selected_pairs))
    for key in all_keys:
        pair = _parse_pair_key(key)
        if pair is None:
            continue
        role_payload = role_by_candidate_pair.get(key, {}) if isinstance(role_by_candidate_pair.get(key), dict) else {}
        domain_relation = str(role_payload.get("domain_relation", "unknown_domain_relation"))
        evidence_class = str(role_payload.get("evidence_class", "unselected_candidate"))
        topology_role = str(role_payload.get("topology_role", "unassigned_dynamic_context"))
        selected = key in selected_pairs
        pair_roles[key] = {
            **_sequence_context(pair, sequence_length=sequence_length, domain_relation=domain_relation),
            "topology_role": topology_role,
            "evidence_class": evidence_class,
            "role_decision": topology_role,
            "selected": selected,
            "support": support.get(key),
            "mean_frequency": mean_freq.get(key) or role_payload.get("tail_frequency_mean"),
            "tail_frequency_min": role_payload.get("tail_frequency_min"),
            "tail_frequency_max": role_payload.get("tail_frequency_max"),
            "chemical_score": chem.get(key) if key in chem else role_payload.get("chemical_score"),
            "dca_score": dca.get(key) if key in dca else role_payload.get("dca_score"),
            "inside_effective_balanced": role_payload.get("inside_effective_balanced"),
            "role_assignment_basis": [
                "declared_domain_boundaries",
                "domain_relation",
                "hierarchical_topology_policy",
                "replica_support_or_candidate_context",
                "DCA_support",
                "adaptive_chemical_policy",
            ],
        }

    selected_domain_core = _pairs_to_keys(selected_n) + _pairs_to_keys(selected_c)
    positive = bool(selected_pairs)
    if selected_domain_core and not _pairs_to_keys(selected_interdomain):
        final_status = "multi_domain_composite_domain_core_signal_found_with_dynamic_separation_context;interdomain_hinge_not_yet_proven;claim_allowed=false"
    elif _pairs_to_keys(selected_interdomain):
        final_status = "multi_domain_composite_interdomain_signal_found_with_dynamic_separation_context;claim_allowed=false"
    elif positive:
        final_status = "multi_domain_composite_support_signal_found_with_dynamic_separation_context;claim_allowed=false"
    else:
        final_status = "multi_domain_composite_clean_abstain_with_dynamic_separation_context;claim_allowed=false"

    return {
        "protein": "1CLL",
        "artifact_status": "present",
        "artifact_path": str(payload.get("_artifact_path", artifact_path)),
        "target_role": payload.get("target_role", "multi_domain_composite"),
        "grammar_policy": "multi_domain_composite_dynamic_separation_hierarchical_topology",
        "topology_policy": payload.get("topology_policy", "hierarchical_domain_core_plus_interdomain"),
        "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
        "fixed_residue_cutoff_used": False,
        "domain_boundaries": payload.get("domain_boundaries"),
        "selected_pairs": selected_pairs,
        "selected_domain_core": selected_domain_core,
        "selected_N_domain_core": _pairs_to_keys(selected_n),
        "selected_C_domain_core": _pairs_to_keys(selected_c),
        "selected_hinge_or_interdomain": _pairs_to_keys(selected_interdomain),
        "selected_local_support": _pairs_to_keys(selected_local),
        "selected_medium_support": _pairs_to_keys(selected_medium),
        "dynamic_pair_roles": pair_roles,
        "replica_support": support,
        "mean_frequency_by_selected_pair": mean_freq,
        "chemical_policy": payload.get("chemical_policy", "adaptive_soft_guard"),
        "chemical_score_by_selected_pair": chem,
        "dca_support": {
            "dca_score_by_selected_pair": dca,
            "dca_mean_selected": band.get("dca_mean_selected"),
            "dca_background_enrichment_ratio": band.get("dca_background_enrichment_ratio"),
            "dca_background_enrichment_pass": band.get("dca_background_enrichment_pass"),
            "dca_absolute_support_pass": band.get("dca_absolute_support_pass"),
        },
        "noise_added": band.get("noise_added"),
        "long_range_evidence_polluted": band.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": band.get("classification_coverage_ratio"),
        "claim_lock_status": payload.get("claim_lock_status") or _lookup(payload, "claim_lock_check", "status"),
        "claim_lock_failed_checks": payload.get("claim_lock_failed_checks") or _lookup(payload, "claim_lock_check", "failed_checks", default=[]),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "positive_evidence_found": positive,
        "final_status": final_status,
    }


def normalize_4ake_dynamic(bridge_payload: Optional[dict[str, Any]] = None, v14_payload: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    # V15 does not synthesize 4AKE evidence from PDB-only or visual-only files.
    # If a V15 4AKE bridge artifact exists, it may enter the panel, but its row
    # remains claim-disabled and positive only when backed by a machine-readable
    # role artifact.
    if isinstance(bridge_payload, dict) and isinstance(bridge_payload.get("protein_row"), dict):
        row = dict(bridge_payload["protein_row"])
        row["protein"] = "4AKE"
        row["claim_allowed"] = False
        row["biological_transfer_claim_allowed"] = False
        row["separation_policy"] = "dynamic_contextual_role_assignment_no_fixed_residue_cutoff"
        row["fixed_residue_cutoff_used"] = False
        row.setdefault("dynamic_pair_roles", {})
        row.setdefault("selected_pairs", [])
        row.setdefault("positive_evidence_found", False)
        return row
    return {
        "protein": "4AKE",
        "artifact_status": "missing_machine_readable_grammar_artifact",
        "target_role": "domain_hinge_object",
        "grammar_policy": "dynamic_separation_not_evaluated_artifact_missing",
        "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
        "fixed_residue_cutoff_used": False,
        "selected_pairs": [],
        "dynamic_pair_roles": {},
        "claim_lock_status": "not_evaluated_machine_readable_artifact_missing",
        "claim_lock_failed_checks": ["machine_readable_4AKE_grammar_artifact_present"],
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "positive_evidence_found": False,
        "final_status": "4AKE_not_evaluated_machine_readable_grammar_artifact_missing_claim_disabled",
    }


def _global_status(rows: list[dict[str, Any]]) -> str:
    positives = [row["protein"] for row in rows if row.get("positive_evidence_found")]
    missing = [row["protein"] for row in rows if str(row.get("artifact_status", "")).startswith("missing")]
    bridge_pending = [
        row["protein"]
        for row in rows
        if "bridge_pending" in str(row.get("artifact_status", ""))
        or "bridge_pending" in str(row.get("claim_lock_status", ""))
    ]
    if len(positives) >= 2 and "4AKE" in bridge_pending:
        return "dynamic_separation_grammar_positive_on_1UBQ_1CLL_4AKE_bridge_pending_claim_disabled"
    if len(positives) >= 2 and "4AKE" in missing:
        return "dynamic_separation_grammar_positive_on_1UBQ_1CLL_4AKE_missing_claim_disabled"
    if len(positives) == 3:
        return "dynamic_separation_grammar_coherent_across_three_object_types_claim_disabled"
    if positives:
        return "partial_dynamic_separation_grammar_signal_found_claim_disabled"
    return "clean_abstain_no_dynamic_separation_grammar_signal_claim_disabled"


def _coherence_checks(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "no_claim_allowed_anywhere": all(row.get("claim_allowed") is False for row in rows),
        "no_fixed_residue_cutoff_used_anywhere": all(row.get("fixed_residue_cutoff_used") is False for row in rows),
        "all_present_rows_have_dynamic_pair_roles": all(
            bool(row.get("dynamic_pair_roles")) for row in rows if row.get("positive_evidence_found")
        ),
        "positive_evidence_proteins": [row["protein"] for row in rows if row.get("positive_evidence_found")],
        "missing_artifacts": [row["protein"] for row in rows if str(row.get("artifact_status", "")).startswith("missing")],
        "bridge_pending_artifacts": [
            row["protein"]
            for row in rows
            if "bridge_pending" in str(row.get("artifact_status", ""))
            or "bridge_pending" in str(row.get("claim_lock_status", ""))
        ],
    }


def _table_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "protein": row.get("protein"),
        "artifact_status": row.get("artifact_status"),
        "target_role": row.get("target_role"),
        "grammar_policy": row.get("grammar_policy"),
        "separation_policy": row.get("separation_policy"),
        "fixed_residue_cutoff_used": row.get("fixed_residue_cutoff_used"),
        "selected_pairs": json.dumps(row.get("selected_pairs", []), sort_keys=True),
        "replica_support": json.dumps(row.get("replica_support", {}), sort_keys=True),
        "chemical_policy": row.get("chemical_policy"),
        "topology_policy": row.get("topology_policy"),
        "noise_added": row.get("noise_added"),
        "long_range_evidence_polluted": row.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": row.get("classification_coverage_ratio"),
        "claim_lock_status": row.get("claim_lock_status"),
        "claim_lock_failed_checks": json.dumps(row.get("claim_lock_failed_checks", []), sort_keys=True),
        "final_status": row.get("final_status"),
        "claim_allowed": row.get("claim_allowed"),
    }


def _write_report(path: Path, *, certificate: dict[str, Any]) -> None:
    rows = certificate["protein_rows"]
    lines = [
        "# V15 Dynamic Separation Grammar Readout",
        "",
        "This report is postprocess-only. It does not rerun molecular dynamics and it does not use native precision to select evidence.",
        "",
        "Sequence separation is descriptive context only. No fixed residue-count cutoff is used as a gate.",
        "",
        f"Global status: `{certificate['global_status']}`",
        f"Claim allowed: `{certificate['claim_allowed']}`",
        "",
        "## Grammar axes",
    ]
    for axis in certificate["grammar_axes"]:
        lines.append(f"- `{axis}`")
    lines.extend(["", "## Protein rows"])
    for row in rows:
        lines.extend([
            f"### {row.get('protein')}",
            f"- Artifact status: `{row.get('artifact_status')}`",
            f"- Target role: `{row.get('target_role')}`",
            f"- Grammar policy: `{row.get('grammar_policy')}`",
            f"- Separation policy: `{row.get('separation_policy')}`",
            f"- Fixed residue cutoff used: `{row.get('fixed_residue_cutoff_used')}`",
            f"- Selected pairs: `{row.get('selected_pairs', [])}`",
            f"- Claim lock: `{row.get('claim_lock_status')}`",
            f"- Final status: `{row.get('final_status')}`",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run V15 dynamic separation grammar readout.")
    parser.add_argument("--1ubq-cert", default=str(DEFAULT_UBQ_CERT))
    parser.add_argument("--1cll-cert", default=str(DEFAULT_CLL_CERT))
    parser.add_argument("--v14-cert", default=str(DEFAULT_V14_CERT))
    parser.add_argument("--4ake-bridge-cert", default=str(DEFAULT_4AKE_BRIDGE_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    ubq_payload = _read_json(Path(args.__dict__["1ubq_cert"]))
    cll_payload = _read_json(Path(args.__dict__["1cll_cert"]))
    v14_payload = _read_json(Path(args.v14_cert))
    bridge_payload = _read_json(Path(args.__dict__["4ake_bridge_cert"]))

    rows = [
        normalize_4ake_dynamic(bridge_payload, v14_payload),
        normalize_1ubq_dynamic(ubq_payload, Path(args.__dict__["1ubq_cert"])),
        normalize_1cll_dynamic(cll_payload, Path(args.__dict__["1cll_cert"])),
    ]
    checks = _coherence_checks(rows)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    certificate = {
        "kind": "V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
        "fixed_residue_cutoff_used": False,
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "global_status": _global_status(rows),
        "grammar_axes": GRAMMAR_AXES,
        "coherence_checks": checks,
        "positive_evidence_proteins": checks["positive_evidence_proteins"],
        "missing_artifacts": checks["missing_artifacts"],
        "protein_rows": rows,
        "interpretation": (
            "V15 retires fixed sequence-distance cutoffs as evidence gates. Pair separation is reported as dynamic "
            "context, while role validity is decided by purpose/topology, DCA, geometry, replica support, chemical policy, "
            "noise/pollution guards, and claim locks."
        ),
    }

    cert_path = out_dir / "v15_dynamic_separation_grammar_readout_certificate.json"
    table_path = out_dir / "v15_dynamic_separation_grammar_table.csv"
    report_path = out_dir / "V15_DYNAMIC_SEPARATION_GRAMMAR_REPORT.md"
    cert_path.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(table_path, [_table_row(row) for row in rows])
    _write_report(report_path, certificate=certificate)

    print(json.dumps({
        "kind": certificate["kind"],
        "certificate": str(cert_path),
        "table": str(table_path),
        "report": str(report_path),
        "global_status": certificate["global_status"],
        "positive_evidence_proteins": certificate["positive_evidence_proteins"],
        "missing_artifacts": certificate["missing_artifacts"],
        "fixed_residue_cutoff_used": certificate["fixed_residue_cutoff_used"],
        "claim_allowed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
