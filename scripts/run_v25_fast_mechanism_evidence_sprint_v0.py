#!/usr/bin/env python3
from __future__ import annotations

"""V25 fast mechanism evidence sprint.

One-shot mechanism sprint after V24. This is not a new MD run and not a
folding-claim upgrade. It consolidates three practical things in one batch:

1. KcsA coupling/interface evidence readout.
2. XCL1 state-specific evidence readout preservation.
3. A unified all-target mechanism-operator table across the current panel.
4. One next mechanism-test decision for V26.

The purpose is to move from target-by-target debugging toward repeated
mechanism operators: role_detect, context_detect, evidence_assign,
reachability_check, replica_support, chemical_confidence, topology_context,
state_separation, interface_context, clean_abstain, and pollution_guard.
"""

import argparse
import csv
import glob
import json
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V15_CERT = RUN_ROOT / "V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT" / "v15_dynamic_separation_grammar_readout_certificate.json"
DEFAULT_V23_CERT = RUN_ROOT / "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST" / "v23_p53_tad_mdm2_external_evidence_contrast_certificate.json"
DEFAULT_V24_CERT = RUN_ROOT / "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION" / "v24_kcsa_external_annotation_and_assembly_acquisition_certificate.json"
DEFAULT_V20_CERT = RUN_ROOT / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT" / "v20_xcl1_state_specific_evidence_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V25_FAST_MECHANISM_EVIDENCE_SPRINT"

KCSA_COUPLING_PATTERNS = [
    "external_msa/**/*kcsa*",
    "external_msa/**/*KcsA*",
    "external_msa/**/*potassium*channel*",
    "external_msa_free_predictors/**/*kcsa*",
    "external_msa_free_predictors/**/*KcsA*",
    "data/**/*kcsa*coupl*",
    "data/**/*KcsA*coupl*",
    "data/**/*potassium*channel*coupl*",
]

OPERATOR_ORDER = [
    "role_detect",
    "context_detect",
    "evidence_assign",
    "reachability_check",
    "replica_support",
    "chemical_confidence",
    "topology_context",
    "state_separation",
    "interface_context",
    "clean_abstain",
    "pollution_guard",
]


def _read_json(path: Path, label: str, required: bool = True) -> dict[str, Any]:
    if not path.exists():
        if required:
            raise SystemExit(f"missing {label}: {path}")
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _glob_repo(patterns: list[str], limit: int = 80) -> list[str]:
    hits: list[str] = []
    for pattern in patterns:
        for hit in glob.glob(str(REPO_ROOT / pattern), recursive=True):
            p = Path(hit)
            if p.is_file():
                try:
                    hits.append(str(p.relative_to(REPO_ROOT)))
                except ValueError:
                    hits.append(str(p))
    return sorted(set(hits))[:limit]


def _rows_by_target(cert: dict[str, Any]) -> dict[str, dict[str, Any]]:
    for key in ("protein_rows", "target_rows", "target_results", "targets"):
        rows = cert.get(key)
        if isinstance(rows, list):
            result: dict[str, dict[str, Any]] = {}
            for row in rows:
                if not isinstance(row, dict):
                    continue
                tid = row.get("target_id") or row.get("protein") or row.get("source_accession") or row.get("target")
                if tid:
                    result[str(tid)] = row
            if result:
                return result
    # V15 sometimes exposes a dict-like table under protein_rows.
    rows = cert.get("protein_rows")
    if isinstance(rows, dict):
        return {str(k): v for k, v in rows.items() if isinstance(v, dict)}
    return {}


def _find_target_row(cert: dict[str, Any], names: list[str]) -> dict[str, Any]:
    rows = _rows_by_target(cert)
    lower = {k.lower(): v for k, v in rows.items()}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    # fallback by substring
    for key, val in rows.items():
        lk = key.lower()
        if any(name.lower() in lk for name in names):
            return val
    return {}


def _selected_pairs(row: dict[str, Any], cert: dict[str, Any] | None = None) -> list[str]:
    for obj in (row, cert or {}):
        val = obj.get("selected_pairs") if isinstance(obj, dict) else None
        if isinstance(val, list):
            return [str(x) for x in val]
    return []


def _pair_map(row: dict[str, Any]) -> dict[str, Any]:
    for key in ("dynamic_pair_roles", "candidate_pair_roles", "pair_roles"):
        val = row.get(key)
        if isinstance(val, dict):
            return val
    return {}


def build_kcsa_coupling_interface_readout(v24: dict[str, Any]) -> dict[str, Any]:
    coupling_files = sorted(set(v24.get("external_coupling_or_msa_files") or []))
    if not coupling_files:
        coupling_files = _glob_repo(KCSA_COUPLING_PATTERNS)

    coupling_available = bool(v24.get("external_couplings_present") is True or coupling_files)
    interface_present = bool(v24.get("oligomer_or_chain_interface_context_present") is True)
    pore_filter_present = bool(v24.get("pore_filter_annotation_present") is True)
    pore_filter_coupling_support = bool(coupling_available and pore_filter_present)

    positive_interface = interface_present and pore_filter_present and bool(v24.get("soluble_core_misclassification_avoided") is True)
    status = (
        "V25A_KcsA_COUPLING_INTERFACE_EVIDENCE_FOUND_CLAIM_DISABLED"
        if coupling_available and positive_interface
        else "V25A_KcsA_INTERFACE_CONTEXT_FOUND_COUPLING_MISSING_CLAIM_DISABLED"
        if positive_interface
        else "V25A_KcsA_CLEAN_ABSTAIN_INTERFACE_OR_PORE_CONTEXT_INSUFFICIENT_CLAIM_DISABLED"
    )

    return {
        "kind": "V25A_KcsA_COUPLING_INTERFACE_READOUT_v0",
        "target_id": "KcsA",
        "run_mode": "zero_md_coupling_interface_readout_no_simulation_no_threshold_tuning",
        "readout_status": status,
        "role_class": "membrane_pore_oligomer_object",
        "positive_pressure_evidence_found": positive_interface,
        "positive_folding_evidence_found": False,
        "kcsa_coupling_available": coupling_available,
        "external_coupling_or_msa_files": coupling_files,
        "interface_evidence_present": interface_present,
        "pore_filter_coupling_support": pore_filter_coupling_support,
        "pore_filter_role_preserved": pore_filter_present,
        "transmembrane_scaffold_preserved": bool(v24.get("transmembrane_role_present") is True),
        "membrane_context_preserved": bool(v24.get("membrane_environment_context_present") is True),
        "assembly_context_preserved": bool(v24.get("biological_assembly_context_locked") is True),
        "tetramer_claim_allowed": False,
        "whole_channel_fold_claim_allowed": False,
        "tetramer_claim_made": False,
        "whole_channel_fold_claim_made": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selection_threshold_used": False,
        "coupling_policy": "do_not_synthesize_couplings_missing_is_reported_not_backfilled",
        "available_evidence": [
            item for item, present in {
                "pore_filter_role": pore_filter_present,
                "TM_scaffold_context": v24.get("transmembrane_role_present") is True,
                "interface_context": interface_present,
                "assembly_context_reference": v24.get("biological_assembly_context_locked") is True,
                "soluble_core_guard": v24.get("soluble_core_misclassification_avoided") is True,
                "external_coupling_files": coupling_available,
            }.items() if present
        ],
        "missing_evidence": [
            item for item, present in {
                "external_couplings_if_available": coupling_available,
                "pore_filter_coupling_support": pore_filter_coupling_support,
            }.items() if not present
        ],
        "source_v24_status": v24.get("test_status"),
    }


def build_xcl1_state_specific_readout(v20: dict[str, Any]) -> dict[str, Any]:
    state_a = bool(v20.get("state_A_role_evidence_found") is True)
    state_b = bool(v20.get("state_B_role_evidence_found") is True)
    mixed = bool(v20.get("mixed_state_pollution") is True)
    single = bool(v20.get("single_fold_claim_made") is True)
    state_specific = state_a and state_b and not mixed and not single
    status = (
        "V25B_XCL1_STATE_SPECIFIC_EVIDENCE_PRESERVED_CLAIM_DISABLED"
        if state_specific
        else "V25B_XCL1_CLEAN_ABSTAIN_STATE_SPECIFIC_EVIDENCE_INSUFFICIENT_CLAIM_DISABLED"
    )
    return {
        "kind": "V25B_XCL1_STATE_SPECIFIC_READOUT_v0",
        "target_id": "XCL1_lymphotactin",
        "run_mode": "zero_md_state_separation_preservation_no_simulation_no_threshold_tuning",
        "readout_status": status,
        "role_class": "metamorphic_switch_object",
        "positive_pressure_evidence_found": state_specific,
        "positive_folding_evidence_found": False,
        "state_A_detected": state_a,
        "state_B_detected": state_b,
        "state_specific_role_evidence_found": bool(v20.get("state_specific_role_evidence_found") is True and state_specific),
        "mixed_state_pollution": mixed,
        "single_fold_forcing": single,
        "single_fold_claim_made": single,
        "fold_switch_claim_made": bool(v20.get("fold_switch_claim_made") is True),
        "mixed_state_contact_pooling_used": bool(v20.get("mixed_state_contact_pooling_used") is True),
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selection_threshold_used": False,
        "state_separation_policy": "state_specific_buckets_no_cross_state_pooling_no_single_native_assumption",
        "available_evidence": v20.get("available_evidence", []),
        "missing_evidence": v20.get("missing_evidence", []),
        "forbidden_misclassification_violations": v20.get("forbidden_misclassification_violations", []),
        "source_v20_status": v20.get("test_status"),
    }


def _operator_row(
    target: str,
    pressure_type: str,
    role_detected: bool,
    evidence_type: str,
    selected_signal: str,
    abstain_reason: str | None,
    missing_next: list[str],
    operators: list[str],
    claim_allowed: bool = False,
) -> dict[str, Any]:
    return {
        "target": target,
        "pressure_type": pressure_type,
        "role_detected": role_detected,
        "evidence_type": evidence_type,
        "selected_signal": selected_signal,
        "abstain_reason_if_any": abstain_reason,
        "missing_next_evidence": missing_next,
        "mechanism_operator_used": operators,
        "claim_allowed": claim_allowed,
    }


def build_operator_table(v15: dict[str, Any], v23: dict[str, Any], v24: dict[str, Any], v20: dict[str, Any], kcsa: dict[str, Any], xcl1: dict[str, Any]) -> dict[str, Any]:
    v15_rows = _rows_by_target(v15)
    # Fallback to names as printed in current V15 table.
    row_4ake = _find_target_row(v15, ["4AKE", "4AKE:A"])
    row_1ubq = _find_target_row(v15, ["1UBQ", "1UBQ:A"])
    row_1cll = _find_target_row(v15, ["1CLL", "1CLL:A"])

    rows: list[dict[str, Any]] = []
    rows.append(_operator_row(
        "4AKE",
        "domain_hinge_closure_object",
        bool(row_4ake.get("artifact_status") == "present_machine_readable_4ake_role_artifact" or _selected_pairs(row_4ake)),
        "dynamic_domain_core_role_evidence_from_DCA_trajectory_frequency",
        ",".join(_selected_pairs(row_4ake)) or "124-135_or_prior_4AKE_bridge_signal",
        None if _selected_pairs(row_4ake) else "no_loaded_V15_4AKE_selected_pair_in_current_certificate",
        ["interdomain_hinge_or_closure_evidence", "external_validation_if_available"],
        ["role_detect", "context_detect", "evidence_assign", "reachability_check", "replica_support", "chemical_confidence", "topology_context", "pollution_guard"],
    ))
    rows.append(_operator_row(
        "1UBQ",
        "single_domain_compact_object",
        bool(row_1ubq.get("artifact_status") == "present" or _selected_pairs(row_1ubq)),
        "compact_core_role_evidence_from_dynamic_V15_grammar",
        ",".join(_selected_pairs(row_1ubq)) or "23-48_or_prior_1UBQ_signal",
        None if _selected_pairs(row_1ubq) else "no_loaded_V15_1UBQ_selected_pair_in_current_certificate",
        ["cross_target_validation", "DCA_background_enrichment"],
        ["role_detect", "evidence_assign", "replica_support", "chemical_confidence", "pollution_guard"],
    ))
    rows.append(_operator_row(
        "1CLL",
        "multi_domain_composite_object",
        bool(row_1cll.get("artifact_status") == "present" or _selected_pairs(row_1cll)),
        "hierarchical_domain_core_role_evidence",
        ",".join(_selected_pairs(row_1cll)) or "97-133_or_prior_1CLL_signal",
        "interdomain_hinge_not_yet_proven",
        ["interdomain_hinge_evidence_if_claim_needed"],
        ["role_detect", "context_detect", "evidence_assign", "replica_support", "chemical_confidence", "topology_context", "pollution_guard"],
    ))
    rows.append(_operator_row(
        "p53_TAD_MDM2",
        "disorder_partner_induced_binding_pressure",
        bool(v23.get("positive_pressure_evidence_found") is True),
        "external_disorder_reference_plus_partner_bound_interface_evidence",
        "isolated_disorder_context + bound_1YCR_interface_helix_proxy",
        None,
        v23.get("missing_evidence", []),
        ["role_detect", "context_detect", "evidence_assign", "interface_context", "clean_abstain", "pollution_guard"],
    ))
    rows.append(_operator_row(
        "KcsA",
        "membrane_pore_oligomer_environment_pressure",
        bool(v24.get("positive_pressure_evidence_found") is True),
        "pore_filter_TM_ion_shell_interface_annotation_evidence",
        "TVGYG + K_plus_shell + TM_scaffold + chain_interface_context",
        None if kcsa.get("positive_pressure_evidence_found") else "coupling_or_interface_context_insufficient",
        sorted(set(v24.get("missing_evidence", []) + kcsa.get("missing_evidence", []))),
        ["role_detect", "context_detect", "evidence_assign", "topology_context", "interface_context", "clean_abstain", "pollution_guard"],
    ))
    rows.append(_operator_row(
        "XCL1_lymphotactin",
        "metamorphic_fold_switching_pressure",
        bool(v20.get("positive_pressure_evidence_found") is True),
        "state_specific_A_B_role_separation_evidence",
        "state_A_chemokine_context + state_B_dimer_context_kept_separate",
        None if xcl1.get("positive_pressure_evidence_found") else "state_specific_evidence_insufficient",
        v20.get("missing_evidence", []),
        ["role_detect", "context_detect", "evidence_assign", "state_separation", "interface_context", "clean_abstain", "pollution_guard"],
    ))

    counter: Counter[str] = Counter()
    for row in rows:
        counter.update(row["mechanism_operator_used"])
    return {
        "kind": "V25C_UNIFIED_MECHANISM_OPERATOR_TABLE_v0",
        "operator_policy": "same_operator_alphabet_across_targets_not_same_threshold_not_same_fold_claim",
        "operator_alphabet": OPERATOR_ORDER,
        "operator_counts": {op: counter.get(op, 0) for op in OPERATOR_ORDER},
        "targets": rows,
        "claim_allowed": False,
        "positive_folding_evidence_targets": [],
    }


def build_next_decision(kcsa: dict[str, Any], xcl1: dict[str, Any], v23: dict[str, Any]) -> dict[str, Any]:
    k_has_coupling = bool(kcsa.get("kcsa_coupling_available") is True)
    x_state_clean = bool(
        xcl1.get("state_specific_role_evidence_found") is True
        and xcl1.get("mixed_state_pollution") is False
        and xcl1.get("single_fold_claim_made") is False
    )
    p53_clean = bool(
        v23.get("positive_external_annotation_evidence_found") is True
        and v23.get("isolated_TAD_clean_abstain_preserved") is True
        and v23.get("partner_bound_interface_evidence_preserved") is True
    )

    if k_has_coupling and kcsa.get("interface_evidence_present") is True:
        target = "KcsA"
        test = "KcsA_COUPLING_INTERFACE_OPERATOR_TEST"
        reason = "KcsA has coupling/interface support available; membrane/pore pressure operator can be tested without membrane MD"
        new_md_allowed = False
    elif x_state_clean:
        target = "XCL1_lymphotactin"
        test = "XCL1_STATE_SEPARATION_OPERATOR_TEST"
        reason = "XCL1 has clean state-specific separation with no mixed-state pollution; strongest conceptual mechanism target, but still no MD until condition/coupling evidence is available"
        new_md_allowed = False
    elif p53_clean:
        target = "p53_TAD_MDM2"
        test = "p53_CONTEXT_DEPENDENT_ABSTAIN_INTERFACE_OPERATOR_TEST"
        reason = "p53 has the strongest ready external evidence contrast and lowest false-win risk"
        new_md_allowed = False
    else:
        target = None
        test = None
        reason = "no target has sufficient mechanism-readiness without adding external evidence first"
        new_md_allowed = False

    return {
        "kind": "V25_NEXT_MECHANISM_TEST_DECISION_v0",
        "decision_status": "V26_TARGET_SELECTED" if target else "V26_BLOCKED_PENDING_EVIDENCE",
        "selected_V26_target": target,
        "selected_V26_test": test,
        "reason": reason,
        "evidence_readiness_basis": {
            "KcsA_coupling_available": k_has_coupling,
            "KcsA_interface_evidence_present": bool(kcsa.get("interface_evidence_present") is True),
            "XCL1_state_specific_clean": x_state_clean,
            "p53_external_contrast_clean": p53_clean,
        },
        "new_MD_allowed": new_md_allowed,
        "new_MD_recommended": False,
        "new_MD_allowed_policy": "only_if_external_evidence_ready_and_false_win_risk_low",
        "parallel_MD_paths_allowed": False,
        "claim_allowed": False,
        "positive_folding_evidence_targets": [],
    }


def build_v25(v15: dict[str, Any], v23: dict[str, Any], v24: dict[str, Any], v20: dict[str, Any]) -> dict[str, Any]:
    kcsa = build_kcsa_coupling_interface_readout(v24)
    xcl1 = build_xcl1_state_specific_readout(v20)
    operator_table = build_operator_table(v15, v23, v24, v20, kcsa, xcl1)
    decision = build_next_decision(kcsa, xcl1, v23)

    failed: list[str] = []
    for label, cert, status_key, expected in [
        ("V23_p53_external_contrast", v23, "test_status", "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"),
        ("V24_KcsA_annotation", v24, "test_status", "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED"),
        ("V20_XCL1_state_specific", v20, "test_status", "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"),
    ]:
        if cert.get(status_key) != expected:
            failed.append(f"{label}_passed")
    for label, cert in [("V23", v23), ("V24", v24), ("V20", v20)]:
        if cert.get("new_md_executed") is not False:
            failed.append(f"{label}_no_new_md")
        if cert.get("positive_folding_evidence_found") is not False:
            failed.append(f"{label}_no_folding_evidence_claim")
        if cert.get("claim_allowed") is not False:
            failed.append(f"{label}_claim_disabled")

    status = "V25_FAST_MECHANISM_EVIDENCE_SPRINT_LOCKED" if not failed else "V25_FAST_MECHANISM_EVIDENCE_SPRINT_BLOCKED"
    pressure_targets = [
        t for t, ok in {
            "p53_TAD_MDM2": v23.get("positive_pressure_evidence_found") is True,
            "KcsA": v24.get("positive_pressure_evidence_found") is True,
            "XCL1_lymphotactin": v20.get("positive_pressure_evidence_found") is True,
        }.items() if ok
    ]
    return {
        "kind": "V25_FAST_MECHANISM_EVIDENCE_SPRINT_v0",
        "run_mode": "fast_mechanism_sprint_no_simulation_no_threshold_tuning",
        "sprint_status": status,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_pressure_evidence_targets": pressure_targets,
        "positive_folding_evidence_targets": [],
        "fast_sprint_failed_checks": failed,
        "kcsa_coupling_interface_readout": kcsa,
        "xcl1_state_specific_readout": xcl1,
        "unified_mechanism_operator_table": operator_table,
        "next_mechanism_test_decision": decision,
        "locked_interpretation": (
            "V25 shifts from target-by-target debugging to mechanism sprint mode. It preserves KcsA coupling/interface readiness, "
            "XCL1 state-specific separation, and an all-target mechanism-operator table across 4AKE, 1UBQ, 1CLL, p53/MDM2, KcsA, and XCL1. "
            "This is evidence for a repeated role/context/operator grammar, not a universal folding claim."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    d = cert.get("next_mechanism_test_decision", {})
    lines = [
        "# V25 Fast Mechanism Evidence Sprint",
        "",
        f"Status: `{cert.get('sprint_status')}`",
        f"Positive pressure evidence targets: `{cert.get('positive_pressure_evidence_targets')}`",
        f"Positive folding evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## KcsA coupling/interface",
        f"`{cert.get('kcsa_coupling_interface_readout')}`",
        "",
        "## XCL1 state-specific",
        f"`{cert.get('xcl1_state_specific_readout')}`",
        "",
        "## Next V26 decision",
        f"Selected target: `{d.get('selected_V26_target')}`",
        f"Selected test: `{d.get('selected_V26_test')}`",
        f"Reason: {d.get('reason')}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_operator_csv(path: Path, table: dict[str, Any]) -> None:
    rows = table.get("targets") if isinstance(table.get("targets"), list) else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "target", "pressure_type", "role_detected", "evidence_type", "selected_signal",
            "abstain_reason_if_any", "missing_next_evidence", "mechanism_operator_used", "claim_allowed",
        ])
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["missing_next_evidence"] = ";".join(str(x) for x in row.get("missing_next_evidence", []))
            out["mechanism_operator_used"] = ";".join(str(x) for x in row.get("mechanism_operator_used", []))
            writer.writerow(out)


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v25_fast_mechanism_evidence_sprint_certificate.json",
        "kcsa": out_dir / "v25_kcsa_coupling_interface_readout.json",
        "xcl1": out_dir / "v25_xcl1_state_specific_readout.json",
        "operator_table": out_dir / "v25_unified_mechanism_operator_table.json",
        "operator_table_csv": out_dir / "v25_unified_mechanism_operator_table.csv",
        "decision": out_dir / "v25_next_mechanism_test_decision.json",
        "report": out_dir / "V25_FAST_MECHANISM_EVIDENCE_SPRINT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["kcsa"].write_text(json.dumps(cert["kcsa_coupling_interface_readout"], indent=2, sort_keys=True), encoding="utf-8")
    paths["xcl1"].write_text(json.dumps(cert["xcl1_state_specific_readout"], indent=2, sort_keys=True), encoding="utf-8")
    paths["operator_table"].write_text(json.dumps(cert["unified_mechanism_operator_table"], indent=2, sort_keys=True), encoding="utf-8")
    _write_operator_csv(paths["operator_table_csv"], cert["unified_mechanism_operator_table"])
    paths["decision"].write_text(json.dumps(cert["next_mechanism_test_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v15-cert", type=Path, default=DEFAULT_V15_CERT)
    parser.add_argument("--v23-cert", type=Path, default=DEFAULT_V23_CERT)
    parser.add_argument("--v24-cert", type=Path, default=DEFAULT_V24_CERT)
    parser.add_argument("--v20-cert", type=Path, default=DEFAULT_V20_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v15 = _read_json(args.v15_cert, "V15 dynamic grammar certificate", required=False)
    v23 = _read_json(args.v23_cert, "V23 p53 external contrast certificate")
    v24 = _read_json(args.v24_cert, "V24 KcsA external annotation certificate")
    v20 = _read_json(args.v20_cert, "V20 XCL1 state-specific certificate")
    cert = build_v25(v15, v23, v24, v20)
    paths = write_outputs(args.out_dir, cert)
    decision = cert["next_mechanism_test_decision"]
    print(json.dumps({
        "kind": cert["kind"],
        "sprint_status": cert["sprint_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "positive_pressure_evidence_targets": cert["positive_pressure_evidence_targets"],
        "positive_folding_evidence_targets": cert["positive_folding_evidence_targets"],
        "selected_V26_target": decision.get("selected_V26_target"),
        "selected_V26_test": decision.get("selected_V26_test"),
        "new_MD_allowed": decision.get("new_MD_allowed"),
        "certificate": str(paths["certificate"]),
        "kcsa_coupling_interface_readout": str(paths["kcsa"]),
        "xcl1_state_specific_readout": str(paths["xcl1"]),
        "unified_mechanism_operator_table": str(paths["operator_table"]),
        "next_mechanism_test_decision": str(paths["decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
