#!/usr/bin/env python3
from __future__ import annotations

"""V34 KcsA discriminative content controls.

Purpose
-------
V33 proved that real external source files could open a claim-disabled
operator readout for KcsA. That is necessary, but not enough: a weak readout
could still pass by merely counting files/rows.

V34 asks a sharper question:
    Does the imported KcsA evidence contain the actual pore/filter and
    tetramer/interface content expected for KcsA, and do adversarially damaged
    versions fail?

Boundary
--------
- No MD.
- No folding claim.
- No native-metric selection.
- No threshold tuning for prediction.
- Coordinate-derived RCSB source is used only as external source-import
  content for controls; this is not de-novo folding prediction.
"""

import argparse
import copy
import csv
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V33_CERT = RUN_ROOT / "V33_CONSTRAINT_BACKED_OPERATOR_READOUT" / "v33_constraint_backed_operator_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS"

PASSED = "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED"
FAILED = "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_FAILED_CLAIM_DISABLED"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _norm_rel(path: str) -> str:
    return path.replace("\\", "/").strip()


def _load_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing imported source CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [dict(row) for row in rows]


def _int_or_none(value: Any) -> int | None:
    try:
        if value in (None, ""):
            return None
        return int(float(str(value).strip()))
    except Exception:
        return None


def _upper_values(rows: list[dict[str, str]], key: str) -> set[str]:
    return {str(row.get(key, "")).strip().upper() for row in rows if str(row.get(key, "")).strip()}


def _residue_numbers(rows: list[dict[str, str]]) -> set[int]:
    vals: set[int] = set()
    for row in rows:
        for key in ("residue_number", "residue_i", "residue_j"):
            val = _int_or_none(row.get(key))
            if val is not None:
                vals.add(val)
    return vals


def _chains(rows: list[dict[str, str]]) -> set[str]:
    vals: set[str] = set()
    for row in rows:
        for key in ("chain", "chain_a", "chain_b"):
            val = str(row.get(key, "")).strip()
            if val:
                vals.add(val)
    return vals


def _interface_pairs(rows: list[dict[str, str]]) -> set[str]:
    pairs: set[str] = set()
    for row in rows:
        a = str(row.get("chain_a", "")).strip()
        b = str(row.get("chain_b", "")).strip()
        if not a or not b or a == b:
            continue
        left, right = sorted([a, b])
        pairs.add(f"{left}-{right}")
    return pairs


def validate_kcsa_content(pore_rows: list[dict[str, str]], interface_rows: list[dict[str, str]]) -> dict[str, Any]:
    failures: list[str] = []

    pore_residues = _residue_numbers(pore_rows)
    pore_chains = _chains(pore_rows)
    pore_motifs = _upper_values(pore_rows, "filter_motif")
    pore_ions = _upper_values(pore_rows, "ion_name") | _upper_values(pore_rows, "ion_element")
    pore_classes = _upper_values(pore_rows, "constraint_class")

    iface_chains = _chains(interface_rows)
    iface_pairs = _interface_pairs(interface_rows)
    iface_classes = _upper_values(interface_rows, "constraint_class")

    if len(pore_rows) < 4:
        failures.append("pore_filter_too_few_rows")
    if len(interface_rows) < 4:
        failures.append("assembly_interface_too_few_rows")

    # KcsA 1BL8 import should carry concrete filter identity, not just a generic file.
    # The builder produced rows around canonical filter residues 75-78 that contact K+.
    if not {75, 76, 77, 78}.issubset(pore_residues):
        failures.append("missing_canonical_filter_residues_75_78")
    if "TVGYG" not in pore_motifs:
        failures.append("missing_TVGYG_motif_label")
    if not ({"K", "K+", "POT"} & pore_ions):
        failures.append("missing_potassium_ion_identity")
    if "PORE_FILTER_POTASSIUM_COORDINATION_CONTACT" not in pore_classes:
        failures.append("missing_pore_filter_potassium_coordination_class")
    if len(pore_chains) < 4:
        failures.append("pore_filter_not_observed_across_four_chains")

    # Tetramer/interface evidence should contain inter-chain contacts across a 4-chain assembly.
    if "ASSEMBLY_INTERFACE_HEAVY_ATOM_CONTACT" not in iface_classes:
        failures.append("missing_assembly_interface_contact_class")
    if len(iface_chains) < 4:
        failures.append("insufficient_oligomer_chains")
    if len(iface_pairs) < 4:
        failures.append("insufficient_distinct_chain_pairs")

    return {
        "valid": not failures,
        "failures": sorted(set(failures)),
        "summary": {
            "pore_row_count": len(pore_rows),
            "interface_row_count": len(interface_rows),
            "pore_residue_numbers": sorted(pore_residues),
            "pore_chains": sorted(pore_chains),
            "pore_motif_labels": sorted(pore_motifs),
            "pore_ion_names": sorted(pore_ions),
            "pore_constraint_classes": sorted(pore_classes),
            "interface_chains": sorted(iface_chains),
            "interface_chain_pair_count": len(iface_pairs),
            "interface_chain_pairs": sorted(iface_pairs),
            "interface_constraint_classes": sorted(iface_classes),
        },
    }


def _source_paths_from_v33(v33: dict[str, Any]) -> tuple[list[str], list[str]]:
    rows = [row for row in _as_list(v33.get("selected_target_import_rows")) if isinstance(row, dict)]
    pore_paths: list[str] = []
    iface_paths: list[str] = []
    for row in rows:
        rel = _norm_rel(str(row.get("file_path", "")))
        if not rel:
            continue
        if row.get("evidence_type") == "pore_filter_coupling":
            pore_paths.append(rel)
        elif row.get("evidence_type") == "assembly_interface_constraint":
            iface_paths.append(rel)
    return pore_paths, iface_paths


def _load_rows(paths: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rel in paths:
        rows.extend(_load_csv(REPO_ROOT / rel))
    return rows


def _control(control_id: str, observed: dict[str, Any], expected_valid: bool, expected_failure_any: list[str], reason: str) -> dict[str, Any]:
    failures = set(_as_list(observed.get("failures")))
    expected_set = set(expected_failure_any)
    pass_failures = True if expected_valid else bool(expected_set & failures)
    passed = observed.get("valid") is expected_valid and pass_failures
    return {
        "control_id": control_id,
        "expected_valid": expected_valid,
        "observed_valid": observed.get("valid"),
        "expected_failure_any": expected_failure_any,
        "observed_failures": observed.get("failures"),
        "passed": passed,
        "reason": reason,
    }


def build_v34(v33: dict[str, Any]) -> dict[str, Any]:
    precondition_failures: list[str] = []
    if v33.get("readout_status") != "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED":
        precondition_failures.append("V33_readout_not_passed")
    if v33.get("selected_V33_target") != "KcsA":
        precondition_failures.append("V33_selected_target_not_KcsA")
    for key in ["claim_allowed", "new_MD_allowed", "positive_folding_evidence_found", "folding_problem_solved"]:
        if v33.get(key) is not False:
            precondition_failures.append(f"V33_{key}_must_be_false")

    pore_paths, iface_paths = _source_paths_from_v33(v33)
    if not pore_paths:
        precondition_failures.append("V33_missing_pore_source_paths")
    if not iface_paths:
        precondition_failures.append("V33_missing_interface_source_paths")

    if precondition_failures:
        return {
            "kind": "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_v0",
            "run_mode": "content_discriminative_controls_no_MD_no_threshold_tuning_claim_disabled",
            "control_status": FAILED,
            "precondition_failures": sorted(set(precondition_failures)),
            "controls": [],
            "control_count": 0,
            "passed_control_count": 0,
            "claim_allowed": False,
            "new_MD_allowed": False,
            "positive_folding_evidence_found": False,
            "folding_problem_solved": False,
            "next_action": "fix_V33_preconditions_before_discriminative_controls",
        }

    pore_rows = _load_rows(pore_paths)
    iface_rows = _load_rows(iface_paths)

    baseline = validate_kcsa_content(pore_rows, iface_rows)
    controls = [
        _control(
            "baseline_real_kcsa_content_signature_valid",
            baseline,
            True,
            [],
            "Real imported KcsA source rows must carry pore/filter identity, K+ coordination identity, and tetramer/interface identity.",
        )
    ]

    shuffled_residue_pore = copy.deepcopy(pore_rows)
    for idx, row in enumerate(shuffled_residue_pore):
        if "residue_number" in row:
            row["residue_number"] = str(900 + idx)
        if "residue_i" in row:
            row["residue_i"] = str(900 + idx)
        if "residue_j" in row:
            row["residue_j"] = str(950 + idx)
    controls.append(_control(
        "pore_residue_shuffle_breaks_canonical_filter_identity",
        validate_kcsa_content(shuffled_residue_pore, iface_rows),
        False,
        ["missing_canonical_filter_residues_75_78"],
        "A residue-number randomized pore file must not still count as the KcsA canonical filter.",
    ))

    fake_motif_pore = copy.deepcopy(pore_rows)
    for row in fake_motif_pore:
        if "filter_motif" in row:
            row["filter_motif"] = "FAKE"
    controls.append(_control(
        "pore_motif_relabel_breaks_TVGYG_identity",
        validate_kcsa_content(fake_motif_pore, iface_rows),
        False,
        ["missing_TVGYG_motif_label"],
        "A pore file without the TVGYG label must not pass as KcsA filter grammar.",
    ))

    fake_ion_pore = copy.deepcopy(pore_rows)
    for row in fake_ion_pore:
        if "ion_name" in row:
            row["ion_name"] = "NA"
        if "ion_element" in row:
            row["ion_element"] = "NA"
    controls.append(_control(
        "pore_ion_relabel_breaks_potassium_identity",
        validate_kcsa_content(fake_ion_pore, iface_rows),
        False,
        ["missing_potassium_ion_identity"],
        "KcsA filter evidence must not survive if K+ identity is removed or relabeled.",
    ))

    fake_class_pore = copy.deepcopy(pore_rows)
    for row in fake_class_pore:
        if "constraint_class" in row:
            row["constraint_class"] = "generic_contact"
    controls.append(_control(
        "pore_constraint_class_relabel_breaks_filter_contact_identity",
        validate_kcsa_content(fake_class_pore, iface_rows),
        False,
        ["missing_pore_filter_potassium_coordination_class"],
        "Generic contacts must not be promoted into pore/filter potassium-coordination evidence.",
    ))

    same_chain_iface = copy.deepcopy(iface_rows)
    for row in same_chain_iface:
        if "chain_a" in row:
            row["chain_a"] = "A"
        if "chain_b" in row:
            row["chain_b"] = "A"
    controls.append(_control(
        "interface_same_chain_collapse_breaks_tetramer_context",
        validate_kcsa_content(pore_rows, same_chain_iface),
        False,
        ["insufficient_oligomer_chains", "insufficient_distinct_chain_pairs"],
        "A collapsed same-chain interface file must not pass as tetramer/interface context.",
    ))

    singleton_pore = copy.deepcopy(pore_rows[:1])
    singleton_iface = copy.deepcopy(iface_rows[:1])
    controls.append(_control(
        "singleton_rows_do_not_define_operator_context",
        validate_kcsa_content(singleton_pore, singleton_iface),
        False,
        ["pore_filter_too_few_rows", "assembly_interface_too_few_rows"],
        "One row is not enough to define a stable operator context for either pore or interface.",
    ))

    swapped = validate_kcsa_content(copy.deepcopy(iface_rows), copy.deepcopy(pore_rows))
    controls.append(_control(
        "swapped_pore_and_interface_schema_fails_content_validation",
        swapped,
        False,
        ["missing_TVGYG_motif_label", "missing_assembly_interface_contact_class", "insufficient_distinct_chain_pairs"],
        "Swapping pore/interface schemas must not pass merely because rows exist.",
    ))

    all_passed = all(row.get("passed") is True for row in controls)
    status = PASSED if all_passed else FAILED
    return {
        "kind": "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_v0",
        "run_mode": "content_discriminative_controls_no_MD_no_threshold_tuning_claim_disabled",
        "control_status": status,
        "precondition_failures": [],
        "source_boundary": "external_coordinate_derived_content_controls_only_not_de_novo_folding_prediction",
        "source_pore_paths": pore_paths,
        "source_interface_paths": iface_paths,
        "baseline_content_signature": baseline.get("summary"),
        "controls": controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row.get("passed") is True),
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "fixed_threshold_policy": "content_invariant_controls_only_no_prediction_threshold_tuning",
        "native_metrics_used_for_selection": False,
        "next_action": "move_to_non_coordinate_external_evolutionary_or_blind_holdout_tests" if status == PASSED else "fix_discriminative_content_controls_before_any_claim_or_MD",
        "locked_interpretation": (
            "Passing V34 means the imported KcsA source files are not merely counted: damaged pore/filter identity, K+ identity, constraint-class identity, tetramer/interface identity, singleton rows, and swapped schemas fail. This is still not de-novo folding prediction and does not solve universal protein folding."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V34 KcsA Discriminative Content Controls",
        "",
        f"Status: `{cert.get('control_status')}`",
        f"Passed controls: `{cert.get('passed_control_count')}` / `{cert.get('control_count')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        f"Folding solved: `{cert.get('folding_problem_solved')}`",
        "",
        "## Controls",
    ]
    for row in _as_list(cert.get("controls")):
        lines.extend([
            f"### {row.get('control_id')}",
            f"Passed: `{row.get('passed')}`",
            f"Observed valid: `{row.get('observed_valid')}`",
            f"Expected valid: `{row.get('expected_valid')}`",
            f"Observed failures: `{row.get('observed_failures')}`",
            f"Expected failure any: `{row.get('expected_failure_any')}`",
            f"Reason: {row.get('reason')}",
            "",
        ])
    lines.extend([
        "## Baseline content signature",
        "```json",
        json.dumps(cert.get("baseline_content_signature", {}), indent=2, sort_keys=True),
        "```",
        "",
        "## Locked interpretation",
        str(cert.get("locked_interpretation", "")),
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v34_kcsa_discriminative_content_controls_certificate.json"
    report_path = out_dir / "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_REPORT.md"
    decision_path = out_dir / "v34_kcsa_discriminative_content_controls_next_decision.json"
    decision = {
        "kind": "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_NEXT_DECISION_v0",
        "decision_status": cert.get("control_status"),
        "next_action": cert.get("next_action"),
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }
    cert = {
        **cert,
        "next_decision": decision,
        "artifacts": {
            "certificate": str(cert_path),
            "report": str(report_path),
            "decision": str(decision_path),
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V34 KcsA discriminative content controls.")
    parser.add_argument("--v33-certificate", type=Path, default=DEFAULT_V33_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v33 = _read_json(args.v33_certificate, "V33 certificate")
    cert = build_v34(v33)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "control_status": cert.get("control_status"),
        "control_count": cert.get("control_count"),
        "passed_control_count": cert.get("passed_control_count"),
        "claim_allowed": cert.get("claim_allowed"),
        "new_MD_allowed": cert.get("new_MD_allowed"),
        "positive_folding_evidence_found": cert.get("positive_folding_evidence_found"),
        "folding_problem_solved": cert.get("folding_problem_solved"),
        "next_action": cert.get("next_action"),
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
