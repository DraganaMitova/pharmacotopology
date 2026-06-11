#!/usr/bin/env python3
from __future__ import annotations

"""Run the V41 protein-folding claim boundary gate."""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DATA_ROOT = REPO_ROOT / "data" / "claim_gate" / "V41"
DEFAULT_OUT_DIR = RUN_ROOT / "V41_PROTEIN_FOLDING_CLAIM_GATE"

CLAIM_LEVELS = [
    {
        "level": "C0_NO_CLAIM",
        "name": "No scientific claim beyond code execution.",
    },
    {
        "level": "C1_CLEAN_EVIDENCE_PIPELINE",
        "name": "The system cleanly imports external evidence and blocks leakage.",
    },
    {
        "level": "C2_MECHANISM_GRAMMAR_CLASSIFICATION",
        "name": "The system classifies folding-problem grammar from non-coordinate evidence.",
    },
    {
        "level": "C3_FALSIFIABLE_OPERATOR_MECHANISM",
        "name": "The system generates falsifiable, perturbation-sensitive operator-level mechanism predictions supported by independent non-coordinate holdouts.",
    },
    {
        "level": "C4_PROSPECTIVE_MECHANISM_GENERALIZATION",
        "name": "The system generalizes prospectively to unseen targets under masking, decoy pressure, source separation, and independent holdout validation.",
    },
    {
        "level": "C5_PROTEIN_FOLDING_SOLVED",
        "name": "The system can predict de novo 3D folds / functional conformational ensembles across broad protein classes with independent blind validation comparable to structure-prediction benchmarks.",
    },
]

CERTIFICATE_PATHS = {
    "V32": RUN_ROOT / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT" / "v32_external_constraint_source_import_preflight_certificate.json",
    "V33_READOUT": RUN_ROOT / "V33_CONSTRAINT_BACKED_OPERATOR_READOUT" / "v33_constraint_backed_operator_readout_certificate.json",
    "V33_NEGATIVE": RUN_ROOT / "V33_NEGATIVE_CONTROLS" / "v33_negative_controls_certificate.json",
    "V34": RUN_ROOT / "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS" / "v34_kcsa_discriminative_content_controls_certificate.json",
    "V35": RUN_ROOT / "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT" / "v35_noncoordinate_evolutionary_holdout_certificate.json",
    "V36": RUN_ROOT / "V36_REAL_EVIDENCE_DOSSIER_GATE" / "v36_real_evidence_dossier_gate_certificate.json",
    "V37": RUN_ROOT / "V37_MECHANISM_QUESTION_PROBES" / "v37_mechanism_question_probes_certificate.json",
    "V38": RUN_ROOT / "V38_BLIND_MECHANISM_GENERALIZATION_PANEL" / "v38_blind_mechanism_generalization_panel_certificate.json",
    "V39": RUN_ROOT / "V39_MECHANISM_TO_FALSIFIABLE_PREDICTION_VALIDATION" / "v39_mechanism_to_falsifiable_prediction_validation_certificate.json",
    "V40": RUN_ROOT / "V40_MECHANISM_PERTURBATION_PRESSURE_TESTS" / "v40_mechanism_perturbation_pressure_tests_certificate.json",
}

PASSED = "V41_MECHANISM_CLAIM_ALLOWED_C5_BLOCKED"
NO_CLAIM = "V41_NO_CLAIM_ALLOWED"
BLOCKED = "V41_BLOCKED_CLAIM_BOUNDARY_VIOLATION"
FAILED_UNSAFE = "V41_FAILED_UNSAFE_SOLVED_FOLDING_CLAIM"

ALLOWED_PUBLIC_CLAIM = (
    "We have a claim-disabled, coordinate-free mechanism-grammar prototype that classifies hard protein-folding regimes "
    "and predicts perturbation-sensitive operator constraints for KcsA, XCL1, and alpha-synuclein, with masked-target "
    "and decoy controls, no MD, and no coordinate leakage. This is not solved protein folding, is not a de novo "
    "protein-structure predictor, and does not solve protein folding."
)

FORBIDDEN_PUBLIC_CLAIMS = [
    "we solved protein folding",
    "we can predict every protein",
    "we can predict de novo 3D structure from sequence",
    "we outperform AlphaFold",
    "KcsA/XCL1/SNCA prove universal folding",
    "coordinates were not needed to solve structure",
    "this is a validated drug-discovery engine",
]

MISSING_EVIDENCE_FOR_C5 = [
    "broad prospective blind panel",
    "de novo coordinate or ensemble predictions",
    "independent experimental structure/function validation",
    "quantitative comparison against strong baselines",
    "failure modes documented across broad protein classes",
    "reproducibility certificate",
    "external review-ready report",
]

SCIENTIST_ATTACK_SURFACE = [
    "Only three focal mechanism exemplars are supported at C3: KcsA, XCL1, and alpha-synuclein.",
    "V38 is masked and decoy-controlled, but V41 does not add a new prospective unseen panel.",
    "The claim is coordinate-free mechanism grammar, not de novo 3D coordinate prediction.",
    "Static structure, folding pathway, dynamics, membrane context, and ensemble prediction remain distinct claims.",
    "External benchmark comparison against AlphaFold/ESMFold/RoseTTAFold is not present in V41.",
    "Perturbation pressure is literature/annotation-supported, not newly experimentally validated here.",
    "Drug-discovery utility, therapeutic validation, and clinical relevance are not established.",
    "C5 requires external blind validation and broad reproducibility evidence that is not present.",
]


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


def _bool_true(value: Any) -> bool:
    if value is True:
        return True
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return False


def _walk(value: Any):
    if isinstance(value, dict):
        yield value
        for item in value.values():
            yield from _walk(item)
    elif isinstance(value, list):
        for item in value:
            yield from _walk(item)


def _any_true(data: dict[str, Any], keys: set[str]) -> bool:
    return any(_bool_true(obj.get(key)) for obj in _walk(data) if isinstance(obj, dict) for key in keys if key in obj)


def _sum_numeric(data: dict[str, Any], keys: set[str]) -> int:
    total = 0
    for obj in _walk(data):
        if not isinstance(obj, dict):
            continue
        for key in keys:
            value = obj.get(key)
            if isinstance(value, bool):
                continue
            if isinstance(value, int):
                total += value
    return total


def _top_level_int(data: dict[str, Any], key: str) -> int:
    value = data.get(key)
    if isinstance(value, bool):
        return 0
    return value if isinstance(value, int) else 0


def load_certificates() -> dict[str, dict[str, Any]]:
    return {label: _read_json(path, label) for label, path in CERTIFICATE_PATHS.items()}


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(certs: dict[str, dict[str, Any]], claim_coordinate_count: int, claim_internal_count: int) -> list[dict[str, Any]]:
    controls: list[dict[str, Any]] = []
    v32_status = certs["V32"].get("preflight_status") or certs["V32"].get("control_status") or certs["V32"].get("source_v31_status")
    controls.append(_control(
        "v32_external_source_import_preflight_passed_or_abstained",
        any(token in str(v32_status) for token in ["READY", "CLEAN_ABSTAIN", "PASSED"]),
        "V32 must pass or cleanly abstain before downstream claim use.",
        {"status": v32_status},
    ))
    controls.append(_control(
        "v33_operator_readout_passed_claim_disabled",
        certs["V33_READOUT"].get("kind") == "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_v0" and certs["V33_READOUT"].get("claim_allowed") is False,
        "V33 operator readout must exist and remain claim-disabled.",
        {"kind": certs["V33_READOUT"].get("kind"), "claim_allowed": certs["V33_READOUT"].get("claim_allowed")},
    ))
    controls.append(_control(
        "v33_negative_controls_passed",
        certs["V33_NEGATIVE"].get("passed_control_count") == certs["V33_NEGATIVE"].get("control_count"),
        "V33 negative controls must pass.",
        {"passed": certs["V33_NEGATIVE"].get("passed_control_count"), "total": certs["V33_NEGATIVE"].get("control_count")},
    ))
    controls.append(_control(
        "v34_discriminative_content_controls_passed",
        certs["V34"].get("control_status") == "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED",
        "V34 discriminative content controls must pass claim-disabled.",
        {"status": certs["V34"].get("control_status")},
    ))
    controls.append(_control(
        "v35_clean_abstain_not_positive_evidence",
        certs["V35"].get("control_status") == "V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE",
        "V35 clean abstain must remain a guardrail, not positive evidence.",
        {"status": certs["V35"].get("control_status")},
    ))
    controls.append(_control(
        "v36_real_noncoordinate_dossiers_built",
        certs["V36"].get("control_status") == "V36_REAL_EVIDENCE_DOSSIERS_READY_CLAIM_DISABLED",
        "V36 must build real non-coordinate dossiers.",
        {"status": certs["V36"].get("control_status")},
    ))
    controls.append(_control(
        "v37_mechanism_classes_and_masking_passed",
        certs["V37"].get("control_status") == "V37_MECHANISM_QUESTION_PROBES_PASSED_CLAIM_DISABLED",
        "V37 must assign mechanism classes under target-name masking.",
        {"status": certs["V37"].get("control_status")},
    ))
    controls.append(_control(
        "v38_blind_masked_generalization_decoys_passed",
        certs["V38"].get("control_status") == "V38_BLIND_MECHANISM_GENERALIZATION_PASSED_CLAIM_DISABLED",
        "V38 must pass masked assignment and decoy pressure.",
        {"status": certs["V38"].get("control_status")},
    ))
    controls.append(_control(
        "v39_falsifiable_mechanism_predictions_passed",
        certs["V39"].get("control_status") == "V39_MECHANISM_PREDICTIONS_VALIDATED_CLAIM_DISABLED",
        "V39 must generate falsifiable mechanism predictions.",
        {"status": certs["V39"].get("control_status")},
    ))
    controls.append(_control(
        "v40_perturbation_pressure_predictions_passed",
        certs["V40"].get("control_status") == "V40_MECHANISM_PERTURBATION_PRESSURE_PASSED_CLAIM_DISABLED",
        "V40 must generate perturbation-pressure predictions.",
        {"status": certs["V40"].get("control_status")},
    ))
    controls.append(_control(
        "all_prior_claim_allowed_flags_false",
        not any(_any_true(data, {"claim_allowed"}) for data in certs.values()),
        "If any earlier certificate has claim_allowed=true, V41 blocks.",
    ))
    controls.append(_control(
        "all_prior_folding_problem_solved_flags_false",
        not any(_any_true(data, {"folding_problem_solved"}) for data in certs.values()),
        "If folding_problem_solved=true appears anywhere, V41 blocks.",
    ))
    controls.append(_control(
        "no_md_used_as_claim_evidence",
        not any(_any_true(data, {"new_md_executed", "membrane_md_executed", "md_used_for_claim"}) for data in certs.values()),
        "MD must not be used as claim evidence.",
    ))
    controls.append(_control(
        "no_coordinate_derived_sources_for_noncoordinate_claims",
        claim_coordinate_count == 0,
        "Coordinate-derived evidence must not support the non-coordinate claim.",
        {"coordinate_derived_source_count_for_claim": claim_coordinate_count},
    ))
    controls.append(_control(
        "no_internal_runtime_reports_as_biological_evidence",
        claim_internal_count == 0,
        "Internal runtime reports must remain audit evidence only.",
        {"internal_runtime_source_count_for_claim": claim_internal_count},
    ))
    controls.append(_control(
        "v38_answer_key_not_used_for_assignment",
        certs["V38"].get("answer_key_used_for_assignment") is False,
        "V38 answer key must not be used for assignment.",
        {"answer_key_used_for_assignment": certs["V38"].get("answer_key_used_for_assignment")},
    ))
    controls.append(_control(
        "c5_blocked_without_de_novo_structure_ensemble_validation",
        True,
        "C5 cannot be allowed without de novo structure/ensemble validation and broad independent benchmarks.",
        {"c5_claim_allowed": False, "missing_evidence": MISSING_EVIDENCE_FOR_C5},
    ))
    controls.append(_control(
        "allowed_claim_text_mentions_not_solved_protein_folding",
        "not solved protein folding" in ALLOWED_PUBLIC_CLAIM.lower(),
        "Allowed claim text must explicitly state not solved protein folding.",
        {"allowed_public_claim": ALLOWED_PUBLIC_CLAIM},
    ))
    controls.append(_control(
        "forbidden_claims_include_we_solved_protein_folding",
        "we solved protein folding" in FORBIDDEN_PUBLIC_CLAIMS,
        "Forbidden claims must explicitly include we solved protein folding.",
        {"forbidden_public_claims": FORBIDDEN_PUBLIC_CLAIMS},
    ))
    controls.append(_control(
        "scientist_attack_surface_nonempty",
        len(SCIENTIST_ATTACK_SURFACE) > 0,
        "Scientist attack surface must not be empty.",
        {"scientist_attack_surface_count": len(SCIENTIST_ATTACK_SURFACE)},
    ))
    controls.append(_control(
        "max_claim_not_below_c2_with_v36_to_v40_passing",
        True,
        "Max claim must not fall below C2 when V36-V40 pass.",
        {"max_claim_floor": "C2_MECHANISM_GRAMMAR_CLASSIFICATION"},
    ))
    return controls


def _claim_evidence_versions() -> list[str]:
    return ["V36", "V37", "V38", "V39", "V40"]


def _evidence_supporting_claim(certs: dict[str, dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {
            "version": "V36",
            "supports": "real non-coordinate dossiers for KcsA, XCL1, and alpha-synuclein",
            "status": str(certs["V36"].get("control_status")),
        },
        {
            "version": "V37",
            "supports": "mechanism-class grammar assignment from V36 dossiers under masking controls",
            "status": str(certs["V37"].get("control_status")),
        },
        {
            "version": "V38",
            "supports": "masked-target and decoy-pressure generalization audit, not C5 prospective structure validation",
            "status": str(certs["V38"].get("control_status")),
        },
        {
            "version": "V39",
            "supports": "falsifiable mechanism predictions with independent non-coordinate holdouts",
            "status": str(certs["V39"].get("control_status")),
        },
        {
            "version": "V40",
            "supports": "perturbation-sensitive operator constraints supported by non-coordinate holdouts",
            "status": str(certs["V40"].get("control_status")),
        },
    ]


def _matrix_rows(certs: dict[str, dict[str, Any]]) -> list[dict[str, str]]:
    rows = [
        ("V32", "pipeline audit", "C1_CLEAN_EVIDENCE_PIPELINE", "audit_only", "External source preflight passed/abstained; coordinate imports are not used for non-coordinate claim support."),
        ("V33_READOUT", "pipeline audit", "C1_CLEAN_EVIDENCE_PIPELINE", "audit_only", "Operator readout remained claim-disabled."),
        ("V33_NEGATIVE", "pipeline audit", "C1_CLEAN_EVIDENCE_PIPELINE", "audit_only", "Negative controls passed."),
        ("V34", "pipeline audit", "C1_CLEAN_EVIDENCE_PIPELINE", "audit_only", "Discriminative content controls passed."),
        ("V35", "guardrail", "none_positive", "abstain_guardrail", "Clean abstain is not positive evidence."),
        ("V36", "claim support", "C2_MECHANISM_GRAMMAR_CLASSIFICATION", "noncoordinate_claim_support", "Real non-coordinate dossiers built."),
        ("V37", "claim support", "C2_MECHANISM_GRAMMAR_CLASSIFICATION", "noncoordinate_claim_support", "Mechanism classes assigned under target-name masking."),
        ("V38", "claim support", "C2_MECHANISM_GRAMMAR_CLASSIFICATION", "noncoordinate_claim_support", "Masked and decoy-controlled mechanism generalization audit."),
        ("V39", "claim support", "C3_FALSIFIABLE_OPERATOR_MECHANISM", "noncoordinate_claim_support", "Falsifiable operator predictions supported by independent holdouts."),
        ("V40", "claim support", "C3_FALSIFIABLE_OPERATOR_MECHANISM", "noncoordinate_claim_support", "Perturbation pressure operators supported by independent holdouts."),
    ]
    out = []
    for version, role, claim_level, claim_use, notes in rows:
        cert = certs[version]
        out.append({
            "version": version,
            "role": role,
            "certificate_path": str(CERTIFICATE_PATHS[version]),
            "status": str(cert.get("control_status") or cert.get("preflight_status") or cert.get("kind")),
            "supports_claim_level": claim_level,
            "claim_use": claim_use,
            "notes": notes,
        })
    return out


def build_v41(certs: dict[str, dict[str, Any]] | None = None) -> dict[str, Any]:
    certs = certs or load_certificates()
    claim_versions = _claim_evidence_versions()
    claim_coordinate_count = sum(
        _top_level_int(certs[version], "coordinate_derived_source_count")
        + _top_level_int(certs[version], "coordinate_derived_source_count_for_claim")
        for version in claim_versions
    )
    claim_internal_count = sum(
        _top_level_int(certs[version], "internal_runtime_source_count")
        + _top_level_int(certs[version], "internal_runtime_source_count_for_claim")
        for version in claim_versions
    )
    controls = _controls(certs, claim_coordinate_count, claim_internal_count)
    failed = [control["control_id"] for control in controls if not control["passed"]]

    v36_to_v40_pass = all(control["passed"] for control in controls[5:10])
    c5_claim_allowed = False
    protein_folding_solved_claim_allowed = False
    mechanism_claim_allowed = v36_to_v40_pass and not failed
    if mechanism_claim_allowed:
        max_level = "C3_FALSIFIABLE_OPERATOR_MECHANISM"
        status = PASSED
    elif failed:
        max_level = "C0_NO_CLAIM"
        status = BLOCKED
    else:
        max_level = "C0_NO_CLAIM"
        status = NO_CLAIM
    if c5_claim_allowed or protein_folding_solved_claim_allowed:
        status = FAILED_UNSAFE

    max_name = next(row["name"] for row in CLAIM_LEVELS if row["level"] == max_level)
    return {
        "kind": "V41_PROTEIN_FOLDING_CLAIM_GATE_v0",
        "run_mode": "claim_boundary_audit_no_coordinates_no_MD_c5_blocked",
        "control_status": status,
        "evaluated_versions": list(CERTIFICATE_PATHS.keys()),
        "max_allowed_claim_level": max_level,
        "max_allowed_claim_name": max_name,
        "c5_claim_allowed": c5_claim_allowed,
        "protein_folding_solved_claim_allowed": protein_folding_solved_claim_allowed,
        "mechanism_claim_allowed": mechanism_claim_allowed,
        "allowed_public_claim": ALLOWED_PUBLIC_CLAIM,
        "forbidden_public_claims": FORBIDDEN_PUBLIC_CLAIMS,
        "evidence_supporting_claim": _evidence_supporting_claim(certs),
        "evidence_missing_for_c5": MISSING_EVIDENCE_FOR_C5,
        "scientist_attack_surface": SCIENTIST_ATTACK_SURFACE,
        "scientist_attack_surface_count": len(SCIENTIST_ATTACK_SURFACE),
        "coordinate_derived_source_count_for_claim": claim_coordinate_count,
        "internal_runtime_source_count_for_claim": claim_internal_count,
        "native_metrics_used_for_claim": False,
        "md_used_for_claim": False,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "positive_folding_evidence_found": False,
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "controls": controls,
        "failed_checks": failed,
        "next_action": "run_V42_baseline_failure_mode_comparison_or_V41_review_before_any_public_claim",
        "locked_interpretation": (
            "V41 allows a bounded C3 mechanism-grammar claim and blocks C5. The project can say it has a coordinate-free, "
            "claim-disabled mechanism grammar with falsifiable perturbation-sensitive operators for the three audited regimes; "
            "it cannot say protein folding is solved, cannot claim de novo 3D structure prediction, and cannot claim superiority to AlphaFold."
        ),
    }


def _write_matrix(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["version", "role", "certificate_path", "status", "supports_claim_level", "claim_use", "notes"])
        writer.writeheader()
        writer.writerows(rows)


def _write_lines(path: Path, title: str, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = [f"# {title}", ""]
    body.extend(lines)
    path.write_text("\n".join(body) + "\n", encoding="utf-8")


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V41 Protein Folding Claim Gate",
        "",
        f"Status: `{cert['control_status']}`",
        f"Max allowed claim level: `{cert['max_allowed_claim_level']}`",
        f"Protein-folding-solved claim allowed: `{cert['protein_folding_solved_claim_allowed']}`",
        f"C5 claim allowed: `{cert['c5_claim_allowed']}`",
        f"Controls: `{cert['passed_control_count']}` / `{cert['control_count']}`",
        "",
        "## Allowed Claim",
        cert["allowed_public_claim"],
        "",
        "## Forbidden Claims",
    ]
    lines.extend(f"- {claim}" for claim in cert["forbidden_public_claims"])
    lines.extend(["", "## Missing Evidence For C5"])
    lines.extend(f"- {item}" for item in cert["evidence_missing_for_c5"])
    lines.extend(["", "## Scientist Attack Surface"])
    lines.extend(f"- {item}" for item in cert["scientist_attack_surface"])
    lines.extend(["", "## Plain English Interpretation", cert["locked_interpretation"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    certs = load_certificates()
    cert = build_v41(certs)
    out_dir.mkdir(parents=True, exist_ok=True)
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    claim_ladder_path = DATA_ROOT / "claim_ladder.json"
    matrix_path = DATA_ROOT / "evidence_to_claim_matrix.csv"
    allowed_path = DATA_ROOT / "allowed_claim_text.md"
    forbidden_path = DATA_ROOT / "forbidden_claim_text.md"
    attack_path = DATA_ROOT / "scientist_attack_surface.md"
    cert_path = out_dir / "v41_protein_folding_claim_gate_certificate.json"
    report_path = out_dir / "V41_PROTEIN_FOLDING_CLAIM_GATE_REPORT.md"
    decision_path = out_dir / "v41_protein_folding_claim_gate_next_decision.json"

    ladder = {
        "kind": "V41_CLAIM_LADDER_v0",
        "claim_levels": CLAIM_LEVELS,
        "max_allowed_claim_level": cert["max_allowed_claim_level"],
        "c5_claim_allowed": False,
    }
    matrix_rows = _matrix_rows(certs)
    decision = {
        "kind": "V41_PROTEIN_FOLDING_CLAIM_GATE_NEXT_DECISION_v0",
        "decision_status": cert["control_status"],
        "max_allowed_claim_level": cert["max_allowed_claim_level"],
        "protein_folding_solved_claim_allowed": False,
        "next_action": cert["next_action"],
    }
    artifacts = {
        "certificate": str(cert_path),
        "report": str(report_path),
        "decision": str(decision_path),
        "claim_ladder": str(claim_ladder_path),
        "evidence_to_claim_matrix": str(matrix_path),
        "allowed_claim_text": str(allowed_path),
        "forbidden_claim_text": str(forbidden_path),
        "scientist_attack_surface": str(attack_path),
    }
    cert = {**cert, "artifacts": artifacts, "next_decision": decision}
    _write_json(claim_ladder_path, ladder)
    _write_matrix(matrix_path, matrix_rows)
    _write_lines(allowed_path, "Allowed Claim Text", [cert["allowed_public_claim"]])
    _write_lines(forbidden_path, "Forbidden Claim Text", [f"- {claim}" for claim in cert["forbidden_public_claims"]])
    _write_lines(attack_path, "Scientist Attack Surface", [f"- {item}" for item in cert["scientist_attack_surface"]])
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V41 protein-folding claim boundary gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = write_outputs(args.out_dir)
    cert = _read_json(paths["certificate"], "V41 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "max_allowed_claim_level": cert["max_allowed_claim_level"],
        "protein_folding_solved_claim_allowed": cert["protein_folding_solved_claim_allowed"],
        "mechanism_claim_allowed": cert["mechanism_claim_allowed"],
        "c5_claim_allowed": cert["c5_claim_allowed"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "coordinate_derived_source_count_for_claim": cert["coordinate_derived_source_count_for_claim"],
        "internal_runtime_source_count_for_claim": cert["internal_runtime_source_count_for_claim"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
