#!/usr/bin/env python3
from __future__ import annotations

"""Run V49 reviewer attack-surface hardening audit."""

import argparse
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V49_REVIEWER_ATTACK_SURFACE_HARDENING"

PASSED = "V49_REVIEWER_ATTACK_SURFACE_HARDENING_PASSED"
FAILED = "V49_REVIEWER_ATTACK_SURFACE_HARDENING_FAILED"

REQUIRED_SOURCE_ROLE_FIELDS = [
    "source_class",
    "source_role",
    "spatial_proxy",
    "coordinate_derived",
    "internal_runtime",
    "allowed_for_prediction",
    "allowed_for_holdout",
    "allowed_for_claim",
    "leakage_risk",
    "rationale",
]

REQUIRED_SCORE_FIELDS = [
    "prediction_id",
    "target",
    "mechanism_class",
    "operator_bucket",
    "region_or_state",
    "predicted_effect",
    "perturbation",
    "expected_direction",
    "confidence",
    "prediction_source_ids",
    "falsification_criteria",
    "holdout_evidence_ids",
    "score_label",
    "score_reason",
    "scoring_pre_registered",
]

REQUIRED_SCORE_LABELS = [
    "supported",
    "partially_supported",
    "contradicted",
    "not_testable",
    "blocked_for_leakage",
]

EVIDENCE_CLASSES = [
    "pure_non_coordinate",
    "spatial_proxy_non_coordinate",
    "coordinate_derived",
    "internal_runtime",
]

V46_CERT_REL = "first_contact_clean_pharmacotopology_layer_run/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK/v46_cftr_f508del_membrane_multidomain_attack_certificate.json"
V46_REPORT_REL = "first_contact_clean_pharmacotopology_layer_run/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK_REPORT.md"


def _read_text(path: Path) -> str:
    if not path.exists():
        raise SystemExit(f"missing required file: {path}")
    return path.read_text(encoding="utf-8")


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {
        "control_id": control_id,
        "passed": bool(passed),
        "reason": reason,
        "observed": observed,
    }


def _section(text: str, heading: str) -> str:
    pattern = re.compile(rf"^## {re.escape(heading)}\s*$", re.MULTILINE)
    match = pattern.search(text)
    if not match:
        return ""
    next_heading = re.search(r"^## ", text[match.end():], re.MULTILINE)
    if not next_heading:
        return text[match.end():]
    return text[match.end(): match.end() + next_heading.start()]


def _has_all(text: str, needles: list[str]) -> bool:
    lowered = text.lower()
    return all(needle.lower() in lowered for needle in needles)


def negative_control_ids_from_readme(readme_text: str) -> list[str]:
    section = _section(readme_text, "Negative / Null Controls")
    return re.findall(r"^- `([^`]+)`", section, flags=re.MULTILINE)


def missing_source_role_fields(source_row: dict[str, Any]) -> list[str]:
    return [field for field in REQUIRED_SOURCE_ROLE_FIELDS if field not in source_row]


def _audit_controls(
    readme_text: str,
    taxonomy_text: str,
    rubric_text: str,
    inventory_text: str,
) -> list[dict[str, Any]]:
    readme_lower = readme_text.lower()
    taxonomy_lower = taxonomy_text.lower()
    rubric_lower = rubric_text.lower()
    inventory_lower = inventory_text.lower()
    negative_controls = negative_control_ids_from_readme(readme_text)
    return [
        _control(
            "readme_evidence_boundary_section_present",
            bool(_section(readme_text, "Evidence Boundary")),
            "README must contain an Evidence Boundary section.",
        ),
        _control(
            "readme_operator_scoring_section_present",
            bool(_section(readme_text, "Operator Scoring")),
            "README must contain an Operator Scoring section.",
        ),
        _control(
            "readme_latest_completed_full_cycle_section_present",
            bool(_section(readme_text, "Latest Completed Full Cycle")),
            "README must contain a Latest Completed Full Cycle section.",
        ),
        _control(
            "readme_negative_null_controls_section_present",
            bool(_section(readme_text, "Negative / Null Controls")),
            "README must contain a Negative / Null Controls section.",
        ),
        _control(
            "readme_distinguishes_pure_from_spatial_proxy",
            _has_all(readme_text, ["pure_non_coordinate", "spatial_proxy_non_coordinate", "must never be hidden inside pure"]),
            "README must distinguish pure non-coordinate evidence from spatial-proxy non-coordinate evidence.",
        ),
        _control(
            "readme_states_spatial_proxy_can_encode_spatial_information",
            "spatial-proxy evidence can encode spatial information" in readme_lower,
            "README must state that spatial-proxy evidence can encode spatial information.",
        ),
        _control(
            "readme_coordinate_derived_blocked_before_sealing",
            "coordinate-derived sources are blocked before sealing" in readme_lower,
            "README must state that coordinate-derived sources are blocked before sealing.",
        ),
        _control(
            "readme_lists_v46_certificate_and_report_paths",
            V46_CERT_REL in readme_text and V46_REPORT_REL in readme_text,
            "README must list V46 certificate and report paths.",
        ),
        _control(
            "readme_states_v46_completed_full_cycle",
            "v46 completed seal -> holdout -> validation" in readme_lower,
            "README must state that V46 completed seal -> holdout -> validation.",
        ),
        _control(
            "readme_lists_at_least_eight_negative_null_controls",
            len(set(negative_controls)) >= 8,
            "README must list at least eight negative/null controls.",
            {"control_count": len(set(negative_controls)), "control_ids": sorted(set(negative_controls))},
        ),
        _control(
            "operator_scoring_rubric_labels_frozen",
            all(label in rubric_lower for label in REQUIRED_SCORE_LABELS),
            "Scoring rubric must contain supported, partially_supported, contradicted, not_testable, and blocked_for_leakage labels.",
        ),
        _control(
            "evidence_taxonomy_forbids_internal_runtime_as_biological_evidence",
            "never biological prediction evidence" in taxonomy_lower and "never biological claim evidence" in taxonomy_lower,
            "Evidence taxonomy must forbid internal runtime artifacts as biological evidence.",
        ),
        _control(
            "readme_keeps_folding_problem_solved_false",
            re.search(r"folding_problem_solved\s*=\s*false", readme_lower) is not None,
            "README must preserve folding_problem_solved = false.",
        ),
        _control(
            "readme_forbids_universal_protein_folding_solved_claim",
            "universal protein folding is solved" in readme_lower and "not allowed" in readme_lower,
            "README must forbid a universal protein folding solved claim.",
        ),
        _control(
            "evidence_taxonomy_categories_present",
            all(category in taxonomy_text for category in EVIDENCE_CLASSES),
            "Evidence taxonomy must define all four evidence classes.",
        ),
        _control(
            "source_role_fields_frozen",
            all(field in taxonomy_text for field in REQUIRED_SOURCE_ROLE_FIELDS),
            "Evidence taxonomy must list the required source role fields.",
        ),
        _control(
            "operator_scoring_fields_frozen",
            all(field in rubric_text for field in REQUIRED_SCORE_FIELDS),
            "Operator scoring rubric must list the required prediction row fields.",
        ),
        _control(
            "negative_control_inventory_present",
            all(control_id in inventory_lower for control_id in [
                "random_sequence_control",
                "shuffled_sequence_control",
                "swapped_evidence_control",
                "wrong_target_control",
                "generic_annotation_only_control",
                "coordinate_leakage_control",
                "internal_runtime_leakage_control",
                "forced_wrong_grammar_control",
            ]),
            "Negative control inventory must contain the required null-control families.",
        ),
    ]


def audit_texts(
    readme_text: str,
    taxonomy_text: str,
    rubric_text: str,
    inventory_text: str,
    out_dir: Path = DEFAULT_OUT_DIR,
) -> dict[str, Any]:
    controls = _audit_controls(readme_text, taxonomy_text, rubric_text, inventory_text)
    failed = [control["control_id"] for control in controls if not control["passed"]]
    status = PASSED if not failed else FAILED
    negative_controls = negative_control_ids_from_readme(readme_text)
    cert_path = out_dir / "v49_reviewer_attack_surface_audit_certificate.json"
    report_path = out_dir / "V49_REVIEWER_ATTACK_SURFACE_HARDENING_REPORT.md"
    cert = {
        "kind": "V49_REVIEWER_ATTACK_SURFACE_HARDENING_CERTIFICATE_v0",
        "status": status,
        "run_mode": "reviewer_attack_surface_documentation_and_protocol_audit_no_MD_no_new_target",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_taxonomy_present": all(category in taxonomy_text for category in EVIDENCE_CLASSES),
        "spatial_proxy_boundary_present": _has_all(taxonomy_text, ["spatial_proxy_non_coordinate", "can encode spatial information", "must not be represented as `pure_non_coordinate`"]),
        "operator_scoring_rubric_present": all(field in rubric_text for field in REQUIRED_SCORE_FIELDS) and all(label in rubric_text for label in REQUIRED_SCORE_LABELS),
        "latest_completed_cycle_documented": all(control["passed"] for control in controls if control["control_id"] in {
            "readme_latest_completed_full_cycle_section_present",
            "readme_lists_v46_certificate_and_report_paths",
            "readme_states_v46_completed_full_cycle",
        }),
        "negative_control_inventory_present": any(control["control_id"] == "negative_control_inventory_present" and control["passed"] for control in controls),
        "v46_full_cycle_visible": all(control["passed"] for control in controls if control["control_id"] in {
            "readme_lists_v46_certificate_and_report_paths",
            "readme_states_v46_completed_full_cycle",
        }),
        "readme_claim_boundary_preserved": all(control["passed"] for control in controls if control["control_id"] in {
            "readme_keeps_folding_problem_solved_false",
            "readme_forbids_universal_protein_folding_solved_claim",
        }),
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "failed_checks": failed,
        "controls": controls,
        "negative_null_control_count": len(set(negative_controls)),
        "negative_null_control_ids": sorted(set(negative_controls)),
        "required_source_role_fields": REQUIRED_SOURCE_ROLE_FIELDS,
        "required_score_fields": REQUIRED_SCORE_FIELDS,
        "required_score_labels": REQUIRED_SCORE_LABELS,
        "folding_problem_solved": False,
        "claim_allowed": status == PASSED,
        "allowed_claim_text": (
            "The repository now makes the evidence boundary, operator scoring rubric, completed V46 cycle, and negative/null controls reviewer-visible. This is not a universal protein-folding solved claim."
            if status == PASSED else ""
        ),
        "forbidden_claims": [
            "universal protein folding is solved",
            "the project is an AlphaFold replacement",
            "coordinates were predicted de novo",
            "spatial-proxy evidence is the same as pure non-coordinate evidence",
            "internal runtime artifacts are biological evidence",
            "failed predictions may be repaired after holdout opening",
        ],
        "next_action": (
            "Use the V49 source-role fields and scoring rubric for the next live packet before opening holdouts."
            if status == PASSED else
            "Fix the failed reviewer audit checks before adding another target."
        ),
        "artifacts": {
            "certificate": str(cert_path),
            "report": str(report_path),
            "v46_certificate": str(REPO_ROOT / V46_CERT_REL),
            "v46_report": str(REPO_ROOT / V46_REPORT_REL),
        },
        "plain_english_interpretation": (
            "V49 does not add biology. It makes the reviewer-facing protocol sharper: evidence classes are explicit, spatial proxies are labeled, operator scoring is row-based, V46's full cycle is visible, and negative controls are required."
            if status == PASSED else
            "V49 found reviewer-facing gaps that must be fixed before the methodology is considered hardened."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V49 Reviewer Attack-Surface Hardening",
        "",
        f"Status: `{cert['status']}`",
        f"folding_problem_solved: `{cert['folding_problem_solved']}`",
        f"claim_allowed: `{cert['claim_allowed']}`",
        "",
        "## Reviewer Questions",
        f"- Evidence taxonomy present: `{cert['evidence_taxonomy_present']}`",
        f"- Spatial-proxy boundary present: `{cert['spatial_proxy_boundary_present']}`",
        f"- Operator scoring rubric present: `{cert['operator_scoring_rubric_present']}`",
        f"- Latest completed cycle documented: `{cert['latest_completed_cycle_documented']}`",
        f"- Negative control inventory present: `{cert['negative_control_inventory_present']}`",
        f"- V46 full cycle visible: `{cert['v46_full_cycle_visible']}`",
        f"- README claim boundary preserved: `{cert['readme_claim_boundary_preserved']}`",
        "",
        "## Controls",
        f"Passed `{cert['passed_control_count']}` / `{cert['control_count']}`.",
    ]
    if cert["failed_checks"]:
        lines.extend(["", "## Failed Checks"])
        for check in cert["failed_checks"]:
            lines.append(f"- `{check}`")
    lines.extend([
        "",
        "## Negative / Null Controls",
    ])
    for control_id in cert["negative_null_control_ids"]:
        lines.append(f"- `{control_id}`")
    lines.extend([
        "",
        "## Artifacts",
        f"- Certificate: `{cert['artifacts']['certificate']}`",
        f"- Report: `{cert['artifacts']['report']}`",
        f"- V46 certificate: `{cert['artifacts']['v46_certificate']}`",
        f"- V46 report: `{cert['artifacts']['v46_report']}`",
        "",
        "## Plain English Interpretation",
        cert["plain_english_interpretation"],
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v49(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    readme_text = _read_text(REPO_ROOT / "README.md")
    taxonomy_text = _read_text(REPO_ROOT / "docs" / "EVIDENCE_TAXONOMY.md")
    rubric_text = _read_text(REPO_ROOT / "docs" / "OPERATOR_SCORING_RUBRIC.md")
    inventory_text = _read_text(REPO_ROOT / "docs" / "NEGATIVE_CONTROL_INVENTORY.md")
    out_dir.mkdir(parents=True, exist_ok=True)
    cert = audit_texts(readme_text, taxonomy_text, rubric_text, inventory_text, out_dir=out_dir)
    cert_path = out_dir / "v49_reviewer_attack_surface_audit_certificate.json"
    report_path = out_dir / "V49_REVIEWER_ATTACK_SURFACE_HARDENING_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V49 reviewer attack-surface hardening audit.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v49(args.out_dir)
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "failed_checks": cert["failed_checks"],
        "folding_problem_solved": cert["folding_problem_solved"],
        "claim_allowed": cert["claim_allowed"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
