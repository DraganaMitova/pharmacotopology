#!/usr/bin/env python3
from __future__ import annotations

"""Run V60 bounded Protein Esperanto solved-claim gate.

V60 does not run new science.  It reads V50-V59 certificates and computes the
allowed claim boundary.
"""

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import stable_hash  # noqa: E402


DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V60"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V60_BOUNDED_PROTEIN_ESPERANTO_SOLVED_CLAIM_GATE"

STRONG_PASS = "V60_BOUNDED_PROCESS_REPLICATION_CLAIM_PASSED_REVIEW_REQUIRED"
BEHAVIOR_ONLY = "V60_REAL_SEQUENCE_BEHAVIOR_SUPPORTED_PROCESS_CLAIM_BLOCKED"
BLOCKED = "V60_SOLVED_CLAIM_BLOCKED_ENGINE_REVISION_REQUIRED"

INPUT_CERTIFICATES = {
    "v50_v56": RUN_ROOT / "V50_V56_PROTEIN_ESPERANTO_ENGINE" / "v50_v56_protein_esperanto_engine_certificate.json",
    "v57": RUN_ROOT / "V57_BLIND_PROTEIN_ESPERANTO_GENERALIZATION_GATE" / "v57_blind_protein_esperanto_generalization_certificate.json",
    "v58": RUN_ROOT / "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_GATE" / "v58_real_sequence_time_blind_folding_replication_certificate.json",
    "v58_audit": RUN_ROOT / "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT" / "v58_error_autopsy_and_calibration_audit_certificate.json",
    "v59": RUN_ROOT / "V59_REAL_FOLDING_PROCESS_REPLICATION_PANEL" / "v59_real_folding_process_replication_certificate.json",
}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    return path


def load_claim_inputs() -> dict[str, dict[str, Any]]:
    return {name: _read_json(path, name) for name, path in INPUT_CERTIFICATES.items()}


def _status(cert: dict[str, Any]) -> str:
    return str(cert.get("status", ""))


def _input_summary(name: str, cert: dict[str, Any], passed: bool, reason: str) -> dict[str, Any]:
    return {
        "input": name,
        "status": _status(cert),
        "passed_for_v60": bool(passed),
        "reason": reason,
        "certificate_hash": cert.get("certificate_hash"),
    }


def _control_observed(cert: dict[str, Any], control_id: str) -> dict[str, Any]:
    for row in cert.get("controls", []):
        if row.get("control_id") == control_id and isinstance(row.get("observed"), dict):
            return row["observed"]
    return {}


def compute_v60(inputs: dict[str, dict[str, Any]]) -> dict[str, Any]:
    v50 = inputs["v50_v56"]
    v57 = inputs["v57"]
    v58 = inputs["v58"]
    v58_audit = inputs["v58_audit"]
    v59 = inputs["v59"]

    v50_ok = (
        _status(v50) == "V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_PASSED_REVIEW_REQUIRED"
        and v50.get("passed_control_count") == v50.get("control_count")
        and v50.get("folding_problem_solved") is False
    )
    v57_ok = (
        _status(v57) == "V57_BLIND_GENERALIZATION_PASSED_REVIEW_REQUIRED"
        and v57.get("passed_control_count") == v57.get("control_count")
        and v57.get("folding_problem_solved") is False
    )
    v58_ok = (
        _status(v58) == "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_PASSED_REVIEW_REQUIRED"
        and v58.get("target_count") == 20
        and v58.get("passed_control_count") == v58.get("control_count")
        and v58.get("failure_cases_reported") is True
        and v58.get("folding_problem_solved") is False
        and v58.get("coordinate_truth_used_before_seal") is False
        and v58.get("alphafold_used_before_seal") is False
    )
    v58_audit_ok = (
        _status(v58_audit) == "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_COMPLETED_REVIEW_REQUIRED"
        and v58_audit.get("accepted_accuracy") == 1.0
        and v58_audit.get("failure_preserved") is True
        and v58_audit.get("engine_biology_modified") is False
    )
    v59_ok = (
        _status(v59) == "V59_REAL_FOLDING_PROCESS_REPLICATION_PASSED_REVIEW_REQUIRED"
        and v59.get("passed_control_count") == v59.get("control_count")
        and v59.get("accepted_process_accuracy") == 1.0
        and v59.get("process_holdouts_opened_after_seal") is True
        and v59.get("failures_preserved") is True
        and v59.get("coordinate_truth_used_before_seal") is False
        and v59.get("alphafold_used_before_seal") is False
    )
    behavior_supported = v50_ok and v57_ok and v58_ok and v58_audit_ok
    process_supported = behavior_supported and v59_ok
    v59_engine_revision = _control_observed(v59, "engine_revision_explicit_and_bounded")
    if process_supported:
        outcome = STRONG_PASS
        claim_allowed = True
        bounded_process_supported = True
        allowed_claim = (
            "The protein folding problem is not universally solved here; however, after the V59-required bounded engine hardening, "
            "the current Protein Esperanto engine supports a bounded solution claim on accepted targets: it predicts folding regime, "
            "important regions, topology/observable behavior, and selected folding-process observables under source-separated validation."
        )
    elif behavior_supported:
        outcome = BEHAVIOR_ONLY
        claim_allowed = True
        bounded_process_supported = False
        allowed_claim = (
            "The engine predicts real-sequence folding regimes, important regions, and topology/observable proxies on accepted targets, "
            "but folding-process replication remains unproven."
        )
    else:
        outcome = BLOCKED
        claim_allowed = False
        bounded_process_supported = False
        allowed_claim = "The current engine is promising but not yet claim-ready beyond internal/restricted validation."

    input_rows = [
        _input_summary("V50-V56 engine suite", v50, v50_ok, "engine suite and controls"),
        _input_summary("V57 blind generalization", v57, v57_ok, "blind generalization and controls"),
        _input_summary("V58 real-sequence gate", v58, v58_ok, "real-sequence behavior and leakage controls"),
        _input_summary("V58 calibration audit", v58_audit, v58_audit_ok, "accepted subset calibration and preserved failures"),
        _input_summary("V59 process panel", v59, v59_ok, "real folding-process evidence and controls"),
    ]
    cert = {
        "kind": "V60_BOUNDED_PROTEIN_ESPERANTO_SOLVED_CLAIM_CERTIFICATE_v0",
        "status": outcome,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_rows": input_rows,
        "protein_esperanto_engine_claim_allowed": claim_allowed,
        "bounded_accepted_target_process_replication_supported": bounded_process_supported,
        "bounded_real_sequence_behavior_supported": behavior_supported,
        "universal_folding_problem_solved": False,
        "atomistic_folding_solved": False,
        "all_sequence_prediction_solved": False,
        "folding_problem_solved": False,
        "alphafold_replaced": False,
        "atomistic_md_performed": False,
        "readme_touched": False,
        "engine_revision_after_v58_required": bool(v59.get("engine_biology_modified")),
        "engine_revision_scope": v59_engine_revision.get("engine_revision_scope", []),
        "engine_operator_set_modified_in_v59": any(
            row.get("control_id") == "engine_revision_explicit_and_bounded"
            and row.get("observed", {}).get("engine_operator_set_modified") is True
            for row in v59.get("controls", [])
        ),
        "engine_mechanism_class_set_modified_in_v59": any(
            row.get("control_id") == "engine_revision_explicit_and_bounded"
            and row.get("observed", {}).get("engine_mechanism_class_set_modified") is True
            for row in v59.get("controls", [])
        ),
        "allowed_claim": allowed_claim,
        "forbidden_claims": [
            "Universal protein folding is solved.",
            "All sequences are solved.",
            "Atomistic folding is solved.",
            "AlphaFold is replaced.",
            "Every unsupported target can be answered instead of abstained.",
            "A V59 engine revision can be hidden under a frozen-engine claim.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V60 Bounded Protein Esperanto Solved-Claim Gate",
        "",
        f"Status: `{cert['status']}`",
        f"Claim allowed: `{cert['protein_esperanto_engine_claim_allowed']}`",
        f"Bounded process replication supported: `{cert['bounded_accepted_target_process_replication_supported']}`",
        f"Universal folding problem solved: `{cert['universal_folding_problem_solved']}`",
        f"Atomistic folding solved: `{cert['atomistic_folding_solved']}`",
        f"All sequence prediction solved: `{cert['all_sequence_prediction_solved']}`",
        f"Engine revision after V58 required: `{cert['engine_revision_after_v58_required']}`",
        "",
        "## Inputs",
    ]
    for row in cert["input_rows"]:
        lines.append(f"- `{row['input']}` status `{row['status']}` passed `{row['passed_for_v60']}`")
    lines.extend([
        "",
        "## Allowed Claim",
        cert["allowed_claim"],
        "",
        "## Forbidden Claims",
    ])
    for claim in cert["forbidden_claims"]:
        lines.append(f"- {claim}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v60(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    inputs = load_claim_inputs()
    cert = compute_v60(inputs)
    data_cert_path = _write_json(DATA_ROOT / "v60_bounded_solved_claim_certificate.json", cert)
    data_report_path = DATA_ROOT / "V60_BOUNDED_PROTEIN_ESPERANTO_SOLVED_CLAIM_GATE_REPORT.md"
    _write_report(data_report_path, cert)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v60_bounded_protein_esperanto_solved_claim_certificate.json"
    report_path = out_dir / "V60_BOUNDED_PROTEIN_ESPERANTO_SOLVED_CLAIM_GATE_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert)
    return {
        "data_certificate": data_cert_path,
        "data_report": data_report_path,
        "certificate": cert_path,
        "report": report_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V60 bounded Protein Esperanto solved-claim gate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v60(args.out_dir)
    cert = _read_json(paths["certificate"], "V60 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "protein_esperanto_engine_claim_allowed": cert["protein_esperanto_engine_claim_allowed"],
        "bounded_accepted_target_process_replication_supported": cert["bounded_accepted_target_process_replication_supported"],
        "universal_folding_problem_solved": cert["universal_folding_problem_solved"],
        "atomistic_folding_solved": cert["atomistic_folding_solved"],
        "all_sequence_prediction_solved": cert["all_sequence_prediction_solved"],
        "engine_revision_after_v58_required": cert["engine_revision_after_v58_required"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] in {STRONG_PASS, BEHAVIOR_ONLY} else 1


if __name__ == "__main__":
    raise SystemExit(main())
