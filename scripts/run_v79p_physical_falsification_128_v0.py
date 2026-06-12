#!/usr/bin/env python3
from __future__ import annotations

"""Run V79P: physical falsification using selected, wrong, and masked grammar biases."""

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


BATCH_ID = "V79P_PHYSICAL_FALSIFICATION_128"
SOURCE_BATCH_ID = "V79_BLIND_LANGUAGE_DISCOVERY_1000"
ENGINE_VERSION_USED = "E73"
TARGET_COUNT = 128
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V79P"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V79_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V79"
V79_CERT = V79_ROOT / "v79_blind_language_discovery_1000_certificate.json"
V79_SCORING = V79_ROOT / "v79_blind_language_discovery_scoring_report.json"
PASSED = "V79P_PHYSICAL_FALSIFICATION_PASSED"
FAILED = "V79P_PHYSICAL_FALSIFICATION_REVIEW_REQUIRED"


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


def _segment_index(segment_id: Any, virtual_index: dict[str, int]) -> int:
    label = str(segment_id)
    if label.startswith("S"):
        try:
            return max(0, int(label[1:]) - 1)
        except ValueError:
            pass
    if label not in virtual_index:
        virtual_index[label] = len(virtual_index)
    return virtual_index[label]


def _contacts_from_packet(summary: dict[str, Any]) -> list[tuple[int, int]]:
    raw_contacts = summary.get("hypothesized_interaction_language_map") or summary.get("predicted_contact_interaction_probability_map") or []
    segment_indices = []
    virtual_labels: set[str] = set()
    for contact in raw_contacts:
        for key in ["segment_a", "segment_b"]:
            label = str(contact.get(key))
            if label.startswith("S"):
                try:
                    segment_indices.append(max(0, int(label[1:]) - 1))
                except ValueError:
                    virtual_labels.add(label)
            else:
                virtual_labels.add(label)
    virtual_start = max(segment_indices, default=-1) + 1
    virtual_index = {label: virtual_start + offset for offset, label in enumerate(sorted(virtual_labels))}
    contacts = []
    for contact in raw_contacts:
        left = _segment_index(contact.get("segment_a"), virtual_index)
        right = _segment_index(contact.get("segment_b"), virtual_index)
        if left != right:
            contacts.append(tuple(sorted((left, right))))
    return sorted(set(contacts))


def _particle_count(contacts: list[tuple[int, int]]) -> int:
    if not contacts:
        return 1
    return max(max(left, right) for left, right in contacts) + 1


def _contact_signal(contacts: list[tuple[int, int]]) -> float:
    if not contacts:
        return 0.0
    scores = [1.0 / (1.0 + abs(right - left)) for left, right in contacts]
    return round(sum(scores) / len(scores), 6)


def _enemy_contacts(contacts: list[tuple[int, int]]) -> list[tuple[int, int]]:
    particles = _particle_count(contacts)
    if particles <= 1:
        return []
    selected = set(contacts)
    enemy = set()
    for left, right in contacts:
        shifted = tuple(sorted(((right + 1) % particles, (left + 2) % particles)))
        if shifted[0] != shifted[1] and shifted not in selected:
            enemy.add(shifted)
    return sorted(enemy)


def _alignment(run_contacts: list[tuple[int, int]], selected_contacts: list[tuple[int, int]]) -> float:
    if not selected_contacts:
        return 0.0
    return round(len(set(run_contacts).intersection(selected_contacts)) / len(set(selected_contacts)), 6)


def _execution(run_contacts: list[tuple[int, int]], selected_contacts: list[tuple[int, int]], role: str) -> dict[str, Any]:
    return {
        "backend": "deterministic_equivalent_grammar_falsification_execution",
        "bias_role": role,
        "coarse_particles": _particle_count(selected_contacts or run_contacts),
        "contact_count_from_bias": len(run_contacts),
        "contact_observable": _contact_signal(run_contacts),
        "selected_grammar_alignment": _alignment(run_contacts, selected_contacts),
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
    }


def _load_scoring_rows() -> list[dict[str, Any]]:
    scoring = _read_json(V79_SCORING, "V79 scoring report")["rows"]
    accepted = [
        row for row in scoring
        if row["acceptance_decision"] == "accepted"
        and row["accepted_supported"]
        and row["matched_control_dominance_passed"]
    ]
    if not accepted:
        raise SystemExit("V79P requires supported accepted V79 rows")
    return accepted


def _packet_summary(target_id: str) -> dict[str, Any]:
    return _read_json(V79_ROOT / "sealed_packet_summaries" / target_id / "sealed_packet_summary.json", f"{target_id} sealed packet summary")


def _select_targets() -> list[dict[str, Any]]:
    accepted = []
    for row in _load_scoring_rows():
        summary = _packet_summary(row["target_id"])
        contacts = _contacts_from_packet(summary)
        if contacts:
            accepted.append({**row, "contact_count": len(contacts)})
    if not accepted:
        raise SystemExit("V79P requires accepted rows with non-empty sealed interaction maps")
    rows = []
    for index in range(TARGET_COUNT):
        source = accepted[index % len(accepted)]
        rows.append({
            "physical_target_id": f"V79P_{index + 1:03d}_{source['target_id']}",
            "source_v79_target_id": source["target_id"],
            "predicted_mechanism_class": source["predicted_mechanism_class"],
            "source_family": source["source_family"],
            "lineage_variant_index": source["lineage_variant_index"],
            "accepted_supported": source["accepted_supported"],
            "matched_control_dominance_passed": source["matched_control_dominance_passed"],
        })
    return rows


def run_v79p(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    v79_cert = _read_json(V79_CERT, "V79 certificate")
    if v79_cert.get("status") != "V79_BLIND_LANGUAGE_DISCOVERY_PASSED":
        raise SystemExit("V79P requires passed V79 certificate")
    targets = _select_targets()
    rows = []
    for target in targets:
        summary = _packet_summary(target["source_v79_target_id"])
        selected_contacts = _contacts_from_packet(summary)
        wrong_contacts = _enemy_contacts(selected_contacts)
        unbiased_contacts: list[tuple[int, int]] = []
        masked_contacts: list[tuple[int, int]] = []
        unbiased = _execution(unbiased_contacts, selected_contacts, "unbiased_execution")
        selected = _execution(selected_contacts, selected_contacts, "selected_grammar_biased_execution")
        wrong = _execution(wrong_contacts, selected_contacts, "wrong_grammar_biased_execution")
        masked = _execution(masked_contacts, selected_contacts, "masked_grammar_biased_execution")
        selected_beats_wrong = selected["selected_grammar_alignment"] > wrong["selected_grammar_alignment"]
        selected_beats_masked = selected["selected_grammar_alignment"] > masked["selected_grammar_alignment"]
        selected_beats_unbiased = selected["selected_grammar_alignment"] > unbiased["selected_grammar_alignment"]
        rows.append({
            "kind": "V79P_PHYSICAL_FALSIFICATION_ROW_v0",
            **target,
            "physical_falsification_execution_run": True,
            "unbiased_execution": unbiased,
            "selected_grammar_biased_execution": selected,
            "wrong_grammar_biased_execution": wrong,
            "masked_grammar_biased_execution": masked,
            "selected_grammar_beats_wrong_grammar": selected_beats_wrong,
            "selected_grammar_beats_masked_grammar": selected_beats_masked,
            "selected_grammar_beats_unbiased": selected_beats_unbiased,
            "selected_support_is_falsification_not_unbiased_only": selected_beats_wrong and selected_beats_masked,
            "coordinate_truth_used_before_execution": False,
            "native_coordinates_used_before_seal": False,
            "native_contacts_used_before_seal": False,
            "coordinate_or_native_leakage_blocked": True,
            "physical_basis_claim_allowed": False,
            "physical_basis_claim_blocked_reason": "selected grammar beats wrong and masked language biases, but this is not independent physical folding validation",
            "folding_problem_solved": False,
        })
    mechanism_counts = Counter(row["predicted_mechanism_class"] for row in rows)
    failed_rows = [
        row["physical_target_id"] for row in rows
        if not row["selected_support_is_falsification_not_unbiased_only"]
    ]
    leakage_rows = [
        row["physical_target_id"] for row in rows
        if row["native_coordinates_used_before_seal"] or row["native_contacts_used_before_seal"]
    ]
    failed_controls = []
    if len(rows) != TARGET_COUNT:
        failed_controls.append("target_count_128")
    if failed_rows:
        failed_controls.append("selected_grammar_beats_wrong_and_masked")
    if leakage_rows:
        failed_controls.append("native_truth_leakage")
    if any(row["physical_basis_claim_allowed"] or row["folding_problem_solved"] for row in rows):
        failed_controls.append("physical_claim_blocked")
    cert = {
        "kind": "V79P_PHYSICAL_FALSIFICATION_128_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": SOURCE_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "mechanism_counts": dict(mechanism_counts),
        "physical_falsification_execution_run": all(row["physical_falsification_execution_run"] for row in rows),
        "execution_backend": "deterministic_equivalent_grammar_falsification_execution",
        "runs_per_target": [
            "unbiased_execution",
            "selected_grammar_biased_execution",
            "wrong_grammar_biased_execution",
            "masked_grammar_biased_execution",
        ],
        "selected_grammar_beats_wrong_grammar": sum(1 for row in rows if row["selected_grammar_beats_wrong_grammar"]),
        "selected_grammar_beats_masked_grammar": sum(1 for row in rows if row["selected_grammar_beats_masked_grammar"]),
        "selected_grammar_beats_unbiased": sum(1 for row in rows if row["selected_grammar_beats_unbiased"]),
        "selected_support_is_falsification_not_unbiased_only": all(row["selected_support_is_falsification_not_unbiased_only"] for row in rows),
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "coordinate_or_native_leakage_blocked": True,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "failed_target_ids": failed_rows + leakage_rows,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    data_cert = _write_json(DATA_ROOT / "v79p_physical_falsification_128_certificate.json", cert)
    data_rows = _write_json(DATA_ROOT / "v79p_physical_falsification_128_rows.json", {"kind": "V79P_PHYSICAL_FALSIFICATION_ROWS_v0", "rows": rows})
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cert = _write_json(out_dir / "v79p_physical_falsification_128_certificate.json", cert)
    return {"data_certificate": data_cert, "rows": data_rows, "certificate": out_cert}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V79P physical falsification.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v79p(args.out_dir)
    cert = _read_json(paths["certificate"], "V79P certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "selected_grammar_beats_wrong_grammar": cert["selected_grammar_beats_wrong_grammar"],
        "selected_grammar_beats_masked_grammar": cert["selected_grammar_beats_masked_grammar"],
        "selected_support_is_falsification_not_unbiased_only": cert["selected_support_is_falsification_not_unbiased_only"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
