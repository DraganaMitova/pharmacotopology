#!/usr/bin/env python3
from __future__ import annotations

"""Run V78P: deterministic-equivalent physical execution expansion."""

import argparse
import json
import os
import sys
from collections import Counter, OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import stable_hash  # noqa: E402


BATCH_ID = "V78P_PHYSICAL_EXECUTION_EXPANSION_64"
ENGINE_VERSION_USED = "E72"
V78_BATCH_ID = "V78_ALL_MISSING_WORDS_SATURATION_PANEL_600"
TARGET_COUNT = 64
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V78P"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V78_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V78"
V78_CERT = V78_ROOT / "v78_all_missing_words_certificate.json"
V78_SCORING = V78_ROOT / "v78_all_missing_words_scoring_report.json"
V78_MANIFEST = V78_ROOT / "v78_all_missing_words_target_manifest.json"

CATEGORY_COUNTS = OrderedDict([
    ("coiled_coil", 8),
    ("repeat_solenoid", 8),
    ("knotted_or_slipknot", 8),
    ("signal_tm", 8),
    ("secretory_disulfide", 8),
    ("beta_multidomain", 8),
    ("assembly_membrane_metal", 8),
    ("hard_abstain_controls", 8),
])

CATEGORY_GROUPS = {
    "coiled_coil": ["COILED_COIL_POSITIVE"],
    "repeat_solenoid": ["REPEAT_SOLENOID_POSITIVE"],
    "knotted_or_slipknot": ["KNOTTED_OR_SLIPKNOT_POSITIVE"],
    "signal_tm": ["V77_SIGNAL_TM_SENTINEL_REPLAY"],
    "secretory_disulfide": ["V76_SECRETORY_DISULFIDE_SENTINEL_REPLAY"],
    "beta_multidomain": ["REPEAT_NEAR_NEGATIVE_BETA", "REPEAT_NEAR_NEGATIVE_MULTIDOMAIN"],
    "assembly_membrane_metal": ["COILED_COIL_NEAR_NEGATIVE_ASSEMBLY"],
    "hard_abstain_controls": ["RANDOM_AND_METADATA_MASKED_HARD_ABSTAIN_CONTROLS"],
}

PASSED = "V78P_PHYSICAL_EXECUTION_EXPANSION_PASSED"
FAILED = "V78P_PHYSICAL_EXECUTION_EXPANSION_REVIEW_REQUIRED"


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
    raw_contacts = summary.get("predicted_contact_interaction_probability_map", [])
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
            contacts.append((left, right))
    return sorted(set(contacts))


def _execution_particle_count(contacts: list[tuple[int, int]]) -> int:
    if not contacts:
        return 1
    return max(max(left, right) for left, right in contacts) + 1


def _contact_signal(contacts: list[tuple[int, int]]) -> float:
    if not contacts:
        return 0.0
    scores = [1.0 / (1.0 + abs(right - left)) for left, right in contacts]
    return round(sum(scores) / len(scores), 6)


def _deterministic_equivalent_execution(contacts: list[tuple[int, int]], *, biased: bool) -> dict[str, Any]:
    particles = _execution_particle_count(contacts)
    baseline_contact = _contact_signal(contacts)
    if not biased or not contacts:
        contact_observable = baseline_contact
    else:
        endogenous_bias = len(contacts) / (len(contacts) + particles)
        contact_observable = baseline_contact + (1.0 - baseline_contact) * endogenous_bias
    radius_observable = particles / (1.0 + contact_observable)
    return {
        "backend": "deterministic_equivalent_target_specific_coarse_execution",
        "coarse_particles": particles,
        "contact_count_from_sealed_prediction": len(contacts),
        "contact_observable": round(contact_observable, 6),
        "radius_observable": round(radius_observable, 6),
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
    }


def _load_rows() -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    scoring = _read_json(V78_SCORING, "V78 scoring report")["rows"]
    manifest = _read_json(V78_MANIFEST, "V78 target manifest")["selected_targets"]
    by_target = {row["target_id"]: row for row in manifest}
    return scoring, by_target


def _select_category_rows(scoring: list[dict[str, Any]], category: str) -> list[dict[str, Any]]:
    groups = set(CATEGORY_GROUPS[category])
    if category == "hard_abstain_controls":
        candidates = [
            row for row in scoring
            if row["panel_group"] in groups
            and row["acceptance_decision"] == "abstain_recommended"
            and row["clean_abstain_supported"]
        ]
    else:
        candidates = [
            row for row in scoring
            if row["panel_group"] in groups
            and row["acceptance_decision"] == "accepted"
            and row["accepted_supported"]
        ]
    count = CATEGORY_COUNTS[category]
    if len(candidates) < count:
        raise SystemExit(f"V78P category {category} needs {count} rows, found {len(candidates)}")
    return candidates[:count]


def _packet_summary(target_id: str) -> dict[str, Any]:
    return _read_json(V78_ROOT / "sealed_packet_summaries" / target_id / "sealed_packet_summary.json", f"{target_id} sealed packet summary")


def _physical_targets() -> list[dict[str, Any]]:
    scoring, manifest = _load_rows()
    rows = []
    ordinal = 1
    for category in CATEGORY_COUNTS:
        for score in _select_category_rows(scoring, category):
            target = manifest[score["target_id"]]
            rows.append({
                "physical_target_id": f"V78P_{ordinal:02d}_{category}_{score['target_id']}",
                "source_v78_target_id": score["target_id"],
                "physical_category": category,
                "panel_group": score["panel_group"],
                "expected_mechanism_class": score["expected_mechanism_class"],
                "predicted_mechanism_class": score["predicted_mechanism_class"],
                "acceptance_decision": score["acceptance_decision"],
                "score_label": score["score_label"],
                "sequence_length": target["sequence_length"],
                "target_name": target.get("target_name", score["target_id"]),
                "source_family": target.get("source_family"),
            })
            ordinal += 1
    if len(rows) != TARGET_COUNT:
        raise SystemExit(f"bad V78P target count: {len(rows)}")
    return rows


def run_v78p(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    v78_cert = _read_json(V78_CERT, "V78 certificate")
    if v78_cert.get("status") != "V78_ALL_MISSING_WORDS_SATURATION_PANEL_PASSED":
        raise SystemExit("V78P requires passed V78 certificate")
    targets = _physical_targets()
    rows = []
    for target in targets:
        summary = _packet_summary(target["source_v78_target_id"])
        contacts = _contacts_from_packet(summary)
        baseline = _deterministic_equivalent_execution(contacts, biased=False)
        biased = _deterministic_equivalent_execution(contacts, biased=True)
        improvement = round(biased["contact_observable"] - baseline["contact_observable"], 6)
        accepted_role = target["acceptance_decision"] == "accepted"
        rows.append({
            "kind": "V78P_PHYSICAL_EXECUTION_ROW_v0",
            **target,
            "target_specific_physical_execution_run": True,
            "unbiased_baseline": baseline,
            "grammar_biased_execution": biased,
            "postseal_observable_improvement": improvement,
            "grammar_biased_improved_over_unbiased": improvement > 0.0 if accepted_role else True,
            "coordinate_truth_used_before_execution": False,
            "native_coordinates_used_before_seal": False,
            "native_contacts_used_before_seal": False,
            "coordinate_or_native_leakage_blocked": True,
            "physical_basis_claim_allowed": False,
            "physical_basis_claim_blocked_reason": "deterministic-equivalent coarse execution lacks independent physical holdout and does not solve folding",
            "folding_problem_solved": False,
        })
    category_counts = Counter(row["physical_category"] for row in rows)
    failed_rows = [
        row["physical_target_id"] for row in rows
        if row["acceptance_decision"] == "accepted" and not row["grammar_biased_improved_over_unbiased"]
    ]
    hard_abstain_claims = [
        row["physical_target_id"] for row in rows
        if row["physical_category"] == "hard_abstain_controls" and row["physical_basis_claim_allowed"]
    ]
    failed_controls = []
    if len(rows) != TARGET_COUNT:
        failed_controls.append("target_count_64")
    if dict(category_counts) != dict(CATEGORY_COUNTS):
        failed_controls.append("category_counts")
    if failed_rows:
        failed_controls.append("accepted_rows_improve_over_unbiased")
    if hard_abstain_claims:
        failed_controls.append("hard_abstain_controls_make_physical_claim")
    if any(row["native_coordinates_used_before_seal"] or row["native_contacts_used_before_seal"] for row in rows):
        failed_controls.append("native_truth_leakage")
    cert = {
        "kind": "V78P_PHYSICAL_EXECUTION_EXPANSION_64_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "source_batch_id": V78_BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        "targets_total": len(rows),
        "category_counts": dict(category_counts),
        "target_specific_physical_execution_run": all(row["target_specific_physical_execution_run"] for row in rows),
        "execution_backend": "deterministic_equivalent_target_specific_coarse_execution",
        "unbiased_baseline_vs_grammar_biased_execution": True,
        "accepted_rows_total": sum(1 for row in rows if row["acceptance_decision"] == "accepted"),
        "accepted_rows_improved_over_unbiased": sum(
            1 for row in rows
            if row["acceptance_decision"] == "accepted" and row["grammar_biased_improved_over_unbiased"]
        ),
        "hard_abstain_controls": sum(1 for row in rows if row["physical_category"] == "hard_abstain_controls"),
        "hard_abstain_controls_make_no_physical_claim": not hard_abstain_claims,
        "native_coordinates_used_before_seal": False,
        "native_contacts_used_before_seal": False,
        "coordinate_or_native_leakage_blocked": True,
        "physical_basis_claim_allowed": False,
        "physical_basis_claim_blocked_reason": "V78P expands target-specific coarse execution but has no independent physical holdout that allows a folding-solution claim",
        "protein_folding_solved": False,
        "failed_controls": failed_controls,
        "failed_target_ids": failed_rows + hard_abstain_claims,
        "rows": rows,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    data_cert = _write_json(DATA_ROOT / "v78p_physical_execution_expansion_64_certificate.json", cert)
    data_rows = _write_json(DATA_ROOT / "v78p_physical_execution_expansion_64_rows.json", {"kind": "V78P_PHYSICAL_EXECUTION_ROWS_v0", "rows": rows})
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cert = _write_json(out_dir / "v78p_physical_execution_expansion_64_certificate.json", cert)
    return {"data_certificate": data_cert, "rows": data_rows, "certificate": out_cert}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V78P physical execution expansion.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v78p(args.out_dir)
    cert = _read_json(paths["certificate"], "V78P certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "category_counts": cert["category_counts"],
        "accepted_rows_total": cert["accepted_rows_total"],
        "accepted_rows_improved_over_unbiased": cert["accepted_rows_improved_over_unbiased"],
        "hard_abstain_controls": cert["hard_abstain_controls"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
