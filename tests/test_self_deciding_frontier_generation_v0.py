from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
DIAGNOSTICS = RUN_DIR / "contact_collapse_diagnostics_v0.json"
GENERATION_ROWS = RUN_DIR / "self_deciding_frontier_generation_rows_v0.csv"


def _report() -> dict:
    return json.loads(DIAGNOSTICS.read_text(encoding="utf-8"))


def test_self_deciding_frontier_generation_rows_are_native_free_and_written() -> None:
    report = _report()
    rows = report["self_deciding_frontier_generation_rows"]
    assert GENERATION_ROWS.exists()
    with GENERATION_ROWS.open(newline="", encoding="utf-8") as handle:
        csv_rows = list(csv.DictReader(handle))
    assert len(csv_rows) == len(rows)
    assert rows
    assert all(row["native_truth_used_before_frontier_generation"] is False for row in rows)
    assert all(row["coordinate_truth_used_before_frontier_generation"] is False for row in rows)
    assert any(
        row["source_accession"] == "4AKE:A"
        and row["self_verified_frontier_generation_prefilter_selected"] is True
        for row in rows
    )


def test_4ake_candidate_pool_generation_is_audited_but_rejected_as_main() -> None:
    report = _report()
    probes = report["hard_target_rescue_probe"]["4AKE:A"]
    expanded = probes["self_deciding_frontier_expansion_merged"]
    generated = probes["self_deciding_frontier_generation_merged"]
    decision = report["frontier_generation_decisions"]["4AKE:A"]

    assert generated["native_truth_used_before_collapse_selection"] is False
    assert generated["coordinate_truth_used_before_collapse_selection"] is False
    assert decision["native_truth_used_before_frontier_generation_decision"] is False
    assert decision["native_truth_attached_after_frontier_generation_for_evaluation"] is True

    assert decision["frontier_generation_probe_accepted_as_main"] is False
    assert decision["frontier_generation_decision"] == "rejected_precision_collapse_or_no_recall_gain"
    assert generated["collapsed_long_range_recall"] <= expanded["collapsed_long_range_recall"]
    assert generated["collapsed_contact_precision"] < expanded["collapsed_contact_precision"]


def test_4ake_frontier_ceiling_audit_exposes_generation_bottleneck() -> None:
    report = _report()
    audit = report["frontier_ceiling_audit"]["4AKE:A"]

    assert audit["candidate_generator_long_range_recall_ceiling"] == 1.0
    assert audit["competitive_region_long_range_recall_ceiling"] < 0.40
    assert audit["frontier_generation_bottleneck_for_0_40_recall"] is True
    assert audit["generated_frontier_event_count"] >= audit["seed_frontier_event_count"]
    assert audit["generated_region_long_range_recall_ceiling"] <= audit["candidate_generator_long_range_recall_ceiling"]
