#!/usr/bin/env python3
"""Build the final no-AlphaFold 4AKE frontier report.

The report compares three bounded native-free branches:
1. iterative distance-geometry diffusion,
2. leave-one-out empirical sequence-contact prior trained on non-4AKE rows,
3. empirical local contacts fused with DG long-range contacts.

4AKE native coordinates are used only after each branch has selected contacts.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_leave_one_out_empirical_contact_prior import build_leave_one_out_empirical_contact_prior  # noqa: E402
from pharmacotopology.folding_native_contact_eval import ContactPair, evaluate_contact_prediction, normalized_contact_pairs  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows  # noqa: E402

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_DG_SELECTED = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_iterative_distance_geometry_diffusion_v0" / "4ake_iterative_distance_geometry_diffusion_selected_contacts.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_noaf_final_frontier_v0"


def _row_by_accession(rows: Sequence[object], accession: str):
    matches = [row for row in rows if getattr(row, "source_accession") == accession]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {accession}, found {len(matches)}")
    return matches[0]


def _read_pairs_csv(path: Path) -> tuple[ContactPair, ...]:
    if not path.exists():
        return ()
    pairs: list[ContactPair] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if str(row.get("selected", "true")).lower() not in ("true", "1", "yes"):
                continue
            i = row.get("i") or row.get("left") or row.get("residue_i")
            j = row.get("j") or row.get("right") or row.get("residue_j")
            if i is None or j is None:
                continue
            pairs.append((int(i), int(j)))
    return normalized_contact_pairs(pairs)


def _metric_row(name: str, native_pairs: Iterable[ContactPair], predicted_pairs: Iterable[ContactPair]) -> dict[str, object]:
    predicted = tuple(normalized_contact_pairs(predicted_pairs))
    metric = evaluate_contact_prediction(native_pairs=native_pairs, predicted_pairs=predicted)
    return {
        "branch": name,
        "predicted_contact_count": metric.predicted_contact_count,
        "predicted_long_range_contact_count": sum(1 for pair in predicted if pair[1] - pair[0] >= 24),
        **metric.to_dict(),
        "folding_problem_solved_after_native_audit": (
            metric.native_contact_precision >= 0.70
            and metric.native_contact_recall >= 0.70
            and metric.long_range_contact_recall >= 0.70
        ),
    }


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="4AKE final no-AF frontier report")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--dg-selected-contacts", default=str(DEFAULT_DG_SELECTED))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    target = _row_by_accession(rows, args.source_accession)
    training = tuple(row for row in rows if row.source_accession != target.source_accession and row.row_id != target.row_id)
    native_pairs = tuple(target.native_contact_pairs())

    empirical = build_leave_one_out_empirical_contact_prior(
        target_row=target,
        training_rows=training,
        budget_fraction=2.06,
        max_degree=6,
    )
    empirical_pairs = tuple(empirical.selected_pairs)
    dg_pairs = _read_pairs_csv(Path(args.dg_selected_contacts))
    fusion_pairs = normalized_contact_pairs(
        list(empirical_pairs) + [pair for pair in dg_pairs if pair[1] - pair[0] >= 24]
    )

    metric_rows = [
        _metric_row("iterative_dg_diffusion", native_pairs, dg_pairs),
        _metric_row("leave_one_out_empirical_sequence_prior", native_pairs, empirical_pairs),
        _metric_row("empirical_local_plus_dg_long_range_fusion", native_pairs, fusion_pairs),
    ]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(out_dir / "4ake_noaf_final_frontier_metrics.csv", metric_rows)

    solved = any(bool(row["folding_problem_solved_after_native_audit"]) for row in metric_rows)
    best_by_f1 = max(metric_rows, key=lambda row: float(row["contact_map_f1"]))
    best_by_recall = max(metric_rows, key=lambda row: float(row["native_contact_recall"]))
    report = {
        "kind": "4ake_no_alphafold_final_frontier_report_v0",
        "source_accession": target.source_accession,
        "sequence_hash": target.sequence_sha256,
        "sequence_length": target.sequence_length,
        "raw_sequence_persisted": False,
        "gif_generation_used": False,
        "internet_required_to_run_project_scripts": False,
        "alphafold_used_as_evidence": False,
        "msa_used_as_evidence": False,
        "target_native_truth_used_before_selection": False,
        "target_coordinate_truth_used_before_selection": False,
        "training_coordinate_truth_used_from_non_target_rows_in_empirical_branch": True,
        "branches": metric_rows,
        "best_branch_by_f1_after_native_audit": best_by_f1["branch"],
        "best_branch_by_recall_after_native_audit": best_by_recall["branch"],
        "folding_problem_solved_after_native_audit": solved,
        "benchmark_claim_allowed_after_native_audit": solved,
        "claim_rejection_reason": "none" if solved else "no_noaf_branch_reached_0_70_precision_0_70_recall_0_70_long_range_recall",
    }
    (out_dir / "4ake_noaf_final_frontier_report.json").write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    summary = """# 4AKE no-AlphaFold final frontier v0

Final audited branches:

| branch | precision | recall | F1 | long-range recall | solved |
|---|---:|---:|---:|---:|---:|
"""
    for row in metric_rows:
        summary += (
            f"| {row['branch']} | {row['native_contact_precision']:.6f} | {row['native_contact_recall']:.6f} | "
            f"{row['contact_map_f1']:.6f} | {row['long_range_contact_recall']:.6f} | "
            f"{str(row['folding_problem_solved_after_native_audit']).lower()} |\n"
        )
    summary += f"""

Verdict:

```text
folding_problem_solved_after_native_audit = {str(solved).lower()}
claim_rejection_reason = {report['claim_rejection_reason']}
```

Interpretation:

The leave-one-out empirical prior is the strongest no-AF branch by F1 and improves overall contact recovery, but it is local-contact heavy and has zero long-range recall. Iterative DG diffusion supplies some long-range signal, but not enough. The fused branch improves recall over DG alone but still fails the 0.70/0.70/0.70 solve boundary.

So the final no-AF conclusion remains: 4AKE is not cracked by hand physics, DG diffusion, or a small empirical prior. A real MSA-free learned structure predictor or much stronger generative model is still missing.
"""
    (out_dir / "4AKE_NOAF_FINAL_FRONTIER_V0_SUMMARY.md").write_text(summary, encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
