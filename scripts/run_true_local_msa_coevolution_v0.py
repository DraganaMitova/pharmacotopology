#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows
from pharmacotopology.folding_true_local_msa_coevolution import TrueLocalMSAPacket, run_true_local_msa_packet

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "true_local_msa_coevolution_v0"


def _write_csv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _flatten_row(row: dict) -> dict:
    flat = dict(row)
    metric = flat.pop("metric_after_native_audit", {})
    if isinstance(metric, dict):
        for key, value in metric.items():
            flat["metric_%s" % key] = value
    return flat


def _write_packet(packet: TrueLocalMSAPacket, out_dir: Path) -> dict:
    payload = packet.to_dict()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "true_local_msa_coevolution_report.json"
    rows_path = out_dir / "true_local_msa_coevolution_rows.csv"
    contacts_path = out_dir / "true_local_msa_coevolution_contacts.csv"
    decisions_path = out_dir / "true_local_msa_coevolution_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_row(row) for row in payload["rows"]])
    _write_csv(contacts_path, payload["contacts"])
    _write_csv(decisions_path, payload["decisions"])
    return {"report": str(report_path), "rows": str(rows_path), "contacts": str(contacts_path), "decisions": str(decisions_path)}


def _summary(packet: TrueLocalMSAPacket) -> dict:
    return {
        "kind": packet.kind,
        "source_mode": packet.source_mode,
        "row_count": packet.row_count,
        "msa_path": packet.msa_path,
        "raw_msa_available_for_true_local_mi": packet.raw_msa_available_for_true_local_mi,
        "top_safe_dca_anchor_count_requested": packet.top_safe_dca_anchor_count_requested,
        "local_window": packet.local_window,
        "threshold": packet.threshold,
        "min_filtered_sequences": packet.min_filtered_sequences,
        "mean_safe_anchor_count": packet.mean_safe_anchor_count,
        "mean_accepted_local_pair_count": packet.mean_accepted_local_pair_count,
        "mean_native_contact_precision_after_audit": packet.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": packet.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": packet.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": packet.mean_contact_map_f1_after_audit,
        "true_local_msa_claim_allowed": packet.true_local_msa_claim_allowed,
        "universal_physical_law_claim_allowed": packet.universal_physical_law_claim_allowed,
        "folding_problem_solved": packet.folding_problem_solved,
        "claim_rejection_reason": packet.claim_rejection_reason,
    }


def _write_action_required(out_dir: Path, msa_path: str, source_accession: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    text = """# ACTION REQUIRED: raw MSA missing for true local co-evolution

The true local co-evolution experiment cannot run without a raw aligned MSA.
This runner intentionally refuses to fall back to the proxy channel.

Requested source accession: {source_accession}
Requested MSA path: `{msa_path}`

Build the MSA locally, then rerun:

```bash
PYTHONPATH=src bash scripts/build_4ake_msa_with_jackhmmer_v0.sh \
  --database /path/to/uniprot_sprot.fasta \
  --out-dir external_msa/4ake_jackhmmer

PYTHONPATH=src python3 scripts/run_true_local_msa_coevolution_v0.py \
  --source-accession 4AKE:A \
  --msa external_msa/4ake_jackhmmer/4ake_jackhmmer.afa \
  --out-dir first_contact_clean_pharmacotopology_layer_run/true_local_msa_coevolution_v0
```

Expected accepted formats: aligned FASTA/A2M/A3M or Stockholm. The first record
should be the query sequence or one record must ungap-match the 4AKE sequence.
""".format(source_accession=source_accession, msa_path=msa_path)
    (out_dir / "ACTION_REQUIRED_RAW_MSA.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run true raw-MSA local co-evolution expansion around DCA anchors.")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--msa", required=True, help="Raw aligned MSA file: FASTA/A2M/A3M/Stockholm")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--top-anchor-count", type=int, default=50)
    parser.add_argument("--window", type=int, default=5)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--degree-cap", type=int, default=5)
    parser.add_argument("--min-filtered-sequences", type=int, default=16)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    if not Path(args.msa).exists():
        _write_action_required(out_dir, args.msa, args.source_accession)
        summary = {
            "kind": "true_local_msa_coevolution_expansion_v0",
            "status": "action_required_raw_msa_missing",
            "source_accession": args.source_accession,
            "msa_path": args.msa,
            "raw_msa_available_for_true_local_mi": False,
            "proxy_fallback_used": False,
            "out_dir": str(out_dir),
            "action_required_file": str(out_dir / "ACTION_REQUIRED_RAW_MSA.md"),
            "full_suite_run": False,
            "gif_generation_used": False,
            "predictor_subprocesses_used": False,
            "native_truth_used_before_selection": False,
            "coordinate_truth_used_before_selection": False,
        }
        (out_dir / "true_local_msa_coevolution_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        print(json.dumps(summary, indent=2, sort_keys=True))
        return

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    constraints = load_coupling_dataset(Path(args.external_coupling_file)).constraints
    packet = run_true_local_msa_packet(
        rows,
        constraints,
        msa_path=Path(args.msa),
        evaluation_source_accessions=(args.source_accession,),
        top_anchor_count=args.top_anchor_count,
        window=args.window,
        threshold=args.threshold,
        degree_cap=args.degree_cap,
        min_filtered_sequences=args.min_filtered_sequences,
    )
    artifacts = _write_packet(packet, out_dir)
    summary = _summary(packet)
    summary.update({
        "out_dir": str(out_dir),
        "artifacts": artifacts,
        "full_suite_run": False,
        "gif_generation_used": False,
        "predictor_subprocesses_used": False,
        "native_truth_used_before_selection": False,
        "coordinate_truth_used_before_selection": False,
    })
    (out_dir / "true_local_msa_coevolution_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
