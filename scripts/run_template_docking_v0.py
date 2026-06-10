#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_template_docking import (
    TEMPLATE_DOCKING_KIND,
    collect_pdb_templates,
    _is_kinase_like,
    parse_benchmark_templates_with_filters,
    run_template_docking_v0,
)

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "template_docking_v0"
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            serializable = {}
            for key, value in row.items():
                if isinstance(value, (tuple, list)):
                    serializable[key] = ",".join(str(v) for v in value)
                else:
                    serializable[key] = value
            writer.writerow(serializable)


def _flatten_report(report: dict[str, object]) -> dict[str, object]:
    payload = dict(report)
    metric = payload.pop("metric_after_native_audit", {})
    if isinstance(metric, dict):
        for key, value in metric.items():
            payload[f"metric_{key}"] = value
    payload["template_rows"] = payload.get("template_rows", [])
    payload["voted_pairs"] = payload.get("voted_pairs", [])
    payload["winner_pairs"] = payload.get("winner_pairs", [])
    payload["predicted_pairs"] = payload.get("predicted_pairs", [])
    return payload


def _read_accession_lines(path: str | None) -> set[str]:
    if not path:
        return set()
    requested = Path(path)
    if not requested.exists():
        raise SystemExit(f"Accession file does not exist: {requested}")
    values = set()
    for line in requested.read_text(encoding="utf-8").splitlines():
        cleaned = line.strip()
        if not cleaned or cleaned.startswith("#"):
            continue
        for token in cleaned.replace(",", " ").split():
            if token:
                values.add(token.strip())
    return values


def _to_csv_rows(values: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for raw in values:
        item = dict(raw)
        if isinstance(item.get("template_ids"), tuple):
            item["template_ids"] = ",".join(item["template_ids"])
        out.append(item)
    return out


def main() -> None:
    def _parse_chemical_threshold(raw: str) -> float | None:
        lower = raw.strip().lower()
        if lower in {"auto", "a", ""}:
            return None
        return float(raw)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--template-dir", default=None)
    parser.add_argument("--template-chain", default=None)
    parser.add_argument("--include-benchmark-templates", action="store_true", default=True)
    parser.add_argument("--no-benchmark-templates", dest="include_benchmark_templates", action="store_false")
    parser.add_argument("--max-template-count", type=int, default=None)
    parser.add_argument("--contact-cutoff", type=float, default=8.0)
    parser.add_argument(
        "--chemical-threshold",
        type=_parse_chemical_threshold,
        default=None,
        help="float threshold or 'auto' (default)",
    )
    parser.add_argument("--no-auto-kinase-filter", action="store_true", default=False)
    parser.add_argument("--match-target-architecture", action="store_true", default=False)
    parser.add_argument("--only-kinase-templates", action="store_true", default=False)
    parser.add_argument("--template-allowlist", default=None, help="Template accession allowlist file")
    parser.add_argument("--template-exclude", default=None, help="Template accession exclude list file")
    parser.add_argument(
        "--domain-scan-step",
        type=int,
        default=1,
        help="Domain-scan step for faster domain-aware alignment (>=1)",
    )
    parser.add_argument(
        "--domain-scan-window",
        type=int,
        default=0,
        help="Window around midpoint for domain-aware scan (0 = full scan)",
    )
    parser.add_argument("--jobs", type=int, default=1, help="Number of parallel workers (0 = all CPUs)")
    parser.add_argument("--min-template-support", type=int, default=1)
    parser.add_argument("--long-range-threshold", type=int, default=24)
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    target_row = None
    for row in rows:
        if row.source_accession == args.source_accession:
            target_row = row
            break
    if target_row is None:
        raise SystemExit(f"No target row for source-accession={args.source_accession}")

    target_reference_fold_class = target_row.reference_fold_class
    target_architecture_axis = target_row.truth_axes.get("architecture_axis", "unknown")
    target_is_kinase_like = (
        _is_kinase_like(target_row.source_id)
        or _is_kinase_like(target_row.source_accession)
        or _is_kinase_like(target_row.reference_fold_class)
    )
    effective_only_kinase_templates = (
        args.only_kinase_templates
        or (not args.no_auto_kinase_filter and target_is_kinase_like)
    )
    match_target_architecture = (
        args.match_target_architecture
        or (effective_only_kinase_templates and not args.no_auto_kinase_filter)
    )

    templates = []
    if args.include_benchmark_templates:
        templates.extend(
            parse_benchmark_templates_with_filters(
                rows,
                args.source_accession,
                only_kinase_templates=effective_only_kinase_templates,
                allowed_source_accessions=_read_accession_lines(args.template_allowlist),
                excluded_source_accessions=_read_accession_lines(args.template_exclude),
                target_reference_fold_class=target_reference_fold_class,
                target_architecture_axis=target_architecture_axis,
                restrict_same_architecture=match_target_architecture,
            )
        )
        if effective_only_kinase_templates and args.max_template_count is None:
            templates = list(templates)[:50]

    template_paths = []
    if args.template_dir:
        template_dir = Path(args.template_dir)
        if not template_dir.exists():
            raise SystemExit(f"Template directory does not exist: {template_dir}")
        template_paths = [
            p
            for p in sorted(template_dir.iterdir())
            if p.is_file() and p.suffix.lower() == ".pdb"
        ]
    if template_paths:
        templates.extend(
            collect_pdb_templates(
                template_paths,
                chain_id=args.template_chain,
                max_template_count=args.max_template_count,
            )
        )

    source_mode_parts = []
    if args.include_benchmark_templates:
        source_mode_parts.append("benchmark")
    if template_paths:
        source_mode_parts.append("pdb")
    if not source_mode_parts:
        source_mode_parts.append("none")

    if args.max_template_count is not None:
        templates = list(templates)[: args.max_template_count]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report, constraint_payload, _ = run_template_docking_v0(
        target_row=target_row,
        templates=tuple(templates),
        long_range_threshold=args.long_range_threshold,
        chemical_threshold=args.chemical_threshold,
        contact_cutoff_angstrom=args.contact_cutoff,
        minimum_sequence_separation=3,
        min_template_support=args.min_template_support,
        domain_scan_step=args.domain_scan_step,
        domain_scan_window=args.domain_scan_window,
        jobs=args.jobs,
        source_mode="+".join(source_mode_parts),
    )

    report_payload = report.to_dict()
    report_path = out_dir / "template_docking_v0_report.json"
    template_rows_path = out_dir / "template_docking_v0_template_rows.csv"
    voted_pairs_path = out_dir / "template_docking_v0_voted_pairs.csv"
    winner_pairs_path = out_dir / "template_docking_v0_winner_pairs.csv"
    predicted_pairs_path = out_dir / "template_docking_v0_predicted_pairs.csv"
    constraints_path = out_dir / "template_docking_v0_constraints.json"
    certificate_path = out_dir / "template_docking_v0_certificate.json"

    report_path.write_text(json.dumps(_flatten_report(report_payload), indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(template_rows_path, report_payload["template_rows"])  # type: ignore[arg-type]
    _write_csv(voted_pairs_path, _to_csv_rows(list(report_payload["voted_pairs"])))  # type: ignore[arg-type]
    _write_csv(winner_pairs_path, _to_csv_rows(list(report_payload["winner_pairs"])))  # type: ignore[arg-type]
    _write_csv(predicted_pairs_path, _to_csv_rows(list(report_payload["predicted_pairs"])))  # type: ignore[arg-type]
    constraints_path.write_text(json.dumps(constraint_payload, indent=2, sort_keys=True), encoding="utf-8")

    summary = {
        "kind": TEMPLATE_DOCKING_KIND,
        "source_accession": target_row.source_accession,
        "out_dir": str(out_dir),
        "template_count": report.template_count,
        "templates_used": report.templates_used,
        "templates_with_mapped_pairs": report.templates_with_mapped_pairs,
        "winner_template_id": report.winner_template_id,
        "winner_template_source": report.winner_template_source,
        "predicted_pair_count": report.predicted_pair_count,
        "native_contact_precision": report.native_contact_precision,
        "native_contact_recall": report.native_contact_recall,
        "long_range_contact_recall": report.long_range_recall,
        "contact_map_f1": report.contact_map_f1,
        "false_contact_rate": report.false_contact_rate,
        "artifacts": {
            "report": str(report_path),
            "template_rows": str(template_rows_path),
            "voted_pairs": str(voted_pairs_path),
            "winner_pairs": str(winner_pairs_path),
            "predicted_pairs": str(predicted_pairs_path),
            "constraints": str(constraints_path),
        },
    }
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(certificate_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
