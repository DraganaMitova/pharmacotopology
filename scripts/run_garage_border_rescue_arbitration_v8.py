#!/usr/bin/env python3
from __future__ import annotations

"""Postprocess-only V8 border-rescue arbitration on existing V5 trajectories."""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Optional, Sequence

from run_openmm_tmd_replicas_v0 import (  # type: ignore[import-not-found]
    DEFAULT_BENCHMARK_FILE,
    DEFAULT_EXTERNAL_COUPLING_FILE,
    _load_row,
    _vote_pairs,
)
from pharmacotopology.folding_template_docking import chemical_score

from run_garage_rescue_role_selector_v7 import (  # type: ignore[import-not-found]
    CORE_LONG_RANGE_SEPARATION,
    _aggregate_pair_status,
    _build_stable_pairs_with_lane_min_separation,
    _classify_contact_role,
    _coerce_pair,
    _collect_replica_lanes,
    _collect_replica_map,
    _collect_replica_reachability,
    _extract_target_pairs,
    _read_json,
    _summarize_pair_reachability,
)

Pair = tuple[int, int]


def _safe_ratio(count: int, total: int) -> float:
    if total <= 0:
        return 0.0
    return round(count / total, 6)


def _resolve_dca_thresholds_after(first_metadata: dict[str, Any]) -> tuple[float, dict[str, float]]:
    lane_thresholds = first_metadata.get("resolved_dca_threshold_by_lane")
    if isinstance(lane_thresholds, dict):
        return lane_thresholds.get("balanced_rescue", 0.0), {
            str(k): float(v) if isinstance(v, (int, float)) else float(v or 0.0)
            for k, v in lane_thresholds.items()
        }
    return 0.0, {}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run GARAGE_BORDER_RESCUE_ARBITRATION_V8 postprocess on existing V5 trajectories.",
    )
    parser.add_argument(
        "--source-run-dir",
        default=str(
            Path(__file__).resolve().parents[1]
            / "first_contact_clean_pharmacotopology_layer_run"
            / "GARAGE_TAIL_REACHABILITY_WITH_STRICT_GUARD_V5"
        ),
    )
    parser.add_argument(
        "--source-certificate",
        default="",
        help="Optional explicit source certificate JSON path.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(
            Path(__file__).resolve().parents[1]
            / "first_contact_clean_pharmacotopology_layer_run"
            / "GARAGE_BORDER_RESCUE_ARBITRATION_V8"
        ),
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--readout-coupling-file", default="")
    parser.add_argument("--anchor-profile-file", default="")
    parser.add_argument("--audit-pairs-json", default="")

    parser.add_argument("--strict-min-separation", type=int, default=24)
    parser.add_argument("--balanced-min-separation", type=int, default=24)
    parser.add_argument("--rescue-min-separation", type=int, default=12)
    parser.add_argument("--monitor-min-separation", type=int, default=24)
    parser.add_argument("--long-range-min-separation", type=int, default=CORE_LONG_RANGE_SEPARATION)
    parser.add_argument(
        "--strict-seq-boundary",
        action="store_true",
        default=False,
        help="Use strict (<=) sequence-separation min-separation boundary checks.",
    )

    parser.add_argument(
        "--rescue-late-vote-threshold",
        type=int,
        default=2,
        help="Minimum replica support under rescue-late rule.",
    )
    parser.add_argument(
        "--rescue-tail-frequency-threshold",
        type=float,
        default=0.5,
        help="Tail frequency threshold for rescue-late arbitration candidates.",
    )
    parser.add_argument("--write-json", default="")

    args = parser.parse_args()

    source_run = Path(args.source_run_dir)
    if not source_run.exists():
        raise SystemExit(f"source run directory missing: {source_run}")

    if args.source_certificate:
        certificate_path = Path(args.source_certificate).expanduser()
        if not certificate_path.exists():
            certificate_path = source_run / "openmm_tmd_replicas_v0_certificate.json"
    else:
        certificate_path = source_run / "openmm_tmd_replicas_v0_certificate.json"
    source_payload = _read_json(certificate_path)

    source_accession = source_payload.get("source_accession", "4AKE:A")
    row = _load_row(Path(args.benchmark_file), source_accession)

    topology_mode = str(source_payload.get("topology_mode", "interdomain"))
    domain_boundaries = str(source_payload.get("domain_boundaries", "1-60,61-110,111-170,171-214"))
    tail_fraction = float(source_payload.get("tail_fraction", 0.20))
    contact_cutoff = float(source_payload.get("contact_cutoff_angstrom", 7.0))

    readout_coupling = Path(args.readout_coupling_file) if args.readout_coupling_file else Path(args.external_coupling_file)
    if not readout_coupling.exists():
        readout_coupling = Path(args.external_coupling_file)
    anchor_profile = Path(args.anchor_profile_file) if args.anchor_profile_file else None

    dca_threshold = source_payload.get("dca_threshold", "auto")
    frequency_threshold = source_payload.get("frequency_threshold", "auto")
    strict_dca_threshold = source_payload.get("strict_dca_threshold")
    balanced_strong_dca_threshold = source_payload.get("balanced_strong_dca_threshold", 0.80)
    balanced_rescue_dca_threshold = source_payload.get("balanced_rescue_dca_threshold", 0.70)
    monitor_dca_threshold = source_payload.get("monitor_dca_threshold")

    chemical_threshold = float(source_payload.get("chemical_threshold", 0.5))
    strict_chemical_threshold = source_payload.get("strict_chemical_threshold")
    balanced_strong_chemical_threshold = source_payload.get("balanced_strong_chemical_threshold")
    balanced_rescue_chemical_threshold = source_payload.get("balanced_rescue_chemical_threshold")
    monitor_chemical_threshold = source_payload.get("monitor_chemical_threshold")

    max_degree = int(source_payload.get("max_degree", 4))
    control_count = int(source_payload.get("control_count", 0))
    max_control_overlap_fraction = float(source_payload.get("max_control_overlap_fraction", 0.25))

    vote_threshold = int(source_payload.get("vote_threshold", 3))
    vote_min_average_chem = float(source_payload.get("vote_min_average_chemical", 0.5))
    strict_vote_threshold = source_payload.get("strict_vote_threshold")
    balanced_vote_threshold = source_payload.get("balanced_vote_threshold")
    rescue_vote_threshold = source_payload.get("balanced_rescue_vote_threshold")
    monitor_vote_threshold = source_payload.get("monitor_vote_threshold")
    strict_vote_min_average_chem = source_payload.get("strict_vote_min_average_chemical")
    balanced_vote_min_average_chem = source_payload.get("balanced_vote_min_average_chemical")
    rescue_vote_min_average_chem = source_payload.get("balanced_rescue_vote_min_average_chemical")
    monitor_vote_min_average_chem = source_payload.get("monitor_vote_min_average_chemical")

    target_pairs = _extract_target_pairs(
        source_payload,
        Path(args.audit_pairs_json) if args.audit_pairs_json else None,
    )
    if args.audit_pairs_json:
        explicit_pairs = _extract_target_pairs(
            source_payload,
            Path(args.audit_pairs_json),
        )
    else:
        explicit_pairs = []

    default_border_pairs = [
        _coerce_pair((93, 117)),
        _coerce_pair((96, 120)),
    ]

    if not target_pairs:
        raise SystemExit("target rescue pairs are empty; provide --audit-pairs-json or source pair_audit")

    if args.audit_pairs_json:
        target_pairs_set = set(explicit_pairs)
    else:
        target_pairs_set = set(
            pair
            for pair in target_pairs
            if pair in default_border_pairs
        )

    if not target_pairs_set:
        # Fallback: keep the two canonical border-long-range candidates if none were provided.
        target_pairs_set = set(default_border_pairs)

    preflight = source_payload.get("preflight")
    forced_lane_pairs: dict[Pair, str] = {}
    effective_anchor_pairs: set[Pair] = set()
    if isinstance(preflight, dict):
        preflight_effective = preflight.get("effective")
        if isinstance(preflight_effective, dict):
            for lane_name, raw_pairs in preflight_effective.items():
                if not isinstance(raw_pairs, list):
                    continue
                for item in raw_pairs:
                    if isinstance(item, (list, tuple)) and len(item) == 2:
                        forced_lane_pairs[_coerce_pair(item)] = str(lane_name)
                        effective_anchor_pairs.add(_coerce_pair(item))

    min_sep_after = {
        "strict": int(args.strict_min_separation),
        "balanced": int(args.balanced_min_separation),
        "balanced_rescue": int(args.rescue_min_separation),
        "monitor": int(args.monitor_min_separation),
        "unknown": int(args.long_range_min_separation),
    }

    replica_dirs = sorted(
        [
            item
            for item in source_run.glob("replica_*")
            if item.is_dir() and (item / "openmm_dca_restrained_trajectory.pdb").exists()
        ]
    )
    if not replica_dirs:
        raise SystemExit(f"no replica trajectories found in {source_run}")

    after_payloads: list[dict[Pair, dict[str, float]]] = []
    after_metadata: list[dict[str, Any]] = []
    after_reasons: list[dict[Pair, str]] = []
    after_lanes: list[dict[Pair, str]] = []
    after_reachability: list[dict[Pair, dict[str, object]]] = []

    for replica_dir in replica_dirs:
        trajectory = replica_dir / "openmm_dca_restrained_trajectory.pdb"
        if not trajectory.exists():
            continue

        after_stable, after_meta, after_reason_map, after_lane_map, after_reach = _build_stable_pairs_with_lane_min_separation(
            trajectory=trajectory,
            row=row,
            sequence=row.sequence,
            tail_fraction=tail_fraction,
            contact_cutoff_angstrom=contact_cutoff,
            min_separation_by_lane=min_sep_after,
            frequency_threshold=frequency_threshold,
            dca_threshold=dca_threshold,
            chemical_threshold=chemical_threshold,
            topology_mode=topology_mode,
            domain_boundaries_raw=domain_boundaries,
            max_degree=max_degree,
            control_count=control_count,
            max_control_overlap_fraction=max_control_overlap_fraction,
            coupling_file=readout_coupling,
            anchor_profile_file=anchor_profile,
            audit_pairs=target_pairs_set,
            forced_lane_pairs=forced_lane_pairs,
            strict_dca_threshold=strict_dca_threshold,
            balanced_strong_dca_threshold=balanced_strong_dca_threshold,
            balanced_rescue_dca_threshold=balanced_rescue_dca_threshold,
            monitor_dca_threshold=monitor_dca_threshold,
            strict_chemical_threshold=strict_chemical_threshold,
            balanced_strong_chemical_threshold=balanced_strong_chemical_threshold,
            balanced_rescue_chemical_threshold=balanced_rescue_chemical_threshold,
            monitor_chemical_threshold=monitor_chemical_threshold,
            min_separation_inclusive=not args.strict_seq_boundary,
        )

        after_payloads.append(after_stable)
        after_metadata.append(after_meta)
        after_reasons.append(after_reason_map)
        after_lanes.append(after_lane_map)
        after_reachability.append(after_reach)

    if not after_payloads:
        raise SystemExit("no successful trajectories could be postprocessed")

    (
        after_voted,
        _,
        after_lane_support,
        after_lane_selected,
        _,
        after_fail_reasons,
    ) = _vote_pairs(
        after_payloads,
        vote_threshold=vote_threshold,
        min_average_chem=vote_min_average_chem,
        strict_vote_threshold=strict_vote_threshold,
        balanced_vote_threshold=balanced_vote_threshold,
        balanced_rescue_vote_threshold=rescue_vote_threshold,
        monitor_vote_threshold=monitor_vote_threshold,
        strict_min_average_chem=strict_vote_min_average_chem,
        balanced_strong_min_average_chem=balanced_vote_min_average_chem,
        balanced_rescue_min_average_chem=rescue_vote_min_average_chem,
        monitor_min_average_chem=monitor_vote_min_average_chem,
    )

    after_selected = {(left, right) for left, right, _, _, _ in after_voted}
    after_support: Counter[Pair] = Counter()
    for payload in after_payloads:
        for pair in payload:
            after_support[pair] += 1

    dca_threshold_by_lane = {}
    dca_default = balanced_rescue_dca_threshold
    if after_metadata:
        dca_default, dca_threshold_by_lane = _resolve_dca_thresholds_after(after_metadata[0])

    from run_openmm_tmd_replicas_v0 import _load_dca_scores

    dca_scores, _dca_classes = _load_dca_scores(
        coupling_file=readout_coupling,
        row=row,
        sequence_length=row.sequence_length,
    )

    border_long_range_candidates: list[dict[str, Any]] = []
    border_rescue_selected_count = 0
    border_rescue_monitor_count = 0
    local_support_count = 0

    for pair in sorted(target_pairs_set):
        left, right = pair
        seq_sep = right - left

        after_rep_status = _collect_replica_map(after_reasons, pair)
        after_lane_rep = _collect_replica_lanes(after_lanes, pair)

        after_status = _aggregate_pair_status(list(after_rep_status.values()))
        status_values = list(after_rep_status.values()) if after_rep_status else ["missing_from_tail"]

        after_lanes_seen = sorted(set(after_lane_rep.values()))
        after_lane = after_lanes_seen[0] if after_lanes_seen else "monitor"

        role = _classify_contact_role(after_lane, seq_sep)
        if role == "local_support":
            local_support_count += 1

        after_pair_reach = _collect_replica_reachability(after_reachability, pair)
        reach_summary = _summarize_pair_reachability(
            after_pair_reach,
            selected=pair in after_selected,
        )
        tail_freq = float(reach_summary.get("tail_contact_frequency", 0.0) or 0.0)
        per_replica_presence = reach_summary.get("per_replica_tail_presence", {})
        replica_presence_count = int(reach_summary.get("tail_presence_count", 0))

        selected_under_current_rule = pair in after_selected

        topology_ok = "topology" not in status_values
        dca_score = float(dca_scores.get(pair, 0.0))
        dca_support = dca_score >= dca_threshold_by_lane.get("balanced_rescue", dca_default)

        geometry_reachable = tail_freq > 0.0 or replica_presence_count > 0

        candidate_eligible = (
            role == "border_long_range_rescue"
            and topology_ok
            and pair in effective_anchor_pairs
            and tail_freq >= args.rescue_tail_frequency_threshold
        )

        selected_under_rescue_late_rule = False
        final_role_decision = "not_eligible"

        if candidate_eligible:
            if dca_support:
                selected_under_rescue_late_rule = after_support.get(pair, 0) >= args.rescue_late_vote_threshold
                if selected_under_rescue_late_rule:
                    final_role_decision = "dca_backed_late_rescue"
                else:
                    final_role_decision = "vote_unstable_contact"
            elif geometry_reachable:
                final_role_decision = "geometry_reachable_but_dca_weak_monitor_contact"
            else:
                final_role_decision = "not_eligible"

        if final_role_decision == "dca_backed_late_rescue":
            border_rescue_selected_count += 1
        elif final_role_decision == "geometry_reachable_but_dca_weak_monitor_contact":
            border_rescue_monitor_count += 1

        chemistry = 0.0
        for payload in after_payloads:
            payload_entry = payload.get(pair)
            if payload_entry and "chemical_score" in payload_entry:
                chemistry = float(payload_entry["chemical_score"])
                break
        if chemistry == 0.0:
            try:
                chemistry = float(chemical_score(row.sequence[left - 1], row.sequence[right - 1]))
            except Exception:
                chemistry = 0.0

        vote_count = int(after_support.get(pair, 0))
        vote_fail_reason = after_fail_reasons.get(pair, "missing_from_tail")

        border_long_range_candidates.append(
            {
                "pair": [left, right],
                "sequence_separation": seq_sep,
                "role": role,
                "DCA_support": bool(dca_support),
                "DCA_score": round(dca_score, 6),
                "DCA_threshold": dca_threshold_by_lane.get("balanced_rescue", dca_default),
                "tail_frequency": tail_freq,
                "replica_presence": per_replica_presence,
                "vote_count": vote_count,
                "chemical_score": chemistry,
                "topology_ok": bool(topology_ok),
                "selected_under_current_rule": bool(selected_under_current_rule),
                "selected_under_rescue_late_rule": bool(selected_under_rescue_late_rule),
                "final_role_decision": final_role_decision,
                "min_sep_status": after_status,
                "vote_failure": vote_fail_reason,
                "lane_admission": after_lane,
            }
        )

    long_range_selected_before = 0
    long_range_selected_after = 0

    selected_by_border_pairs = [
        item["pair"]
        for item in border_long_range_candidates
        if item.get("selected_under_rescue_late_rule", False)
    ]

    summary: dict[str, Any] = {
        "run_type": "GARAGE_BORDER_RESCUE_ARBITRATION_V8",
        "run_mode": "postprocess_no_new_simulation",
        "source_run_dir": str(source_run),
        "source_certificate": str(certificate_path),
        "claim_allowed": False,
        "long_range_min_separation": args.long_range_min_separation,
        "strict_seq_boundary": args.strict_seq_boundary,
        "target_pairs": sorted([list(pair) for pair in sorted(target_pairs_set)]),
        "rescue_tail_frequency_threshold": args.rescue_tail_frequency_threshold,
        "rescue_late_vote_threshold": args.rescue_late_vote_threshold,
        "border_rescue_pairs": border_long_range_candidates,
        "border_rescue_selected_count": int(border_rescue_selected_count),
        "border_rescue_monitor_count": int(border_rescue_monitor_count),
        "local_support_count": int(local_support_count),
        "selected_by_border_pairs": selected_by_border_pairs,
        "selected_by_rule": selected_by_border_pairs,
        "long_range_selected_before": long_range_selected_before,
        "long_range_selected_after": long_range_selected_after,
        "long_range_evidence_polluted": False,
        "noise_added": int(after_lane_selected.get("monitor", 0) + after_lane_selected.get("unknown", 0)),
        "lane_vote_support_after": dict(after_lane_support),
        "lane_vote_support_current_rule": dict(after_lane_selected),
        "pair_in_blueprint": _safe_ratio(len(target_pairs_set & effective_anchor_pairs), len(target_pairs_set)),
        "preflight": {
            "source_preflight_ready": bool(source_payload.get("preflight", {}).get("v5_ready", False)),
            "v5_requirements": source_payload.get("preflight", {}).get("v5_requirements", {}),
            "failure_types": source_payload.get("preflight", {}).get("failure_type")
            or source_payload.get("preflight", {}).get("notes", []),
        },
    }

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    out_file = out_root / "border_rescue_arbitration_v8_postprocess.json"
    if args.write_json:
        out_file = Path(args.write_json)
        out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
