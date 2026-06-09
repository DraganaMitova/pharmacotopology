#!/usr/bin/env python3
from __future__ import annotations

"""Runtime-oriented postprocess-only V9 role-aware rescue selector.

Promotes V8 border-rescue arbitration into a role-aware final selected set.
This keeps the same V5 trajectories as input and emits category-separated
selection outputs:
- strict scaffold
- balanced core
- border long-range rescue
- monitor-only border rescue
- local support
- medium support
- diagnostic shell
"""

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from run_openmm_tmd_replicas_v0 import (  # type: ignore[import-not-found]
    DEFAULT_BENCHMARK_FILE,
    DEFAULT_EXTERNAL_COUPLING_FILE,
    _load_dca_scores,
    _load_row,
    _vote_pairs,
)
from pharmacotopology.folding_template_docking import chemical_score

from run_garage_rescue_role_selector_v7 import (  # type: ignore[import-not-found]
    _aggregate_pair_status,
    _build_stable_pairs_with_lane_min_separation,
    _collect_replica_lanes,
    _collect_replica_map,
    _collect_replica_reachability,
    _classify_contact_role,
    _coerce_pair,
    _extract_target_pairs,
    _read_json,
    _safe_ratio,
    _summarize_pair_reachability,
)

Pair = tuple[int, int]


def _resolve_dca_thresholds_after(first_metadata: dict[str, Any]) -> tuple[float, dict[str, float]]:
    lane_thresholds = first_metadata.get("resolved_dca_threshold_by_lane")
    if isinstance(lane_thresholds, dict):
        return lane_thresholds.get("balanced_rescue", 0.0), {
            str(k): float(v) if isinstance(v, (int, float)) else float(v or 0.0)
            for k, v in lane_thresholds.items()
        }
    return 0.0, {}


def _category_bucket(after_lane: str, role: str, final_decision: str) -> str:
    if role == "border_long_range_rescue" and final_decision == "selected_border_rescue":
        return "selected_border_rescue"
    if role == "border_long_range_rescue" and final_decision == "monitor_only_border_rescue":
        return "monitor_only_border_rescue"
    if role == "local_support":
        return "selected_local_support"
    if role == "medium_support":
        return "selected_medium_support"
    if role == "diagnostic_shell":
        return "diagnostic_shell"
    if role == "true_long_range_core":
        if after_lane == "strict":
            return "selected_strict_scaffold"
        if after_lane == "balanced":
            return "selected_balanced_core"
    return "unassigned"


def _long_range_readout(
    selected_strict_scaffold: set[Pair],
    selected_balanced_core: set[Pair],
    selected_border_rescue: set[Pair],
) -> list[list[int]]:
    selected = sorted(selected_strict_scaffold | selected_balanced_core | selected_border_rescue)
    return [list(pair) for pair in selected]


def _support_readout(
    selected_local_support: set[Pair],
    selected_medium_support: set[Pair],
) -> list[list[int]]:
    selected = sorted(selected_local_support | selected_medium_support)
    return [list(pair) for pair in selected]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run GARAGE_ROLE_AWARE_RESCUE_SELECTOR_V9 postprocess on existing V5 trajectories.",
    )
    parser.add_argument(
        "--source-run-dir",
        default=str(
            REPO_ROOT
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
            REPO_ROOT
            / "first_contact_clean_pharmacotopology_layer_run"
            / "GARAGE_ROLE_AWARE_RESCUE_SELECTOR_V9"
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
    parser.add_argument("--long-range-min-separation", type=int, default=24)
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
    if not target_pairs:
        raise SystemExit("target pairs are empty; provide --audit-pairs-json or source pair_audit")
    target_pairs_set = set(target_pairs)

    forced_lane_pairs: dict[Pair, str] = {}
    preflight = source_payload.get("preflight")
    if isinstance(preflight, dict):
        preflight_effective = preflight.get("effective")
        if isinstance(preflight_effective, dict):
            for lane_name, raw_pairs in preflight_effective.items():
                if not isinstance(raw_pairs, list):
                    continue
                for item in raw_pairs:
                    if isinstance(item, (list, tuple)) and len(item) == 2:
                        forced_lane_pairs[_coerce_pair(item)] = str(lane_name)

    effective_anchor_pairs: set[Pair] = set()
    if isinstance(preflight, dict):
        preflight_effective = preflight.get("effective")
        if isinstance(preflight_effective, dict):
            for raw_pairs in preflight_effective.values():
                if not isinstance(raw_pairs, list):
                    continue
                for item in raw_pairs:
                    if isinstance(item, (list, tuple)) and len(item) == 2:
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

        stable, meta, reason_map, lane_map, reach = _build_stable_pairs_with_lane_min_separation(
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

        after_payloads.append(stable)
        after_metadata.append(meta)
        after_reasons.append(reason_map)
        after_lanes.append(lane_map)
        after_reachability.append(reach)

    if not after_payloads:
        raise SystemExit("no successful trajectories could be postprocessed")

    after_voted, _, after_lane_support, _, _, after_fail_reasons = _vote_pairs(
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
    after_selected_by_lane: dict[str, set[Pair]] = {
        "strict": set(),
        "balanced": set(),
        "balanced_rescue": set(),
        "monitor": set(),
        "unknown": set(),
    }
    for left, right, _, _, lane in after_voted:
        after_selected_by_lane.setdefault(lane, set()).add((left, right))

    after_support: Counter[Pair] = Counter()
    for payload in after_payloads:
        for pair in payload:
            after_support[pair] += 1

    dca_default, dca_threshold_by_lane = _resolve_dca_thresholds_after(after_metadata[0])
    dca_scores, _ = _load_dca_scores(
        coupling_file=readout_coupling,
        row=row,
        sequence_length=row.sequence_length,
    )

    selected_strict_scaffold: set[Pair] = set()
    selected_balanced_core: set[Pair] = set()
    selected_border_rescue: set[Pair] = set()
    monitor_only_border_rescue: set[Pair] = set()
    selected_local_support: set[Pair] = set()
    selected_medium_support: set[Pair] = set()
    diagnostic_shell_selected: set[Pair] = set()

    border_rescue_selected_count = 0
    border_rescue_monitor_count = 0
    local_support_count = 0

    pair_reports: list[dict[str, Any]] = []

    for pair in sorted(target_pairs_set):
        left, right = pair
        seq_sep = right - left

        after_rep_status = _collect_replica_map(after_reasons, pair)
        after_lane_rep = _collect_replica_lanes(after_lanes, pair)
        after_status = _aggregate_pair_status(list(after_rep_status.values()))

        after_lanes_seen = sorted(set(after_lane_rep.values()))
        after_lane = after_lanes_seen[0] if after_lanes_seen else "monitor"

        role = _classify_contact_role(after_lane, seq_sep)
        selected_under_current_rule = pair in after_selected

        after_pair_reach = _collect_replica_reachability(after_reachability, pair)
        reach_summary = _summarize_pair_reachability(after_pair_reach, selected=selected_under_current_rule)
        tail_freq = float(reach_summary.get("tail_contact_frequency", 0.0) or 0.0)
        replica_presence_count = int(reach_summary.get("tail_presence_count", 0))
        replica_presence = reach_summary.get("per_replica_tail_presence", {})

        topology_ok = "topology" not in after_status
        dca_threshold_used = float(dca_threshold_by_lane.get("balanced_rescue", dca_default))
        dca_score = float(dca_scores.get(pair, 0.0))
        dca_support = dca_score >= dca_threshold_used

        geometry_reachable = tail_freq > 0.0 or replica_presence_count > 0

        final_decision = "not_eligible"
        selected_under_rescue_late_rule = False

        if role == "border_long_range_rescue":
            candidate_eligible = (
                topology_ok
                and pair in effective_anchor_pairs
                and tail_freq >= args.rescue_tail_frequency_threshold
            )

            if candidate_eligible and dca_support:
                selected_under_rescue_late_rule = after_support.get(pair, 0) >= args.rescue_late_vote_threshold
                if selected_under_rescue_late_rule:
                    final_decision = "selected_border_rescue"
                    selected_border_rescue.add(pair)
                    border_rescue_selected_count += 1
                else:
                    final_decision = "vote_unstable_border_rescue"
                    border_rescue_monitor_count += 1
            elif candidate_eligible and geometry_reachable:
                final_decision = "monitor_only_border_rescue"
                monitor_only_border_rescue.add(pair)
                border_rescue_monitor_count += 1

        if role in {"local_support", "medium_support", "diagnostic_shell", "true_long_range_core"}:
            if selected_under_current_rule:
                bucket = _category_bucket(after_lane, role, final_decision)
                if bucket == "selected_local_support":
                    selected_local_support.add(pair)
                    local_support_count += 1
                elif bucket == "selected_medium_support":
                    selected_medium_support.add(pair)
                elif bucket == "diagnostic_shell":
                    diagnostic_shell_selected.add(pair)
                elif bucket == "selected_strict_scaffold":
                    selected_strict_scaffold.add(pair)
                elif bucket == "selected_balanced_core":
                    selected_balanced_core.add(pair)

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

        pair_reports.append(
            {
                "pair": [left, right],
                "sequence_separation": seq_sep,
                "role": role,
                "lane_admission": after_lane,
                "min_sep_status": after_status,
                "DCA_support": bool(dca_support),
                "DCA_score": round(dca_score, 6),
                "DCA_threshold": dca_threshold_used,
                "tail_frequency": tail_freq,
                "replica_presence": replica_presence,
                "vote_count": vote_count,
                "chemical_score": chemistry,
                "topology_ok": bool(topology_ok),
                "selected_under_current_rule": bool(selected_under_current_rule),
                "selected_under_rescue_late_rule": bool(selected_under_rescue_late_rule),
                "final_role_decision": final_decision,
                "vote_failure": vote_fail_reason,
            }
        )

    selected_pairs: set[Pair] = set().union(
        selected_strict_scaffold,
        selected_balanced_core,
        selected_border_rescue,
        selected_local_support,
        selected_medium_support,
        monitor_only_border_rescue,
        diagnostic_shell_selected,
    )

    long_range_selected = _long_range_readout(
        selected_strict_scaffold,
        selected_balanced_core,
        selected_border_rescue,
    )
    support_selected = _support_readout(selected_local_support, selected_medium_support)

    summary: dict[str, Any] = {
        "run_type": "GARAGE_ROLE_AWARE_RESCUE_SELECTOR_V9",
        "run_mode": "postprocess_no_new_simulation",
        "claim_allowed": False,
        "source_run_dir": str(source_run),
        "source_certificate": str(certificate_path),
        "target_pairs": sorted([list(pair) for pair in sorted(target_pairs_set)]),
        "rescue_tail_frequency_threshold": args.rescue_tail_frequency_threshold,
        "rescue_late_vote_threshold": args.rescue_late_vote_threshold,
        "selected_strict_scaffold": sorted([list(pair) for pair in sorted(selected_strict_scaffold)]),
        "selected_balanced_core": sorted([list(pair) for pair in sorted(selected_balanced_core)]),
        "selected_border_rescue": sorted([list(pair) for pair in sorted(selected_border_rescue)]),
        "monitor_only_border_rescue": sorted([list(pair) for pair in sorted(monitor_only_border_rescue)]),
        "selected_local_support": sorted([list(pair) for pair in sorted(selected_local_support)]),
        "selected_medium_support": sorted([list(pair) for pair in sorted(selected_medium_support)]),
        "selected_diagnostic_shell": sorted([list(pair) for pair in sorted(diagnostic_shell_selected)]),
        "long_range_evidence": {
            "pairs": long_range_selected,
            "count": len(long_range_selected),
        },
        "support_evidence": {
            "pairs": support_selected,
            "count": len(support_selected),
        },
        "monitor_only": {
            "pairs": sorted([list(pair) for pair in sorted(monitor_only_border_rescue)]),
            "count": len(monitor_only_border_rescue),
        },
        "diagnostic_only": {
            "pairs": sorted([list(pair) for pair in sorted(diagnostic_shell_selected)]),
            "count": len(diagnostic_shell_selected),
        },
        "long_range_evidence_polluted": False,
        "noise_added": int(after_lane_support.get("monitor", 0) + after_lane_support.get("unknown", 0)),
        "pair_in_blueprint": _safe_ratio(len(selected_pairs & effective_anchor_pairs), len(selected_pairs)),
        "selected_pair_count": len(selected_pairs),
        "selected_pairs": sorted([list(pair) for pair in sorted(selected_pairs)]),
        "selected_by_lane_after": {
            lane: sorted([list(pair) for pair in sorted(pairs)])
            for lane, pairs in after_selected_by_lane.items()
        },
        "selected_by_role": {
            "strict_scaffold": sorted([list(pair) for pair in sorted(selected_strict_scaffold)]),
            "balanced_core": sorted([list(pair) for pair in sorted(selected_balanced_core)]),
            "border_rescue": sorted([list(pair) for pair in sorted(selected_border_rescue)]),
            "monitor_border_rescue": sorted([list(pair) for pair in sorted(monitor_only_border_rescue)]),
            "local_support": sorted([list(pair) for pair in sorted(selected_local_support)]),
            "medium_support": sorted([list(pair) for pair in sorted(selected_medium_support)]),
            "diagnostic_shell": sorted([list(pair) for pair in sorted(diagnostic_shell_selected)]),
        },
        "metrics": {
            "border_rescue_selected_count": int(border_rescue_selected_count),
            "border_rescue_monitor_count": int(border_rescue_monitor_count),
            "local_support_count": int(local_support_count),
            "strict_balanced_long_range_count": int(len(selected_strict_scaffold) + len(selected_balanced_core)),
        },
        "rescue_pairs": pair_reports,
        "lane_vote_support_after": dict(after_lane_support),
        "preflight": {
            "source_preflight_ready": bool(source_payload.get("preflight", {}).get("v5_ready", False)),
            "v5_requirements": source_payload.get("preflight", {}).get("v5_requirements", {}),
            "failure_types": source_payload.get("preflight", {}).get("failure_type")
            or source_payload.get("preflight", {}).get("notes", []),
        },
    }

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    out_file = out_root / "role_aware_rescue_selector_v9_postprocess.json"
    if args.write_json:
        out_file = Path(args.write_json)
        out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
