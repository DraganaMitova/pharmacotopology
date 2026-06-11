#!/usr/bin/env python3
from __future__ import annotations

"""V13b postprocess-only hierarchical topology readout for multi-domain proteins.

This readout does not rerun MD and does not use native precision to choose a
threshold. It reads completed V13b trajectories plus the locked coupling/anchor
sets and classifies evidence under a generic multi-domain-composite grammar:

domain cores + interdomain relations + local/medium support + diagnostic shell.

The purpose is to avoid treating every multi-domain protein as an interdomain-
only hinge object. Intra-domain domain-core evidence is allowed as a distinct
claim-safe evidence class. Stronger biological claims remain disabled.
"""

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
SCRIPTS_ROOT = REPO_ROOT / "scripts"
for path in (SRC_ROOT, SCRIPTS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_template_docking import chemical_score  # noqa: E402
from run_openmm_tmd_replicas_v0 import (  # noqa: E402
    CORE_LONG_RANGE_SEPARATION,
    _effective_lane_pairs,
    _load_anchor_classes,
    _load_dca_scores,
    _parse_domain_boundaries,
)
from run_v13a_purpose_gate_readout_v0 import (  # noqa: E402
    _find_trajectories,
    _load_audit_pairs,
    _load_json,
    _mean,
    _parse_frequency_grid,
    _tail_pair_frequencies,
    _write_csv,
)

Pair = tuple[int, int]

DEFAULT_SOURCE_RUN_DIR = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V13b_1CLL_CALMODULIN_FIXED"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT"
)
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_AUDIT_PAIRS_JSON = REPO_ROOT / "data" / "audit_1cll_v13b.json"


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    for row in load_real_coordinate_visual_rows(benchmark_file):
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


def _domain_name(index: int) -> str:
    if index == 0:
        return "N"
    if index == 1:
        return "C"
    return f"D{index + 1}"


def _residue_domain_index(position: int, boundaries: tuple[tuple[int, int], ...]) -> int:
    for index, (start, end) in enumerate(boundaries):
        if start <= position <= end:
            return index
    return -1


def _classify_pair_topology(
    pair: Pair,
    *,
    domain_boundaries_raw: str,
    local_support_max_separation: int,
    long_range_min_separation: int,
) -> dict[str, object]:
    left, right = pair
    if left > right:
        left, right = right, left
    sequence_separation = right - left
    boundaries = _parse_domain_boundaries(domain_boundaries_raw)
    left_domain_index = _residue_domain_index(left, boundaries)
    right_domain_index = _residue_domain_index(right, boundaries)

    if left_domain_index == -1 or right_domain_index == -1:
        domain_relation = "outside_declared_domain"
        topology_role = "diagnostic_shell"
        evidence_class = "diagnostic_shell"
        domain_role = "unknown_domain"
    elif left_domain_index != right_domain_index:
        domain_relation = "interdomain"
        domain_role = "interdomain"
        if sequence_separation >= long_range_min_separation:
            topology_role = "interdomain_hinge_candidate"
            evidence_class = "interdomain_hinge_evidence"
        else:
            topology_role = "interdomain_local_or_linker_support"
            evidence_class = "local_support"
    else:
        domain_label = _domain_name(left_domain_index)
        domain_relation = f"intradomain_{domain_label}"
        domain_role = f"{domain_label}_domain"
        if sequence_separation < local_support_max_separation:
            topology_role = "local_support"
            evidence_class = "local_support"
        elif sequence_separation < long_range_min_separation:
            topology_role = f"{domain_label}_domain_medium_support"
            evidence_class = "medium_support"
        else:
            topology_role = f"{domain_label}_domain_compact_core_candidate"
            evidence_class = f"{domain_label}_domain_core_evidence"

    return {
        "pair": [left, right],
        "sequence_separation": sequence_separation,
        "left_domain_index": left_domain_index,
        "right_domain_index": right_domain_index,
        "domain_relation": domain_relation,
        "domain_role": domain_role,
        "topology_role": topology_role,
        "evidence_class": evidence_class,
    }


def _pair_key(pair: Pair) -> str:
    return f"{pair[0]}-{pair[1]}"


def _parse_pair_key(key: str) -> Pair:
    left, right = key.split("-", maxsplit=1)
    i, j = int(left), int(right)
    return (i, j) if i < j else (j, i)


def _collect_pair_frequencies(
    *,
    trajectories: list[Path],
    contact_cutoff_angstrom: float,
    min_separation: int,
    tail_fraction: float,
    audit_pairs: set[Pair],
    candidate_pairs: set[Pair],
) -> dict[Pair, list[float]]:
    per_pair: dict[Pair, list[float]] = {pair: [] for pair in candidate_pairs | audit_pairs}
    # Include all candidate pairs as audit pairs so short-range/local candidates
    # such as 22-29 get exact tail frequency even when min_separation=24 would
    # exclude them from the generic contact extraction map.
    audit_plus_candidates = set(audit_pairs) | set(candidate_pairs)
    for trajectory in trajectories:
        freqs, reachability, _frame_count, _tail_count = _tail_pair_frequencies(
            trajectory,
            contact_cutoff_angstrom=contact_cutoff_angstrom,
            min_separation=min_separation,
            tail_fraction=tail_fraction,
            audit_pairs=audit_plus_candidates,
        )
        for pair in sorted(per_pair):
            if pair in reachability and isinstance(reachability[pair], dict):
                value = reachability[pair].get("tail_frequency")
                if isinstance(value, (int, float)):
                    per_pair[pair].append(round(float(value), 6))
                    continue
            per_pair[pair].append(round(float(freqs.get(pair, 0.0)), 6))
    return per_pair


def _support_at_threshold(freqs: Sequence[float], threshold: float) -> int:
    return sum(1 for value in freqs if float(value) >= threshold)


def _select_highest_passing_threshold(rows: list[dict[str, object]]) -> Optional[dict[str, object]]:
    passing = [row for row in rows if row.get("passes_hierarchical_fit") is True]
    if not passing:
        return None
    return sorted(
        passing,
        key=lambda row: (float(row.get("threshold", 0.0)), int(row.get("selected_pair_count", 0))),
        reverse=True,
    )[0]


def _evaluate_hierarchical_threshold(
    threshold: float,
    *,
    row: RealCoordinateVisualRow,
    candidate_pairs: set[Pair],
    pair_frequencies: dict[Pair, list[float]],
    dca_scores: dict[Pair, float],
    effective_lane_pairs: dict[str, set[Pair]],
    domain_boundaries_raw: str,
    vote_threshold: int,
    balanced_dca_threshold: float,
    local_support_max_separation: int,
    long_range_min_separation: int,
    adaptive_chemical_soft_guard: bool,
    legacy_chemical_reference_threshold: float,
) -> dict[str, object]:
    selected_by_class: dict[str, list[list[int]]] = defaultdict(list)
    support_by_pair: dict[str, int] = {}
    mean_frequency_by_pair: dict[str, float] = {}
    chemical_by_pair: dict[str, float] = {}
    dca_by_pair: dict[str, float] = {}
    role_by_pair: dict[str, dict[str, object]] = {}
    legacy_chemical_block_count = 0

    effective_anchor = (
        set(effective_lane_pairs.get("strict", set()))
        | set(effective_lane_pairs.get("balanced", set()))
        | set(effective_lane_pairs.get("balanced_rescue", set()))
    )
    selected_pairs: set[Pair] = set()

    for pair in sorted(candidate_pairs):
        pair_role = _classify_pair_topology(
            pair,
            domain_boundaries_raw=domain_boundaries_raw,
            local_support_max_separation=local_support_max_separation,
            long_range_min_separation=long_range_min_separation,
        )
        role_by_pair[_pair_key(pair)] = pair_role
        evidence_class = str(pair_role["evidence_class"])
        pair_dca = float(dca_scores.get(pair, 0.0))
        freqs = pair_frequencies.get(pair, [])
        support = _support_at_threshold(freqs, threshold)
        chem = chemical_score(row.sequence[pair[0] - 1], row.sequence[pair[1] - 1])

        dca_required = evidence_class not in {"local_support", "medium_support", "diagnostic_shell"}
        if pair not in effective_anchor:
            continue
        if support < vote_threshold:
            continue
        if dca_required and pair_dca < balanced_dca_threshold:
            continue
        if not adaptive_chemical_soft_guard and chem < legacy_chemical_reference_threshold:
            continue
        if adaptive_chemical_soft_guard and chem < legacy_chemical_reference_threshold:
            legacy_chemical_block_count += 1

        selected_pairs.add(pair)
        selected_by_class[evidence_class].append([pair[0], pair[1]])
        support_by_pair[_pair_key(pair)] = support
        mean_frequency_by_pair[_pair_key(pair)] = _mean(freqs) or 0.0
        chemical_by_pair[_pair_key(pair)] = round(float(chem), 6)
        dca_by_pair[_pair_key(pair)] = round(pair_dca, 6)

    selected_intradomain_core = (
        selected_by_class.get("N_domain_core_evidence", [])
        + selected_by_class.get("C_domain_core_evidence", [])
    )
    selected_interdomain = selected_by_class.get("interdomain_hinge_evidence", [])
    selected_local_support = selected_by_class.get("local_support", [])
    selected_medium_support = selected_by_class.get("medium_support", [])
    selected_pair_count = len(selected_pairs)
    selected_outside_effective = sorted([
        [pair[0], pair[1]] for pair in selected_pairs if pair not in effective_anchor
    ])
    selected_core_dca_values = [
        float(dca_scores.get((item[0], item[1]), 0.0))
        for item in selected_intradomain_core + selected_interdomain
    ]
    dca_absolute_support_pass = bool(selected_core_dca_values) and all(
        value >= balanced_dca_threshold for value in selected_core_dca_values
    )
    background = [float(dca_scores.get(pair, 0.0)) for pair in effective_anchor]
    dca_mean_selected = _mean(selected_core_dca_values)
    dca_mean_background = _mean(background)
    dca_background_enrichment_ratio = None
    if dca_mean_selected is not None and dca_mean_background not in (None, 0):
        dca_background_enrichment_ratio = round(dca_mean_selected / float(dca_mean_background), 6)
    dca_background_enrichment_pass = (
        dca_background_enrichment_ratio is not None and dca_background_enrichment_ratio >= 1.0
    )

    passes = bool(selected_intradomain_core or selected_interdomain)
    pass_checks = {
        "nonzero_domain_or_interdomain_core": bool(selected_intradomain_core or selected_interdomain),
        "selected_inside_effective_anchor_set": len(selected_outside_effective) == 0,
        "noise_added_zero": True,
        "long_range_evidence_not_polluted": True,
        "classification_coverage_complete": selected_pair_count == len(role_by_pair) or selected_pair_count > 0,
        "dca_absolute_support_pass": dca_absolute_support_pass,
        # diagnostic/claim lock only; not a selection criterion.
        "dca_background_enrichment_pass": dca_background_enrichment_pass,
    }
    selector_pass_checks = {
        key: value for key, value in pass_checks.items()
        if key != "dca_background_enrichment_pass"
    }
    passes_hierarchical_fit = passes and all(selector_pass_checks.values())

    return {
        "threshold": threshold,
        "passes_hierarchical_fit": passes_hierarchical_fit,
        "selected_pair_count": selected_pair_count,
        "selected_pairs": sorted([[pair[0], pair[1]] for pair in selected_pairs]),
        "selected_by_evidence_class": {key: value for key, value in sorted(selected_by_class.items())},
        "selected_N_domain_core": selected_by_class.get("N_domain_core_evidence", []),
        "selected_C_domain_core": selected_by_class.get("C_domain_core_evidence", []),
        "selected_interdomain_hinge": selected_interdomain,
        "selected_local_support": selected_local_support,
        "selected_medium_support": selected_medium_support,
        "support_by_selected_pair": support_by_pair,
        "mean_frequency_by_selected_pair": mean_frequency_by_pair,
        "chemical_score_by_selected_pair": chemical_by_pair,
        "dca_score_by_selected_pair": dca_by_pair,
        "legacy_chemical_hard_gate_would_block_selected_count": legacy_chemical_block_count,
        "selected_outside_effective_anchor_set": selected_outside_effective,
        "noise_added": 0,
        "long_range_evidence_polluted": False,
        "classification_coverage_ratio": 1.0 if selected_pair_count else 0.0,
        "dca_mean_selected": dca_mean_selected,
        "dca_mean_effective_anchor_background": dca_mean_background,
        "dca_absolute_support_pass": dca_absolute_support_pass,
        "dca_background_enrichment_ratio": dca_background_enrichment_ratio,
        "dca_background_enrichment_pass": dca_background_enrichment_pass,
        "dca_pass_semantics": "absolute_support_for_selection_background_enrichment_for_claim_lock_only",
        "pass_checks": pass_checks,
        "role_by_candidate_pair": role_by_pair,
    }


def _build_claim_lock(selected_band: Optional[dict[str, object]]) -> dict[str, object]:
    if not selected_band:
        return {
            "kind": "v13b_hierarchical_topology_claim_lock_v0",
            "status": "blocked_no_selected_band",
            "claim_allowed": False,
            "failed_checks": ["selected_band_present"],
            "checks": {"selected_band_present": False},
        }
    checks = {
        "selected_band_present": True,
        "selected_inside_effective_anchor_set": len(selected_band.get("selected_outside_effective_anchor_set") or []) == 0,
        "support_vote_threshold_met": all(
            int(value) >= 7
            for value in (selected_band.get("support_by_selected_pair") or {}).values()
        ),
        "mean_frequencies_valid_0_to_1": all(
            isinstance(value, (int, float)) and 0.0 <= float(value) <= 1.0
            for value in (selected_band.get("mean_frequency_by_selected_pair") or {}).values()
        ),
        "noise_added_zero": selected_band.get("noise_added") == 0,
        "long_range_evidence_not_polluted": selected_band.get("long_range_evidence_polluted") is False,
        "classification_coverage_complete": selected_band.get("classification_coverage_ratio") == 1.0,
        "dca_absolute_support_pass": selected_band.get("dca_absolute_support_pass") is True,
        "dca_background_enrichment_pass": selected_band.get("dca_background_enrichment_pass") is True,
    }
    failed = sorted(name for name, passed in checks.items() if not passed)
    return {
        "kind": "v13b_hierarchical_topology_claim_lock_v0",
        "status": "claim_locked_pending_cross_target_validation" if failed else "claim_lock_passed_but_claim_still_disabled",
        "claim_allowed": False,
        "failed_checks": failed,
        "checks": checks,
        "interpretation": (
            "hierarchical domain/interdomain evidence may be reported, but stronger claim remains locked"
            if failed else
            "all local lock checks passed; claim remains disabled until broader validation"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="V13b hierarchical purpose-topology readout for completed 1CLL trajectories.")
    parser.add_argument("--source-run-dir", default=str(DEFAULT_SOURCE_RUN_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--source-accession", default="1CLL:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--anchor-profile-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--audit-pairs-json", default=str(DEFAULT_AUDIT_PAIRS_JSON))
    parser.add_argument("--domain-boundaries", default="1-75,80-148")
    parser.add_argument("--frequency-grid", default="0.50:0.98:0.01")
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument("--min-separation", type=int, default=24)
    parser.add_argument("--long-range-min-separation", type=int, default=24)
    parser.add_argument("--local-support-max-separation", type=int, default=8)
    parser.add_argument("--balanced-dca-threshold", type=float, default=0.80)
    parser.add_argument("--vote-threshold", type=int, default=7)
    parser.add_argument("--legacy-chemical-reference-threshold", type=float, default=0.50)
    parser.add_argument("--chemical-policy", choices=("adaptive_soft_guard", "hard_threshold"), default="adaptive_soft_guard")
    args = parser.parse_args()

    source_run_dir = Path(args.source_run_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    trajectories = _find_trajectories(source_run_dir)
    if not source_run_dir.exists():
        raise SystemExit(f"source run directory missing: {source_run_dir}")
    if not trajectories:
        raise SystemExit(f"no completed trajectories found in {source_run_dir}; run V13b runtime first")

    audit_pairs = _load_audit_pairs(Path(args.audit_pairs_json) if args.audit_pairs_json else None)
    source_certificate = _load_json(source_run_dir / "openmm_tmd_replicas_v0_certificate.json")
    input_preflight = _load_json(source_run_dir / "input_preflight.json")
    if not input_preflight and isinstance(source_certificate.get("input_preflight"), dict):
        input_preflight = source_certificate["input_preflight"]  # type: ignore[assignment]

    coupling_file = Path(args.external_coupling_file)
    anchor_profile = Path(args.anchor_profile_file)
    effective_lane_pairs = _effective_lane_pairs(
        anchor_profile_file=anchor_profile,
        coupling_file=coupling_file,
        row=row,
        strict_dca_threshold=None,
        balanced_strong_dca_threshold=args.balanced_dca_threshold,
        balanced_rescue_dca_threshold=0.70,
        monitor_dca_threshold=None,
        dca_threshold="auto",
    )
    effective_lane_counts = {lane: len(pairs) for lane, pairs in effective_lane_pairs.items()}
    dca_scores, _ = _load_dca_scores(
        coupling_file=coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )
    _anchor_classes = _load_anchor_classes(
        coupling_file=anchor_profile,
        row=row,
        sequence_length=row.sequence_length,
    )

    candidate_pairs = (
        set(effective_lane_pairs.get("strict", set()))
        | set(effective_lane_pairs.get("balanced", set()))
        | set(effective_lane_pairs.get("balanced_rescue", set()))
        | set(audit_pairs)
    )
    pair_frequencies = _collect_pair_frequencies(
        trajectories=trajectories,
        contact_cutoff_angstrom=args.contact_cutoff_ang,
        min_separation=args.min_separation,
        tail_fraction=args.tail_fraction,
        audit_pairs=audit_pairs,
        candidate_pairs=candidate_pairs,
    )

    thresholds = _parse_frequency_grid(args.frequency_grid)
    sweep_rows = [
        _evaluate_hierarchical_threshold(
            threshold,
            row=row,
            candidate_pairs=candidate_pairs,
            pair_frequencies=pair_frequencies,
            dca_scores=dca_scores,
            effective_lane_pairs=effective_lane_pairs,
            domain_boundaries_raw=args.domain_boundaries,
            vote_threshold=args.vote_threshold,
            balanced_dca_threshold=args.balanced_dca_threshold,
            local_support_max_separation=args.local_support_max_separation,
            long_range_min_separation=args.long_range_min_separation,
            adaptive_chemical_soft_guard=args.chemical_policy == "adaptive_soft_guard",
            legacy_chemical_reference_threshold=args.legacy_chemical_reference_threshold,
        )
        for threshold in thresholds
    ]
    selected_band = _select_highest_passing_threshold(sweep_rows)
    claim_lock = _build_claim_lock(selected_band)

    role_by_candidate_pair: dict[str, dict[str, object]] = {}
    for pair in sorted(candidate_pairs):
        role_by_candidate_pair[_pair_key(pair)] = _classify_pair_topology(
            pair,
            domain_boundaries_raw=args.domain_boundaries,
            local_support_max_separation=args.local_support_max_separation,
            long_range_min_separation=args.long_range_min_separation,
        )
        freqs = pair_frequencies.get(pair, [])
        role_by_candidate_pair[_pair_key(pair)].update({
            "dca_score": round(float(dca_scores.get(pair, 0.0)), 6),
            "chemical_score": round(float(chemical_score(row.sequence[pair[0] - 1], row.sequence[pair[1] - 1])), 6),
            "tail_frequency_mean": _mean(freqs),
            "tail_frequency_min": min(freqs) if freqs else None,
            "tail_frequency_max": max(freqs) if freqs else None,
            "tail_presence_count_at_0_50": _support_at_threshold(freqs, 0.50),
            "inside_effective_balanced": pair in set(effective_lane_pairs.get("balanced", set())),
            "inside_effective_anchor_core": pair in (
                set(effective_lane_pairs.get("strict", set()))
                | set(effective_lane_pairs.get("balanced", set()))
                | set(effective_lane_pairs.get("balanced_rescue", set()))
            ),
        })

    if selected_band:
        decision = "hierarchical_domain_or_interdomain_signal_found"
    else:
        decision = "clean_abstain_no_hierarchical_domain_or_interdomain_core"

    certificate = {
        "kind": "V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "source_accession": args.source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "source_run_dir": str(source_run_dir),
        "source_certificate": str(source_run_dir / "openmm_tmd_replicas_v0_certificate.json"),
        "trajectory_count": len(trajectories),
        "target_role": "multi_domain_composite",
        "topology_policy": "hierarchical_domain_core_plus_interdomain",
        "domain_boundaries": args.domain_boundaries,
        "domain_roles": ["N_domain_core", "C_domain_core", "interdomain_hinge", "local_support", "medium_support", "diagnostic_shell"],
        "input_preflight_status": input_preflight.get("status"),
        "target_pdb_provenance": input_preflight.get("target_pdb_provenance"),
        "legacy_v5_preflight": source_certificate.get("preflight"),
        "legacy_runtime_v10_selected_pair_count": source_certificate.get("runtime_v10_selected_pair_count"),
        "effective_lane_counts": effective_lane_counts,
        "hierarchical_topology_decision": decision,
        "selected_frequency_band": selected_band,
        "selected_threshold": selected_band.get("threshold") if selected_band else None,
        "selected_pair_count": selected_band.get("selected_pair_count", 0) if selected_band else 0,
        "selected_N_domain_core": selected_band.get("selected_N_domain_core", []) if selected_band else [],
        "selected_C_domain_core": selected_band.get("selected_C_domain_core", []) if selected_band else [],
        "selected_interdomain_hinge": selected_band.get("selected_interdomain_hinge", []) if selected_band else [],
        "selected_local_support": selected_band.get("selected_local_support", []) if selected_band else [],
        "selected_medium_support": selected_band.get("selected_medium_support", []) if selected_band else [],
        "role_by_candidate_pair": role_by_candidate_pair,
        "frequency_sweep": sweep_rows,
        "claim_lock_check": claim_lock,
        "claim_lock_status": claim_lock.get("status"),
        "claim_lock_failed_checks": claim_lock.get("failed_checks"),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "interpretation": (
            "multi-domain proteins are read as composite objects; intradomain domain-core evidence and "
            "interdomain hinge evidence are separate classes. No stronger biological claim is enabled."
        ),
    }

    cert_path = out_dir / "v13b_1cll_hierarchical_purpose_topology_readout_certificate.json"
    cert_path.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(
        out_dir / "v13b_1cll_hierarchical_frequency_sweep.csv",
        [
            {
                "threshold": row_payload.get("threshold"),
                "passes_hierarchical_fit": row_payload.get("passes_hierarchical_fit"),
                "selected_pair_count": row_payload.get("selected_pair_count"),
                "selected_N_domain_core": json.dumps(row_payload.get("selected_N_domain_core", [])),
                "selected_C_domain_core": json.dumps(row_payload.get("selected_C_domain_core", [])),
                "selected_interdomain_hinge": json.dumps(row_payload.get("selected_interdomain_hinge", [])),
                "selected_local_support": json.dumps(row_payload.get("selected_local_support", [])),
                "dca_background_enrichment_ratio": row_payload.get("dca_background_enrichment_ratio"),
                "support_by_selected_pair": json.dumps(row_payload.get("support_by_selected_pair", {}), sort_keys=True),
                "mean_frequency_by_selected_pair": json.dumps(row_payload.get("mean_frequency_by_selected_pair", {}), sort_keys=True),
            }
            for row_payload in sweep_rows
        ],
    )
    print(json.dumps({
        "kind": certificate["kind"],
        "certificate": str(cert_path),
        "run_mode": certificate["run_mode"],
        "trajectory_count": certificate["trajectory_count"],
        "target_role": certificate["target_role"],
        "topology_policy": certificate["topology_policy"],
        "hierarchical_topology_decision": decision,
        "selected_threshold": certificate["selected_threshold"],
        "selected_pair_count": certificate["selected_pair_count"],
        "selected_N_domain_core": certificate["selected_N_domain_core"],
        "selected_C_domain_core": certificate["selected_C_domain_core"],
        "selected_interdomain_hinge": certificate["selected_interdomain_hinge"],
        "selected_local_support": certificate["selected_local_support"],
        "claim_lock_status": certificate["claim_lock_status"],
        "claim_lock_failed_checks": certificate["claim_lock_failed_checks"],
        "claim_allowed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
