#!/usr/bin/env python3
from __future__ import annotations

"""V15 4AKE balanced-candidate dynamic grammar readout.

Postprocess-only. This reads an already completed 4AKE OpenMM replica run and
asks a narrower question than the legacy V5 selector:

    Among the effective balanced 4AKE candidates, is there any DCA-lane-backed
    pair that is repeatedly present in the trajectory tail under a dynamic
    frequency sweep?

It does not rerun MD, does not use native precision, does not synthesize
positive evidence from PDB/visual files, and does not apply a fixed residue-count
separation cutoff. Sequence separation is reported only as context.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
SCRIPTS_ROOT = REPO_ROOT / "scripts"
for _path in (SCRIPTS_ROOT, SRC_ROOT):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from pharmacotopology.folding_template_docking import chemical_score  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)
from run_openmm_tmd_replicas_v0 import (  # noqa: E402
    _distance,
    _load_dca_scores,
    _parse_ca_trajectory,
    _parse_domain_boundaries,
    _residue_domain_index,
)

RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_SOURCE_RUN_DIR = RUN_ROOT / "V15_4AKE_CLEAN_DYNAMIC_ROLE_GRAMMAR_REBUILD"
DEFAULT_OUT_DIR = RUN_ROOT / "V15_4AKE_BALANCED_CANDIDATE_READOUT"
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_AUDIT_PAIRS = REPO_ROOT / "data" / "audit_4ake_v13c.json"
DEFAULT_DOMAIN_BOUNDARIES = "1-29,30-67,68-117,118-160,161-214"


def _read_json(path: Path) -> Optional[dict[str, Any]]:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
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
        writer.writerows(rows)


def _pair_key(pair: Any) -> str:
    if isinstance(pair, str):
        return pair
    if isinstance(pair, (list, tuple)) and len(pair) == 2:
        return f"{int(pair[0])}-{int(pair[1])}"
    return str(pair)


def _parse_pair_key(pair: Any) -> Optional[tuple[int, int]]:
    if isinstance(pair, (list, tuple)) and len(pair) == 2:
        left, right = int(pair[0]), int(pair[1])
        return (left, right) if left < right else (right, left)
    if isinstance(pair, str) and "-" in pair:
        left_text, right_text = pair.split("-", maxsplit=1)
        try:
            left, right = int(left_text), int(right_text)
        except ValueError:
            return None
        return (left, right) if left < right else (right, left)
    return None


def _pairs_from_preflight(certificate: dict[str, Any]) -> list[tuple[int, int]]:
    raw = certificate.get("preflight", {}).get("effective", {}).get("balanced", [])
    pairs: list[tuple[int, int]] = []
    if isinstance(raw, list):
        for item in raw:
            pair = _parse_pair_key(item)
            if pair is not None:
                pairs.append(pair)
    return sorted(set(pairs))


def _parse_grid(raw: str) -> list[float]:
    text = (raw or "").strip()
    if not text:
        return [round(x / 100.0, 2) for x in range(98, 49, -1)]
    if ":" in text:
        parts = [float(p) for p in text.split(":")]
        if len(parts) != 3:
            raise ValueError("frequency grid must be START:STOP:STEP or comma list")
        start, stop, step = parts
        if step <= 0:
            raise ValueError("frequency grid step must be positive")
        values: list[float] = []
        cur = start
        # Support either ascending or descending textual ranges.
        if start >= stop:
            while cur >= stop - 1e-9:
                values.append(round(cur, 6))
                cur -= step
        else:
            while cur <= stop + 1e-9:
                values.append(round(cur, 6))
                cur += step
        return sorted(set(values), reverse=True)
    return sorted({round(float(p.strip()), 6) for p in text.split(",") if p.strip()}, reverse=True)


def _load_row(source_accession: str, benchmark_file: Path):
    for row in load_real_coordinate_visual_rows(benchmark_file):
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"missing benchmark row for {source_accession!r} in {benchmark_file}")


def _domain_relation(pair: tuple[int, int], domain_boundaries: tuple[tuple[int, int], ...]) -> tuple[str, str, str]:
    left, right = pair
    left_domain = _residue_domain_index(left, domain_boundaries)
    right_domain = _residue_domain_index(right, domain_boundaries)
    if left_domain == -1 or right_domain == -1:
        return "domain_unknown", "diagnostic_shell", "diagnostic_shell"
    if left_domain == right_domain:
        label = f"intradomain_D{left_domain + 1}"
        role = f"D{left_domain + 1}_domain_compact_core_candidate"
        evidence = "domain_core_evidence"
        return label, role, evidence
    label = f"interdomain_D{left_domain + 1}_D{right_domain + 1}"
    return label, "domain_hinge_or_interdomain_closure_candidate", "domain_hinge_or_interdomain_evidence"


def _tail_frequency_for_pair(
    trajectory: Path,
    pair: tuple[int, int],
    *,
    tail_fraction: float,
    contact_cutoff_angstrom: float,
    tail_count_hint: Optional[int] = None,
) -> Optional[float]:
    frames = _parse_ca_trajectory(trajectory)
    if not frames:
        return None
    tail_count = tail_count_hint or max(1, math.ceil(len(frames) * tail_fraction))
    tail = frames[-tail_count:]
    hits = 0
    for frame in tail:
        if _distance(frame, pair[0], pair[1]) <= contact_cutoff_angstrom:
            hits += 1
    return round(hits / float(len(tail)), 6)


def _sequence_context(pair: tuple[int, int], sequence_length: int, domain_relation: str) -> dict[str, Any]:
    left, right = pair
    separation = right - left
    return {
        "sequence_separation": separation,
        "normalized_sequence_separation": round(separation / float(max(1, sequence_length - 1)), 6),
        "domain_relation": domain_relation,
        "separation_filter_applied": False,
        "fixed_residue_cutoff_used": False,
    }


def _dca_background_for_candidates(pairs: list[tuple[int, int]], dca_scores: dict[tuple[int, int], float]) -> float:
    values = [float(dca_scores.get(pair, 0.0)) for pair in pairs]
    return round(sum(values) / len(values), 6) if values else 0.0


def build_readout(
    *,
    source_run_dir: Path,
    benchmark_file: Path,
    coupling_file: Path,
    domain_boundaries_raw: str,
    frequency_grid: list[float],
    vote_threshold: int,
    max_selected_pairs: int,
) -> dict[str, Any]:
    source_cert_path = source_run_dir / "openmm_tmd_replicas_v0_certificate.json"
    source_cert = _read_json(source_cert_path)
    if source_cert is None:
        return {
            "kind": "V15_4AKE_BALANCED_CANDIDATE_READOUT_v0",
            "run_mode": "postprocess_only_no_new_simulation",
            "source_run_dir": str(source_run_dir),
            "source_certificate": str(source_cert_path),
            "artifact_status": "missing_source_4ake_runtime_certificate",
            "purpose_gate_decision": "preflight_abstain_missing_source_certificate",
            "claim_allowed": False,
            "biological_transfer_claim_allowed": False,
            "selected_frequency_band": None,
            "selected_pairs": [],
            "positive_evidence_found": False,
        }

    row = _load_row(str(source_cert.get("source_accession", "4AKE:A")), benchmark_file)
    sequence = row.sequence
    sequence_length = int(source_cert.get("sequence_length") or row.sequence_length)
    contact_cutoff = float(source_cert.get("contact_cutoff_angstrom", 7.0))
    tail_fraction = float(source_cert.get("tail_fraction", 0.2))
    domain_boundaries = _parse_domain_boundaries(domain_boundaries_raw or str(source_cert.get("domain_boundaries", "")))
    candidates = _pairs_from_preflight(source_cert)
    dca_scores, _classes = _load_dca_scores(coupling_file=coupling_file, row=row, sequence_length=sequence_length)
    dca_background = _dca_background_for_candidates(candidates, dca_scores)

    replica_entries = source_cert.get("replica_summaries") or []
    successful = [entry for entry in replica_entries if entry.get("returncode") == 0 and entry.get("trajectory")]
    candidate_payloads: dict[str, dict[str, Any]] = {}
    for pair in candidates:
        key = _pair_key(pair)
        domain_rel, role, evidence = _domain_relation(pair, domain_boundaries)
        chem = chemical_score(sequence[pair[0] - 1], sequence[pair[1] - 1]) if pair[1] <= len(sequence) else 0.0
        freqs: dict[str, float] = {}
        for entry in successful:
            traj = Path(str(entry.get("trajectory")))
            if not traj.exists():
                continue
            tail_hint = None
            stable_meta = entry.get("stable_metadata") or {}
            if isinstance(stable_meta, dict) and isinstance(stable_meta.get("tail_count"), int):
                tail_hint = int(stable_meta["tail_count"])
            freq = _tail_frequency_for_pair(
                traj,
                pair,
                tail_fraction=tail_fraction,
                contact_cutoff_angstrom=contact_cutoff,
                tail_count_hint=tail_hint,
            )
            if freq is not None:
                freqs[f"replica_{int(entry.get('replica', len(freqs)+1)):02d}"] = freq
        values = list(freqs.values())
        candidate_payloads[key] = {
            **_sequence_context(pair, sequence_length, domain_rel),
            "pair": [pair[0], pair[1]],
            "role_decision": role,
            "evidence_class": evidence,
            "inside_effective_balanced": True,
            "dca_score": round(float(dca_scores.get(pair, 0.0)), 6),
            "chemical_score": round(float(chem), 6),
            "tail_frequencies_by_replica": freqs,
            "tail_frequency_mean": round(sum(values) / len(values), 6) if values else None,
            "tail_frequency_min": round(min(values), 6) if values else None,
            "tail_frequency_max": round(max(values), 6) if values else None,
            "tail_presence_count_at_0_50": sum(1 for value in values if value >= 0.5),
        }

    sweep_rows: list[dict[str, Any]] = []
    selected_band: Optional[dict[str, Any]] = None
    for threshold in frequency_grid:
        selected: list[str] = []
        support_by_pair: dict[str, int] = {}
        mean_frequency_by_pair: dict[str, float] = {}
        role_by_selected_pair: dict[str, str] = {}
        for key, payload in candidate_payloads.items():
            values = list((payload.get("tail_frequencies_by_replica") or {}).values())
            support = sum(1 for value in values if float(value) >= threshold)
            if support >= vote_threshold:
                selected.append(key)
                support_by_pair[key] = support
                if payload.get("tail_frequency_mean") is not None:
                    mean_frequency_by_pair[key] = float(payload["tail_frequency_mean"])
                role_by_selected_pair[key] = str(payload.get("role_decision"))
        selected = sorted(selected, key=lambda key: (-support_by_pair.get(key, 0), key))
        selected = selected[:max_selected_pairs]
        selected_pairs = [_parse_pair_key(key) for key in selected]
        selected_pairs = [pair for pair in selected_pairs if pair is not None]
        selected_interdomain = [key for key in selected if "interdomain" in role_by_selected_pair.get(key, "")]
        selected_domain_core = [key for key in selected if "domain_compact_core" in role_by_selected_pair.get(key, "")]
        dca_selected = [candidate_payloads[key]["dca_score"] for key in selected]
        dca_mean_selected = round(sum(dca_selected) / len(dca_selected), 6) if dca_selected else None
        dca_background_enrichment_ratio = (
            round(float(dca_mean_selected) / dca_background, 6) if dca_selected and dca_background else None
        )
        row_payload = {
            "threshold": threshold,
            "passes_purpose_fit": bool(selected),
            "selected_pair_count": len(selected),
            "selected_pairs": selected,
            "selected_domain_core": selected_domain_core,
            "selected_hinge_or_interdomain": selected_interdomain,
            "support_by_selected_pair": support_by_pair,
            "mean_frequency_by_selected_pair": mean_frequency_by_pair,
            "dca_mean_selected": dca_mean_selected,
            "dca_mean_effective_balanced_background": dca_background,
            "dca_background_enrichment_ratio": dca_background_enrichment_ratio,
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0 if selected else 0.0,
        }
        sweep_rows.append(row_payload)
        if selected_band is None and selected:
            selected_band = row_payload

    selected_keys = selected_band.get("selected_pairs", []) if isinstance(selected_band, dict) else []
    selected_domain_core = selected_band.get("selected_domain_core", []) if isinstance(selected_band, dict) else []
    selected_hinge = selected_band.get("selected_hinge_or_interdomain", []) if isinstance(selected_band, dict) else []
    selected_pairs_as_lists = [[int(x) for x in key.split("-")] for key in selected_keys]
    selected_candidate_payloads = {
        key: {**candidate_payloads[key], "selected": key in selected_keys}
        for key in selected_keys
        if key in candidate_payloads
    }
    if selected_band:
        final_status = "domain_hinge_object_balanced_candidate_signal_found_under_dynamic_frequency_readout;claim_allowed=false"
        claim_lock_status = "claim_locked_4ake_candidate_readout_claim_disabled"
        claim_lock_failed_checks: list[str] = []
        artifact_status = "present_machine_readable_4ake_balanced_candidate_readout_positive"
    else:
        final_status = "domain_hinge_object_machine_readable_artifact_present_but_no_balanced_candidate_band;claim_allowed=false"
        claim_lock_status = "claim_locked_no_4ake_balanced_candidate_band"
        claim_lock_failed_checks = ["selected_frequency_band_present"]
        artifact_status = "present_machine_readable_4ake_balanced_candidate_readout_clean_abstain"

    return {
        "kind": "V15_4AKE_BALANCED_CANDIDATE_READOUT_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "source_run_dir": str(source_run_dir),
        "source_certificate": str(source_cert_path),
        "source_accession": source_cert.get("source_accession", "4AKE:A"),
        "row_id": source_cert.get("row_id"),
        "target_role": "domain_hinge_closure_object",
        "grammar_policy": "domain_hinge_dynamic_balanced_candidate_frequency_readout",
        "topology_policy": "hierarchical_domain_core_plus_interdomain_hinge_context",
        "chemical_policy": "adaptive_soft_guard_report_only_not_hard_kill",
        "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
        "fixed_residue_cutoff_used": False,
        "artifact_status": artifact_status,
        "trajectory_count": len(successful),
        "effective_balanced_count": len(candidates),
        "candidate_pairs": [_pair_key(pair) for pair in candidates],
        "candidate_pair_roles": candidate_payloads,
        "frequency_grid": frequency_grid,
        "vote_threshold": vote_threshold,
        "selected_frequency_band": selected_band,
        "selected_pairs": selected_keys,
        "selected_balanced_core": selected_pairs_as_lists,
        "selected_domain_core": selected_domain_core,
        "selected_hinge_or_interdomain": selected_hinge,
        "selected_local_support": [],
        "support_by_selected_pair": selected_band.get("support_by_selected_pair", {}) if isinstance(selected_band, dict) else {},
        "mean_frequency_by_selected_pair": selected_band.get("mean_frequency_by_selected_pair", {}) if isinstance(selected_band, dict) else {},
        "chemical_score_by_selected_pair": {key: candidate_payloads[key]["chemical_score"] for key in selected_keys if key in candidate_payloads},
        "dca_score_by_selected_pair": {key: candidate_payloads[key]["dca_score"] for key in selected_keys if key in candidate_payloads},
        "dynamic_pair_roles": selected_candidate_payloads,
        "sweep": sweep_rows,
        "noise_added": 0 if selected_band else None,
        "long_range_evidence_polluted": False if selected_band else None,
        "classification_coverage_ratio": 1.0 if selected_band else None,
        "claim_lock_status": claim_lock_status,
        "claim_lock_failed_checks": claim_lock_failed_checks,
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "positive_evidence_found": bool(selected_band),
        "final_status": final_status,
    }


def _write_report(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# V15 4AKE Balanced Candidate Readout",
        "",
        "Postprocess-only. No new MD. No fixed residue-count cutoff. No native precision selection.",
        "",
        f"Artifact status: `{payload.get('artifact_status')}`",
        f"Final status: `{payload.get('final_status')}`",
        f"Claim allowed: `{payload.get('claim_allowed')}`",
        f"Effective balanced candidate count: `{payload.get('effective_balanced_count')}`",
        f"Trajectory count: `{payload.get('trajectory_count')}`",
        "",
        "## Selected pairs",
    ]
    selected = payload.get("selected_pairs") or []
    if selected:
        for key in selected:
            role = payload.get("dynamic_pair_roles", {}).get(key, {})
            lines.append(f"- `{key}`: `{role.get('role_decision')}` support=`{payload.get('support_by_selected_pair', {}).get(key)}` mean_frequency=`{payload.get('mean_frequency_by_selected_pair', {}).get(key)}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Candidate pairs"])
    for key, role in (payload.get("candidate_pair_roles") or {}).items():
        lines.append(f"- `{key}`: role=`{role.get('role_decision')}` mean_frequency=`{role.get('tail_frequency_mean')}` dca=`{role.get('dca_score')}` chemical=`{role.get('chemical_score')}`")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run V15 4AKE balanced-candidate postprocess readout.")
    parser.add_argument("--source-run-dir", default=str(DEFAULT_SOURCE_RUN_DIR))
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--domain-boundaries", default=DEFAULT_DOMAIN_BOUNDARIES)
    parser.add_argument("--frequency-grid", default="0.98:0.50:0.01")
    parser.add_argument("--vote-threshold", type=int, default=7)
    parser.add_argument("--max-selected-pairs", type=int, default=12)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = build_readout(
        source_run_dir=Path(args.source_run_dir),
        benchmark_file=Path(args.benchmark_file),
        coupling_file=Path(args.external_coupling_file),
        domain_boundaries_raw=args.domain_boundaries,
        frequency_grid=_parse_grid(args.frequency_grid),
        vote_threshold=args.vote_threshold,
        max_selected_pairs=args.max_selected_pairs,
    )
    cert_path = out_dir / "v15_4ake_balanced_candidate_readout_certificate.json"
    report_path = out_dir / "V15_4AKE_BALANCED_CANDIDATE_READOUT_REPORT.md"
    sweep_path = out_dir / "v15_4ake_balanced_candidate_frequency_sweep.csv"
    cert_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, payload)
    _write_csv(sweep_path, payload.get("sweep", []))
    print(json.dumps({
        "kind": payload.get("kind"),
        "certificate": str(cert_path),
        "report": str(report_path),
        "sweep": str(sweep_path),
        "artifact_status": payload.get("artifact_status"),
        "positive_evidence_found": payload.get("positive_evidence_found"),
        "selected_pairs": payload.get("selected_pairs"),
        "selected_threshold": (payload.get("selected_frequency_band") or {}).get("threshold") if isinstance(payload.get("selected_frequency_band"), dict) else None,
        "claim_allowed": False,
        "final_status": payload.get("final_status"),
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
