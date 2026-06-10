#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _read_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"read_error": str(exc)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize V13 runtime provenance for a run dir.")
    parser.add_argument("run_dir")
    args = parser.parse_args()
    run_dir = Path(args.run_dir)
    cert = _read_json(run_dir / "openmm_tmd_replicas_v0_certificate.json")
    input_preflight = _read_json(run_dir / "input_preflight.json") or _read_json(run_dir / "input_preflight_failed.json")
    trajectories = sorted(run_dir.glob("replica_*/openmm_dca_restrained_trajectory.pdb"))
    logs = sorted(run_dir.glob("replica_*/run_openmm_tmd_replica.log"))
    error_hits = []
    for log in logs:
        text = log.read_text(encoding="utf-8", errors="replace")
        for marker in ("FileNotFoundError", "target_pdb_missing", "external_coupling_file_missing", "no anchors selected"):
            if marker in text:
                error_hits.append({"log": str(log), "marker": marker})
                break
    summary = {
        "run_dir": str(run_dir),
        "certificate_present": bool(cert),
        "input_preflight_status": input_preflight.get("status"),
        "input_preflight_failure_type": input_preflight.get("failure_type"),
        "target_pdb": cert.get("target_pdb") or input_preflight.get("target_pdb"),
        "successful_replicas": cert.get("successful_replicas"),
        "replica_count_in_certificate": len(cert.get("replica_summaries", [])) if cert else 0,
        "trajectory_count": len(trajectories),
        "runtime_v10_selected_pair_count": cert.get("runtime_v10_selected_pair_count"),
        "preflight_v5_ready": (cert.get("preflight") or {}).get("v5_ready") if cert else None,
        "preflight_failed_checks": ((cert.get("preflight") or {}).get("v5_requirements") or {}).get("failed_checks") if cert else None,
        "replica_error_hits": error_hits[:20],
        "claim_allowed": cert.get("claim_allowed", False) if cert else False,
        "physics_interpretation_allowed": cert.get("physics_interpretation_allowed", False) if cert else False,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
