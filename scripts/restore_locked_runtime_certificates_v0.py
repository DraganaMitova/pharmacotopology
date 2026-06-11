#!/usr/bin/env python3
"""Restore exported locked certificates back into the runtime artifact tree.

This restores JSON certificates only. It cannot restore heavy trajectories or MD
outputs, and it deliberately refuses to invent missing science artifacts.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

ROOT = Path(os.environ.get("REPO_ROOT", Path.cwd())).resolve()
RUN_ROOT = ROOT / "first_contact_clean_pharmacotopology_layer_run"
IN = ROOT / "data" / "locked_runtime_certificates"


def main() -> None:
    index_path = IN / "locked_runtime_certificates_index.json"
    if not index_path.exists():
        raise SystemExit(f"missing export index: {index_path}")
    index = json.loads(index_path.read_text(encoding="utf-8"))
    restored = []
    missing = []
    for rel in index.get("exported", []):
        src = ROOT / rel
        if not src.exists():
            missing.append(str(src))
            continue
        # rel is data/locked_runtime_certificates/<DIR>/<certificate>
        parts = Path(rel).parts
        if len(parts) < 4 or parts[0] != "data" or parts[1] != "locked_runtime_certificates":
            continue
        runtime_dir = parts[2]
        filename = parts[3]
        dst_dir = RUN_ROOT / runtime_dir
        dst = dst_dir / filename
        json.loads(src.read_text(encoding="utf-8"))
        dst_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        restored.append(str(dst.relative_to(ROOT)))
    report = {
        "kind": "LOCKED_RUNTIME_CERTIFICATE_RESTORE_REPORT_v0",
        "claim_allowed": False,
        "restored_count": len(restored),
        "missing_count": len(missing),
        "restored": restored,
        "missing": missing,
        "policy": "certificates_only_no_trajectories_no_md_outputs",
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
