#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = ["*.sto", "*.stockholm", "*.a2m", "*.a3m", "*.aln", "*.afa", "*.msa", "*msa*.fasta", "*msa*.fa"]
EXCLUDE_SUBSTRINGS = ["__pycache__", ".pytest_cache"]

candidates = []
for pattern in PATTERNS:
    for path in ROOT.rglob(pattern):
        rel = str(path.relative_to(ROOT))
        if any(part in rel for part in EXCLUDE_SUBSTRINGS):
            continue
        if path.is_file():
            candidates.append({"path": rel, "size_bytes": path.stat().st_size})

likely_raw_msa = [c for c in candidates if c["size_bytes"] > 1000 and not c["path"].endswith(".md")]
payload = {
    "repo_root": str(ROOT),
    "candidate_count": len(candidates),
    "likely_raw_msa_count": len(likely_raw_msa),
    "candidates": sorted(candidates, key=lambda x: x["path"]),
    "likely_raw_msa": sorted(likely_raw_msa, key=lambda x: x["path"]),
    "raw_msa_for_4ake_found": any("4ake" in c["path"].lower() for c in likely_raw_msa),
}
out_dir = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "true_local_msa_coevolution_v0"
out_dir.mkdir(parents=True, exist_ok=True)
(out_dir / "raw_msa_presence_scan.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
print(json.dumps(payload, indent=2, sort_keys=True))
