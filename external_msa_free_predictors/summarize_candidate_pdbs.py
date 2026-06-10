#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path
root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
items=[]
for p in sorted(root.rglob("*.pdb")):
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
        atoms = sum(1 for line in text.splitlines() if line.startswith("ATOM"))
        size = p.stat().st_size
    except Exception:
        atoms=0; size=0
    if size > 0:
        items.append({"path": str(p), "size_bytes": size, "atom_lines": atoms, "usable": atoms > 0})
print(json.dumps({"kind":"candidate_pdb_scan_v1","root":str(root),"candidate_count":len(items),"usable_count":sum(i['usable'] for i in items),"candidates":items}, indent=2))
