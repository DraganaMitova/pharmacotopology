#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else ".")
rows=[]
for rp in sorted(root.rglob("msa_free_single_sequence_structure_probe_report.json")):
    try:
        rep=json.loads(rp.read_text())
    except Exception as e:
        rows.append({"report_path":str(rp),"error":str(e)})
        continue
    d=rep.get("direct_structure_metric",{})
    r=rep.get("report",{})
    s=rep.get("safety",{})
    rows.append({
        "report_path":str(rp),
        "source_accession":rep.get("source_accession"),
        "pdb_path":rep.get("single_sequence_source_run",{}).get("predicted_pdb_path"),
        "status":rep.get("single_sequence_source_run",{}).get("status"),
        "direct_precision":d.get("native_contact_precision"),
        "direct_recall":d.get("native_contact_recall"),
        "ensemble_precision":r.get("contact_precision"),
        "ensemble_recall":r.get("contact_recall"),
        "independent_structure_pair_count":r.get("independent_structure_pair_count"),
        "folding_problem_solved":s.get("folding_problem_solved"),
        "script_safety_rejection":s.get("script_safety_rejection"),
    })
solved=[x for x in rows if x.get("folding_problem_solved")]
best=sorted(rows, key=lambda x: ((x.get("direct_recall") or 0)+(x.get("ensemble_recall") or 0), (x.get("direct_precision") or 0)+(x.get("ensemble_precision") or 0)), reverse=True)[:5]
print(json.dumps({"kind":"tryhard_final_verdict_v1","report_root":str(root),"tested_reports":len(rows),"solved_count":len(solved),"folding_problem_solved":bool(solved),"best":best,"all":rows}, indent=2, sort_keys=True))
