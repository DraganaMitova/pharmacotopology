#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

def source_id_for(path: Path) -> str:
    name = path.stem.lower()
    if "omega" in str(path).lower():
        return "omegafold_single_sequence_4ake"
    if "esmfold2" in str(path).lower():
        return "esmfold2_single_sequence_4ake"
    if "esmfold" in str(path).lower() or "esm" in str(path).lower():
        return "esmfold_single_sequence_4ake"
    if "spired" in str(path).lower():
        return "spired_single_sequence_4ake"
    return "custom_single_sequence_4ake"

def atom_count(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith("ATOM"))

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--candidate-root", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--min-votes", type=int, default=1)
    args = ap.parse_args()
    project = Path(args.project).resolve()
    cand_root = Path(args.candidate_root).resolve()
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    results=[]
    for pdb in sorted(cand_root.rglob("*.pdb")):
        atoms = atom_count(pdb)
        if atoms <= 0:
            continue
        sid = source_id_for(pdb)
        od = out_root / pdb.stem
        cmd = [
            sys.executable,
            str(project / "scripts" / "run_msa_free_single_sequence_structure_probe_v0.py"),
            "--source-accession", "4AKE:A",
            "--predicted-pdb", str(pdb),
            "--predicted-source-id", sid,
            "--predicted-pdb-chain", "A",
            "--out-dir", str(od),
            "--min-votes", str(args.min_votes),
            "--include-sequence-physical-priors",
        ]
        env = dict(**__import__("os").environ)
        env["PYTHONPATH"] = str(project / "src")
        completed = subprocess.run(cmd, cwd=str(project), env=env, text=True, capture_output=True, timeout=240)
        report_path = od / "msa_free_single_sequence_structure_probe_report.json"
        record = {"pdb": str(pdb), "atom_lines": atoms, "source_id": sid, "returncode": completed.returncode, "report_path": str(report_path), "stdout_tail": completed.stdout[-2000:], "stderr_tail": completed.stderr[-2000:]}
        if report_path.exists():
            try:
                rep = json.loads(report_path.read_text())
                safety = rep.get("safety", {})
                record["folding_problem_solved"] = safety.get("folding_problem_solved")
                record["folding_solution_mode"] = safety.get("folding_solution_mode")
                record["direct_structure_solved"] = safety.get("direct_structure_solved")
                record["ensemble_contact_collapse_solved"] = safety.get("ensemble_contact_collapse_solved")
                record["direct_precision"] = rep.get("direct_structure_metric", {}).get("native_contact_precision")
                record["direct_recall"] = rep.get("direct_structure_metric", {}).get("native_contact_recall")
                record["ensemble_precision"] = rep.get("report", {}).get("contact_precision")
                record["ensemble_recall"] = rep.get("report", {}).get("contact_recall")
                record["script_safety_rejection"] = rep.get("safety", {}).get("script_safety_rejection")
            except Exception as e:
                record["parse_error"] = str(e)
        results.append(record)
    payload={"kind":"probe_candidate_pdbs_v1","candidate_root":str(cand_root),"out_root":str(out_root),"tested_count":len(results),"results":results}
    print(json.dumps(payload, indent=2, sort_keys=True))
    (out_root / "candidate_probe_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True))
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
