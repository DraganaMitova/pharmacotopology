#!/usr/bin/env python3
"""Local ESMFold v1 runner.

This requires fair-esm with ESMFold dependencies and model weights. It is bounded
by the caller's timeout (the all-sources shell runner wraps it with timeout).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> str:
    seq = "".join(line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith(">"))
    seq = "".join(ch for ch in seq.upper() if ch.isalpha())
    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        raise ValueError(f"invalid amino acid letters: {invalid}")
    return seq


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", nargs="?", default="4ake.fasta")
    parser.add_argument("output_pdb", nargs="?", default="4ake_esmfold_local.pdb")
    parser.add_argument("--cpu", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--report", default=None)
    args = parser.parse_args()
    started = time.time()
    out = Path(args.output_pdb)
    report = Path(args.report) if args.report else out.with_suffix(".local_esmfold_report.json")
    payload = {
        "kind": "local_esmfold_v1_runner",
        "status": "not_started",
        "fasta": args.fasta,
        "output_pdb": str(out),
        "device": "cpu" if args.cpu else "auto_cuda_if_available",
        "chunk_size": args.chunk_size,
        "error": "",
        "elapsed_seconds": None,
        "pdb_exists": False,
        "pdb_atom_lines": 0,
    }
    try:
        seq = read_fasta(Path(args.fasta))
        payload["sequence_length"] = len(seq)
        import torch
        import esm
        model = esm.pretrained.esmfold_v1()
        model = model.eval()
        if not args.cpu and torch.cuda.is_available():
            model = model.cuda()
            payload["device"] = "cuda"
        else:
            payload["device"] = "cpu"
        if hasattr(model, "set_chunk_size"):
            model.set_chunk_size(args.chunk_size)
        with torch.no_grad():
            pdb = model.infer_pdb(seq)
        atom_lines = sum(1 for line in pdb.splitlines() if line.startswith("ATOM"))
        payload["pdb_atom_lines"] = atom_lines
        if atom_lines <= 0:
            payload["status"] = "invalid_pdb_returned"
            return_code = 3
        else:
            out.write_text(pdb, encoding="utf-8")
            payload["status"] = "success"
            return_code = 0
    except Exception as exc:  # noqa: BLE001
        payload["status"] = "exception"
        payload["error"] = f"{type(exc).__name__}: {exc}"
        return_code = 10
    finally:
        payload["elapsed_seconds"] = round(time.time() - started, 3)
        payload["pdb_exists"] = out.exists() and out.stat().st_size > 0
        if payload["pdb_exists"]:
            payload["pdb_atom_lines"] = sum(1 for line in out.read_text(errors="replace").splitlines() if line.startswith("ATOM"))
        report.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        print(json.dumps(payload, indent=2, sort_keys=True))
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
