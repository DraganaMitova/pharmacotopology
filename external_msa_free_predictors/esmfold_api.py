#!/usr/bin/env python3
"""Bounded ESMFold/ESM Atlas API runner for 4AKE or any FASTA.

Important details:
- Uses raw text/plain POST body, matching the documented ESM Atlas API usage.
- Defaults work with no arguments: ./esmfold_api.py
- Does not persist any extra raw-sequence artifacts beyond the user-supplied FASTA.
- Writes a JSON status report and exits nonzero if no valid PDB is produced.
"""
from __future__ import annotations

import argparse
import json
import socket
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    seq = "".join(line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">"))
    seq = "".join(ch for ch in seq.upper() if ch.isalpha())
    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        raise ValueError(f"invalid amino acid letters: {invalid}")
    if len(seq) < 10:
        raise ValueError("sequence too short")
    return seq


def write_report(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run ESMFold API on a FASTA and save PDB.")
    parser.add_argument("fasta", nargs="?", default="4ake.fasta")
    parser.add_argument("output_pdb", nargs="?", default="4ake_esmfold_api.pdb")
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--api-url", default=API_URL)
    parser.add_argument("--report", default=None)
    args = parser.parse_args()

    fasta = Path(args.fasta)
    output_pdb = Path(args.output_pdb)
    report_path = Path(args.report) if args.report else output_pdb.with_suffix(".esmfold_api_report.json")
    started = time.time()
    payload = {
        "kind": "esmfold_api_runner_v1",
        "api_url": args.api_url,
        "fasta": str(fasta),
        "output_pdb": str(output_pdb),
        "timeout_seconds": args.timeout_seconds,
        "status": "not_started",
        "pdb_exists": False,
        "pdb_atom_lines": 0,
        "elapsed_seconds": None,
        "error": "",
        "http_status": None,
    }
    try:
        seq = read_fasta_sequence(fasta)
        payload["sequence_length"] = len(seq)
        req = urllib.request.Request(
            args.api_url,
            data=seq.encode("utf-8"),
            headers={"Content-Type": "text/plain"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=args.timeout_seconds) as response:
            payload["http_status"] = getattr(response, "status", None)
            pdb_text = response.read().decode("utf-8", errors="replace")
        atom_lines = sum(1 for line in pdb_text.splitlines() if line.startswith("ATOM"))
        payload["pdb_atom_lines"] = atom_lines
        if atom_lines <= 0:
            payload["status"] = "invalid_pdb_returned"
            output_pdb.with_suffix(".raw_response.txt").write_text(pdb_text[:200000], encoding="utf-8")
            return_code = 3
        else:
            output_pdb.write_text(pdb_text, encoding="utf-8")
            payload["pdb_exists"] = output_pdb.exists() and output_pdb.stat().st_size > 0
            payload["status"] = "success"
            return_code = 0
    except urllib.error.HTTPError as exc:
        body = ""
        try:
            body = exc.read().decode("utf-8", errors="replace")[:4000]
        except Exception:
            body = ""
        payload.update(status="http_error", http_status=exc.code, error=f"HTTPError: {exc.reason}; body={body}")
        return_code = 4
    except urllib.error.URLError as exc:
        payload.update(status="url_error", error=f"URLError: {exc.reason}")
        return_code = 5
    except socket.gaierror as exc:
        payload.update(status="dns_error", error=f"gaierror: {exc}")
        return_code = 6
    except Exception as exc:  # noqa: BLE001 - command-line runner should report everything
        payload.update(status="exception", error=f"{type(exc).__name__}: {exc}")
        return_code = 10
    finally:
        payload["elapsed_seconds"] = round(time.time() - started, 3)
        payload["pdb_exists"] = output_pdb.exists() and output_pdb.stat().st_size > 0
        if payload["pdb_exists"] and payload["pdb_atom_lines"] == 0:
            try:
                payload["pdb_atom_lines"] = sum(
                    1 for line in output_pdb.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith("ATOM")
                )
            except Exception:
                pass
        write_report(report_path, payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
