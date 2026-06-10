#!/usr/bin/env python3
from __future__ import annotations

"""Terminate runaway project-local test/script workers before a bounded run.

The script intentionally refuses to kill core sandbox/system services. It only
matches Python/pytest processes whose command line points at this repository and
excludes its own PID/ancestors.
"""

import argparse
import json
import os
import signal
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PROJECT_MARKERS = (str(REPO_ROOT), "pharmacotopology", "scripts/run_", "pytest")
PROCESS_MARKERS = ("python", "pytest")


@dataclass(frozen=True)
class ProcessCandidate:
    pid: int
    ppid: int
    pcpu: float
    command: str
    args: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _ancestor_pids() -> set[int]:
    ancestors = {os.getpid()}
    current = os.getppid()
    while current and current not in ancestors:
        ancestors.add(current)
        try:
            with open(f"/proc/{current}/stat", "r", encoding="utf-8") as handle:
                fields = handle.read().split()
            current = int(fields[3])
        except OSError:
            break
    return ancestors


def _process_rows() -> tuple[ProcessCandidate, ...]:
    output = subprocess.check_output(
        ["ps", "-eo", "pid=,ppid=,pcpu=,comm=,args="],
        text=True,
    )
    rows: list[ProcessCandidate] = []
    for raw_line in output.splitlines():
        parts = raw_line.strip().split(None, 4)
        if len(parts) < 5:
            continue
        pid_raw, ppid_raw, pcpu_raw, command, args = parts
        try:
            rows.append(
                ProcessCandidate(
                    pid=int(pid_raw),
                    ppid=int(ppid_raw),
                    pcpu=float(pcpu_raw),
                    command=command,
                    args=args,
                )
            )
        except ValueError:
            continue
    return tuple(rows)


def find_project_cpu_workers(*, min_cpu: float) -> tuple[ProcessCandidate, ...]:
    ancestors = _ancestor_pids()
    matches: list[ProcessCandidate] = []
    for row in _process_rows():
        haystack = f"{row.command} {row.args}".lower()
        if row.pid in ancestors or row.ppid in ancestors:
            continue
        if row.pcpu < min_cpu:
            continue
        if not any(marker in haystack for marker in PROCESS_MARKERS):
            continue
        if not any(marker.lower() in haystack for marker in PROJECT_MARKERS):
            continue
        if "kill_project_cpu_workers_v0.py" in haystack:
            continue
        matches.append(row)
    return tuple(sorted(matches, key=lambda item: (-item.pcpu, item.pid)))


def terminate(candidates: tuple[ProcessCandidate, ...], *, dry_run: bool) -> dict[str, object]:
    terminated: list[int] = []
    if not dry_run:
        for candidate in candidates:
            try:
                os.kill(candidate.pid, signal.SIGTERM)
                terminated.append(candidate.pid)
            except ProcessLookupError:
                continue
        time.sleep(0.7)
        for candidate in candidates:
            try:
                os.kill(candidate.pid, 0)
            except ProcessLookupError:
                continue
            try:
                os.kill(candidate.pid, signal.SIGKILL)
            except ProcessLookupError:
                continue
    return {
        "repo_root": str(REPO_ROOT),
        "dry_run": dry_run,
        "matched_count": len(candidates),
        "matched": [candidate.to_dict() for candidate in candidates],
        "terminated_pids": terminated,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-cpu", type=float, default=5.0)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    candidates = find_project_cpu_workers(min_cpu=args.min_cpu)
    print(json.dumps(terminate(candidates, dry_run=args.dry_run), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
