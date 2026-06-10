#!/usr/bin/env python3
"""Cross-platform timeout wrapper for macOS/Linux.
Usage: python3 run_with_timeout.py SECONDS COMMAND [ARGS...]
"""
from __future__ import annotations
import subprocess, sys, time, json


def main() -> int:
    if len(sys.argv) < 3:
        print('usage: run_with_timeout.py SECONDS COMMAND [ARGS...]', file=sys.stderr)
        return 2
    try:
        seconds = float(sys.argv[1])
    except ValueError:
        print(f'invalid timeout seconds: {sys.argv[1]}', file=sys.stderr)
        return 2
    cmd = sys.argv[2:]
    started = time.time()
    try:
        proc = subprocess.run(cmd, timeout=seconds)
        return int(proc.returncode)
    except subprocess.TimeoutExpired:
        elapsed = round(time.time() - started, 3)
        print(json.dumps({
            'kind': 'cross_platform_timeout_v1',
            'status': 'timed_out',
            'timeout_seconds': seconds,
            'elapsed_seconds': elapsed,
            'command': cmd,
        }, indent=2), file=sys.stderr)
        return 124

if __name__ == '__main__':
    raise SystemExit(main())
