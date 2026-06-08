#!/usr/bin/env python3
from __future__ import annotations

"""Run the repo's default non-visual test suite without external pytest plugins."""

import os
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _terminate_process_group(process: subprocess.Popen[object]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    except PermissionError:
        process.terminate()
    deadline = time.monotonic() + 2.0
    while process.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            return
        except PermissionError:
            process.kill()


def main() -> None:
    env = dict(os.environ)
    env.setdefault("PYTHONPATH", "src")
    env.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    env.setdefault("PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS", "60")
    timeout_seconds = int(env["PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS"])
    command = [sys.executable, "-m", "pytest", "-q"]
    process: subprocess.Popen[object] = subprocess.Popen(
        command,
        cwd=REPO_ROOT,
        env=env,
        start_new_session=True,
    )
    try:
        returncode = process.wait(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        _terminate_process_group(process)
        print(
            f"fast suite timed out after {timeout_seconds}s; "
            "pytest process group was killed",
            file=sys.stderr,
        )
        raise SystemExit(124)
    raise SystemExit(returncode)


if __name__ == "__main__":
    main()
