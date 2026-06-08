#!/usr/bin/env python3
from __future__ import annotations

"""Run the repo's default non-visual test suite without external pytest plugins."""

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    env = dict(os.environ)
    env.setdefault("PYTHONPATH", "src")
    env.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    env.setdefault("PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS", "60")
    command = [sys.executable, "-m", "pytest", "-q"]
    raise SystemExit(subprocess.run(command, cwd=REPO_ROOT, env=env).returncode)


if __name__ == "__main__":
    main()
