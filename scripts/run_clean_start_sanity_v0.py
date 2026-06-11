#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
V32_CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT" / "v32_external_constraint_source_import_preflight_certificate.json"


def _env() -> dict[str, str]:
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    prefix = f"{REPO_ROOT / 'src'}{os.pathsep}{REPO_ROOT / 'scripts'}"
    env["PYTHONPATH"] = prefix if not existing else f"{prefix}{os.pathsep}{existing}"
    env.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    env.setdefault("REPO_ROOT", str(REPO_ROOT))
    return env


def _run(label: str, cmd: list[str], timeout_s: int) -> dict[str, Any]:
    print(f"\n=== {label} ===")
    print(" ".join(cmd))
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            env=_env(),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired as exc:
        out = exc.stdout or ""
        if isinstance(out, bytes):
            out = out.decode("utf-8", errors="replace")
        print(out[-4000:])
        return {"label": label, "ok": False, "timeout": True, "returncode": None}

    print(proc.stdout[-8000:])
    return {"label": label, "ok": proc.returncode == 0, "timeout": False, "returncode": proc.returncode}


def _verify_v32_certificate() -> dict[str, Any]:
    print("\n=== verify V32 certificate ===")
    cert = json.loads(V32_CERT.read_text(encoding="utf-8"))
    summary = {
        "preflight_status": cert.get("preflight_status"),
        "selected_V33_target": cert.get("selected_V33_target"),
        "selected_V33_panel": cert.get("selected_V33_panel"),
        "claim_allowed": cert.get("claim_allowed"),
        "new_MD_allowed": cert.get("new_MD_allowed"),
        "new_MD_recommended": cert.get("new_MD_recommended"),
        "provenance_clean": cert.get("provenance_clean"),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    ok = (
        summary["preflight_status"] in {
            "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED",
            "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED",
        }
        and summary["claim_allowed"] is False
        and summary["new_MD_allowed"] is False
        and summary["new_MD_recommended"] is False
        and summary["provenance_clean"] is True
    )
    return {"label": "verify V32 certificate", "ok": ok, "timeout": False, "returncode": 0 if ok else 1, "summary": summary}


def main() -> int:
    results: list[dict[str, Any]] = []
    results.append(_run("compile source/scripts", [sys.executable, "-m", "compileall", "-q", "src", "scripts"], 60))
    results.append(_run(
        "targeted V30-V32 tests only",
        [
            sys.executable, "-m", "pytest", "-q",
            "tests/test_v30_external_constraint_and_coupling_acquisition_sprint.py",
            "tests/test_v31_constraint_backed_operator_readout_preflight.py",
            "tests/test_v32_external_constraint_source_import_preflight.py",
        ],
        60,
    ))
    results.append(_run("rerun V32 import preflight", [sys.executable, "scripts/run_v32_external_constraint_source_import_preflight_v0.py"], 60))
    results.append(_verify_v32_certificate())

    ok = all(item.get("ok") for item in results)
    print("\n=== clean-start sanity summary ===")
    print(json.dumps({"ok": ok, "results": results}, indent=2, sort_keys=True))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
