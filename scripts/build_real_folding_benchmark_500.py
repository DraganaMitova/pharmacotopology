from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_reference_loader import (  # noqa: E402
    load_folding_reference_dataset,
)
from pharmacotopology.folding_structure_features import (  # noqa: E402
    build_locked_benchmark_payload,
)
from pharmacotopology.folding_topology import FoldingReferenceExample  # noqa: E402


DEFAULT_OUTPUT_PATH = Path("data/folding_benchmarks_real_500.locked.json")


def _git_commit_hash() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    return result.stdout.strip()


def build_real_folding_benchmark_payload(
    *,
    input_path: Path | None,
    target_size: int,
    lock_requested: bool,
    recipe_commit_hash: str = "",
) -> dict[str, object]:
    references: Sequence[FoldingReferenceExample] = ()
    if input_path is not None:
        dataset = load_folding_reference_dataset(input_path, require_external=True)
        references = dataset.references

    payload = build_locked_benchmark_payload(
        references,
        target_size=target_size,
        lock_requested=lock_requested,
        recipe_commit_hash=recipe_commit_hash,
    )
    payload["input_benchmark_file"] = str(input_path) if input_path else ""
    if input_path is None:
        payload["status"] = "empty_locked_dataset_shell_no_external_rows_attached"
        payload["builder_note"] = (
            "No input file was supplied, so this output is a target shell. "
            "It is not evidence and will not pass --require-external benchmark runs."
        )
    else:
        payload["status"] = (
            "locked_external_benchmark"
            if payload["locked_after_generation"]
            else "external_benchmark_loaded_but_not_locked"
        )
    return payload


def write_real_folding_benchmark_payload(
    payload: dict[str, object],
    output_path: Path,
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a locked real folding benchmark file from externally sourced "
            "reference rows. Without --input, writes an honest empty target shell."
        )
    )
    parser.add_argument(
        "--input",
        default="",
        help="Optional externally sourced benchmark JSON to validate and lock.",
    )
    parser.add_argument(
        "--size",
        type=int,
        default=500,
        help="Target benchmark size. Default: 500.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Output benchmark JSON path.",
    )
    parser.add_argument(
        "--lock",
        action="store_true",
        help="Request locked/no-retuning certificate when blockers are clear.",
    )
    parser.add_argument(
        "--recipe-commit",
        default="",
        help="Optional commit hash to record. Defaults to current git HEAD.",
    )
    args = parser.parse_args()

    if args.size <= 0:
        parser.error("--size must be positive")

    input_path = Path(args.input) if args.input else None
    recipe_commit_hash = args.recipe_commit or _git_commit_hash()
    try:
        payload = build_real_folding_benchmark_payload(
            input_path=input_path,
            target_size=args.size,
            lock_requested=args.lock,
            recipe_commit_hash=recipe_commit_hash,
        )
    except ValueError as exc:
        parser.error(str(exc))

    output = write_real_folding_benchmark_payload(payload, Path(args.output))
    print(output)
    blockers = payload.get("lock_certificate", {}).get("lock_blockers", [])
    if blockers:
        print("lock_blockers=" + ",".join(str(item) for item in blockers))


if __name__ == "__main__":
    main()
