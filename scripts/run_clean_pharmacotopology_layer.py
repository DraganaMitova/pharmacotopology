from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.layer import (  # noqa: E402
    DEFAULT_PHARMACOTOPOLOGY_SURFACE,
    run_clean_pharmacotopology_layer,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the clean KNOT pharmacotopology layer simulation."
    )
    parser.add_argument(
        "--run-dir",
        default="first_contact_clean_pharmacotopology_layer_run",
        help="Directory for clean pharmacotopology layer artifacts.",
    )
    parser.add_argument(
        "--surface",
        default=DEFAULT_PHARMACOTOPOLOGY_SURFACE,
        help="Quarantined operator surface for the simulation request.",
    )
    args = parser.parse_args()

    report = run_clean_pharmacotopology_layer(
        Path(args.run_dir),
        surface=args.surface,
    )
    print(json.dumps(asdict(report), ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
