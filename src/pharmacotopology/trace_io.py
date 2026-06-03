from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pharmacotopology.field_alphabet import FieldEvent, FieldStatus


def append_jsonl(path: Path, item: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as file:
        file.write(json.dumps(item, ensure_ascii=False, sort_keys=True) + "\n")


def event_group(run_dir: Path, cycle: int) -> None:
    for event in (
        (FieldEvent.ACTIVATED, FieldStatus.CONTRACT_ACTIVE),
        (FieldEvent.CANDIDATE, FieldStatus.ACTION_CANDIDATE),
        (FieldEvent.EXECUTED, FieldStatus.ACTION_RETURNED),
        (FieldEvent.MEMORY_WRITE, FieldStatus.MEMORY_WRITTEN),
        (FieldEvent.STOPPED, FieldStatus.EXTERNAL_STOP),
    ):
        append_jsonl(
            run_dir / "audit.jsonl",
            {
                "cycle": cycle if event[0] != FieldEvent.ACTIVATED else cycle - 1,
                "kind": event[0],
                "message": event[1],
            },
        )
