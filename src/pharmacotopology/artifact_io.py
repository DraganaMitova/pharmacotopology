from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Mapping, Sequence


def stable_float_text(value: float, *, digits: int = 6) -> str:
    if not math.isfinite(value):
        raise ValueError(f"non-finite float cannot be serialized: {value!r}")
    if abs(value) < 0.5 * 10**-digits:
        value = 0.0
    text = f"{value:.{digits}f}".rstrip("0").rstrip(".")
    if "." not in text:
        text += ".0"
    return text


def _stable_csv_value(value: object) -> object:
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        return stable_float_text(value)
    if isinstance(value, Mapping):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    if isinstance(value, (list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def stable_csv_row(
    row: Mapping[str, object],
    fieldnames: Sequence[str],
) -> dict[str, object]:
    return {key: _stable_csv_value(row.get(key, "")) for key in fieldnames}


def write_csv_rows(rows: Sequence[Mapping[str, object]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(stable_csv_row(row, fieldnames))
    return path
