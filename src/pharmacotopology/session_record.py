from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List

from pharmacotopology.field_alphabet import FIELD_CLOSED_CHANNELS, FieldResponse


@dataclass(frozen=True)
class BoundaryStatus:
    psi_0: bool = False
    psi_1: bool = False
    psi_2: bool = False
    psi_3: bool = False
    infinity_0: bool = False
    kappa_delta: bool = False
    phi_omega: bool = False
    omega_ext: bool = True
    closed_channels: List[str] = field(
        default_factory=lambda: list(FIELD_CLOSED_CHANNELS)
    )


@dataclass(frozen=True)
class ManualSessionRecord:
    session_id: str
    turn_index: int
    iota: Dict[str, object]
    psi: str
    mu: str
    rho_n: int
    boundary_status: BoundaryStatus
    omega: str
    created_at: float
    psi_kind: str = FieldResponse.READOUT


class JsonlSessionRecordStore:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def read_records(self) -> List[ManualSessionRecord]:
        if not self.path.exists():
            return []

        records: List[ManualSessionRecord] = []
        for line in self.path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            item = json.loads(line)
            item.setdefault("psi_kind", FieldResponse.READOUT)
            item["boundary_status"] = BoundaryStatus(**item["boundary_status"])
            records.append(ManualSessionRecord(**item))
        return records

    def write_record(self, record: ManualSessionRecord) -> None:
        with self.path.open("a", encoding="utf-8") as file:
            file.write(
                json.dumps(asdict(record), ensure_ascii=False, sort_keys=True) + "\n"
            )


def default_boundary_status() -> BoundaryStatus:
    return BoundaryStatus()


def boundary_status_dict() -> Dict[str, object]:
    return asdict(default_boundary_status())


def next_turn_index(records: List[ManualSessionRecord]) -> int:
    return len(records) + 1


def now_timestamp() -> float:
    return time.time()
