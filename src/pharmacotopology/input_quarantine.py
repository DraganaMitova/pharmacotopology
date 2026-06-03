from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path


@dataclass(frozen=True)
class InputAtom:
    ref: str
    h: str
    n: int
    kind: str

    def field_packet(self) -> dict[str, object]:
        return {
            "ι.ref": self.ref,
            "ι.h": self.h,
            "ι.n": self.n,
            "ι.kind": self.kind,
        }


@dataclass(frozen=True)
class InputSurfaceRecord:
    ref: str
    h: str
    n: int
    kind: str
    text: str
    created_at: float


class InputQuarantine:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def ingest(self, text: str, kind: str = "operator") -> InputAtom:
        existing = self._count_existing()
        ref = f"ι.{existing + 1:06d}"
        digest = "sha256:" + sha256(text.encode("utf-8")).hexdigest()
        atom = InputAtom(ref=ref, h=digest, n=len(text), kind=kind)
        record = InputSurfaceRecord(
            ref=atom.ref,
            h=atom.h,
            n=atom.n,
            kind=atom.kind,
            text=text,
            created_at=time.time(),
        )

        with self.path.open("a", encoding="utf-8") as file:
            file.write(
                json.dumps(asdict(record), ensure_ascii=False, sort_keys=True) + "\n"
            )

        return atom

    def _count_existing(self) -> int:
        if not self.path.exists():
            return 0
        return sum(
            1
            for line in self.path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
