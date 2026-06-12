from __future__ import annotations

import io
import tokenize
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"


def test_engine_has_no_arbitrary_decimal_threshold_or_coefficient_literals() -> None:
    tokens = tokenize.generate_tokens(io.StringIO(ENGINE.read_text(encoding="utf-8")).readline)
    forbidden: list[tuple[int, str]] = []
    for token in tokens:
        if token.type != tokenize.NUMBER or "." not in token.string:
            continue
        if token.string in {"0.0", "1.0"}:
            continue
        forbidden.append((token.start[0], token.string))

    assert forbidden == []
