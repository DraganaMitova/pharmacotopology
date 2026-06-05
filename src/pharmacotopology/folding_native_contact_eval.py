from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from typing import Iterable, Sequence


ContactPair = tuple[int, int]


@dataclass(frozen=True)
class ContactMetricPacket:
    native_contact_count: int
    predicted_contact_count: int
    true_positive_contacts: int
    false_positive_contacts: int
    false_negative_contacts: int
    native_contact_recall: float
    native_contact_precision: float
    contact_map_f1: float
    long_range_contact_recall: float
    short_range_contact_recall: float
    false_contact_rate: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def normalized_contact_pair(pair: Sequence[int]) -> ContactPair:
    if len(pair) != 2:
        raise ValueError("contact pair must contain exactly two indexes")
    left = int(pair[0])
    right = int(pair[1])
    if left == right:
        raise ValueError("contact pair cannot self-contact")
    return (left, right) if left < right else (right, left)


def normalized_contact_pairs(pairs: Iterable[Sequence[int]]) -> tuple[ContactPair, ...]:
    return tuple(sorted({normalized_contact_pair(pair) for pair in pairs}))


def contact_map_hash(pairs: Iterable[Sequence[int]]) -> str:
    normalized = normalized_contact_pairs(pairs)
    encoded = ";".join(f"{left}-{right}" for left, right in normalized)
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def evaluate_contact_prediction(
    *,
    native_pairs: Iterable[Sequence[int]],
    predicted_pairs: Iterable[Sequence[int]],
    long_range_threshold: int = 24,
    short_range_threshold: int = 12,
) -> ContactMetricPacket:
    native = set(normalized_contact_pairs(native_pairs))
    predicted = set(normalized_contact_pairs(predicted_pairs))
    true_positive = native & predicted
    false_positive = predicted - native
    false_negative = native - predicted
    native_long = {
        pair for pair in native if abs(pair[1] - pair[0]) >= long_range_threshold
    }
    native_short = {
        pair for pair in native if abs(pair[1] - pair[0]) <= short_range_threshold
    }
    true_positive_long = true_positive & native_long
    true_positive_short = true_positive & native_short

    recall = len(true_positive) / len(native) if native else (1.0 if not predicted else 0.0)
    precision = (
        len(true_positive) / len(predicted)
        if predicted
        else (1.0 if not native else 0.0)
    )
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    long_recall = (
        len(true_positive_long) / len(native_long)
        if native_long
        else 1.0
    )
    short_recall = (
        len(true_positive_short) / len(native_short)
        if native_short
        else 1.0
    )
    false_rate = len(false_positive) / len(predicted) if predicted else 0.0
    return ContactMetricPacket(
        native_contact_count=len(native),
        predicted_contact_count=len(predicted),
        true_positive_contacts=len(true_positive),
        false_positive_contacts=len(false_positive),
        false_negative_contacts=len(false_negative),
        native_contact_recall=_rounded(recall),
        native_contact_precision=_rounded(precision),
        contact_map_f1=_rounded(f1),
        long_range_contact_recall=_rounded(long_recall),
        short_range_contact_recall=_rounded(short_recall),
        false_contact_rate=_rounded(false_rate),
    )
