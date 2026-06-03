from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional, Sequence, Tuple

from pharmacotopology.field_alphabet import FieldEvent, FieldResponse


START = "START"

ALLOWED_TRANSITIONS = {
    START: frozenset({FieldEvent.ACTIVATED, FieldEvent.DENIED}),
    FieldEvent.ACTIVATED: frozenset({FieldEvent.CANDIDATE}),
    FieldEvent.CANDIDATE: frozenset(
        {FieldEvent.EXECUTED, FieldEvent.DENIED, FieldEvent.STOPPED}
    ),
    FieldEvent.EXECUTED: frozenset({FieldEvent.MEMORY_WRITE, FieldEvent.STOPPED}),
    FieldEvent.MEMORY_WRITE: frozenset({FieldEvent.STOPPED}),
    FieldEvent.DENIED: frozenset({FieldEvent.STOPPED}),
    FieldEvent.STOPPED: frozenset(),
}

KNOWN_EVENT_SYMBOLS = frozenset(
    {
        FieldEvent.ACTIVATED,
        FieldEvent.CANDIDATE,
        FieldEvent.EXECUTED,
        FieldEvent.MEMORY_WRITE,
        FieldEvent.STOPPED,
        FieldEvent.DENIED,
    }
)


@dataclass(frozen=True)
class FieldRuleViolation:
    index: int
    previous: str
    current: str
    message: str


@dataclass(frozen=True)
class FieldRuleResult:
    valid: bool
    terminal: Optional[str]
    events: Tuple[str, ...]
    violations: Tuple[FieldRuleViolation, ...]

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class FieldTraceRuleResult:
    valid: bool
    groups: Tuple[FieldRuleResult, ...]
    violations: Tuple[FieldRuleViolation, ...]

    def to_dict(self) -> dict:
        return asdict(self)


def validate_event_sequence(
    events: Sequence[str],
    *,
    speech_gate_open: bool = False,
) -> FieldRuleResult:
    previous = START
    violations: list[FieldRuleViolation] = []

    for index, current in enumerate(events):
        if current == FieldResponse.VOICE and not speech_gate_open:
            violations.append(
                FieldRuleViolation(
                    index=index,
                    previous=previous,
                    current=current,
                    message="ψ.voice is closed",
                )
            )
            previous = current
            continue

        if current not in KNOWN_EVENT_SYMBOLS:
            violations.append(
                FieldRuleViolation(
                    index=index,
                    previous=previous,
                    current=current,
                    message="unknown event symbol",
                )
            )
            previous = current
            continue

        allowed = ALLOWED_TRANSITIONS.get(previous, frozenset())
        if current not in allowed:
            violations.append(
                FieldRuleViolation(
                    index=index,
                    previous=previous,
                    current=current,
                    message=f"{current} cannot follow {previous}",
                )
            )

        previous = current

    terminal = events[-1] if events else None
    if events and terminal != FieldEvent.STOPPED:
        violations.append(
            FieldRuleViolation(
                index=len(events) - 1,
                previous=events[-2] if len(events) > 1 else START,
                current=str(terminal),
                message="active path did not terminate at ev.ω",
            )
        )

    return FieldRuleResult(
        valid=not violations,
        terminal=terminal,
        events=tuple(events),
        violations=tuple(violations),
    )


def split_event_groups(events: Sequence[str]) -> Tuple[Tuple[str, ...], ...]:
    groups: list[tuple[str, ...]] = []
    current: list[str] = []

    for event in events:
        if (
            event == FieldEvent.ACTIVATED
            and current
            and current[-1] == FieldEvent.STOPPED
        ):
            groups.append(tuple(current))
            current = []
        current.append(event)

    if current:
        groups.append(tuple(current))

    return tuple(groups)


def validate_event_groups(
    events: Sequence[str],
    *,
    speech_gate_open: bool = False,
) -> FieldTraceRuleResult:
    groups = split_event_groups(events)
    results = tuple(
        validate_event_sequence(group, speech_gate_open=speech_gate_open)
        for group in groups
    )
    violations = tuple(
        violation
        for result in results
        for violation in result.violations
    )
    return FieldTraceRuleResult(
        valid=not violations,
        groups=results,
        violations=violations,
    )
