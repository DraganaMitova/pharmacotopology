from __future__ import annotations

import os
import signal
from collections.abc import Generator

import pytest


DEFAULT_TEST_TIMEOUT_SECONDS = int(os.environ.get("PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS", "60"))


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-visual",
        action="store_true",
        default=False,
        help="run opt-in GIF/visual-proof tests",
    )
    parser.addoption(
        "--run-full-suite",
        action="store_true",
        default=False,
        help="run opt-in slow/full-suite tests",
    )


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "visual: opt-in GIF/visual proof test")
    config.addinivalue_line("markers", "slow: opt-in slow test")
    config.addinivalue_line("markers", "full_suite: opt-in full-suite regression")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    skip_visual = pytest.mark.skip(reason="GIF/visual proof generation is opt-in; pass --run-visual to run")
    skip_full = pytest.mark.skip(reason="slow/full-suite tests are opt-in; pass --run-full-suite to run")
    run_visual = bool(config.getoption("--run-visual"))
    run_full = bool(config.getoption("--run-full-suite"))
    for item in items:
        keywords = set(item.keywords)
        if "visual" in keywords and not run_visual:
            item.add_marker(skip_visual)
        if ({"slow", "full_suite"} & keywords) and not run_full:
            item.add_marker(skip_full)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item: pytest.Item) -> Generator[None, None, None]:
    """Hard-stop tests that accidentally return to hanging full-suite behavior."""

    if DEFAULT_TEST_TIMEOUT_SECONDS <= 0 or not hasattr(signal, "SIGALRM"):
        yield
        return

    def _timeout_handler(signum: int, frame: object) -> None:  # pragma: no cover - only fires on hang
        raise TimeoutError(
            f"pytest item exceeded {DEFAULT_TEST_TIMEOUT_SECONDS}s hard timeout: {item.nodeid}"
        )

    previous_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, DEFAULT_TEST_TIMEOUT_SECONDS)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
