from __future__ import annotations

import os
import signal
from collections.abc import Generator

import pytest


DEFAULT_TEST_TIMEOUT_SECONDS = int(os.environ.get("PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS", "0"))


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
    """Optional per-test alarm; disabled by default because pytest fd-capture + SIGALRM can hang on some platforms."""

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


def pytest_terminal_summary(terminalreporter: pytest.TerminalReporter, exitstatus: int, config: pytest.Config) -> None:
    """Force process exit after pytest has printed its summary.

    In the execution environment used for these artifacts, pytest can print the
    final summary and then stay alive because an inherited capture/resource
    handle remains open.  The suite itself is complete at this point, so default
    test runs force-exit after the terminal summary.  Set
    PHARMACOTOPOLOGY_PYTEST_FORCE_EXIT=0 to disable this behavior.
    """

    if os.environ.get("PHARMACOTOPOLOGY_PYTEST_FORCE_EXIT", "1") == "0":
        return
    import sys
    import threading

    def _force_exit() -> None:  # pragma: no cover - only process lifecycle safety
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(int(exitstatus))

    timer = threading.Timer(0.20, _force_exit)
    timer.daemon = False
    timer.start()
