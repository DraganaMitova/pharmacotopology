"""Repository-local Python startup guards.

The CI/sandbox Python environment can have unrelated third-party pytest plugins
installed globally.  Some of those plugins keep background hooks/threads alive
after this project's tiny test suite has already finished, which makes ordinary
`pytest` look like it is hanging.  Keep this repo hermetic by disabling external
plugin autoload unless the caller explicitly overrides it before startup.
"""
from __future__ import annotations

import os

os.environ.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
