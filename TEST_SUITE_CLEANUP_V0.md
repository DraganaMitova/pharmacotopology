# Test Suite Cleanup V0

The repository has been reduced to one active regression test:

```text
tests/test_4ake_alphafold_real_source_end_to_end_v0.py
```

Removed from the active test suite:

- old 1CLL recovery/forecast/perfect tests
- old adaptive-floor and threshold-probe tests
- old blind/holdout tests
- old root-level `test_*.py` smoke files
- earlier partial independent-ensemble unit tests that were superseded by the real-source end-to-end probe

The remaining test is the one that currently matters: it runs the real 4AKE AlphaFold independent-source ensemble probe from the actual script into a temporary output directory, checks the PDB provenance/hash lock, verifies the leakage gates, and asserts the locked 4AKE rescue metrics.

It includes a 30-second subprocess timeout so the suite fails clearly instead of hanging.
