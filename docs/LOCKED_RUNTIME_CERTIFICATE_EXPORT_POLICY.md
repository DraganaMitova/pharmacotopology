# Locked Runtime Certificate Export Policy

`first_contact_clean_pharmacotopology_layer_run/` is a generated runtime output directory and may be ignored by Git because it can contain heavy trajectories and other regenerated artifacts.

However, small locked milestone certificates are reproducibility-critical. After a successful milestone run, export them into:

```text
data/locked_runtime_certificates/
```

using:

```bash
python3 scripts/export_locked_runtime_certificates_v0.py
```

The export copies JSON certificates only. It does not copy trajectories, PDB trajectories, MD outputs, or heavy data.

If runtime artifacts are accidentally deleted, restore exported certificates with:

```bash
python3 scripts/restore_locked_runtime_certificates_v0.py
```

This cannot restore missing MD trajectories or runtime-generated evidence that was never exported; it only restores previously exported locked certificates.
