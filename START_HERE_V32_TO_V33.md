# START HERE — V32 clean start toward V33

This package is the clean continuation point from the uploaded project.

## What is clean in this zip

Kept:

- `src/`, `scripts/`, `tests/`, `docs/`
- `data/` source/evidence material and locked certificates
- `external_msa/` because it is real source data, not cache
- small V30/V31/V32 reports/certificates under `first_contact_clean_pharmacotopology_layer_run/`

Removed:

- `.git/`
- `.venv/`
- `__MACOSX/`, `.DS_Store`, AppleDouble files
- Python caches: `__pycache__/`, `.pytest_cache/`, `*.pyc`

Repaired:

- `.gitignore` is now a real ignore file, not a pasted shell block.
- Added `scripts/run_clean_start_sanity_v0.py` so you do not need to run a hanging full suite.

## Current scientific state

The correct state is still:

```text
V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED
selected_V33_target = None
selected_V33_panel = None
claim_allowed = False
new_MD_allowed = False
new_MD_recommended = False
```

This is not a failure. It means the anti-leakage gate is working. V32 is refusing to use generated internal reports or placeholder rows as real external constraints.

## First test after unzipping

Run only the clean-start sanity runner:

```bash
python3 scripts/run_clean_start_sanity_v0.py
```

Expected result:

```text
compile source/scripts: OK
targeted V30-V32 tests only: 13 passed
rerun V32 import preflight: OK
verify V32 certificate: OK
```

Do not use the full test suite as your first check. Some older tests are not the right current gate and can waste time or appear to hang. The clean-start runner is the safe gate for the current V32→V33 transition.

## What to do next scientifically

Do not run V33 yet.

The next real move is KcsA external-source import, not MD and not operator readout.

Fill this active manifest only after you have real external files:

```text
data/external_constraints/v32_external_constraint_source_import_manifest.json
```

For KcsA, V32 requires two real external constraint buckets:

```text
data/external_constraints/KcsA/pore_filter/
data/external_constraints/KcsA/assembly_interface/
```

The two manifest rows must look like this in meaning:

```text
target = KcsA
evidence_type = pore_filter_coupling
state_or_context = TVGYG_selectivity_filter_or_pore_filter
file_path = data/external_constraints/KcsA/pore_filter/<real file>
source_name = real database / paper / external pipeline name
source_url_or_citation = stable source or citation
```

and:

```text
target = KcsA
evidence_type = assembly_interface_constraint
state_or_context = tetramer_or_chain_interface_context
file_path = data/external_constraints/KcsA/assembly_interface/<real file>
source_name = real database / paper / external pipeline name
source_url_or_citation = stable source or citation
```

Rules:

- Do not point the manifest to `first_contact_clean_pharmacotopology_layer_run/`.
- Do not point the manifest to `data/locked_runtime_certificates/`.
- Do not use annotation-only rows as constraint-backed evidence.
- Do not pool XCL1 state A and state B unless both have separate real external constraints.
- Do not tune thresholds after seeing native/contact metrics.

## Success gate before V33

Rerun V32 after importing real files. The only acceptable transition state is:

```text
preflight_status = V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED
selected_V33_target = KcsA
claim_allowed = False
new_MD_allowed = False
new_MD_recommended = False
```

Only then should the next script be a V33 no-MD, constraint-backed operator readout. Even then, the claim remains disabled until a locked, pre-registered panel beats controls without leakage.

## Path toward the larger goal

The near-term goal is not “solve all protein folding” in one jump. The clean path is:

1. Prove real external source import for KcsA without leakage.
2. Run a no-MD V33 operator readout only after V32 selects KcsA.
3. Compare against matched/randomized/adversarial controls before any claim.
4. Only after no-MD readout passes should MD/pressure tests be considered.
5. Then repeat on a harder protein class, not before the KcsA gate is clean.
