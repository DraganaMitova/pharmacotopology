# COMMANDS_RUN_NOW — V32 → V33 clean continuation

Use these commands from Terminal on macOS. No git commands. No full-suite runs. No MD.

## 0) Unzip and enter the clean project

If the ZIP is in Downloads:

```bash
cd "$HOME/Downloads" || exit 1
unzip -q pharmacotopology_clean_v33_start_COMMANDS.zip -d "$HOME/Desktop"
cd "$HOME/Desktop/pharmacotopology_clean_v33_start" || exit 1
```

If Safari renamed it, check the exact name first:

```bash
ls -lh "$HOME/Downloads" | grep pharmacotopology
```

## 1) First sanity test — expected to PASS and ABSTAIN

```bash
python3 scripts/run_clean_start_sanity_v0.py
```

Expected scientific result:

```text
13 passed
preflight_status = V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED
selected_V33_target = None
claim_allowed = False
new_MD_allowed = False
```

This is correct. It means the system refuses placeholders/internal evidence.

## 2) Manual V32 preflight rerun

```bash
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
```

Expected before real KcsA source import:

```text
V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED
```

## 3) Only targeted V30–V32 tests

```bash
python3 -m pytest -q \
  tests/test_v30_external_constraint_and_coupling_acquisition_sprint.py \
  tests/test_v31_constraint_backed_operator_readout_preflight.py \
  tests/test_v32_external_constraint_source_import_preflight.py
```

Expected:

```text
13 passed
```

## 4) Prepare KcsA real external source import

Do not use generated runtime reports as evidence. Put real external KcsA evidence files here:

```bash
mkdir -p data/external_constraints/KcsA/pore_filter
mkdir -p data/external_constraints/KcsA/assembly_interface
open data/external_constraints/KcsA
open -e data/external_constraints/v32_external_constraint_source_import_manifest.json
```

Edit the manifest so KcsA rows point to real files, not `<real_external...>` placeholders.

Good target state:

```text
KcsA pore/filter evidence exists
KcsA assembly/interface evidence exists
source_name is real database/paper/pipeline
source_url_or_citation is real
file_path exists
not under first_contact_clean_pharmacotopology_layer_run
not under data/locked_runtime_certificates
claim_allowed stays false
```

## 5) After real KcsA import, rerun the gate

```bash
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
```

Target result after real external KcsA source files are imported:

```text
preflight_status = V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED
selected_V33_target = KcsA
claim_allowed = False
new_MD_allowed = False
new_MD_recommended = False
```

## 6) Do not run yet

Do not run V33 operator readout until step 5 gives KcsA selected.
Do not run full suite.
Do not run MD.
Do not claim protein folding solved.

The next scientific win is not a folding claim. It is this:

```text
Real external KcsA pore/filter + assembly/interface evidence opens V33 while placeholders, annotations, and internal reports stay blocked.
```
