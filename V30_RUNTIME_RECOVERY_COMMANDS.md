# V30 Runtime Recovery

Use this after accidentally deleting/replacing `first_contact_clean_pharmacotopology_layer_run`.

This recovery restores small locked V15 summary certificates from prior locked outputs, then rebuilds V15→V30 runtime artifacts in order. It does **not** restore raw MD trajectories and does **not** upgrade any claims.

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

bash "$REPO_ROOT/scripts/rebuild_runtime_chain_v15_to_v30_after_recovery.sh"
```

Expected end status:

```text
V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED
claim_allowed=false
new_MD_allowed=false
positive_folding_evidence_targets=[]
```

Important: this recovery uses `scripts/recover_v15_locked_source_snapshots_v0.py`, which writes small claim-disabled summary certificates for 1UBQ, 1CLL, and 4AKE from previously locked outputs. It is not a raw trajectory restore.
MD rerun is still forbidden unless raw trajectories are explicitly required later.
