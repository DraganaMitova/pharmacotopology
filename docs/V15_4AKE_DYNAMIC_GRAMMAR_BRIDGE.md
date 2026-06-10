# V15_4AKE_DYNAMIC_GRAMMAR_BRIDGE

This layer creates a claim-disabled 4AKE bridge artifact for the V15 dynamic grammar panel.

It does **not** synthesize evidence from PDB-only or visual-only files. If no machine-readable 4AKE role-aware steering artifact exists, the bridge records:

```text
4AKE_legacy_visual_input_material_present_but_machine_readable_role_artifact_missing_claim_disabled
```

If a machine-readable role-aware 4AKE artifact is later restored or produced, the bridge can normalize selected strict/balanced/hinge/rescue pairs into the dynamic grammar row.
