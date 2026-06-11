# Pharmacotopology

Mechanism-language experiments for hard protein-folding regimes.

This repository tests a modest but useful idea: before claiming structure-level success, can a system produce sealed, falsifiable mechanism packets from non-coordinate evidence, then survive post-seal holdout validation?

It is not a universal protein-folding solver. It is not an AlphaFold replacement. It is a small, audit-heavy research workspace for asking whether operator-level biological explanations can be made reproducible, source-separated, and hard to overclaim.

## Current Frontier

Recent live-target packets now cover three different hard classes:

| Version | Target | Class | Status |
| --- | --- | --- | --- |
| V44 | FUS-LC residues 1-214 | IDP / phase separation | sealed packet, review-required |
| V45 | TDP-43 LCD residues 274-414 | IDP / phase separation replication | sealed packet, review-required |
| V46 | CFTR F508del | membrane multidomain folding / proteostasis | committed and synced |
| V47 | RfaH-CTD residues 101-162 | metamorphic alpha/beta fold switching | committed and synced |
| V48 | SARS-CoV-2 ORF6 | viral accessory host hijacking | generated, review-required |

Every packet keeps:

```text
folding_problem_solved = false
```

The strongest current claim is narrower:

```text
sealed mechanism-language packets can generate and validate operator-level predictions
across several hard folding regimes without coordinate leakage before prediction
```

## How The Runs Work

Each live-target attack follows the same discipline:

1. Acquire allowed non-coordinate sources.
2. Write a sealed prediction packet and hash it.
3. Keep coordinates, native contacts, AlphaFold-style models, internal runtime artifacts, and holdouts out of prediction.
4. Open independent holdout evidence only after sealing.
5. Score operator buckets, perturbation predictions, leakage controls, and forbidden-claim guards.
6. Emit a certificate and a short report.

The useful output is not a picture of a fold. It is a claim-controlled packet: what mechanism is proposed, what would falsify it, which shortcuts were blocked, and which post-seal evidence supported or contradicted it.

## Evidence Boundary

V49 makes the source boundary explicit. Evidence is not just "coordinate" or "non-coordinate"; every source must declare which class and role it plays.

- `pure_non_coordinate`: sequence, UniProt features, Pfam/InterPro, DisProt, PTM annotations, function annotations, motif labels, composition, charge, hydrophobicity, low-complexity, and secondary-structure propensity.
- `spatial_proxy_non_coordinate`: DCA/evolutionary couplings, FRET, crosslinking, NMR chemical shifts, SAXS, EPR/DEER, distance restraints, and contact-like biochemical constraints.
- `coordinate_derived`: PDB/mmCIF coordinates, contact maps from structures, AlphaFold/ESMFold/RoseTTAFold coordinates, native distance matrices, structural alignments, and PDB-derived interface contacts.
- `internal_runtime`: certificates, generated summaries, previous prediction outputs, and other `first_contact_clean_pharmacotopology_layer_run` artifacts.

The sharp boundary is between `pure_non_coordinate` and `spatial_proxy_non_coordinate`: spatial-proxy evidence can encode spatial information, so it must never be hidden inside pure non-coordinate evidence. Spatial-proxy sources must be explicitly tagged and assigned one role: `prediction_input`, `holdout_validation`, `baseline_only`, or `blocked`.

Coordinate-derived sources are blocked before sealing. They may appear only after sealing as `holdout_validation` sources. Internal runtime artifacts are audit metadata only and are never biological prediction evidence.

All new source rows must include `source_class`, `source_role`, `spatial_proxy`, `coordinate_derived`, `internal_runtime`, `allowed_for_prediction`, `allowed_for_holdout`, `allowed_for_claim`, `leakage_risk`, and `rationale`.

## Operator Scoring

Operator-level scoring is frozen as a pre-registered row schema, not a moving narrative. Each prediction row must include:

```text
prediction_id
target
mechanism_class
operator_bucket
region_or_state
predicted_effect
perturbation
expected_direction
confidence
prediction_source_ids
falsification_criteria
holdout_evidence_ids
score_label
score_reason
scoring_pre_registered
```

Allowed `score_label` values are:

```text
supported
partially_supported
contradicted
not_testable
blocked_for_leakage
```

For CFTR F508del, this means the packet cannot simply say "validated." It must map each operator-level claim to post-seal evidence: NBD1 stability, NBD1-only partial rescue, interdomain/interface correction, trafficking/proteostasis context, and forbidden shortcut rejection.

## Latest Completed Full Cycle

V46 completed seal -> holdout -> validation as the latest README-visible full cycle.

- Sealed packet: `data/live_unsolved_targets/V46/CFTR_F508del/prediction_sealed/sealed_prediction_packet.json`
- Post-seal holdout manifest: `data/live_unsolved_targets/V46/CFTR_F508del/holdouts_postseal/holdout_manifest.json`
- Validation result: `data/live_unsolved_targets/V46/CFTR_F508del/validation/validation_result.json`
- Certificate: `first_contact_clean_pharmacotopology_layer_run/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK/v46_cftr_f508del_membrane_multidomain_attack_certificate.json`
- Report: `first_contact_clean_pharmacotopology_layer_run/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK/V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FOLDING_RESCUE_ATTACK_REPORT.md`
- Validated predictions: `12`
- Contradicted predictions: `0`
- Controls passed: `17 / 17`
- Leakage counts before prediction: coordinate-derived `0`, internal runtime `0`, holdout leakage `False`

V46's supported claims remain bounded: F508del weakens NBD1 stability, NBD1-only rescue is partial, interdomain/interface correction is needed, and trafficking/proteostasis is maturation context rather than atomic proof.

## Negative / Null Controls

Future live packets must expose these controls in the README-visible protocol, certificate, or both:

- `random_sequence_control`: random sequence input must not create a target-specific packet.
- `shuffled_sequence_control`: shuffled sequence input must weaken or block operator assignment.
- `swapped_evidence_control`: evidence from another target must not validate target-specific predictions.
- `wrong_target_control`: wrong target labels must not pass as the live target.
- `generic_annotation_only_control`: generic function annotation alone must not create a full packet.
- `coordinate_leakage_control`: coordinate-derived sources before sealing must block prediction use.
- `internal_runtime_leakage_control`: internal reports and certificates must not become biological evidence.
- `forced_wrong_grammar_control`: forcing the wrong mechanism grammar must fail or cleanly abstain.
- `failed_prediction_not_repaired_after_holdout`: failed predictions must remain failed after holdout opening.
- `holdout_opened_before_seal_control`: holdout evidence opened before sealing must block the run.

## Known Weaknesses / Reviewer Attack Surface

The next hardening layer is V49 reviewer attack-surface hardening. It does not add a target and does not strengthen the biological claim. It makes four reviewer questions mechanically auditable:

- What counts as non-coordinate evidence?
- How are operator-level predictions scored?
- Where is a complete seal -> holdout -> report cycle visible?
- Which negative/null controls are required?

## Reproduce The Latest Packets

From the repository root:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/build_v46_cftr_f508del_membrane_multidomain_sources_v0.py
python3 scripts/run_v46_cftr_f508del_membrane_multidomain_attack_v0.py
python3 -m pytest -q tests/test_v46_cftr_f508del_membrane_multidomain_attack.py
```

Run the reviewer hardening audit:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/run_v49_reviewer_attack_surface_audit_v0.py
python3 scripts/print_v49_reviewer_attack_surface_audit.py
python3 -m pytest -q tests/test_v49_reviewer_attack_surface_audit.py
```

## Repository Map

- `scripts/` builds source manifests, sealed packets, validation artifacts, and printed reports.
- `tests/` contains focused regression tests for each packet's scientific and leakage gates.
- `docs/` records the protocol and claim boundary for each version.
- `data/live_unsolved_targets/` stores source manifests, sealed predictions, post-seal holdouts, and validation JSON.
- `first_contact_clean_pharmacotopology_layer_run/` stores run certificates and human-readable reports.

## Scientific Boundary

Useful claims are welcome here. Grand claims are not.

Allowed:

```text
This run produced a sealed, source-separated mechanism packet that passed its
defined holdout and leakage controls.
```

Not allowed:

```text
Universal protein folding is solved.
Coordinates were predicted de novo.
External review is unnecessary.
```

The work is most interesting when it stays humble: each packet should make the next experiment clearer, not make the conclusion louder.
