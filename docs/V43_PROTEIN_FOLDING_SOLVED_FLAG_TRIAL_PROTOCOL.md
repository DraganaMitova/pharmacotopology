# V43 Protein Folding Solved Flag Trial Protocol

## Purpose

V43 is the first solved-candidate trial. It is allowed to set:

`protein_folding_solved_candidate = true`

only if all hard thresholds pass.

It must never set:

`folding_problem_solved = true`

The absolute solved flag remains false pending external blind review, broad benchmarking, and independent replication.

## Panel

V43 uses a minimum 24-target panel, imported from the V42 de novo mechanism-language challenge and augmented with structural, ensemble, function, and contact-scoring requirements.

Required target classes:

- compact soluble folded proteins
- membrane proteins/channels/transporters
- weak or shallow evolutionary-information proteins
- disordered or partially disordered proteins
- metamorphic/multistate/allosteric proteins
- folding-upon-binding or contextual-ordering proteins

## Sealing Rule

All prediction packets are written and hashed before post-seal holdout validation is opened.

Prediction packets include:

- `prediction_timestamp`
- `prediction_hash`
- `prediction_inputs_manifest`
- `blocked_inputs_manifest`
- `no_holdout_access_before_hash: true`
- `mechanism_class`
- `operator_regions`
- `predicted_contact_or_region_constraints`
- `predicted_low_resolution_structure_or_ensemble`
- `predicted_perturbation_effects`
- `falsification_criteria`

## Blocked Inputs Before Sealing

- PDB/mmCIF coordinates
- PDB-derived contacts
- native contact maps
- AlphaFold/ESMFold/RoseTTAFold predicted coordinates
- V33/V34 coordinate-derived files
- runtime reports as biological evidence
- holdout files
- answer keys

## Post-Seal Holdouts

After sealing only, V43 may use:

- structure/contact scoring summaries
- contact maps derived after sealing for scoring only
- experimental disorder/ensemble annotations
- state-specific function evidence
- oligomer/interface/ligand annotations
- mutagenesis/function annotations
- literature validation sources

## Solved-Candidate Thresholds

`protein_folding_solved_candidate` can be true only if all thresholds pass:

- prediction sealed before holdout
- panel target count at least 24
- holdout validated target count at least 20
- mechanism accuracy at least 0.80
- hard-class precision at least 0.80
- hard-class recall at least 0.75
- operator support at least 0.65
- low-resolution structure/ensemble support at least 0.65
- perturbation support at least 0.60
- contact precision at least 0.40 or contact enrichment at least 2.0
- false hard grammar promotions at most 2
- false single-fold promotions at most 2
- beats random, keyword, majority, and simple sequence baselines
- no coordinate/internal/holdout leakage before prediction
- no native metrics before prediction
- all controls pass
- failures are documented

## Safe Pass Status

`V43_PROTEIN_FOLDING_SOLVED_CANDIDATE_PASSED_REVIEW_REQUIRED`

This status allows solved-candidate wording only. It does not allow the universal solved claim.
