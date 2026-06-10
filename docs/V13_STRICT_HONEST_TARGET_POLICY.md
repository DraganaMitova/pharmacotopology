# V13 strict honest target policy

Official V13 transfer/provenance tests must use only real downloaded target structures that contain an exact full-length benchmark sequence segment.

Disallowed for official results:

- synthetic/generated target coordinates
- native/debug fallback targets
- partial-prefix or truncated targets
- selector/role-rule tuning during target repair

For the 1UBQ transfer repair lane, the default independent target candidates are RCSB `1D3Z` and `5DK8`, both treated as external real ubiquitin structures. The target preparer accepts a candidate only if the full 76-residue benchmark sequence can be found and renumbered exactly.

For 1CLL, the default target source is AFDB only. The script intentionally does not fall back to the native RCSB 1CLL structure.

All generated certificates must keep `claim_allowed=false` until the run has complete runtime provenance and passes the locked interpretation gates.
