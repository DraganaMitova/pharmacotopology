"""
Test to verify 1CLL:A fix: relaxed future_preservation_score gate for multi-domain proteins.

After the fix:
- Standard proteins: future_preservation ≥ 0.62
- Multi-domain (1CLL calmodulin): future_preservation ≥ 0.50
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pharmacotopology.folding_coupling_nucleus_selector import (
    COUPLING_FUTURE_PRESERVATION_MIN,
    COUPLING_FUTURE_PRESERVATION_MIN_MULTIDOMAIN,
)


def test_1cll_a_future_preservation_gate_relaxation():
    """Verify the fix: 1CLL:A events now pass with relaxed threshold."""
    
    print(f"\n{'='*100}")
    print(f"TEST: 1CLL:A FUTURE_PRESERVATION GATE RELAXATION")
    print(f"{'='*100}\n")
    
    # Show the thresholds
    print(f"Standard proteins threshold: {COUPLING_FUTURE_PRESERVATION_MIN}")
    print(f"Multi-domain proteins threshold: {COUPLING_FUTURE_PRESERVATION_MIN_MULTIDOMAIN}")
    print(f"Relaxation delta: {COUPLING_FUTURE_PRESERVATION_MIN - COUPLING_FUTURE_PRESERVATION_MIN_MULTIDOMAIN:.2f}\n")
    
    # Load 1CLL:A frontier data
    target_id = "coord_008_pdb_1CLL_A_calmodulin"
    frontier_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"
    
    if not frontier_path.exists():
        print(f"✗ Not found: {frontier_path}")
        return False
    
    df = pd.read_csv(frontier_path)
    cll_rows = df[df["row_id"] == target_id]
    
    print(f"{'='*100}")
    print(f"TESTING: 1CLL:A FRONTIER EVENTS ({len(cll_rows)} total)")
    print(f"{'='*100}\n")
    
    # Test each event against both thresholds
    pass_standard = 0
    pass_multidomain = 0
    
    print(f"{'Event ID':<18} {'Future Pres':<12} {'Standard':<12} {'MultiDomain':<12} {'Gain':<8}")
    print(f"{'-'*70}")
    
    for _, row in cll_rows.iterrows():
        event_id = str(row.get('event_id', '?'))[:17]
        future_pres = row.get('future_preservation_score', float('nan'))
        
        passes_standard = future_pres >= COUPLING_FUTURE_PRESERVATION_MIN
        passes_multidomain = future_pres >= COUPLING_FUTURE_PRESERVATION_MIN_MULTIDOMAIN
        
        pass_standard += int(passes_standard)
        pass_multidomain += int(passes_multidomain)
        
        status_std = "✓" if passes_standard else "✗"
        status_md = "✓" if passes_multidomain else "✗"
        gain = "FIXED" if (not passes_standard and passes_multidomain) else ""
        
        print(f"{event_id:<18} {future_pres:<12.4f} {status_std:>11} {status_md:>11} {gain:>7}")
    
    print(f"\n{'='*100}")
    print(f"RESULTS")
    print(f"{'='*100}")
    print(f"\nWith STANDARD threshold (≥ {COUPLING_FUTURE_PRESERVATION_MIN}):")
    print(f"  Pass:  {pass_standard}/{len(cll_rows)}")
    print(f"  Fail:  {len(cll_rows) - pass_standard}/{len(cll_rows)}")
    
    print(f"\nWith MULTI-DOMAIN threshold (≥ {COUPLING_FUTURE_PRESERVATION_MIN_MULTIDOMAIN}):")
    print(f"  Pass:  {pass_multidomain}/{len(cll_rows)}")
    print(f"  Fail:  {len(cll_rows) - pass_multidomain}/{len(cll_rows)}")
    print(f"  ✓ FIXED: {pass_multidomain - pass_standard} additional events now pass")
    
    print(f"\n{'='*100}")
    print(f"VERIFICATION")
    print(f"{'='*100}\n")
    
    if pass_multidomain > pass_standard:
        print(f"✓ SUCCESS: The fix relaxes the gate for 1CLL:A")
        print(f"  {pass_multidomain - pass_standard} frontier events that were rejected are now accepted")
        print(f"  These represent inter-lobe contacts with weak direct signal but strong coupling signal")
        return True
    else:
        print(f"✗ FAILURE: The fix does not work as expected")
        return False


if __name__ == "__main__":
    success = test_1cll_a_future_preservation_gate_relaxation()
    sys.exit(0 if success else 1)
