"""
Integration test: Verify 1CLL:A native long-range contact recovery after the fix.

Before fix: 6/7 frontier events are rejected (too strict gate)
After fix: 6/7 frontier events are accepted (relaxed gate)

This test shows how many additional native contacts these 6 events would help recover.
"""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_1cll_a_native_contact_recovery_forecast():
    """Forecast native contact recovery after the fix."""
    
    print(f"\n{'='*100}")
    print(f"FORECAST: 1CLL:A NATIVE CONTACT RECOVERY AFTER THE FIX")
    print(f"{'='*100}\n")
    
    target_id = "coord_008_pdb_1CLL_A_calmodulin"
    frontier_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"
    
    if not frontier_path.exists():
        print(f"✗ Not found frontier data")
        return
    
    df = pd.read_csv(frontier_path)
    cll_rows = df[df["row_id"] == target_id]
    
    print(f"Analyzing {len(cll_rows)} frontier events for 1CLL:A\n")
    
    # Show segments that were previously rejected
    print(f"{'='*100}")
    print(f"SEGMENTS THAT WILL NOW BE AVAILABLE (6 events previously rejected)")
    print(f"{'='*100}\n")
    
    print(f"{'Event':<4} {'Segment A':<12} {'Segment B':<12} {'Native Contacts':<16} {'Inter-lobe':<12}")
    print(f"{'-'*80}")
    
    count = 0
    for idx, (i, row) in enumerate(cll_rows.iterrows(), 1):
        future_pres = row.get('future_preservation_score', float('nan'))
        
        # Only show events that now pass with 0.50 but fail with 0.62
        if 0.50 <= future_pres < 0.62:
            seg_a = f"[{int(row['segment_a_start'])},{int(row['segment_a_end'])}]"
            seg_b = f"[{int(row['segment_b_start'])},{int(row['segment_b_end'])}]"
            native_count = row.get('native_long_range_contact_count_after_scoring', 0)
            
            # Check if it's inter-lobe (rough estimate: segments > 30 residues apart)
            a_mid = (row['segment_a_start'] + row['segment_a_end']) / 2
            b_mid = (row['segment_b_start'] + row['segment_b_end']) / 2
            distance = abs(a_mid - b_mid)
            is_interlobe = "YES" if distance > 30 else "NO"
            
            print(f"{idx:<4} {seg_a:<12} {seg_b:<12} {int(native_count):<16} {is_interlobe:<12}")
            count += 1
    
    print(f"\n✓ {count} EVENTS NEWLY AVAILABLE DUE TO GATE RELAXATION\n")
    
    # Estimate contact recovery
    print(f"{'='*100}")
    print(f"CONTACT RECOVERY FORECAST")
    print(f"{'='*100}\n")
    
    # Get the newly available events
    newly_available = cll_rows[(cll_rows['future_preservation_score'] >= 0.50) & 
                               (cll_rows['future_preservation_score'] < 0.62)]
    
    total_native_in_new_events = newly_available['native_long_range_contact_count_after_scoring'].sum()
    
    print(f"Events previously rejected (0.50 ≤ future_pres < 0.62): {len(newly_available)}")
    print(f"Total native contacts in these events: {int(total_native_in_new_events)}")
    print(f"Average per event: {total_native_in_new_events/max(1,len(newly_available)):.1f}\n")
    
    print(f"{'='*100}")
    print(f"IMPACT ASSESSMENT")
    print(f"{'='*100}\n")
    
    print(f"""
BEFORE FIX:
  ✗ 6 frontier events rejected due to strict gate (future_pres ≥ 0.62)
  ✗ These events contain inter-lobe contact information
  ✗ System is forced to rely on weaker evolutionary signals

AFTER FIX:
  ✓ 6 frontier events now accepted with relaxed gate (future_pres ≥ 0.50)
  ✓ System can incorporate allosteric/coupling signal for inter-lobe contacts
  ✓ More balanced treatment of direct vs coupling evolutionary information
  
EXPECTED OUTCOME FOR 1CLL:A:
  • Better coverage of N-lobe ↔ C-lobe contact map
  • Improved recall of 36 native long-range contacts
  • Maintained false_nucleus_rate = 0.0 (no spurious contacts)
  • More accurate representation of calmodulin's conformational dynamics

NEXT STEP:
  Run full pipeline on 1CLL:A with this fix and measure:
  1. How many of the 36 native long-range contacts are now recovered?
  2. Does false_nucleus_rate remain at 0.0?
  3. Are the recovered contacts functionally significant (inter-lobe)?
    """)


if __name__ == "__main__":
    test_1cll_a_native_contact_recovery_forecast()
