"""
Surgical fix for 1CLL:A collapse: Lower future_preservation_score gate
for calmodulin's inter-lobe contacts.

The diagnostic revealed:
- frontier_events = 7 (high-quality replacement candidates)
- selected_events = 671 (but many duplicates)
- future_preservation_score = 0.4935-0.5667 (mean 0.5452)
- GATE THRESHOLD = 0.62 (too strict!)

Hypothesis: Calmodulin's inter-lobe contacts have weak *direct* evolutionary
signal in FUTURE positions (because lobes move independently), but strong
*coupling* signal from multi-body coevolution.

Solution: For proteins like calmodulin (multi-domain, inter-lobe contacts),
use relaxed threshold: future_preservation ≥ 0.50 instead of 0.62

Test this on 1CLL:A only.
"""

from pathlib import Path
import sys
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def analyze_1cll_a_future_scores():
    """Dump all 1CLL:A frontier events with future_preservation_score."""
    
    target_id = "coord_008_pdb_1CLL_A_calmodulin"
    frontier_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"
    
    if not frontier_path.exists():
        print(f"✗ Not found: {frontier_path}")
        return
    
    df = pd.read_csv(frontier_path)
    cll_rows = df[df["row_id"] == target_id]
    
    print(f"\n{'='*100}")
    print(f"1CLL:A FRONTIER EVENTS - FUTURE PRESERVATION ANALYSIS")
    print(f"{'='*100}")
    print(f"Total frontier events: {len(cll_rows)}\n")
    
    # Show all events with key scoring columns
    cols_to_show = [
        "frontier_kind", "source_selector", "event_id",
        "segment_a_start", "segment_a_end",
        "segment_b_start", "segment_b_end",
        "future_preservation_score",  # ← KEY COLUMN
        "blocked_future_pressure",
        "direct_support_score",
        "contact_cluster_gain",
        "coupling_nucleus_score",
        "exclusion_reasons"
    ]
    
    existing_cols = [c for c in cols_to_show if c in cll_rows.columns]
    
    # Display with formatting
    for idx, (i, row) in enumerate(cll_rows.iterrows(), 1):
        print(f"Event {idx}:")
        print(f"  ID: {row.get('event_id', '?')}")
        print(f"  Selector: {row.get('source_selector', '?')} → {row.get('target_selector', '?')}")
        print(f"  Segment: [{int(row['segment_a_start'])},{int(row['segment_a_end'])}] -- [{int(row['segment_b_start'])},{int(row['segment_b_end'])}]")
        print(f"  FUTURE_PRESERVATION: {row.get('future_preservation_score', float('nan')):.4f} ⚠️ (gate: ≥0.62)")
        print(f"  blocked_future: {row.get('blocked_future_pressure', float('nan')):.4f} (gate: ≤0.16)")
        print(f"  direct_support: {row.get('direct_support_score', float('nan')):.4f} (gate: ≥0.22)")
        print(f"  contact_cluster_gain: {row.get('contact_cluster_gain', float('nan')):.4f}")
        print(f"  coupling_nucleus_score: {row.get('coupling_nucleus_score', float('nan')):.4f}")
        
        reasons = row.get('exclusion_reasons', 'SELECTED')
        if pd.isna(reasons) or reasons == '':
            print(f"  Status: ✓ SELECTED")
        else:
            print(f"  Status: ✗ REJECTED")
            print(f"  Reasons: {reasons}")
        print()
    
    print(f"{'='*100}")
    print(f"ANALYSIS")
    print(f"{'='*100}")
    
    # Check how many would pass with relaxed threshold
    below_062 = (cll_rows["future_preservation_score"] < 0.62).sum()
    above_050 = (cll_rows["future_preservation_score"] >= 0.50).sum()
    
    print(f"\nWith current gate (future_preservation ≥ 0.62):")
    print(f"  Events passing: {len(cll_rows) - below_062}/{len(cll_rows)}")
    print(f"  Events failing: {below_062}/{len(cll_rows)}")
    
    print(f"\nWith relaxed gate (future_preservation ≥ 0.50):")
    print(f"  Events passing: {above_050}/{len(cll_rows)}")
    print(f"  Events failing: {len(cll_rows) - above_050}/{len(cll_rows)}")
    
    print(f"\n{'='*100}")
    print(f"RECOMMENDATION")
    print(f"{'='*100}")
    print("""
For multi-domain proteins like CALMODULIN (1CLL:A):
- Inter-lobe contacts have WEAK direct evolutionary signal in future positions
- But STRONG coupling signal from allosteric coevolution
- Current threshold (≥0.62) filters out these allosteric contacts

FIX: Create protein-specific gate thresholds:
  - Standard proteins: future_preservation ≥ 0.62
  - Multi-domain (CalM, etc): future_preservation ≥ 0.50
  - OR: Use ensemble approach combining direct + coupling signals

Test: Temporarily lower threshold to 0.50 for 1CLL:A and measure:
  1. How many additional events get selected?
  2. Do they correspond to missing native long-range contacts?
  3. Does false_nucleus_rate stay at 0.0?
    """)


if __name__ == "__main__":
    analyze_1cll_a_future_scores()
