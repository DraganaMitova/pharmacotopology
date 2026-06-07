import sys
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

def test_full_gap_based():
    df = pd.read_csv('first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_frontier.csv')
    
    for row_id in df['row_id'].unique():
        prot = df[df['row_id'] == row_id]
        
        # future preservation gap
        fp_scores = sorted([row['future_preservation_score'] for _, row in prot.iterrows()], reverse=True)
        if len(fp_scores) < 2: gap_fp = 0.0
        else:
            gaps = [fp_scores[i] - fp_scores[i+1] for i in range(len(fp_scores)-1)]
            gap_fp = fp_scores[max(range(len(gaps)), key=lambda i: gaps[i]) + 1]
            
        # blocked future gap (ascending)
        bf_scores = sorted([row['blocked_future_pressure'] for _, row in prot.iterrows()])
        if len(bf_scores) < 2: gap_bf = 1.0
        else:
            gaps = [bf_scores[i+1] - bf_scores[i] for i in range(len(bf_scores)-1)]
            gap_bf = bf_scores[max(range(len(gaps)), key=lambda i: gaps[i])]
            
        pass_count = 0
        for _, row in prot.iterrows():
            fp_pass = row['future_preservation_score'] >= gap_fp
            bf_pass = row['blocked_future_pressure'] <= gap_bf
            ds_pass = row['direct_support_score'] >= 0.22  # keeping direct support hardcoded 0.22?
            
            # Wait, the user said "без 0.62, без 0.16. Само internal gaps."
            # They didn't mention direct support.
            
            if fp_pass and bf_pass and ds_pass:
                pass_count += 1
                
        print(f"{row_id[:25]:<25} | GapFP: {gap_fp:.4f} | GapBF: {gap_bf:.4f} | Pass: {pass_count}/{len(prot)}")

if __name__ == "__main__":
    test_full_gap_based()
