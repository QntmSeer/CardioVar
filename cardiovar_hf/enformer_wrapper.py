import torch
import numpy as np
from enformer_pytorch import Enformer
from api_integrations import fetch_genomic_sequence

# Global model instance to avoid reloading
_MODEL = None

def get_model():
    global _MODEL
    if _MODEL is None:
        print(">> Loading Enformer model (this may take a moment)...")
        try:
            # Load pre-trained Enformer
            _MODEL = Enformer.from_pretrained('EleutherAI/enformer-official-rough')
            _MODEL.eval()  # Set to evaluation mode
            print(">> Enformer model loaded successfully")
        except Exception as e:
            print(f">> Failed to load Enformer: {e}")
            return None
    return _MODEL

# ... (one_hot_encode function remains same) ...

def predict_variant_impact_dl(chrom, pos, ref, alt):
    # ... (args docstring) ...
    model = get_model()
    if model is None:
        return None
        
    # Enformer requires 196,608 bp context
    SEQUENCE_LENGTH = 196_608
    start = pos - (SEQUENCE_LENGTH // 2)
    end = start + SEQUENCE_LENGTH
    
    # 1. Fetch Reference Sequence
    print(f">> Fetching {SEQUENCE_LENGTH}bp sequence for {chrom}:{pos}...")
    ref_seq = fetch_genomic_sequence(chrom, start, end)
    
    if not ref_seq or len(ref_seq) != SEQUENCE_LENGTH:
        print(f">> Failed to fetch correct sequence length. Got {len(ref_seq) if ref_seq else 0}")
        return None
        
    # 2. Create Alternate Sequence
    # ... (logic remains same) ...
    
    # Verify reference matches
    fetched_ref = ref_seq[rel_pos:rel_pos+len(ref)]
    if fetched_ref != ref:
        print(f">> Reference mismatch! Expected {ref}, got {fetched_ref}")
        # Continue anyway for demo, but warn
    
    # ... (alt seq construction) ...
    if len(ref) == len(alt):
        alt_seq = ref_seq[:rel_pos] + alt + ref_seq[rel_pos+len(ref):]
    else:
        print(">> Indels not fully supported in this demo version")
        return None

    # ... (encoding logic) ...
    
    # 4. Run Prediction
    print(">> Running Enformer inference...")
    with torch.no_grad():
        ref_pred = model(ref_tensor)
        alt_pred = model(alt_tensor)
        
    # Enformer returns dictionary with 'human' and 'mouse' heads
    # We want 'human' head
    # Shape: (batch, seq_len, tracks) -> (1, 896, 5313)
    
    ref_human = ref_pred['human'].numpy()[0]
    alt_human = alt_pred['human'].numpy()[0]
    
    # 5. Calculate Impact (L1 norm of difference across all tracks)
    # This gives us a profile of change across the 896 bins (each bin ~128bp)
    delta_profile = np.abs(alt_human - ref_human).mean(axis=1) # Mean across tracks
    
    # Center is the variant position
    center_idx = len(delta_profile) // 2
    
    # Extract a window around the center for visualization
    # The model output is lower resolution (128bp bins)
    # We'll interpolate to get a smooth curve for the UI
    
    return {
        "raw_delta": delta_profile,
        "center_idx": center_idx,
        "max_impact": np.max(delta_profile)
    }
