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
            
            # Move to GPU if available
            if torch.cuda.is_available():
                _MODEL = _MODEL.cuda()
                print(f">> Using GPU: {torch.cuda.get_device_name(0)}")
            else:
                print(">> Using CPU (GPU not available)")
            
            _MODEL.eval()  # Set to evaluation mode
            print(">> Enformer model loaded successfully")
        except Exception as e:
            print(f">> Failed to load Enformer: {e}")
            return None
    return _MODEL



def one_hot_encode(seq):
    """Convert DNA sequence to one-hot encoding for Enformer."""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    seq = seq.upper()
    
    # Create one-hot matrix
    one_hot = np.zeros((len(seq), 4), dtype=np.float32)
    for i, nucleotide in enumerate(seq):
        if nucleotide in mapping and mapping[nucleotide] < 4:
            one_hot[i, mapping[nucleotide]] = 1.0
        # N or unknown nucleotides get all zeros
    
    return one_hot


def predict_variant_impact_dl(chrom, pos, ref, alt):
    """
    Predict variant impact using Enformer deep learning model.
    
    Args:
        chrom: Chromosome (e.g., "chr22")
        pos: Genomic position
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        Dict with raw_delta, center_idx, max_impact or None if failed
    """
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
    # Calculate relative position of variant within the fetched sequence
    rel_pos = SEQUENCE_LENGTH // 2  # Variant is at center of sequence
    
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


    # 3. One-hot encode sequences
    ref_encoded = one_hot_encode(ref_seq)
    alt_encoded = one_hot_encode(alt_seq)
    
    # Convert to torch tensors and add batch dimension
    ref_tensor = torch.from_numpy(ref_encoded).unsqueeze(0)  # (1, seq_len, 4)
    alt_tensor = torch.from_numpy(alt_encoded).unsqueeze(0)
    
    # Move tensors to same device as model (CPU or GPU)
    device = next(model.parameters()).device
    ref_tensor = ref_tensor.to(device)
    alt_tensor = alt_tensor.to(device)
    
    # 4. Run Prediction
    print(">> Running Enformer inference...")
    with torch.no_grad():
        ref_pred = model(ref_tensor)
        alt_pred = model(alt_tensor)
        
    # Enformer returns dictionary with 'human' and 'mouse' heads
    # We want 'human' head
    # Shape: (batch, seq_len, tracks) -> (1, 896, 5313)
    
    # Move to CPU for numpy operations
    ref_human = ref_pred['human'].cpu().numpy()[0]
    alt_human = alt_pred['human'].cpu().numpy()[0]
    
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
