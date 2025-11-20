import numpy as np

def compute_variant_impact(chrom, pos, ref, alt, window_size=100):
    """
    Core logic for variant impact prediction.
    Returns synthetic but realistic data structure.
    """
    # 1. Generate Curve Data
    x = np.arange(-window_size, window_size + 1)
    np.random.seed(pos % 10000) # Deterministic seed based on pos
    
    # Signal shape
    direction = 1 if (pos % 2 == 0) else -1
    signal = 3.5 * direction * np.exp(-0.02 * (x)**2)
    noise = np.random.normal(0, 0.3, len(x))
    delta_rna = signal + noise
    
    # 2. Calculate Metrics
    max_idx = np.argmax(np.abs(delta_rna))
    max_delta = float(delta_rna[max_idx])
    max_pos_rel = int(x[max_idx])
    
    # 3. Mock gnomAD Frequency
    freq = 0.00005 if ref == 'A' else 0.05
    
    # 4. Mock Gene Symbol
    gene_symbol = "MYH9" if chrom == "chr22" else "PCSK9"
    
    # 5. Mock Tracks
    # Gene Structure (Exons)
    exons = [{"start": -50, "end": 20, "label": "Exon 1"}]
    
    # Conservation (PhyloP)
    cons_scores = np.random.normal(0.5, 1.0, len(x))
    cons_scores[window_size-10:window_size+10] += 2.0 # Conserved peak
    
    # 6. Tissue-Specific Effects (mock)
    tissues = ["Heart LV", "Heart RA", "Aorta", "Coronary Artery", "Liver", "Brain", "Kidney"]
    tissue_effects = []
    for tissue in tissues:
        # Heart tissues have higher impact
        if "Heart" in tissue or "Aorta" in tissue or "Coronary" in tissue:
            effect = abs(max_delta) * np.random.uniform(0.7, 1.2)
        else:
            effect = abs(max_delta) * np.random.uniform(0.1, 0.4)
        tissue_effects.append({"tissue": tissue, "delta": round(effect, 2)})
    
    # 7. Background Distribution (mock - simulate other variants in this gene)
    np.random.seed(pos % 100)
    background_deltas = np.random.normal(0, 1.5, 200)  # 200 background variants
    background_deltas = background_deltas.tolist()
    
    # Calculate percentile
    all_deltas = background_deltas + [abs(max_delta)]
    percentile = (np.sum(np.array(all_deltas) < abs(max_delta)) / len(all_deltas)) * 100
    
    return {
        "variant_id": f"{chrom}:{pos}:{ref}:{alt}",
        "metrics": {
            "max_delta": round(max_delta, 4),
            "max_pos_rel": max_pos_rel,
            "gnomad_freq": freq,
            "gene_symbol": gene_symbol,
            "percentile": round(percentile, 1)
        },
        "curve": {
            "x": x.tolist(),
            "y": delta_rna.tolist()
        },
        "tracks": {
            "exons": exons,
            "conservation": cons_scores.tolist()
        },
        "tissue_effects": tissue_effects,
        "background_distribution": {
            "background_deltas": background_deltas,
            "variant_delta": abs(max_delta)
        }
    }
