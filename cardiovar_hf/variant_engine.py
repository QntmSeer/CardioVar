import numpy as np

def compute_variant_impact(chrom, pos, ref, alt, assembly="GRCh38", window_size=100):
    """
    Core logic for variant impact prediction.
    Returns synthetic but realistic data structure.
    
    Args:
        chrom: Chromosome (e.g., "chr22")
        pos: Position
        ref: Reference allele
        alt: Alternate allele
        assembly: Genome build ("GRCh38" or "GRCh37")
        window_size: Window size around variant
    
    Raises:
        ValueError: If assembly is not GRCh38
    """
    # Validate assembly
    if assembly not in ["GRCh38", "GRCh37"]:
        raise ValueError(f"Invalid assembly: {assembly}. Must be 'GRCh38' or 'GRCh37'")
    
    if assembly == "GRCh37":
        raise ValueError("Currently only GRCh38 coordinates are supported in this demo. Please switch to GRCh38.")
    
    # 1. Gene Symbol Mapping (based on position ranges for GRCh38)
    gene_map = {
        "chr1": [
            (156100000, 156200000, "LMNA"),      # LMNA region
            (55000000, 55100000, "PCSK9"),       # PCSK9 region
            (236700000, 236750000, "ACTN2"),     # ACTN2 region
        ],
        "chr2": [
            (21000000, 21100000, "APOB"),        # APOB region
            (178500000, 178600000, "TTN"),       # TTN region
        ],
        "chr3": [
            (46850000, 46900000, "MYL3"),        # MYL3 region
        ],
        "chr22": [
            (36100000, 36400000, "MYH9"),        # MYH9 region
        ]
    }
    
    # Find gene symbol based on position
    gene_symbol = "UNKNOWN"
    if chrom in gene_map:
        for start, end, gene in gene_map[chrom]:
            if start <= pos <= end:
                gene_symbol = gene
                break
    
    # Fallback to chromosome-based if no match
    if gene_symbol == "UNKNOWN":
        gene_symbol = {"chr22": "MYH9", "chr1": "PCSK9", "chr2": "APOB", "chr3": "MYL3"}.get(chrom, "GENE_X")
    
    # 2. Generate Curve Data with Realistic Model
    x = np.arange(-window_size, window_size + 1)
    np.random.seed(pos % 10000)  # Deterministic seed based on pos
    
    # Determine variant type based on position and alleles
    is_transition = (ref in ['A', 'G'] and alt in ['A', 'G']) or (ref in ['C', 'T'] and alt in ['C', 'T'])
    is_splice_region = (pos % 100) < 10  # Mock: ~10% are near splice sites
    is_regulatory = (pos % 50) < 5  # Mock: ~10% are in regulatory regions
    
    # Base effect magnitude (depends on variant type)
    if is_splice_region:
        # Splice site variants: sharp, localized effect
        base_magnitude = np.random.uniform(3.5, 5.5)
        spread = 15  # Narrow effect
        shape = 'sharp'
    elif is_regulatory:
        # Regulatory variants: broader, moderate effect
        base_magnitude = np.random.uniform(2.0, 4.0)
        spread = 40  # Broad effect
        shape = 'broad'
    else:
        # Coding variants: moderate, intermediate spread
        base_magnitude = np.random.uniform(1.5, 3.5)
        spread = 25
        shape = 'moderate'
    
    # Direction (gain or loss of function)
    direction = 1 if (pos % 2 == 0) else -1
    
    # Generate signal based on variant type
    if shape == 'sharp':
        # Splice site: sharp peak with exponential decay
        signal = base_magnitude * direction * np.exp(-0.15 * (x)**2)
        # Add asymmetry (splice sites affect downstream more)
        signal[x > 0] *= 0.7
    elif shape == 'broad':
        # Regulatory: broad Gaussian with long tail
        signal = base_magnitude * direction * np.exp(-0.01 * (x)**2)
        # Add secondary peak (distal regulatory effect)
        secondary = 0.3 * base_magnitude * direction * np.exp(-0.02 * (x - 30)**2)
        signal += secondary
    else:
        # Coding: intermediate with slight asymmetry
        signal = base_magnitude * direction * np.exp(-0.04 * (x)**2)
        # Transitions (A<->G, C<->T) have slightly different patterns
        if is_transition:
            signal *= 0.85  # Transitions often have milder effects
    
    # Add realistic noise (heteroscedastic - more noise at extremes)
    noise_level = 0.2 + 0.1 * np.abs(x) / window_size
    noise = np.random.normal(0, noise_level, len(x))
    
    # Add occasional outliers (biological variability)
    outlier_mask = np.random.random(len(x)) < 0.05
    noise[outlier_mask] += np.random.normal(0, 0.8, np.sum(outlier_mask))
    
    delta_rna = signal + noise
    
    # Ensure realistic bounds (RNA-seq changes rarely exceed ±10)
    delta_rna = np.clip(delta_rna, -8, 8)
    
    # 3. Calculate Metrics
    max_idx = np.argmax(np.abs(delta_rna))
    max_delta = float(delta_rna[max_idx])
    max_pos_rel = int(x[max_idx])
    
    # 4. Real gnomAD Frequency (with fallback)
    from api_integrations import fetch_gnomad_frequency
    
    freq = fetch_gnomad_frequency(chrom, pos, ref, alt)
    if freq is None:
        # Fallback to mock if API fails
        freq = np.random.uniform(0.00001, 0.0001)
        print(f"Using mock frequency for {chrom}:{pos} (gnomAD API unavailable)")
    else:
        print(f"Real gnomAD frequency for {chrom}:{pos}: {freq}")
    
    # 5. Conservation Track - Real PhyloP from UCSC (with fallback)
    from api_integrations import fetch_ucsc_phylop
    
    # Calculate genomic coordinates for conservation window
    cons_start = max(0, pos - window_size)
    cons_end = pos + window_size + 1
    
    # Try to fetch real PhyloP scores
    cons_scores = fetch_ucsc_phylop(chrom, cons_start, cons_end)
    
    if cons_scores is not None and len(cons_scores) == len(x):
        # Successfully got real data
        print(f"✅ Using real PhyloP scores for {chrom}:{pos}")
        cons_scores = np.array(cons_scores)
    else:
        # Fallback to synthetic
        print(f"⚠️ Using synthetic conservation for {chrom}:{pos} (UCSC API unavailable)")
        cons_scores = np.random.normal(0.5, 1.0, len(x))
        cons_scores[window_size-10:window_size+10] += 2.0  # Conserved peak
    
    # Gene Structure (Exons) - Mock for now
    exons = [{"start": -50, "end": 20, "label": "Exon 1"}]
    
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
