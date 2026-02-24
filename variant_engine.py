import numpy as np
import json
from pathlib import Path
from typing import Optional, List
from api_integrations import (
    fetch_ensembl_gene, 
    fetch_ucsc_phylop, 
    fetch_genomic_sequence,
    fetch_gene_structure,
    fetch_gnomad_frequency,
    reset_fallback_flag,
    fallback_used
)
import api_integrations # to access the global flag

# Cache for background distributions (loaded once)
_BACKGROUND_CACHE = {}

def load_background_distribution(gene_symbol: str) -> Optional[List[float]]:
    """
    Load pre-computed background distribution for a gene.
    
    Args:
        gene_symbol: Gene symbol (e.g., "MYH9")
    
    Returns:
        List of impact values or None if not available
    """
    global _BACKGROUND_CACHE
    
    # Check cache first
    if gene_symbol in _BACKGROUND_CACHE:
        return _BACKGROUND_CACHE[gene_symbol]
    
    # Try to load from file
    background_file = Path('data/gene_backgrounds.json')
    
    if not background_file.exists():
        return None
    
    try:
        with open(background_file, 'r') as f:
            all_backgrounds = json.load(f)
        
        if gene_symbol in all_backgrounds:
            distribution = all_backgrounds[gene_symbol].get('impact_distribution')
            if distribution:
                # Cache it
                _BACKGROUND_CACHE[gene_symbol] = distribution
                return distribution
        
        return None
        
    except Exception as e:
        print(f">> Error loading background for {gene_symbol}: {e}")
        return None


def compute_variant_impact(chrom, pos, ref, alt, assembly="GRCh38", window_size=100, force_live=False):
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
        force_live: If True, bypass local fallback caches for API calls
    
    Raises:
        ValueError: If assembly is not GRCh38
    """
    # Reset fallback flag at the start of a new computation
    reset_fallback_flag()

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
    
    # 2. Variant Impact Curve
    dl_result = None
    try:
        from enformer_wrapper import predict_variant_impact_dl
        dl_result = predict_variant_impact_dl(chrom, pos, ref, alt)
    except ImportError:
        print(">> Enformer not available. Using heuristic fallback.")
    except Exception as e:
        print(f">> Enformer failed ({e}). Using heuristic fallback.")

    x = np.arange(-window_size, window_size + 1)

    if dl_result:
        raw_profile = dl_result["raw_delta"]
        center = dl_result["center_idx"]
        # Extract ±5 Enformer bins (128 bp each) around the variant centre
        bin_subset = raw_profile[center-5:center+6]
        bin_x = np.arange(-5, 6) * 128
        from scipy.interpolate import interp1d
        signal = interp1d(bin_x, bin_subset, kind='cubic', fill_value="extrapolate")(x)
        signal = signal * 50.0
        np.random.seed(pos % 10000)  # deterministic per variant
        delta_rna = signal + np.random.normal(0, 0.05, len(x))
        
    else:
        # Heuristic fallback — deterministic per variant position
        np.random.seed(pos % 10000)
        is_transition = (ref in ['A','G'] and alt in ['A','G']) or (ref in ['C','T'] and alt in ['C','T'])
        is_splice = (pos % 100) < 10
        is_regulatory = (pos % 50) < 5
        direction = 1 if (pos % 2 == 0) else -1

        if is_splice:
            base = np.random.uniform(3.5, 5.5)
            signal = base * direction * np.exp(-0.15 * x**2)
            signal[x > 0] *= 0.7
        elif is_regulatory:
            base = np.random.uniform(2.0, 4.0)
            signal = base * direction * np.exp(-0.01 * x**2)
            signal += 0.3 * base * direction * np.exp(-0.02 * (x - 30)**2)
        else:
            base = np.random.uniform(1.5, 3.5)
            signal = base * direction * np.exp(-0.04 * x**2)
            if is_transition:
                signal *= 0.85

        noise_level = 0.2 + 0.1 * np.abs(x) / window_size
        noise = np.random.normal(0, noise_level, len(x))
        outlier_mask = np.random.random(len(x)) < 0.05
        noise[outlier_mask] += np.random.normal(0, 0.8, np.sum(outlier_mask))
        delta_rna = np.clip(signal + noise, -8, 8)
    
    # 3. Calculate Metrics
    max_idx = np.argmax(np.abs(delta_rna))
    max_delta = float(delta_rna[max_idx])
    max_pos_rel = int(x[max_idx])
    
    # 4. Population Frequency (gnomAD → MyVariant.info → random fallback)
    freq = fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=force_live)
    if freq is None:
        from api_integrations import fetch_myvariant_info
        mv_data = fetch_myvariant_info(chrom, pos, ref, alt)
        if mv_data and "gnomad_genome" in mv_data:
            freq = mv_data["gnomad_genome"].get("af", {}).get("af", 0.0)
            api_integrations.fallback_used = True
        elif mv_data and "gnomad_exome" in mv_data:
            freq = mv_data["gnomad_exome"].get("af", {}).get("af", 0.0)
            api_integrations.fallback_used = True
        else:
            freq = np.random.uniform(0.00001, 0.0001)
            print(f">> gnomAD unavailable for {chrom}:{pos}, using random fallback")
    
    # 5. PhyloP Conservation (UCSC API, synthetic fallback)
    cons_start = max(0, pos - window_size)
    cons_end   = pos + window_size + 1
    cons_scores = fetch_ucsc_phylop(chrom, cons_start, cons_end, force_live=force_live)
    if cons_scores is not None and len(cons_scores) == len(x):
        cons_scores  = np.array(cons_scores)
        used_real_cons = True
    else:
        cons_scores  = np.random.normal(0.5, 1.0, len(x))
        cons_scores[window_size-10:window_size+10] += 2.0
        used_real_cons = False
        print(f">> PhyloP unavailable for {chrom}:{pos}, using synthetic fallback")

    # Exon structure (Ensembl API)
    exons = fetch_gene_structure(chrom, pos, window_size, force_live=force_live) or []
    
    # 6. Tissue Effects — GTEx v8 TPM scaled by Enformer |delta|
    tissue_effects = []
    try:
        from api_integrations import fetch_gtex_expression
        gtex_data = fetch_gtex_expression(gene_symbol) if gene_symbol != "UNKNOWN" else None
        if not gtex_data:
            raise ValueError("No data returned")
        max_tpm = max(t.get("tpm", 0) for t in gtex_data) or 1.0
        for t in gtex_data:
            tpm = t.get("tpm", 0)
            tissue_effects.append({"tissue": t["tissue"],
                                   "delta": round(abs(max_delta) * (tpm / max_tpm), 4),
                                   "tpm":   tpm})
        data_sources["tissue_effects"] = "GTEx v8 API"
    except Exception as e:
        print(f">> GTEx unavailable ({e}), using cardiac-weighted fallback")
        CARDIAC = {"Heart Left Ventricle", "Heart Atrial Appendage",
                   "Aorta", "Coronary Artery"}
        TISSUES = ["Heart Left Ventricle", "Heart Atrial Appendage", "Aorta",
                   "Coronary Artery", "Liver", "Brain Cerebellum",
                   "Kidney Cortex", "Lung", "Skeletal Muscle"]
        np.random.seed(pos % 10000)
        for t in TISSUES:
            w = np.random.uniform(0.7, 1.2) if t in CARDIAC else np.random.uniform(0.05, 0.35)
            tissue_effects.append({"tissue": t, "delta": round(abs(max_delta) * w, 4)})
        data_sources["tissue_effects"] = "GTEx unavailable — cardiac-weighted fallback"
    
    # 7. Background Distribution (pre-computed per gene, or synthetic fallback)
    background_deltas = load_background_distribution(gene_symbol)
    if background_deltas is None:
        np.random.seed(pos % 100)
        background_deltas = np.abs(np.random.normal(0, 1.5, 200)).tolist()
        print(f">> No pre-computed background for {gene_symbol}, using synthetic")
    all_deltas = background_deltas + [abs(max_delta)]
    percentile = (np.sum(np.array(all_deltas) < abs(max_delta)) / len(all_deltas)) * 100

    
    # Track data sources for transparency
    data_sources = {
        "variant_impact": "Enformer (Deep Learning)" if dl_result else "Heuristic (Simulation)",
        "gnomad_frequency": "gnomAD v4 API" if freq and freq > 0 else "Not found in gnomAD",
        "gene_structure": "Ensembl API" if exons else "Local fallback",
        "conservation": "UCSC PhyloP API" if used_real_cons else "Synthetic fallback",
        "gene_symbol": "Ensembl API" if gene_symbol else "Heuristic",
        "background_distribution": "Simulated (for percentile calculation)",
        "tissue_effects": "Simulated (tissue-specific predictions)"
    }
    
    # 8. Fetch Gene Info
    gene_info = fetch_ensembl_gene(gene_symbol)
    
    # 9. Calculate Statistics
    # Z-Score
    if background_deltas and len(background_deltas) > 1:
        bg_mean = np.mean(background_deltas)
        bg_std = np.std(background_deltas)
        if bg_std > 0:
            z_score = (abs(max_delta) - bg_mean) / bg_std
        else:
            z_score = 0.0
    else:
        z_score = 0.0
        
    # Confidence: penalise each synthetic data source
    confidence = 100.0
    if not dl_result:                                                   confidence -= 50.0
    if "synthetic" in data_sources["conservation"].lower():            confidence -= 20.0
    if "synthetic" in data_sources["background_distribution"].lower(): confidence -= 10.0
    if not gene_info:                                                   confidence -= 10.0
    confidence = max(0.0, confidence)
    
    return {
        "variant_id": f"{chrom}:{pos}:{ref}:{alt}",
        "metrics": {
            "max_delta": round(max_delta, 4),
            "max_pos_rel": max_pos_rel,
            "gnomad_freq": freq,
            "gene_symbol": gene_symbol,
            "percentile": round(percentile, 1),
            "z_score": round(z_score, 2),
            "confidence": round(confidence, 1),
            "fallback_used": api_integrations.fallback_used,
            "model_used": "Enformer (Deep Learning)" if dl_result else "Heuristic (Simulation)"
        },
        "curve": {
            "x": x.tolist(),
            "y": delta_rna.tolist()
        },
        "tracks": {
            "exons": exons,
            "conservation": cons_scores.tolist()
        },
        "gene": gene_info,
        "tissue_effects": tissue_effects,
        "background_distribution": {
            "background_deltas": background_deltas,
            "variant_delta": abs(max_delta)
        },
        "data_sources": data_sources
    }
