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
    Core logic for variant impact prediction using Enformer.
    
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
        RuntimeError: If Enformer is unavailable or fails to produce a prediction
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
    gene_symbol = "UNANNOTATED"
    if chrom in gene_map:
        for start, end, gene in gene_map[chrom]:
            if start <= pos <= end:
                gene_symbol = gene
                break
    
    # 2. Variant Impact Curve (Enformer required — no synthetic fallback)
    try:
        from enformer_wrapper import predict_variant_impact_dl
    except ImportError as e:
        raise RuntimeError(
            "Enformer is not installed or could not be imported. "
            "Install Enformer dependencies before running variant impact prediction."
        ) from e

    try:
        dl_result = predict_variant_impact_dl(chrom, pos, ref, alt)
    except Exception as e:
        raise RuntimeError(
            f"Enformer prediction failed for {chrom}:{pos}:{ref}>{alt}: {e}"
        ) from e

    if not dl_result or "raw_delta" not in dl_result or "center_idx" not in dl_result:
        raise RuntimeError(
            f"Enformer returned no valid prediction for {chrom}:{pos}:{ref}>{alt}."
        )

    x = np.arange(-window_size, window_size + 1)

    raw_profile = dl_result["raw_delta"]
    center = dl_result["center_idx"]
    # Extract ±5 Enformer bins (128 bp each) around the variant centre
    bin_subset = raw_profile[center-5:center+6]
    bin_x = np.arange(-5, 6) * 128
    from scipy.interpolate import interp1d
    delta_rna = interp1d(bin_x, bin_subset, kind='cubic', fill_value="extrapolate")(x)
    
    # 3. Calculate Metrics
    max_idx = np.argmax(np.abs(delta_rna))
    max_delta = float(delta_rna[max_idx])
    max_pos_rel = int(x[max_idx])
    
    # 4. Population Frequency (gnomAD → MyVariant.info → None/Unavailable)
    gnomad_source = "gnomAD API"
    freq = fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=force_live)
    if freq is None:
        from api_integrations import fetch_myvariant_info
        mv_data = fetch_myvariant_info(chrom, pos, ref, alt)
        if mv_data and "gnomad_genome" in mv_data:
            freq = mv_data["gnomad_genome"].get("af", {}).get("af")
            gnomad_source = "MyVariant.info (gnomAD genome)"
        elif mv_data and "gnomad_exome" in mv_data:
            freq = mv_data["gnomad_exome"].get("af", {}).get("af")
            gnomad_source = "MyVariant.info (gnomAD exome)"
        else:
            freq = None
            gnomad_source = "Unavailable (gnomAD & MyVariant.info lookup failed)"
            print(f">> gnomAD unavailable for {chrom}:{pos}; reporting Unavailable")
    
    # 5. PhyloP Conservation (UCSC API — no synthetic fallback)
    cons_start = max(0, pos - window_size)
    cons_end   = pos + window_size + 1
    cons_scores = fetch_ucsc_phylop(chrom, cons_start, cons_end, force_live=force_live)
    if cons_scores is not None and len(cons_scores) == len(x):
        cons_scores = np.array(cons_scores)
        used_real_cons = True
    else:
        cons_scores = None
        used_real_cons = False
        print(f">> PhyloP unavailable for {chrom}:{pos}; reporting Unavailable")

    # Exon structure (Ensembl API)
    exons = fetch_gene_structure(chrom, pos, window_size, force_live=force_live) or []
    
    # 6. Tissue Effects — GTEx v8 TPM scaled by Enformer |delta| (no synthetic fallback)
    tissue_effects = []
    tissue_effects_source = "GTEx v8 API"
    try:
        from api_integrations import fetch_gtex_expression
        gtex_data = fetch_gtex_expression(gene_symbol) if gene_symbol != "UNANNOTATED" else None
        if not gtex_data:
            raise ValueError("No GTEx expression data returned")
        max_tpm = max(t.get("tpm", 0) for t in gtex_data) or 1.0
        for t in gtex_data:
            tpm = t.get("tpm", 0)
            tissue_effects.append({"tissue": t["tissue"],
                                   "delta": round(abs(max_delta) * (tpm / max_tpm), 4),
                                   "tpm":   tpm})
    except Exception as e:
        print(f">> GTEx unavailable ({e}); reporting Unavailable")
        tissue_effects = []
        tissue_effects_source = "Unavailable (GTEx API failed or unannotated gene)"
    
    # 7. Background Distribution (pre-computed per gene — no synthetic fallback)
    background_deltas = load_background_distribution(gene_symbol) if gene_symbol != "UNANNOTATED" else None
    if background_deltas and len(background_deltas) > 1:
        all_deltas = background_deltas + [abs(max_delta)]
        percentile = round(float((np.sum(np.array(all_deltas) < abs(max_delta)) / len(all_deltas)) * 100), 1)
        bg_mean = np.mean(background_deltas)
        bg_std = np.std(background_deltas)
        z_score = round(float((abs(max_delta) - bg_mean) / bg_std), 2) if bg_std > 0 else 0.0
        bg_source = "Pre-computed gene background"
    else:
        background_deltas = None
        percentile = None
        z_score = None
        bg_source = "Unavailable (no pre-computed gene background)"
        print(f">> No pre-computed background for {gene_symbol}; reporting Unavailable")

    # 8. Real ClinVar neighbouring variants (±50 kb, GRCh38)
    try:
        from api_integrations import fetch_clinvar_region_variants
        clinvar_region = fetch_clinvar_region_variants(chrom, pos, window=50000, force_live=force_live)
    except Exception as e:
        print(f">> ClinVar region lookup failed ({e}); reporting Unavailable")
        clinvar_region = None
    clinvar_source = (
        f"NCBI ClinVar API ({len(clinvar_region)} variants in ±50 kb)"
        if clinvar_region is not None
        else "Unavailable (NCBI ClinVar API failed)"
    )

    # Track data sources for transparency
    enformer_top_tracks = dl_result.get("top_tracks")
    data_sources = {
        "variant_impact": "Enformer (Deep Learning)",
        "enformer_tracks": "Enformer (top 20 human tracks)" if enformer_top_tracks else "Unavailable",
        "gnomad_frequency": gnomad_source,
        "gene_structure": "Ensembl API" if exons else "Unavailable (no exons returned)",
        "conservation": "UCSC PhyloP API" if used_real_cons else "Unavailable (UCSC PhyloP API failed)",
        "gene_symbol": "Position-mapped (GRCh38)" if gene_symbol != "UNANNOTATED" else "Unannotated",
        "background_distribution": bg_source,
        "tissue_effects": tissue_effects_source,
        "clinvar": clinvar_source,
    }
    
    # 9. Fetch Gene Info
    gene_info = fetch_ensembl_gene(gene_symbol) if gene_symbol != "UNANNOTATED" else None
        
    # Confidence: penalise each unavailable data source
    confidence = 100.0
    if freq is None:                  confidence -= 20.0
    if not used_real_cons:            confidence -= 20.0
    if background_deltas is None:     confidence -= 10.0
    if not tissue_effects:            confidence -= 10.0
    if not gene_info:                 confidence -= 10.0
    if clinvar_region is None:        confidence -= 10.0
    confidence = max(0.0, confidence)
    
    return {
        "variant_id": f"{chrom}:{pos}:{ref}:{alt}",
        "metrics": {
            "max_delta": round(max_delta, 4),
            "max_pos_rel": max_pos_rel,
            "gnomad_freq": freq,
            "gene_symbol": gene_symbol,
            "percentile": percentile,
            "z_score": z_score,
            "confidence": round(confidence, 1),
            "fallback_used": api_integrations.fallback_used,
            "model_used": "Enformer (Deep Learning)"
        },
        "curve": {
            "x": x.tolist(),
            "y": delta_rna.tolist()
        },
        "tracks": {
            "exons": exons,
            "conservation": cons_scores.tolist() if cons_scores is not None else None
        },
        "enformer_top_tracks": enformer_top_tracks,
        "clinvar_region": clinvar_region,
        "gene": gene_info,
        "tissue_effects": tissue_effects,
        "background_distribution": {
            "background_deltas": background_deltas,
            "variant_delta": abs(max_delta)
        },
        "data_sources": data_sources
    }
