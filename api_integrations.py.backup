"""
Real API integrations for CardioVar.
Fetches data from gnomAD, Ensembl, and other public databases.
"""
import requests
import json
from typing import Optional, Dict, Any, List

# API Endpoints
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_API = "https://rest.ensembl.org"

def fetch_gnomad_frequency(chrom: str, pos: int, ref: str, alt: str) -> Optional[float]:
    """
    Fetch real allele frequency from gnomAD v4.
    
    Args:
        chrom: Chromosome (e.g., "chr22" or "22")
        pos: Position
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        Allele frequency or None if not found
    """
    try:
        # Remove 'chr' prefix if present (gnomAD uses 1, 2, 3... not chr1, chr2...)
        chrom_clean = chrom.replace("chr", "")
        
        # GraphQL query for gnomAD
        query = """
        query VariantQuery($variantId: String!) {
          variant(variantId: $variantId, dataset: gnomad_r4) {
            genome {
              ac
              an
              af
            }
          }
        }
        """
        
        # Format: chrom-pos-ref-alt
        variant_id = f"{chrom_clean}-{pos}-{ref}-{alt}"
        
        variables = {"variantId": variant_id}
        
        response = requests.post(
            GNOMAD_API,
            json={"query": query, "variables": variables},
            timeout=10
        )
        
        if response.status_code == 200:
            data = response.json()
            if data.get("data") and data["data"].get("variant"):
                genome = data["data"]["variant"].get("genome")
                if genome and genome.get("af") is not None:
                    return float(genome["af"])
        
        return None
        
    except Exception as e:
        print(f"gnomAD API error: {e}")
        return None


def fetch_ensembl_gene(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """
    Fetch gene information from Ensembl REST API.
    
    Args:
        gene_symbol: Gene symbol (e.g., "MYH9")
    
    Returns:
        Gene data dictionary or None
    """
    try:
        # Lookup gene by symbol
        url = f"{ENSEMBL_API}/lookup/symbol/homo_sapiens/{gene_symbol}"
        headers = {"Content-Type": "application/json"}
        
        response = requests.get(url, headers=headers, timeout=10)
        
        if response.status_code == 200:
            gene_data = response.json()
            
            # Get additional details
            gene_id = gene_data.get("id")
            
            # Fetch protein domains (if coding gene)
            protein_data = None
            if gene_data.get("biotype") == "protein_coding":
                protein_url = f"{ENSEMBL_API}/overlap/id/{gene_id}?feature=protein_feature"
                protein_response = requests.get(protein_url, headers=headers, timeout=10)
                if protein_response.status_code == 200:
                    protein_data = protein_response.json()
            
            return {
                "symbol": gene_symbol,
                "name": gene_data.get("description", ""),
                "ensembl_id": gene_id,
                "chromosome": gene_data.get("seq_region_name"),
                "start": gene_data.get("start"),
                "end": gene_data.get("end"),
                "strand": gene_data.get("strand"),
                "biotype": gene_data.get("biotype"),
                "protein_features": protein_data
            }
        
        return None
        
    except Exception as e:
        print(f"Ensembl API error: {e}")
        return None


def fetch_clinvar_data(chrom: str, pos: int, ref: str, alt: str) -> list:
    """
    Fetch ClinVar data for a variant.
    Note: This is a simplified version. Full implementation would use NCBI E-utilities.
    
    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        List of ClinVar entries
    """
    try:
        # For now, return empty list
        # Full implementation would query:
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar
        
        # This requires more complex parsing of XML responses
        # Placeholder for future implementation
        return []
        
    except Exception as e:
        print(f"ClinVar API error: {e}")
        return []


def fetch_ucsc_phylop(chrom: str, start: int, end: int) -> Optional[list]:
    """
    Fetch real PhyloP conservation scores from UCSC Genome Browser.
    
    PhyloP scores measure evolutionary conservation:
    - Positive scores: conserved (slower evolution)
    - Negative scores: fast-evolving
    - Range typically -14 to +6
    
    Args:
        chrom: Chromosome (e.g., "chr22")
        start: Start position (0-based)
        end: End position (exclusive)
    
    Returns:
        List of PhyloP scores or None if fetch fails
    """
    try:
        # UCSC API endpoint for PhyloP100way (100 vertebrate species)
        base_url = "https://api.genome.ucsc.edu/getData/track"
        
        params = {
            "genome": "hg38",
            "track": "phyloP100way",
            "chrom": chrom,
            "start": start,
            "end": end
        }
        
        response = requests.get(base_url, params=params, timeout=15)
        
        if response.status_code == 200:
            data = response.json()
            
            # Extract scores from response
            if "phyloP100way" in data:
                track_data = data["phyloP100way"]
                
                # UCSC returns data in various formats, handle the most common
                if isinstance(track_data, dict) and "data" in track_data:
                    scores = track_data["data"]
                elif isinstance(track_data, list):
                    scores = track_data
                else:
                    return None
                
                # Ensure we have the right number of scores
                expected_length = end - start
                if len(scores) == expected_length:
                    return [float(s) if s is not None else 0.0 for s in scores]
        
        return None
        
    except Exception as e:
        print(f"UCSC PhyloP API error: {e}")
        return None


def fetch_gtex_expression(gene_symbol: str) -> Optional[List[Dict[str, Any]]]:
    """
    Fetch real tissue expression from GTEx API.
    Returns median TPM values for key tissues.
    
    Args:
        gene_symbol: Gene symbol (e.g., "MYH9")
        
    Returns:
        List of dictionaries with 'tissue' and 'tpm' keys
    """
    try:
        # GTEx API v2
        url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        
        params = {
            "geneId": gene_symbol,
            "tissueSiteDetailId": [
                "Heart_Left_Ventricle",
                "Heart_Atrial_Appendage", 
                "Artery_Aorta",
                "Artery_Coronary",
                "Liver",
                "Brain_Cortex",
                "Kidney_Cortex",
                "Muscle_Skeletal"
            ]
        }
        
        headers = {"Content-Type": "application/json"}
        
        response = requests.get(url, params=params, headers=headers, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            if "medianGeneExpression" in data:
                results = []
                # Map GTEx tissue names to our display names
                name_map = {
                    "Heart_Left_Ventricle": "Heart LV",
                    "Heart_Atrial_Appendage": "Heart RA",
                    "Artery_Aorta": "Aorta",
                    "Artery_Coronary": "Coronary Artery",
                    "Liver": "Liver",
                    "Brain_Cortex": "Brain",
                    "Kidney_Cortex": "Kidney",
                    "Muscle_Skeletal": "Skeletal Muscle"
                }
                
                for item in data["medianGeneExpression"]:
                    tissue_id = item.get("tissueSiteDetailId")
                    tpm = item.get("median")
                    
                    if tissue_id in name_map and tpm is not None:
                        results.append({
                            "tissue": name_map[tissue_id],
                            "tpm": float(tpm)
                        })
                
                if results:
                    return results
        
        return None
        
    except Exception as e:
        print(f"GTEx API error: {e}")
        return None


def fetch_genomic_sequence(chrom: str, start: int, end: int) -> Optional[str]:
    """
    Fetch genomic sequence from UCSC API.
    Needed for deep learning models (Enformer).
    
    Args:
        chrom: Chromosome (e.g., "chr1")
        start: Start position (0-based)
        end: End position (exclusive)
        
    Returns:
        DNA sequence string (upper case) or None
    """
    try:
        # UCSC API limit is often 10MB, so 200kb is fine
        base_url = "https://api.genome.ucsc.edu/getData/sequence"
        
        params = {
            "genome": "hg38",
            "chrom": chrom,
            "start": start,
            "end": end
        }
        
        response = requests.get(base_url, params=params, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            if "dna" in data:
                return data["dna"].upper()
                
        return None
        
    except Exception as e:
        print(f"UCSC Sequence API error: {e}")
        return None


def load_fallback_gene_data(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """
    Load gene data from local JSON file as fallback.
    
    Args:
        gene_symbol: Gene symbol
    
    Returns:
        Gene data or None
    """
    try:
        with open("data/gene_annotations.json", "r") as f:
            genes = json.load(f)
            for gene in genes:
                if gene["symbol"].upper() == gene_symbol.upper():
                    return gene
        return None
    except Exception as e:
        print(f"Fallback data error: {e}")
        return None


def load_fallback_related_data(chrom: str, pos: int, ref: str, alt: str) -> list:
    """
    Load related variant data from local JSON file as fallback.
    
    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        List of related data entries
    """
    try:
        with open("data/related_variants.json", "r") as f:
            data = json.load(f)
            key = f"{chrom}:{pos}:{ref}:{alt}"
            return data.get(key, [])
    except Exception as e:
        print(f"Fallback related data error: {e}")
        return []
