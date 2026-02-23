import io
import pandas as pd

def parse_vcf(file_content: str) -> pd.DataFrame:
    """
    Parse a VCF file string into a pandas DataFrame with chrom, pos, ref, alt.
    Handles standard VCF format (tab-separated, header lines start with #).
    """
    lines = file_content.splitlines()
    data = []
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
            
        parts = line.split("\t")
        if len(parts) < 5:
            continue
            
        chrom = parts[0]
        # Ensure chrom starts with 'chr' if not present (optional, but good for consistency)
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
            
        try:
            pos = int(parts[1])
        except ValueError:
            continue
            
        ref = parts[3]
        alts = parts[4].split(",")  # Handle multiple ALTs
        
        for alt in alts:
            data.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt
            })
            
    return pd.DataFrame(data)
