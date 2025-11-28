import json
import numpy as np
import pandas as pd
import os

# Define paths
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

GNOMAD_FILE = os.path.join(DATA_DIR, "gnomad_fallback.json")
PHYLOP_FILE = os.path.join(DATA_DIR, "phylop_fallback.npy")
GTEX_FILE = os.path.join(DATA_DIR, "gtex_expression.tsv")
DOMAINS_FILE = os.path.join(DATA_DIR, "protein_domains_fallback.json")
GENE_ANNOTATIONS_FILE = os.path.join(DATA_DIR, "gene_annotations.json")

# --- Genes and Regions to Mock ---
# Based on variant_engine.py gene_map
GENES = {
    "MYH9": {"chrom": "chr22", "start": 36100000, "end": 36400000, "protein_len": 1960},
    "LMNA": {"chrom": "chr1", "start": 156100000, "end": 156200000, "protein_len": 664},
    "PCSK9": {"chrom": "chr1", "start": 55000000, "end": 55100000, "protein_len": 692},
    "ACTN2": {"chrom": "chr1", "start": 236700000, "end": 236750000, "protein_len": 894},
    "APOB": {"chrom": "chr2", "start": 21000000, "end": 21100000, "protein_len": 4563},
    "TTN": {"chrom": "chr2", "start": 178500000, "end": 178600000, "protein_len": 34350},
    "MYL3": {"chrom": "chr3", "start": 46850000, "end": 46900000, "protein_len": 195},
}

def create_gnomad_fallback():
    print("Creating gnomAD fallback data...")
    # Create a few mock variants for each gene region
    data = {}
    for gene, info in GENES.items():
        chrom = info["chrom"].replace("chr", "")
        # Create 5 random variants per gene
        for _ in range(5):
            pos = np.random.randint(info["start"], info["end"])
            ref = np.random.choice(["A", "C", "G", "T"])
            alt = np.random.choice([b for b in ["A", "C", "G", "T"] if b != ref])
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
            data[variant_id] = float(np.random.uniform(0.00001, 0.01))
            
    # Add the specific test variant we use often
    data["22-36191400-A-C"] = 0.00012
    
    with open(GNOMAD_FILE, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved {len(data)} variants to {GNOMAD_FILE}")

def create_phylop_fallback():
    print("Creating PhyloP fallback data...")
    # We can't easily mock every base, but we can mock the specific test regions
    # For the demo, we'll just store a generic "random" array that the API integration
    # can use if it doesn't find an exact match, OR we can rely on the API integration's
    # logic to handle missing keys? 
    # Actually, the current api_integrations.py checks for exact key (chrom, start, end).
    # To make this robust for *any* query in the demo regions, we might need a different strategy.
    # But for now, let's just add the specific test window for MYH9.
    
    data = {}
    # MYH9 test variant window
    key = ("chr22", 36191300, 36191500) # 200bp window around 36191400
    # Generate realistic-looking conservation scores (positive = conserved)
    scores = np.random.normal(0.5, 1.0, 200)
    scores[90:110] += 2.0 # Peak conservation at variant
    data[key] = scores.tolist()
    
    np.save(PHYLOP_FILE, data)
    print(f"Saved PhyloP data to {PHYLOP_FILE}")

def create_gtex_fallback():
    print("Creating GTEx fallback data...")
    tissues = [
        "Heart - Left Ventricle", "Heart - Atrial Appendage", "Artery - Aorta", 
        "Artery - Coronary", "Liver", "Brain - Cortex", "Kidney - Cortex"
    ]
    
    rows = []
    for gene in GENES:
        # Heart genes get higher expression in heart
        is_heart_gene = gene in ["MYH9", "MYL3", "TTN", "ACTN2", "LMNA"]
        
        for tissue in tissues:
            if is_heart_gene and "Heart" in tissue:
                tpm = np.random.uniform(50, 200)
            elif is_heart_gene and "Artery" in tissue:
                tpm = np.random.uniform(20, 80)
            elif gene == "PCSK9" and "Liver" in tissue:
                tpm = np.random.uniform(100, 300)
            elif gene == "APOB" and ("Liver" in tissue or "Intestine" in tissue):
                tpm = np.random.uniform(200, 500)
            else:
                tpm = np.random.uniform(0.1, 10)
                
            rows.append({
                "gene_symbol": gene,
                "tissue": tissue,
                "tpm": round(tpm, 2)
            })
            
    df = pd.DataFrame(rows)
    df.to_csv(GTEX_FILE, sep="\t", index=False)
    print(f"Saved GTEx data for {len(GENES)} genes to {GTEX_FILE}")

def create_protein_domains_fallback():
    print("Creating Protein Domains fallback data...")
    data = {}
    
    for gene, info in GENES.items():
        # Mock domains
        domains = []
        length = info["protein_len"]
        num_domains = np.random.randint(2, 6)
        
        for i in range(num_domains):
            start = np.random.randint(1, length - 50)
            end = start + np.random.randint(20, 100)
            if end > length: end = length
            domains.append({
                "name": f"Domain_{i+1}",
                "start": start,
                "end": end,
                "type": "Pfam"
            })
            
        data[gene] = {
            "protein_length": length,
            "protein_domains": domains
        }
        
    with open(DOMAINS_FILE, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved domains for {len(GENES)} genes to {DOMAINS_FILE}")

def create_gene_annotations_fallback():
    print("Creating Gene Annotations fallback data...")
    data = []
    for gene, info in GENES.items():
        data.append({
            "symbol": gene,
            "name": f"{gene} mock description",
            "ensembl_id": f"ENSG00000{np.random.randint(100000, 999999)}",
            "chromosome": info["chrom"],
            "start": info["start"],
            "end": info["end"],
            "biotype": "protein_coding",
            "pathways": ["Pathway A", "Pathway B"],
            "disease_associations": ["Cardiomyopathy", "Arrhythmia"] if gene in ["MYH9", "TTN", "LMNA"] else ["Other"],
            "links": {"Genecards": "http://...", "OMIM": "http://..."}
        })
        
    with open(GENE_ANNOTATIONS_FILE, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved annotations for {len(GENES)} genes to {GENE_ANNOTATIONS_FILE}")

if __name__ == "__main__":
    create_gnomad_fallback()
    create_phylop_fallback()
    create_gtex_fallback()
    create_protein_domains_fallback()
    create_gene_annotations_fallback()
    print("All fallback data generated successfully.")
