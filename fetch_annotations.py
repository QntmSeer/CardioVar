import random
import pandas as pd

def fetch_clinvar(chrom, pos, ref, alt):
    """
    Mock function to fetch ClinVar significance.
    """
    # Mock logic: Randomly assign significance based on position
    seed = pos % 10
    if seed < 2:
        return "Pathogenic"
    elif seed < 4:
        return "Likely Pathogenic"
    elif seed < 7:
        return "Benign"
    else:
        return "Uncertain Significance"

def fetch_gwas(chrom, pos):
    """
    Mock function to fetch nearby GWAS traits.
    """
    traits = [
        "Atrial Fibrillation",
        "Coronary Artery Disease",
        "Hypertension",
        "Heart Failure",
        "QT Interval",
        "None"
    ]
    # Deterministic mock
    idx = (pos // 1000) % len(traits)
    return traits[idx]

def fetch_heart_expression(gene):
    """
    Mock function to return expression levels in heart tissues.
    """
    # Mock expression data (TPM)
    return {
        "Heart - Left Ventricle": round(random.uniform(10, 100), 2),
        "Heart - Atrial Appendage": round(random.uniform(5, 80), 2),
        "Artery - Coronary": round(random.uniform(2, 50), 2)
    }

def get_heart_assays():
    """
    Returns a list of available heart-specific assays.
    """
    return pd.DataFrame([
        {"Assay ID": "A001", "Tissue": "Heart Left Ventricle", "Ontology": "UBERON:0002084", "Source": "ENCODE"},
        {"Assay ID": "A002", "Tissue": "Heart Atrial Appendage", "Ontology": "UBERON:0006618", "Source": "GTEx"},
        {"Assay ID": "A003", "Tissue": "Coronary Artery", "Ontology": "UBERON:0002111", "Source": "Roadmap Epigenomics"},
        {"Assay ID": "A004", "Tissue": "Cardiomyocyte", "Ontology": "CL:0000746", "Source": "AlphaGenome Internal"}
    ])

def fetch_gnomad_freq(chrom, pos, ref, alt):
    """
    Mock function to fetch allele frequency from gnomAD.
    """
    # Mock logic: Rare variants for 'A', common for others
    if ref == 'A':
        return 0.00005 # Very rare
    else:
        return 0.05 # Common (5%)
