"""
Clinical Smoke Test Suite for CardioVar
Using real pathogenic variants that should be in ClinVar
"""
# Fix Windows console encoding
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests
import json
import time

API_URL = "http://localhost:8000"

# High-impact pathogenic variants from your clinical test suite
CLINICAL_VARIANTS = [
    {
        "name": "LMNA (DCM) - pathogenic missense",
        "chrom": "1",
        "pos": 156137666,
        "ref": "C",
        "alt": "T",
        "gene": "LMNA",
        "hgvs": "c.1621C>T p.Arg541Cys",
        "expected": "High percentile, ClinVar pathogenic, cardiac tissue relevance"
    },
    {
        "name": "SCN5A (Brugada/LQTS) - pathogenic missense",
        "chrom": "3",
        "pos": 38606709,
        "ref": "C",
        "alt": "T",
        "gene": "SCN5A",
        "hgvs": "c.1100G>A p.Arg367His",
        "expected": "High percentile, electrical/arrhythmia disease links"
    },
    {
        "name": "SCO2 (cardio-encephalomyopathy) - pathogenic",
        "chrom": "22",
        "pos": 50523994,
        "ref": "G",
        "alt": "A",
        "gene": "SCO2",
        "hgvs": "c.418G>A p.Glu140Lys",
        "expected": "High percentile, chr22 fallback, mito/cardiac disease"
    },
    {
        "name": "MYH9 - pathogenic missense",
        "chrom": "22",
        "pos": 36305975,
        "ref": "G",
        "alt": "A",
        "gene": "MYH9",
        "hgvs": "c.2114G>A p.Arg705His",
        "expected": "Gene annotation rendering, chr22 context plots"
    },
    {
        "name": "PCSK9 - protective/benign LoF (control)",
        "chrom": "1",
        "pos": 55039974,
        "ref": "G",
        "alt": "T",
        "gene": "PCSK9",
        "hgvs": "c.137G>T p.Arg46Leu",
        "expected": "Lower percentile, benign/protective ClinVar"
    }
]

def test_variant(variant):
    print(f"\n{'='*70}")
    print(f"🧬 {variant['name']}")
    print(f"{'='*70}")
    print(f"Gene: {variant['gene']}")
    print(f"Variant: chr{variant['chrom']}:{variant['pos']}:{variant['ref']}>{variant['alt']}")
    print(f"HGVS: {variant['hgvs']}")
    print(f"Expected: {variant['expected']}")
    print(f"{'-'*70}")
    
    try:
        # Test /related-data endpoint
        print(f"\n📡 Testing /related-data endpoint...")
        resp = requests.get(
            f"{API_URL}/related-data",
            params={
                "chrom": variant["chrom"],
                "pos": variant["pos"],
                "ref": variant["ref"],
                "alt": variant["alt"]
            },
            timeout=30
        )
        
        if resp.status_code == 200:
            data = resp.json()
            
            print(f"✅ Status: {resp.status_code} OK")
            
            # ClinVar Results
            print(f"\n🏥 ClinVar Results:")
            if data.get("clinvar"):
                clinvar = data["clinvar"]
                print(f"  ✅ ClinVar data FOUND!")
                print(f"  📋 Variant ID: {clinvar.get('variant_id')}")
                print(f"  🆔 Variation ID: {clinvar.get('variation_id')}")
                print(f"  ⚕️  Clinical Significance: {clinvar.get('clinical_significance')}")
                print(f"  ⭐ Review Status: {clinvar.get('review_status')}")
                print(f"  🌟 Gold Stars: {clinvar.get('gold_stars')}")
                print(f"  📅 Release Date: {clinvar.get('release_date')}")
                print(f"  📝 Submissions: {len(clinvar.get('submissions', []))}")
                
                if clinvar.get('submissions'):
                    print(f"\n  📄 First Submission Details:")
                    sub = clinvar['submissions'][0]
                    print(f"     Submitter: {sub.get('submitter')}")
                    print(f"     Significance: {sub.get('clinical_significance')}")
                    if sub.get('conditions'):
                        print(f"     Conditions:")
                        for cond in sub['conditions'][:3]:  # Show first 3
                            print(f"       - {cond.get('name')} (MedGen: {cond.get('medgen_id')})")
            else:
                print(f"  ⚠️  NO ClinVar data found")
                print(f"  ℹ️  This variant may not be in ClinVar database")
            
            # dbSNP Results
            print(f"\n🧬 dbSNP Results:")
            if data.get("dbsnp"):
                dbsnp = data["dbsnp"]
                print(f"  ✅ dbSNP data found!")
                print(f"  🆔 rsID: {dbsnp.get('rsid')}")
                print(f"  🧬 Gene: {dbsnp.get('gene')}")
                print(f"  📊 Global MAF: {dbsnp.get('global_maf')}")
            else:
                print(f"  ⚠️  No dbSNP data")
            
            # Summary
            has_clinvar = bool(data.get("clinvar"))
            has_dbsnp = bool(data.get("dbsnp"))
            
            print(f"\n{'='*70}")
            if has_clinvar:
                print(f"✅ SUCCESS: ClinVar data retrieved via gnomAD API")
            else:
                print(f"⚠️  WARNING: No ClinVar data (variant may not be in database)")
            print(f"{'='*70}")
                
        else:
            print(f"❌ API Error: {resp.status_code}")
            print(f"Response: {resp.text[:200]}")
            
    except Exception as e:
        print(f"❌ Exception: {type(e).__name__}: {e}")

def test_gene_annotations():
    """Test gene annotations endpoint"""
    print(f"\n{'='*70}")
    print(f"🧬 TESTING GENE ANNOTATIONS")
    print(f"{'='*70}")
    
    test_genes = ["LMNA", "MYH9", "SCN5A", "PCSK9"]
    
    for gene in test_genes:
        print(f"\n📋 Testing gene: {gene}")
        try:
            resp = requests.get(
                f"{API_URL}/gene-annotations",
                params={"gene": gene},
                timeout=15
            )
            if resp.status_code == 200:
                data = resp.json()
                print(f"  ✅ Success!")
                print(f"  Source: {data.get('source')}")
                print(f"  Protein Domains: {len(data.get('protein_domains', []))}")
                print(f"  Expression Data: {len(data.get('expression', []))}")
            else:
                print(f"  ❌ Error: {resp.status_code}")
        except Exception as e:
            print(f"  ❌ Exception: {e}")

if __name__ == "__main__":
    print("\n" + "="*70)
    print("🏥 CARDIOVAR CLINICAL SMOKE TEST SUITE")
    print("="*70)
    print("Testing ClinVar integration with known pathogenic variants")
    print("="*70)
    
    time.sleep(2)  # Wait for API
    
    # Test each clinical variant
    for i, variant in enumerate(CLINICAL_VARIANTS, 1):
        print(f"\n\n🧪 TEST {i}/{len(CLINICAL_VARIANTS)}")
        test_variant(variant)
        if i < len(CLINICAL_VARIANTS):
            time.sleep(3)  # Rate limiting between tests
    
    # Test gene annotations
    test_gene_annotations()
    
    print("\n\n" + "="*70)
    print("✅ TESTING COMPLETE")
    print("="*70)
    print("\nNext steps:")
    print("1. Review ClinVar results above")
    print("2. Test in dashboard UI with these variants")
    print("3. Check batch processing with batch_smoketest.csv")
    print("="*70)
