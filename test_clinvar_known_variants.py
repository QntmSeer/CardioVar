"""
Test ClinVar integration with known pathogenic variants from gnomAD
"""
import requests
import json
import time

API_URL = "http://localhost:8000"

# Known pathogenic variants that should be in ClinVar
TEST_VARIANTS = [
    {
        "name": "BRCA1 pathogenic",
        "chrom": "17",
        "pos": 43094464,
        "ref": "A",
        "alt": "T",
        "expected": "Should have ClinVar data"
    },
    {
        "name": "APOE e4 allele",
        "chrom": "19",
        "pos": 44908684,
        "ref": "T",
        "alt": "C",
        "expected": "Common variant, likely in ClinVar"
    },
    {
        "name": "LDLR pathogenic",
        "chrom": "19",
        "pos": 11091630,
        "ref": "G",
        "alt": "A",
        "expected": "Should have ClinVar data"
    },
    {
        "name": "MYH9 test variant",
        "chrom": "22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C",
        "expected": "May or may not be in ClinVar"
    }
]

def test_variant(variant):
    print(f"\n{'='*60}")
    print(f"Testing: {variant['name']}")
    print(f"Variant: chr{variant['chrom']}:{variant['pos']}:{variant['ref']}>{variant['alt']}")
    print(f"Expected: {variant['expected']}")
    print(f"{'='*60}")
    
    try:
        # Test the /related-data endpoint
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
            
            print(f"\n✅ API Response: {resp.status_code}")
            print(f"\nClinVar Data:")
            if data.get("clinvar"):
                clinvar = data["clinvar"]
                print(f"  ✅ ClinVar data found!")
                print(f"  Variant ID: {clinvar.get('variant_id')}")
                print(f"  Variation ID: {clinvar.get('variation_id')}")
                print(f"  Clinical Significance: {clinvar.get('clinical_significance')}")
                print(f"  Review Status: {clinvar.get('review_status')}")
                print(f"  Gold Stars: {clinvar.get('gold_stars')}")
                print(f"  Release Date: {clinvar.get('release_date')}")
                print(f"  Submissions: {len(clinvar.get('submissions', []))}")
                
                if clinvar.get('submissions'):
                    print(f"\n  First Submission:")
                    sub = clinvar['submissions'][0]
                    print(f"    Submitter: {sub.get('submitter')}")
                    print(f"    Significance: {sub.get('clinical_significance')}")
                    print(f"    Conditions: {len(sub.get('conditions', []))}")
            else:
                print(f"  ⚠️  No ClinVar data (variant may not be in ClinVar)")
            
            print(f"\ndbSNP Data:")
            if data.get("dbsnp"):
                dbsnp = data["dbsnp"]
                print(f"  ✅ dbSNP data found!")
                print(f"  rsID: {dbsnp.get('rsid')}")
                print(f"  Gene: {dbsnp.get('gene')}")
                print(f"  Global MAF: {dbsnp.get('global_maf')}")
            else:
                print(f"  ⚠️  No dbSNP data")
                
        else:
            print(f"\n❌ API Error: {resp.status_code}")
            print(f"Response: {resp.text}")
            
    except Exception as e:
        print(f"\n❌ Exception: {e}")

if __name__ == "__main__":
    print("\n" + "="*60)
    print("TESTING CLINVAR INTEGRATION VIA GNOMAD")
    print("="*60)
    
    time.sleep(2)  # Wait for API to be ready
    
    for variant in TEST_VARIANTS:
        test_variant(variant)
        time.sleep(2)  # Rate limiting
    
    print("\n" + "="*60)
    print("TESTING COMPLETE")
    print("="*60)
