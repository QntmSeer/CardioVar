import requests
import time
import json

API_URL = "http://localhost:8000"

def test_clinvar_gnomad():
    print("\n--- Testing ClinVar via gnomAD (chr22-36191400-A-C) ---")
    # This is a known pathogenic variant in MYH9
    try:
        resp = requests.get(f"{API_URL}/related-data", params={
            "chrom": "22",
            "pos": 36191400,
            "ref": "A",
            "alt": "C"
        })
        if resp.status_code == 200:
            data = resp.json()
            print(f"Success! Found related data")
            print(f"ClinVar: {data.get('clinvar')}")
            print(f"dbSNP: {data.get('dbsnp')}")
            
            if data.get("clinvar"):
                clinvar = data["clinvar"]
                print(f"\nClinVar Details:")
                print(f"  Clinical Significance: {clinvar.get('clinical_significance')}")
                print(f"  Review Status: {clinvar.get('review_status')}")
                print(f"  Gold Stars: {clinvar.get('gold_stars')}")
                print(f"  Submissions: {len(clinvar.get('submissions', []))}")
            else:
                print("WARNING: No ClinVar data found")
                
        else:
            print(f"Failed: {resp.status_code} - {resp.text}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    time.sleep(3) # Wait for API
    test_clinvar_gnomad()
