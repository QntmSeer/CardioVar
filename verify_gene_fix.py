import requests
import time
import json

API_URL = "http://localhost:8000"

def test_gene_annotations():
    print("\n--- Testing Gene Annotations (MYH9) ---")
    # MYH9 is in our fallback data
    gene = "MYH9"
    try:
        resp = requests.get(f"{API_URL}/gene-annotations", params={"gene": gene})
        if resp.status_code == 200:
            data = resp.json()
            print(f"Success! Found data for {gene}")
            print(f"Source: {data.get('source')}")
            print(f"Protein Domains: {len(data.get('protein_domains', []))}")
            print(f"Expression Data: {len(data.get('expression', []))}")
            
            if "protein_domains" in data and len(data["protein_domains"]) > 0:
                print("Protein domains found (OK)")
            else:
                print("WARNING: No protein domains found")
                
        else:
            print(f"Failed: {resp.status_code} - {resp.text}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    time.sleep(3) # Wait for API
    test_gene_annotations()
