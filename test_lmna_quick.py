# Fix Windows console encoding
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests

API_URL = "http://localhost:8000"

# Test just the LMNA variant
print("Testing LMNA pathogenic variant...")
print("chr1:156137666:C>T (c.1621C>T p.Arg541Cys)")
print("="*60)

resp = requests.get(
    f"{API_URL}/related-data",
    params={
        "chrom": "1",
        "pos": 156137666,
        "ref": "C",
        "alt": "T"
    },
    timeout=30
)

print(f"Status: {resp.status_code}")

if resp.status_code == 200:
    data = resp.json()
    
    print("\n✅ API Response received")
    
    if data.get("clinvar"):
        print("\n🎉 ClinVar DATA FOUND!")
        print(f"Data: {data['clinvar']}")
    else:
        print("\n⚠️  No ClinVar data")
    
    if data.get("dbsnp"):
        print(f"\ndbSNP: {data['dbsnp']}")
else:
    print(f"Error: {resp.text}")
