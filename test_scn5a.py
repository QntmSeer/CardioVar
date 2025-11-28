# Fix Windows console encoding
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests

API_URL = "http://localhost:8000"

print("Testing SCN5A variant (chr3:38606709:C>T)")
print("="*60)

# Test gene annotations
print("\n1. Testing /gene-annotations for SCN5A...")
try:
    r = requests.get(f"{API_URL}/gene-annotations", params={"gene": "SCN5A"}, timeout=10)
    print(f"Status: {r.status_code}")
    if r.status_code == 200:
        data = r.json()
        print(f"Source: {data.get('source')}")
        print(f"Protein Domains: {len(data.get('protein_domains', []))}")
        print(f"Expression Data: {len(data.get('expression', []))}")
    else:
        print(f"Error: {r.text}")
except Exception as e:
    print(f"Exception: {e}")

# Test related data
print("\n2. Testing /related-data for chr3:38606709:C>T...")
try:
    r = requests.get(
        f"{API_URL}/related-data",
        params={"chrom": "3", "pos": 38606709, "ref": "C", "alt": "T"},
        timeout=10
    )
    print(f"Status: {r.status_code}")
    if r.status_code == 200:
        data = r.json()
        print(f"ClinVar: {'✅ Found' if data.get('clinvar') else '❌ Not found'}")
        if data.get('clinvar'):
            print(f"  Classification: {data['clinvar'].get('germline_classification')}")
            print(f"  Gene: {data['clinvar'].get('gene_symbol')}")
        print(f"dbSNP: {'✅ Found' if data.get('dbsnp') else '❌ Not found'}")
        if data.get('dbsnp'):
            print(f"  rsID: {data['dbsnp'].get('rsid')}")
    else:
        print(f"Error: {r.text}")
except Exception as e:
    print(f"Exception: {e}")

# Test variant impact
print("\n3. Testing /variant-impact for chr3:38606709:C>T...")
try:
    r = requests.post(
        f"{API_URL}/variant-impact",
        json={"chrom": "chr3", "pos": 38606709, "ref": "C", "alt": "T"},
        timeout=10
    )
    print(f"Status: {r.status_code}")
    if r.status_code == 200:
        data = r.json()
        print(f"Variant ID: {data.get('variant_id')}")
        print(f"Gene: {data.get('metrics', {}).get('gene_symbol')}")
        print(f"Max Delta: {data.get('metrics', {}).get('max_delta')}")
    else:
        print(f"Error: {r.text[:200]}")
except Exception as e:
    print(f"Exception: {e}")

print("\n" + "="*60)
print("Testing complete!")
