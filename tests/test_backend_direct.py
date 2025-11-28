
from fastapi.testclient import TestClient
import sys
import os

# Add parent directory to path to import api
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api import app

client = TestClient(app)

def test_variant_impact_endpoint():
    """Test the /variant-impact endpoint directly using TestClient."""
    payload = {
        "chrom": "chr22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C",
        "window_size": 100,
        "force_live": False
    }
    
    print(f"Testing /variant-impact with payload: {payload}")
    response = client.post("/variant-impact", json=payload)
    
    assert response.status_code == 200
    data = response.json()
    
    # Check structure
    assert "variant_id" in data
    assert "metrics" in data
    assert "curve" in data
    assert "tracks" in data
    
    # Check metrics
    metrics = data["metrics"]
    assert "max_delta" in metrics
    assert "gene_symbol" in metrics
    
    print("✅ /variant-impact endpoint passed")

def test_gene_annotations_endpoint():
    """Test the /gene-annotations endpoint directly."""
    gene = "MYH9"
    print(f"Testing /gene-annotations with gene: {gene}")
    
    response = client.get(f"/gene-annotations?gene={gene}")
    
    assert response.status_code == 200
    data = response.json()
    
    assert data["symbol"] == gene
    assert "name" in data
    assert "expression" in data or "source" in data
    
    print("✅ /gene-annotations endpoint passed")

def test_system_status_endpoint():
    """Test the /system-status endpoint."""
    response = client.get("/system-status")
    assert response.status_code == 200
    data = response.json()
    assert "cpu_percent" in data
    assert "memory_percent" in data
    print("✅ /system-status endpoint passed")

if __name__ == "__main__":
    try:
        test_variant_impact_endpoint()
        test_gene_annotations_endpoint()
        test_system_status_endpoint()
        print("\nALL DIRECT BACKEND TESTS PASSED!")
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        exit(1)
