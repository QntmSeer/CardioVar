import os
import requests
import pytest

BASE_URL = os.getenv("CARDIOVAR_BASE_URL", "http://localhost:8000")

def test_variant_impact_smoke():
    payload = {
        "chrom": "chr22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C",
        "window_size": 100
    }
    print(f"Testing {BASE_URL}/variant-impact with payload: {payload}")
    try:
        # Increased timeout to 60s because Enformer loading/inference on CPU is slow
        response = requests.post(f"{BASE_URL}/variant-impact", json=payload, timeout=60)
        response.raise_for_status()
        data = response.json()
        
        # Check for expected top-level keys based on variant_engine.py
        # Expected: variant_id, metrics, curve, tracks
        assert "variant_id" in data
        assert "metrics" in data
        assert "curve" in data
        assert "tracks" in data
        
        # Check specific metrics
        metrics = data["metrics"]
        for key in ["max_delta", "max_pos_rel", "gnomad_freq", "gene_symbol"]:
            assert key in metrics
            
        print("✅ Backend Smoke Test Passed!")
        
    except Exception as e:
        pytest.fail(f"Backend smoke test failed: {e}")

if __name__ == "__main__":
    test_variant_impact_smoke()
