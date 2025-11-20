import os
import requests
import pytest

BASE_URL = os.getenv("CARDIOVAR_BASE_URL", "http://localhost:8000")

def test_variant_impact_grch38():
    """Test variant impact with GRCh38 assembly (should succeed)"""
    payload = {
        "assembly": "GRCh38",
        "chrom": "chr22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C",
        "window_size": 100
    }
    print(f"Testing {BASE_URL}/variant-impact with GRCh38")
    try:
        r = requests.post(f"{BASE_URL}/variant-impact", json=payload, timeout=10)
        r.raise_for_status()
        data = r.json()
        
        # Check for expected top-level keys
        assert "variant_id" in data
        assert "metrics" in data
        assert "curve" in data
        assert "tracks" in data
        
        # Check specific metrics
        metrics = data["metrics"]
        for key in ["max_delta", "max_pos_rel", "gnomad_freq", "gene_symbol", "percentile"]:
            assert key in metrics
            
        print("✅ GRCh38 Test Passed!")
        
    except Exception as e:
        pytest.fail(f"GRCh38 test failed: {e}")

def test_variant_impact_grch37():
    """Test variant impact with GRCh37 assembly (should fail with 400)"""
    payload = {
        "assembly": "GRCh37",
        "chrom": "chr22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C",
        "window_size": 100
    }
    print(f"Testing {BASE_URL}/variant-impact with GRCh37 (expecting 400 error)")
    try:
        r = requests.post(f"{BASE_URL}/variant-impact", json=payload, timeout=10)
        
        # Should get 400 error
        assert r.status_code == 400, f"Expected 400, got {r.status_code}"
        
        # Check error message
        error_data = r.json()
        assert "detail" in error_data
        assert "GRCh38" in error_data["detail"], "Error message should mention GRCh38"
        
        print("✅ GRCh37 Test Passed (correctly rejected)!")
        
    except AssertionError:
        raise
    except Exception as e:
        pytest.fail(f"GRCh37 test failed: {e}")

def test_invalid_assembly():
    """Test with invalid assembly name (should fail with 400)"""
    payload = {
        "assembly": "hg19",  # Invalid format
        "chrom": "chr22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C"
    }
    print(f"Testing {BASE_URL}/variant-impact with invalid assembly")
    try:
        r = requests.post(f"{BASE_URL}/variant-impact", json=payload, timeout=10)
        
        # Should get 400 error
        assert r.status_code == 400, f"Expected 400, got {r.status_code}"
        
        print("✅ Invalid Assembly Test Passed!")
        
    except AssertionError:
        raise
    except Exception as e:
        pytest.fail(f"Invalid assembly test failed: {e}")

if __name__ == "__main__":
    print("🧪 Running Assembly Tests...\n")
    test_variant_impact_grch38()
    test_variant_impact_grch37()
    test_invalid_assembly()
    print("\n✨ All assembly tests completed!")
