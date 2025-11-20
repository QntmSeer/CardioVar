import time
import sys
from utils import get_delta_rna_seq, get_gene_info
from fetch_annotations import fetch_clinvar

def run_test(name, check_func):
    sys.stdout.write(f"TEST: {name} ... ")
    sys.stdout.flush()
    time.sleep(0.5) # Simulate work
    try:
        check_func()
        print("PASS")
        return True
    except Exception as e:
        print(f"FAIL ({e})")
        return False

def test_pipeline():
    print("Starting CardioVar Live Autotest...\n")
    
    # Test 1: Data Generation
    def check_data():
        x, y = get_delta_rna_seq("chr22", 36191400, "A", "C")
        assert len(x) == len(y)
        assert len(x) > 0
    run_test("Synthetic Data Generation", check_data)
    
    # Test 2: Gene Lookup
    def check_genes():
        hits = get_gene_info("chr22", 36191400)
        assert any(h['Gene'] == 'MYH9' for h in hits)
    run_test("Gene Annotation Lookup", check_genes)
    
    # Test 3: ClinVar Fetch
    def check_clinvar():
        sig = fetch_clinvar("chr22", 36191400, "A", "C")
        assert isinstance(sig, str)
    run_test("External API (ClinVar) Mock", check_clinvar)
    
    # Test 4: Edge Case (Negative Position)
    def check_edge():
        x, y = get_delta_rna_seq("chr1", -100, "A", "T")
        # Should run without crashing, even if pos is weird
        assert len(x) > 0
    run_test("Edge Case Handling", check_edge)
    
    print("\nAll tests completed.")

if __name__ == "__main__":
    test_pipeline()
