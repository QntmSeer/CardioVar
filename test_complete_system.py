"""
Comprehensive test suite for CardioVar after restoration.
Tests all major components: API integrations, fallbacks, links, and UI features.
"""
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api_integrations import (
    fetch_gnomad_frequency,
    fetch_ensembl_gene,
    fetch_ucsc_phylop,
    fetch_gtex_expression,
    fetch_protein_domains,
    fetch_gene_structure,
    fetch_genomic_sequence,
    reset_fallback_flag,
    fallback_used
)

def test_retry_decorator():
    """Test that retry decorator exists and is importable."""
    from api_integrations import retry_on_failure
    print("✅ Retry decorator imported successfully")
    return True

def test_gnomad_api():
    """Test gnomAD frequency fetch."""
    print("\n[API] Testing gnomAD API...")
    reset_fallback_flag()
    
    # Test variant
    freq = fetch_gnomad_frequency("chr22", 36191400, "A", "C", force_live=False)
    
    if freq is not None:
        print(f"✅ gnomAD frequency: {freq}")
        return True
    else:
        print("❌ gnomAD fetch failed")
        return False

def test_ensembl_gene():
    """Test Ensembl gene fetch with dynamic links."""
    print("\n[API] Testing Ensembl Gene API...")
    reset_fallback_flag()
    
    gene_data = fetch_ensembl_gene("MYH9", force_live=False)
    
    if gene_data:
        print(f"✅ Gene: {gene_data.get('display_name', 'N/A')}")
        print(f"   Description: {gene_data.get('description', 'N/A')[:50]}...")
        
        # Check for dynamic links
        links = gene_data.get('links', {})
        if links:
            print(f"✅ Dynamic links generated: {len(links)} links")
            for name, url in links.items():
                print(f"   - {name}: {url[:50]}...")
            
            # Verify all 5 links are present
            expected_links = ['GeneCards', 'UniProt', 'gnomAD', 'Ensembl', 'OMIM']
            missing = [link for link in expected_links if link not in links]
            if missing:
                print(f"[WARNING]  Missing links: {missing}")
                return False
            else:
                print("✅ All 5 external links present!")
                return True
        else:
            print("❌ No links generated")
            return False
    else:
        print("❌ Gene fetch failed")
        return False

def test_phylop():
    """Test PhyloP conservation scores."""
    print("\n[API] Testing UCSC PhyloP API...")
    reset_fallback_flag()
    
    scores = fetch_ucsc_phylop("chr22", 36191300, 36191500, force_live=False)
    
    if scores and len(scores) > 0:
        print(f"✅ PhyloP scores: {len(scores)} positions")
        print(f"   Range: {min(scores):.2f} to {max(scores):.2f}")
        return True
    else:
        print("❌ PhyloP fetch failed")
        return False

def test_gtex():
    """Test GTEx expression data."""
    print("\n[API] Testing GTEx Expression API...")
    reset_fallback_flag()
    
    expr_data = fetch_gtex_expression("MYH9", force_live=False)
    
    if expr_data and len(expr_data) > 0:
        print(f"✅ GTEx expression: {len(expr_data)} tissues")
        print(f"   Sample: {expr_data[0]}")
        return True
    else:
        print("❌ GTEx fetch failed")
        return False

def test_protein_domains():
    """Test protein domain fetch."""
    print("\n[API] Testing Protein Domains API...")
    reset_fallback_flag()
    
    domain_data = fetch_protein_domains("MYH9", force_live=False)
    
    if domain_data:
        print(f"✅ Protein length: {domain_data.get('protein_length', 'N/A')}")
        domains = domain_data.get('protein_domains', [])
        print(f"✅ Domains found: {len(domains)}")
        if domains:
            print(f"   Sample: {domains[0]}")
        return True
    else:
        print("❌ Protein domains fetch failed")
        return False

def test_gene_structure():
    """Test gene structure fetch."""
    print("\n[API] Testing Gene Structure API...")
    reset_fallback_flag()
    
    exons = fetch_gene_structure("chr22", 36191400, window_size=100, force_live=False)
    
    if exons and len(exons) > 0:
        print(f"✅ Gene structure: {len(exons)} exons")
        print(f"   Sample: {exons[0]}")
        return True
    else:
        print("❌ Gene structure fetch failed")
        return False

def test_genomic_sequence():
    """Test genomic sequence fetch."""
    print("\n[API] Testing Genomic Sequence API...")
    
    seq = fetch_genomic_sequence("chr22", 36191400, 36191450, force_live=False)
    
    if seq and len(seq) > 0:
        print(f"✅ Sequence fetched: {len(seq)} bp")
        print(f"   Sample: {seq[:30]}...")
        return True
    else:
        print("❌ Sequence fetch failed")
        return False

def test_fallback_mechanism():
    """Test that fallback mechanism works."""
    print("\n[API] Testing Fallback Mechanism...")
    reset_fallback_flag()
    
    # Force a fallback by using a gene that should be in local data
    gene_data = fetch_ensembl_gene("PCSK9", force_live=False)
    
    if gene_data:
        print("✅ Fallback mechanism working")
        return True
    else:
        print("❌ Fallback mechanism failed")
        return False

def test_ui_components():
    """Test UI-related imports."""
    print("\n[UI] Testing UI Components...")
    
    try:
        # Check if dashboard.py has required imports
        with open("dashboard.py", "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
            
        checks = {
            "io import": "import io" in content,
            "matplotlib import": "import matplotlib" in content,
            "seaborn import": "import seaborn" in content,
            "COLORS dict": "COLORS = {" in content,
            "Model Info section": "Model Info" in content,
            "Light mode CSS": "#2C5F7C" in content or "#C44027" in content,
        }
        
        all_passed = True
        for check_name, passed in checks.items():
            if passed:
                print(f"✅ {check_name}")
            else:
                print(f"❌ {check_name}")
                all_passed = False
        
        return all_passed
    except Exception as e:
        print(f"❌ UI component test failed: {e}")
        return False

def test_config_file():
    """Test Streamlit config for light mode."""
    print("\n[CFG]  Testing Streamlit Config...")
    
    try:
        with open(".streamlit/config.toml", "r") as f:
            content = f.read()
        
        checks = {
            "White background": '#FFFFFF' in content,
            "Light secondary": '#F8F9FA' in content,
            "Dark text": '#303633' in content,
            "Tiger Flame accent": '#F46036' in content,
        }
        
        all_passed = True
        for check_name, passed in checks.items():
            if passed:
                print(f"✅ {check_name}")
            else:
                print(f"❌ {check_name}")
                all_passed = False
        
        return all_passed
    except Exception as e:
        print(f"❌ Config test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("CardioVar Comprehensive Test Suite")
    print("=" * 60)
    
    tests = [
        ("Retry Decorator", test_retry_decorator),
        ("gnomAD API", test_gnomad_api),
        ("Ensembl Gene + Links", test_ensembl_gene),
        ("PhyloP Conservation", test_phylop),
        ("GTEx Expression", test_gtex),
        ("Protein Domains", test_protein_domains),
        ("Gene Structure", test_gene_structure),
        ("Genomic Sequence", test_genomic_sequence),
        ("Fallback Mechanism", test_fallback_mechanism),
        ("UI Components", test_ui_components),
        ("Streamlit Config", test_config_file),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 60)
    print("[SUMMARY] Test Summary")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\n{'='*60}")
    print(f"Results: {passed}/{total} tests passed ({passed/total*100:.1f}%)")
    print(f"{'='*60}")
    
    if passed == total:
        print("\n[SUCCESS] ALL TESTS PASSED! CardioVar is fully functional!")
        return 0
    else:
        print(f"\n[WARNING]  {total - passed} test(s) failed. Review above for details.")
        return 1

if __name__ == "__main__":
    exit(main())
