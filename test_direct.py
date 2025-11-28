"""Direct test of variant_engine to see actual error"""
import sys
sys.path.insert(0, '.')

try:
    from variant_engine import compute_variant_impact
    
    print("Testing variant_engine.compute_variant_impact...")
    result = compute_variant_impact(
        chrom="chr22",
        pos=36191400,
        ref="A",
        alt="C",
        window_size=100
    )
    
    print("[SUCCESS] Variant engine test passed!")
    print(f"Variant ID: {result['variant_id']}")
    print(f"Gene: {result['metrics']['gene_symbol']}")
    print(f"Max Delta: {result['metrics']['max_delta']}")
    print(f"Exons found: {len(result['tracks']['exons'])}")
    if result['tracks']['exons']:
        print(f"First exon: {result['tracks']['exons'][0]}")
    
except Exception as e:
    print(f"[ERROR] {e}")
    import traceback
    traceback.print_exc()
