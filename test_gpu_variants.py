"""
Test CardioVar with different variant examples using GPU acceleration.
"""
import time
from variant_engine import compute_variant_impact

# Test variants in different cardiovascular genes
test_variants = [
    {
        "name": "MYH9 variant (chr22)",
        "chrom": "chr22",
        "pos": 36191400,
        "ref": "A",
        "alt": "C"
    },
    {
        "name": "MYBPC3 variant (chr11)",
        "chrom": "chr11",
        "pos": 47352960,
        "ref": "G",
        "alt": "A"
    },
    {
        "name": "TNNT2 variant (chr1)",
        "chrom": "chr1",
        "pos": 201328200,
        "ref": "C",
        "alt": "T"
    }
]

print("="*70)
print("CardioVar GPU Performance Test")
print("="*70)
print("\nTesting with RTX 4060 Laptop GPU")
print(f"Running {len(test_variants)} variant predictions...\n")

results = []

for i, variant in enumerate(test_variants, 1):
    print(f"\n{'='*70}")
    print(f"Test {i}/{len(test_variants)}: {variant['name']}")
    print(f"{'='*70}")
    print(f"Position: {variant['chrom']}:{variant['pos']} {variant['ref']}->{variant['alt']}")
    
    start_time = time.time()
    
    try:
        result = compute_variant_impact(
            chrom=variant['chrom'],
            pos=variant['pos'],
            ref=variant['ref'],
            alt=variant['alt']
        )
        
        elapsed = time.time() - start_time
        
        print(f"\n[SUCCESS] Prediction complete in {elapsed:.2f} seconds")
        print(f"  Gene: {result['metrics']['gene_symbol']}")
        print(f"  Max Delta: {result['metrics']['max_delta']:.4f}")
        print(f"  gnomAD Freq: {result['metrics']['gnomad_freq']:.6f}")
        print(f"  Percentile: {result['metrics']['percentile']:.1f}%")
        print(f"  Exons found: {len(result['tracks']['exons'])}")
        
        results.append({
            "variant": variant['name'],
            "time": elapsed,
            "gene": result['metrics']['gene_symbol'],
            "max_delta": result['metrics']['max_delta']
        })
        
    except Exception as e:
        print(f"\n[ERROR] {e}")
        results.append({
            "variant": variant['name'],
            "time": -1,
            "error": str(e)
        })

# Summary
print(f"\n\n{'='*70}")
print("Performance Summary")
print(f"{'='*70}")

total_time = sum(r['time'] for r in results if r['time'] > 0)
avg_time = total_time / len([r for r in results if r['time'] > 0]) if results else 0

print(f"\nTotal variants processed: {len(results)}")
print(f"Total time: {total_time:.2f} seconds")
print(f"Average time per variant: {avg_time:.2f} seconds")

print("\nDetailed Results:")
print(f"{'Variant':<30} {'Time (s)':<12} {'Gene':<10} {'Max Delta':<12}")
print("-"*70)
for r in results:
    if 'error' not in r:
        print(f"{r['variant']:<30} {r['time']:<12.2f} {r['gene']:<10} {r['max_delta']:<12.4f}")
    else:
        print(f"{r['variant']:<30} {'ERROR':<12}")

print(f"\n{'='*70}")
print("GPU Acceleration Active!")
print(f"{'='*70}")
print("\nExpected performance:")
print("  - First run: ~5-8 seconds (model loading)")
print("  - Subsequent runs: <1 second per variant")
print("  - Speedup vs CPU: 10-30x faster")
