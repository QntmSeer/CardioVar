from msaexplorer.explore import MSA
import io

# Dummy FASTA
fasta_data = """>seq1
ATGCGTACGTTAG
>seq2
ATGCGTACGT-AG
>seq3
ATGCGTACGTTAG
"""

# Try passing string directly (if it supports file path, maybe I need to save it first)
try:
    msa = MSA(fasta_data)
    print("Initialized with string")
except Exception as e:
    print(f"Failed with string: {e}")
    # Try with file
    with open("test.fasta", "w") as f:
        f.write(fasta_data)
    try:
        msa = MSA("test.fasta")
        print("Initialized with file")
        print(f"Length: {msa.length}")
    except Exception as e:
        print(f"Failed with file: {e}")
