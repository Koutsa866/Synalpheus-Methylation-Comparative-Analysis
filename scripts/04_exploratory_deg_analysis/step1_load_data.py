#!/usr/bin/env python3
import pandas as pd

# Step 1: Load DESeq2 wholetable data
print("Loading S. brooksi data...")
brooksi = pd.read_csv('Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular', sep='\t')
print(f"  Loaded {len(brooksi)} transcripts")
print(f"  Columns: {list(brooksi.columns)}")
print(f"\nFirst 3 rows:")
print(brooksi.head(3))

print("\n" + "="*80 + "\n")

print("Loading S. elizabethae data...")
elizabethae = pd.read_csv('Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular', sep='\t')
print(f"  Loaded {len(elizabethae)} transcripts")
print(f"  Columns: {list(elizabethae.columns)}")
print(f"\nFirst 3 rows:")
print(elizabethae.head(3))

print("\n✓ Step 1 complete")
