#!/usr/bin/env python3
"""
ID Debugger: Diagnose why DIAMOND hits don't match DEA IDs
"""

import pandas as pd

print("="*60)
print("ID MISMATCH DEBUGGER")
print("="*60)

# 1. Check DIAMOND output IDs
print("\n1. DIAMOND Hit Subject IDs (first 5):")
with open("Results/05_ortholog_analysis/chacei_vs_brooksi.tsv", 'r') as f:
    for i in range(5):
        line = f.readline().strip().split('\t')
        print(f"   [{i+1}] sseqid: '{line[1]}'")

# 2. Check DEA table IDs
print("\n2. DEA Table Transcript IDs (first 5):")
dea_df = pd.read_csv("Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular", 
                     sep='\t', nrows=5)
for i in range(5):
    print(f"   [{i+1}] ID: '{dea_df.iloc[i, 0]}'")

# 3. Check protein file IDs
print("\n3. Protein File IDs (first 5):")
with open("Data/transcriptomics/transcriptomes/S.brooksi_proteins.faa", 'r') as f:
    count = 0
    for line in f:
        if line.startswith('>'):
            print(f"   [{count+1}] Header: '{line.strip()}'")
            count += 1
            if count >= 5:
                break

# 4. Check DNA transcriptome IDs
print("\n4. DNA Transcriptome IDs (first 5):")
with open("Data/transcriptomics/transcriptomes/S.brooksi_brain_filtered_transcriptome.fasta", 'r') as f:
    count = 0
    for line in f:
        if line.startswith('>'):
            print(f"   [{count+1}] Header: '{line.strip()}'")
            count += 1
            if count >= 5:
                break

print("\n" + "="*60)
print("DIAGNOSIS:")
print("Compare the ID formats above to identify the mismatch.")
print("="*60)
