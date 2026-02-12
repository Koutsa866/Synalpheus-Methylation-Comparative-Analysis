#!/usr/bin/env python3
import pandas as pd

# Load data
brooksi = pd.read_csv('Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular', sep='\t')
elizabethae = pd.read_csv('Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular', sep='\t')

# Step 2: Classify DEGs
def classify_degs(df, species_name):
    df['regulation_direction'] = 'Non-DEG'
    df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] > 0), 'regulation_direction'] = 'Queen-upregulated'
    df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] < 0), 'regulation_direction'] = 'Worker-upregulated'
    
    counts = df['regulation_direction'].value_counts()
    print(f"\n{species_name}:")
    print(f"  Queen-upregulated: {counts.get('Queen-upregulated', 0)}")
    print(f"  Worker-upregulated: {counts.get('Worker-upregulated', 0)}")
    print(f"  Non-DEG: {counts.get('Non-DEG', 0)}")
    print(f"  Total: {len(df)}")
    
    return df

brooksi = classify_degs(brooksi, "S. brooksi")
elizabethae = classify_degs(elizabethae, "S. elizabethae")

print("\n✓ Step 2 complete")
