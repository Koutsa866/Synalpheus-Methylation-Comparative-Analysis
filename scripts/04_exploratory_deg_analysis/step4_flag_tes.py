#!/usr/bin/env python3
import pandas as pd
import re

print("Step 4: Flagging Transposable Elements and merging data")

# Load DESeq2 data with classifications
brooksi = pd.read_csv('Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular', sep='\t')
elizabethae = pd.read_csv('Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular', sep='\t')

# Classify DEGs
brooksi['regulation_direction'] = 'Non-DEG'
brooksi.loc[(brooksi['padj'] < 0.05) & (brooksi['log2FoldChange'] > 0), 'regulation_direction'] = 'Queen-upregulated'
brooksi.loc[(brooksi['padj'] < 0.05) & (brooksi['log2FoldChange'] < 0), 'regulation_direction'] = 'Worker-upregulated'

elizabethae['regulation_direction'] = 'Non-DEG'
elizabethae.loc[(elizabethae['padj'] < 0.05) & (elizabethae['log2FoldChange'] > 0), 'regulation_direction'] = 'Queen-upregulated'
elizabethae.loc[(elizabethae['padj'] < 0.05) & (elizabethae['log2FoldChange'] < 0), 'regulation_direction'] = 'Worker-upregulated'

# Load UniProt annotations
brooksi_annot = pd.read_csv('Results/05_transcriptomics/brooksi_uniprot.tsv', sep='\t', header=None,
                            names=['transcript_id', 'uniprot_id', 'pident', 'length', 'evalue', 'bitscore', 'description'])
elizabethae_annot = pd.read_csv('Results/05_transcriptomics/elizabethae_uniprot.tsv', sep='\t', header=None,
                                names=['transcript_id', 'uniprot_id', 'pident', 'length', 'evalue', 'bitscore', 'description'])

# Flag TEs
te_keywords = r'transpos|retrotrans|gag|pol polyprotein|reverse transcriptase|tigger|mariner|line|sine|ltr|retrovirus|copia|gypsy|helitron|integrase'

def flag_tes(df):
    df['is_te'] = df['description'].str.contains(te_keywords, case=False, na=False)
    return df

brooksi_annot = flag_tes(brooksi_annot)
elizabethae_annot = flag_tes(elizabethae_annot)

print(f"\nS. brooksi:")
print(f"  Total annotated: {len(brooksi_annot)}")
print(f"  TEs: {brooksi_annot['is_te'].sum()}")
print(f"  Non-TEs: {(~brooksi_annot['is_te']).sum()}")

print(f"\nS. elizabethae:")
print(f"  Total annotated: {len(elizabethae_annot)}")
print(f"  TEs: {elizabethae_annot['is_te'].sum()}")
print(f"  Non-TEs: {(~elizabethae_annot['is_te']).sum()}")

# Merge with DESeq2 data
brooksi_merged = brooksi.merge(brooksi_annot, left_index=True, right_on='transcript_id', how='inner')
elizabethae_merged = elizabethae.merge(elizabethae_annot, left_index=True, right_on='transcript_id', how='inner')

print(f"\nMerged data:")
print(f"  S. brooksi: {len(brooksi_merged)} annotated transcripts with DEG info")
print(f"  S. elizabethae: {len(elizabethae_merged)} annotated transcripts with DEG info")

# Save intermediate results
brooksi_merged.to_csv('Results/05_transcriptomics/brooksi_merged.csv', index=False)
elizabethae_merged.to_csv('Results/05_transcriptomics/elizabethae_merged.csv', index=False)

print("\n✓ Step 4 complete")
