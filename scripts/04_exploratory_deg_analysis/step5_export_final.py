#!/usr/bin/env python3
import pandas as pd

print("Step 5: Combining data and exporting final dataset")

# Load merged data
brooksi = pd.read_csv('Results/05_transcriptomics/brooksi_merged.csv')
elizabethae = pd.read_csv('Results/05_transcriptomics/elizabethae_merged.csv')

# Add species column
brooksi['species'] = 'brooksi'
elizabethae['species'] = 'elizabethae'

# Select and reorder columns
columns = ['species', 'transcript_id', 'uniprot_id', 'pident', 'evalue', 'description', 'is_te', 'log2FoldChange', 'regulation_direction']

brooksi_final = brooksi[columns]
elizabethae_final = elizabethae[columns]

# Combine both species
combined = pd.concat([brooksi_final, elizabethae_final], ignore_index=True)

print(f"\nFinal dataset:")
print(f"  Total transcripts: {len(combined)}")
print(f"  S. brooksi: {len(brooksi_final)}")
print(f"  S. elizabethae: {len(elizabethae_final)}")

print(f"\nBreakdown by regulation:")
print(combined.groupby(['species', 'regulation_direction']).size())

print(f"\nTransposable elements:")
print(f"  Total TEs: {combined['is_te'].sum()}")
print(f"  Total non-TEs: {(~combined['is_te']).sum()}")

# Export
output_file = 'Results/05_transcriptomics/gene_functions_with_regulation_wholetable.csv'
combined.to_csv(output_file, index=False)
print(f"\n✓ Exported to: {output_file}")

print("\n✓ Step 5 complete - Pipeline finished!")
