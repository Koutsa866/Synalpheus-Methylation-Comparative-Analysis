#!/usr/bin/env python3
import pandas as pd

# Load the combined DEG file
deg_meth = pd.read_csv('Results/05_transcriptomics/deg_direction_methylation.csv')

# Load mapping files to identify species
brooksi_map = pd.read_csv('Results/05_transcriptomics/brooksi_deg_to_chacei.tsv', sep='\t', 
                          names=['transcript_id','chacei_gene','pident','evalue','bitscore'])
elizabethae_map = pd.read_csv('Results/05_transcriptomics/elizabethae_deg_to_chacei.tsv', sep='\t',
                              names=['transcript_id','chacei_gene','pident','evalue','bitscore'])

# Tag species
brooksi_genes = set(brooksi_map['chacei_gene'])
elizabethae_genes = set(elizabethae_map['chacei_gene'])

deg_meth['species'] = deg_meth['chacei_gene'].apply(
    lambda x: 'brooksi' if x in brooksi_genes else ('elizabethae' if x in elizabethae_genes else 'unknown')
)

# Separate by species
brooksi_degs = deg_meth[deg_meth['species'] == 'brooksi'].copy()
elizabethae_degs = deg_meth[deg_meth['species'] == 'elizabethae'].copy()

print("="*80)
print("SPECIES SEPARATION ANALYSIS")
print("="*80)

for species_name, df in [('S. brooksi', brooksi_degs), ('S. elizabethae', elizabethae_degs)]:
    print(f"\n{species_name}:")
    print(f"  Total DEGs: {len(df)}")
    
    in_data = df[df['methylation_status'] != 'not_in_data']
    print(f"  In methylation data: {len(in_data)}")
    
    methylated = in_data[in_data['methylation_status'] == 'methylated']
    print(f"  Methylated: {len(methylated)} ({len(methylated)/len(in_data)*100:.1f}%)")
    
    queen = in_data[in_data['direction'] == 'Queen_up']
    worker = in_data[in_data['direction'] == 'Worker_up']
    queen_meth = queen[queen['methylation_status'] == 'methylated']
    worker_meth = worker[worker['methylation_status'] == 'methylated']
    
    print(f"  Queen-up: {len(queen_meth)}/{len(queen)} methylated ({len(queen_meth)/len(queen)*100:.1f}%)")
    print(f"  Worker-up: {len(worker_meth)}/{len(worker)} methylated ({len(worker_meth)/len(worker)*100:.1f}%)")

# Export separated files
brooksi_degs.to_csv('Results/06_deg_methylation/brooksi_deg_methylation.csv', index=False)
elizabethae_degs.to_csv('Results/06_deg_methylation/elizabethae_deg_methylation.csv', index=False)

print(f"\n✓ Exported species-separated files")
