#!/usr/bin/env python3
"""
HIGH-QUALITY DEG METHYLATION ANALYSIS
Only uses high-confidence orthologs (E-value < 1e-20, identity > 60%)
"""
import pandas as pd

# Load mappings
brooksi_map = pd.read_csv('Results/05_transcriptomics/brooksi_deg_to_chacei.tsv', sep='\t',
                          names=['transcript_id','chacei_gene','pident','evalue','bitscore'])
elizabethae_map = pd.read_csv('Results/05_transcriptomics/elizabethae_deg_to_chacei.tsv', sep='\t',
                              names=['transcript_id','chacei_gene','pident','evalue','bitscore'])

# FILTER TO HIGH QUALITY ONLY
print("="*80)
print("HIGH-QUALITY ORTHOLOG FILTERING")
print("="*80)

quality_threshold_evalue = 1e-20
quality_threshold_pident = 60

brooksi_hq = brooksi_map[(brooksi_map['evalue'] < quality_threshold_evalue) & 
                          (brooksi_map['pident'] > quality_threshold_pident)]
elizabethae_hq = elizabethae_map[(elizabethae_map['evalue'] < quality_threshold_evalue) & 
                                  (elizabethae_map['pident'] > quality_threshold_pident)]

print(f"\nS. brooksi:")
print(f"  Total mappings: {len(brooksi_map)}")
print(f"  High-quality (E<1e-20, ID>60%): {len(brooksi_hq)} ({len(brooksi_hq)/len(brooksi_map)*100:.1f}%)")
print(f"  Median E-value: {brooksi_hq['evalue'].median():.2e}")
print(f"  Median % identity: {brooksi_hq['pident'].median():.1f}%")

print(f"\nS. elizabethae:")
print(f"  Total mappings: {len(elizabethae_map)}")
print(f"  High-quality (E<1e-20, ID>60%): {len(elizabethae_hq)} ({len(elizabethae_hq)/len(elizabethae_map)*100:.1f}%)")
print(f"  Median E-value: {elizabethae_hq['evalue'].median():.2e}")
print(f"  Median % identity: {elizabethae_hq['pident'].median():.1f}%")

# Combine high-quality mappings
all_hq = pd.concat([brooksi_hq, elizabethae_hq])
hq_genes = set(all_hq['chacei_gene'])

# Load methylation data
deg_meth = pd.read_csv('Results/05_transcriptomics/deg_direction_methylation.csv')

# Filter to high-quality orthologs only
deg_meth_hq = deg_meth[deg_meth['chacei_gene'].isin(hq_genes)].copy()

print(f"\n" + "="*80)
print("HIGH-QUALITY DEG METHYLATION ANALYSIS")
print("="*80)

print(f"\nTotal DEGs with high-quality orthologs: {len(deg_meth_hq)}")

in_data = deg_meth_hq[deg_meth_hq['methylation_status'] != 'not_in_data']
print(f"In methylation data: {len(in_data)}")

if len(in_data) > 0:
    methylated = in_data[in_data['methylation_status'] == 'methylated']
    unmethylated = in_data[in_data['methylation_status'] == 'unmethylated']
    
    print(f"Methylated: {len(methylated)} ({len(methylated)/len(in_data)*100:.1f}%)")
    print(f"Unmethylated: {len(unmethylated)} ({len(unmethylated)/len(in_data)*100:.1f}%)")
    
    print(f"\nGenome-wide: 16.1% methylated")
    print(f"High-quality DEGs: {len(methylated)/len(in_data)*100:.1f}% methylated")
    
    if len(methylated) > 0:
        depletion = 16.1 / (len(methylated)/len(in_data)*100)
        print(f"Depletion factor: {depletion:.2f}x")
    
    # Queen vs Worker
    queen = in_data[in_data['direction'] == 'Queen_up']
    worker = in_data[in_data['direction'] == 'Worker_up']
    
    if len(queen) > 0 and len(worker) > 0:
        queen_meth = queen[queen['methylation_status'] == 'methylated']
        worker_meth = worker[worker['methylation_status'] == 'methylated']
        
        print(f"\nQueen-up: {len(queen_meth)}/{len(queen)} methylated ({len(queen_meth)/len(queen)*100:.1f}%)")
        print(f"Worker-up: {len(worker_meth)}/{len(worker)} methylated ({len(worker_meth)/len(worker)*100:.1f}%)")
    
    # Export
    deg_meth_hq.to_csv('Results/06_deg_methylation/high_quality_deg_methylation.csv', index=False)
    print(f"\n✓ Exported to Results/06_deg_methylation/high_quality_deg_methylation.csv")
    
    # Show which genes these are
    if len(methylated) > 0:
        print(f"\nMethylated high-quality DEGs:")
        for _, row in methylated.iterrows():
            print(f"  {row['chacei_gene']} | {row['direction']}")
else:
    print("\n⚠️  No high-quality DEGs found in methylation data")
    print("Sample size too small for analysis")
