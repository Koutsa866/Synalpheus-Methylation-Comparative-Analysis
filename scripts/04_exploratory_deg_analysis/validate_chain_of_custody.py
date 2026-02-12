#!/usr/bin/env python3
"""
DATA VALIDATION: Chain of Custody for DEG Methylation Analysis
Implements 4-step validation checklist before any analysis
"""
import pandas as pd
import numpy as np

print("="*80)
print("CHAIN OF CUSTODY VALIDATION")
print("="*80)

# ============================================================================
# STEP 1: SPOT-CHECK VALIDATION (Manual Verification)
# ============================================================================
print("\n[STEP 1] SPOT-CHECK VALIDATION")
print("-"*80)

deg_meth = pd.read_csv('Results/05_transcriptomics/deg_direction_methylation.csv')
print(f"Total entries in deg_direction_methylation.csv: {len(deg_meth)}")

# Pick 3 random genes for manual verification
sample_genes = deg_meth[deg_meth['methylation_status'] != 'not_in_data'].sample(min(3, len(deg_meth)))
print("\n3 GENES FOR MANUAL VERIFICATION:")
for idx, row in sample_genes.iterrows():
    print(f"\n  Gene: {row['chacei_gene']}")
    print(f"    Direction: {row['direction']}")
    print(f"    Methylation: {row['methylation_status']}")
    print(f"    → Verify in: Data/methylation/bedmethyl_50p_mC_cov30.bed")
    print(f"    → Verify in: Data/transcriptomics/deseq2_results/")

# ============================================================================
# STEP 2: DIAMOND ALIGNMENT QUALITY
# ============================================================================
print("\n\n[STEP 2] DIAMOND ALIGNMENT QUALITY")
print("-"*80)

brooksi_map = pd.read_csv('Results/05_transcriptomics/brooksi_deg_to_chacei.tsv', sep='\t',
                          names=['transcript_id','chacei_gene','pident','evalue','bitscore'])
elizabethae_map = pd.read_csv('Results/05_transcriptomics/elizabethae_deg_to_chacei.tsv', sep='\t',
                              names=['transcript_id','chacei_gene','pident','evalue','bitscore'])

for species, df in [('S. brooksi', brooksi_map), ('S. elizabethae', elizabethae_map)]:
    print(f"\n{species}:")
    print(f"  Total mappings: {len(df)}")
    print(f"  E-value < 1e-10: {(df['evalue'] < 1e-10).sum()} ({(df['evalue'] < 1e-10).sum()/len(df)*100:.1f}%)")
    print(f"  Percent identity > 70%: {(df['pident'] > 70).sum()} ({(df['pident'] > 70).sum()/len(df)*100:.1f}%)")
    print(f"  Median E-value: {df['evalue'].median():.2e}")
    print(f"  Median % identity: {df['pident'].median():.1f}%")
    
    # Flag low-quality mappings
    low_quality = df[(df['evalue'] > 1e-10) | (df['pident'] < 70)]
    if len(low_quality) > 0:
        print(f"  ⚠️  WARNING: {len(low_quality)} low-quality mappings detected")

# ============================================================================
# STEP 3: MAPPING BIAS CHECK
# ============================================================================
print("\n\n[STEP 3] MAPPING BIAS CHECK")
print("-"*80)

# Load full transcriptome sizes
brooksi_total = 411  # from brooksi_deg_ids.txt
elizabethae_total = 275  # from elizabethae_deg_ids.txt

brooksi_mapped = len(brooksi_map)
elizabethae_mapped = len(elizabethae_map)

print(f"\nS. brooksi:")
print(f"  Total DEGs: {brooksi_total}")
print(f"  Mapped to S. chacei: {brooksi_mapped} ({brooksi_mapped/brooksi_total*100:.1f}%)")

print(f"\nS. elizabethae:")
print(f"  Total DEGs: {elizabethae_total}")
print(f"  Mapped to S. chacei: {elizabethae_mapped} ({elizabethae_mapped/elizabethae_total*100:.1f}%)")

# Check if methylated genes map differently
in_meth_data = deg_meth[deg_meth['methylation_status'] != 'not_in_data']
methylated = in_meth_data[in_meth_data['methylation_status'] == 'methylated']
unmethylated = in_meth_data[in_meth_data['methylation_status'] == 'unmethylated']

print(f"\nMapping to methylation data:")
print(f"  Total mapped DEGs: {len(deg_meth)}")
print(f"  In methylation dataset: {len(in_meth_data)} ({len(in_meth_data)/len(deg_meth)*100:.1f}%)")
print(f"  Methylated: {len(methylated)}")
print(f"  Unmethylated: {len(unmethylated)}")

# ============================================================================
# STEP 4: SANITY CHECK TABLE
# ============================================================================
print("\n\n[STEP 4] SANITY CHECK TABLE")
print("-"*80)

# Load genome-wide methylation data
genes_meth = pd.read_csv('Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv', sep='\t')
genes_unmeth = pd.read_csv('Results/03_methylation_analysis/genes_without_methylated_promoters.tsv', sep='\t')
total_genes = len(genes_meth) + len(genes_unmeth)

sanity_table = pd.DataFrame({
    'Metric': [
        'Input DEGs (combined)',
        'S. brooksi DEGs',
        'S. elizabethae DEGs',
        'Successfully mapped to S. chacei',
        'In S. chacei methylation data',
        'Methylated DEGs found',
        'Unmethylated DEGs found',
        'Total S. chacei genes in genome',
        'Total S. chacei genes methylated',
        'Genome-wide methylation rate',
        'DEG methylation rate',
        'Depletion factor'
    ],
    'Value': [
        f"{len(deg_meth)}",
        f"{brooksi_mapped} ({brooksi_mapped/brooksi_total*100:.1f}%)",
        f"{elizabethae_mapped} ({elizabethae_mapped/elizabethae_total*100:.1f}%)",
        f"{len(deg_meth)} ({len(deg_meth)/(brooksi_total+elizabethae_total)*100:.1f}%)",
        f"{len(in_meth_data)} ({len(in_meth_data)/len(deg_meth)*100:.1f}%)",
        f"{len(methylated)}",
        f"{len(unmethylated)}",
        f"{total_genes:,}",
        f"{len(genes_meth):,}",
        f"{len(genes_meth)/total_genes*100:.1f}%",
        f"{len(methylated)/len(in_meth_data)*100:.1f}%",
        f"{(len(genes_meth)/total_genes) / (len(methylated)/len(in_meth_data)):.2f}x"
    ]
})

print(sanity_table.to_string(index=False))

# ============================================================================
# METADATA LOG
# ============================================================================
print("\n\n[METADATA LOG]")
print("-"*80)
print("Files used in this analysis:")
print("  Transcriptome DEGs:")
print("    - Data/transcriptomics/deseq2_results/S.brooksi_brain_DEseq2_result.tabular")
print("    - Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEseq2_result.tabular")
print("  Methylome:")
print("    - Data/methylation/bedmethyl_50p_mC_cov30.bed")
print("    - Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv")
print("    - Results/03_methylation_analysis/genes_without_methylated_promoters.tsv")
print("  Mapping:")
print("    - Results/05_transcriptomics/brooksi_deg_to_chacei.tsv")
print("    - Results/05_transcriptomics/elizabethae_deg_to_chacei.tsv")
print("  Integration:")
print("    - Results/05_transcriptomics/deg_direction_methylation.csv")

# ============================================================================
# VALIDATION SUMMARY
# ============================================================================
print("\n\n" + "="*80)
print("VALIDATION SUMMARY")
print("="*80)

warnings = []

# Check alignment quality
if (brooksi_map['evalue'] > 1e-10).sum() > len(brooksi_map) * 0.1:
    warnings.append("⚠️  >10% of S. brooksi mappings have E-value > 1e-10")
if (elizabethae_map['evalue'] > 1e-10).sum() > len(elizabethae_map) * 0.1:
    warnings.append("⚠️  >10% of S. elizabethae mappings have E-value > 1e-10")

# Check mapping rate
if len(in_meth_data)/len(deg_meth) < 0.3:
    warnings.append(f"⚠️  Only {len(in_meth_data)/len(deg_meth)*100:.1f}% of DEGs in methylation data")

# Check sample sizes
queen = in_meth_data[in_meth_data['direction'] == 'Queen_up']
worker = in_meth_data[in_meth_data['direction'] == 'Worker_up']
if len(queen) < 10:
    warnings.append(f"⚠️  Small Queen sample size (n={len(queen)})")
if len(worker) < 20:
    warnings.append(f"⚠️  Small Worker sample size (n={len(worker)})")

if warnings:
    print("\nWARNINGS:")
    for w in warnings:
        print(f"  {w}")
else:
    print("\n✓ All validation checks passed")

print("\n" + "="*80)
print("NEXT STEPS:")
print("  1. Manually verify the 3 spot-check genes listed above")
print("  2. Review alignment quality statistics")
print("  3. If all checks pass, proceed with figure generation")
print("="*80)
