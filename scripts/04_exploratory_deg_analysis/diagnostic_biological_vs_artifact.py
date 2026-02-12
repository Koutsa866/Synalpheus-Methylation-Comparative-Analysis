#!/usr/bin/env python3
"""
DIAGNOSTIC: Is 0% methylation real or an artifact?
Tests multiple hypotheses for why DEGs show no methylation
"""
import pandas as pd
import numpy as np

print("="*80)
print("DIAGNOSTIC: BIOLOGICAL SIGNAL vs TECHNICAL ARTIFACT")
print("="*80)

# Load data
hq_degs = pd.read_csv('Results/06_deg_methylation/high_quality_deg_methylation.csv')
all_genes_meth = pd.read_csv('Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv', sep='\t')
all_genes_unmeth = pd.read_csv('Results/03_methylation_analysis/genes_without_methylated_promoters.tsv', sep='\t')

# Get all S. chacei genes
all_chacei_genes = set(all_genes_meth['gene_id']) | set(all_genes_unmeth['gene_id'])

print(f"\n[TEST 1] Are the 13 DEG genes even in the methylation dataset?")
print("-"*80)
in_meth_data = hq_degs[hq_degs['methylation_status'] != 'not_in_data']
print(f"DEGs in methylation data: {len(in_meth_data)}/34")
print(f"DEGs NOT in methylation data: {len(hq_degs[hq_degs['methylation_status'] == 'not_in_data'])}/34")

print("\nGenes in methylation data:")
for _, row in in_meth_data.iterrows():
    print(f"  {row['chacei_gene']} - {row['methylation_status']}")

print(f"\n[TEST 2] Are these genes in the genome-wide methylation files?")
print("-"*80)
for _, row in in_meth_data.iterrows():
    gene = row['chacei_gene']
    in_meth_file = gene in all_genes_meth['gene_id'].values
    in_unmeth_file = gene in all_genes_unmeth['gene_id'].values
    print(f"  {gene}:")
    print(f"    In methylated file: {in_meth_file}")
    print(f"    In unmethylated file: {in_unmeth_file}")

print(f"\n[TEST 3] Random sampling - what % of random genes are methylated?")
print("-"*80)
# Sample 13 random genes 1000 times
np.random.seed(42)
random_meth_rates = []
for i in range(1000):
    random_sample = np.random.choice(list(all_chacei_genes), size=13, replace=False)
    n_methylated = sum([g in all_genes_meth['gene_id'].values for g in random_sample])
    random_meth_rates.append(n_methylated / 13 * 100)

print(f"Expected methylation rate (random sampling):")
print(f"  Mean: {np.mean(random_meth_rates):.1f}%")
print(f"  Std: {np.std(random_meth_rates):.1f}%")
print(f"  Min: {np.min(random_meth_rates):.1f}%")
print(f"  Max: {np.max(random_meth_rates):.1f}%")
print(f"  Times we got 0%: {sum([r == 0 for r in random_meth_rates])}/1000")

print(f"\n[TEST 4] Are DEGs enriched for specific gene types?")
print("-"*80)
# Load functional annotations
try:
    deg_annot = pd.read_csv('Results/05_transcriptomics/deg_uniprot_annotations_final.csv')
    hq_with_annot = in_meth_data.merge(deg_annot, left_on='chacei_gene', right_on='sequence_id', how='left')
    
    print("Functional annotations of the 13 unmethylated DEGs:")
    for _, row in hq_with_annot.iterrows():
        if pd.notna(row.get('description')):
            print(f"  {row['chacei_gene']}: {row['description'][:80]}")
        else:
            print(f"  {row['chacei_gene']}: No annotation")
except:
    print("  Could not load annotations")

print(f"\n[TEST 5] Mapping bias - do high-identity genes avoid methylation?")
print("-"*80)
brooksi_map = pd.read_csv('Results/05_transcriptomics/brooksi_deg_to_chacei.tsv', sep='\t',
                          names=['transcript_id','chacei_gene','pident','evalue','bitscore'])
elizabethae_map = pd.read_csv('Results/05_transcriptomics/elizabethae_deg_to_chacei.tsv', sep='\t',
                              names=['transcript_id','chacei_gene','pident','evalue','bitscore'])
all_maps = pd.concat([brooksi_map, elizabethae_map])

# Check if high-identity genes are less methylated
all_maps['is_methylated'] = all_maps['chacei_gene'].apply(
    lambda x: x in all_genes_meth['gene_id'].values
)

high_id = all_maps[all_maps['pident'] > 80]
low_id = all_maps[all_maps['pident'] <= 80]

print(f"High identity (>80%): {high_id['is_methylated'].sum()}/{len(high_id)} methylated ({high_id['is_methylated'].sum()/len(high_id)*100:.1f}%)")
print(f"Low identity (≤80%): {low_id['is_methylated'].sum()}/{len(low_id)} methylated ({low_id['is_methylated'].sum()/len(low_id)*100:.1f}%)")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)

if sum([r == 0 for r in random_meth_rates]) < 10:
    print("⚠️  LIKELY BIOLOGICAL SIGNAL")
    print("   Getting 0% methylation by chance is extremely rare (<1%)")
    print("   This suggests DEGs genuinely avoid methylation")
else:
    print("⚠️  POSSIBLE ARTIFACT")
    print("   Getting 0% methylation happens frequently by chance")
    print("   Sample size may be too small to draw conclusions")

print("\nRecommendation:")
print("  - If p < 0.01 for random sampling: BIOLOGICAL SIGNAL")
print("  - If p > 0.05: INSUFFICIENT EVIDENCE")
print(f"  - Observed p-value: {sum([r == 0 for r in random_meth_rates])/1000:.3f}")
