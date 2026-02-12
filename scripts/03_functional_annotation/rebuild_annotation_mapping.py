#!/usr/bin/env python3
"""
REBUILD ANNOTATION MAPPING
Use BED file + protein FASTA order to correctly map protein IDs to gene IDs
"""

import pandas as pd
from Bio import SeqIO
from collections import defaultdict

print("="*70)
print("REBUILDING ANNOTATION MAPPING")
print("="*70)

# Step 1: Load BED file to get gene coordinates
print("\n1. Loading gene coordinates from BED file...")
bed_file = "archive/old_results/results (1)/gene_promoters_fixed.bed"
bed = pd.read_csv(bed_file, sep='\t', header=None,
                  names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand'])

# Create full gene IDs: contig_ID=gX
bed['full_gene_id'] = bed['chrom'] + '_' + bed['gene_id']
print(f"   Loaded {len(bed)} gene entries")
print(f"   Sample: {bed['full_gene_id'].head(3).tolist()}")

# Step 2: Load protein FASTA to get protein order
print("\n2. Loading protein sequences...")
protein_fasta = "archive/old_results/results (1)/filtered_genes_proteins_unique.faa"
proteins = []
for rec in SeqIO.parse(protein_fasta, "fasta"):
    proteins.append(rec.id)

print(f"   Loaded {len(proteins)} proteins")
print(f"   Sample: {proteins[:5]}")

# Step 3: Build mapping based on gene numbering
print("\n3. Building protein→gene mapping...")

# Group BED entries by gene number
gene_groups = defaultdict(list)
for _, row in bed.iterrows():
    gene_num = row['gene_id'].replace('ID=', '')  # g1, g2, etc.
    gene_groups[gene_num].append(row['full_gene_id'])

print(f"   Found {len(gene_groups)} gene numbers")
print(f"   g1 has {len(gene_groups['g1'])} instances")

# Map proteins to genes
# Protein ID format: g1_1, g1_2, g2_1, etc.
# This means: gene_number_isoform
protein_to_gene = {}
gene_counters = defaultdict(int)

for protein_id in proteins:
    # Parse protein ID: g1_1 → gene_num=g1, isoform=1
    parts = protein_id.split('_')
    if len(parts) == 2:
        gene_num = parts[0]  # g1
        isoform = int(parts[1])  # 1
        
        # Get the list of genes with this number
        if gene_num in gene_groups:
            gene_list = gene_groups[gene_num]
            # Use counter to assign to correct gene instance
            idx = gene_counters[gene_num]
            if idx < len(gene_list):
                protein_to_gene[protein_id] = gene_list[idx]
                # Increment counter when we've seen all isoforms for this gene
                # (This assumes isoforms are grouped together)
                gene_counters[gene_num] += 1

print(f"   Mapped {len(protein_to_gene)} proteins to genes")

# Step 4: Load DIAMOND results and apply mapping
print("\n4. Loading DIAMOND results...")
diamond = pd.read_csv("Results/04_functional_annotation/diamond_results.tsv",
                      sep='\t', header=None,
                      names=['protein_id', 'uniprot_id', 'pident', 'evalue', 'bitscore', 'description'])
print(f"   {len(diamond)} annotations")

# Apply mapping
diamond['gene_id'] = diamond['protein_id'].map(protein_to_gene)

# Step 5: Add methylation status
print("\n5. Adding methylation status...")
promoter_meth = pd.read_csv("Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv", sep='\t')
meth_genes = set(promoter_meth['gene_id'])

diamond['methylation_status'] = diamond['gene_id'].apply(
    lambda x: 'methylated' if x in meth_genes else 'unmethylated' if pd.notna(x) else 'unknown'
)

# Step 6: Save
print("\n6. Saving corrected mapping...")
output_file = "Results/04_functional_annotation/gene_uniprot_mapping_FIXED.csv"
diamond.to_csv(output_file, index=False)

# Statistics
print("\n" + "="*70)
print("RESULTS")
print("="*70)
print(f"Total annotations: {len(diamond)}")
print(f"Unique genes annotated: {diamond['gene_id'].nunique()}")
print(f"Successfully mapped: {diamond['gene_id'].notna().sum()}")
print(f"Failed to map: {diamond['gene_id'].isna().sum()}")
print(f"Methylated: {(diamond['methylation_status']=='methylated').sum()}")
print(f"Unmethylated: {(diamond['methylation_status']=='unmethylated').sum()}")

print(f"\nSaved to: {output_file}")
print("\nSample mappings:")
print(diamond[['protein_id', 'gene_id', 'methylation_status']].head(10))

print("\n✓ Mapping rebuilt!")
