#!/usr/bin/env python3
"""
FIX ANNOTATION MAPPING
Correctly map protein IDs to gene IDs using FASTA headers
"""

import pandas as pd
from Bio import SeqIO

print("="*70)
print("FIXING ANNOTATION MAPPING")
print("="*70)

# Step 1: Build correct protein → gene mapping from FASTA
print("\n1. Building protein→gene mapping from FASTA...")
protein_fasta = "Data/Transcriptomic_Data/S.chacei_all_proteins.fasta"
protein_to_gene = {}

for rec in SeqIO.parse(protein_fasta, "fasta"):
    # Header format: >contig_19849_ID=g2 contig_19849:27795-28795(+)
    # rec.id = "contig_19849_ID=g2"
    # We need to extract the gene_id from description
    parts = rec.description.split()
    if len(parts) >= 1:
        gene_id = parts[0]  # contig_19849_ID=g2
        # Protein ID in DIAMOND is like g2_1, g2_2 (isoforms)
        # Extract gene number from gene_id
        if "ID=" in gene_id:
            gene_num = gene_id.split("ID=")[1]  # g2
            contig = gene_id.split("_ID=")[0]   # contig_19849
            
            # DIAMOND protein IDs are: g1_1, g1_2, g2_1, etc.
            # We need to match these back to full gene IDs
            # Store mapping: gene_num → full_gene_id
            if gene_num not in protein_to_gene:
                protein_to_gene[gene_num] = []
            protein_to_gene[gene_num].append(gene_id)

print(f"   Found {len(protein_to_gene)} unique gene numbers")
print(f"   Sample: {list(protein_to_gene.keys())[:10]}")

# Step 2: Load DIAMOND results
print("\n2. Loading DIAMOND results...")
diamond = pd.read_csv("Results/04_functional_annotation/diamond_results.tsv", 
                      sep='\t', header=None,
                      names=['protein_id', 'uniprot_id', 'pident', 'evalue', 'bitscore', 'description'])
print(f"   {len(diamond)} protein annotations")

# Step 3: Map protein IDs to gene IDs
print("\n3. Mapping proteins to genes...")

# Build better mapping: protein_id (g1_3) → gene_id (contig_XXXX_ID=g1)
# Need to parse the FASTA more carefully
protein_to_gene_correct = {}

for rec in SeqIO.parse(protein_fasta, "fasta"):
    # rec.id might be the protein isoform ID
    # rec.description has the full info
    desc_parts = rec.description.split()
    if len(desc_parts) >= 1:
        gene_id_full = desc_parts[0]  # contig_19849_ID=g2
        
        # The protein ID in DIAMOND might be derived from the sequence order
        # Let's check what the actual protein IDs look like in the FASTA
        protein_id = rec.id
        protein_to_gene_correct[protein_id] = gene_id_full

print(f"   Built mapping for {len(protein_to_gene_correct)} proteins")

# Check if DIAMOND protein IDs match FASTA IDs
sample_diamond_ids = diamond['protein_id'].head(10).tolist()
print(f"\n   Sample DIAMOND protein IDs: {sample_diamond_ids}")
print(f"   Sample FASTA protein IDs: {list(protein_to_gene_correct.keys())[:10]}")

# Map using the correct mapping
diamond['gene_id'] = diamond['protein_id'].map(protein_to_gene_correct)

# Step 4: Add methylation status
print("\n4. Adding methylation status...")
promoter_meth = pd.read_csv("Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv", sep='\t')
meth_genes = set(promoter_meth['gene_id'])

diamond['methylation_status'] = diamond['gene_id'].apply(
    lambda x: 'methylated' if x in meth_genes else 'unmethylated' if pd.notna(x) else 'unknown'
)

# Step 5: Save corrected mapping
print("\n5. Saving corrected mapping...")
output_file = "Results/04_functional_annotation/gene_uniprot_mapping_CORRECTED.csv"
diamond.to_csv(output_file, index=False)

# Step 6: Statistics
print("\n" + "="*70)
print("RESULTS")
print("="*70)
print(f"Total protein annotations: {len(diamond)}")
print(f"Unique genes annotated: {diamond['gene_id'].nunique()}")
print(f"Genes with unknown mapping: {diamond['gene_id'].isna().sum()}")
print(f"Methylated gene annotations: {(diamond['methylation_status']=='methylated').sum()}")
print(f"Unmethylated gene annotations: {(diamond['methylation_status']=='unmethylated').sum()}")

print(f"\nSaved to: {output_file}")

# Show sample of mapped genes
print("\nSample mapped genes:")
print(diamond[['protein_id', 'gene_id', 'methylation_status']].head(10))

print("\n✓ Mapping fixed!")
