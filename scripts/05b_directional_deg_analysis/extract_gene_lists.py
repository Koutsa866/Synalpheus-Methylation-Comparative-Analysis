#!/usr/bin/env python3
"""
Extract Gene IDs for Functional Analysis
Extracts S. chacei gene IDs for up-regulated and down-regulated methylated genes.

Purpose:
  Generate gene lists for downstream functional annotation, BLAST searches,
  or sequence extraction.

Inputs:
  - Results/05_ortholog_analysis/chacei_vs_brooksi.tsv
  - Results/05_ortholog_analysis/chacei_vs_elizabethae.tsv
  - Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular
  - Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular

Outputs:
  - Results/05b_directional_analysis/gene_lists/elizabethae_upregulated_methylated_genes.csv
  - Results/05b_directional_analysis/gene_lists/elizabethae_downregulated_methylated_genes.csv
  - Results/05b_directional_analysis/gene_lists/brooksi_upregulated_methylated_genes.csv
  - Results/05b_directional_analysis/gene_lists/brooksi_downregulated_methylated_genes.csv

Author: Philip Koutsaftis
Advisor: Dr. Solomon Chak
Institution: Denison University
Date: 2025
"""

import pandas as pd
import os

def extract_directional_genes(species, ortho_path, dea_path):
    """
    Extract up-regulated and down-regulated methylated gene IDs.
    
    Args:
        species (str): Species name ('brooksi' or 'elizabethae')
        ortho_path (str): Path to DIAMOND orthology results
        dea_path (str): Path to DESeq2 whole-table results
    
    Returns:
        tuple: (up_regulated_df, down_regulated_df)
    """
    print(f"\nExtracting genes for S. {species}...")
    
    # Load DEA data
    dea_df = pd.read_csv(dea_path, sep='\t', index_col=0, low_memory=False)
    padj_col = dea_df.columns[8]
    lfc_col = dea_df.columns[5]
    
    # Load orthology data
    ortho_df = pd.read_csv(ortho_path, sep='\t', 
                           names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore'])
    
    # Extract methylation status
    ortho_df['Methylation'] = ortho_df['qseqid'].apply(lambda x: x.split('|')[1])
    methylated = ortho_df[ortho_df['Methylation'] == 'Methylated'].copy()
    
    # Map expression data
    methylated['padj'] = methylated['sseqid'].map(dea_df[padj_col])
    methylated['log2FC'] = methylated['sseqid'].map(dea_df[lfc_col])
    
    # Filter for up-regulated
    up_regulated = methylated[(methylated['padj'] < 0.05) & (methylated['log2FC'] > 0)].copy()
    up_regulated['chacei_gene_id'] = up_regulated['qseqid'].apply(lambda x: x.split('|')[0])
    
    # Filter for down-regulated
    down_regulated = methylated[(methylated['padj'] < 0.05) & (methylated['log2FC'] < 0)].copy()
    down_regulated['chacei_gene_id'] = down_regulated['qseqid'].apply(lambda x: x.split('|')[0])
    
    print(f"  Up-regulated methylated genes: {len(up_regulated)}")
    print(f"  Down-regulated methylated genes: {len(down_regulated)}")
    
    return up_regulated, down_regulated

# --- EXECUTE ANALYSIS ---
if __name__ == "__main__":
    print("\n" + "="*60)
    print("GENE LIST EXTRACTION")
    print("="*60)
    
    # Create output directory
    os.makedirs('Results/05b_directional_analysis/gene_lists', exist_ok=True)
    
    # Extract for S. elizabethae
    eliz_up, eliz_down = extract_directional_genes(
        "elizabethae",
        "Results/05_ortholog_analysis/chacei_vs_elizabethae.tsv",
        "Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular"
    )
    
    # Extract for S. brooksi
    brook_up, brook_down = extract_directional_genes(
        "brooksi",
        "Results/05_ortholog_analysis/chacei_vs_brooksi.tsv",
        "Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular"
    )
    
    # Save results
    print("\n" + "="*60)
    print("Saving gene lists...")
    
    # S. elizabethae
    eliz_up[['chacei_gene_id', 'sseqid', 'log2FC', 'padj', 'pident']].to_csv(
        'Results/05b_directional_analysis/gene_lists/elizabethae_upregulated_methylated_genes.csv',
        index=False
    )
    print("  ✓ elizabethae_upregulated_methylated_genes.csv")
    
    eliz_down[['chacei_gene_id', 'sseqid', 'log2FC', 'padj', 'pident']].to_csv(
        'Results/05b_directional_analysis/gene_lists/elizabethae_downregulated_methylated_genes.csv',
        index=False
    )
    print("  ✓ elizabethae_downregulated_methylated_genes.csv")
    
    # S. brooksi
    brook_up[['chacei_gene_id', 'sseqid', 'log2FC', 'padj', 'pident']].to_csv(
        'Results/05b_directional_analysis/gene_lists/brooksi_upregulated_methylated_genes.csv',
        index=False
    )
    print("  ✓ brooksi_upregulated_methylated_genes.csv")
    
    brook_down[['chacei_gene_id', 'sseqid', 'log2FC', 'padj', 'pident']].to_csv(
        'Results/05b_directional_analysis/gene_lists/brooksi_downregulated_methylated_genes.csv',
        index=False
    )
    print("  ✓ brooksi_downregulated_methylated_genes.csv")
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"S. elizabethae: {len(eliz_up)} up, {len(eliz_down)} down")
    print(f"S. brooksi: {len(brook_up)} up, {len(brook_down)} down")
    print("\n✓ Gene list extraction complete!")
