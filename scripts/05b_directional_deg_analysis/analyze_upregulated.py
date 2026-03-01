#!/usr/bin/env python3
"""
Directional DEG Analysis: Up-Regulated Genes
Tests if S. chacei promoter methylation predicts UP-REGULATION in target species.

Biological Question:
  Does promoter methylation in S. chacei predict which genes are ACTIVATED 
  (up-regulated) in the social brains of S. elizabethae and S. brooksi?

Method:
  - Classification: UP = padj < 0.05 AND log2FoldChange > 0
  - Background: All other genes (down-regulated + non-significant)
  - Statistical Test: Chi-square (or Fisher's exact if cell count < 5)

Inputs:
  - Results/05_ortholog_analysis/chacei_vs_brooksi.tsv
  - Results/05_ortholog_analysis/chacei_vs_elizabethae.tsv
  - Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular
  - Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular

Outputs:
  - Results/05b_directional_analysis/upregulated/brooksi_upregulated_contingency.csv
  - Results/05b_directional_analysis/upregulated/elizabethae_upregulated_contingency.csv
  - Results/05b_directional_analysis/upregulated/upregulated_summary_statistics.csv

Author: Philip Koutsaftis
Advisor: Dr. Solomon Chak
Institution: Denison University
Date: 2025
"""

import pandas as pd
import os
from scipy.stats import chi2_contingency, fisher_exact

def run_upregulated_analysis(species, ortho_path, dea_path):
    """
    Analyze enrichment of methylated genes in UP-REGULATED genes.
    
    Args:
        species (str): Species name ('brooksi' or 'elizabethae')
        ortho_path (str): Path to DIAMOND orthology results
        dea_path (str): Path to DESeq2 whole-table results
    
    Returns:
        dict: Analysis results including p-value, odds ratio, contingency table
    """
    print(f"\n{'='*60}")
    print(f"ANALYSIS: S. {species.upper()}")
    print(f"{'='*60}")

    # 1. Load DEA Table
    dea_df = pd.read_csv(dea_path, sep='\t', index_col=0, low_memory=False)
    
    # Column indices after index_col=0
    padj_col = dea_df.columns[8]   # Column 9 (padj)
    lfc_col = dea_df.columns[5]    # Column 6 (log2FoldChange)

    # Convert to numeric, handle NaN
    dea_df[padj_col] = pd.to_numeric(dea_df[padj_col], errors='coerce').fillna(1.0)
    dea_df[lfc_col] = pd.to_numeric(dea_df[lfc_col], errors='coerce').fillna(0.0)

    print(f"Loaded DEA data: {len(dea_df)} transcripts")
    print(f"  DEGs (padj < 0.05): {(dea_df[padj_col] < 0.05).sum()}")
    print(f"  Up-regulated (padj < 0.05 & log2FC > 0): {((dea_df[padj_col] < 0.05) & (dea_df[lfc_col] > 0)).sum()}")

    # 2. Load DIAMOND Orthology Results
    ortho_cols = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
    ortho_df = pd.read_csv(ortho_path, sep='\t', names=ortho_cols)
    print(f"Loaded DIAMOND hits: {len(ortho_df)} orthologs")

    # 3. Filter for "Successful Mappers Only"
    mapped_ids = set(dea_df.index.tolist())
    final_df = ortho_df[ortho_df['sseqid'].isin(mapped_ids)].copy()
    print(f"Successful mappers: {len(final_df)} genes")

    # 4. Extract Methylation Status from FASTA headers (gene_id|Status)
    final_df['Methylation'] = final_df['qseqid'].apply(lambda x: x.split('|')[1])

    # 5. Map Up-regulated Status
    deg_status_map = dict(zip(dea_df.index, zip(dea_df[padj_col], dea_df[lfc_col])))
    
    def classify_expression(transcript_id):
        padj, lfc = deg_status_map.get(transcript_id, (1.0, 0.0))
        if padj < 0.05 and lfc > 0:  # UP-REGULATED
            return 'UP'
        else:
            return 'Background'
    
    final_df['Exp_Status'] = final_df['sseqid'].apply(classify_expression)

    # 6. Build 2x2 Contingency Table
    ctable = pd.crosstab(final_df['Methylation'], final_df['Exp_Status'])
    
    # Ensure 2x2 structure
    for row in ['Methylated', 'Unmethylated']:
        if row not in ctable.index: 
            ctable.loc[row] = [0, 0]
    for col in ['UP', 'Background']:
        if col not in ctable.columns: 
            ctable[col] = 0
    
    ctable = ctable.loc[['Methylated', 'Unmethylated'], ['UP', 'Background']]

    print("\nContingency Table:")
    print(ctable)
    print(f"\nTotal genes analyzed: {ctable.sum().sum()}")

    # 7. Statistical Test
    if (ctable < 5).any().any():
        print("\n⚠ Small cell count detected. Using Fisher's Exact Test.")
        odds_ratio, p_val = fisher_exact(ctable)
    else:
        print("\nUsing Chi-Square Test of Independence.")
        chi2, p_val, dof, expected = chi2_contingency(ctable)
        
        # Calculate Odds Ratio: (a*d)/(b*c)
        a, b = ctable.iloc[0, 0], ctable.iloc[0, 1]  # Meth-UP, Meth-Background
        c, d = ctable.iloc[1, 0], ctable.iloc[1, 1]  # Unmeth-UP, Unmeth-Background
        odds_ratio = (a * d) / (b * c) if (b * c) != 0 else float('inf')
        print(f"  Chi-Square statistic: {chi2:.4f}")

    print(f"\n{'='*60}")
    print(f"RESULTS:")
    print(f"  P-value: {p_val:.4e}")
    print(f"  Odds Ratio: {odds_ratio:.4f}")
    
    if p_val < 0.05:
        direction = "enriched" if odds_ratio > 1 else "depleted"
        print(f"  ✓ SIGNIFICANT: Methylated genes are {direction} in up-regulated genes (p < 0.05)")
    else:
        print(f"  ✗ NOT SIGNIFICANT: No association between methylation and up-regulation")
    print(f"{'='*60}\n")
    
    return {
        'species': species,
        'n_genes': ctable.sum().sum(),
        'p_value': p_val,
        'odds_ratio': odds_ratio,
        'significant': p_val < 0.05,
        'contingency_table': ctable
    }

# --- EXECUTE ANALYSIS ---
if __name__ == "__main__":
    print("\n" + "="*60)
    print("UP-REGULATED GENE ANALYSIS")
    print("="*60)
    
    # Create output directory
    os.makedirs('Results/05b_directional_analysis/upregulated', exist_ok=True)
    
    results = []
    
    # Analyze S. brooksi
    results.append(run_upregulated_analysis(
        "brooksi", 
        "Results/05_ortholog_analysis/chacei_vs_brooksi.tsv", 
        "Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular"
    ))
    
    # Analyze S. elizabethae
    results.append(run_upregulated_analysis(
        "elizabethae", 
        "Results/05_ortholog_analysis/chacei_vs_elizabethae.tsv", 
        "Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular"
    ))
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for r in results:
        sig = "✓ SIGNIFICANT" if r['significant'] else "✗ NOT SIGNIFICANT"
        print(f"S. {r['species']:12s}: {sig:20s} (p={r['p_value']:.4e}, OR={r['odds_ratio']:.2f}, n={r['n_genes']})")
    print("="*60)
    
    # Save results
    print("\nSaving results...")
    
    # 1. Save summary statistics
    summary_df = pd.DataFrame([{
        'species': r['species'],
        'n_genes_analyzed': r['n_genes'],
        'p_value': r['p_value'],
        'odds_ratio': r['odds_ratio'],
        'significant': r['significant']
    } for r in results])
    summary_df.to_csv('Results/05b_directional_analysis/upregulated/upregulated_summary_statistics.csv', index=False)
    print("  ✓ Results/05b_directional_analysis/upregulated/upregulated_summary_statistics.csv")
    
    # 2. Save contingency tables
    for r in results:
        ctable = r['contingency_table']
        ctable.to_csv(f"Results/05b_directional_analysis/upregulated/{r['species']}_upregulated_contingency.csv")
        print(f"  ✓ Results/05b_directional_analysis/upregulated/{r['species']}_upregulated_contingency.csv")
    
    print("\n✓ Up-regulated analysis complete!")

