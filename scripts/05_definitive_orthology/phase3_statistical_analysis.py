#!/usr/bin/env python3
"""
Pipeline: 05_definitive_orthology
Phase: 3 - Statistical Integration
Goal: Integrate methylation and DEG status, perform Chi-square tests on 2×2 contingency tables
Note: Implements Dr. Chak's requirement for whole-transcriptome background normalization
Inputs:
  - Results/05_ortholog_analysis/brooksi_diamond_results.tsv
  - Results/05_ortholog_analysis/elizabethae_diamond_results.tsv
  - Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular
  - Data/transcriptomics/deseq2_results/S.elizabethae_brain_DEA_wholetable.tabular
Outputs:
  - Results/05_ortholog_analysis/brooksi_contingency_table.csv
  - Results/05_ortholog_analysis/elizabethae_contingency_table.csv
  - Results/05_ortholog_analysis/summary_statistics.csv
Author: Philip Koutsaftis
Date: 2025
"""

"""
Phase 3: Statistical Integration and Analysis
Integrates methylation status with DEG status to test correlation
"""

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

def run_integration(species, ortho_path, dea_path):
    print(f"\n{'='*60}")
    print(f"ANALYSIS: S. {species.upper()}")
    print(f"{'='*60}")

    # 1. Load DEA Table (Transcript IDs are row indices, padj is column 9)
    dea_df = pd.read_csv(dea_path, sep='\t', index_col=0, low_memory=False)
    
    padj_col = dea_df.columns[8]  # Column 9 becomes index 8 after index_col=0

    # Convert padj to numeric, handle NaN (low-count genes)
    dea_df[padj_col] = pd.to_numeric(dea_df[padj_col], errors='coerce')
    dea_df[padj_col] = dea_df[padj_col].fillna(1.0)  # Treat NaN as non-significant

    print(f"Loaded DEA data: {len(dea_df)} transcripts")
    print(f"  DEGs (padj < 0.05): {(dea_df[padj_col] < 0.05).sum()}")

    # 2. Load DIAMOND Hits
    ortho_cols = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
    ortho_df = pd.read_csv(ortho_path, sep='\t', names=ortho_cols)
    print(f"Loaded DIAMOND hits: {len(ortho_df)} orthologs")

    # 3. Filter for "Successful Mappers Only" (PI requirement)
    mapped_ids = set(dea_df.index.tolist())  # Use index, not column
    final_df = ortho_df[ortho_df['sseqid'].isin(mapped_ids)].copy()
    print(f"Successful mappers: {len(final_df)} genes")

    # 4. Extract Methylation Status from FASTA headers (gene_id|Status)
    final_df['Methylation'] = final_df['qseqid'].apply(lambda x: x.split('|')[1])

    # 5. Map DEG Status from DEA table
    deg_status_map = dict(zip(dea_df.index, dea_df[padj_col]))  # Use index
    final_df['Exp_Status'] = final_df['sseqid'].apply(
        lambda x: 'DEG' if deg_status_map.get(x, 1.0) < 0.05 else 'Non-DEG'
    )

    # 6. Build 2x2 Contingency Table
    ctable = pd.crosstab(final_df['Methylation'], final_df['Exp_Status'])
    
    # Ensure 2x2 structure (handle missing categories)
    for row in ['Methylated', 'Unmethylated']:
        if row not in ctable.index: 
            ctable.loc[row] = [0, 0]
    for col in ['DEG', 'Non-DEG']:
        if col not in ctable.columns: 
            ctable[col] = 0
    
    ctable = ctable.loc[['Methylated', 'Unmethylated'], ['DEG', 'Non-DEG']]

    print("\nContingency Table:")
    print(ctable)
    print(f"\nTotal genes analyzed: {ctable.sum().sum()}")

    # 7. Statistical Test
    # Use Fisher's Exact if any cell < 5, otherwise Chi-Square
    if (ctable < 5).any().any():
        print("\n⚠ Small cell count detected. Using Fisher's Exact Test.")
        odds_ratio, p_val = fisher_exact(ctable)
    else:
        print("\nUsing Chi-Square Test of Independence.")
        chi2, p_val, dof, expected = chi2_contingency(ctable)
        
        # Calculate Odds Ratio: (a*d)/(b*c)
        a, b = ctable.iloc[0, 0], ctable.iloc[0, 1]  # Meth-DEG, Meth-NonDEG
        c, d = ctable.iloc[1, 0], ctable.iloc[1, 1]  # Unmeth-DEG, Unmeth-NonDEG
        odds_ratio = (a * d) / (b * c) if (b * c) != 0 else float('inf')
        print(f"  Chi-Square statistic: {chi2:.4f}")

    print(f"\n{'='*60}")
    print(f"RESULTS:")
    print(f"  P-value: {p_val:.4e}")
    print(f"  Odds Ratio: {odds_ratio:.4f}")
    
    if p_val < 0.05:
        direction = "enriched" if odds_ratio > 1 else "depleted"
        print(f"  ✓ SIGNIFICANT: Methylated genes are {direction} in DEGs (p < 0.05)")
    else:
        print(f"  ✗ NOT SIGNIFICANT: No association between methylation and DEG status")
    print(f"{'='*60}\n")
    
    # Return results for summary
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
    print("PHASE 3: METHYLATION-DEG CORRELATION ANALYSIS")
    print("="*60)
    
    results = []
    
    # Analyze S. brooksi
    results.append(run_integration(
        "brooksi", 
        "Results/05_ortholog_analysis/chacei_vs_brooksi.tsv", 
        "Data/transcriptomics/deseq2_results/S.brooksi_brain_DEA_wholetable.tabular"
    ))
    
    # Analyze S. elizabethae
    results.append(run_integration(
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
    
    # Save results to files
    print("\nSaving results...")
    
    # 1. Save summary statistics
    summary_df = pd.DataFrame([{
        'species': r['species'],
        'n_genes_analyzed': r['n_genes'],
        'p_value': r['p_value'],
        'odds_ratio': r['odds_ratio'],
        'significant': r['significant']
    } for r in results])
    summary_df.to_csv('Results/05_ortholog_analysis/summary_statistics.csv', index=False)
    print("  ✓ Results/05_ortholog_analysis/summary_statistics.csv")
    
    # 2. Save contingency tables
    for r in results:
        ctable = r['contingency_table']
        ctable.to_csv(f"Results/05_ortholog_analysis/{r['species']}_contingency_table.csv")
        print(f"  ✓ Results/05_ortholog_analysis/{r['species']}_contingency_table.csv")
    
    print("\n✓ All results saved to Results/05_ortholog_analysis/")
