#!/usr/bin/env python3
"""
Directional Comparison Summary
Compares Any DEG, UP, and DOWN results to identify mechanistic patterns.

Biological Question:
  Is the methylation-DEG correlation driven by gene ACTIVATION or REPRESSION?

Method:
  - Loads results from all three analyses (Any DEG, UP, DOWN)
  - Generates comparison table
  - Identifies dominant mechanism per species

Inputs:
  - Results/05_ortholog_analysis/summary_statistics.csv (Any DEG)
  - Results/05b_directional_analysis/upregulated/upregulated_summary_statistics.csv
  - Results/05b_directional_analysis/downregulated/downregulated_summary_statistics.csv

Outputs:
  - Results/05b_directional_analysis/directional_comparison_summary.csv

Author: Philip Koutsaftis
Advisor: Dr. Solomon Chak
Institution: Denison University
Date: 2025
"""

import pandas as pd
import os

def load_results():
    """Load all three analysis results."""
    
    # Load Any DEG results
    any_deg = pd.read_csv('Results/05_ortholog_analysis/summary_statistics.csv')
    any_deg['analysis_type'] = 'Any_DEG'
    
    # Load Up-regulated results
    up_reg = pd.read_csv('Results/05b_directional_analysis/upregulated/upregulated_summary_statistics.csv')
    up_reg['analysis_type'] = 'Up_Regulated'
    
    # Load Down-regulated results
    down_reg = pd.read_csv('Results/05b_directional_analysis/downregulated/downregulated_summary_statistics.csv')
    down_reg['analysis_type'] = 'Down_Regulated'
    
    return any_deg, up_reg, down_reg

def create_comparison_table(any_deg, up_reg, down_reg):
    """Create side-by-side comparison table."""
    
    # Combine all results
    all_results = pd.concat([any_deg, up_reg, down_reg], ignore_index=True)
    
    # Pivot to wide format
    comparison = all_results.pivot(
        index='species',
        columns='analysis_type',
        values=['p_value', 'odds_ratio', 'significant']
    )
    
    # Flatten column names
    comparison.columns = ['_'.join(col).strip() for col in comparison.columns.values]
    comparison = comparison.reset_index()
    
    return comparison

def interpret_results(comparison):
    """Interpret the directional pattern for each species."""
    
    interpretations = []
    
    for _, row in comparison.iterrows():
        species = row['species']
        
        # Check significance
        any_sig = row['significant_Any_DEG']
        up_sig = row['significant_Up_Regulated']
        down_sig = row['significant_Down_Regulated']
        
        # Determine mechanism
        if any_sig and down_sig and not up_sig:
            mechanism = "Repression-driven (methylation predicts silencing)"
        elif any_sig and up_sig and not down_sig:
            mechanism = "Activation-driven (methylation predicts activation)"
        elif any_sig and up_sig and down_sig:
            mechanism = "Bidirectional (both activation and repression)"
        elif not any_sig:
            mechanism = "No correlation (lineage-specific loss)"
        else:
            mechanism = "Unclear pattern"
        
        interpretations.append({
            'species': species,
            'mechanism': mechanism,
            'any_deg_p': row['p_value_Any_DEG'],
            'up_p': row['p_value_Up_Regulated'],
            'down_p': row['p_value_Down_Regulated']
        })
    
    return pd.DataFrame(interpretations)

# --- EXECUTE ANALYSIS ---
if __name__ == "__main__":
    print("\n" + "="*60)
    print("DIRECTIONAL COMPARISON ANALYSIS")
    print("="*60)
    
    # Create output directory
    os.makedirs('Results/05b_directional_analysis', exist_ok=True)
    
    # Load results
    print("\nLoading results...")
    any_deg, up_reg, down_reg = load_results()
    print("  ✓ Any DEG results loaded")
    print("  ✓ Up-regulated results loaded")
    print("  ✓ Down-regulated results loaded")
    
    # Create comparison table
    print("\nGenerating comparison table...")
    comparison = create_comparison_table(any_deg, up_reg, down_reg)
    
    # Interpret results
    print("\nInterpreting directional patterns...")
    interpretation = interpret_results(comparison)
    
    # Display results
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)
    print("\nP-values by Analysis Type:")
    print(comparison[['species', 'p_value_Any_DEG', 'p_value_Up_Regulated', 'p_value_Down_Regulated']].to_string(index=False))
    
    print("\n" + "="*60)
    print("MECHANISTIC INTERPRETATION")
    print("="*60)
    for _, row in interpretation.iterrows():
        print(f"\nS. {row['species']}:")
        print(f"  {row['mechanism']}")
        print(f"  Any DEG: p={row['any_deg_p']:.4f}")
        print(f"  Up-regulated: p={row['up_p']:.4f}")
        print(f"  Down-regulated: p={row['down_p']:.4f}")
    
    # Save results
    print("\n" + "="*60)
    print("Saving results...")
    comparison.to_csv('Results/05b_directional_analysis/directional_comparison_summary.csv', index=False)
    print("  ✓ Results/05b_directional_analysis/directional_comparison_summary.csv")
    
    interpretation.to_csv('Results/05b_directional_analysis/mechanistic_interpretation.csv', index=False)
    print("  ✓ Results/05b_directional_analysis/mechanistic_interpretation.csv")
    
    print("\n✓ Comparison analysis complete!")
