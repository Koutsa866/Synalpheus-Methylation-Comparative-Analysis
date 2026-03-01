#!/bin/bash
"""
Master Script: Complete Directional DEG Analysis
Runs all four analysis scripts in sequence.

Usage:
  bash run_complete_directional_analysis.sh

Workflow:
  1. analyze_upregulated.py - Test for activation signal
  2. analyze_downregulated.py - Test for repression signal
  3. compare_directions.py - Identify dominant mechanism
  4. extract_gene_lists.py - Export gene IDs for annotation

Author: Philip Koutsaftis
Advisor: Dr. Solomon Chak
Institution: Denison University
Date: 2025
"""

set -e  # Exit on error

echo "========================================"
echo "DIRECTIONAL DEG ANALYSIS PIPELINE"
echo "========================================"
echo ""

# Activate conda environment
echo "Activating shrimp_local environment..."
eval "$(conda shell.bash hook)"
conda activate shrimp_local

# Step 1: Up-regulated analysis
echo ""
echo "Step 1/4: Running up-regulated gene analysis..."
python scripts/05b_directional_deg_analysis/analyze_upregulated.py

# Step 2: Down-regulated analysis
echo ""
echo "Step 2/4: Running down-regulated gene analysis..."
python scripts/05b_directional_deg_analysis/analyze_downregulated.py

# Step 3: Comparison
echo ""
echo "Step 3/4: Comparing directional results..."
python scripts/05b_directional_deg_analysis/compare_directions.py

# Step 4: Gene list extraction
echo ""
echo "Step 4/4: Extracting gene lists..."
python scripts/05b_directional_deg_analysis/extract_gene_lists.py

# Summary
echo ""
echo "========================================"
echo "PIPELINE COMPLETE"
echo "========================================"
echo ""
echo "Results saved to:"
echo "  - Results/05b_directional_analysis/upregulated/"
echo "  - Results/05b_directional_analysis/downregulated/"
echo "  - Results/05b_directional_analysis/gene_lists/"
echo "  - Results/05b_directional_analysis/directional_comparison_summary.csv"
echo ""
echo "✓ All analyses complete!"
