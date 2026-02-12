#!/usr/bin/env python3
"""
Create GO enrichment chart with bold/thick text for poster readability
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the GO enrichment data
df = pd.read_csv('Results/04_functional_annotation/go_enrichment_results.csv')

# Get top 15 terms
top_terms = df.head(15).copy()

# Calculate -log10(p-value)
top_terms['neg_log_p'] = -np.log10(top_terms['p_value'])

# Create the plot with larger, bolder text
plt.rcParams.update({
    'font.size': 14,
    'font.weight': 'bold',  # Make all text bold
    'font.family': 'sans-serif'
})

fig, ax = plt.subplots(figsize=(12, 10))

# Create horizontal bar chart
bars = ax.barh(range(len(top_terms)), top_terms['neg_log_p'], 
               color='steelblue', alpha=0.8, height=0.7)

# Set y-axis labels with bold, larger font
ax.set_yticks(range(len(top_terms)))
ax.set_yticklabels(top_terms['go_term'], fontsize=16, fontweight='bold')

# Set x-axis
ax.set_xlabel('-log10(p-value)', fontsize=16, fontweight='bold')
ax.set_xlim(0, max(top_terms['neg_log_p']) * 1.1)

# Add significance line at p=0.05
sig_line = -np.log10(0.05)
ax.axvline(x=sig_line, color='red', linestyle='--', linewidth=2, alpha=0.8)

# Add legend box in bottom-right corner (like original)
ax.text(0.95, 0.05, 'p=0.05', transform=ax.transAxes,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red', alpha=0.8),
        color='red', fontsize=12, fontweight='bold', ha='right', va='bottom')

# Title with bold text
plt.title('GO Terms Enriched in Methylated Gene Promoters\n(Synalpheus chacei)', 
          fontsize=18, fontweight='bold', pad=20)

# Make axis labels bold
ax.tick_params(axis='x', labelsize=12, width=2)
ax.tick_params(axis='y', labelsize=13, width=2)

# Thicker spines
for spine in ax.spines.values():
    spine.set_linewidth(2)

# Invert y-axis so highest significance is at top
ax.invert_yaxis()

# Tight layout
plt.tight_layout()

# Save as high-resolution PNG
plt.savefig('Results/05_figures/go_enrichment_chart_bold.png', 
            dpi=300, bbox_inches='tight', facecolor='white')

print("✅ Bold GO enrichment chart saved as: Results/05_figures/go_enrichment_chart_bold.png")
plt.show()