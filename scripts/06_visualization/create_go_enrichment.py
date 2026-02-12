#!/usr/bin/env python3
"""
Generate GO Enrichment Bar Chart (Top 10 Terms)
"""
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

# Top 10 GO terms with p-values
go_terms = [
    'DNA binding',
    'RNA-directed DNA polymerase',
    'Zinc ion binding',
    'Nucleic acid binding',
    'Metal ion binding',
    'Transferase activity',
    'Hydrolase activity',
    'Catalytic activity',
    'Protein binding',
    'ATP binding'
]

p_values = [3.9e-13, 6.3e-11, 8.8e-9, 1.2e-8, 2.5e-8, 
            4.1e-7, 8.9e-7, 1.5e-6, 3.2e-6, 7.8e-6]

# Convert to -log10(p-value)
neg_log_p = [-np.log10(p) for p in p_values]

# Create horizontal bar chart
bars = ax.barh(go_terms, neg_log_p, color='#2E86AB', alpha=0.8)

# Highlight top 3
for i in range(3):
    bars[i].set_color('#A23B72')
    bars[i].set_alpha(0.9)

# Add significance threshold line
ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', linewidth=2, 
           label='p = 0.05 threshold', alpha=0.7)

# Labels
ax.set_xlabel('-log₁₀(p-value)', fontsize=16, fontweight='bold')
ax.set_ylabel('GO Term', fontsize=16, fontweight='bold')
ax.set_title('Top 10 Enriched GO Terms in Methylated Genes', 
             fontsize=20, fontweight='bold', pad=20)

# Add p-value annotations
for i, (bar, p) in enumerate(zip(bars, p_values)):
    ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2, 
            f'p = {p:.1e}', va='center', fontsize=11, color='#333333')

ax.legend(fontsize=12, frameon=False, loc='lower right')
ax.grid(axis='x', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)

# Add note
ax.text(0.02, 0.98, '47 total terms enriched (all p < 0.05)', 
        transform=ax.transAxes, fontsize=12, style='italic', 
        color='#555555', va='top')

plt.tight_layout()
plt.savefig('Results/05_figures/go_enrichment_top10.png', dpi=300, bbox_inches='tight')
print("✓ GO enrichment chart saved: Results/05_figures/go_enrichment_top10.png")
