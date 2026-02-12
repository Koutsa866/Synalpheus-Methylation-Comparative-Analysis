#!/usr/bin/env python3
"""
Generate Methylation Landscape Bar Chart
"""
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

categories = ['Methylated\nPromoters', 'Unmethylated\nPromoters']
gene_counts = [1103, 5745]
annotation_rates = [2.8, 94.0]

x = np.arange(len(categories))
width = 0.35

# Gene counts
ax1 = ax
bars1 = ax1.bar(x - width/2, gene_counts, width, label='Gene Count', 
                color='#2E86AB', alpha=0.8)

# Annotation rates on secondary axis
ax2 = ax1.twinx()
bars2 = ax2.bar(x + width/2, annotation_rates, width, label='Annotation Rate (%)', 
                color='#A23B72', alpha=0.8)

# Labels
ax1.set_ylabel('Number of Genes', fontsize=16, fontweight='bold', color='#2E86AB')
ax2.set_ylabel('Annotation Success Rate (%)', fontsize=16, fontweight='bold', color='#A23B72')
ax1.set_xlabel('Promoter Methylation Status', fontsize=16, fontweight='bold')
ax1.set_title('Methylation Landscape in S. chacei Genome', fontsize=20, fontweight='bold', pad=20)
ax1.set_xticks(x)
ax1.set_xticklabels(categories, fontsize=14)

# Add value labels
for i, (bar, val) in enumerate(zip(bars1, gene_counts)):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 100, 
             f'{val:,}', ha='center', va='bottom', fontsize=14, fontweight='bold', color='#2E86AB')

for i, (bar, val) in enumerate(zip(bars2, annotation_rates)):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, 
             f'{val}%', ha='center', va='bottom', fontsize=14, fontweight='bold', color='#A23B72')

# Legends
ax1.legend(loc='upper left', fontsize=12, frameon=False)
ax2.legend(loc='upper right', fontsize=12, frameon=False)

# Grid
ax1.yaxis.grid(True, alpha=0.3, linestyle='--')
ax1.set_axisbelow(True)

# Total annotation
ax1.text(0.5, 0.02, 'Total: 6,848 genes | 80% annotated overall', 
         transform=ax.transAxes, ha='center', fontsize=12, 
         style='italic', color='#555555')

plt.tight_layout()
plt.savefig('Results/05_figures/methylation_landscape.png', dpi=300, bbox_inches='tight')
print("✓ Methylation landscape saved: Results/05_figures/methylation_landscape.png")
