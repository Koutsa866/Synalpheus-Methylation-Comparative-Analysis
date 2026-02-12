#!/usr/bin/env python3
"""
Create TE enrichment figure for poster
Shows 77% TE enrichment in annotated methylated genes
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Data: 24 TEs out of 31 annotated methylated genes
te_count = 24
functional_count = 7
total = 31

# Calculate percentages
te_pct = (te_count / total) * 100
func_pct = (functional_count / total) * 100

# Create figure with poster-appropriate styling
fig, ax = plt.subplots(figsize=(8, 6), facecolor='white')

# Create stacked horizontal bar
categories = ['Annotated\nMethylated Genes\n(n=31)']
te_bar = ax.barh(categories, te_count, color='#d62728', edgecolor='#206082', linewidth=2)
func_bar = ax.barh(categories, functional_count, left=te_count, color='#4472C4', 
                   edgecolor='#206082', linewidth=2)

# Add value labels inside bars
ax.text(te_count/2, 0, f'TEs\n24 (77%)', ha='center', va='center', 
        fontsize=16, fontweight='bold', color='white')
ax.text(te_count + functional_count/2, 0, f'Functional\n7 (23%)', ha='center', va='center',
        fontsize=14, fontweight='bold', color='white')

# Styling
ax.set_xlabel('Number of Genes', fontsize=14, fontweight='bold')
ax.set_xlim(0, 35)
ax.set_ylim(-0.5, 0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params(axis='both', labelsize=12, width=2)

# Title
ax.set_title('TE Enrichment in Methylated Gene Promoters', 
             fontsize=16, fontweight='bold', pad=20)

# Add statistical annotation
ax.text(0.98, 0.95, 'p = 3.9 × 10⁻¹³', transform=ax.transAxes,
        fontsize=12, fontweight='bold', ha='right', va='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('Results/05_figures/te_enrichment_bar.png', dpi=300, bbox_inches='tight')
print("✓ Created: Results/05_figures/te_enrichment_bar.png")

# Also create a pie chart version
fig2, ax2 = plt.subplots(figsize=(7, 7), facecolor='white')

colors = ['#d62728', '#4472C4']
explode = (0.05, 0)
wedges, texts, autotexts = ax2.pie([te_count, functional_count], 
                                     explode=explode,
                                     labels=['Transposable\nElements', 'Functional\nGenes'],
                                     colors=colors,
                                     autopct='%1.0f%%',
                                     startangle=90,
                                     textprops={'fontsize': 14, 'fontweight': 'bold'},
                                     wedgeprops={'edgecolor': '#206082', 'linewidth': 2})

# Make percentage text white and larger
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontsize(18)
    autotext.set_fontweight('bold')

# Add center text
ax2.text(0, 0, f'n = {total}\ngenes', ha='center', va='center',
         fontsize=14, fontweight='bold')

ax2.set_title('TE Composition of Annotated Methylated Genes\n(p = 3.9 × 10⁻¹³)', 
              fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('Results/05_figures/te_enrichment_pie.png', dpi=300, bbox_inches='tight')
print("✓ Created: Results/05_figures/te_enrichment_pie.png")

print("\nBoth figures created successfully!")
print("- Bar chart: te_enrichment_bar.png (recommended for poster)")
print("- Pie chart: te_enrichment_pie.png (alternative)")