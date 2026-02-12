#!/usr/bin/env python3
"""
Create 77% TE Impact Slide with Pie Chart
Combines large percentage display with visual pie chart
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Data
te_count = 24
functional_count = 7
total = 31
te_percentage = 77.4

# Create figure with two subplots side by side
fig = plt.figure(figsize=(16, 9), facecolor='white', dpi=300)

# Left side: Giant 77%
ax1 = plt.subplot(1, 2, 1)
ax1.text(0.5, 0.5, '77%', 
         ha='center', va='center',
         fontsize=180, fontweight='bold',
         color='#A23B72',
         transform=ax1.transAxes)
ax1.text(0.5, 0.25, 'of methylated genes\nare transposable elements',
         ha='center', va='center',
         fontsize=28, fontweight='bold',
         color='#1a1a1a',
         transform=ax1.transAxes)
ax1.axis('off')

# Right side: Pie chart
ax2 = plt.subplot(1, 2, 2)
colors = ['#A23B72', '#2E86AB']
explode = (0.1, 0)

wedges, texts, autotexts = ax2.pie(
    [te_count, functional_count],
    explode=explode,
    labels=['Transposable\nElements', 'Other\nGenes'],
    colors=colors,
    autopct='%1.0f%%',
    startangle=90,
    textprops={'fontsize': 20, 'fontweight': 'bold'},
    wedgeprops={'edgecolor': '#1a1a1a', 'linewidth': 3}
)

# Make percentage text white and larger
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontsize(32)
    autotext.set_fontweight('bold')

# Add center annotation
ax2.text(0, 0, f'n = {total}', 
         ha='center', va='center',
         fontsize=18, fontweight='bold',
         color='#1a1a1a')

# Add p-value annotation
ax2.text(0, -1.5, 'p = 3.9 × 10⁻¹³',
         ha='center', va='center',
         fontsize=20, fontweight='bold',
         color='#1a1a1a',
         bbox=dict(boxstyle='round,pad=0.5', 
                  facecolor='#FFE6F0', 
                  edgecolor='#A23B72', 
                  linewidth=2))

plt.tight_layout()
plt.savefig('Results/05_figures/te_77_impact_with_pie.png', 
            dpi=300, bbox_inches='tight', facecolor='white')
print("✓ Created: Results/05_figures/te_77_impact_with_pie.png")
print("  Giant 77% on left, pie chart on right")
