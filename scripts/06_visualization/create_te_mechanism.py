#!/usr/bin/env python3
"""
Generate TE Silencing Conceptual Diagram
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Wedge

fig, ax = plt.subplots(figsize=(16, 9), dpi=300)

# Left side: Active TE (unmethylated)
ax.text(3.5, 8, 'Unmethylated TE', ha='center', fontsize=24, fontweight='bold', color='#A23B72')

# DNA strand (unmethylated)
rect1 = FancyBboxPatch((2, 6), 3, 0.4, boxstyle="round,pad=0.05",
                       edgecolor='#333', facecolor='#E8E8E8', linewidth=2)
ax.add_patch(rect1)
ax.text(3.5, 6.2, 'TE Promoter', ha='center', fontsize=14, color='#333')

# Transcription arrow
arrow1 = FancyArrowPatch((3.5, 5.8), (3.5, 4.5), arrowstyle='->', 
                        mutation_scale=40, linewidth=4, color='#A23B72')
ax.add_patch(arrow1)

# mRNA
rect2 = FancyBboxPatch((2.5, 3.8), 2, 0.3, boxstyle="round,pad=0.05",
                       edgecolor='#A23B72', facecolor='#FFE6F0', linewidth=2)
ax.add_patch(rect2)
ax.text(3.5, 4, 'mRNA', ha='center', fontsize=14, fontweight='bold', color='#A23B72')

# Protein
circle1 = mpatches.Circle((3.5, 2.5), 0.5, edgecolor='#A23B72', 
                         facecolor='#FFB3D9', linewidth=2)
ax.add_patch(circle1)
ax.text(3.5, 2.5, 'TE\nProtein', ha='center', va='center', 
        fontsize=12, fontweight='bold', color='#A23B72')

# Genome instability
ax.text(3.5, 1.2, '⚠ Genome Instability', ha='center', fontsize=16, 
        fontweight='bold', color='#D32F2F')

# Right side: Silenced TE (methylated)
ax.text(12.5, 8, 'Methylated TE', ha='center', fontsize=24, fontweight='bold', color='#2E86AB')

# DNA strand with methylation marks
rect3 = FancyBboxPatch((11, 6), 3, 0.4, boxstyle="round,pad=0.05",
                       edgecolor='#333', facecolor='#E8E8E8', linewidth=2)
ax.add_patch(rect3)
ax.text(12.5, 6.2, 'TE Promoter', ha='center', fontsize=14, color='#333')

# Methylation marks (CH3)
for x_pos in [11.5, 12, 12.5, 13, 13.5]:
    ax.text(x_pos, 6.6, 'CH₃', ha='center', fontsize=12, 
            fontweight='bold', color='#2E86AB')

# Blocked transcription
ax.text(12.5, 5, '✖', ha='center', fontsize=60, color='#D32F2F', alpha=0.7)
ax.text(12.5, 4.2, 'Transcription\nBlocked', ha='center', fontsize=14, 
        fontweight='bold', color='#D32F2F')

# Checkmark for stability
ax.text(12.5, 2.5, '✓', ha='center', fontsize=80, color='#4CAF50', alpha=0.7)
ax.text(12.5, 1.2, 'Genome Stability', ha='center', fontsize=16, 
        fontweight='bold', color='#4CAF50')

# Central arrow showing the process
arrow_main = FancyArrowPatch((6, 6.2), (10, 6.2), arrowstyle='->', 
                            mutation_scale=50, linewidth=5, color='#2E86AB')
ax.add_patch(arrow_main)
ax.text(8, 6.8, 'DNA Methylation', ha='center', fontsize=18, 
        fontweight='bold', color='#2E86AB')

# Title
ax.text(8, 9.5, 'Epigenetic Genome Defense Mechanism', ha='center', 
        fontsize=28, fontweight='bold', color='#1a1a1a')

ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis('off')

plt.tight_layout(pad=0)
plt.savefig('Results/05_figures/te_silencing_mechanism.png', dpi=300, 
            bbox_inches='tight', facecolor='white')
print("✓ TE silencing diagram saved: Results/05_figures/te_silencing_mechanism.png")
