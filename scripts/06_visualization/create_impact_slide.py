#!/usr/bin/env python3
"""
Generate Impact Slide: 77% TE Statistic
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(16, 9), dpi=300)

# Large percentage
ax.text(0.5, 0.65, '77%', ha='center', va='center', 
        fontsize=180, fontweight='bold', color='#2E86AB',
        transform=ax.transAxes)

# Subtitle
ax.text(0.5, 0.40, 'of annotated methylated genes', ha='center', va='center',
        fontsize=42, color='#333333', transform=ax.transAxes)
ax.text(0.5, 0.32, 'are transposable elements', ha='center', va='center',
        fontsize=42, fontweight='bold', color='#A23B72', transform=ax.transAxes)

# Statistical significance
ax.text(0.5, 0.18, 'p = 3.9 × 10⁻¹³', ha='center', va='center',
        fontsize=36, color='#555555', style='italic', transform=ax.transAxes)

# Add decorative elements
circle = mpatches.Circle((0.5, 0.65), 0.25, transform=ax.transAxes,
                        fill=False, edgecolor='#2E86AB', linewidth=8, alpha=0.3)
ax.add_patch(circle)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout(pad=0)
plt.savefig('Results/05_figures/te_percentage_impact.png', dpi=300, 
            bbox_inches='tight', facecolor='white')
print("✓ Impact slide saved: Results/05_figures/te_percentage_impact.png")
