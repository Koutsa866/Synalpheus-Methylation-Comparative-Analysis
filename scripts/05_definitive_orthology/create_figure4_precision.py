#!/usr/bin/env python3
"""
Figure 4 - Precision Schematic Version
Exact statistics with L-shaped drop lines to x-axis
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch

# --- CONFIGURATION ---
sns.set_style("white") 
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300
DEG_COLOR = '#e74c3c'  # Social DEG
BG_COLOR = '#3498db'   # Background Blue

def create_publication_figure():
    # Use a large layout to ensure no overlap
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, :])

    # 1. Exact Statistics
    # Brooksi: 20/842 (2.4%), 133/4673 (2.8%)
    # Elizabethae: 84/783 (10.7%), 349/4229 (8.3%)
    species_list = [
        ('S. brooksi', [20, 822], [133, 4540]), 
        ('S. elizabethae', [84, 699], [349, 3880])
    ]

    for idx, (name, meth_counts, unmeth_counts) in enumerate(species_list):
        # Calculate Proportions accurately
        m_total = sum(meth_counts)
        u_total = sum(unmeth_counts)
        m_deg_p = (meth_counts[0] / m_total) * 100
        u_deg_p = (unmeth_counts[0] / u_total) * 100
        
        off = idx * 3
        x_pos = [off, off + 1]
        heights = [m_deg_p, u_deg_p]
        
        # Plot 100% Background Bars (Blue)
        ax1.bar(x_pos, [100, 100], color=BG_COLOR, alpha=0.6, edgecolor='black', width=0.7, zorder=1)
        # Plot DEG Bars (Red)
        ax1.bar(x_pos, heights, color=DEG_COLOR, alpha=0.9, edgecolor='black', width=0.7, zorder=2)

        # 2. PRECISION CONNECTORS (L-shaped drop line to X-axis)
        for i, h in enumerate(heights):
            x_coord = x_pos[i]
            
            # Direction of the callout (pointing towards the center of each species group)
            callout_dir = 0.45 if i == 0 else -0.45
            callout_x = x_coord + callout_dir

            # Draw the "Drop Line": 
            # (Bar transition -> Side extension -> Drop to X-axis)
            line_coords_x = [x_coord, callout_x, callout_x]
            line_coords_y = [h, h, 0]
            
            ax1.plot(line_coords_x, line_coords_y, color='black', 
                     linestyle=':', linewidth=0.8, alpha=0.6, zorder=3)
            
            # 3. PERCENTAGE TEXT (Small, Black, Not Bold)
            # Placed slightly above the horizontal segment of the dotted line
            ax1.text(callout_x, h + 0.5, f'{h:.1f}%', 
                     ha='center', va='bottom', fontsize=9, 
                     fontweight='normal', color='black')
            
            # Non-DEG label at the top 100% mark (Subtle)
            ax1.text(x_coord, 101.5, f'{(100-h):.1f}%', ha='center', fontsize=8, color='#7f8c8d')

        # 4. Species and Global Methylation Status Labels
        ax1.text(off + 0.5, -15, rf"$\mathit{{{name}}}$", ha='center', fontsize=14, fontweight='bold')
        meth_total_rate = (m_total / (m_total + u_total)) * 100
        ax1.text(off + 0.5, -21, f"({meth_total_rate:.1f}% Methylated)", 
                 ha='center', fontsize=10, color='#555555', style='italic')

    # --- FORMATTING PANEL A ---
    ax1.set_ylim(0, 125)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_xticks([0, 1, 3, 4])
    ax1.set_xticklabels(['Methylated', 'Unmethylated', 'Methylated', 'Unmethylated'], fontweight='bold')
    ax1.set_ylabel('Percentage of Genes (%)', fontweight='bold', fontsize=12)
    ax1.set_title('A. Differential Gene Recruitment by Methylation Status', loc='left', fontsize=15, fontweight='bold')
    
    # Clean Legend
    legend_elements = [Patch(facecolor=DEG_COLOR, edgecolor='black', label='Social DEG'),
                       Patch(facecolor=BG_COLOR, alpha=0.6, edgecolor='black', label='Non-DEG Background')]
    ax1.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)

    sns.despine(ax=ax1, offset=10)
    
    # Save output
    plt.savefig('Results/05_ortholog_analysis/figure4_combined_summary.png', bbox_inches='tight', dpi=300)
    print("✓ Figure 4 (Precision Schematic) saved.")

if __name__ == "__main__":
    create_publication_figure()
