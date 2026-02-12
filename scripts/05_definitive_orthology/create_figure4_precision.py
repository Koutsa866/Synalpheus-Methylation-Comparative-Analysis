#!/usr/bin/env python3
"""
Figure 4 - Precision Version
Exact statistics with dotted line annotations
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# --- CONFIGURATION ---
sns.set_style("white") 
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300
DEG_COLOR = '#e74c3c'  # Social DEG
BG_COLOR = '#3498db'   # Background Blue

def create_publication_figure():
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, :])

    # 1. Exact Statistics from your previous successful run
    # Brooksi: 20/842 (2.4%), 133/4673 (2.8%)
    # Elizabethae: 84/783 (10.7%), 349/4229 (8.3%)
    
    species_list = [
        ('S. brooksi', [20, 822], [133, 4540]), 
        ('S. elizabethae', [84, 699], [349, 3880])
    ]

    for idx, (name, meth_counts, unmeth_counts) in enumerate(species_list):
        # Calculate Proportions
        m_total = sum(meth_counts)
        u_total = sum(unmeth_counts)
        m_deg_p = (meth_counts[0] / m_total) * 100
        u_deg_p = (unmeth_counts[0] / u_total) * 100
        
        off = idx * 3
        x_pos = [off, off + 1]
        heights = [m_deg_p, u_deg_p]
        
        # Plot Bars (Bottom Red, Top Blue)
        # We plot the full 100% background first, then overlay the DEG
        ax1.bar(x_pos, [100, 100], color=BG_COLOR, alpha=0.6, edgecolor='black', width=0.7)
        ax1.bar(x_pos, heights, color=DEG_COLOR, alpha=0.9, edgecolor='black', width=0.7)

        # 2. ADD PRECISION ANNOTATIONS (The Dotted Line + Out-of-the-way Text)
        for i, h in enumerate(heights):
            x = x_pos[i]
            # Thin dotted line extending from the bar toward the species center
            line_dir = 0.4 if i == 0 else -0.4
            ax1.plot([x, x + line_dir], [h, h], color='black', linestyle=':', linewidth=1, alpha=0.7)
            
            # Text placement: On the end of the dotted line, slightly above
            ax1.text(x + line_dir, h + 1, f'{h:.1f}%', ha='center', va='bottom', 
                     fontsize=10, fontweight='bold', color='black')
            
            # Non-DEG label (at the top 100% mark)
            ax1.text(x, 102, f'{(100-h):.1f}%', ha='center', fontsize=8, color='gray')

        # 3. Labeling and Stats
        ax1.text(off + 0.5, -15, rf"$\mathit{{{name}}}$", ha='center', fontsize=14, fontweight='bold')
        meth_rate = (m_total / (m_total + u_total)) * 100
        ax1.text(off + 0.5, -22, f"({meth_rate:.1f}% Methylated)", ha='center', fontsize=10, color='gray', style='italic')

    # Formatting Panel A
    ax1.set_ylim(0, 125)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_xticks([0, 1, 3, 4])
    ax1.set_xticklabels(['Methylated', 'Unmethylated', 'Methylated', 'Unmethylated'], fontweight='bold')
    ax1.set_ylabel('Percentage of Genes (%)', fontweight='bold', fontsize=12)
    ax1.set_title('A. Differential Gene Recruitment by Methylation Status', loc='left', fontsize=15, fontweight='bold')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=DEG_COLOR, edgecolor='black', label='Social DEG'),
                       Patch(facecolor=BG_COLOR, alpha=0.6, edgecolor='black', label='Non-DEG Background')]
    ax1.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

    sns.despine(ax=ax1, offset=10)
    plt.savefig('Results/05_ortholog_analysis/figure4_combined_summary.png', bbox_inches='tight', dpi=300)
    print("✓ Precision Figure 4 saved.")

if __name__ == "__main__":
    create_publication_figure()
