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
            # Placed directly above the red DEG segment
            ax1.text(x_coord, h + 2, f'{h:.1f}%', 
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
    
    # --- PANEL B: FOREST PLOT ---
    ax2 = fig.add_subplot(gs[1, 0])
    
    # Calculate OR and CI for each species
    def calc_or_ci(meth, unmeth):
        a, b = meth[0], meth[1]  # DEG, Non-DEG in methylated
        c, d = unmeth[0], unmeth[1]  # DEG, Non-DEG in unmethylated
        or_val = (a * d) / (b * c) if (b * c) != 0 else 0
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        log_or = np.log(or_val)
        ci_low = np.exp(log_or - 1.96 * se)
        ci_high = np.exp(log_or + 1.96 * se)
        return or_val, ci_low, ci_high
    
    ors = []
    cis_low = []
    cis_high = []
    p_values = [0.507, 0.028]  # From your statistical results
    
    for meth, unmeth in [(species_list[0][1], species_list[0][2]), (species_list[1][1], species_list[1][2])]:
        o, l, h = calc_or_ci(meth, unmeth)
        ors.append(o)
        cis_low.append(l)
        cis_high.append(h)
    
    y_pos = [1, 0]
    colors = ['#95a5a6', '#27ae60']  # Gray for brooksi, green for elizabethae
    
    for i in range(2):
        ax2.plot([cis_low[i], cis_high[i]], [y_pos[i], y_pos[i]], color=colors[i], lw=3, zorder=1)
        ax2.scatter(ors[i], y_pos[i], s=200, color=colors[i], edgecolor='black', lw=2, zorder=2)
        if p_values[i] < 0.05:
            ax2.text(cis_high[i] + 0.05, y_pos[i], '*', fontsize=20, fontweight='bold', va='center', color='#27ae60')
    
    ax2.axvline(x=1, color='red', linestyle='--', lw=2)
    ax2.text(0.5, 1.3, r'$\leftarrow$ Depleted', color='gray', ha='center', fontweight='bold', fontsize=9)
    ax2.text(1.5, 1.3, r'Enriched $\rightarrow$', color='#27ae60', ha='center', fontweight='bold', fontsize=9)
    
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([r'$\mathit{S.\ brooksi}$', r'$\mathit{S.\ elizabethae}$'], fontsize=12)
    ax2.set_xlabel('Odds Ratio (95% CI)', fontweight='bold')
    ax2.set_xlim(0, 2.0)
    ax2.set_title('B. Epigenetic Effect Size', loc='left', fontsize=15, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    # --- PANEL C: SUMMARY TABLE ---
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    table_data = []
    for i, (name, meth, unmeth) in enumerate(species_list):
        n_genes = sum(meth) + sum(unmeth)
        meth_pct = (sum(meth) / n_genes) * 100
        sig = "Significant" if p_values[i] < 0.05 else "Not Sig."
        table_data.append([name, n_genes, f"{meth_pct:.1f}%", f"{ors[i]:.2f}", f"{p_values[i]:.3f}", sig])
    
    table = ax3.table(cellText=table_data, 
                      colLabels=['Species', 'N Genes', '% Meth', 'OR', 'P-value', 'Result'],
                      loc='center', cellLoc='center', bbox=[0, 0.2, 1, 0.6])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    
    # Style header
    for i in range(6):
        table[(0, i)].set_facecolor('#34495e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style rows
    for i in range(1, 3):
        table[(i, 0)].get_text().set_fontstyle('italic')
        res_cell = table[(i, 5)]
        res_cell.set_facecolor('#d5f4e6' if "Significant" in res_cell.get_text().get_text() else '#f2f4f4')
    
    ax3.set_title('C. Statistical Summary', loc='left', fontsize=15, fontweight='bold')
    
    # Save output
    plt.savefig('Results/05_ortholog_analysis/figure4_combined_summary.png', bbox_inches='tight', dpi=300)
    print("✓ Figure 4 (Complete 3-Panel) saved.")

if __name__ == "__main__":
    create_publication_figure()
