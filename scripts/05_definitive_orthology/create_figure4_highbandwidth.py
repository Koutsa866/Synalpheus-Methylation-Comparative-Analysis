#!/usr/bin/env python3
"""
High-Bandwidth Figure 4 - Publication Ready
Adds directional anchors, moves legend, italicizes species names
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
import seaborn as sns

# --- STYLE & CONFIGURATION ---
sns.set_style("white")  # Cleaner background to make red pop
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 300
SIG_COLOR = '#27ae60'  # Professional Green
NS_COLOR = '#95a5a6'   # Neutral Gray
DEG_COLOR = '#e74c3c'  # Social DEG (High Contrast Red)
BG_COLOR = '#3498db'   # Background

# Load Data
brooksi = pd.read_csv('Results/05_ortholog_analysis/brooksi_contingency_table.csv', index_col=0)
elizabethae = pd.read_csv('Results/05_ortholog_analysis/elizabethae_contingency_table.csv', index_col=0)
summary = pd.read_csv('Results/05_ortholog_analysis/summary_statistics.csv')

def get_stats(ctable):
    a, b = ctable.iloc[0, 0], ctable.iloc[0, 1]
    c, d = ctable.iloc[1, 0], ctable.iloc[1, 1]
    or_val = (a * d) / (b * c) if (b * c) != 0 else 0
    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    low = np.exp(np.log(or_val) - 1.96 * se)
    high = np.exp(np.log(or_val) + 1.96 * se)
    return or_val, low, high

# Build Figure
fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)

# Calculate Global Methylation Rates
brooksi_meth_rate = (brooksi.loc['Methylated'].sum() / brooksi.sum().sum()) * 100
eliz_meth_rate = (elizabethae.loc['Methylated'].sum() / elizabethae.sum().sum()) * 100
meth_rates = [brooksi_meth_rate, eliz_meth_rate]

# --- PANEL A: ENRICHMENT BARS ---
ax1 = fig.add_subplot(gs[0, :])
species_data = [('S. brooksi', brooksi), ('S. elizabethae', elizabethae)]

for idx, (name, ct) in enumerate(species_data):
    props = ct.div(ct.sum(axis=1), axis=0) * 100
    off = idx * 3
    x_pos = np.array([0, 1]) + off
    
    # Plot Bars
    bars_deg = ax1.bar(x_pos, props['DEG'], color=DEG_COLOR, alpha=0.9, edgecolor='black', linewidth=1.2, label='Social DEG' if idx==0 else "")
    bars_bg = ax1.bar(x_pos, props['Non-DEG'], bottom=props['DEG'], color=BG_COLOR, alpha=0.6, edgecolor='black', linewidth=1, label='Non-DEG Background' if idx==0 else "")
    
    # ADD PERCENTAGE LABELS (High-Contrast Upgrade)
    for i, bar in enumerate(bars_deg):
        val = props['DEG'].iloc[i]
        # Bold label for DEG portion
        ax1.text(bar.get_x() + bar.get_width()/2, val + 1, f'{val:.1f}%', 
                 ha='center', va='bottom', fontweight='bold', fontsize=12, color=DEG_COLOR)
        
        # Subdued label for Non-DEG portion
        bg_val = props['Non-DEG'].iloc[i]
        ax1.text(bar.get_x() + bar.get_width()/2, 98, f'{bg_val:.1f}%', 
                 ha='center', va='top', fontsize=9, color='white', alpha=0.9)
    
    # Species and Global Methylation Labels
    ax1.text(off + 0.5, -12, rf"$\mathit{{{name}}}$", ha='center', fontweight='bold', fontsize=13)
    ax1.text(off + 0.5, -18, f"({meth_rates[idx]:.1f}% Methylated)", ha='center', fontsize=10, color='gray', style='italic')
    
    # Stat Boxes
    p = summary.iloc[idx]['p_value']
    color = SIG_COLOR if p < 0.05 else NS_COLOR
    ax1.text(off + 0.5, 112, f"OR={summary.iloc[idx]['odds_ratio']:.2f}\np={p:.4f}", 
             ha='center', va='center', fontweight='bold', color=color,
             bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=color, lw=2))

ax1.set_ylim(-25, 135)
ax1.set_yticks([0, 20, 40, 60, 80, 100])
ax1.set_xticks([0, 1, 3, 4])
ax1.set_xticklabels(['Methylated', 'Unmethylated', 'Methylated', 'Unmethylated'], fontweight='bold')
ax1.set_ylabel('Percentage of Genes (%)', fontweight='bold')
ax1.set_title('A. Differential Gene Recruitment by Methylation Status', loc='left', fontsize=14, fontweight='bold')

# Directional Arrows
ax1.annotate('', xy=(0.5, -22), xytext=(0, -22), arrowprops=dict(arrowstyle='<-', color='gray', lw=1.5))
ax1.text(0.25, -24, 'Depleted', ha='center', fontsize=9, color='gray')
ax1.annotate('', xy=(3.5, -22), xytext=(3, -22), arrowprops=dict(arrowstyle='->', color=SIG_COLOR, lw=1.5))
ax1.text(3.25, -24, 'Enriched', ha='center', fontsize=9, color=SIG_COLOR, fontweight='bold')

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=True)
sns.despine(ax=ax1, offset=10)

# --- PANEL B: FOREST PLOT (WITH DIRECTIONAL ANCHORS) ---
ax2 = fig.add_subplot(gs[1, 0])
y_pos = [1, 0]
for i in range(2):
    p = summary.iloc[i]['p_value']
    color = SIG_COLOR if p < 0.05 else NS_COLOR
    or_v, lo, hi = get_stats([brooksi, elizabethae][i])
    ax2.plot([lo, hi], [y_pos[i], y_pos[i]], color=color, lw=3)
    ax2.scatter(or_v, y_pos[i], s=200, color=color, edgecolor='black', zorder=3)
    if p < 0.05: ax2.text(hi + 0.05, y_pos[i], '*', fontsize=20, fontweight='bold', va='center')

ax2.axvline(x=1, color='red', linestyle='--', lw=2)
# Interpretive Anchors
ax2.text(0.5, 1.3, r'$\leftarrow$ Depleted in DEGs', color='gray', ha='center', fontweight='bold')
ax2.text(1.5, 1.3, r'Enriched in DEGs $\rightarrow$', color=SIG_COLOR, ha='center', fontweight='bold')

ax2.set_yticks(y_pos)
ax2.set_yticklabels([r'$\mathit{S.\ brooksi}$', r'$\mathit{S.\ elizabethae}$'], fontsize=12)
ax2.set_xlabel('Odds Ratio (95% CI)', fontweight='bold')
ax2.set_xlim(0, 2.0)
ax2.set_title('B. Evolutionary Conservation of Epigenetic Effect', loc='left', fontsize=14, fontweight='bold')
ax2.grid(axis='x', alpha=0.3)

# --- PANEL C: CLEAN JOURNAL TABLE ---
ax3 = fig.add_subplot(gs[1, 1])
ax3.axis('off')
table_data = []
for i, row in summary.iterrows():
    res = "Significant" if row['p_value'] < 0.05 else "Not Sig."
    table_data.append([f"S. {row['species']}", row['n_genes_analyzed'], f"{meth_rates[i]:.1f}%", f"{row['odds_ratio']:.2f}", f"{row['p_value']:.4f}", res])

table = ax3.table(cellText=table_data, colLabels=['Species', 'N Genes', '% Meth', 'OR', 'P-value', 'Result'], 
                  loc='center', cellLoc='center', bbox=[0, 0, 1, 0.7])
table.auto_set_font_size(False)
table.set_fontsize(10)

# Style header
for i in range(6):
    table[(0, i)].set_facecolor('#34495e')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Style Table
for i in range(1, 3):
    table[(i, 0)].get_text().set_fontstyle('italic')
    res_cell = table[(i, 5)]
    res_cell.set_facecolor('#d5f4e6' if "Significant" in res_cell.get_text().get_text() else '#f2f4f4')

ax3.set_title('C. Statistical Summary Table', loc='left', fontsize=14, fontweight='bold')

plt.savefig('Results/05_ortholog_analysis/figure4_combined_summary.png', bbox_inches='tight', dpi=300)
print("✓ High-bandwidth Figure 4 saved.")
