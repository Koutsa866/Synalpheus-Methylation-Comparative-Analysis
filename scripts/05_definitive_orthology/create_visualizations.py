#!/usr/bin/env python3
"""
Pipeline: 05_definitive_orthology
Phase: 4 - Visualization
Goal: Generate publication-quality figures for thesis defense
Note: Figure 2 (Forest Plot) is the primary result showing lineage-specific epigenetic conservation
Inputs:
  - Results/05_ortholog_analysis/brooksi_contingency_table.csv
  - Results/05_ortholog_analysis/elizabethae_contingency_table.csv
  - Results/05_ortholog_analysis/summary_statistics.csv
Outputs:
  - Results/05_ortholog_analysis/figure1_enrichment_plot.png
  - Results/05_ortholog_analysis/figure2_forest_plot.png (PRIMARY RESULT)
  - Results/05_ortholog_analysis/figure3_mosaic_plot.png
  - Results/05_ortholog_analysis/figure4_combined_summary.png
Author: Philip Koutsaftis
Date: 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
import seaborn as sns

# --- STYLE CONFIGURATION ---
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 300
# High-bandwidth color palette
DEG_COLOR = '#e74c3c'      # Salmon/Red
NON_DEG_COLOR = '#3498db'  # Blue
SIG_COLOR = '#27ae60'      # Green
NS_COLOR = '#95a5a6'       # Gray

# Load contingency tables
brooksi = pd.read_csv('Results/05_ortholog_analysis/brooksi_contingency_table.csv', index_col=0)
elizabethae = pd.read_csv('Results/05_ortholog_analysis/elizabethae_contingency_table.csv', index_col=0)

# Load summary stats
summary = pd.read_csv('Results/05_ortholog_analysis/summary_statistics.csv')

def calculate_or_ci(ctable):
    """Calculate Odds Ratio and 95% CI"""
    a, b = ctable.iloc[0, 0], ctable.iloc[0, 1]
    c, d = ctable.iloc[1, 0], ctable.iloc[1, 1]
    or_val = (a * d) / (b * c) if (b * c) != 0 else float('inf')
    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
    log_or = np.log(or_val)
    return or_val, np.exp(log_or - 1.96 * se_log_or), np.exp(log_or + 1.96 * se_log_or)

# Prep stats
ors, lows, highs = [], [], []
for ct in [brooksi, elizabethae]:
    o, l, h = calculate_or_ci(ct)
    ors.append(o); lows.append(l); highs.append(h)

# ============================================================
# FIGURE 1: Normalized Stacked Bar Chart (Enrichment Plot)
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for idx, (species, ctable) in enumerate([('S. brooksi', brooksi), ('S. elizabethae', elizabethae)]):
    # Convert to percentages
    prop_table = ctable.div(ctable.sum(axis=1), axis=0) * 100
    
    # Plot stacked bar
    prop_table.plot(kind='bar', stacked=True, ax=axes[idx], 
                    color=['#e74c3c', '#3498db'], width=0.6)
    
    axes[idx].set_title(f'$\mathit{{{species}}}$', fontsize=14, fontweight='bold')
    axes[idx].set_xlabel('Methylation Status', fontsize=12)
    axes[idx].set_ylabel('Percentage of Genes (%)', fontsize=12)
    axes[idx].set_xticklabels(axes[idx].get_xticklabels(), rotation=0)
    axes[idx].legend(title='Expression Status', loc='upper left', bbox_to_anchor=(0, 1))
    axes[idx].set_ylim(0, 100)
    
    # Add percentage labels
    for container in axes[idx].containers:
        axes[idx].bar_label(container, fmt='%.1f%%', label_type='center')

plt.tight_layout()
plt.savefig('Results/05_ortholog_analysis/figure1_enrichment_plot.png', dpi=300, bbox_inches='tight')
print("✓ Figure 1: Enrichment plot saved")

# ============================================================
# FIGURE 2: Forest Plot (Odds Ratio with 95% CI)
# ============================================================

def calculate_or_ci(ctable):
    """Calculate Odds Ratio and 95% CI"""
    a, b = ctable.iloc[0, 0], ctable.iloc[0, 1]  # Meth-DEG, Meth-NonDEG
    c, d = ctable.iloc[1, 0], ctable.iloc[1, 1]  # Unmeth-DEG, Unmeth-NonDEG
    
    # Odds Ratio
    or_val = (a * d) / (b * c) if (b * c) != 0 else float('inf')
    
    # Standard Error of log(OR)
    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
    
    # 95% CI
    log_or = np.log(or_val)
    ci_lower = np.exp(log_or - 1.96 * se_log_or)
    ci_upper = np.exp(log_or + 1.96 * se_log_or)
    
    return or_val, ci_lower, ci_upper

# Calculate ORs and CIs
brooksi_or, brooksi_ci_low, brooksi_ci_high = calculate_or_ci(brooksi)
eliz_or, eliz_ci_low, eliz_ci_high = calculate_or_ci(elizabethae)

# Create forest plot
fig, ax = plt.subplots(figsize=(10, 6))

species_names = ['$\mathit{S. brooksi}$', '$\mathit{S. elizabethae}$']
ors = [brooksi_or, eliz_or]
ci_lows = [brooksi_ci_low, eliz_ci_low]
ci_highs = [brooksi_ci_high, eliz_ci_high]
p_values = summary['p_value'].tolist()

# Plot points and error bars
y_positions = [1, 0]
colors = ['#95a5a6' if p > 0.05 else '#27ae60' for p in p_values]

for i, (y, or_val, ci_low, ci_high, color, p_val) in enumerate(zip(y_positions, ors, ci_lows, ci_highs, colors, p_values)):
    # Error bar
    ax.plot([ci_low, ci_high], [y, y], color=color, linewidth=2)
    # Point
    ax.scatter(or_val, y, s=200, color=color, zorder=3, edgecolors='black', linewidth=1.5)
    # Add OR and p-value text
    sig = "*" if p_val < 0.05 else "ns"
    ax.text(ci_high + 0.15, y, f'OR={or_val:.2f}\np={p_val:.4f} {sig}', 
            va='center', fontsize=10)

# Reference line at OR=1
ax.axvline(x=1, color='red', linestyle='--', linewidth=1.5, label='OR = 1 (No Effect)')

# Formatting
ax.set_yticks(y_positions)
ax.set_yticklabels(species_names, fontsize=12)
ax.set_xlabel('Odds Ratio (95% CI)', fontsize=12, fontweight='bold')
ax.set_title('Methylation Enrichment in DEGs Across Species', fontsize=14, fontweight='bold')
ax.set_xlim(0, 2.0)  # Fixed range for clarity
ax.legend(loc='upper right')
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig('Results/05_ortholog_analysis/figure2_forest_plot.png', dpi=300, bbox_inches='tight')
print("✓ Figure 2: Forest plot saved")

# ============================================================
# FIGURE 3: Mosaic Plot (Chi-Square Visualization)
# ============================================================

from statsmodels.graphics.mosaicplot import mosaic

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, (species, ctable) in enumerate([('$\mathit{S. brooksi}$', brooksi), ('$\mathit{S. elizabethae}$', elizabethae)]):
    # Calculate standardized residuals for coloring
    chi2, p, dof, expected = chi2_contingency(ctable)
    residuals = (ctable - expected) / np.sqrt(expected)
    
    # Prepare data for mosaic
    data = {}
    for meth_status in ctable.index:
        for exp_status in ctable.columns:
            data[(meth_status, exp_status)] = ctable.loc[meth_status, exp_status]
    
    # Create mosaic
    mosaic(data, ax=axes[idx], title=f'{species}\n(p={summary.iloc[idx]["p_value"]:.4f})',
           labelizer=lambda k: f'{k[0]}\n{k[1]}',
           properties=lambda key: {'color': 'lightblue' if residuals.loc[key] > 0 else 'lightcoral'})
    
    axes[idx].set_xlabel('')
    axes[idx].set_ylabel('')

plt.tight_layout()
plt.savefig('Results/05_ortholog_analysis/figure3_mosaic_plot.png', dpi=300, bbox_inches='tight')
print("✓ Figure 3: Mosaic plot saved")

# ============================================================
# BONUS: Combined Summary Figure
# ============================================================

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# Panel A: Enrichment bars
ax1 = fig.add_subplot(gs[0, :])
for idx, (species, ctable) in enumerate([('S. brooksi', brooksi), ('S. elizabethae', elizabethae)]):
    prop_table = ctable.div(ctable.sum(axis=1), axis=0) * 100
    x_offset = idx * 3
    x_pos = np.array([0, 1]) + x_offset
    
    ax1.bar(x_pos, prop_table['DEG'], width=0.6, label=f'{species} - DEG' if idx == 0 else None,
            color=['#e74c3c', '#e74c3c'], alpha=0.7, edgecolor='black')
    ax1.bar(x_pos, prop_table['Non-DEG'], width=0.6, bottom=prop_table['DEG'],
            label=f'{species} - Non-DEG' if idx == 0 else None,
            color=['#3498db', '#3498db'], alpha=0.7, edgecolor='black')
    
    # Add species labels
    ax1.text(x_offset + 0.5, -8, f'$\mathit{{{species}}}$', ha='center', fontsize=12, fontweight='bold')
    
    # Add statistics annotation
    or_val = ors[idx]
    p_val = p_values[idx]
    sig_marker = "*" if p_val < 0.05 else "ns"
    stats_text = f'OR={or_val:.2f}, p={p_val:.3f} {sig_marker}'
    text_color = '#27ae60' if p_val < 0.05 else '#95a5a6'
    ax1.text(x_offset + 0.5, 105, stats_text, ha='center', fontsize=10, 
             fontweight='bold', color=text_color, bbox=dict(boxstyle='round', 
             facecolor='white', edgecolor=text_color, linewidth=1.5))

ax1.set_xticks([0, 1, 3, 4])
ax1.set_xticklabels(['Methylated', 'Unmethylated', 'Methylated', 'Unmethylated'])
ax1.set_ylabel('Percentage (%)', fontsize=12)
ax1.set_title('A. DEG Enrichment by Methylation Status', fontsize=14, fontweight='bold', loc='left')
ax1.set_ylim(0, 115)  # Increased to accommodate statistics text
ax1.legend(loc='upper right')

# Panel B: Forest plot
ax2 = fig.add_subplot(gs[1, 0])
for i, (y, or_val, ci_low, ci_high, color, p_val) in enumerate(zip(y_positions, ors, ci_lows, ci_highs, colors, p_values)):
    ax2.plot([ci_low, ci_high], [y, y], color=color, linewidth=2)
    ax2.scatter(or_val, y, s=150, color=color, zorder=3, edgecolors='black', linewidth=1.5)
ax2.axvline(x=1, color='red', linestyle='--', linewidth=1.5)
ax2.set_yticks(y_positions)
ax2.set_yticklabels(species_names)
ax2.set_xlabel('Odds Ratio (95% CI)', fontsize=11)
ax2.set_title('B. Effect Size Comparison', fontsize=14, fontweight='bold', loc='left')
ax2.grid(axis='x', alpha=0.3)

# Panel C: Summary statistics table
ax3 = fig.add_subplot(gs[1, 1])
ax3.axis('off')

table_data = []
for _, row in summary.iterrows():
    sig = "✓ Significant" if row['significant'] else "✗ Not Significant"
    table_data.append([
        row['species'].replace('_', ' ').title(),
        f"{row['n_genes_analyzed']}",
        f"{row['odds_ratio']:.2f}",
        f"{row['p_value']:.4f}",
        sig
    ])

table = ax3.table(cellText=table_data,
                  colLabels=['Species', 'N Genes', 'Odds Ratio', 'P-value', 'Result'],
                  cellLoc='center',
                  loc='center',
                  bbox=[0, 0, 1, 1])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# Style header
for i in range(5):
    table[(0, i)].set_facecolor('#34495e')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Color-code results
for i in range(1, 3):
    if table_data[i-1][4] == "✓ Significant":
        table[(i, 4)].set_facecolor('#d5f4e6')
    else:
        table[(i, 4)].set_facecolor('#fadbd8')

ax3.set_title('C. Statistical Summary', fontsize=14, fontweight='bold', loc='left', pad=20)

plt.savefig('Results/05_ortholog_analysis/figure4_combined_summary.png', dpi=300, bbox_inches='tight')
print("✓ Figure 4: Combined summary figure saved")

print("\n✓ All visualizations complete!")
print("  Saved to: Results/05_ortholog_analysis/")
