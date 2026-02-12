#!/usr/bin/env python3
"""
Generate CpG island distribution plots for poster
1. CpG island size distribution
2. Methylation level distribution across CpG sites
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print("Loading data...")

# Load CpG island data
cpg_islands = pd.read_csv('Results/03_methylation_analysis/high_methylated_cpg_islands.bed', 
                          sep='\t', header=None, names=['chrom', 'start', 'end'])
cpg_islands['length'] = cpg_islands['end'] - cpg_islands['start']

# Load methylation data
meth_data = pd.read_csv('Results/03_methylation_analysis/filtered_bedmethyl_50p_cov30.csv')
# Column is 'percent_modified' not 'methylation_pct'
meth_data['methylation_pct'] = meth_data['percent_modified']

print(f"CpG islands: {len(cpg_islands):,}")
print(f"Methylation sites: {len(meth_data):,}")

# Set up plotting style
plt.rcParams.update({
    'font.size': 14,
    'font.weight': 'bold',
    'font.family': 'sans-serif'
})

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# ===== PLOT 1: CpG Island Size Distribution =====
ax1.hist(cpg_islands['length'], bins=50, color='steelblue', alpha=0.8, edgecolor='black')
ax1.set_xlabel('CpG Island Length (bp)', fontsize=16, fontweight='bold')
ax1.set_ylabel('Count', fontsize=16, fontweight='bold')
ax1.set_title('CpG Island Size Distribution\n(Methylated Islands)', 
              fontsize=18, fontweight='bold', pad=15)

# Add statistics
median_len = cpg_islands['length'].median()
mean_len = cpg_islands['length'].mean()
ax1.axvline(median_len, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax1.text(0.65, 0.95, f'Median: {median_len:.0f} bp\nMean: {mean_len:.0f} bp\nn = {len(cpg_islands):,}',
         transform=ax1.transAxes, fontsize=12, fontweight='bold',
         verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.tick_params(axis='both', labelsize=12, width=2)
for spine in ax1.spines.values():
    spine.set_linewidth(2)

# ===== PLOT 2: Methylation Level Distribution =====
ax2.hist(meth_data['methylation_pct'], bins=50, color='darkgreen', alpha=0.8, edgecolor='black')
ax2.set_xlabel('Methylation Level (%)', fontsize=16, fontweight='bold')
ax2.set_ylabel('Count', fontsize=16, fontweight='bold')
ax2.set_title('Methylation Level Distribution\n(≥50% Methylation, ≥30× Coverage)', 
              fontsize=18, fontweight='bold', pad=15)

# Add statistics
median_meth = meth_data['methylation_pct'].median()
mean_meth = meth_data['methylation_pct'].mean()
ax2.axvline(median_meth, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax2.text(0.05, 0.95, f'Median: {median_meth:.1f}%\nMean: {mean_meth:.1f}%\nn = {len(meth_data):,}',
         transform=ax2.transAxes, fontsize=12, fontweight='bold',
         verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax2.tick_params(axis='both', labelsize=12, width=2)
for spine in ax2.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()

# Save figure
output_file = 'Results/05_figures/cpg_methylation_distributions.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✅ Distribution plots saved: {output_file}")

# Also save as individual plots for flexibility
fig1, ax = plt.subplots(figsize=(8, 6))
ax.hist(cpg_islands['length'], bins=50, color='steelblue', alpha=0.8, edgecolor='black')
ax.set_xlabel('CpG Island Length (bp)', fontsize=16, fontweight='bold')
ax.set_ylabel('Count', fontsize=16, fontweight='bold')
ax.set_title('CpG Island Size Distribution\n(Methylated Islands)', 
             fontsize=18, fontweight='bold', pad=15)
ax.axvline(median_len, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.text(0.65, 0.95, f'Median: {median_len:.0f} bp\nMean: {mean_len:.0f} bp\nn = {len(cpg_islands):,}',
        transform=ax.transAxes, fontsize=12, fontweight='bold',
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.tick_params(axis='both', labelsize=12, width=2)
for spine in ax.spines.values():
    spine.set_linewidth(2)
plt.tight_layout()
plt.savefig('Results/05_figures/cpg_island_size_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
print("✅ Individual plot saved: cpg_island_size_distribution.png")

fig2, ax = plt.subplots(figsize=(8, 6))
ax.hist(meth_data['methylation_pct'], bins=50, color='darkgreen', alpha=0.8, edgecolor='black')
ax.set_xlabel('Methylation Level (%)', fontsize=16, fontweight='bold')
ax.set_ylabel('Count', fontsize=16, fontweight='bold')
ax.set_title('Methylation Level Distribution\n(≥50% Methylation, ≥30× Coverage)', 
             fontsize=18, fontweight='bold', pad=15)
ax.axvline(median_meth, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax.text(0.05, 0.95, f'Median: {median_meth:.1f}%\nMean: {mean_meth:.1f}%\nn = {len(meth_data):,}',
        transform=ax.transAxes, fontsize=12, fontweight='bold',
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.tick_params(axis='both', labelsize=12, width=2)
for spine in ax.spines.values():
    spine.set_linewidth(2)
plt.tight_layout()
plt.savefig('Results/05_figures/methylation_level_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
print("✅ Individual plot saved: methylation_level_distribution.png")

plt.show()

print("\n📊 Summary Statistics:")
print(f"\nCpG Islands:")
print(f"  Total: {len(cpg_islands):,}")
print(f"  Median length: {median_len:.0f} bp")
print(f"  Mean length: {mean_len:.0f} bp")
print(f"  Min length: {cpg_islands['length'].min():.0f} bp")
print(f"  Max length: {cpg_islands['length'].max():.0f} bp")

print(f"\nMethylation Levels:")
print(f"  Total sites: {len(meth_data):,}")
print(f"  Median: {median_meth:.1f}%")
print(f"  Mean: {mean_meth:.1f}%")
print(f"  Min: {meth_data['percent_modified'].min():.1f}%")
print(f"  Max: {meth_data['percent_modified'].max():.1f}%")
