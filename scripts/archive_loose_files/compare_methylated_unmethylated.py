#!/usr/bin/env python3
"""
Functional comparison: Methylated vs Unmethylated genes
Shows what genes are "ON" (unmethylated) vs "OFF" (methylated)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print("Loading annotation data...")

# Load gene annotations
annotations = pd.read_csv('Results/04_functional_annotation/gene_uniprot_mapping.csv')

print(f"Total annotated genes: {len(annotations):,}")
print(f"Methylated: {len(annotations[annotations['methylation_status']=='methylated']):,}")
print(f"Unmethylated: {len(annotations[annotations['methylation_status']=='unmethylated']):,}")

# Separate methylated and unmethylated
methylated = annotations[annotations['methylation_status'] == 'methylated']
unmethylated = annotations[annotations['methylation_status'] == 'unmethylated']

# Count TEs in each group
def is_te(description):
    """Check if gene is a TE based on description"""
    te_keywords = ['transpos', 'Transpos', 'retrotrans', 'Retrotrans', 
                   'Gag', 'Pol polyprotein', 'reverse transcriptase',
                   'Tigger', 'Mariner', 'LINE', 'SINE', 'LTR', 
                   'Retrovirus', 'Copia', 'Gypsy', 'Helitron', 'integrase']
    return any(keyword in str(description) for keyword in te_keywords)

methylated['is_te'] = methylated['description'].apply(is_te)
unmethylated['is_te'] = unmethylated['description'].apply(is_te)

meth_te_count = methylated['is_te'].sum()
meth_non_te_count = len(methylated) - meth_te_count
unmeth_te_count = unmethylated['is_te'].sum()
unmeth_non_te_count = len(unmethylated) - unmeth_te_count

print(f"\n=== TE COMPOSITION ===")
print(f"Methylated genes:")
print(f"  TEs: {meth_te_count} ({meth_te_count/len(methylated)*100:.1f}%)")
print(f"  Non-TEs: {meth_non_te_count} ({meth_non_te_count/len(methylated)*100:.1f}%)")
print(f"\nUnmethylated genes:")
print(f"  TEs: {unmeth_te_count} ({unmeth_te_count/len(unmethylated)*100:.1f}%)")
print(f"  Non-TEs: {unmeth_non_te_count} ({unmeth_non_te_count/len(unmethylated)*100:.1f}%)")

# Create comparison plot
plt.rcParams.update({
    'font.size': 14,
    'font.weight': 'bold',
    'font.family': 'sans-serif'
})

fig, ax = plt.subplots(figsize=(10, 7))

# Data for grouped bar chart
categories = ['Methylated\n(OFF)', 'Unmethylated\n(ON)']
te_counts = [meth_te_count/len(methylated)*100, unmeth_te_count/len(unmethylated)*100]
non_te_counts = [meth_non_te_count/len(methylated)*100, unmeth_non_te_count/len(unmethylated)*100]

x = np.arange(len(categories))
width = 0.35

bars1 = ax.bar(x - width/2, te_counts, width, label='Transposable Elements', 
               color='#8B0000', alpha=0.8, edgecolor='black', linewidth=2)
bars2 = ax.bar(x + width/2, non_te_counts, width, label='Functional Genes', 
               color='#006400', alpha=0.8, edgecolor='black', linewidth=2)

# Add value labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%',
                ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Percentage of Genes', fontsize=16, fontweight='bold')
ax.set_xlabel('Gene Methylation Status', fontsize=16, fontweight='bold')
ax.set_title('Functional Composition: Methylated vs Unmethylated Genes', 
             fontsize=18, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=14, fontweight='bold')
ax.legend(fontsize=13, loc='upper right', framealpha=0.9)
ax.set_ylim(0, 100)

# Add sample sizes
ax.text(0, -15, f'n={len(methylated)}', ha='center', fontsize=11, 
        transform=ax.transData, fontweight='bold')
ax.text(1, -15, f'n={len(unmethylated)}', ha='center', fontsize=11, 
        transform=ax.transData, fontweight='bold')

ax.tick_params(axis='both', labelsize=12, width=2)
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()
plt.savefig('Results/05_figures/methylated_vs_unmethylated_comparison.png', 
            dpi=300, bbox_inches='tight', facecolor='white')
print("\n✅ Comparison plot saved: methylated_vs_unmethylated_comparison.png")

# Get top functional categories for unmethylated genes (non-TEs)
print("\n=== TOP FUNCTIONS IN UNMETHYLATED GENES (Non-TEs) ===")
unmeth_non_te = unmethylated[~unmethylated['is_te']]

# Extract general function categories from descriptions
def extract_function(desc):
    """Extract general function from description"""
    desc_str = str(desc).lower()
    if 'kinase' in desc_str:
        return 'Kinase'
    elif 'ribosom' in desc_str:
        return 'Ribosomal protein'
    elif 'transcription' in desc_str:
        return 'Transcription factor'
    elif 'metabol' in desc_str or 'synthase' in desc_str or 'dehydrogenase' in desc_str:
        return 'Metabolism'
    elif 'transport' in desc_str or 'channel' in desc_str:
        return 'Transport'
    elif 'signal' in desc_str:
        return 'Signaling'
    elif 'structural' in desc_str or 'cytoskeleton' in desc_str or 'actin' in desc_str or 'tubulin' in desc_str:
        return 'Structural/Cytoskeleton'
    elif 'dna repair' in desc_str or 'helicase' in desc_str:
        return 'DNA repair'
    elif 'protease' in desc_str or 'peptidase' in desc_str:
        return 'Protease'
    else:
        return 'Other'

unmeth_non_te['function_category'] = unmeth_non_te['description'].apply(extract_function)
function_counts = unmeth_non_te['function_category'].value_counts().head(10)

print("\nTop 10 functional categories in unmethylated genes:")
for func, count in function_counts.items():
    pct = count / len(unmeth_non_te) * 100
    print(f"  {func}: {count} ({pct:.1f}%)")

# Save summary statistics
summary = pd.DataFrame({
    'Category': ['Methylated (OFF)', 'Unmethylated (ON)'],
    'Total_Genes': [len(methylated), len(unmethylated)],
    'TE_Count': [meth_te_count, unmeth_te_count],
    'TE_Percentage': [meth_te_count/len(methylated)*100, unmeth_te_count/len(unmethylated)*100],
    'Functional_Count': [meth_non_te_count, unmeth_non_te_count],
    'Functional_Percentage': [meth_non_te_count/len(methylated)*100, unmeth_non_te_count/len(unmethylated)*100]
})

summary.to_csv('Results/04_functional_annotation/methylated_vs_unmethylated_summary.csv', index=False)
print("\n✅ Summary statistics saved: methylated_vs_unmethylated_summary.csv")

plt.show()

print("\n" + "="*60)
print("KEY FINDING:")
print(f"Methylated genes: {meth_te_count/len(methylated)*100:.1f}% TEs (SILENCED)")
print(f"Unmethylated genes: {unmeth_non_te_count/len(unmethylated)*100:.1f}% functional (ACTIVE)")
print("="*60)
