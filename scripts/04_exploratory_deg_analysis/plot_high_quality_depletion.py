#!/usr/bin/env python3
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Set the style
sns.set_theme(style="white")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# --- PLOT 1: The "Complete Depletion" Bar ---
categories = ['Genome-wide\n(n=6,848)', 'High-Quality DEGs\n(n=13)']
values = [16.1, 0.0]
colors = ['#34495e', '#e74c3c']

sns.barplot(x=categories, y=values, palette=colors, ax=ax1, edgecolor='black', linewidth=1.5)
ax1.set_ylim(0, 20)
ax1.set_ylabel('Methylation Rate (%)', fontsize=12, fontweight='bold')
ax1.set_title('A. Genomic Methylation Depletion', fontsize=14, fontweight='bold')

# Add genome baseline reference
ax1.axhline(y=16.1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(1.5, 16.5, 'Genome baseline', fontsize=9, style='italic', color='gray')

# Add value labels
ax1.text(0, 16.5, '16.1%', ha='center', va='bottom', fontsize=12, fontweight='bold')
ax1.text(1, 0.5, '0.0%', ha='center', va='bottom', fontsize=12, fontweight='bold', color='red')

# Add depletion annotation
ax1.text(0.5, 18, 'Complete\nDepletion', ha='center', fontsize=11, 
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

# --- PLOT 2: Individual Gene Dots ---
# Creating a dummy dataframe for the 13 genes
data = {
    'Caste': ['Queen-up']*2 + ['Worker-up']*11,
    'Status': [0]*13  # 0 for Unmethylated
}
df_dots = pd.DataFrame(data)

sns.stripplot(data=df_dots, x='Caste', y='Status', jitter=0.2, size=12, 
              palette=['#f39c12', '#3498db'], ax=ax2, linewidth=1, edgecolor='black')

# Only show unmethylated region (no genes are methylated)
ax2.set_yticks([0])
ax2.set_yticklabels(['Unmethylated\n(13/13 genes)'], fontweight='bold')
ax2.set_ylim(-0.3, 0.5)
ax2.set_title('B. Status of Individual DEGs', fontsize=14, fontweight='bold')
ax2.set_xlabel('Regulation Direction', fontsize=12)
ax2.set_ylabel('')

# Annotate n values
ax2.text(0, -0.2, 'n=2', ha='center', fontsize=10, style='italic')
ax2.text(1, -0.2, 'n=11', ha='center', fontsize=10, style='italic')

# Add note about methylated genes
ax2.text(0.5, 0.4, 'No methylated DEGs detected', ha='center', fontsize=9, 
         style='italic', color='red')

plt.tight_layout()
plt.savefig('Results/05_figures/deg_methylation_depletion_high_quality.png', dpi=300, bbox_inches='tight')
plt.savefig('Results/05_figures/deg_methylation_depletion_high_quality.pdf', bbox_inches='tight')
print("✓ Saved to Results/05_figures/deg_methylation_depletion_high_quality.png")
print("✓ Saved to Results/05_figures/deg_methylation_depletion_high_quality.pdf")
plt.show()
