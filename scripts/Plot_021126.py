import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# Ensure directory exists
output_dir = 'Results/05_figures/'
os.makedirs(output_dir, exist_ok=True)

# Set the style
sns.set_theme(style="white")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# --- PLOT 1: The "Complete Depletion" Bar ---
categories = ['Genome-wide\n(n=6,848)', 'High-Quality DEGs\n(n=13)']
values = [16.1, 0.0]
colors = ['#34495e', '#e74c3c']

sns.barplot(x=categories, y=values, palette=colors, ax=ax1, edgecolor='black', linewidth=1.5)

# Add baseline reference line
ax1.axhline(y=16.1, color='gray', linestyle='--', linewidth=1.5, alpha=0.6, label='Genome Baseline (16.1%)')

ax1.set_ylim(0, 22) # Extra room for labels
ax1.set_ylabel('Methylation Rate (%)', fontsize=12, fontweight='bold')
ax1.set_title('A. Genomic Methylation Depletion', fontsize=14, fontweight='bold')
ax1.legend(loc='upper right', frameon=False)

# Add value labels and callout
ax1.text(0, 16.5, '16.1%', ha='center', va='bottom', fontsize=12, fontweight='bold')
ax1.text(1, 1.0, 'COMPLETE\nDEPLETION', ha='center', va='bottom', fontsize=10, fontweight='bold', color='#e74c3c')
ax1.annotate('0.0%', xy=(1, 0), xytext=(1, 4), arrowprops=dict(arrowstyle='->', color='black'), ha='center', fontweight='bold')

# --- PLOT 2: Individual Gene Dots (Standard Binary Format) ---
# Creating dataframe for the 13 genes
data = {
    'Caste': ['Queen-up']*2 + ['Worker-up']*11,
    'Status': ['Unmethylated']*13  # All unmethylated
}
df_dots = pd.DataFrame(data)

sns.stripplot(data=df_dots, x='Caste', y='Status', jitter=0.25, size=14, 
              palette=['#f39c12', '#3498db'], ax=ax2, linewidth=1.5, edgecolor='black', alpha=0.8)

# Show both categories (standard binary visualization)
ax2.set_ylim(-0.5, 1.5)
ax2.set_yticks([0, 1])
ax2.set_yticklabels(['Unmethylated\n(13/13 genes)', 'Methylated\n(0/13 genes)'], fontweight='bold')

ax2.set_title('B. Methylation Status of High-Quality DEGs', fontsize=14, fontweight='bold')
ax2.set_xlabel('Regulation Direction', fontsize=12, fontweight='bold')
ax2.set_ylabel('')

# Add shaded region for "no data" zone
ax2.axhspan(0.5, 1.5, alpha=0.1, color='red', zorder=0)
ax2.text(0.5, 1.0, 'No genes in this category', ha='center', fontsize=10, 
         color='#e74c3c', style='italic', fontweight='bold')

# Annotate n values
ax2.text(0, -0.35, 'n=2', ha='center', fontsize=11, fontweight='bold')
ax2.text(1, -0.35, 'n=11', ha='center', fontsize=11, fontweight='bold')

plt.tight_layout()

# Add figure-level caption about ortholog quality
fig.text(0.5, 0.02, 'High-confidence orthologs only (E-value < 1e-20, identity > 60%). Median identity: 92.8% (S. brooksi), 84.1% (S. elizabethae)',
         ha='center', fontsize=9, style='italic', color='gray')

# Save in both formats
save_path = os.path.join(output_dir, 'deg_methylation_depletion_high_quality')
plt.savefig(f'{save_path}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{save_path}.pdf', bbox_inches='tight')

print(f"Success: Figures saved to {output_dir}")
plt.show()