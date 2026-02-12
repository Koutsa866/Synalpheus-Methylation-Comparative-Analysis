#!/usr/bin/env python3
"""
Create professional donut chart for PAG33 poster
Shows promoter methylation distribution
"""

import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent.parent / "Results" / "05_figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Data from results
labels = ['Methylated Promoters', 'Unmethylated Promoters']
counts = [1103, 5745]
colors = ['#003366', '#D3D3D3']  # Navy Blue and Light Gray
explode = (0.1, 0)  # "Pop out" the methylated slice

# Styling
plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif'})

fig, ax = plt.subplots(figsize=(8, 8))

# Create the pie chart
wedges, texts, autotexts = ax.pie(
    counts, 
    labels=labels, 
    autopct='%1.1f%%', 
    startangle=90, 
    colors=colors, 
    explode=explode, 
    pctdistance=0.85,
    textprops={'weight': 'bold'}
)

# Draw a circle in the center to make it a donut
centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig.gca().add_artist(centre_circle)

# Add Total N in the center
plt.text(0, 0, f"Total Genes\n{sum(counts):,}", 
         ha='center', va='center', fontsize=18, weight='bold')

plt.title("Genomic Distribution of Promoter Methylation", 
          fontsize=20, pad=20, weight='bold')

# Save as high-res PNG and PDF for PAG
output_png = OUTPUT_DIR / 'promoter_methylation_donut.png'
output_pdf = OUTPUT_DIR / 'promoter_methylation_donut.pdf'

plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.savefig(output_pdf, bbox_inches='tight')

print(f"✅ Donut chart saved:")
print(f"   PNG: {output_png}")
print(f"   PDF: {output_pdf}")

plt.show()
