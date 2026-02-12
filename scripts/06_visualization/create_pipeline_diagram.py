#!/usr/bin/env python3
"""
Generate Methods Pipeline Diagram - Detailed Version
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(16, 9), dpi=300)

# Pipeline boxes with tools/methods
boxes = [
    {'x': 1.5, 'y': 6.5, 'text': 'ONT\nSequencing', 'stat': 'Draft Genome', 'tool': 'MinION'},
    {'x': 4.5, 'y': 6.5, 'text': 'Gene\nPrediction', 'stat': '6,848 genes', 'tool': 'AUGUSTUS'},
    {'x': 7.5, 'y': 6.5, 'text': 'Methylation\nCalling', 'stat': '2.6M CpG sites', 'tool': '≥30× cov, ≥50% meth'},
    {'x': 10.5, 'y': 6.5, 'text': 'Promoter\nAnalysis', 'stat': '1,103 methylated (16%)', 'tool': 'Upstream regions'},
    {'x': 13.5, 'y': 6.5, 'text': 'Functional\nAnnotation', 'stat': '5,504 annotated (80%)', 'tool': 'DIAMOND + SwissProt'},
    {'x': 7.5, 'y': 2.5, 'text': 'GO\nEnrichment', 'stat': '47 terms enriched', 'tool': "Fisher's exact test"},
]

# Draw boxes
for box in boxes:
    rect = FancyBboxPatch((box['x']-1.2, box['y']-0.7), 2.4, 1.4,
                          boxstyle="round,pad=0.1", 
                          edgecolor='#2E86AB', facecolor='#E8F4F8',
                          linewidth=3)
    ax.add_patch(rect)
    # Step name
    ax.text(box['x'], box['y']+0.3, box['text'], ha='center', va='center',
            fontsize=16, fontweight='bold', color='#1a1a1a')
    # Tool/method
    ax.text(box['x'], box['y']-0.05, box['tool'], ha='center', va='center',
            fontsize=11, color='#555555', style='italic')
    # Result/stat
    ax.text(box['x'], box['y']-0.45, box['stat'], ha='center', va='center',
            fontsize=13, color='#A23B72', fontweight='bold')

# Draw arrows (top row)
for i in range(4):
    arrow = FancyArrowPatch((boxes[i]['x']+1.3, boxes[i]['y']), 
                           (boxes[i+1]['x']-1.3, boxes[i+1]['y']),
                           arrowstyle='->', mutation_scale=35, 
                           linewidth=4, color='#2E86AB')
    ax.add_patch(arrow)

# Draw arrow from Annotation to GO Enrichment
arrow_down = FancyArrowPatch((boxes[4]['x'], boxes[4]['y']-0.8), 
                            (boxes[5]['x'], boxes[5]['y']+0.8),
                            arrowstyle='->', mutation_scale=35, 
                            linewidth=4, color='#2E86AB')
ax.add_patch(arrow_down)

# Title
ax.text(7.5, 8.2, 'Genome-Wide Methylation Analysis Pipeline', ha='center', va='center',
        fontsize=28, fontweight='bold', color='#1a1a1a')

# Key finding box (LARGER)
finding_text = 'Key Finding:\n77% of annotated\nmethylated genes\nare TEs'
ax.text(1.5, 1.2, finding_text, ha='center', va='center',
        fontsize=20, color='#1a1a1a', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.6', facecolor='#FFE6F0', 
                 edgecolor='#A23B72', linewidth=4))

ax.set_xlim(0, 15)
ax.set_ylim(0, 9)
ax.axis('off')

plt.tight_layout(pad=0)
plt.savefig('Results/05_figures/pipeline_diagram.png', dpi=300, 
            bbox_inches='tight', facecolor='white')
print("✓ Pipeline diagram saved: Results/05_figures/pipeline_diagram.png")
