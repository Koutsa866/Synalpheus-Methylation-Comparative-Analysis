#!/usr/bin/env python3
"""
Generate Eusociality Evolution Diagram
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

# Species groups with eusociality
species = [
    {'name': 'Ants', 'y': 6, 'icon': '🐜', 'group': 'Insects'},
    {'name': 'Bees', 'y': 5, 'icon': '🐝', 'group': 'Insects'},
    {'name': 'Termites', 'y': 4, 'icon': '🪲', 'group': 'Insects'},
    {'name': 'Naked Mole Rats', 'y': 2.5, 'icon': '🐭', 'group': 'Mammals'},
    {'name': 'Synalpheus Shrimp', 'y': 1, 'icon': '🦐', 'group': 'Marine'},
]

# Draw species boxes
for sp in species:
    # Color by group
    if sp['group'] == 'Insects':
        color = '#FFD700'
    elif sp['group'] == 'Mammals':
        color = '#FF6B6B'
    else:
        color = '#2E86AB'
    
    rect = mpatches.FancyBboxPatch((8, sp['y']-0.3), 5, 0.6,
                                   boxstyle="round,pad=0.1",
                                   edgecolor='#333', facecolor=color,
                                   linewidth=2, alpha=0.7)
    ax.add_patch(rect)
    ax.text(9, sp['y'], sp['icon'], ha='center', va='center', fontsize=32)
    ax.text(10.5, sp['y'], sp['name'], ha='center', va='center', 
            fontsize=18, fontweight='bold', color='#1a1a1a')

# Convergent evolution arrows
for sp in species:
    ax.annotate('', xy=(7.8, sp['y']), xytext=(3, 3.5),
                arrowprops=dict(arrowstyle='->', lw=2, color='#666', alpha=0.5))

# Central point
circle = mpatches.Circle((3, 3.5), 0.4, edgecolor='#A23B72', 
                        facecolor='#FFE6F0', linewidth=3)
ax.add_patch(circle)
ax.text(3, 3.5, '?', ha='center', va='center', fontsize=32, 
        fontweight='bold', color='#A23B72')

# Labels
ax.text(3, 2.5, 'Convergent\nEvolution', ha='center', fontsize=16, 
        fontweight='bold', color='#A23B72')

# Title
ax.text(7, 7.5, 'Independent Origins of Eusociality', ha='center', 
        fontsize=26, fontweight='bold', color='#1a1a1a')

# Key features box
features = '• Reproductive division of labor\n• Cooperative brood care\n• Overlapping generations'
ax.text(1, 6.5, features, fontsize=14, color='#333',
        bbox=dict(boxstyle='round', facecolor='#F0F0F0', 
                 edgecolor='#2E86AB', linewidth=2))

ax.set_xlim(0, 14)
ax.set_ylim(0, 8)
ax.axis('off')

plt.tight_layout(pad=0)
plt.savefig('Results/05_figures/eusociality_evolution.png', dpi=300, 
            bbox_inches='tight', facecolor='white')
print("✓ Eusociality diagram saved: Results/05_figures/eusociality_evolution.png")
