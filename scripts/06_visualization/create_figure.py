#!/usr/bin/env python3
"""Create GO Enrichment Bar Chart"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT = Path(__file__).parent / "analysis_output"

# Load results
df = pd.read_csv(OUTPUT / "go_enrichment_results.csv")
sig = df[df['significant'] & (df['enriched_in'] == 'methylated')].head(15)

# Calculate -log10(p-value)
sig['neg_log_p'] = -np.log10(sig['p_value'])

# Create bar chart
fig, ax = plt.subplots(figsize=(10, 8))
y_pos = range(len(sig))
ax.barh(y_pos, sig['neg_log_p'], color='steelblue')
ax.set_yticks(y_pos)
ax.set_yticklabels(sig['go_term'], fontsize=10)
ax.set_xlabel('-log10(p-value)', fontsize=12)
ax.set_title('GO Terms Enriched in Methylated Gene Promoters\n(Synalpheus chacei)', fontsize=14, fontweight='bold')
ax.axvline(x=1.3, color='red', linestyle='--', linewidth=1, label='p=0.05')
ax.legend()
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(OUTPUT / "go_enrichment_chart.png", dpi=300, bbox_inches='tight')
print(f"✅ Chart saved: {OUTPUT / 'go_enrichment_chart.png'}")
