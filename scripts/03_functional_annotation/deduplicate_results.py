#!/usr/bin/env python3
"""Deduplicate DIAMOND results - keep best hit per gene"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]

# Load new results
df = pd.read_csv(BASE / "Results/04_functional_annotation/gene_uniprot_mapping_NEW.csv")

print(f"Before deduplication: {len(df)} rows, {df['gene_id'].nunique()} unique genes")

# Keep best hit per gene (lowest e-value, highest bitscore)
df_dedup = df.sort_values(['gene_id', 'evalue', 'bitscore'], ascending=[True, True, False])
df_dedup = df_dedup.drop_duplicates(subset='gene_id', keep='first')

print(f"After deduplication: {len(df_dedup)} rows, {df_dedup['gene_id'].nunique()} unique genes")

# Save deduplicated results
output = BASE / "Results/04_functional_annotation/gene_uniprot_mapping_DEDUP.csv"
df_dedup.to_csv(output, index=False)

# Report statistics
meth = df_dedup[df_dedup['methylation_status'] == 'methylated']
unmeth = df_dedup[df_dedup['methylation_status'] == 'unmethylated']

print(f"\n✅ DEDUPLICATED RESULTS:")
print(f"   Total unique genes: {len(df_dedup)}")
print(f"   Methylated: {len(meth)}")
print(f"   Unmethylated: {len(unmeth)}")

# Check TE enrichment
te_pattern = r'[Tt]ranspos|[Rr]etrotrans|Gag|Pol polyprotein|reverse transcriptase|Tigger|Mariner|LINE|SINE|LTR|Retrovirus|Copia|Gypsy|Helitron|integrase'
te_count = meth['description'].str.contains(te_pattern, case=False, na=False).sum()
print(f"   TEs in methylated: {te_count}/{len(meth)} = {te_count/len(meth)*100:.1f}%")

print(f"\n   Saved to: {output}")
