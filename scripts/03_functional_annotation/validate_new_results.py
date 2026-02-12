#!/usr/bin/env python3
"""Validate integrity of new DIAMOND results"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]

# Load data
old = pd.read_csv(BASE / "Results/04_functional_annotation/gene_uniprot_mapping.csv")
new = pd.read_csv(BASE / "Results/04_functional_annotation/gene_uniprot_mapping_NEW.csv")
meth = pd.read_csv(BASE / "Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv", sep="\t")
unmeth = pd.read_csv(BASE / "Results/03_methylation_analysis/genes_without_methylated_promoters.tsv", sep="\t")

print("=" * 70)
print("DATA INTEGRITY VALIDATION")
print("=" * 70)

# 1. Check total gene counts
print("\n1. GENE COUNT VERIFICATION")
print(f"   Methylated genes (source): {len(meth)}")
print(f"   Unmethylated genes (source): {len(unmeth)}")
print(f"   Total: {len(meth) + len(unmeth)}")

# 2. Compare old vs new
print("\n2. OLD vs NEW COMPARISON")
print(f"   Old total hits: {len(old)}")
print(f"   New total hits: {len(new)}")
print(f"   Old methylated: {(old['methylation_status']=='methylated').sum()}")
print(f"   New methylated: {(new['methylation_status']=='methylated').sum()}")

# 3. Check for duplicates
print("\n3. DUPLICATE CHECK")
print(f"   Old unique genes: {old['gene_id'].nunique()} / {len(old)}")
print(f"   New unique genes: {new['gene_id'].nunique()} / {len(new)}")
if new['gene_id'].nunique() != len(new):
    print(f"   ⚠️  WARNING: {len(new) - new['gene_id'].nunique()} duplicate genes!")

# 4. Identity distribution
print("\n4. ALIGNMENT QUALITY (% identity)")
print(f"   Old mean: {old['pident'].mean():.1f}%")
print(f"   New mean: {new['pident'].mean():.1f}%")
print(f"   Old median: {old['pident'].median():.1f}%")
print(f"   New median: {new['pident'].median():.1f}%")
print(f"   New min: {new['pident'].min():.1f}%")
print(f"   New <40% identity: {(new['pident'] < 40).sum()} ({(new['pident'] < 40).sum()/len(new)*100:.1f}%)")

# 5. E-value distribution
print("\n5. E-VALUE DISTRIBUTION")
print(f"   Old mean: {old['evalue'].mean():.2e}")
print(f"   New mean: {new['evalue'].mean():.2e}")
print(f"   New >1e-5: {(new['evalue'] > 1e-5).sum()} ({(new['evalue'] > 1e-5).sum()/len(new)*100:.1f}%)")
print(f"   New >1e-3: {(new['evalue'] > 1e-3).sum()}")

# 6. TE content check
te_pattern = r'[Tt]ranspos|[Rr]etrotrans|Gag|Pol polyprotein|reverse transcriptase|Tigger|Mariner|LINE|SINE|LTR|Retrovirus|Copia|Gypsy|Helitron|integrase'
old_meth = old[old['methylation_status'] == 'methylated']
new_meth = new[new['methylation_status'] == 'methylated']

old_te = old_meth['description'].str.contains(te_pattern, case=False, na=False).sum()
new_te = new_meth['description'].str.contains(te_pattern, case=False, na=False).sum()

print("\n6. TE ENRICHMENT")
print(f"   Old: {old_te}/{len(old_meth)} = {old_te/len(old_meth)*100:.1f}%")
print(f"   New: {new_te}/{len(new_meth)} = {new_te/len(new_meth)*100:.1f}%")

# 7. Sample annotations
print("\n7. SAMPLE NEW ANNOTATIONS (first 10 methylated)")
sample = new_meth.head(10)[['gene_id', 'pident', 'evalue', 'description']]
for _, row in sample.iterrows():
    is_te = "🔴 TE" if pd.notna(row['description']) and pd.Series([row['description']]).str.contains(te_pattern, case=False).iloc[0] else "🔵 Func"
    print(f"   {is_te} | {row['pident']:.0f}% | {row['evalue']:.1e} | {row['description'][:50]}")

# 8. Check for suspicious patterns
print("\n8. QUALITY FLAGS")
low_id = (new['pident'] < 30).sum()
high_eval = (new['evalue'] > 1e-3).sum()
print(f"   ⚠️  Alignments <30% identity: {low_id}")
print(f"   ⚠️  E-values >1e-3: {high_eval}")
if low_id > 0 or high_eval > 0:
    print("   ⚠️  Some hits may be false positives")

print("\n" + "=" * 70)
print("VALIDATION COMPLETE")
print("=" * 70)
