#!/usr/bin/env python3
import pandas as pd
import re

# Load data
deg_meth = pd.read_csv('Results/05_transcriptomics/deg_direction_methylation.csv')
deg_annot = pd.read_csv('Results/05_transcriptomics/deg_uniprot_annotations_final.csv')

# Merge
deg_meth['sequence_id'] = deg_meth['chacei_gene'].str.replace('_ID=', '_') + '.t1'
merged = deg_meth.merge(deg_annot[['sequence_id', 'description']], on='sequence_id', how='left')

# Flag TEs
te_pattern = r'[Tt]ranspos|[Rr]etrotrans|Gag|Pol polyprotein|reverse transcriptase|Tigger|Mariner|LINE|SINE|LTR|Retrovirus|Copia|Gypsy|Helitron|integrase|ZN566|C2H2-type'
merged['is_te'] = merged['description'].fillna('').str.contains(te_pattern, regex=True)

# Filter to genes in methylation data
in_data = merged[merged['methylation_status'] != 'not_in_data'].copy()

print("=" * 80)
print("DEG METHYLATION ANALYSIS (136 DEGs from S. elizabethae)")
print("=" * 80)

print(f"\nTotal DEGs: {len(deg_meth)}")
print(f"DEGs in S. chacei methylation data: {len(in_data)} ({len(in_data)/len(deg_meth)*100:.1f}%)")
print(f"DEGs with functional annotation: {merged['description'].notna().sum()}")

print("\n" + "=" * 80)
print("1. METHYLATION BY DIRECTION")
print("=" * 80)
ct = pd.crosstab(in_data['direction'], in_data['methylation_status'], margins=True)
print(ct)
print(f"\nQueen-up methylated: {len(in_data[(in_data['direction']=='Queen_up') & (in_data['methylation_status']=='methylated')])}/{len(in_data[in_data['direction']=='Queen_up'])} = {len(in_data[(in_data['direction']=='Queen_up') & (in_data['methylation_status']=='methylated')])/len(in_data[in_data['direction']=='Queen_up'])*100:.1f}%")
print(f"Worker-up methylated: {len(in_data[(in_data['direction']=='Worker_up') & (in_data['methylation_status']=='methylated')])}/{len(in_data[in_data['direction']=='Worker_up'])} = {len(in_data[(in_data['direction']=='Worker_up') & (in_data['methylation_status']=='methylated')])/len(in_data[in_data['direction']=='Worker_up'])*100:.1f}%")

print("\n" + "=" * 80)
print("2. TE ENRICHMENT IN DEGs")
print("=" * 80)
annot = merged[merged['description'].notna()]
print(f"Annotated DEGs: {len(annot)}")
print(f"TEs in annotated DEGs: {annot['is_te'].sum()} ({annot['is_te'].sum()/len(annot)*100:.1f}%)")
print(f"\nTE DEGs by methylation:")
te_meth = annot[annot['is_te']].groupby('methylation_status').size()
print(te_meth)

print("\n" + "=" * 80)
print("3. METHYLATED DEGs")
print("=" * 80)
meth_degs = in_data[in_data['methylation_status'] == 'methylated']
print(f"Total methylated DEGs: {len(meth_degs)}")
for _, row in meth_degs.iterrows():
    desc = row['description'] if pd.notna(row['description']) else 'No annotation'
    te_flag = ' [TE]' if row['is_te'] else ''
    print(f"  {row['chacei_gene']} | {row['direction']} | {desc}{te_flag}")

print("\n" + "=" * 80)
print("4. FUNCTIONAL COMPARISON")
print("=" * 80)
print("\nMethylated DEG functions:")
meth_annot = merged[(merged['methylation_status']=='methylated') & (merged['description'].notna())]
for desc in meth_annot['description'].unique():
    print(f"  - {desc}")

print("\nUnmethylated DEG functions (sample):")
unmeth_annot = merged[(merged['methylation_status']=='unmethylated') & (merged['description'].notna())]
for desc in unmeth_annot['description'].head(10).unique():
    print(f"  - {desc}")

# Export
merged.to_csv('Results/06_deg_methylation/deg_methylation_functional_analysis.csv', index=False)
print(f"\n✓ Exported to Results/06_deg_methylation/deg_methylation_functional_analysis.csv")
