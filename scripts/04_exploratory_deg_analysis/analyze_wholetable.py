#!/usr/bin/env python3
import pandas as pd
import re

deg = pd.read_csv('Results/05_transcriptomics/gene_functions_with_regulation_wholetable.csv')
meth = pd.read_csv('Results/05_transcriptomics/deg_direction_methylation.csv')

# Map transcript IDs to chacei genes
brooksi_map = pd.read_csv('Results/05_transcriptomics/brooksi_deg_to_chacei.tsv', sep='\t', names=['transcript_id','chacei_gene','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
deg_mapped = deg.merge(brooksi_map[['transcript_id','chacei_gene']], on='transcript_id', how='left')
deg_mapped = deg_mapped.merge(meth[['chacei_gene','methylation_status']], on='chacei_gene', how='left')
deg_mapped['methylation_status'] = deg_mapped['methylation_status'].fillna('not_mapped')

in_data = deg_mapped[deg_mapped['methylation_status'].isin(['methylated','unmethylated'])]

print("="*80)
print("DEG METHYLATION ANALYSIS (8,022 transcripts)")
print("="*80)
print(f"Total transcripts: {len(deg)}")
print(f"Mapped to S. chacei: {(deg_mapped['chacei_gene'].notna()).sum()}")
print(f"In methylation data: {len(in_data)}")

print("\n" + "="*80)
print("1. METHYLATION BY DIRECTION")
print("="*80)
ct = pd.crosstab(in_data['regulation_direction'], in_data['methylation_status'], margins=True)
print(ct)
q_meth = len(in_data[(in_data['regulation_direction']=='Queen-upregulated') & (in_data['methylation_status']=='methylated')])
q_tot = len(in_data[in_data['regulation_direction']=='Queen-upregulated'])
w_meth = len(in_data[(in_data['regulation_direction']=='Worker-upregulated') & (in_data['methylation_status']=='methylated')])
w_tot = len(in_data[in_data['regulation_direction']=='Worker-upregulated'])
print(f"\nQueen-up methylated: {q_meth}/{q_tot} = {q_meth/q_tot*100:.1f}%")
print(f"Worker-up methylated: {w_meth}/{w_tot} = {w_meth/w_tot*100:.1f}%")

print("\n" + "="*80)
print("2. TE ENRICHMENT")
print("="*80)
te_ct = pd.crosstab(in_data['is_te'], in_data['methylation_status'], margins=True)
print(te_ct)
te_meth = len(in_data[(in_data['is_te']==True) & (in_data['methylation_status']=='methylated')])
te_tot = len(in_data[in_data['is_te']==True])
if te_tot > 0:
    print(f"\nTEs methylated: {te_meth}/{te_tot} = {te_meth/te_tot*100:.1f}%")
else:
    print(f"\nNo TEs found in mapped DEGs")

print("\n" + "="*80)
print("3. METHYLATED DEGs")
print("="*80)
meth_degs = in_data[in_data['methylation_status']=='methylated'].sort_values('regulation_direction')
print(f"Total: {len(meth_degs)}")
for _, r in meth_degs.iterrows():
    te = '[TE]' if r['is_te'] else ''
    print(f"  {r['transcript_id']} | {r['regulation_direction']} | {r['description'][:60]} {te}")

deg_mapped.to_csv('Results/06_deg_methylation/wholetable_deg_methylation_analysis.csv', index=False)
print(f"\n✓ Exported to Results/06_deg_methylation/wholetable_deg_methylation_analysis.csv")
