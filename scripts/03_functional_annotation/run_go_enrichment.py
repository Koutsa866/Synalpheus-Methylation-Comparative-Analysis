#!/usr/bin/env python3
"""GO Enrichment Analysis - Methylated vs Unmethylated Genes"""

import pandas as pd
import requests
from scipy.stats import fisher_exact
from pathlib import Path
import time

BASE = Path(__file__).parent
OUTPUT = BASE / "analysis_output"

# Load gene mapping
print("Loading gene annotations...")
df = pd.read_csv(OUTPUT / "gene_uniprot_mapping.csv")

meth = df[df['methylation_status'] == 'methylated']['uniprot_id'].str.split('|').str[1].unique()
unmeth = df[df['methylation_status'] == 'unmethylated']['uniprot_id'].str.split('|').str[1].unique()

print(f"Methylated genes: {len(meth)}")
print(f"Unmethylated genes: {len(unmeth)}")

# Fetch GO terms from UniProt
print("\nFetching GO annotations from UniProt...")
def get_go_terms(uniprot_ids, batch_size=500):
    go_map = {}
    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i+batch_size]
        query = " OR ".join([f"accession:{uid}" for uid in batch])
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,go_p,go_f,go_c&format=tsv"
        
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')[1:]
                for line in lines:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        acc = parts[0]
                        go_terms = []
                        for col in parts[1:]:
                            if col:
                                go_terms.extend([t.split(' [')[0] for t in col.split('; ') if t])
                        go_map[acc] = go_terms
            time.sleep(0.5)
        except Exception as e:
            print(f"  Error fetching batch: {e}")
    return go_map

meth_go = get_go_terms(list(meth))
unmeth_go = get_go_terms(list(unmeth))

# Build GO term counts
print("\nCounting GO terms...")
all_go_terms = set()
for terms in list(meth_go.values()) + list(unmeth_go.values()):
    all_go_terms.update(terms)

print(f"Total unique GO terms: {len(all_go_terms)}")

# Fisher's exact test for each GO term
print("\nRunning enrichment analysis...")
results = []
for go_term in all_go_terms:
    meth_with = sum(1 for g in meth if g in meth_go and go_term in meth_go[g])
    meth_without = len(meth) - meth_with
    unmeth_with = sum(1 for g in unmeth if g in unmeth_go and go_term in unmeth_go[g])
    unmeth_without = len(unmeth) - unmeth_with
    
    if meth_with > 0 or unmeth_with > 0:
        table = [[meth_with, meth_without], [unmeth_with, unmeth_without]]
        odds_ratio, p_value = fisher_exact(table)
        
        results.append({
            'go_term': go_term,
            'methylated_count': meth_with,
            'unmethylated_count': unmeth_with,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'enriched_in': 'methylated' if odds_ratio > 1 else 'unmethylated'
        })

# Save results
results_df = pd.DataFrame(results).sort_values('p_value')
results_df['significant'] = results_df['p_value'] < 0.05

output_file = OUTPUT / "go_enrichment_results.csv"
results_df.to_csv(output_file, index=False)

# Summary
sig = results_df[results_df['significant']]
print(f"\n✅ COMPLETE!")
print(f"  Results: {output_file}")
print(f"  Total GO terms tested: {len(results_df)}")
print(f"  Significant (p<0.05): {len(sig)}")
print(f"  Enriched in methylated: {len(sig[sig['enriched_in']=='methylated'])}")
print(f"  Enriched in unmethylated: {len(sig[sig['enriched_in']=='unmethylated'])}")

# Top 10 enriched in methylated genes
print("\nTop 10 GO terms enriched in METHYLATED genes:")
top_meth = sig[sig['enriched_in']=='methylated'].head(10)
for _, row in top_meth.iterrows():
    print(f"  {row['go_term']} (p={row['p_value']:.2e}, OR={row['odds_ratio']:.2f})")
