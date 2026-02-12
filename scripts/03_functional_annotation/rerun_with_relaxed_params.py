#!/usr/bin/env python3
"""
Re-run DIAMOND annotation with relaxed parameters
Uses existing gene_uniprot_mapping.csv as template for correct gene IDs
"""

import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO

BASE = Path(__file__).resolve().parents[2]

# Load existing annotation to get correct gene→protein mapping
print("Loading existing annotation mapping...")
existing = pd.read_csv(BASE / "Results/04_functional_annotation/gene_uniprot_mapping.csv")
print(f"  Found {existing['gene_id'].nunique()} unique genes")
print(f"  Found {existing['protein_id'].nunique()} unique proteins")

# Get all gene IDs we need to annotate
meth = pd.read_csv(BASE / "Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv", sep="\t")
unmeth = pd.read_csv(BASE / "Results/03_methylation_analysis/genes_without_methylated_promoters.tsv", sep="\t")
all_genes = set(meth['gene_id']) | set(unmeth['gene_id'])
print(f"\nTotal genes to annotate: {len(all_genes)}")

# Create protein→gene mapping from existing data
protein_to_gene = dict(zip(existing['protein_id'], existing['gene_id']))
gene_to_meth_status = {}
for g in meth['gene_id']:
    gene_to_meth_status[g] = 'methylated'
for g in unmeth['gene_id']:
    gene_to_meth_status[g] = 'unmethylated'

# Find the protein FASTA file that was actually used
# Check multiple possible locations
protein_files = [
    BASE / "Data/proteins/predicted_proteins.faa",
    BASE / "Results/01_gene_prediction/predicted_proteins.faa",
    BASE / "Data/filtered_genes_proteins.faa"
]

protein_fasta = None
for pf in protein_files:
    if pf.exists():
        protein_fasta = pf
        break

if not protein_fasta:
    print("\n❌ ERROR: Cannot find protein FASTA file")
    print("Searched:")
    for pf in protein_files:
        print(f"  - {pf}")
    exit(1)

print(f"\nUsing protein file: {protein_fasta}")

# Extract proteins that match our genes
print("\nExtracting proteins for annotation...")
query_proteins = []
protein_count = 0

for record in SeqIO.parse(protein_fasta, "fasta"):
    if record.id in protein_to_gene:
        query_proteins.append(record)
        protein_count += 1

if len(query_proteins) == 0:
    print("❌ ERROR: No proteins matched! Protein IDs don't match existing annotation.")
    print(f"Sample protein IDs from FASTA: {[r.id for r in list(SeqIO.parse(protein_fasta, 'fasta'))[:5]]}")
    print(f"Sample protein IDs from annotation: {list(protein_to_gene.keys())[:5]}")
    exit(1)

# Save query proteins
query_faa = Path("/tmp/query_proteins_relaxed.faa")
SeqIO.write(query_proteins, query_faa, "fasta")
print(f"  Extracted {len(query_proteins)} proteins")

# Setup DIAMOND database
print("\nSetting up DIAMOND database...")
db_file = BASE / "uniprot_sprot.fasta"
if not db_file.exists():
    print(f"❌ ERROR: SwissProt database not found: {db_file}")
    exit(1)

dmnd_db = Path("/tmp/sprot_relaxed")
if not (Path("/tmp/sprot_relaxed.dmnd")).exists():
    subprocess.run(["diamond", "makedb", "--in", str(db_file), "--db", str(dmnd_db)], check=True)

# Run DIAMOND with RELAXED parameters
print("\nRunning DIAMOND BLASTP with relaxed parameters...")
print("  E-value: 1e-3")
print("  Query coverage: 50%")
print("  Subject coverage: 50%")
print("  Identity: 30%")

diamond_out = Path("/tmp/diamond_relaxed.tsv")
subprocess.run([
    "diamond", "blastp",
    "--query", str(query_faa),
    "--db", str(dmnd_db),
    "--out", str(diamond_out),
    "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore", "stitle",
    "--max-target-seqs", "1",
    "--evalue", "1e-3",
    "--query-cover", "50",
    "--subject-cover", "50",
    "--id", "30",
    "--threads", "4"
], check=True)

# Parse results
print("\nParsing results...")
df = pd.read_csv(diamond_out, sep="\t", names=["protein_id", "uniprot_id", "pident", "evalue", "bitscore", "description"])
df['gene_id'] = df['protein_id'].map(protein_to_gene)
df['methylation_status'] = df['gene_id'].map(gene_to_meth_status)

# Deduplicate by gene (keep best hit)
df = df.sort_values(['gene_id', 'evalue', 'bitscore'], ascending=[True, True, False])
df = df.drop_duplicates(subset='gene_id', keep='first')

# Save results
output = BASE / "Results/04_functional_annotation/gene_uniprot_mapping_RELAXED.csv"
df.to_csv(output, index=False)

# Report statistics
meth_genes = df[df['methylation_status'] == 'methylated']
unmeth_genes = df[df['methylation_status'] == 'unmethylated']

print(f"\n✅ COMPLETE!")
print(f"   Total genes annotated: {len(df)}")
print(f"   Methylated: {len(meth_genes)}")
print(f"   Unmethylated: {len(unmeth_genes)}")

# Check TE enrichment
te_pattern = r'[Tt]ranspos|[Rr]etrotrans|Gag|Pol polyprotein|reverse transcriptase|Tigger|Mariner|LINE|SINE|LTR|Retrovirus|Copia|Gypsy|Helitron|integrase'
te_count = meth_genes['description'].str.contains(te_pattern, case=False, na=False).sum()
print(f"   TEs in methylated: {te_count}/{len(meth_genes)} = {te_count/len(meth_genes)*100:.1f}%")

print(f"\n   Saved to: {output}")
print(f"\nCompare with original:")
print(f"   Original: {(existing['methylation_status']=='methylated').sum()} methylated genes")
print(f"   Relaxed:  {len(meth_genes)} methylated genes")
