#!/usr/bin/env python3
"""
CORRECTED: Re-run DIAMOND with relaxed parameters
Fixes:
1. Uses correct protein→gene mapping from existing results
2. Deduplicates by gene (not protein)
3. Validates all inputs before running
4. Uses /tmp to avoid Google Drive sync issues
"""

import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO

BASE = Path(__file__).resolve().parents[2]

print("=" * 70)
print("DIAMOND RE-ANNOTATION WITH RELAXED PARAMETERS")
print("=" * 70)

# STEP 1: Load and validate existing annotation
print("\n[1/7] Loading existing annotation...")
existing_file = BASE / "Results/04_functional_annotation/gene_uniprot_mapping.csv"
if not existing_file.exists():
    print(f"❌ ERROR: {existing_file} not found")
    exit(1)

existing = pd.read_csv(existing_file)
print(f"  ✓ Loaded {len(existing)} annotations")
print(f"  ✓ Unique genes: {existing['gene_id'].nunique()}")
print(f"  ✓ Unique proteins: {existing['protein_id'].nunique()}")

# STEP 2: Load methylation data
print("\n[2/7] Loading methylation data...")
meth = pd.read_csv(BASE / "Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv", sep="\t")
unmeth = pd.read_csv(BASE / "Results/03_methylation_analysis/genes_without_methylated_promoters.tsv", sep="\t")
print(f"  ✓ Methylated genes: {len(meth)}")
print(f"  ✓ Unmethylated genes: {len(unmeth)}")
print(f"  ✓ Total: {len(meth) + len(unmeth)}")

# Create gene→methylation mapping
gene_to_status = {}
for g in meth['gene_id']:
    gene_to_status[g] = 'methylated'
for g in unmeth['gene_id']:
    gene_to_status[g] = 'unmethylated'

# STEP 3: Build protein→gene mapping from existing results
print("\n[3/7] Building protein→gene mapping...")
protein_to_gene = dict(zip(existing['protein_id'], existing['gene_id']))
print(f"  ✓ Mapped {len(protein_to_gene)} proteins to genes")

# STEP 4: Find and validate protein FASTA
print("\n[4/7] Locating protein FASTA file...")
possible_locations = [
    BASE / "archive/old_results/results (1)/filtered_genes_proteins_unique.faa",
    BASE / "archive/old_results/results (1)/archive/filtered_genes_proteins.faa",
    BASE / "Data/proteins/predicted_proteins.faa",
    BASE / "Results/01_gene_prediction/predicted_proteins.faa"
]

protein_fasta = None
for loc in possible_locations:
    if loc.exists():
        # Verify it has the right proteins
        sample_ids = [rec.id for rec in list(SeqIO.parse(loc, "fasta"))[:10]]
        matches = sum(1 for pid in sample_ids if pid in protein_to_gene)
        if matches > 0:
            protein_fasta = loc
            print(f"  ✓ Found: {loc}")
            print(f"  ✓ Validated: {matches}/10 sample proteins match")
            break

if not protein_fasta:
    print("\n❌ ERROR: Cannot find valid protein FASTA file")
    print("Searched locations:")
    for loc in possible_locations:
        print(f"  - {loc} {'(exists)' if loc.exists() else '(not found)'}")
    print("\nProtein IDs needed (from annotation):")
    print(f"  {list(protein_to_gene.keys())[:5]}")
    exit(1)

# STEP 5: Extract proteins and deduplicate by gene
print("\n[5/7] Extracting proteins (keeping best per gene)...")
gene_to_best_protein = {}

for record in SeqIO.parse(protein_fasta, "fasta"):
    if record.id in protein_to_gene:
        gene_id = protein_to_gene[record.id]
        # Keep first protein per gene (assumes sorted by quality)
        if gene_id not in gene_to_best_protein:
            gene_to_best_protein[gene_id] = record

print(f"  ✓ Extracted {len(gene_to_best_protein)} unique genes")
print(f"  ✓ Methylated: {sum(1 for g in gene_to_best_protein if gene_to_status.get(g)=='methylated')}")
print(f"  ✓ Unmethylated: {sum(1 for g in gene_to_best_protein if gene_to_status.get(g)=='unmethylated')}")

# Save query proteins
query_faa = Path("/tmp/query_proteins_corrected.faa")
SeqIO.write(list(gene_to_best_protein.values()), query_faa, "fasta")

# STEP 6: Setup and run DIAMOND
print("\n[6/7] Running DIAMOND BLASTP...")
db_file = BASE / "uniprot_sprot.fasta"
if not db_file.exists():
    print(f"❌ ERROR: {db_file} not found")
    exit(1)

dmnd_db = Path("/tmp/sprot_corrected")
if not Path("/tmp/sprot_corrected.dmnd").exists():
    print("  Building database...")
    subprocess.run(["diamond", "makedb", "--in", str(db_file), "--db", str(dmnd_db)], 
                   check=True, capture_output=True)

print("  Parameters:")
print("    E-value: 1e-3 (relaxed)")
print("    Query coverage: 50%")
print("    Subject coverage: 50%")
print("    Min identity: 30%")

diamond_out = Path("/tmp/diamond_corrected.tsv")
result = subprocess.run([
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
], capture_output=True, text=True)

if result.returncode != 0:
    print(f"❌ DIAMOND failed: {result.stderr}")
    exit(1)

# STEP 7: Parse and validate results
print("\n[7/7] Processing results...")
df = pd.read_csv(diamond_out, sep="\t", 
                 names=["protein_id", "uniprot_id", "pident", "evalue", "bitscore", "description"])

# Map protein→gene→methylation
df['gene_id'] = df['protein_id'].map(protein_to_gene)
df['methylation_status'] = df['gene_id'].map(gene_to_status)

# Verify no duplicates
if df['gene_id'].nunique() != len(df):
    print(f"  ⚠️  WARNING: Found duplicates, deduplicating...")
    df = df.sort_values(['gene_id', 'evalue'], ascending=[True, True])
    df = df.drop_duplicates(subset='gene_id', keep='first')

# Save to project directory
output = BASE / "Results/04_functional_annotation/gene_uniprot_mapping_RELAXED.csv"
df.to_csv(output, index=False)

# Calculate statistics
meth_df = df[df['methylation_status'] == 'methylated']
unmeth_df = df[df['methylation_status'] == 'unmethylated']

te_pattern = r'[Tt]ranspos|[Rr]etrotrans|Gag|Pol polyprotein|reverse transcriptase|Tigger|Mariner|LINE|SINE|LTR|Retrovirus|Copia|Gypsy|Helitron|integrase'
te_count = meth_df['description'].str.contains(te_pattern, case=False, na=False).sum()

print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)
print(f"\nTotal genes annotated: {len(df)} / {len(gene_to_best_protein)}")
print(f"  Methylated: {len(meth_df)} / {len(meth)}")
print(f"  Unmethylated: {len(unmeth_df)} / {len(unmeth)}")
print(f"\nTE enrichment in methylated genes:")
print(f"  {te_count} / {len(meth_df)} = {te_count/len(meth_df)*100:.1f}%")
print(f"\nAlignment quality:")
print(f"  Mean identity: {df['pident'].mean():.1f}%")
print(f"  Median identity: {df['pident'].median():.1f}%")
print(f"  Mean e-value: {df['evalue'].mean():.2e}")

print(f"\n✅ Saved to: {output}")

# Compare with original
orig_meth = (existing['methylation_status']=='methylated').sum()
print(f"\nComparison with original (1e-5 parameters):")
print(f"  Original methylated: {orig_meth}")
print(f"  Relaxed methylated:  {len(meth_df)}")
print(f"  Difference: +{len(meth_df) - orig_meth} genes")

print("\n" + "=" * 70)
