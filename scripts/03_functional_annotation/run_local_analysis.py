#!/usr/bin/env python3
"""Local Diamond BLASTP + GO Enrichment Pipeline"""

import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE = Path(__file__).resolve().parents[2]
RESULTS = BASE / "Results" / "03_methylation_analysis"
RESULTS_OLD = BASE / "archive" / "old_results" / "results (1)"
OUTPUT = Path("/tmp/diamond_analysis")
OUTPUT.mkdir(exist_ok=True)
# Step 1: Extract gene IDs
print("Step 1: Extracting gene IDs...")
meth = pd.read_csv(RESULTS / "genes_with_high_methylated_promoters.tsv", sep="\t")
unmeth = pd.read_csv(RESULTS / "genes_without_methylated_promoters.tsv", sep="\t")

meth_genes = set(meth['gene_id'].str.strip())
unmeth_genes = set(unmeth['gene_id'].str.strip())

print(f"  Methylated: {len(meth_genes)}, Unmethylated: {len(unmeth_genes)}")

# Step 2: Extract target proteins from existing file
print("Step 2: Extracting protein sequences...")
all_proteins = RESULTS_OLD / "filtered_genes_proteins_unique.faa"
proteins = []
protein_to_gene = {}  # Map protein ID to full gene ID

for rec in SeqIO.parse(all_proteins, "fasta"):
    # Match protein ID (g1_1) to gene ID (contig_XXXXX_ID=g1)
    gene_num = rec.id.split("_")[0]  # g1, g2, etc.
    
    # Find matching gene in methylated or unmethylated lists
    matched_gene = None
    for gene_id in meth_genes | unmeth_genes:
        if f"ID={gene_num}" in gene_id:
            matched_gene = gene_id
            break
    
    if matched_gene:
        clean_seq = ''.join(c for c in str(rec.seq) if c.isalpha())
        if clean_seq:
            new_rec = SeqRecord(Seq(clean_seq), id=rec.id, description="")
            proteins.append(new_rec)
            protein_to_gene[rec.id] = matched_gene

query_faa = OUTPUT / "query_proteins.faa"
SeqIO.write(proteins, query_faa, "fasta")
print(f"  Extracted {len(proteins)} protein sequences")

# Step 3: Build Diamond database
print("Step 3: Building Diamond database...")
db_file = BASE / "uniprot_sprot.fasta"
dmnd_db = OUTPUT / "sprot"
if not (OUTPUT / "sprot.dmnd").exists():
    subprocess.run(["diamond", "makedb", "--in", str(db_file), "--db", str(dmnd_db)], check=True)

# Step 4: Run Diamond BLASTP
print("Step 4: Running Diamond BLASTP...")
diamond_out = OUTPUT / "diamond_results.tsv"
subprocess.run([
    "diamond", "blastp",
    "--query", str(query_faa),
    "--db", str(dmnd_db),
    "--out", str(diamond_out),
    "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore", "stitle",
    "--max-target-seqs", "1",
    "--evalue", "1e-3",        # Relaxed for cross-species
    "--query-cover", "50",     # 50% query coverage
    "--subject-cover", "50",   # 50% subject coverage
    "--id", "30",              # 30% minimum identity
    "--threads", "4"
], check=True)

# Step 5: Create gene-UniProt mapping
print("Step 5: Creating gene mapping...")
df = pd.read_csv(diamond_out, sep="\t", names=["protein_id", "uniprot_id", "pident", "evalue", "bitscore", "description"])
df['gene_id'] = df['protein_id'].map(protein_to_gene)
df['methylation_status'] = df['gene_id'].apply(lambda x: 'methylated' if x in meth_genes else 'unmethylated')

mapping_file = OUTPUT / "gene_uniprot_mapping.csv"
df.to_csv(mapping_file, index=False)

# Copy results back to project
final_output = BASE / "Results" / "04_functional_annotation" / "gene_uniprot_mapping_NEW.csv"
df.to_csv(final_output, index=False)

print(f"\n✅ COMPLETE!")
print(f"  Results: {final_output}")
print(f"  Total hits: {len(df)}")
print(f"  Methylated hits: {(df['methylation_status']=='methylated').sum()}")
print(f"  Unmethylated hits: {(df['methylation_status']=='unmethylated').sum()}")
