#!/usr/bin/env python3
"""
BLASTP Analysis Pipeline
Prepares protein sequences for BLAST analysis against SwissProt database.
"""

import subprocess
import sys
import tempfile
import shutil
from pathlib import Path
from Bio import SeqIO


def deduplicate_fasta(input_path, output_path):
    """Remove duplicate sequences from FASTA file."""
    seen_seqs = set()
    unique_records = []
    
    for record in SeqIO.parse(input_path, "fasta"):
        seq_str = str(record.seq)
        if seq_str not in seen_seqs:
            seen_seqs.add(seq_str)
            unique_records.append(record)
    
    SeqIO.write(unique_records, output_path, "fasta")
    print(f"Saved {len(unique_records)} unique sequences to: {output_path}")
    return len(unique_records)


def count_sequences(fasta_path):
    """Count total sequences in FASTA file."""
    count = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
    print(f"Total protein sequences: {count}")
    return count


def extract_sample(input_path, output_path, n=5):
    """Extract first n sequences from FASTA file."""
    records = list(SeqIO.parse(input_path, "fasta"))[:n]
    SeqIO.write(records, output_path, "fasta")
    print(f"Sample FASTA written: {output_path}")
    return len(records)


def setup_blast_db(db_url, db_dir, db_name="swissprot_db"):
    """Download and create BLAST database."""
    db_dir = Path(db_dir)
    db_dir.mkdir(exist_ok=True)
    
    if subprocess.run(["which", "makeblastdb"], capture_output=True).returncode != 0:
        print("makeblastdb not found. Installing BLAST+...")
        try:
            subprocess.run(["brew", "install", "blast"], check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Failed to install BLAST+. Install manually: brew install blast")
            return
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        fasta_gz = tmpdir / "uniprot_sprot.fasta.gz"
        fasta = tmpdir / "uniprot_sprot.fasta"
        
        print("Downloading SwissProt database...")
        subprocess.run(["curl", "-o", str(fasta_gz), db_url], check=True)
        
        print("Extracting database...")
        subprocess.run(["gunzip", "-f", str(fasta_gz)], check=True)
        
        print("Creating BLAST database...")
        subprocess.run(
            ["makeblastdb", "-in", "uniprot_sprot.fasta", "-dbtype", "prot", "-out", db_name],
            cwd=str(tmpdir),
            check=True
        )
        
        print("Moving database files...")
        for file in tmpdir.glob(f"{db_name}*"):
            shutil.move(str(file), str(db_dir / file.name))
        
        print(f"BLAST database created: {db_dir / db_name}")


def main():
    script_dir = Path(__file__).parent.parent
    base_dir = script_dir / "results (1)"
    
    input_fasta = base_dir / "filtered_genes_proteins_cleaned.faa"
    dedup_fasta = base_dir / "filtered_genes_proteins_cleaned_dedup.faa"
    sample_fasta = base_dir / "test_5_proteins.faa"
    
    print("=== BLASTP Analysis Pipeline ===\n")
    
    print("Step 1: Counting sequences...")
    count_sequences(input_fasta)
    
    print("\nStep 2: Deduplicating sequences...")
    deduplicate_fasta(input_fasta, dedup_fasta)
    
    print("\nStep 3: Extracting sample...")
    extract_sample(dedup_fasta, sample_fasta, n=5)
    
    print("\nStep 4: Setting up BLAST database...")
    db_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    db_dir = base_dir / "blast_db"
    setup_blast_db(db_url, db_dir)
    
    print("\n=== Pipeline Complete ===")


if __name__ == "__main__":
    main()
