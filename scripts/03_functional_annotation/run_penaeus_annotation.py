#!/usr/bin/env python3
"""
Run Diamond BLASTP against Penaeus vannamei proteome.
"""

import subprocess
import pandas as pd
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DB_DIR = BASE_DIR / "Data" / "databases"
QUERY_FILE = BASE_DIR / "Data" / "unannotated_proteins" / "unannotated_proteins.fasta"
OUTPUT_DIR = BASE_DIR / "Results" / "04_functional_annotation"
OUTPUT_FILE = OUTPUT_DIR / "penaeus_results.tsv"

def run_diamond_blast():
    """Run Diamond BLASTP against Penaeus database."""
    print("=" * 60)
    print("RUNNING DIAMOND BLASTP - PENAEUS VANNAMEI")
    print("=" * 60)
    
    cmd = [
        "diamond", "blastp",
        "--db", str(DB_DIR / "penaeus_db.dmnd"),
        "--query", str(QUERY_FILE),
        "--out", str(OUTPUT_FILE),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore", "stitle",
        "--evalue", "1e-3",        # Relaxed for cross-species
        "--query-cover", "50",     # 50% query coverage
        "--subject-cover", "50",   # 50% subject coverage
        "--id", "30",              # 30% minimum identity
        "--max-target-seqs", "1",
        "--threads", "8",
        "--block-size", "2.0",
        "--index-chunks", "1"
    ]
    
    print(f"\nCommand: {' '.join(cmd)}\n")
    print("Running Diamond BLAST (this may take 15-30 minutes)...")
    print("Progress will be shown below:\n")
    
    subprocess.run(cmd, check=True)
    
    print("\n" + "=" * 60)
    print("DIAMOND BLAST COMPLETE")
    print("=" * 60)

def parse_results():
    """Parse Diamond results and report statistics."""
    df = pd.read_csv(OUTPUT_FILE, sep='\t', header=None,
                     names=['protein_id', 'uniprot_id', 'pident', 'evalue', 'bitscore', 'description'])
    
    print(f"\nResults saved to: {OUTPUT_FILE}")
    print(f"\nNew annotations found: {len(df)}")
    print(f"\nTop 10 hits:")
    print(df[['protein_id', 'pident', 'evalue', 'description']].head(10).to_string(index=False))
    
    print(f"\n\nIdentity distribution:")
    print(df['pident'].describe())
    
    return len(df)

def main():
    # Check if query file exists
    if not QUERY_FILE.exists():
        print(f"ERROR: Query file not found: {QUERY_FILE}")
        print("Please run extract_unannotated_proteins.py first")
        return
    
    # Check if database exists
    if not (DB_DIR / "penaeus_db.dmnd").exists():
        print(f"ERROR: Penaeus database not found")
        print("Please run setup_penaeus_db.sh first")
        return
    
    # Run Diamond BLAST
    run_diamond_blast()
    
    # Parse and report results
    new_annotations = parse_results()
    
    print("\n" + "=" * 60)
    print(f"PHASE 1 COMPLETE: +{new_annotations} new annotations")
    print("=" * 60)

if __name__ == "__main__":
    main()
