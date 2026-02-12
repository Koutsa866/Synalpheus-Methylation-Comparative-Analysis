#!/usr/bin/env python3
"""
Pipeline: 05_definitive_orthology
Phase: 1 - Gene Body Extraction
Goal: Extract gene body sequences (start-to-stop, NO promoters) with methylation labels
Note: Implements Dr. Chak's requirement for gene-body extraction (not promoters)
Inputs: 
  - Data/genome/S_chacei_genome.fasta
  - Data/genome/S_chacei_genes.gff3
  - Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv
  - Results/03_methylation_analysis/genes_without_methylated_promoters.tsv
Outputs:
  - Results/05_ortholog_analysis/S_chacei_gene_bodies.fasta
Author: Philip Koutsaftis
Date: 2025
"""

"""
Phase 1: Extract S. chacei Gene Body Sequences with Methylation Labels
Goal: Extract Gene Body (Start-to-Stop) sequences and label headers with 
methylation status for downstream statistical comparison.
"""

from Bio import SeqIO
import sys
import os

# --- CONFIGURATION ---
GENOME = "Data/genome/assembly.fasta"
GFF = "Results/01_gene_prediction/augustus_complete_merged.gff"
METHYLATED = "Results/03_methylation_analysis/genes_with_high_methylated_promoters.tsv"
UNMETHYLATED = "Results/03_methylation_analysis/genes_without_methylated_promoters.tsv"
OUTPUT = "Results/05_ortholog_analysis/S_chacei_gene_bodies.fasta"

def load_methylation_status():
    """
    Load both Methylated and Unmethylated IDs into a single dictionary.
    This creates the 'Background' necessary for valid statistics.
    """
    status_map = {}
    print("Loading methylation data...")
    
    # Helper to process files
    def process_file(filepath, label, col=0):
        count = 0
        if not os.path.exists(filepath):
            print(f"Warning: File not found: {filepath}")
            return 0
            
        with open(filepath, 'r') as f:
            try:
                next(f) # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    # Safety check: Ensure the column index exists for this line
                    if len(parts) > col:
                        gene_id = parts[col]
                        status_map[gene_id] = label
                        count += 1
            except StopIteration:
                pass # Empty file
        return count

    m_count = process_file(METHYLATED, "Methylated")
    u_count = process_file(UNMETHYLATED, "Unmethylated", col=3)  # Unmethylated uses column 3
    
    print(f"  Loaded {m_count} Methylated and {u_count} Unmethylated genes.")
    return status_map

def extract_sequences():
    """Main extraction logic using memory-efficient indexing"""
    # 1. Load Methylation Data
    meth_status = load_methylation_status()
    if not meth_status:
        print("Error: No methylation data loaded (or files are empty). Exiting.")
        sys.exit(1)

    # 2. Index Genome (Memory Safe)
    # SeqIO.index creates a SQLite-like index; it doesn't load sequences into RAM.
    print(f"Indexing genome: {GENOME} (Memory-safe mode)...")
    if not os.path.exists(GENOME):
        print(f"Error: Genome file not found at {GENOME}")
        sys.exit(1)
        
    try:
        genome_idx = SeqIO.index(GENOME, "fasta")
    except ValueError as e:
        print(f"Error indexing FASTA: {e}")
        sys.exit(1)

    # 3. Parse GFF and Extract
    print(f"Parsing GFF and extracting gene bodies to {OUTPUT}...")
    
    count_extracted = 0
    count_skipped = 0
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    
    with open(GFF, 'r') as gff_handle, open(OUTPUT, 'w') as out_handle:
        for line in gff_handle:
            if line.startswith("#") or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            
            # GOLDEN RULE: Extract 'gene' features only (Start-to-Stop / Gene body)
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            
            # Define chrom FIRST to avoid NameError when building gene_id
            chrom = parts[0]
            start = int(parts[3]) - 1 # Convert GFF 1-based to Python 0-based
            end = int(parts[4])
            strand = parts[6]
            
            # Robust ID Parsing from Column 9
            attributes_str = parts[8]
            short_id = None
            for attr in attributes_str.split(';'):
                if attr.strip().startswith("ID="):
                    short_id = attr.strip().split("=")[1]
                    break
            
            if not short_id:
                continue

            # Construct the Long ID: 'contig_NAME_ID=gX'
            gene_id = f"{chrom}_ID={short_id}"
            
            # Check if this gene is in our study (Methylated OR Unmethylated)
            if gene_id not in meth_status:
                continue
            
            # Extract sequence from indexed genome
            if chrom in genome_idx:
                # Fetch sequence directly from disk (via index)
                full_chrom_seq = genome_idx[chrom].seq
                gene_seq = full_chrom_seq[start:end]
                
                # Strand Awareness: Reverse Complement if the gene is on '-' strand
                if strand == '-':
                    gene_seq = gene_seq.reverse_complement()
                
                # Write to Output with Label
                # Header: >contig_1_ID=g1|Methylated
                label = meth_status[gene_id]
                header = f">{gene_id}|{label}"
                
                out_handle.write(f"{header}\n{str(gene_seq)}\n")
                count_extracted += 1
            else:
                count_skipped += 1

    print(f"Success.")
    print(f"  Total Extracted: {count_extracted}")
    print(f"  Missing from Genome: {count_skipped}")

if __name__ == "__main__":
    extract_sequences()
