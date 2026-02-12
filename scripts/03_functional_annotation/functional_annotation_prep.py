#!/usr/bin/env python3
"""
Functional Annotation Preparation
Extracts gene IDs from filtered promoter-CpG pairs and prepares for BLASTP analysis.
"""

import re
import subprocess
from pathlib import Path
import pandas as pd


def extract_gene_ids(filtered_csv, output_txt):
    """Extract unique gene IDs from filtered results."""
    df = pd.read_csv(filtered_csv)
    df['gene_id_clean'] = df['gene_id'].str.replace("ID=", "", regex=False)
    unique_genes = df['gene_id_clean'].dropna().unique()
    
    with open(output_txt, "w") as f:
        for gene_id in unique_genes:
            f.write(f"{gene_id}\n")
    
    print(f"Gene ID list saved to: {output_txt}")
    return len(unique_genes)


def count_unique_genes(parsed_csv):
    """Count unique genes near CpG islands."""
    df = pd.read_csv(parsed_csv)
    gene_ids = df['gene_id'].unique()
    count = len(gene_ids)
    print(f"Total unique genes near CpG islands: {count}")
    return count


def merge_gff_files(augustus_dir, output_gff):
    """Merge all AUGUSTUS GFF files into one."""
    augustus_path = Path(augustus_dir)
    gff_files = list(augustus_path.glob("contig_*_augustus.gff"))
    
    with open(output_gff, "w") as outfile:
        for gff_file in gff_files:
            with open(gff_file) as infile:
                outfile.write(infile.read())
    
    print(f"Merged {len(gff_files)} GFF files to: {output_gff}")


def extract_gff_gene_ids(gff_path):
    """Extract gene IDs from merged GFF file."""
    gene_ids = set()
    
    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("# start gene"):
                match = re.search(r"# start gene (\S+)", line)
                if match:
                    gene_ids.add(match.group(1))
    
    print(f"Sample gene IDs from GFF: {list(gene_ids)[:10]}")
    print(f"Total unique gene IDs in GFF: {len(gene_ids)}")
    return gene_ids


def install_gffread():
    """Install gffread if not available."""
    if subprocess.run(["which", "gffread"], capture_output=True).returncode != 0:
        print("gffread not found. Installing...")
        try:
            subprocess.run(["brew", "install", "gffread"], check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Failed to install gffread. Install manually: brew install gffread")
            return False
    return True


def main():
    script_dir = Path(__file__).parent.parent
    results_dir = script_dir / "results (1)"
    
    filtered_csv = results_dir / "filtered_promoter_cpg_pairs.csv"
    parsed_csv = results_dir / "promoter_cpg_parsed.csv"
    blastp_output = results_dir / "filtered_gene_ids_for_blastp.txt"
    augustus_dir = results_dir / "augustus_output"
    merged_gff = results_dir / "merged_augustus_predictions.gff"
    
    print("=== Functional Annotation Preparation ===\n")
    
    print("Step 1: Extracting gene IDs for BLASTP...")
    gene_count = extract_gene_ids(filtered_csv, blastp_output)
    
    print(f"\nStep 2: Counting unique genes...")
    count_unique_genes(parsed_csv)
    
    print("\nStep 3: Installing gffread...")
    if not install_gffread():
        print("Skipping GFF processing")
        return
    
    print("\nStep 4: Merging AUGUSTUS GFF files...")
    merge_gff_files(augustus_dir, merged_gff)
    
    print("\nStep 5: Extracting gene IDs from GFF...")
    extract_gff_gene_ids(merged_gff)
    
    print("\n=== Preparation Complete ===")


if __name__ == "__main__":
    main()
