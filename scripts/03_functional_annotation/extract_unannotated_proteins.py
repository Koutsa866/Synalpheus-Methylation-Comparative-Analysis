#!/usr/bin/env python3
"""
Extract protein sequences for genes without functional annotations.
"""

import os
from pathlib import Path
from Bio import SeqIO

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
RESULTS_DIR = BASE_DIR / "Results" / "04_functional_annotation"
AUGUSTUS_DIR = BASE_DIR / "Results" / "01_gene_prediction" / "augustus_output"
OUTPUT_DIR = BASE_DIR / "Data" / "unannotated_proteins"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def get_annotated_proteins():
    """Get set of protein IDs that already have annotations."""
    annotated = set()
    mapping_file = RESULTS_DIR / "gene_uniprot_mapping.csv"
    
    with open(mapping_file) as f:
        next(f)  # Skip header
        for line in f:
            protein_id = line.split(',')[0]
            annotated.add(protein_id)
    
    print(f"Found {len(annotated)} annotated proteins")
    return annotated

def extract_all_proteins():
    """Extract all protein sequences from AUGUSTUS GFF3 files."""
    all_proteins = {}
    
    # Search for all GFF3 files in augustus_output subdirectories
    for gff_file in AUGUSTUS_DIR.rglob("*.gff3"):
        print(f"Processing {gff_file.name}...")
        
        with open(gff_file) as f:
            current_protein = None
            sequence_lines = []
            
            for line in f:
                if line.startswith("# protein sequence = ["):
                    # Start of protein sequence
                    sequence_lines = [line.split("[")[1].rstrip("]\n")]
                elif line.startswith("# ") and current_protein and sequence_lines:
                    # Continuation of sequence
                    sequence_lines.append(line[2:].rstrip("]\n"))
                elif line.startswith("# end gene"):
                    # End of gene - save protein
                    if current_protein and sequence_lines:
                        seq = ''.join(sequence_lines).replace(' ', '')
                        all_proteins[current_protein] = seq
                        current_protein = None
                        sequence_lines = []
                elif not line.startswith("#"):
                    # GFF3 feature line
                    fields = line.strip().split('\t')
                    if len(fields) >= 9 and fields[2] == "CDS":
                        # Extract protein ID from attributes
                        attrs = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                        if 'Parent' in attrs:
                            current_protein = attrs['Parent']
    
    print(f"Extracted {len(all_proteins)} total proteins")
    return all_proteins

def write_unannotated_fasta(all_proteins, annotated_proteins):
    """Write FASTA file of unannotated proteins."""
    output_file = OUTPUT_DIR / "unannotated_proteins.fasta"
    unannotated_count = 0
    
    with open(output_file, 'w') as f:
        for protein_id, sequence in all_proteins.items():
            if protein_id not in annotated_proteins:
                f.write(f">{protein_id}\n{sequence}\n")
                unannotated_count += 1
    
    print(f"\nWrote {unannotated_count} unannotated proteins to {output_file}")
    return unannotated_count

def main():
    print("=" * 60)
    print("EXTRACTING UNANNOTATED PROTEINS")
    print("=" * 60)
    
    # Get annotated proteins
    annotated = get_annotated_proteins()
    
    # Extract all proteins from AUGUSTUS
    all_proteins = extract_all_proteins()
    
    # Write unannotated proteins
    unannotated_count = write_unannotated_fasta(all_proteins, annotated)
    
    print("\n" + "=" * 60)
    print(f"SUMMARY:")
    print(f"  Total proteins: {len(all_proteins)}")
    print(f"  Annotated: {len(annotated)}")
    print(f"  Unannotated: {unannotated_count}")
    print(f"  Annotation rate: {len(annotated)/len(all_proteins)*100:.1f}%")
    print("=" * 60)

if __name__ == "__main__":
    main()
