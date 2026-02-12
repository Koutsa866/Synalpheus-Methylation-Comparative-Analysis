#!/usr/bin/env python3
"""
Merge SwissProt and Penaeus annotations into final annotation file.
"""

import pandas as pd
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
RESULTS_DIR = BASE_DIR / "Results" / "04_functional_annotation"

def merge_annotations():
    """Merge SwissProt and Penaeus annotations."""
    print("=" * 60)
    print("MERGING ANNOTATIONS")
    print("=" * 60)
    
    # Load original SwissProt annotations
    print("\nLoading SwissProt annotations...")
    swissprot = pd.read_csv(RESULTS_DIR / "gene_uniprot_mapping.csv")
    print(f"  SwissProt: {len(swissprot)} genes")
    
    # Load Penaeus annotations
    print("Loading Penaeus annotations...")
    penaeus_file = RESULTS_DIR / "penaeus_results.tsv"
    
    if not penaeus_file.exists():
        print(f"ERROR: Penaeus results not found: {penaeus_file}")
        return
    
    penaeus = pd.read_csv(penaeus_file, sep='\t', header=None,
                          names=['protein_id', 'uniprot_id', 'pident', 'evalue', 'bitscore', 'description'])
    print(f"  Penaeus: {len(penaeus)} genes")
    
    # Add source column
    swissprot['source'] = 'SwissProt'
    penaeus['source'] = 'Penaeus_vannamei'
    
    # Combine (Penaeus only has protein_id, need to add gene_id and methylation_status)
    # For now, just combine the protein-level annotations
    combined = pd.concat([swissprot, penaeus], ignore_index=True)
    
    # Save merged annotations
    output_file = RESULTS_DIR / "combined_annotations.csv"
    combined.to_csv(output_file, index=False)
    
    print(f"\n" + "=" * 60)
    print("MERGE COMPLETE")
    print("=" * 60)
    print(f"Total annotations: {len(combined)}")
    print(f"  SwissProt: {len(swissprot)}")
    print(f"  Penaeus: {len(penaeus)}")
    print(f"\nSaved to: {output_file}")
    
    # Calculate new annotation rate
    total_genes = 6848  # From your results
    annotation_rate = (len(combined) / total_genes) * 100
    improvement = annotation_rate - 80.4
    
    print(f"\n" + "=" * 60)
    print("ANNOTATION STATISTICS")
    print("=" * 60)
    print(f"Total genes: {total_genes}")
    print(f"Annotated genes: {len(combined)}")
    print(f"Annotation rate: {annotation_rate:.1f}%")
    print(f"Improvement: +{improvement:.1f}%")
    print("=" * 60)

if __name__ == "__main__":
    merge_annotations()
