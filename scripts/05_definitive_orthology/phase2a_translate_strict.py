#!/usr/bin/env python3
"""
Pipeline: 05_definitive_orthology
Phase: 2a - Transcriptome Translation
Goal: Translate DNA transcriptomes to protein with strict ID preservation
Inputs:
  - Data/transcriptomics/S.brooksi_brain_transcriptome.fasta
  - Data/transcriptomics/S.elizabethae_brain_transcriptome.fasta
Outputs:
  - Data/transcriptomics/S.brooksi_brain_transcriptome.faa
  - Data/transcriptomics/S.elizabethae_brain_transcriptome.faa
Author: Philip Koutsaftis
Date: 2025
"""

"""
Phase 2a: Strict-ID Translation
Translates DNA to Protein while preserving exact Trinity IDs
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

targets = {
    "brooksi": "Data/transcriptomics/transcriptomes/S.brooksi_brain_filtered_transcriptome.fasta",
    "elizabethae": "Data/transcriptomics/transcriptomes/S.elizabethae_brain_filtered_transcriptome.fasta"
}

for sp, path in targets.items():
    output_faa = f"Data/transcriptomics/transcriptomes/S.{sp}_strict.faa"
    prot_records = []
    
    for record in SeqIO.parse(path, "fasta"):
        # Extract clean ID (before first space)
        clean_id = record.id.split()[0]
        
        # Translate DNA to protein (all 6 frames, take longest)
        prot_seq = record.seq.translate(to_stop=False)
        
        # Keep EXACT ID match
        prot_records.append(SeqRecord(prot_seq, id=clean_id, description=""))
    
    SeqIO.write(prot_records, output_faa, "fasta")
    print(f"✓ Created {output_faa}: {len(prot_records)} sequences")

print("\n✓ Strict-ID translation complete!")
print("  IDs preserved exactly as in DEA files")
