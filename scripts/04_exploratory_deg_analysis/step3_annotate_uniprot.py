#!/usr/bin/env python3
from Bio import SeqIO
import subprocess
import os

print("Step 3: Annotating transcripts with UniProt")

# Translate transcripts to proteins
def translate_fasta(input_fasta, output_fasta, species):
    print(f"\nTranslating {species} transcripts...")
    with open(output_fasta, 'w') as out:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            best_prot = ""
            for frame in range(3):
                seq = record.seq[frame:]
                prot = seq.translate(to_stop=True)
                if len(prot) > len(best_prot):
                    best_prot = str(prot)
            for frame in range(3):
                seq = record.seq.reverse_complement()[frame:]
                prot = seq.translate(to_stop=True)
                if len(prot) > len(best_prot):
                    best_prot = str(prot)
            if best_prot:
                out.write(f">{record.id}\n{best_prot}\n")
    print(f"  ✓ Translated to {output_fasta}")

translate_fasta('Data/transcriptomics/transcriptomes/S.brooksi_brain_filtered_transcriptome.fasta',
                'Data/transcriptomics/transcriptomes/brooksi_proteins.faa', 'S. brooksi')

translate_fasta('Data/transcriptomics/transcriptomes/S.elizabethae_brain_filtered_transcriptome.fasta',
                'Data/transcriptomics/transcriptomes/elizabethae_proteins.faa', 'S. elizabethae')

# Make DIAMOND database
print("\nCreating DIAMOND database from UniProt...")
if not os.path.exists('Data/databases/uniprot_sprot.dmnd'):
    subprocess.run(['diamond', 'makedb', '--in', 'Data/databases/uniprot_sprot.fasta',
                    '-d', 'Data/databases/uniprot_sprot'], check=True)
    print("  ✓ Database created")
else:
    print("  ✓ Database already exists")

# Run DIAMOND alignment
print("\nRunning DIAMOND alignment for S. brooksi...")
subprocess.run(['diamond', 'blastp',
                '--query', 'Data/transcriptomics/transcriptomes/brooksi_proteins.faa',
                '--db', 'Data/databases/uniprot_sprot',
                '--out', 'Results/05_transcriptomics/brooksi_uniprot.tsv',
                '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle',
                '--max-target-seqs', '1', '--evalue', '1e-5', '--threads', '4'], check=True)
print("  ✓ Alignment complete")

print("\nRunning DIAMOND alignment for S. elizabethae...")
subprocess.run(['diamond', 'blastp',
                '--query', 'Data/transcriptomics/transcriptomes/elizabethae_proteins.faa',
                '--db', 'Data/databases/uniprot_sprot',
                '--out', 'Results/05_transcriptomics/elizabethae_uniprot.tsv',
                '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle',
                '--max-target-seqs', '1', '--evalue', '1e-5', '--threads', '4'], check=True)
print("  ✓ Alignment complete")

print("\n✓ Step 3 complete")
