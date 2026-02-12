#!/bin/bash
# Run DIAMOND alignment with CORRECT filtered transcriptomes

BASE="/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir"
OUT="${BASE}/Results/06_deg_methylation"

echo "Running DIAMOND blastx for S. elizabethae (filtered transcriptome)..."
diamond blastx \
    --query "${BASE}/Data/transcriptomics/transcriptomes/S.elizabethae_brain_filtered_transcriptome.fasta" \
    --db "${BASE}/Data/databases/chacei_proteins.dmnd" \
    --out "${OUT}/elizabethae_filtered_to_chacei.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --threads 4

echo "Running DIAMOND blastx for S. brooksi (filtered transcriptome)..."
diamond blastx \
    --query "${BASE}/Data/transcriptomics/transcriptomes/S.brooksi_brain_filtered_transcriptome.fasta" \
    --db "${BASE}/Data/databases/chacei_proteins.dmnd" \
    --out "${OUT}/brooksi_filtered_to_chacei.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --threads 4

echo "DIAMOND alignment complete!"
echo "Elizabethae results: ${OUT}/elizabethae_filtered_to_chacei.tsv"
echo "Brooksi results: ${OUT}/brooksi_filtered_to_chacei.tsv"
