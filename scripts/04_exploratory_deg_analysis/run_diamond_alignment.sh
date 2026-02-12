#!/bin/bash
# DIAMOND alignment of transcriptomes to S. chacei proteins

BASE="/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir"
OUT_DIR="${BASE}/Results/06_deg_methylation"

echo "Running DIAMOND blastx for S. elizabethae..."
diamond blastx \
    --query "${BASE}/Data/Elizabethae_Assembly_Trinity.fasta" \
    --db "${BASE}/Data/databases/chacei_proteins.dmnd" \
    --out "${OUT_DIR}/elizabethae_to_chacei.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --threads 4

echo "Running DIAMOND blastx for S. brooksi..."
diamond blastx \
    --query "${BASE}/Data/Brooksi_Assembly_Trinity.fasta" \
    --db "${BASE}/Data/databases/chacei_proteins.dmnd" \
    --out "${OUT_DIR}/brooksi_to_chacei.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --threads 4

echo "DIAMOND alignment complete!"
