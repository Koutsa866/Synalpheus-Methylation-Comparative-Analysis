#!/bin/bash
set -e

# Translate transcripts to proteins
transeq -sequence Data/transcriptomics/transcriptomes/S.brooksi_brain_filtered_transcriptome.fasta \
        -outseq Data/transcriptomics/transcriptomes/S.brooksi_proteins.faa -frame 6 -clean

transeq -sequence Data/transcriptomics/transcriptomes/S.elizabethae_brain_filtered_transcriptome.fasta \
        -outseq Data/transcriptomics/transcriptomes/S.elizabethae_proteins.faa -frame 6 -clean

# Get S. chacei proteins
grep -A1 "^>" Results/01_gene_prediction/predicted_proteins.faa | grep -v "^--$" > Data/chacei_proteins.faa

# Make DIAMOND database
diamond makedb --in Data/chacei_proteins.faa -d Data/chacei_db

# Run DIAMOND alignment
diamond blastp --query Data/transcriptomics/transcriptomes/S.brooksi_proteins.faa \
               --db Data/chacei_db \
               --out Results/06_deg_methylation/brooksi_to_chacei_new.tsv \
               --outfmt 6 qseqid sseqid pident length evalue bitscore \
               --max-target-seqs 1 --evalue 1e-5

diamond blastp --query Data/transcriptomics/transcriptomes/S.elizabethae_proteins.faa \
               --db Data/chacei_db \
               --out Results/06_deg_methylation/elizabethae_to_chacei_new.tsv \
               --outfmt 6 qseqid sseqid pident length evalue bitscore \
               --max-target-seqs 1 --evalue 1e-5

echo "✓ DIAMOND alignment complete"
