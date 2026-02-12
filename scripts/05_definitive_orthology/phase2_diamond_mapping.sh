#!/bin/bash
# Phase 2: DIAMOND Orthology Mapping (Cloud-Safe Version)
set -e

# --- PATHS ---
QUERY="Results/05_ortholog_analysis/S_chacei_gene_bodies.fasta"
OUTDIR="Results/05_ortholog_analysis"
LOCAL_TMP="/tmp/diamond_run_$(date +%s)"

mkdir -p "$OUTDIR"
mkdir -p "$LOCAL_TMP"

# Cleanup function: Ensures /tmp is cleaned even if the script crashes
trap "rm -rf $LOCAL_TMP; exit" INT TERM EXIT

echo "=== Phase 2: DIAMOND Orthology Mapping (Local I/O) ==="

# 1. Build DBs and Run Alignments
for SP in "brooksi" "elizabethae"; do
    echo "Building $SP database in $LOCAL_TMP..."
    diamond makedb \
      --in Data/transcriptomics/transcriptomes/S.${SP}_strict.faa \
      -d "${LOCAL_TMP}/${SP}_db"

    echo "Mapping S. chacei to $SP..."
    # Write the .tsv to LOCAL_TMP first to avoid Google Drive timeouts
    diamond blastx \
      --query "$QUERY" \
      --db "${LOCAL_TMP}/${SP}_db" \
      --out "${LOCAL_TMP}/chacei_vs_${SP}.tsv" \
      --max-target-seqs 1 \
      --sensitive \
      --evalue 1e-10 \
      --outfmt 6 qseqid sseqid pident length evalue bitscore

    # 2. Move the final file back to the Google Drive folder
    echo "Moving results to $OUTDIR..."
    cp "${LOCAL_TMP}/chacei_vs_${SP}.tsv" "${OUTDIR}/chacei_vs_${SP}.tsv"
done

echo ""
echo "✓ Phase 2 complete!"
echo "  S. brooksi hits: $(wc -l < ${OUTDIR}/chacei_vs_brooksi.tsv)"
echo "  S. elizabethae hits: $(wc -l < ${OUTDIR}/chacei_vs_elizabethae.tsv)"
