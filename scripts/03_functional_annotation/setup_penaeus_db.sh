#!/bin/bash
# Download and prepare Penaeus vannamei proteome for Diamond BLAST

set -e  # Exit on error

BASE_DIR="/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir"
DB_DIR="$BASE_DIR/Data/databases"

echo "=========================================="
echo "DOWNLOADING PENAEUS VANNAMEI PROTEOME"
echo "=========================================="

cd "$DB_DIR"

# Download Penaeus vannamei reference proteome
echo "Downloading from UniProt..."
wget -O penaeus_vannamei.fasta.gz \
  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000283509/UP000283509_6689.fasta.gz"

echo "Extracting..."
gunzip penaeus_vannamei.fasta.gz

echo "Building Diamond database..."
diamond makedb \
  --in penaeus_vannamei.fasta \
  --db penaeus_db \
  --threads 8

echo ""
echo "=========================================="
echo "DATABASE READY"
echo "Location: $DB_DIR/penaeus_db.dmnd"
echo "=========================================="
