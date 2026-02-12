#!/bin/bash
# Diagnostic script to check data compatibility
# Run on OSC: bash check_data_compatibility.sh

TABLES_DIR="/users/PDNU0014/usrname2/shrimp_data_tables"
BRAKER_DIR="/users/PDNU0014/usrname2/chacei_braker_rnaseq_V12"

echo "=== DATA COMPATIBILITY CHECK ==="
echo ""

echo "1. Checking methylation table contigs:"
awk 'NR>1 {print $3}' "$TABLES_DIR/genes_methylated_promoters.tsv" | sort -u | head -10
echo "   Total unique contigs in methylation data: $(awk 'NR>1 {print $3}' "$TABLES_DIR/genes_methylated_promoters.tsv" | sort -u | wc -l)"
echo ""

echo "2. Checking BRAKER GTF contigs:"
cut -f1 "$BRAKER_DIR/braker.gtf" | grep -v "^#" | sort -u | head -10
echo "   Total unique contigs in BRAKER GTF: $(cut -f1 "$BRAKER_DIR/braker.gtf" | grep -v "^#" | sort -u | wc -l)"
echo ""

echo "3. Checking for overlap:"
comm -12 \
  <(awk 'NR>1 {print $3}' "$TABLES_DIR/genes_methylated_promoters.tsv" | sort -u) \
  <(cut -f1 "$BRAKER_DIR/braker.gtf" | grep -v "^#" | sort -u) | head -10
echo "   Overlapping contigs: $(comm -12 <(awk 'NR>1 {print $3}' "$TABLES_DIR/genes_methylated_promoters.tsv" | sort -u) <(cut -f1 "$BRAKER_DIR/braker.gtf" | grep -v "^#" | sort -u) | wc -l)"
echo ""

echo "4. Sample gene IDs from methylation table:"
awk 'NR>1 && NR<6 {print $1, $3}' "$TABLES_DIR/genes_methylated_promoters.tsv"
echo ""

echo "5. Sample gene IDs from BRAKER GTF:"
awk -F'\t' '$3=="CDS" {print $1, $9; exit}' "$BRAKER_DIR/braker.gtf"
echo ""

echo "=== DIAGNOSIS ==="
OVERLAP=$(comm -12 <(awk 'NR>1 {print $3}' "$TABLES_DIR/genes_methylated_promoters.tsv" | sort -u) <(cut -f1 "$BRAKER_DIR/braker.gtf" | grep -v "^#" | sort -u) | wc -l)

if [ "$OVERLAP" -gt 100 ]; then
    echo "✓ Data is compatible - contigs match between files"
    echo "  You can proceed with the analysis"
else
    echo "✗ Data incompatibility detected!"
    echo "  The methylation tables and BRAKER files use different contig names"
    echo "  Possible causes:"
    echo "  - Different genome assemblies were used"
    echo "  - Contigs were renamed between analyses"
    echo "  - Files are from different versions of the data"
    echo ""
    echo "  ACTION NEEDED: Check which assembly was used for:"
    echo "  - Methylation analysis: $TABLES_DIR"
    echo "  - BRAKER annotation: $BRAKER_DIR"
fi