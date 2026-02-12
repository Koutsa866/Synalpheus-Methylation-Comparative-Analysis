#!/bin/bash
# Master script to run Phase 1: Penaeus vannamei re-annotation

set -e  # Exit on error

SCRIPT_DIR="/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir/Scripts/03_functional_annotation"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  PHASE 1: PENAEUS VANNAMEI RE-ANNOTATION                  ║"
echo "║  Estimated time: 30-45 minutes                            ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Step 1: Extract unannotated proteins
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1/4: Extracting unannotated proteins"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
python3 "$SCRIPT_DIR/extract_unannotated_proteins.py"

# Step 2: Setup Penaeus database (skip if already exists)
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2/4: Setting up Penaeus database"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [ ! -f "$SCRIPT_DIR/../../Data/databases/penaeus_db.dmnd" ]; then
    bash "$SCRIPT_DIR/setup_penaeus_db.sh"
else
    echo "Database already exists, skipping download..."
fi

# Step 3: Run Diamond BLAST
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 3/4: Running Diamond BLAST (15-30 min)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
python3 "$SCRIPT_DIR/run_penaeus_annotation.py"

# Step 4: Merge annotations
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 4/4: Merging annotations"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
python3 "$SCRIPT_DIR/merge_annotations.py"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  ✓ PHASE 1 COMPLETE                                       ║"
echo "║                                                            ║"
echo "║  Next steps:                                               ║"
echo "║  1. Review Results/04_functional_annotation/               ║"
echo "║     combined_annotations.csv                               ║"
echo "║  2. Re-run GO enrichment with new annotations              ║"
echo "║  3. (Optional) Run Phase 2: InterProScan                   ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
