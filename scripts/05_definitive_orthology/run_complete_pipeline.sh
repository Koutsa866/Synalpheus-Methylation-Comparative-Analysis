#!/bin/bash
# Master Script: Run Complete Ortholog Pipeline
# Executes Phases 2a, 2, and 3 in sequence

set -e

echo "============================================================"
echo "ORTHOLOG PIPELINE: COMPLETE EXECUTION"
echo "============================================================"
echo ""

# Phase 2a: Translate DNA to Protein (Strict ID Preservation)
echo "--- PHASE 2a: Strict-ID Translation ---"
python scripts/06_ortholog_pipeline/phase2a_translate_strict.py
echo ""

# Phase 2: DIAMOND Orthology Mapping
echo "--- PHASE 2: DIAMOND Mapping ---"
bash scripts/06_ortholog_pipeline/phase2_diamond_mapping.sh
echo ""

# Phase 3: Statistical Analysis
echo "--- PHASE 3: Statistical Analysis ---"
python scripts/06_ortholog_pipeline/phase3_statistical_analysis.py

echo ""
echo "============================================================"
echo "✓ PIPELINE COMPLETE"
echo "============================================================"
