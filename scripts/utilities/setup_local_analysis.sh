#!/bin/bash
# Local Setup Script for Shrimp Methylation Analysis
# Run this script to set up your local environment

set -e

echo "=========================================="
echo "Shrimp Analysis Local Setup"
echo "=========================================="

# Step 1: Create directory
echo "Step 1: Creating local data directory..."
mkdir -p ~/shrimp_local_data
cd ~/shrimp_local_data

# Step 2: Check conda installation
echo "Step 2: Checking conda installation..."
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda not found. Please install miniconda:"
    echo "brew install --cask miniconda"
    exit 1
fi

# Step 3: Create conda environment
echo "Step 3: Creating conda environment..."
if conda env list | grep -q "shrimp_local"; then
    echo "Environment 'shrimp_local' already exists. Skipping creation."
else
    conda create -n shrimp_local python=3.10 -y
fi

# Step 4: Install dependencies
echo "Step 4: Installing dependencies..."
eval "$(conda shell.bash hook)"
conda activate shrimp_local
conda install -c bioconda biopython diamond -y

# Step 5: Download SwissProt database
echo "Step 5: Downloading SwissProt database (~90 MB)..."
if [ ! -f "uniprot_sprot.fasta.gz" ]; then
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
else
    echo "SwissProt database already exists. Skipping download."
fi

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Download data from OSC:"
echo "   cd ~/shrimp_local_data"
echo "   scp -r YOUR_OSC_USERNAME@ascend.osc.edu:/users/PDNU0014/usrname2/shrimp_data_tables ./"
echo "   scp -r YOUR_OSC_USERNAME@ascend.osc.edu:/users/PDNU0014/usrname2/chacei_braker_rnaseq_V12 ./"
echo ""
echo "2. Run the analysis:"
echo "   conda activate shrimp_local"
echo "   python3 ~/path/to/shrimp_local_analysis.py"
