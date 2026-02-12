#!/usr/bin/env python3
"""
Upload methylation analysis project to GitHub
Prepares and syncs the complete methylation study
"""

import os
import shutil
from pathlib import Path
import subprocess

# Paths
BASE_DIR = Path("/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir")
GITHUB_DIR = BASE_DIR / "Refactored_Scripts"

def prepare_methylation_project():
    """Organize methylation analysis for GitHub upload"""
    
    # Create project structure
    project_dirs = [
        "methylation_analysis/scripts",
        "methylation_analysis/results/figures", 
        "methylation_analysis/docs",
        "methylation_analysis/sample_data"
    ]
    
    for dir_path in project_dirs:
        (GITHUB_DIR / dir_path).mkdir(parents=True, exist_ok=True)
    
    # Copy key scripts
    scripts_to_copy = [
        ("Scripts/01_gene_prediction/augustus_run_s2.py", "methylation_analysis/scripts/01_gene_prediction.py"),
        ("Scripts/03_methylation_analysis/methylation_analysis.py", "methylation_analysis/scripts/02_methylation_analysis.py"),
        ("Scripts/04_functional_annotation/functional_annotation.py", "methylation_analysis/scripts/03_functional_annotation.py"),
        ("Scripts/04_visualization/create_te_enrichment_figure.py", "methylation_analysis/scripts/04_visualization.py")
    ]
    
    for src, dst in scripts_to_copy:
        src_path = BASE_DIR / src
        dst_path = GITHUB_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {src}")
    
    # Copy key results/figures
    figures_to_copy = [
        "Results/05_figures/te_enrichment_bar.png",
        "Results/05_figures/te_enrichment_pie.png"
    ]
    
    for fig in figures_to_copy:
        src_path = BASE_DIR / fig
        dst_path = GITHUB_DIR / f"methylation_analysis/results/figures/{Path(fig).name}"
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {fig}")
    
    # Create README
    readme_content = """# DNA Methylation Analysis in *Synalpheus chacei*

## Abstract
Synalpheus chacei is a eusocial snapping shrimp with reproductive division of labor, making it a rare marine model for studying epigenetic regulation. We performed assembly-wide CpG methylation (5mC) analysis to test whether promoter methylation regulates social behavior genes.

## Key Findings
- **6,848 genes** predicted via AUGUSTUS from Oxford Nanopore assembly
- **1,103 genes (16.1%)** exhibit highly methylated promoters
- **77.4%** of annotated methylated genes are transposable elements
- **p = 3.9×10⁻¹³** for DNA binding enrichment in methylated genes

## Results Summary
Instead of regulating social behavior genes, promoter methylation in *S. chacei* serves as a **genome defense mechanism**, silencing transposable elements to maintain chromosomal stability essential for long-lived reproductive castes.

## Repository Contents
- `scripts/`: Complete bioinformatics pipeline
- `results/figures/`: Publication-ready figures
- `docs/`: Methods and documentation

## Methods
- **Sequencing**: Oxford Nanopore long-read sequencing
- **Gene Prediction**: AUGUSTUS with Drosophila training set
- **Methylation Calling**: ≥30× coverage, ≥50% methylation threshold
- **Functional Annotation**: DIAMOND + SwissProt database
- **Statistical Analysis**: GO enrichment with Fisher's exact test

## Citation
Koutsaftis, P. (2025). DNA methylation targets transposable elements in the eusocial shrimp *Synalpheus chacei*. PAG33 Conference.

## Contact
Philip Koutsaftis - koutsa_p1@denison.edu
"""
    
    readme_path = GITHUB_DIR / "methylation_analysis/README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    print("✓ Created README.md")
    
    print("\n🎉 Project prepared for GitHub upload!")
    print(f"📁 Files organized in: {GITHUB_DIR}/methylation_analysis/")

def sync_to_github():
    """Run the existing GitHub sync script"""
    sync_script = BASE_DIR / "Scripts/utilities/auto_github_sync.py"
    
    print("\n🚀 Starting GitHub sync...")
    subprocess.run([
        "python3", str(sync_script), 
        "Add complete methylation analysis project - PAG33 ready"
    ])

if __name__ == "__main__":
    prepare_methylation_project()
    
    # Ask user if they want to sync now
    response = input("\nSync to GitHub now? (y/n): ").lower()
    if response == 'y':
        sync_to_github()
    else:
        print("Run the sync later with: python3 Scripts/utilities/auto_github_sync.py")