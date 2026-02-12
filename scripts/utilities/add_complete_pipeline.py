#!/usr/bin/env python3
"""
Add complete methylation analysis pipeline to GitHub
Includes all scripts, notebooks, and poster figures
"""

import os
import shutil
from pathlib import Path
import subprocess

# Paths
BASE_DIR = Path("/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir")
GITHUB_DIR = BASE_DIR / "Refactored_Scripts"

def copy_complete_pipeline():
    """Copy all key analysis scripts and organize properly"""
    
    # Create comprehensive directory structure
    dirs_to_create = [
        "methylation_analysis/scripts/gene_prediction",
        "methylation_analysis/scripts/feature_extraction", 
        "methylation_analysis/scripts/functional_annotation",
        "methylation_analysis/scripts/visualization",
        "methylation_analysis/scripts/utilities",
        "methylation_analysis/notebooks",
        "methylation_analysis/results/figures",
        "methylation_analysis/results/poster_figures",
        "methylation_analysis/docs"
    ]
    
    for dir_path in dirs_to_create:
        (GITHUB_DIR / dir_path).mkdir(parents=True, exist_ok=True)
    
    # Complete script collection
    scripts_to_copy = [
        # Gene prediction
        ("Scripts/01_gene_prediction/augustus_run_s2.py", "methylation_analysis/scripts/gene_prediction/augustus_gene_prediction.py"),
        
        # Feature extraction
        ("Scripts/02_feature_extraction/bed_files_processing.py", "methylation_analysis/scripts/feature_extraction/bed_processing.py"),
        
        # Functional annotation (complete pipeline)
        ("Scripts/03_functional_annotation/run_go_enrichment.py", "methylation_analysis/scripts/functional_annotation/go_enrichment_analysis.py"),
        ("Scripts/03_functional_annotation/merge_annotations.py", "methylation_analysis/scripts/functional_annotation/annotation_merger.py"),
        ("Scripts/03_functional_annotation/blastp_analysis.py", "methylation_analysis/scripts/functional_annotation/blastp_analysis.py"),
        ("Scripts/03_functional_annotation/extract_unannotated_proteins.py", "methylation_analysis/scripts/functional_annotation/extract_unannotated.py"),
        ("Scripts/03_functional_annotation/functional_annotation_prep.py", "methylation_analysis/scripts/functional_annotation/annotation_prep.py"),
        ("Scripts/03_functional_annotation/run_local_analysis.py", "methylation_analysis/scripts/functional_annotation/local_analysis.py"),
        
        # Visualization (complete set)
        ("Scripts/04_visualization/create_te_enrichment_figure.py", "methylation_analysis/scripts/visualization/te_enrichment_plots.py"),
        ("Scripts/04_visualization/create_donut_chart.py", "methylation_analysis/scripts/visualization/donut_charts.py"),
        ("Scripts/04_visualization/create_figure.py", "methylation_analysis/scripts/visualization/general_figures.py"),
        ("Scripts/create_bold_go_chart.py", "methylation_analysis/scripts/visualization/go_enrichment_charts.py"),
        ("Scripts/create_cpg_distribution_plots.py", "methylation_analysis/scripts/visualization/cpg_distribution_plots.py"),
        
        # Core analysis
        ("Scripts/compare_methylated_unmethylated.py", "methylation_analysis/scripts/methylation_comparison.py"),
        
        # Utilities
        ("Scripts/utilities/ContigStats_refactored.py", "methylation_analysis/scripts/utilities/contig_statistics.py"),
        ("cleanup_duplicates.py", "methylation_analysis/scripts/utilities/cleanup_duplicates.py"),
        ("extract_sequences.py", "methylation_analysis/scripts/utilities/sequence_extraction.py"),
        ("reorganize_project.py", "methylation_analysis/scripts/utilities/project_reorganizer.py")
    ]
    
    # Copy notebooks
    notebooks_to_copy = [
        ("Scripts/promoter_cpgIsland_file_processing.ipynb", "methylation_analysis/notebooks/promoter_cpg_analysis.ipynb"),
        ("Scripts/S_elizabithae_AUGUSTUS.ipynb", "methylation_analysis/notebooks/augustus_analysis.ipynb")
    ]
    
    # Copy poster figures and materials
    poster_files = [
        ("Poster_Materials/Methods_1.0.png", "methylation_analysis/results/poster_figures/methods_flowchart_v1.png"),
        ("Poster_Materials/Methods_2.0.png", "methylation_analysis/results/poster_figures/methods_flowchart_v2.png"),
        ("Poster_Materials/Methods_3.0.png", "methylation_analysis/results/poster_figures/methods_flowchart_v3.png"),
        ("Poster_Materials/PAG33_PhilipK_122425550/Slide1.png", "methylation_analysis/results/poster_figures/final_poster_slide.png"),
        ("Poster_Materials/PAG33_PhilipK_122425550.pdf", "methylation_analysis/docs/PAG33_final_poster.pdf"),
        ("Poster_Materials/Philip_Koutsaftis_PAG33Poster (1).pdf", "methylation_analysis/docs/PAG33_poster_alt.pdf"),
        ("Poster_Materials/unnamed.png", "methylation_analysis/results/poster_figures/poster_preview.png")
    ]
    
    # Copy extracted poster media
    poster_media = [
        ("Poster_Materials/pptx_extracted/ppt/media/image1.png", "methylation_analysis/results/poster_figures/poster_figure_1.png"),
        ("Poster_Materials/pptx_extracted/ppt/media/image2.png", "methylation_analysis/results/poster_figures/poster_figure_2.png"),
        ("Poster_Materials/pptx_extracted/ppt/media/image4.png", "methylation_analysis/results/poster_figures/poster_figure_4.png")
    ]
    
    # Copy existing result figures
    result_figures = [
        ("Results/05_figures/te_enrichment_bar.png", "methylation_analysis/results/figures/te_enrichment_bar_chart.png"),
        ("Results/05_figures/te_enrichment_pie.png", "methylation_analysis/results/figures/te_enrichment_pie_chart.png"),
        ("promoter_methylation_profile.png", "methylation_analysis/results/figures/promoter_methylation_profile.png"),
        ("promoter_methylation_profile.pdf", "methylation_analysis/results/figures/promoter_methylation_profile.pdf")
    ]
    
    # Copy all files
    all_files = scripts_to_copy + notebooks_to_copy + poster_files + poster_media + result_figures
    
    copied_count = 0
    for src, dst in all_files:
        src_path = BASE_DIR / src
        dst_path = GITHUB_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {src_path.name}")
            copied_count += 1
        else:
            print(f"⚠️  Not found: {src}")
    
    # Copy additional important files
    additional_files = [
        ("requirements.txt", "methylation_analysis/requirements.txt"),
        ("setup_local_analysis.sh", "methylation_analysis/setup_analysis.sh")
    ]
    
    for src, dst in additional_files:
        src_path = BASE_DIR / src
        dst_path = GITHUB_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {src_path.name}")
            copied_count += 1
    
    print(f"\n✅ Copied {copied_count} files total")
    return copied_count > 0

def create_project_readme():
    """Create a README for the methylation_analysis subfolder"""
    
    readme_content = """# Methylation Analysis Pipeline

Complete bioinformatics pipeline for DNA methylation analysis in *Synalpheus chacei*.

## Directory Structure

### `/scripts/`
- **`gene_prediction/`** - AUGUSTUS gene calling
- **`feature_extraction/`** - BED file processing and genomic features
- **`functional_annotation/`** - DIAMOND, BLAST, and GO enrichment
- **`visualization/`** - All plotting and figure generation scripts
- **`utilities/`** - Helper scripts and data processing tools

### `/notebooks/`
- **`promoter_cpg_analysis.ipynb`** - Interactive CpG island analysis
- **`augustus_analysis.ipynb`** - Gene prediction workflow

### `/results/`
- **`figures/`** - Publication-ready analysis figures
- **`poster_figures/`** - PAG33 conference poster materials

### `/docs/`
- **`PAG33_final_poster.pdf`** - Conference presentation
- **`PAG33_poster_alt.pdf`** - Alternative poster version

## Key Results
- 6,848 genes predicted via AUGUSTUS
- 1,103 genes (16.1%) with methylated promoters  
- 77.4% of annotated methylated genes are transposable elements
- p = 3.9×10⁻¹³ for DNA binding enrichment

## Usage
1. Install requirements: `pip install -r requirements.txt`
2. Run setup: `bash setup_analysis.sh`
3. Execute pipeline scripts in order
4. Generate figures using visualization scripts
"""
    
    readme_path = GITHUB_DIR / "methylation_analysis/README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    print("✓ Created methylation_analysis/README.md")

def commit_and_push():
    """Commit and push all changes to GitHub"""
    
    os.chdir(GITHUB_DIR)
    
    # Add all files
    subprocess.run(["git", "add", "."], check=True)
    
    # Check status
    result = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True)
    
    if result.stdout.strip():
        print(f"\n📝 Changes detected:")
        print(result.stdout)
        
        # Commit changes
        subprocess.run([
            "git", "commit", "-m", 
            "Add complete methylation analysis pipeline - scripts, notebooks, and poster figures"
        ], check=True)
        
        # Push to GitHub
        print("\n🚀 Pushing to GitHub...")
        subprocess.run(["git", "push"], check=True)
        print("✅ Successfully pushed complete pipeline to GitHub!")
        
    else:
        print("\n⚠️  No new changes detected")

if __name__ == "__main__":
    print("📦 Adding complete methylation analysis pipeline to GitHub...")
    
    # Copy all files
    changes_made = copy_complete_pipeline()
    
    # Create project README
    create_project_readme()
    
    if changes_made:
        try:
            commit_and_push()
        except subprocess.CalledProcessError as e:
            print(f"❌ Git error: {e}")
            print("Try running manually: cd Refactored_Scripts && git add . && git commit -m 'Add complete pipeline' && git push")
    
    print("\n🎉 Complete pipeline setup finished!")
    print("Your GitHub repository now contains:")
    print("- Complete Python analysis pipeline")
    print("- Jupyter notebooks for interactive analysis") 
    print("- All poster figures and materials")
    print("- Publication-ready result figures")
    print("- Comprehensive documentation")