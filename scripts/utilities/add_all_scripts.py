#!/usr/bin/env python3
"""
Add all key Python scripts to GitHub repository
"""

import os
import shutil
from pathlib import Path
import subprocess

# Paths
BASE_DIR = Path("/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir")
GITHUB_DIR = BASE_DIR / "Refactored_Scripts"

def copy_all_scripts():
    """Copy all key analysis scripts"""
    
    # Ensure directories exist
    script_dirs = [
        "methylation_analysis/scripts/gene_prediction",
        "methylation_analysis/scripts/methylation_analysis", 
        "methylation_analysis/scripts/functional_annotation",
        "methylation_analysis/scripts/visualization",
        "methylation_analysis/scripts/utilities"
    ]
    
    for dir_path in script_dirs:
        (GITHUB_DIR / dir_path).mkdir(parents=True, exist_ok=True)
    
    # Comprehensive script copying
    scripts_to_copy = [
        # Gene prediction
        ("Scripts/01_gene_prediction/augustus_run_s2.py", "methylation_analysis/scripts/gene_prediction/augustus_gene_prediction.py"),
        
        # Methylation analysis
        ("Scripts/03_methylation_analysis/methylation_analysis.py", "methylation_analysis/scripts/methylation_analysis/methylation_calling.py"),
        ("Scripts/compare_methylated_unmethylated.py", "methylation_analysis/scripts/methylation_analysis/compare_methylation_groups.py"),
        
        # Functional annotation
        ("Scripts/03_functional_annotation/functional_annotation.py", "methylation_analysis/scripts/functional_annotation/diamond_annotation.py"),
        
        # Visualization
        ("Scripts/04_visualization/create_te_enrichment_figure.py", "methylation_analysis/scripts/visualization/te_enrichment_plots.py"),
        ("Scripts/create_bold_go_chart.py", "methylation_analysis/scripts/visualization/go_enrichment_chart.py"),
        ("Scripts/create_cpg_distribution_plots.py", "methylation_analysis/scripts/visualization/cpg_distribution.py"),
        ("Scripts/promoter_methylation_chart.py", "methylation_analysis/scripts/visualization/promoter_methylation.py"),
        
        # Utilities
        ("Scripts/utilities/ContigStats_refactored.py", "methylation_analysis/scripts/utilities/contig_statistics.py"),
        ("cleanup_duplicates.py", "methylation_analysis/scripts/utilities/cleanup_duplicates.py"),
        ("extract_sequences.py", "methylation_analysis/scripts/utilities/extract_sequences.py")
    ]
    
    copied_count = 0
    for src, dst in scripts_to_copy:
        src_path = BASE_DIR / src
        dst_path = GITHUB_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {src_path.name}")
            copied_count += 1
        else:
            print(f"⚠️  Not found: {src}")
    
    print(f"\n✅ Copied {copied_count} Python scripts")
    return copied_count > 0

def copy_additional_files():
    """Copy additional important files"""
    
    additional_files = [
        ("requirements.txt", "methylation_analysis/requirements.txt"),
        ("Poster_Materials/PAG33_PhilipK_122425550.pdf", "methylation_analysis/docs/PAG33_poster.pdf")
    ]
    
    for src, dst in additional_files:
        src_path = BASE_DIR / src
        dst_path = GITHUB_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {src_path.name}")

def force_git_update():
    """Force git to recognize and push changes"""
    
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
            "Add complete Python analysis pipeline - methylation scripts"
        ], check=True)
        
        # Push to GitHub
        print("\n🚀 Pushing to GitHub...")
        subprocess.run(["git", "push"], check=True)
        print("✅ Successfully pushed to GitHub!")
        
    else:
        print("\n⚠️  No new changes detected")

if __name__ == "__main__":
    print("📁 Adding all Python scripts to GitHub repository...")
    
    # Copy all scripts
    changes_made = copy_all_scripts()
    
    # Copy additional files
    copy_additional_files()
    
    if changes_made:
        # Force git update
        try:
            force_git_update()
        except subprocess.CalledProcessError as e:
            print(f"❌ Git error: {e}")
            print("Try running manually: cd Refactored_Scripts && git add . && git commit -m 'Add scripts' && git push")
    
    print("\n🎉 Script copying complete!")
    print("Check your GitHub repository for the new scripts.")