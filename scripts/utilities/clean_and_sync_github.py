#!/usr/bin/env python3
"""
Clean and organize repository for GitHub sync
Removes old README files and organizes project structure
"""

import os
import shutil
from pathlib import Path
import subprocess

# Paths
BASE_DIR = Path("/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir")
GITHUB_DIR = BASE_DIR / "Refactored_Scripts"

def cleanup_old_readmes():
    """Remove old README files to avoid confusion"""
    
    old_readmes = [
        GITHUB_DIR / "README.md",  # Any existing README in sync folder
        BASE_DIR / "docs/README.md",  # Old docs README
        BASE_DIR / "docs/README_updated.md",  # Updated version we created
        GITHUB_DIR / "methylation_analysis/README.md"  # Any existing project README
    ]
    
    for readme in old_readmes:
        if readme.exists():
            readme.unlink()
            print(f"✓ Removed old README: {readme.name}")

def copy_main_readme():
    """Copy the main README to the sync directory"""
    
    main_readme = BASE_DIR / "README.md"
    target_readme = GITHUB_DIR / "README.md"
    
    if main_readme.exists():
        shutil.copy2(main_readme, target_readme)
        print(f"✓ Copied main README to sync directory")
    else:
        print("❌ Main README.md not found!")

def organize_methylation_project():
    """Organize the methylation analysis project"""
    
    # Create clean project structure
    project_dirs = [
        "methylation_analysis/scripts",
        "methylation_analysis/results/figures", 
        "methylation_analysis/docs"
    ]
    
    for dir_path in project_dirs:
        (GITHUB_DIR / dir_path).mkdir(parents=True, exist_ok=True)
    
    # Copy key scripts (clean versions without personal paths)
    scripts_to_copy = [
        ("Scripts/01_gene_prediction/augustus_run_s2.py", "methylation_analysis/scripts/01_gene_prediction.py"),
        ("Scripts/04_visualization/create_te_enrichment_figure.py", "methylation_analysis/scripts/02_visualization.py")
    ]
    
    for src, dst in scripts_to_copy:
        src_path = BASE_DIR / src
        dst_path = GITHUB_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied {src}")
    
    # Copy key figures
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

def create_gitignore():
    """Create comprehensive .gitignore"""
    
    gitignore_content = """# Data files
*.fasta
*.fastq
*.csv
*.tsv
*.gz
*.zip
*.bed
*.bam
*.sam

# System files
.DS_Store
Thumbs.db

# Large directories
Data/
Results/
archive/
notebooks/
Poster_Materials/

# Temporary files
*.tmp
*.log
*~

# Python
__pycache__/
*.pyc
*.pyo
"""
    
    gitignore_path = GITHUB_DIR / ".gitignore"
    with open(gitignore_path, 'w') as f:
        f.write(gitignore_content)
    print("✓ Created comprehensive .gitignore")

def sync_to_github():
    """Run the GitHub sync"""
    sync_script = BASE_DIR / "Scripts/utilities/auto_github_sync.py"
    
    print("\n🚀 Starting GitHub sync...")
    subprocess.run([
        "python3", str(sync_script), 
        "Clean repository structure - PAG33 ready with main README"
    ])

if __name__ == "__main__":
    print("🧹 Cleaning up repository for GitHub sync...")
    
    # Step 1: Clean up old READMEs
    cleanup_old_readmes()
    
    # Step 2: Copy main README
    copy_main_readme()
    
    # Step 3: Organize project
    organize_methylation_project()
    
    # Step 4: Create .gitignore
    create_gitignore()
    
    print("\n✅ Repository cleaned and organized!")
    print(f"📁 Sync directory: {GITHUB_DIR}")
    print("📄 Main README.md will be the repository homepage")
    
    # Ask user if they want to sync now
    response = input("\nSync to GitHub now? (y/n): ").lower()
    if response == 'y':
        sync_to_github()
    else:
        print("Run sync later with: python3 Scripts/utilities/auto_github_sync.py")