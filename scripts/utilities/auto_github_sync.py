#!/usr/bin/env python3
"""
Auto-sync Google Drive files to GitHub
Run this script whenever you want to push updates to GitHub
"""

import subprocess
import os
from datetime import datetime

# Configuration
REPO_DIR = os.path.expanduser("~/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir/Refactored_Scripts")
GITHUB_REPO = "https://github.com/Koutsa866/Epigenetic-Research-Synalpheus-Shrimp---2025.git"
GIT_EMAIL = "koutsa_p1@denison.edu"
GIT_NAME = "Philip Koutsaftis"

def run_command(cmd, cwd=None, capture=True):
    """Run shell command. Capture=True for background tasks, False for interactive (passwords)"""
    if capture:
        result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
        return result.returncode == 0, result.stdout, result.stderr
    else:
        # This version allows you to see "Username:" and "Password:" prompts
        result = subprocess.run(cmd, shell=True, cwd=cwd)
        return result.returncode == 0, "", ""

def ensure_gitignore():
    """SAFETY SHIELD: Ensures large data files aren't uploaded to GitHub"""
    gitignore_path = os.path.join(REPO_DIR, ".gitignore")
    if not os.path.exists(gitignore_path):
        print("🛡️ Creating .gitignore to protect large data files...")
        content = "*.fasta\n*.fastq\n*.csv\n*.tsv\n*.gz\n.DS_Store\nData/\nResults/\n"
        with open(gitignore_path, "w") as f:
            f.write(content)

def sync_to_github(commit_message=None):
    """Sync local changes to GitHub"""
    if not os.path.exists(REPO_DIR):
        print(f"❌ Error: Folder not found at {REPO_DIR}")
        return

    os.chdir(REPO_DIR)
    
    # Run the safety shield first
    ensure_gitignore()
    
    # Check if git repo exists
    if not os.path.exists(".git"):
        print("Initializing Git repository...")
        run_command("git init")
        run_command(f"git config --local user.email '{GIT_EMAIL}'")
        run_command(f"git config --local user.name '{GIT_NAME}'")
        # Try to add remote, ignore if it already exists
        subprocess.run(f"git remote add origin {GITHUB_REPO}", shell=True, capture_output=True)
        run_command("git branch -M main")
    
    # Add all changes
    print("Adding files...")
    run_command("git add .")
    
    # Check if there are changes
    success, stdout, _ = run_command("git status --porcelain")
    if not stdout.strip():
        print("✅ No changes to commit")
        return
    
    # Commit
    if not commit_message:
        commit_message = f"Auto-update: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    
    print(f"Committing: {commit_message}")
    run_command(f'git commit -m "{commit_message}"')
    
    # Push - INTERACTIVE MODE
    print("\nPushing to GitHub...")
    print("👉 When prompted, enter your Username (Koutsa866)")
    print("👉 When prompted for 'Password', PASTE your Personal Access Token.\n")
    
    # We set capture=False here so you can see the login prompts!
    success, _, _ = run_command("git push -u origin main", capture=False)
    
    if success:
        print("\n✅ Successfully synced to GitHub!")
    else:
        print("\n❌ Push failed. Check your Token or Internet connection.")

if __name__ == "__main__":
    import sys
    message = " ".join(sys.argv[1:]) if len(sys.argv) > 1 else None
    sync_to_github(message)
