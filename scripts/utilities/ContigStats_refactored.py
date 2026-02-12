"""
Contig Statistics and Methylation Analysis Pipeline
Analyzes contig lengths, methylation data, and CpG islands
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
from Bio import SeqIO

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Central configuration for all paths and parameters"""
    # Detect environment (Colab vs local)
    IS_COLAB = os.path.exists("/content")
    
    # Base paths
    if IS_COLAB:
        DRIVE_BASE = "/content/drive/MyDrive/2025_SnrProj_Dir"
        OUTDIR = "/content/contig_trim_v1"
    else:
        # Local macOS/Google Drive
        DRIVE_BASE = os.path.expanduser("~/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir")
        # Write to local temp directory instead of Google Drive
        OUTDIR = os.path.expanduser("~/Desktop/contig_trim_v1")
    
    DATA_DIR = f"{DRIVE_BASE}/Data"
    RESULTS_DIR = f"{DRIVE_BASE}/Results"
    FASTA_DIR = DATA_DIR  # All FASTA files location
    
    # Input files
    CONTIG_LENGTHS_CSV = f"{DATA_DIR}/Contig_Lengths_Sheets2.csv"
    BEDMETHYL_FILE = f"{RESULTS_DIR}/filtered_bedmethyl_50p_cov30.csv"
    ASSEMBLY_FASTA = f"{FASTA_DIR}/assembly.fasta"
    CPG_ISLANDS_BED = f"{DATA_DIR}/cpg_islands_sorted.bed"
    GENE_PROMOTERS_BED = f"{DRIVE_BASE}/results (1)/gene_promoters_unique.bed"
    
    # Thresholds
    MIN_CONTIG_LENGTH = 5000  # bp
    METHYLATION_THRESHOLD = 50  # percent
    COVERAGE_THRESHOLD = 30
    HIGH_METHYLATION_THRESHOLD = 70  # percent
    LOW_METHYLATION_THRESHOLD = 30  # percent
    PROMOTER_UPSTREAM = 1000  # bp
    
    # Optional: GFF file for exon/intron analysis
    GENE_GFF = f"{DRIVE_BASE}/results (1)/genes.gff3"  # Update path as needed

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def log(msg):
    """Print timestamped log message"""
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}")
    sys.stdout.flush()

def ensure_dir(path):
    """Create directory if it doesn't exist"""
    os.makedirs(path, exist_ok=True)
    return path

def normalize_contig_ids(series):
    """Normalize contig IDs by stripping quotes and whitespace"""
    return series.astype(str).str.strip().str.replace('"', '', regex=False)

# ============================================================================
# CONTIG LENGTH ANALYSIS
# ============================================================================

def analyze_contig_lengths(csv_path, min_length=5000):
    """Analyze contig length statistics and create histogram"""
    df = pd.read_csv(csv_path)
    
    total_length = df['Contig_Lengths'].sum()
    num_contigs = len(df)
    mean_length = total_length / num_contigs
    std_dev = df['Contig_Lengths'].std()
    
    log(f"Total contigs: {num_contigs:,}")
    log(f"Total length: {total_length:,} bp")
    log(f"Mean length: {mean_length:.2f} bp")
    log(f"Std deviation: {std_dev:.2f} bp")
    
    # Filter and plot
    filtered = df[df['Contig_Lengths'] >= min_length]
    log(f"Contigs ≥{min_length} bp: {len(filtered):,}")
    
    plt.figure(figsize=(10, 6))
    plt.hist(filtered['Contig_Lengths'], bins=50)
    plt.xlabel(f'Contig Lengths (≥ {min_length:,} bp)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Contig Lengths (Cutoff {min_length:,} bp)')
    plt.grid(True)
    plt.show()
    
    return df, filtered

# ============================================================================
# METHYLATION FILTERING
# ============================================================================

def filter_methylation_data(bedmethyl_path, meth_thresh=50, cov_thresh=30):
    """Filter bedmethyl data by methylation and coverage thresholds"""
    df = pd.read_csv(bedmethyl_path)
    df.columns = [c.strip().strip('"').strip("'") for c in df.columns]
    
    required = ["chrom", "start_position", "end_position", "percent_modified", "Nvalid_cov"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    for c in ["percent_modified", "Nvalid_cov"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    
    df["chrom"] = normalize_contig_ids(df["chrom"])
    
    passed = df[(df["percent_modified"] > meth_thresh) & (df["Nvalid_cov"] > cov_thresh)]
    log(f"Methylation sites passing thresholds: {len(passed):,}")
    
    return passed

# ============================================================================
# CONTIG FILTERING
# ============================================================================

def create_allowlist(bedmethyl_df, output_path):
    """Create allowlist of contigs from filtered methylation data"""
    allowset = set(bedmethyl_df["chrom"].dropna().unique())
    
    with open(output_path, "w") as f:
        for cid in sorted(allowset):
            f.write(f"{cid}\n")
    
    log(f"Allowlist contigs: {len(allowset)} → {output_path}")
    return allowset

def filter_assembly(assembly_path, allowset, min_length, output_path):
    """Filter assembly FASTA by allowlist and minimum length"""
    retained = []
    n_total = n_kept = total_bp = 0
    
    with open(output_path, "w") as fout:
        for rec in SeqIO.parse(assembly_path, "fasta"):
            n_total += 1
            rec_id = rec.id.strip()
            length = len(rec.seq)
            
            if rec_id in allowset and length >= min_length:
                SeqIO.write(rec, fout, "fasta")
                retained.append((rec_id, length))
                n_kept += 1
                total_bp += length
    
    log(f"Assembly: {n_total:,} total → {n_kept:,} retained ({total_bp:,} bp)")
    
    return pd.DataFrame(retained, columns=["contig_id", "length_bp"])

# ============================================================================
# CPG ISLAND ANALYSIS
# ============================================================================

def filter_cpg_islands(cpg_bed_path, allowset, output_path):
    """Filter CpG islands to retained contigs"""
    cpg = pd.read_csv(cpg_bed_path, sep="\t", header=None, 
                      names=["chrom", "start", "end"])
    cpg["chrom"] = normalize_contig_ids(cpg["chrom"])
    
    cpg_filtered = cpg[cpg["chrom"].isin(allowset)]
    cpg_filtered.to_csv(output_path, sep="\t", header=False, index=False)
    
    log(f"CpG islands: {len(cpg):,} → {len(cpg_filtered):,} retained")
    return cpg_filtered

def summarize_island_methylation(islands_bed, bedmethyl_csv, output_path):
    """Summarize methylation levels per CpG island after merging adjacent windows"""
    cpg = pd.read_csv(islands_bed, sep="\t", header=None, 
                      names=["chrom", "start", "end"])
    cpg["chrom"] = normalize_contig_ids(cpg["chrom"])
    cpg = cpg.sort_values(["chrom", "start"]).reset_index(drop=True)
    
    # Merge adjacent 100bp windows into biological islands
    log("Merging adjacent CpG windows into biological islands...")
    merged_islands = []
    
    for chrom in sorted(cpg["chrom"].unique()):
        chrom_islands = cpg[cpg["chrom"] == chrom].sort_values("start")
        if chrom_islands.empty:
            continue
        
        current_start = chrom_islands.iloc[0]["start"]
        current_end = chrom_islands.iloc[0]["end"]
        
        for _, row in chrom_islands.iloc[1:].iterrows():
            if row["start"] <= current_end + 1:  # Adjacent or overlapping
                current_end = max(current_end, row["end"])
            else:
                merged_islands.append({"chrom": chrom, "start": current_start, "end": current_end})
                current_start = row["start"]
                current_end = row["end"]
        
        merged_islands.append({"chrom": chrom, "start": current_start, "end": current_end})
    
    cpg = pd.DataFrame(merged_islands)
    log(f"Merged into {len(cpg):,} biological islands")
    
    bm = pd.read_csv(bedmethyl_csv, usecols=["chrom", "start_position", 
                     "end_position", "percent_modified", "Nvalid_cov"])
    bm["chrom"] = normalize_contig_ids(bm["chrom"])
    
    # Restrict to island contigs
    island_contigs = set(cpg["chrom"].unique())
    bm = bm[bm["chrom"].isin(island_contigs)]
    
    cpg = cpg.sort_values(["chrom", "start"]).reset_index(drop=True)
    bm = bm.sort_values(["chrom", "start_position"]).reset_index(drop=True)
    
    records = []
    for chrom in sorted(island_contigs):
        islands_chr = cpg[cpg["chrom"] == chrom]
        sites_chr = bm[bm["chrom"] == chrom]
        
        if islands_chr.empty or sites_chr.empty:
            for _, row in islands_chr.iterrows():
                records.append({
                    "chrom": chrom, "start": int(row["start"]), 
                    "end": int(row["end"]), "n_sites": 0,
                    "mean_percent": np.nan, "max_percent": np.nan
                })
            continue
        
        islands_arr = islands_chr[["start", "end"]].to_numpy()
        sites_arr = sites_chr[["start_position", "end_position", "percent_modified"]].to_numpy()
        
        j = 0
        for istart, iend in islands_arr:
            while j < len(sites_arr) and sites_arr[j, 1] < istart:
                j += 1
            
            pm_vals = []
            k = j
            while k < len(sites_arr) and sites_arr[k, 0] <= iend:
                if not (sites_arr[k, 1] < istart or sites_arr[k, 0] > iend):
                    pm_vals.append(sites_arr[k, 2])
                k += 1
            
            records.append({
                "chrom": chrom, "start": int(istart), "end": int(iend),
                "n_sites": len(pm_vals),
                "mean_percent": float(np.mean(pm_vals)) if pm_vals else np.nan,
                "max_percent": float(np.max(pm_vals)) if pm_vals else np.nan
            })
    
    summary_df = pd.DataFrame(records)
    summary_df.to_csv(output_path, sep="\t", index=False)
    log(f"Island methylation summary: {len(summary_df):,} islands → {output_path}")
    
    return summary_df

# ============================================================================
# CPG ISLAND VISUALIZATION
# ============================================================================

def plot_cpg_island_distributions(cpg_summary, output_dir):
    """Plot CpG island size and methylation distributions with publication-quality styling"""
    cpg = pd.read_csv(cpg_summary, sep="\t")
    cpg['island_size'] = cpg['end'] - cpg['start']
    
    # Classify islands
    cpg['has_methylation'] = cpg['n_sites'] >= 1
    methylated = cpg[cpg['has_methylation']]
    unmethylated = cpg[~cpg['has_methylation']]
    
    # Set publication-quality style
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14
    })
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Size distribution with KDE
    if len(methylated) > 0:
        sns.histplot(methylated['island_size'], bins=50, kde=True, color='#E74C3C', 
                    label='With methylation', alpha=0.6, ax=ax1, stat='density')
    if len(unmethylated) > 0:
        sns.histplot(unmethylated['island_size'], bins=50, kde=True, color='#3498DB',
                    label='Without methylation', alpha=0.6, ax=ax1, stat='density')
    
    ax1.set_xlabel('CpG Island Size (bp)', fontweight='bold')
    ax1.set_ylabel('Density', fontweight='bold')
    ax1.set_title('CpG Island Size Distribution', fontweight='bold', pad=20)
    ax1.legend(frameon=True, fancybox=True, shadow=True)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Percent methylation distribution with KDE
    meth_data = methylated['mean_percent'].dropna()
    if len(meth_data) > 0:
        sns.histplot(meth_data, bins=50, kde=True, color='#F39C12', 
                    alpha=0.7, ax=ax2, stat='density')
    
    ax2.set_xlabel('Mean % Methylation', fontweight='bold')
    ax2.set_ylabel('Density', fontweight='bold')
    ax2.set_title('Distribution of % Methylation Across CpG Islands', fontweight='bold', pad=20)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cpg_island_distributions.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/cpg_island_distributions.pdf", bbox_inches='tight')  # PDF for posters
    plt.show()
    
    log(f"CpG islands with methylation: {len(methylated):,}")
    log(f"CpG islands without methylation: {len(unmethylated):,}")
    log(f"Mean island size (methylated): {methylated['island_size'].mean():.1f} bp")
    log(f"Mean island size (unmethylated): {unmethylated['island_size'].mean():.1f} bp")
    
    return cpg

# ============================================================================
# PROMOTER ANALYSIS
# ============================================================================

def filter_promoters_by_upstream_space(promoters_bed, assembly_fasta, min_upstream=1000, output_path=None):
    """Filter promoters to ensure sufficient upstream space is available"""
    prom = pd.read_csv(promoters_bed, sep="\t", header=None,
                       names=["chrom", "start", "end", "gene_id", "score", "strand"])
    
    # Deduplicate by gene_id first
    n_before = len(prom)
    prom = prom.drop_duplicates(subset=['gene_id'])
    n_after = len(prom)
    if n_before != n_after:
        log(f"Deduplicated: {n_before:,} → {n_after:,} unique genes")
    
    # Get contig lengths
    contig_lengths = {rec.id.strip(): len(rec.seq) for rec in SeqIO.parse(assembly_fasta, "fasta")}
    
    valid_promoters = []
    dropped_edge = 0
    dropped_missing = 0
    dropped_short = 0
    
    for _, row in prom.iterrows():
        contig_len = contig_lengths.get(row['chrom'], 0)
        
        if contig_len == 0:
            dropped_missing += 1
            continue
        
        # Check if promoter region fits within contig bounds
        if row['start'] < 0 or row['end'] > contig_len:
            dropped_edge += 1
            continue
        
        # Check if promoter is at least min_upstream bp
        if (row['end'] - row['start']) < min_upstream:
            dropped_short += 1
            continue
        
        valid_promoters.append(row)
    
    filtered_prom = pd.DataFrame(valid_promoters)
    
    if output_path:
        filtered_prom.to_csv(output_path, sep="\t", header=False, index=False)
    
    log(f"Promoters: {len(prom):,} total")
    log(f"  Retained: {len(filtered_prom):,} with ≥{min_upstream}bp upstream space")
    log(f"  Dropped: {dropped_edge:,} (edge of contig), {dropped_short:,} (too short), {dropped_missing:,} (contig not found)")
    
    return filtered_prom

def find_promoter_cpg_overlaps(promoters_bed, cpg_summary, high_meth_thresh, output_path):
    """Find promoters overlapping highly methylated CpG islands using % overlap"""
    if isinstance(promoters_bed, str):
        prom = pd.read_csv(promoters_bed, sep="\t", header=None,
                           names=["chrom", "start", "end", "gene_id", "score", "strand"])
    else:
        prom = promoters_bed
    
    cpg = pd.read_csv(cpg_summary, sep="\t")
    cpg_high = cpg[(cpg["n_sites"] >= 1) & (cpg["max_percent"] >= high_meth_thresh)]
    
    log(f"High-methyl CpG islands (≥{high_meth_thresh}%): {len(cpg_high):,}")
    
    overlaps = []
    for chrom in sorted(set(prom["chrom"]).intersection(cpg_high["chrom"])):
        p_chr = prom[prom["chrom"] == chrom]
        c_chr = cpg_high[cpg_high["chrom"] == chrom]
        
        merged = p_chr.merge(c_chr, on="chrom", suffixes=("_prom", "_cpg"))
        
        # Calculate overlap
        overlap_start = np.maximum(merged["start_prom"], merged["start_cpg"])
        overlap_end = np.minimum(merged["end_prom"], merged["end_cpg"])
        overlap_len = np.maximum(0, overlap_end - overlap_start)
        
        cpg_len = merged["end_cpg"] - merged["start_cpg"]
        prom_len = merged["end_prom"] - merged["start_prom"]
        
        merged["overlap_bp"] = overlap_len
        merged["overlap_pct_cpg"] = (overlap_len / cpg_len) * 100
        merged["overlap_pct_prom"] = (overlap_len / prom_len) * 100
        
        # Filter by >25% overlap
        mask = (merged["overlap_pct_cpg"] > 25) | (merged["overlap_pct_prom"] > 25)
        overlaps.append(merged[mask])
    
    if overlaps:
        ov = pd.concat(overlaps, ignore_index=True)
        ov.to_csv(output_path, sep="\t", index=False)
        log(f"Promoter–CpG overlaps (>25%): {len(ov):,} → {output_path}")
        return ov
    
    log("No promoter–CpG overlaps found")
    return pd.DataFrame()

def summarize_genes_with_methylated_promoters(overlap_df, output_path):
    """Summarize genes with highly methylated promoters"""
    gene_summary = (
        overlap_df.groupby("gene_id")
        .agg(
            n_high_methyl_cpg=("start_cpg", "nunique"),
            chrom=("chrom", "first"),
            strand=("strand", "first"),
            prom_start=("start_prom", "min"),
            prom_end=("end_prom", "max"),
            mean_overlap_pct=("overlap_pct_cpg", "mean")
        )
        .reset_index()
        .sort_values("n_high_methyl_cpg", ascending=False)
    )
    
    gene_summary.to_csv(output_path, sep="\t", index=False)
    log(f"Genes with methylated promoters: {len(gene_summary):,} → {output_path}")
    
    return gene_summary

def identify_unmethylated_genes(all_genes_bed, methylated_genes_df, output_path):
    """Identify genes without methylated promoters"""
    all_genes = pd.read_csv(all_genes_bed, sep="\t", header=None,
                            names=["chrom", "start", "end", "gene_id", "score", "strand"])
    
    methylated_gene_ids = set(methylated_genes_df["gene_id"])
    unmethylated = all_genes[~all_genes["gene_id"].isin(methylated_gene_ids)]
    
    unmethylated.to_csv(output_path, sep="\t", index=False)
    log(f"Unmethylated genes: {len(unmethylated):,} → {output_path}")
    
    return unmethylated

# ============================================================================
# EXON/INTRON ANALYSIS
# ============================================================================

def extract_exons_introns_from_gff(gff_path, output_exons, output_introns):
    """Extract exon and intron coordinates from GFF file"""
    exons = []
    genes = {}
    
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand = parts[0], parts[1], parts[2], int(parts[3]), int(parts[4]), parts[5], parts[6]
            
            if feature == 'gene':
                gene_id = parts[8].split('ID=')[1].split(';')[0] if 'ID=' in parts[8] else None
                if gene_id:
                    genes[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand, 'exons': []}
            
            elif feature in ['exon', 'CDS']:
                parent = parts[8].split('Parent=')[1].split(';')[0] if 'Parent=' in parts[8] else None
                exons.append({'chrom': chrom, 'start': start, 'end': end, 'strand': strand, 'gene_id': parent})
                if parent and parent in genes:
                    genes[parent]['exons'].append((start, end))
    
    # Write exons
    exon_df = pd.DataFrame(exons)
    exon_df.to_csv(output_exons, sep='\t', index=False, header=False)
    log(f"Exons extracted: {len(exon_df):,} → {output_exons}")
    
    # Calculate introns
    introns = []
    for gene_id, gene_data in genes.items():
        if len(gene_data['exons']) < 2:
            continue
        sorted_exons = sorted(gene_data['exons'])
        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i][1] + 1
            intron_end = sorted_exons[i + 1][0] - 1
            if intron_end > intron_start:
                introns.append({
                    'chrom': gene_data['chrom'],
                    'start': intron_start,
                    'end': intron_end,
                    'gene_id': gene_id,
                    'strand': gene_data['strand']
                })
    
    intron_df = pd.DataFrame(introns)
    intron_df.to_csv(output_introns, sep='\t', index=False, header=False)
    log(f"Introns calculated: {len(intron_df):,} → {output_introns}")
    
    return exon_df, intron_df

def analyze_cpg_in_features(feature_bed, cpg_summary, high_meth_thresh, low_meth_thresh, output_high, output_low):
    """Analyze CpG islands within exons or introns"""
    features = pd.read_csv(feature_bed, sep='\t', header=None,
                          names=['chrom', 'start', 'end', 'gene_id', 'strand'])
    cpg = pd.read_csv(cpg_summary, sep='\t')
    
    overlaps = []
    for chrom in sorted(set(features['chrom']).intersection(cpg['chrom'])):
        f_chr = features[features['chrom'] == chrom]
        c_chr = cpg[cpg['chrom'] == chrom]
        
        merged = f_chr.merge(c_chr, on='chrom', suffixes=('_feat', '_cpg'))
        overlap_start = np.maximum(merged['start_feat'], merged['start_cpg'])
        overlap_end = np.minimum(merged['end_feat'], merged['end_cpg'])
        overlap_len = np.maximum(0, overlap_end - overlap_start)
        
        merged['overlap_bp'] = overlap_len
        merged = merged[merged['overlap_bp'] > 0]
        overlaps.append(merged)
    
    if not overlaps:
        log("No CpG-feature overlaps found")
        return pd.DataFrame(), pd.DataFrame()
    
    all_overlaps = pd.concat(overlaps, ignore_index=True)
    
    # Classify by methylation
    high_meth = all_overlaps[all_overlaps['max_percent'] >= high_meth_thresh]
    low_meth = all_overlaps[all_overlaps['max_percent'] < low_meth_thresh]
    
    high_meth.to_csv(output_high, sep='\t', index=False)
    low_meth.to_csv(output_low, sep='\t', index=False)
    
    log(f"High methylation: {len(high_meth):,} → {output_high}")
    log(f"Low/no methylation: {len(low_meth):,} → {output_low}")
    
    return high_meth, low_meth

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def run_pipeline():
    """Execute the complete analysis pipeline"""
    cfg = Config()
    ensure_dir(cfg.OUTDIR)
    ensure_dir(cfg.RESULTS_DIR)
    
    log("=== STEP 1: Contig Length Analysis ===")
    contig_df, filtered_contigs = analyze_contig_lengths(
        cfg.CONTIG_LENGTHS_CSV, cfg.MIN_CONTIG_LENGTH
    )
    
    log("\n=== STEP 2: Filter Methylation Data ===")
    meth_filtered = filter_methylation_data(
        cfg.BEDMETHYL_FILE, cfg.METHYLATION_THRESHOLD, cfg.COVERAGE_THRESHOLD
    )
    
    log("\n=== STEP 3: Create Contig Allowlist ===")
    allowlist_path = f"{cfg.OUTDIR}/allowlist_contigs.txt"
    allowset = create_allowlist(meth_filtered, allowlist_path)
    
    log("\n=== STEP 4: Filter Assembly ===")
    trimmed_assembly = f"{cfg.OUTDIR}/trimmed_assembly_{cfg.MIN_CONTIG_LENGTH//1000}kb.fa"
    summary_df = filter_assembly(cfg.ASSEMBLY_FASTA, allowset, 
                                  cfg.MIN_CONTIG_LENGTH, trimmed_assembly)
    summary_df.to_csv(f"{cfg.OUTDIR}/retained_contigs_summary.tsv", sep="\t", index=False)
    
    log("\n=== STEP 5: Filter CpG Islands ===")
    cpg_trimmed = f"{cfg.OUTDIR}/cpg_islands_trimmed.bed"
    cpg_filtered = filter_cpg_islands(cfg.CPG_ISLANDS_BED, allowset, cpg_trimmed)
    
    log("\n=== STEP 6: Summarize Island Methylation ===")
    island_summary = f"{cfg.OUTDIR}/cpg_islands_with_methylation_summary.tsv"
    meth_summary = summarize_island_methylation(cpg_trimmed, cfg.BEDMETHYL_FILE, island_summary)
    
    log("\n=== STEP 6b: Plot CpG Island Distributions ===")
    plot_cpg_island_distributions(island_summary, cfg.OUTDIR)
    
    log("\n=== STEP 7: Filter Promoters by Upstream Space ===")
    filtered_promoters_bed = f"{cfg.OUTDIR}/gene_promoters_filtered.bed"
    filtered_promoters = filter_promoters_by_upstream_space(
        cfg.GENE_PROMOTERS_BED, trimmed_assembly, cfg.PROMOTER_UPSTREAM, filtered_promoters_bed
    )
    
    log("\n=== STEP 8: Find Promoter–CpG Overlaps (>25% overlap) ===")
    overlap_file = f"{cfg.RESULTS_DIR}/promoter_high_methyl_cpg_overlap.tsv"
    overlaps = find_promoter_cpg_overlaps(filtered_promoters, island_summary,
                                          cfg.HIGH_METHYLATION_THRESHOLD, overlap_file)
    
    if not overlaps.empty:
        log("\n=== STEP 9: Summarize Genes ===")
        gene_summary_file = f"{cfg.RESULTS_DIR}/genes_with_high_methylated_promoters.tsv"
        gene_summary = summarize_genes_with_methylated_promoters(overlaps, gene_summary_file)
        print("\n[TOP GENES BY PROMOTER METHYLATION]")
        print(gene_summary.head(20).to_string(index=False))
        
        log("\n=== STEP 10: Identify Unmethylated Genes ===")
        unmethylated_file = f"{cfg.RESULTS_DIR}/genes_without_methylated_promoters.tsv"
        unmethylated_genes = identify_unmethylated_genes(filtered_promoters_bed, gene_summary, unmethylated_file)
        print(f"\nMethylated genes: {len(gene_summary):,}")
        print(f"Unmethylated genes: {len(unmethylated_genes):,}")
    
    # Optional: Exon/Intron analysis
    if os.path.exists(cfg.GENE_GFF):
        log("\n=== STEP 11: Extract Exons and Introns ===")
        exons_bed = f"{cfg.OUTDIR}/exons.bed"
        introns_bed = f"{cfg.OUTDIR}/introns.bed"
        extract_exons_introns_from_gff(cfg.GENE_GFF, exons_bed, introns_bed)
        
        log("\n=== STEP 12: Analyze CpG in Exons ===")
        exon_high = f"{cfg.RESULTS_DIR}/exons_high_methylation.tsv"
        exon_low = f"{cfg.RESULTS_DIR}/exons_low_methylation.tsv"
        analyze_cpg_in_features(exons_bed, island_summary, 
                               cfg.HIGH_METHYLATION_THRESHOLD, cfg.LOW_METHYLATION_THRESHOLD,
                               exon_high, exon_low)
        
        log("\n=== STEP 13: Analyze CpG in Introns ===")
        intron_high = f"{cfg.RESULTS_DIR}/introns_high_methylation.tsv"
        intron_low = f"{cfg.RESULTS_DIR}/introns_low_methylation.tsv"
        analyze_cpg_in_features(introns_bed, island_summary,
                               cfg.HIGH_METHYLATION_THRESHOLD, cfg.LOW_METHYLATION_THRESHOLD,
                               intron_high, intron_low)
    
    log("\n=== PIPELINE COMPLETE ===")

if __name__ == "__main__":
    run_pipeline()
