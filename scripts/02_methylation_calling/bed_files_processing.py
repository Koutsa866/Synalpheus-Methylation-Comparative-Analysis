#!/usr/bin/env python3
"""Process GFF files to BED format and perform genomic analyses."""

import os
import pandas as pd
from pathlib import Path
from glob import glob


BASE_DIR = Path('/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir')
GFF_DIR = BASE_DIR / 'results/augustus_output'
BED_DIR = BASE_DIR / 'results/augustus_gene_bed_files'
STRAND_DIR = BASE_DIR / 'results/strand_aware_bed'


def gff_to_bed(gff_path, bed_path):
    """Convert GFF gene entries to BED format."""
    with open(gff_path) as f, open(bed_path, 'w') as out:
        for line in f:
            if '\tgene\t' in line and not line.startswith('#'):
                parts = line.strip().split('\t')
                out.write(f"{parts[0]}\t{int(parts[3])-1}\t{parts[4]}\tgene\n")


def gff_to_strand_bed(gff_path, out):
    """Convert GFF to strand-aware BED format."""
    with open(gff_path) as f:
        for line in f:
            if '\tgene\t' in line:
                fields = line.strip().split('\t')
                gene_id = fields[8].split(';')[0]
                out.write(f"{fields[0]}\t{int(fields[3])-1}\t{fields[4]}\t{gene_id}\t.\t{fields[6]}\n")


def convert_gff_to_bed():
    """Convert all GFF files to BED format."""
    BED_DIR.mkdir(parents=True, exist_ok=True)
    gff_files = glob(str(GFF_DIR / 'augustus_batch*/*.gff'))
    
    for gff in gff_files:
        bed = BED_DIR / f"{Path(gff).stem}.bed"
        if not bed.exists():
            gff_to_bed(gff, bed)
    
    print(f"✅ Converted {len(gff_files)} GFF files")


def merge_bed_files():
    """Merge all BED files into one."""
    output = BED_DIR / 'merged_augustus_genes.bed'
    bed_files = [f for f in glob(str(BED_DIR / '*.bed')) if 'merged' not in f]
    
    with open(output, 'w') as out:
        for bed in bed_files:
            with open(bed) as f:
                out.writelines(f)
    
    print(f"✅ Merged {len(bed_files)} BED files")


def create_strand_aware_bed():
    """Create strand-aware BED file from all GFF files."""
    STRAND_DIR.mkdir(parents=True, exist_ok=True)
    output = STRAND_DIR / 'all_genes_strand_aware.bed'
    gff_files = glob(str(GFF_DIR / '**/*.gff'), recursive=True)
    
    with open(output, 'w') as out:
        for gff in gff_files:
            gff_to_strand_bed(gff, out)
    
    print(f"✅ Created strand-aware BED from {len(gff_files)} files")


def extract_tss():
    """Extract TSS positions from strand-aware BED."""
    bed = STRAND_DIR / 'all_genes_strand_aware.bed'
    df = pd.read_csv(bed, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand'])
    
    df['tss_start'] = df.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'] - 1, axis=1)
    df['tss_end'] = df['tss_start'] + 1
    
    output = STRAND_DIR / 'gene_tss.bed'
    df[['chrom', 'tss_start', 'tss_end', 'gene_id', 'score', 'strand']].to_csv(
        output, sep='\t', header=False, index=False)
    
    print(f"✅ Extracted {len(df)} TSS positions")


def extract_promoters():
    """Extract promoter regions (1kb upstream of TSS)."""
    tss = STRAND_DIR / 'gene_tss.bed'
    df = pd.read_csv(tss, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand'])
    
    prom_start = []
    prom_end = []
    for _, row in df.iterrows():
        if row['strand'] == '+':
            prom_start.append(max(0, row['start'] - 1000))
            prom_end.append(row['start'])
        else:
            prom_start.append(row['end'])
            prom_end.append(row['end'] + 1000)
    
    df['prom_start'] = prom_start
    df['prom_end'] = prom_end
    
    output = STRAND_DIR / 'gene_promoters.bed'
    df[['chrom', 'prom_start', 'prom_end', 'gene_id', 'score', 'strand']].to_csv(
        output, sep='\t', header=False, index=False)
    
    print(f"✅ Extracted {len(df)} promoter regions")


if __name__ == '__main__':
    convert_gff_to_bed()
    merge_bed_files()
    create_strand_aware_bed()
    extract_tss()
    extract_promoters()
