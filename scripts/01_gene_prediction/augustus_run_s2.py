#!/usr/bin/env python3
"""AUGUSTUS batch processing script for gene prediction."""

import subprocess
from pathlib import Path


# Configuration
BASE_DIR = Path('/Users/philip_koutsaftis/Library/CloudStorage/GoogleDrive-philipkoutsaftis@gmail.com/My Drive/2025_SnrProj_Dir')
INPUT_BASE = BASE_DIR / 'Data/contigs'
OUTPUT_BASE = BASE_DIR / 'results/augustus_output'
BATCH_RANGE = range(11, 17)  # batches 11-16
SPECIES = 'fly'


def run_augustus(fasta_path, output_path, species='fly'):
    """Run AUGUSTUS on a single FASTA file."""
    cmd = [
        'augustus',
        '--species', species,
        '--gff3', 'on',
        '--UTR', 'on',
        '--print_utr', 'on',
        '--protein', 'on',
        '--outfile', str(output_path),
        str(fasta_path)
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def process_batch(batch_num, input_base, output_base, species='fly'):
    """Process all FASTA files in a batch."""
    input_batch = input_base / f'contigs_batch_{batch_num}'
    output_batch = output_base / f'augustus_batch{batch_num}'
    output_batch.mkdir(parents=True, exist_ok=True)
    
    fasta_files = list(input_batch.glob('contig_*/contig_*.fasta'))
    
    if not fasta_files:
        print(f'⚠️  No FASTA files in batch {batch_num}')
        return 0
    
    print(f'📦 Batch {batch_num}: {len(fasta_files)} files')
    processed = 0
    
    for fasta in fasta_files:
        output = output_batch / f'{fasta.stem}_augustus.gff'
        
        if output.exists():
            continue
        
        try:
            run_augustus(fasta, output, species)
            processed += 1
            print(f'✅ {fasta.stem}')
        except subprocess.CalledProcessError as e:
            print(f'❌ {fasta.stem}: {e}')
    
    return processed


if __name__ == '__main__':
    print('🚀 Starting AUGUSTUS batch processing (batches 11–16)...')
    total_processed = 0
    
    for batch in BATCH_RANGE:
        total_processed += process_batch(batch, INPUT_BASE, OUTPUT_BASE, SPECIES)
    
    print(f'🎉 Completed! Processed {total_processed} files')
