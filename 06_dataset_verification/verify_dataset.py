#!/usr/bin/env python3
"""
Verify the final dataset to ensure quality and balance.
"""
import os
import pandas as pd
import numpy as np
import argparse
from collections import Counter
import re
import matplotlib.pyplot as plt
from Bio.SeqUtils import gc_fraction

def parse_args():
    parser = argparse.ArgumentParser(description='Verify dataset quality')
    parser.add_argument('--dataset_dir', type=str, default='final_dataset',
                        help='Directory containing the final datasets')
    parser.add_argument('--output_dir', type=str, default='verification_results',
                        help='Output directory for verification results')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load datasets
    datasets = {}
    for split in ['train', 'val', 'test']:
        file_path = os.path.join(args.dataset_dir, f"{split}.csv")
        if os.path.exists(file_path):
            datasets[split] = pd.read_csv(file_path)
    
    # Verify each dataset
    verification_results = {}
    for split, data in datasets.items():
        print(f"Verifying {split} dataset...")
        verification_results[split] = verify_dataset(data, split, args.output_dir)
    
    # Generate summary report
    generate_summary(verification_results, args.output_dir)
    
    print("Verification complete!")

def verify_dataset(data, split, output_dir):
    """
    Verify a dataset and return quality metrics.
    """
    results = {}
    
    # Check class balance
    bacteria_count = sum(data['label'] == 0)
    phage_count = sum(data['label'] == 1)
    total = len(data)
    
    results['class_balance'] = {
        'bacteria': bacteria_count,
        'phage': phage_count,
        'total': total,
        'bacteria_percent': bacteria_count / total * 100,
        'phage_percent': phage_count / total * 100
    }
    
    # Check sequence lengths
    bacteria_lengths = data[data['label'] == 0]['sequence'].apply(len)
    phage_lengths = data[data['label'] == 1]['sequence'].apply(len)
    
    results['sequence_lengths'] = {
        'bacteria_min': bacteria_lengths.min(),
        'bacteria_max': bacteria_lengths.max(),
        'bacteria_mean': bacteria_lengths.mean(),
        'phage_min': phage_lengths.min(),
        'phage_max': phage_lengths.max(),
        'phage_mean': phage_lengths.mean()
    }
    
    # Check for ambiguous nucleotides
    def count_ambiguous(seq):
        return len(re.findall('[^ACGT]', seq.upper()))
    
    bacteria_ambiguous = data[data['label'] == 0]['sequence'].apply(count_ambiguous).sum()
    phage_ambiguous = data[data['label'] == 1]['sequence'].apply(count_ambiguous).sum()
    
    results['ambiguous_nucleotides'] = {
        'bacteria': bacteria_ambiguous,
        'phage': phage_ambiguous,
        'total': bacteria_ambiguous + phage_ambiguous
    }
    
    # Calculate GC content
    bacteria_gc = data[data['label'] == 0]['sequence'].apply(gc_fraction)
    phage_gc = data[data['label'] == 1]['sequence'].apply(gc_fraction)
    
    results['gc_content'] = {
        'bacteria_mean': bacteria_gc.mean(),
        'bacteria_std': bacteria_gc.std(),
        'phage_mean': phage_gc.mean(),
        'phage_std': phage_gc.std()
    }
    
    # Plot GC content distribution
    plt.figure(figsize=(10, 6))
    plt.hist(bacteria_gc, alpha=0.5, bins=50, label='Bacteria')
    plt.hist(phage_gc, alpha=0.5, bins=50, label='Phage')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Count')
    plt.legend()
    plt.title(f'GC Content Distribution - {split.capitalize()} Dataset')
    plt.savefig(os.path.join(output_dir, f"{split}_gc_content.png"))
    plt.close()
    
    # Calculate nucleotide frequencies
    bacteria_seqs = data[data['label'] == 0]['sequence'].str.upper()
    phage_seqs = data[data['label'] == 1]['sequence'].str.upper()
    
    bacteria_nuc_count = Counter(''.join(bacteria_seqs))
    phage_nuc_count = Counter(''.join(phage_seqs))
    
    bacteria_total = sum(bacteria_nuc_count.values())
    phage_total = sum(phage_nuc_count.values())
    
    results['nucleotide_frequencies'] = {
        'bacteria': {nuc: count / bacteria_total * 100 for nuc, count in bacteria_nuc_count.items()},
        'phage': {nuc: count / phage_total * 100 for nuc, count in phage_nuc_count.items()}
    }
    
    # Plot nucleotide frequencies
    nucs = ['A', 'C', 'G', 'T']
    bacteria_freqs = [results['nucleotide_frequencies']['bacteria'].get(nuc, 0) for nuc in nucs]
    phage_freqs = [results['nucleotide_frequencies']['phage'].get(nuc, 0) for nuc in nucs]
    
    plt.figure(figsize=(8, 6))
    x = np.arange(len(nucs))
    width = 0.35
    
    plt.bar(x - width/2, bacteria_freqs, width, label='Bacteria')
    plt.bar(x + width/2, phage_freqs, width, label='Phage')
    
    plt.xlabel('Nucleotide')
    plt.ylabel('Frequency (%)')
    plt.title(f'Nucleotide Frequencies - {split.capitalize()} Dataset')
    plt.xticks(x, nucs)
    plt.legend()
    
    plt.savefig(os.path.join(output_dir, f"{split}_nucleotide_freq.png"))
    plt.close()
    
    return results

def generate_summary(verification_results, output_dir):
    """
    Generate a summary report of verification results.
    """
    summary_file = os.path.join(output_dir, "verification_summary.txt")
    
    with open(summary_file, 'w') as f:
        f.write("Dataset Verification Summary\n")
        f.write("===========================\n\n")
        
        for split, results in verification_results.items():
            f.write(f"{split.upper()} DATASET\n")
            f.write("-" * len(f"{split.upper()} DATASET") + "\n\n")
            
            # Class balance
            f.write("Class Balance:\n")
            f.write(f"  Bacteria: {results['class_balance']['bacteria']} ({results['class_balance']['bacteria_percent']:.2f}%)\n")
            f.write(f"  Phage: {results['class_balance']['phage']} ({results['class_balance']['phage_percent']:.2f}%)\n")
            f.write(f"  Total: {results['class_balance']['total']}\n\n")
            
            # Sequence lengths
            f.write("Sequence Lengths:\n")
            f.write(f"  Bacteria: min={results['sequence_lengths']['bacteria_min']}, max={results['sequence_lengths']['bacteria_max']}, mean={results['sequence_lengths']['bacteria_mean']:.2f}\n")
            f.write(f"  Phage: min={results['sequence_lengths']['phage_min']}, max={results['sequence_lengths']['phage_max']}, mean={results['sequence_lengths']['phage_mean']:.2f}\n\n")
            
            # Ambiguous nucleotides
            f.write("Ambiguous Nucleotides:\n")
            f.write(f"  Bacteria: {results['ambiguous_nucleotides']['bacteria']}\n")
            f.write(f"  Phage: {results['ambiguous_nucleotides']['phage']}\n")
            f.write(f"  Total: {results['ambiguous_nucleotides']['total']}\n\n")
            
            # GC content
            f.write("GC Content:\n")
            f.write(f"  Bacteria: {results['gc_content']['bacteria_mean']:.2f}% ± {results['gc_content']['bacteria_std']:.2f}%\n")
            f.write(f"  Phage: {results['gc_content']['phage_mean']:.2f}% ± {results['gc_content']['phage_std']:.2f}%\n\n")
            
            # Nucleotide frequencies
            f.write("Nucleotide Frequencies:\n")
            f.write("  Bacteria:\n")
            for nuc, freq in sorted(results['nucleotide_frequencies']['bacteria'].items()):
                if nuc in 'ACGT':
                    f.write(f"    {nuc}: {freq:.2f}%\n")
            
            f.write("  Phage:\n")
            for nuc, freq in sorted(results['nucleotide_frequencies']['phage'].items()):
                if nuc in 'ACGT':
                    f.write(f"    {nuc}: {freq:.2f}%\n")
            
            f.write("\n")
    
    print(f"Summary report saved to {summary_file}")

if __name__ == "__main__":
    main()
