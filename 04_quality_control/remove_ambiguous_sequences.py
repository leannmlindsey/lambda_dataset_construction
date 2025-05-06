#!/usr/bin/env python3
"""
Remove sequences containing ambiguous nucleotides (N's).
"""
import os
import re
from Bio import SeqIO
import argparse
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description='Remove sequences with ambiguous nucleotides')
    parser.add_argument('--bacteria_fasta', type=str, default='filtered_data/bacteria_filtered.fasta',
                        help='Filtered bacteria FASTA file')
    parser.add_argument('--phage_fasta', type=str, default='filtered_data/phage_filtered.fasta',
                        help='Filtered phage FASTA file')
    parser.add_argument('--output_dir', type=str, default='clean_data',
                        help='Output directory for cleaned sequences')
    parser.add_argument('--max_n_percent', type=float, default=0.0,
                        help='Maximum allowed percentage of N bases (0.0 = no Ns allowed)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output files
    bacteria_output = os.path.join(args.output_dir, "bacteria_clean.fasta")
    phage_output = os.path.join(args.output_dir, "phage_clean.fasta")
    
    # Process bacteria sequences
    print("Processing bacteria sequences...")
    process_sequences(args.bacteria_fasta, bacteria_output, args.max_n_percent)
    
    # Process phage sequences
    print("Processing phage sequences...")
    process_sequences(args.phage_fasta, phage_output, args.max_n_percent)
    
    print("Cleaning complete!")

def process_sequences(input_fasta, output_fasta, max_n_percent):
    """
    Process sequences to remove those with ambiguous nucleotides.
    """
    cleaned_records = []
    total_records = 0
    removed_records = 0
    
    for record in tqdm(SeqIO.parse(input_fasta, "fasta")):
        total_records += 1
        seq_str = str(record.seq).upper()
        
        # Count Ns
        n_count = seq_str.count('N')
        n_percent = n_count / len(seq_str) * 100
        
        # Keep sequence if it has fewer Ns than the threshold
        if n_percent <= max_n_percent:
            cleaned_records.append(record)
        else:
            removed_records += 1
    
    # Write cleaned sequences
    SeqIO.write(cleaned_records, output_fasta, "fasta")
    
    print(f"  Total sequences: {total_records}")
    print(f"  Removed sequences: {removed_records} ({removed_records/total_records*100:.2f}%)")
    print(f"  Remaining sequences: {len(cleaned_records)} ({len(cleaned_records)/total_records*100:.2f}%)")

if __name__ == "__main__":
    main()
