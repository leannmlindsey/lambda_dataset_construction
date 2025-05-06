#!/usr/bin/env python3
"""
Remove bacteria-like regions from phage sequences based on BLAST results.
"""
import os
import pandas as pd
from Bio import SeqIO
import argparse
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description='Remove bacteria-like regions from phage')
    parser.add_argument('--blast_results', type=str, default='blast_results/phage_vs_bacteria.tsv',
                        help='BLAST results file (phage vs bacteria)')
    parser.add_argument('--phage_fasta', type=str, default='blast_db/combined_phage.fasta',
                        help='Combined phage FASTA file')
    parser.add_argument('--output_fasta', type=str, default='filtered_data/phage_filtered.fasta',
                        help='Output filtered phage FASTA file')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(os.path.dirname(args.output_fasta), exist_ok=True)
    
    # Read BLAST results
    print("Reading BLAST results...")
    blast_cols = [
        'qseqid', 'sseqid', 'pident', 'length', 'qlen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    blast_results = pd.read_csv(args.blast_results, sep='\t', names=blast_cols)
    
    # Group BLAST hits by query sequence
    print("Processing BLAST hits...")
    seq_hits = {}
    for _, row in blast_results.iterrows():
        qseqid = row['qseqid']
        if qseqid not in seq_hits:
            seq_hits[qseqid] = []
        seq_hits[qseqid].append((row['qstart'], row['qend']))
    
    # Merge overlapping regions
    for qseqid in seq_hits:
        # Sort by start position
        regions = sorted(seq_hits[qseqid])
        merged = []
        for region in regions:
            # If merged is empty or current region doesn't overlap with previous
            if not merged or region[0] > merged[-1][1]:
                merged.append(region)
            # If there is overlap, merge with the previous region
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], region[1]))
        seq_hits[qseqid] = merged
    
    # Filter sequences
    print("Filtering phage sequences...")
    filtered_records = []
    
    for record in tqdm(SeqIO.parse(args.phage_fasta, "fasta")):
        seq_id = record.id
        seq = record.seq
        
        # If this sequence has bacteria-like regions, remove them
        if seq_id in seq_hits:
            # Create a mask for the sequence (1 = keep, 0 = remove)
            mask = [1] * len(seq)
            
            # Mark regions to remove
            for start, end in seq_hits[seq_id]:
                for i in range(start-1, end):  # Convert to 0-based
                    if i < len(mask):
                        mask[i] = 0
            
            # Filter sequence using mask
            filtered_seq = ''.join([seq[i] for i in range(len(seq)) if mask[i] == 1])
            
            # Only keep if there's enough sequence left
            if len(filtered_seq) >= 4000:  # Minimum sequence length requirement
                record.seq = filtered_seq
                filtered_records.append(record)
        else:
            # No bacteria-like regions found, keep the sequence as is
            filtered_records.append(record)
    
    # Write filtered sequences
    print(f"Writing {len(filtered_records)} filtered sequences...")
    SeqIO.write(filtered_records, args.output_fasta, "fasta")
    
    print("Filtering complete!")

if __name__ == "__main__":
    main()
