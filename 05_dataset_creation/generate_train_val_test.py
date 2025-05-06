#!/usr/bin/env python3
"""
Generate train, validation, and test datasets from clustered sequences.
"""
import os
import pickle
import random
import argparse
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description='Generate datasets')
    parser.add_argument('--bacteria_fasta', type=str, default='clean_data/bacteria_clean.fasta',
                        help='Clean bacteria FASTA file')
    parser.add_argument('--phage_fasta', type=str, default='clean_data/phage_clean.fasta',
                        help='Clean phage FASTA file')
    parser.add_argument('--bacteria_clusters', type=str, default='clusters/bacteria_clusters.pkl',
                        help='Bacteria clusters pickle file')
    parser.add_argument('--phage_clusters', type=str, default='clusters/phage_clusters.pkl',
                        help='Phage clusters pickle file')
    parser.add_argument('--output_dir', type=str, default='final_dataset',
                        help='Output directory for final datasets')
    parser.add_argument('--fragment_length', type=int, default=4000,
                        help='Length of sequence fragments')
    parser.add_argument('--train_size', type=int, default=40000,
                        help='Number of sequences per class in training set')
    parser.add_argument('--val_size', type=int, default=10000,
                        help='Number of sequences per class in validation set')
    parser.add_argument('--test_size', type=int, default=10000,
                        help='Number of sequences per class in test set')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Set random seed
    random.seed(args.seed)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load clusters
    with open(args.bacteria_clusters, 'rb') as f:
        bacteria_clusters = pickle.load(f)
    
    with open(args.phage_clusters, 'rb') as f:
        phage_clusters = pickle.load(f)
    
    # Assign clusters to datasets
    # Assuming 3 clusters: 0=train, 1=val, 2=test
    dataset_mapping = {
        0: "train",
        1: "val",
        2: "test"
    }
    
    # Load and index sequences
    print("Loading and indexing sequences...")
    bacteria_seqs = index_sequences(args.bacteria_fasta)
    phage_seqs = index_sequences(args.phage_fasta)
    
    # Extract fragments
    print("Extracting fragments...")
    bacteria_fragments = extract_fragments(bacteria_seqs, bacteria_clusters, 
                                           dataset_mapping, args.fragment_length)
    phage_fragments = extract_fragments(phage_seqs, phage_clusters, 
                                        dataset_mapping, args.fragment_length)
    
    # Sample and create datasets
    print("Creating datasets...")
    create_datasets(bacteria_fragments, phage_fragments, args.output_dir, 
                    args.train_size, args.val_size, args.test_size)
    
    print("Dataset generation complete!")

def index_sequences(fasta_file):
    """
    Load sequences from FASTA file and index them by ID.
    """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def extract_fragments(sequences, clusters, dataset_mapping, fragment_length):
    """
    Extract fixed-length fragments from sequences based on clusters.
    """
    fragments = {dataset: [] for dataset in dataset_mapping.values()}
    
    for cluster_id, seq_ids in tqdm(clusters.items()):
        dataset = dataset_mapping.get(cluster_id, "train")  # Default to train if cluster not mapped
        
        for seq_id in seq_ids:
            if seq_id in sequences:
                seq = sequences[seq_id]
                
                # Only process sequences that are long enough
                if len(seq) >= fragment_length:
                    # Calculate number of non-overlapping fragments
                    n_fragments = len(seq) // fragment_length
                    
                    # Extract fragments
                    for i in range(n_fragments):
                        start = i * fragment_length
                        fragment = seq[start:start+fragment_length]
                        
                        # Check that fragment has no ambiguous nucleotides
                        if 'N' not in fragment.upper():
                            fragments[dataset].append({
                                'accession': f"{seq_id}_{start}-{start+fragment_length}",
                                'sequence': fragment
                            })
    
    return fragments

def create_datasets(bacteria_fragments, phage_fragments, output_dir, train_size, val_size, test_size):
    """
    Create balanced datasets for training, validation, and testing.
    """
    for dataset in ['train', 'val', 'test']:
        # Determine sample size
        if dataset == 'train':
            size = train_size
        elif dataset == 'val':
            size = val_size
        else:  # test
            size = test_size
        
        # Sample fragments
        bacteria_sample = sample_fragments(bacteria_fragments[dataset], size)
        phage_sample = sample_fragments(phage_fragments[dataset], size)
        
        # Create dataset with labels (0 for bacteria, 1 for phage)
        data = []
        for fragment in bacteria_sample:
            data.append({
                'accession': fragment['accession'],
                'sequence': fragment['sequence'],
                'label': 0
            })
        
        for fragment in phage_sample:
            data.append({
                'accession': fragment['accession'],
                'sequence': fragment['sequence'],
                'label': 1
            })
        
        # Shuffle data
        random.shuffle(data)
        
        # Save to CSV
        output_file = os.path.join(output_dir, f"{dataset}.csv")
        pd.DataFrame(data).to_csv(output_file, index=False)
        
        print(f"Created {dataset} dataset with {len(data)} sequences ({len(bacteria_sample)} bacteria, {len(phage_sample)} phage)")

def sample_fragments(fragments, size):
    """
    Sample a specified number of fragments, or all if fewer available.
    """
    if len(fragments) <= size:
        return fragments
    return random.sample(fragments, size)

if __name__ == "__main__":
    main()
