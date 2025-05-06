#!/usr/bin/env python3
"""
Create clusters of genomes to ensure no overlap between train, val, and test sets.
"""
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
from sklearn.cluster import DBSCAN
import pickle
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Create genome clusters')
    parser.add_argument('--bacteria_fasta', type=str, default='clean_data/bacteria_clean.fasta',
                        help='Clean bacteria FASTA file')
    parser.add_argument('--phage_fasta', type=str, default='clean_data/phage_clean.fasta',
                        help='Clean phage FASTA file')
    parser.add_argument('--bacteria_metadata', type=str, default='gtdb_data/metadata/bac120_metadata.tsv',
                        help='Bacteria metadata file from GTDB')
    parser.add_argument('--phage_metadata', type=str, default='phagescope_data/combined_metadata.csv',
                        help='Combined phage metadata file')
    parser.add_argument('--output_dir', type=str, default='clusters',
                        help='Output directory for cluster information')
    parser.add_argument('--n_clusters', type=int, default=3,
                        help='Target number of clusters (for train/val/test)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process bacteria
    print("Creating bacteria clusters...")
    bacteria_clusters = create_bacteria_clusters(args.bacteria_fasta, args.bacteria_metadata, args.n_clusters)
    
    # Process phage
    print("Creating phage clusters...")
    phage_clusters = create_phage_clusters(args.phage_fasta, args.phage_metadata, args.n_clusters)
    
    # Save cluster information
    save_clusters(bacteria_clusters, phage_clusters, args.output_dir)
    
    print("Clustering complete!")

def create_bacteria_clusters(fasta_file, metadata_file, n_clusters):
    """
    Create clusters of bacterial genomes based on taxonomy.
    """
    # Read metadata
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # Create mapping from sequence ID to genome
    seq_to_genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract genome accession from record ID
        genome = record.id.split('|')[0] if '|' in record.id else record.id.split('.')[0]
        seq_to_genome[record.id] = genome
    
    # Group by taxonomic information
    if 'phylum' in metadata.columns:
        # Create clusters based on phylum
        genome_to_cluster = {}
        phyla = metadata['phylum'].unique()
        
        # Assign each phylum to a cluster (distribute evenly)
        phylum_to_cluster = {}
        for i, phylum in enumerate(phyla):
            phylum_to_cluster[phylum] = i % n_clusters
        
        # Assign genomes to clusters
        for _, row in metadata.iterrows():
            genome = row['accession']
            phylum = row['phylum']
            genome_to_cluster[genome] = phylum_to_cluster[phylum]
    else:
        # If no taxonomy info, use random assignment
        genomes = list(set(seq_to_genome.values()))
        genome_to_cluster = {}
        for i, genome in enumerate(genomes):
            genome_to_cluster[genome] = i % n_clusters
    
    # Map sequences to clusters
    seq_clusters = {}
    for seq_id, genome in seq_to_genome.items():
        if genome in genome_to_cluster:
            seq_clusters[seq_id] = genome_to_cluster[genome]
        else:
            # If no cluster for genome, assign randomly
            seq_clusters[seq_id] = np.random.randint(0, n_clusters)
    
    # Group by cluster
    clusters = defaultdict(list)
    for seq_id, cluster in seq_clusters.items():
        clusters[cluster].append(seq_id)
    
    return dict(clusters)

def create_phage_clusters(fasta_file, metadata_file, n_clusters):
    """
    Create clusters of phage genomes based on host information.
    """
    # Read metadata if available
    try:
        metadata = pd.read_csv(metadata_file)
        has_metadata = True
    except:
        has_metadata = False
    
    if has_metadata and 'host' in metadata.columns:
        # Create mapping from phage ID to host
        phage_to_host = {}
        for _, row in metadata.iterrows():
            phage_to_host[row['phage_id']] = row['host']
        
        # Create clusters based on host
        host_to_cluster = {}
        hosts = metadata['host'].unique()
        
        # Assign each host to a cluster (distribute evenly)
        for i, host in enumerate(hosts):
            host_to_cluster[host] = i % n_clusters
        
        # Assign phages to clusters based on host
        seq_clusters = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            phage_id = record.id
            if phage_id in phage_to_host and phage_to_host[phage_id] in host_to_cluster:
                seq_clusters[phage_id] = host_to_cluster[phage_to_host[phage_id]]
            else:
                # If no host info, assign randomly
                seq_clusters[phage_id] = np.random.randint(0, n_clusters)
    else:
        # If no metadata or host info, use random assignment
        seq_clusters = {}
        for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
            seq_clusters[record.id] = i % n_clusters
    
    # Group by cluster
    clusters = defaultdict(list)
    for seq_id, cluster in seq_clusters.items():
        clusters[cluster].append(seq_id)
    
    return dict(clusters)

def save_clusters(bacteria_clusters, phage_clusters, output_dir):
    """
    Save cluster information.
    """
    # Save as pickle for programmatic use
    with open(os.path.join(output_dir, "bacteria_clusters.pkl"), "wb") as f:
        pickle.dump(bacteria_clusters, f)
    
    with open(os.path.join(output_dir, "phage_clusters.pkl"), "wb") as f:
        pickle.dump(phage_clusters, f)
    
    # Also save as text for human readability
    with open(os.path.join(output_dir, "bacteria_clusters.txt"), "w") as f:
        for cluster, sequences in bacteria_clusters.items():
            f.write(f"Cluster {cluster}: {len(sequences)} sequences\n")
            for seq in sequences[:10]:  # Show only first 10 for brevity
                f.write(f"  {seq}\n")
            if len(sequences) > 10:
                f.write(f"  ... and {len(sequences)-10} more\n")
    
    with open(os.path.join(output_dir, "phage_clusters.txt"), "w") as f:
        for cluster, sequences in phage_clusters.items():
            f.write(f"Cluster {cluster}: {len(sequences)} sequences\n")
            for seq in sequences[:10]:  # Show only first 10 for brevity
                f.write(f"  {seq}\n")
            if len(sequences) > 10:
                f.write(f"  ... and {len(sequences)-10} more\n")

if __name__ == "__main__":
    main()
