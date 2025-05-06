#!/usr/bin/env python3
"""
Create visualizations to verify dataset quality and biases.
"""
import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import GC
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser(description='Create dataset visualizations')
    parser.add_argument('--dataset_dir', type=str, default='final_dataset',
                        help='Directory containing the final datasets')
    parser.add_argument('--output_dir', type=str, default='visualizations',
                        help='Output directory for visualizations')
    parser.add_argument('--sample_size', type=int, default=1000,
                        help='Number of sequences to sample for k-mer analysis')
    parser.add_argument('--kmer_size', type=int, default=6,
                        help='Size of k-mers for analysis')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load datasets
    train_data = pd.read_csv(os.path.join(args.dataset_dir, "train.csv"))
    val_data = pd.read_csv(os.path.join(args.dataset_dir, "val.csv"))
    test_data = pd.read_csv(os.path.join(args.dataset_dir, "test.csv"))
    
    # 1. Class Distribution
    print("Generating class distribution plots...")
    plot_class_distribution(train_data, val_data, test_data, args.output_dir)
    
    # 2. GC Content Distribution
    print("Generating GC content plots...")
    plot_gc_content(train_data, val_data, test_data, args.output_dir)
    
    # 3. K-mer PCA Visualization
    print("Generating k-mer PCA visualizations...")
    plot_kmer_pca(train_data, args.output_dir, args.sample_size, args.kmer_size)
    
    # 4. Sequence Length Distribution
    print("Generating sequence length distribution...")
    plot_sequence_length(train_data, val_data, test_data, args.output_dir)
    
    # 5. Nucleotide Distribution
    print("Generating nucleotide distribution plots...")
    plot_nucleotide_distribution(train_data, args.output_dir)
    
    print("Visualization complete!")

def plot_class_distribution(train_data, val_data, test_data, output_dir):
    """
    Plot class distribution for all datasets.
    """
    datasets = {
        'Train': train_data,
        'Validation': val_data,
        'Test': test_data
    }
    
    # Prepare data
    labels = []
    counts = []
    dataset_names = []
    
    for name, data in datasets.items():
        for label in [0, 1]:
            count = sum(data['label'] == label)
            labels.append("Bacteria" if label == 0 else "Phage")
            counts.append(count)
            dataset_names.append(name)
    
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        'Dataset': dataset_names,
        'Class': labels,
        'Count': counts
    })
    
    # Plot
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Dataset', y='Count', hue='Class', data=plot_data)
    plt.title('Class Distribution Across Datasets', fontsize=16)
    plt.ylabel('Number of Sequences', fontsize=14)
    plt.xlabel('Dataset', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend(title='Class', fontsize=12)
    
    # Add count labels on top of bars
    for i, p in enumerate(plt.gca().patches):
        plt.gca().annotate(f"{p.get_height():,}", 
                          (p.get_x() + p.get_width() / 2., p.get_height()), 
                          ha = 'center', va = 'bottom', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "class_distribution.png"), dpi=300)
    plt.close()

def plot_gc_content(train_data, val_data, test_data, output_dir):
    """
    Plot GC content distribution for bacteria and phage.
    """
    # Calculate GC content
    train_data['gc'] = train_data['sequence'].apply(GC)
    val_data['gc'] = val_data['sequence'].apply(GC)
    test_data['gc'] = test_data['sequence'].apply(GC)
    
    # Combine datasets
    train_data['dataset'] = 'Train'
    val_data['dataset'] = 'Validation'
    test_data['dataset'] = 'Test'
    
    combined_data = pd.concat([train_data, val_data, test_data])
    combined_data['class'] = combined_data['label'].apply(lambda x: "Bacteria" if x == 0 else "Phage")
    
    # Plot
    plt.figure(figsize=(14, 10))
    
    # Overall distribution
    plt.subplot(2, 2, 1)
    sns.histplot(data=combined_data, x='gc', hue='class', bins=50, kde=True, element='step')
    plt.title('GC Content Distribution - All Datasets', fontsize=14)
    plt.xlabel('GC Content (%)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    
    # Train dataset
    plt.subplot(2, 2, 2)
    sns.histplot(data=combined_data[combined_data['dataset'] == 'Train'], 
                x='gc', hue='class', bins=50, kde=True, element='step')
    plt.title('GC Content - Train Dataset', fontsize=14)
    plt.xlabel('GC Content (%)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    
    # Validation dataset
    plt.subplot(2, 2, 3)
    sns.histplot(data=combined_data[combined_data['dataset'] == 'Validation'], 
                x='gc', hue='class', bins=50, kde=True, element='step')
    plt.title('GC Content - Validation Dataset', fontsize=14)
    plt.xlabel('GC Content (%)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    
    # Test dataset
    plt.subplot(2, 2, 4)
    sns.histplot(data=combined_data[combined_data['dataset'] == 'Test'], 
                x='gc', hue='class', bins=50, kde=True, element='step')
    plt.title('GC Content - Test Dataset', fontsize=14)
    plt.xlabel('GC Content (%)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_distribution.png"), dpi=300)
    plt.close()
    
    # Box plot for GC content comparison
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='dataset', y='gc', hue='class', data=combined_data)
    plt.title('GC Content Comparison Across Datasets', fontsize=16)
    plt.xlabel('Dataset', fontsize=14)
    plt.ylabel('GC Content (%)', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add statistical annotations
    for ds in ['Train', 'Validation', 'Test']:
        bact = combined_data[(combined_data['dataset'] == ds) & (combined_data['class'] == 'Bacteria')]['gc']
        phage = combined_data[(combined_data['dataset'] == ds) & (combined_data['class'] == 'Phage')]['gc']
        
        plt.annotate(f"Bacteria: {bact.mean():.2f}% ± {bact.std():.2f}%", 
                    xy=(0, 0), xytext=(0.1, 0.85 - 0.1*['Train', 'Validation', 'Test'].index(ds)), 
                    xycoords='figure fraction', fontsize=12)
        plt.annotate(f"Phage: {phage.mean():.2f}% ± {phage.std():.2f}%", 
                    xy=(0, 0), xytext=(0.6, 0.85 - 0.1*['Train', 'Validation', 'Test'].index(ds)), 
                    xycoords='figure fraction', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_boxplot.png"), dpi=300)
    plt.close()

def plot_kmer_pca(train_data, output_dir, sample_size, kmer_size):
    """
    Perform PCA on k-mer frequencies and visualize results.
    """
    # Sample data for faster processing
    bacteria_sample = train_data[train_data['label'] == 0].sample(n=min(sample_size, sum(train_data['label'] == 0)))
    phage_sample = train_data[train_data['label'] == 1].sample(n=min(sample_size, sum(train_data['label'] == 1)))

    # Combine samples
    combined_sample = pd.concat([bacteria_sample, phage_sample])

    # Function to convert sequences to k-mers
    def seq_to_kmers(seq, k):
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        return ' '.join(kmers)

    # Convert sequences to k-mer representation
    combined_sample['kmers'] = combined_sample['sequence'].apply(lambda x: seq_to_kmers(x, kmer_size))

    # Vectorize k-mers
    vectorizer = CountVectorizer(max_features=5000)  # Limit features for memory
    X = vectorizer.fit_transform(combined_sample['kmers'])

    # Perform PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X.toarray())

    # Create DataFrame for plotting
    pca_df = pd.DataFrame({
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1],
        'Class': combined_sample['label'].apply(lambda x: "Bacteria" if x == 0 else "Phage")
    })

    # Plot
    plt.figure(figsize=(12, 10))

    # Scatter plot
    sns.scatterplot(x='PC1', y='PC2', hue='Class', data=pca_df, alpha=0.7)
    plt.title(f'{kmer_size}-mer PCA Visualization', fontsize=16)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)', fontsize=14)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)', fontsize=14)
    plt.grid(linestyle='--', alpha=0.7)

    # Add annotation for total explained variance
    total_var = sum(pca.explained_variance_ratio_) * 100
    plt.annotate(f'Total explained variance: {total_var:.2f}%',
                xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{kmer_size}mer_pca.png"), dpi=300)
    plt.close()

    # Find most discriminative k-mers
    feature_names = vectorizer.get_feature_names_out()
    bacteria_freqs = X[combined_sample['label'] == 0].mean(axis=0)
    phage_freqs = X[combined_sample['label'] == 1].mean(axis=0)

    # Calculate ratio of frequencies
    eps = 1e-10  # To avoid division by zero
    bacteria_ratio = bacteria_freqs / (phage_freqs + eps)
    phage_ratio = phage_freqs / (bacteria_freqs + eps)

    # Find top discriminative k-mers
    top_bacteria_idx = np.argsort(bacteria_ratio.A1)[-20:]
    top_phage_idx = np.argsort(phage_ratio.A1)[-20:]

    # Save top discriminative k-mers
    with open(os.path.join(output_dir, f"top_{kmer_size}mers.txt"), 'w') as f:
        f.write(f"Top {kmer_size}-mers enriched in Bacteria:\n")
        for idx in top_bacteria_idx:
            f.write(f"{feature_names[idx]}: {bacteria_ratio.A1[idx]:.2f}x\n")

        f.write(f"\nTop {kmer_size}-mers enriched in Phage:\n")
        for idx in top_phage_idx:
            f.write(f"{feature_names[idx]}: {phage_ratio.A1[idx]:.2f}x\n")

def plot_sequence_length(train_data, val_data, test_data, output_dir):
    """
    Plot sequence length distribution.
    """
    # Calculate sequence lengths
    train_data['length'] = train_data['sequence'].apply(len)
    val_data['length'] = val_data['sequence'].apply(len)
    test_data['length'] = test_data['sequence'].apply(len)

    # Combine datasets
    train_data['dataset'] = 'Train'
    val_data['dataset'] = 'Validation'
    test_data['dataset'] = 'Test'

    combined_data = pd.concat([train_data, val_data, test_data])
    combined_data['class'] = combined_data['label'].apply(lambda x: "Bacteria" if x == 0 else "Phage")

    # Plot
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='dataset', y='length', hue='class', data=combined_data)
    plt.title('Sequence Length Distribution', fontsize=16)
    plt.xlabel('Dataset', fontsize=14)
    plt.ylabel('Sequence Length (bp)', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Add statistical annotations
    for ds in ['Train', 'Validation', 'Test']:
        bact = combined_data[(combined_data['dataset'] == ds) & (combined_data['class'] == 'Bacteria')]['length']
        phage = combined_data[(combined_data['dataset'] == ds) & (combined_data['class'] == 'Phage')]['length']

        plt.annotate(f"Bacteria: {bact.mean():.1f} ± {bact.std():.1f}",
                    xy=(0, 0), xytext=(0.1, 0.85 - 0.1*['Train', 'Validation', 'Test'].index(ds)),
                    xycoords='figure fraction', fontsize=12)
        plt.annotate(f"Phage: {phage.mean():.1f} ± {phage.std():.1f}",
                    xy=(0, 0), xytext=(0.6, 0.85 - 0.1*['Train', 'Validation', 'Test'].index(ds)),
                    xycoords='figure fraction', fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sequence_length_distribution.png"), dpi=300)
    plt.close()

def plot_nucleotide_distribution(train_data, output_dir):
    """
    Plot nucleotide distribution for bacteria and phage.
    """
    # Function to count nucleotides
    def count_nucleotides(seq):
        seq = seq.upper()
        return {
            'A': seq.count('A') / len(seq) * 100,
            'C': seq.count('C') / len(seq) * 100,
            'G': seq.count('G') / len(seq) * 100,
            'T': seq.count('T') / len(seq) * 100
        }

    # Calculate nucleotide frequencies
    bacteria_seqs = train_data[train_data['label'] == 0]['sequence']
    phage_seqs = train_data[train_data['label'] == 1]['sequence']

    bacteria_counts = [count_nucleotides(seq) for seq in bacteria_seqs]
    phage_counts = [count_nucleotides(seq) for seq in phage_seqs]

    # Create DataFrame for plotting
    bacteria_df = pd.DataFrame(bacteria_counts)
    bacteria_df['class'] = 'Bacteria'
    phage_df = pd.DataFrame(phage_counts)
    phage_df['class'] = 'Phage'

    nuc_df = pd.concat([bacteria_df, phage_df])
    nuc_long = pd.melt(nuc_df, id_vars=['class'], value_vars=['A', 'C', 'G', 'T'],
                      var_name='Nucleotide', value_name='Frequency')

    # Plot
    plt.figure(figsize=(14, 10))

    # Boxplot
    plt.subplot(2, 1, 1)
    sns.boxplot(x='Nucleotide', y='Frequency', hue='class', data=nuc_long)
    plt.title('Nucleotide Frequency Distribution', fontsize=16)
    plt.xlabel('Nucleotide', fontsize=14)
    plt.ylabel('Frequency (%)', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Bar plot with averages
    plt.subplot(2, 1, 2)
    nuc_avg = nuc_long.groupby(['class', 'Nucleotide'])['Frequency'].mean().reset_index()
    sns.barplot(x='Nucleotide', y='Frequency', hue='class', data=nuc_avg)
    plt.title('Average Nucleotide Frequencies', fontsize=16)
    plt.xlabel('Nucleotide', fontsize=14)
    plt.ylabel('Average Frequency (%)', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Add values on top of bars
    for i, p in enumerate(plt.gca().patches):
        plt.gca().annotate(f"{p.get_height():.1f}%",
                          (p.get_x() + p.get_width() / 2., p.get_height()),
                          ha = 'center', va = 'bottom', fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "nucleotide_distribution.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
