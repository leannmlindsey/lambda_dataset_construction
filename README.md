# Bacteria vs. Phage DNA Classification Dataset

This repository provides a comprehensive pipeline for creating, validating, and analyzing a high-quality dataset for bacterial and phage DNA classification tasks. The code ensures robust, reproducible dataset construction that avoids common pitfalls like data leakage and sampling bias.

## Background

Distinguishing between bacterial and phage DNA sequences is a fundamental task in metagenomics and microbial genomics. However, creating reliable datasets for machine learning approaches requires careful consideration of:

1. **Data Leakage**: Avoiding features that trivially distinguish the classes
2. **Cross-Contamination**: Ensuring bacterial sequences are free of phage elements and vice versa
3. **Sequence Quality**: Handling ambiguous nucleotides consistently
4. **Train/Test Separation**: Ensuring no taxonomic overlap between splits

This pipeline addresses all these challenges through a series of carefully designed steps.

## Pipeline Overview

![Pipeline Overview](pipeline_overview.png)

The dataset creation process follows these steps:

1. **Data Download**: Retrieve bacterial genomes from GTDB and phage genomes from PhageScope
2. **Database Creation**: Build BLAST databases for cross-contamination filtering
3. **Cross-Filtering**: Remove phage-like regions from bacteria and bacteria-like regions from phages
4. **Quality Control**: Remove sequences with ambiguous nucleotides (N's)
5. **Dataset Creation**: Generate train/validation/test splits based on taxonomic clustering
6. **Verification**: Validate dataset quality and visualize key characteristics

## Prerequisites

- Python 3.8+
- BioPython
- Pandas, NumPy, Matplotlib, Seaborn
- scikit-learn
- BLAST+ command-line tools

Install the required Python packages:

```bash
# Create and activate conda environment
conda create -n dna_classification python=3.8
conda activate dna_classification

# Install dependencies
pip install -r requirements.txt

# BLAST+ tools need to be installed separately
# On Ubuntu/Debian:
# apt-get install ncbi-blast+
# On macOS:
# brew install blast
