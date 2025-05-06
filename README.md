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
```

## Usage
Follow these steps to create and validate the dataset:
### Download Data
```
# Download phage data
python 01_data_download/download_phagescope.py --output_dir data/phagescope

# Download bacterial data (only representative genomes)
python 01_data_download/download_gtdb.py --output_dir data/gtdb --only_representative

# Combine phage metadata
python 01_data_download/combine_phage_metadata.py --input_dir data/phagescope/metadata --output_file data/phagescope/combined_metadata.csv
```
### Create BLAST Databases
```
# Create phage BLAST database
bash 02_database_creation/create_phage_blast_db.sh --phage_dir data/phagescope/fna --output_dir data/blast_db --db_name phage_db

# Create bacteria BLAST database
bash 02_database_creation/create_bacteria_blast_db.sh --bacteria_dir data/gtdb/fna --output_dir data/blast_db --db_name bacteria_db
```
### Cross-Filter Sequences
```
# BLAST bacteria against phage
bash 03_cross_filtering/blast_bacteria_against_phage.sh --bacteria_fasta data/blast_db/combined_bacteria.fasta --phage_db data/blast_db/phage_db --output_dir data/blast_results

# Remove phage regions from bacteria
python 03_cross_filtering/remove_phage_from_bacteria.py --blast_results data/blast_results/bacteria_vs_phage.tsv --bacteria_fasta data/blast_db/combined_bacteria.fasta --output_fasta data/filtered_data/bacteria_filtered.fasta

# BLAST phage against bacteria
bash 03_cross_filtering/blast_phage_against_bacteria.sh --phage_fasta data/blast_db/combined_phage.fasta --bacteria_db data/blast_db/bacteria_db --output_dir data/blast_results

# Remove bacteria regions from phage
python 03_cross_filtering/remove_bacteria_from_phage.py --blast_results data/blast_results/phage_vs_bacteria.tsv --phage_fasta data/blast_db/combined_phage.fasta --output_fasta data/filtered_data/phage_filtered.fasta
```

### Remove Ambiguous Sequences
```
python 04_quality_control/remove_ambiguous_sequences.py --bacteria_fasta data/filtered_data/bacteria_filtered.fasta --phage_fasta data/filtered_data/phage_filtered.fasta --output_dir data/clean_data --max_n_percent 0.0
```
### Create Dataset Splits
```
# Create genome clusters
python 05_dataset_creation/create_genome_clusters.py --bacteria_fasta data/clean_data/bacteria_clean.fasta --phage_fasta data/clean_data/phage_clean.fasta --bacteria_metadata data/gtdb/metadata/bac120_metadata.tsv --phage_metadata data/phagescope/combined_metadata.csv --output_dir data/clusters

# Generate train/val/test datasets
python 05_dataset_creation/generate_train_val_test.py --bacteria_fasta data/clean_data/bacteria_clean.fasta --phage_fasta data/clean_data/phage_clean.fasta --bacteria_clusters data/clusters/bacteria_clusters.pkl --phage_clusters data/clusters/phage_clusters.pkl --output_dir data/final_dataset --fragment_length 4000 --train_size 40000 --val_size 10000 --test_size 10000
```
### Verify and Visualize Dataset
```
# Verify dataset quality
python 06_dataset_verification/verify_dataset.py --dataset_dir data/final_dataset --output_dir data/verification_results

# Generate visualizations
python 06_dataset_verification/visualization.py --dataset_dir data/final_dataset --output_dir data/visualizations
```
## Dataset Structure
The final dataset consists of three CSV files:

- train.csv: 80,000 sequences (40,000 bacteria, 40,000 phage)
- val.csv: 20,000 sequences (10,000 bacteria, 10,000 phage)
- test.csv: 20,000 sequences (10,000 bacteria, 10,000 phage)

Each file has the following columns:

- accession: Unique identifier with genome source and coordinates
- sequence: 4,000 bp DNA sequence
- label: 0 for bacteria, 1 for phage


## Verification Results
The verification process checks for:

Class balance
- Sequence length distribution
- Absence of ambiguous nucleotides
- GC content distribution
- Nucleotide frequency patterns
- k-mer distribution via PCA visualization

Detailed results and visualizations can be found in the verification_results and visualizations directories after running the pipeline.

## Citation

@misc{bacteria_phage_classification_dataset,
  author = {LeAnn M. Lindsey, Nicole L. Pershing, Keith Dufault-Thompson, Anisa Habib, Aaron Schindler, June Round, Xiaofang Jiang, W. Zac Stephens, Anne J. Blaschke, Hari Sundar},
  title = {LAMBDA: A Prophage Detection Benchmark for Genomic Language Models},
  year = {2025},
  publisher = {bioRxiv},
  url = {https://github.com/leannmlindsey/lambda_dataset_construction}
}
