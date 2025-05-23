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
conda env create -f environment.yml
conda activate dna_classification

# BLAST+ tools need to be installed separately
# On Ubuntu/Debian:
# apt-get install ncbi-blast+
# On macOS:
# brew install blast
```

## Usage
Follow these steps to create and validate the dataset:
### Download Data
First download all of the phage sequences and metadata from [PhageScope](https://phagescope.deepomics.org/download) and the bacterial sequences from the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/).

There is no FTP or API for PhageScope so you must download all of the links in the [Phage Meta Data Download](https://phagescope.deepomics.org/download#meta) section as well as the [Phage FASTA File Download](https://phagescope.deepomics.org/download#fasta) section. The scripts below will notify you if you are missing any of the files.
```
# Download phage data
mkdir -P data/phagescope
# Move all of the files downloaded from the above link into this directory.

# Download all representative bacterial sequences from GTDB
mkdir -P data/gtdb 
sh 01_data_download/download_gtdb.sh data/gtdb

# Combine phage metadata
python 01_data_download/combine_phage_metadata.py --input_dir data/phagescope/metadata --output_file data/phagescope/combined_metadata.csv
```
### Create BLAST Databases
Since bacteriophage are present in bacterial genomes, we must remove as many of the prophage sequences as possible to avoid contamination of phage sequences in our labeled bacteria set. Similarly, bacteriophage sometimes incorporate bacterial genes into their genomes via horizontal gene transfer and those also must be removed from the phage dataset.

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
Genomes sometimes contain [IUPAC](https://www.bioinformatics.org/sms/iupac.html) codes for ambiguous nucleotide sequences and because of differences in how the databases use these codes, this can mislead the model, so we will remove these sequences prior to sampling.

```
python 04_quality_control/remove_ambiguous_sequences.py --bacteria_fasta data/filtered_data/bacteria_filtered.fasta --phage_fasta data/filtered_data/phage_filtered.fasta --output_dir data/clean_data --max_n_percent 0.0
```
### Create Dataset Splits
The PhageScope database has clustered each phage genome by sequence similarity and we will select our samples so that we do not have overlap between clusters in our training, validation and test datasets. We will do the same for the bacterial sequences.

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

- Class balance
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
