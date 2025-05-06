#!/bin/bash

gtdb_dir=$1

cd $gtdb_dir

echo "Starting download from https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
echo "This will take approximately 250 GB of space after the fna files are decompressed. Please verify that you have enough space before you begin downloading."
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_taxonomy.tsv.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz

mkdir metadata
mkdir taxonomy
mkdir fna

mv bac120_metadata.tsv.gz metadata
mv bac120_taxonomy.tsv.gz taxonomy
mv gtdb_genomes_reps.tar.gz fna

cd metadata
gunzip bac120_metadata.tsv.gz

cd ../taxonomy
gunzip bac120_taxonomy.tsv.gz

cd ../fna
tar -xvzf gtdb_genomes_reps.tar.gz

echo "Download complete!"
