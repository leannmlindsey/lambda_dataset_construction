#!/usr/bin/env python3
"""
Download bacterial genome data and metadata from GTDB.
"""
import os
import subprocess
import argparse
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import gzip
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description='Download bacterial data from GTDB')
    parser.add_argument('--output_dir', type=str, default='gtdb_data',
                        help='Directory to store downloaded files')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of parallel downloads')
    parser.add_argument('--only_representative', action='store_true',
                        help='Download only representative genomes')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directories
    os.makedirs(f"{args.output_dir}/fna", exist_ok=True)
    os.makedirs(f"{args.output_dir}/metadata", exist_ok=True)
    
    # Download GTDB metadata
    metadata_url = "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv"
    metadata_file = f"{args.output_dir}/metadata/bac120_metadata.tsv"
    
    print(f"Downloading GTDB metadata to {metadata_file}")
    if not os.path.exists(metadata_file):
        subprocess.run(["wget", "-O", metadata_file, metadata_url])
    
    # Read metadata
    print("Reading metadata...")
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # Filter for representative genomes if requested
    if args.only_representative:
        metadata = metadata[metadata['gtdb_representative'] == 't']
    
    print(f"Downloading {len(metadata)} bacterial genomes")
    
    # Function to extract accession from GTDB format
    def extract_accession(gtdb_accession):
        # GTDB format: GB_GCA_000123456.1 or RS_GCF_000123456.1
        # Extract accession part: GCA_000123456.1 or GCF_000123456.1
        return "_".join(gtdb_accession.split('_')[1:3])
    
    # Function to download a single genome
    def download_genome(row):
        accession = extract_accession(row['accession'])
        output_file = f"{args.output_dir}/fna/{accession}.fna"
        
        if os.path.exists(output_file):
            return  # Skip if already downloaded
        
        # Create NCBI download URL
        # Format: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/123/456/GCA_000123456.1/GCA_000123456.1_genomic.fna.gz
        gca_gcf = accession[:3]  # GCA or GCF
        digits = accession.split('_')[1].split('.')[0]  # 000123456
        
        # Split into directory parts
        dir_parts = [digits[i:i+3] for i in range(0, len(digits), 3)]
        
        url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{gca_gcf}/{'/'.join(dir_parts)}/{accession}/{accession}_genomic.fna.gz"
        
        tmp_file = f"{output_file}.gz"
        try:
            subprocess.run(["wget", "-q", "-O", tmp_file, url])
            
            # Decompress
            with gzip.open(tmp_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            # Remove temporary file
            os.remove(tmp_file)
        except Exception as e:
            print(f"Error downloading {accession}: {e}")
    
    # Download genomes in parallel
    print(f"Downloading bacterial genomes using {args.threads} threads...")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        list(tqdm(executor.map(download_genome, [row for _, row in metadata.iterrows()]), 
                  total=len(metadata)))
    
    print("Download complete!")

if __name__ == "__main__":
    main()
