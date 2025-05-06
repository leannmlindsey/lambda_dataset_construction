#!/usr/bin/env python3
"""
Download all phage FNA files and metadata from PhageScope.
"""
import os
import requests
import pandas as pd
import subprocess
import gzip
import shutil
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Download phage data from PhageScope')
    parser.add_argument('--output_dir', type=str, default='phagescope_data',
                        help='Directory to store downloaded files')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of parallel downloads')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directories
    os.makedirs(f"{args.output_dir}/fna", exist_ok=True)
    os.makedirs(f"{args.output_dir}/metadata", exist_ok=True)
    
    # PhageScope API for phage genomes and metadata
    # Note: You might need to adapt this based on PhageScope's actual API
    base_url = "http://phagescope.org/phagedb/api"
    
    print("Fetching PhageScope genome list...")
    # Get list of available phages
    response = requests.get(f"{base_url}/phages")
    phage_list = response.json()
    
    print(f"Found {len(phage_list)} phage genomes")
    
    # Download metadata
    metadata_url = "http://phagescope.org/phagedb/download/metadata"
    metadata_file = f"{args.output_dir}/metadata/phage_metadata.csv"
    
    print(f"Downloading complete metadata file to {metadata_file}")
    subprocess.run(["wget", "-O", metadata_file, metadata_url])
    
    # Function to download and process a single phage genome
    def download_phage(phage_id):
        fna_url = f"{base_url}/phage/{phage_id}/fna"
        output_file = f"{args.output_dir}/fna/{phage_id}.fna"
        
        if os.path.exists(output_file):
            return  # Skip if already downloaded
            
        # Download gzipped fna file
        tmp_file = f"{output_file}.gz"
        try:
            subprocess.run(["wget", "-q", "-O", tmp_file, fna_url])
            
            # Decompress
            with gzip.open(tmp_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            # Remove temporary file
            os.remove(tmp_file)
        except Exception as e:
            print(f"Error downloading {phage_id}: {e}")
    
    # Download phage genomes in parallel
    print(f"Downloading phage genomes using {args.threads} threads...")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        list(tqdm(executor.map(download_phage, [p['id'] for p in phage_list]), 
                  total=len(phage_list)))
    
    print("Download complete!")

if __name__ == "__main__":
    main()
