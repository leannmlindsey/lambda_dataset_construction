#!/usr/bin/env python3
"""
Combine all downloaded phage metadata into a single file.
"""
import os
import pandas as pd
import glob
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Combine phage metadata files')
    parser.add_argument('--input_dir', type=str, default='phagescope_data/metadata',
                        help='Directory containing metadata files')
    parser.add_argument('--output_file', type=str, default='phagescope_data/combined_metadata.csv',
                        help='Output combined metadata file')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Get all metadata files
    metadata_files = glob.glob(f"{args.input_dir}/*.csv")
    
    print(f"Found {len(metadata_files)} metadata files")
    
    # Combine all metadata
    all_metadata = []
    for file in metadata_files:
        try:
            metadata = pd.read_csv(file)
            all_metadata.append(metadata)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    # Concatenate all metadata
    combined_metadata = pd.concat(all_metadata, ignore_index=True)
    
    # Remove duplicates based on phage ID
    combined_metadata = combined_metadata.drop_duplicates(subset=['phage_id'])
    
    print(f"Combined metadata contains {len(combined_metadata)} unique phage entries")
    
    # Save combined metadata
    combined_metadata.to_csv(args.output_file, index=False)
    print(f"Combined metadata saved to {args.output_file}")

if __name__ == "__main__":
    main()
