#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def trim_sequences(input_file, output_file, target_length=4000):
    """
    Trim phage sequences (label=1) to the target length.
    For phage sequences longer than target_length, trims equally from both ends.
    Leaves bacteria sequences (label=0) unchanged.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)
    
    # Count sequences before trimming
    total_before = len(df)
    phage_before = sum(df['label'] == 1)
    bacteria_before = sum(df['label'] == 0)
    
    print(f"Before trimming:")
    print(f"  Total sequences: {total_before}")
    print(f"  Bacteria (label=0): {bacteria_before}")
    print(f"  Phage (label=1): {phage_before}")
    
    # Process each row
    modified_count = 0
    for idx, row in df.iterrows():
        seq = row['sequence']
        label = row['label']
        
        # Only trim phage sequences (label=1) that are longer than target_length
        if label == 1 and len(seq) > target_length:
            seq_length = len(seq)
            
            # Calculate how much to trim from each end
            total_to_trim = seq_length - target_length
            trim_from_start = total_to_trim // 2
            trim_from_end = total_to_trim - trim_from_start
            
            # Trim the sequence
            trimmed_seq = seq[trim_from_start:seq_length-trim_from_end]
            
            # Update the sequence in the dataframe
            df.at[idx, 'sequence'] = trimmed_seq
            modified_count += 1
    
    # Verify sequence lengths after trimming
    phage_lengths = [len(seq) for idx, seq in enumerate(df['sequence']) if df['label'].iloc[idx] == 1]
    bacteria_lengths = [len(seq) for idx, seq in enumerate(df['sequence']) if df['label'].iloc[idx] == 0]
    
    # Write the modified dataframe to the output file
    df.to_csv(output_file, index=False)
    
    # Print a summary
    print(f"Modified {modified_count} phage sequences")
    print(f"After trimming:")
    print(f"  Phage sequences (label=1): min={min(phage_lengths)}, max={max(phage_lengths)}, mean={sum(phage_lengths)/len(phage_lengths):.2f}")
    print(f"  Bacteria sequences (label=0): min={min(bacteria_lengths)}, max={max(bacteria_lengths)}, mean={sum(bacteria_lengths)/len(bacteria_lengths):.2f}")
    print(f"Saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Trim phage sequences to match bacteria sequence length')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    parser.add_argument('-l', '--length', type=int, default=4000, 
                        help='Target sequence length (default: 4000)')
    
    args = parser.parse_args()
    
    trim_sequences(args.input, args.output, args.length)

if __name__ == "__main__":
    main()
