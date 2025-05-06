#!/bin/bash
# Create BLAST database for phage sequences

set -e

# Default parameters
PHAGE_DIR="phagescope_data/fna"
OUTPUT_DIR="blast_db"
DB_NAME="phage_db"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --phage_dir)
        PHAGE_DIR="$2"
        shift
        shift
        ;;
        --output_dir)
        OUTPUT_DIR="$2"
        shift
        shift
        ;;
        --db_name)
        DB_NAME="$2"
        shift
        shift
        ;;
        *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Combine all phage sequences into a single FASTA file
echo "Combining phage sequences..."
COMBINED_FASTA="$OUTPUT_DIR/combined_phage.fasta"

# Check if combined file exists and remove if it does
if [ -f "$COMBINED_FASTA" ]; then
    rm "$COMBINED_FASTA"
fi

# Concatenate all FASTA files
for fasta in "$PHAGE_DIR"/*.fna; do
    cat "$fasta" >> "$COMBINED_FASTA"
done

# Build BLAST database
echo "Creating BLAST database..."
makeblastdb -in "$COMBINED_FASTA" -dbtype nucl -out "$OUTPUT_DIR/$DB_NAME" -title "$DB_NAME" -parse_seqids

echo "BLAST database created: $OUTPUT_DIR/$DB_NAME"
