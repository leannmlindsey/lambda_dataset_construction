#!/bin/bash
# BLAST phage sequences against bacteria database

set -e

# Default parameters
PHAGE_FASTA="blast_db/combined_phage.fasta"
BACTERIA_DB="blast_db/bacteria_db"
OUTPUT_DIR="blast_results"
THREADS=8
IDENTITY=90
COVERAGE=80

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --phage_fasta)
        PHAGE_FASTA="$2"
        shift
        shift
        ;;
        --bacteria_db)
        BACTERIA_DB="$2"
        shift
        shift
        ;;
        --output_dir)
        OUTPUT_DIR="$2"
        shift
        shift
        ;;
        --threads)
        THREADS="$2"
        shift
        shift
        ;;
        --identity)
        IDENTITY="$2"
        shift
        shift
        ;;
        --coverage)
        COVERAGE="$2"
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

# Run BLAST
echo "Running BLAST: phage sequences against bacteria database..."
OUTPUT_FILE="$OUTPUT_DIR/phage_vs_bacteria.tsv"

blastn -query "$PHAGE_FASTA" \
       -db "$BACTERIA_DB" \
       -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
       -perc_identity "$IDENTITY" \
       -qcov_hsp_perc "$COVERAGE" \
       -num_threads "$THREADS" \
       -out "$OUTPUT_FILE"

echo "BLAST completed. Results saved to: $OUTPUT_FILE"
echo "Found $(wc -l < "$OUTPUT_FILE") hits"
