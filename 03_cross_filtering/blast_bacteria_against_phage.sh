#!/bin/bash
# BLAST bacteria sequences against phage database

set -e

# Default parameters
BACTERIA_FASTA="blast_db/combined_bacteria.fasta"
PHAGE_DB="blast_db/phage_db"
OUTPUT_DIR="blast_results"
THREADS=8
IDENTITY=90
COVERAGE=80

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --bacteria_fasta)
        BACTERIA_FASTA="$2"
        shift
        shift
        ;;
        --phage_db)
        PHAGE_DB="$2"
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
echo "Running BLAST: bacteria sequences against phage database..."
OUTPUT_FILE="$OUTPUT_DIR/bacteria_vs_phage.tsv"

blastn -query "$BACTERIA_FASTA" \
       -db "$PHAGE_DB" \
       -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
       -perc_identity "$IDENTITY" \
       -qcov_hsp_perc "$COVERAGE" \
       -num_threads "$THREADS" \
       -out "$OUTPUT_FILE"

echo "BLAST completed. Results saved to: $OUTPUT_FILE"
echo "Found $(wc -l < "$OUTPUT_FILE") hits"
