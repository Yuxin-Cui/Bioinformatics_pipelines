#!/bin/bash

set -e  

#--------------------------------------------------------
# Bioinformatics Pipeline for the process of ONT data
# This script automates the steps for base calling, quality control, and data cleaning
# Requirements: Chiron, pycoQC, R, Porechop, pysam, tensorflow (1.15.0) installed.
#--------------------------------------------------------

conda activate chiron_env

# Set default paths (can be overridden by user)
INPUT_DIR="./ont_data/practicals"
OUTPUT_DIR="./output"
CHIRON_MODEL_DNA="chiron/chiron/model/DNA_default"
CHIRON_MODEL_RNA="chiron/chiron/model/RNA_default"

# Log function for better status tracking
log() {
    echo "[INFO] $(date): $1"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# 1. Check if required tools are available
log "Checking for required tools..."

# Check if chiron, pycoQC, R, and porechop are installed
required_commands=("chiron" "pycoQC" "R" "porechop")
for cmd in "${required_commands[@]}"; do
    if ! command_exists "$cmd"; then
        echo "[ERROR] $cmd is not installed. Exiting."
        exit 1
    fi
done

# 2. Install Chiron (if not already installed)
log "Installing Chiron..."
pip show chiron || pip install chiron || { echo "[ERROR] Failed to install Chiron. Exiting."; exit 1; }

# 3. Perform basecalling with Chiron for DNA (default)
log "Running Chiron basecalling for DNA..."
chiron call -i "$INPUT_DIR/basecalling_practical/fast5/" \
            -o "$OUTPUT_DIR/fastq2" \
            -m "$CHIRON_MODEL_DNA" \
            --preset dna-pre

# Uncomment and use the following line if RNA basecalling is needed
# log "Running Chiron basecalling for RNA..."
# chiron call -i "./chiron/chiron/example_data/RNA" \
#             -o "$OUTPUT_DIR/fastq" \
#             -m "$CHIRON_MODEL_RNA" \
#             --preset rna-pre --mode RNA

# 4. Run PycoQC for quality control
log "Running PycoQC on sequencing summary..."
pycoQC -f "$INPUT_DIR/qc_practical/summaries/run_1/sequencing_summary.txt" \
       -o "$OUTPUT_DIR/run_1.html"

# 5. Find all qc_practical summaries and display them
log "Searching for qc_practical summaries..."
find "$INPUT_DIR" -path "*/qc_practical/summaries"

# 6. Run MinIONQC (R script for further QC)
log "Running MinIONQC R script..."
Rscript MinIONQC.R -i "$INPUT_DIR/qc_practical/summaries" \
                  -o "$OUTPUT_DIR/minion_qc"

# 7. Run Porechop for adapter trimming (example only, change path as needed)
log "Running Porechop to remove adapters..."
porechop -i "$OUTPUT_DIR/all_guppy.fastq" \
         -o "$OUTPUT_DIR/porechopped.fastq" \
         --discard_middle

log "Pipeline completed successfully!"
