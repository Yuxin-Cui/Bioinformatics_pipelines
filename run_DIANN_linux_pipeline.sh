#!/bin/bash

# Install mono
# Guide: https://www.mono-project.com/download/stable/#download-lin
# sudo apt install ca-certificates gnupg
# sudo gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver # hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
# echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] https://download.mono-project.com/repo/ubuntu stable-focal main" | sudo tee /etc/apt/sources.list.d/mono-official-stable.list
# sudo apt update
# sudo apt install mono-devel

# Convert raw to mzWL if necessary
# mono /mnt/e/bioinfo_hub/immunopeptidomics/ThermoRawFileParser1.4.5/ThermoRawFileParser.exe -d=/mnt/e/bioinfo_hub/immunopeptidomics/data/

# Download uniprot database and universal contamination fasta (2022) directly from raw URL
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
# wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta # optional
# source to download: https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/blob/main/Universal%20protein%20contaminant%20FASTA/0602_Universal%20Contaminants.fasta
# # mv 0602_Universal\ Contaminants.fasta 0602_Universal_Contaminants.fasta
# cat UP000005640_9606.fasta 0602_Universal_Contaminants.fasta >  human_protein_and_universal_contaminants.fasta

# Install diann-linux
# wget https://github.com/vdemichev/DiaNN/releases/download/1.9.2/diann-1.9.2.Linux_update_2024-10-31.zip
# cd diann-1.9.2/
# pwd
# echo 'export PATH=$PATH:/mnt/e/bioinfo_hub/immunopeptidomics/diann-1.9.2' >> ~/.bashrc && source ~/.bashrc

# Test it
diann-linux

##################################################################################
# Library-free search mode (DIA)

BASE_DIR="$(pwd)"

# Parameters (same as Docker template) - modified paths
INPUT_FILES=(
    "$BASE_DIR/data/20100611_Velos1_TaGe_SA_Hela_2.mzML"
)

OUTPUT_FILE="$BASE_DIR/result/report.tsv"
TEMP_DIR="$BASE_DIR/quant/"
FASTA_FILE="$BASE_DIR/human_protein_and_universal_contaminants.fasta"

mkdir -p "$BASE_DIR/result/" "$TEMP_DIR"

# Run diann-linux with the same parameters as the Docker template
diann-linux \
    --f "${INPUT_FILES[@]}" \
    --lib "" \
    --threads 12 \
    --verbose 1 \
    --out "$OUTPUT_FILE" \
    --predictor \
    --qvalue 0.01 \
    --matrices \
    --temp "$TEMP_DIR" \
    --fasta "$FASTA_FILE" \
    --fasta-search \
    --met-excision \
    --cut "K*,R*" \
    --missed-cleavages 1 \
    --min-pep-len 5 \
    --max-pep-len 30 \
    --min-pr-mz 400 \
    --max-pr-mz 1200 \
    --min-pr-charge 1 \
    --max-pr-charge 4 \
    --unimod4 \
    --var-mods 1 \
    --var-mod "UniMod:35,15.994915,M" \
    --reanalyse \
    --smart-profiling
