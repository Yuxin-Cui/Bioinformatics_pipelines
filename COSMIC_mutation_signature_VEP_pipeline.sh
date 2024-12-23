#!/bin/bash

# Activate the conda environment
conda activate variant-calling

# A simplified pipeline to predict COSMIC mutation signature using VEP

# Check VEP installation
which vep

# Set working directory dynamically
WORKDIR=$(pwd)

# Define input/output paths using the current working directory
INPUT="$WORKDIR/input/input.vcf.gz"
CACHE="$WORKDIR/cache"
OUTPUT_PATH="$WORKDIR/output/cosmic"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_PATH

# Run VEP to generate a VCF file (no --filter_common for COSMIC analysis)
OUTPUTVCF="$OUTPUT_PATH/vep_out.vcf"
echo "Running VEP to generate VCF output..."
vep -i $INPUT \
    -o $OUTPUTVCF \
    --cache \
    --dir_cache $CACHE \
    --everything \
    --fork 8 \
    --format vcf \
    --per_gene \
    --total_length \
    --force_overwrite

# Run VEP to generate a TXT file (no --filter_common for COSMIC analysis)
OUTPUT="$OUTPUT_PATH/vep_out.txt"
echo "Running VEP to generate TXT output..."
vep -i $INPUT \
    -o $OUTPUT \
    --cache \
    --dir_cache $CACHE \
    --everything \
    --fork 8 \
    --format vcf \
    --per_gene \
    --stats_file "$OUTPUT_PATH/vep_summary.html" \
    --total_length \
    --force_overwrite

# Filter VEP output for COSMIC annotations
echo "Filtering VEP output for COSMIC mutations..."
filter_vep -i $OUTPUT -o "$OUTPUT_PATH/vep_out_filtered.txt" -filter "SOMATIC" --force_overwrite

# Run R script to clean and visualize the results
echo "Processing results with R..."
Rscript -e "
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr)

# Load VEP filtered output
vep_data <- read.table('$OUTPUT_PATH/vep_out_filtered.txt', sep = '\t', header = F, stringsAsFactors = FALSE)

# Define column names
colnames(vep_data) <- c('Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 
                        'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 
                        'Amino_acids', 'Codons', 'Existing_variation', 'Extra')

# Parse the 'Extra' column into key-value pairs
extra_key_value_lists <- lapply(vep_data\$Extra, function(x) {
  key_value_pairs <- strsplit(x, ';')[[1]]
  key_value_pairs <- strsplit(key_value_pairs, '=')
  key_value_pairs
})

# Extract unique keys
all_keys <- unique(unlist(lapply(extra_key_value_lists, function(x) {
  sapply(x, function(pair) pair[1])
})))

# Create a data frame with all keys as columns
extra_df <- do.call(rbind, lapply(extra_key_value_lists, function(x) {
  key_value_list <- setNames(as.list(rep(NA, length(all_keys))), all_keys)
  for (pair in x) {
    key_value_list[[pair[1]]] <- pair[2]
  }
  as.data.frame(key_value_list, stringsAsFactors = FALSE)
}))

# Merge parsed data back into original
vep_data <- cbind(vep_data, extra_df)

# Clean up columns
vep_data <- vep_data %>%
  mutate(across(starts_with('IMPACT') | starts_with('SYMBOL') |
                  starts_with('AF') | starts_with('PHENO') | starts_with('CLIN_SIG'), 
                ~ gsub('=', '', .)))

# Create a plot for consequence frequency
plot <- ggplot(vep_data, aes(x = Consequence)) +
  geom_bar(fill = 'steelblue', color = 'black') +
  labs(title = 'Frequency of Consequence Categories',
       x = 'Consequence',
       y = 'Count') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot and cleaned data
ggsave('$OUTPUT_PATH/consequence_frequency_plot.png', plot = plot, width = 8, height = 6)
write.csv(vep_data, '$OUTPUT_PATH/vep_out_filtered.csv', row.names = FALSE)
"

# Verify pipeline completion
if [ -f "$OUTPUT_PATH/vep_out_filtered.csv" ]; then
    echo "Pipeline completed successfully!"
else
    echo "Error: COSMIC signature output not found. Pipeline failed."
fi
