#!/bin/bash

# Set the input
INPUT=$(pwd)/input/input.vcf.gz

# Check whether the tools have been installed
conda activate variant-calling
which gatk
java --version
which snpEff

# Set path if needed

# Create output directories if they don't exist
mkdir -p "$(pwd)/output/snpeff"

# Run snpEff
java -Xmx6g -jar snpEff/snpEff.jar hg38 "$INPUT" > "$(pwd)/output/snpeff/variants.snpeff.vcf"

gatk VariantsToTable -V "$(pwd)/output/snpeff/variants.snpeff.vcf" -F CHROM -F POS -F TYPE -F ID -F ANN -F LOF -F NMD -GF AD -GF DP -GF GQ -GF GT -O "$(pwd)/output/snpeff/variants.snpeff.txt"

# Download clinvar database (vcf.gz and tbi)
# mkdir -p "$(pwd)/clinvar"
# cd "$(pwd)/clinvar"
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

# Run SnpSift
java -jar snpEff/SnpSift.jar annotate "$(pwd)/clinvar/clinvar.vcf.gz" "$INPUT" > "$(pwd)/output/snpeff/variants.clinvar.vcf"

gatk VariantsToTable -V "$(pwd)/output/snpeff/variants.clinvar.vcf" -F CHROM -F POS -F TYPE -F ID -F ALLELEID -F CLNDN -F CLNSIG -F CLNSIGCONF -F CLNSIGINCL -F CLNVC -F GENEINFO -GF AD -GF GQ -GF GT -O "$(pwd)/output/snpeff/variants.clinvar.txt"
