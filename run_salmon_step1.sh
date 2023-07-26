#!/bin/bash

salmon --version
mkdir -p ref_genome/ref_index

# URL for the reference file
url="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz"

# Check if the reference file already exists
if file_exists "ref_genome/GRCh38_latest_rna.fna.gz"; then
  echo "Reference file already exists. Skipping download."
else
  # Download the reference file
  echo "Downloading reference file..."
  curl -o ref_genome/GRCh38_latest_rna.fna.gz $url
fi


salmon index -t $(pwd)/ref_genome/GRCh38_latest_rna.fna.gz -i $(pwd)/ref_genome/ref_index -p 8


#index=/mnt/e/salmon_test/lv_sig_index

#for fn in /mnt/e/rnaseq_data/clean/Sample_{1..15}A;
#do
#samp=`basename ${fn}`
#echo "Processing sample ${samp}"
#/home/yc/salmon-latest_linux_x86_64/bin/salmon quant -i $index -l A \
#         -1 ${fn}/${samp}_1.fq.gz \
#         -2 ${fn}/${samp}_2.fq.gz \
#         -p 10 --validateMappings -o quants/${samp}_quant
#done 
