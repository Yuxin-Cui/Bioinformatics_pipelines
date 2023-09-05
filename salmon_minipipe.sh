#!/bin/bash

# Check if the Salmon index exists, and create it if it doesn't
if [ ! -d "/mnt/e/genome/index_salmon" ]; then
    salmon index --threads 16 -t /mnt/e/genome/GRCh38_latest_rna.fna.gz -i /mnt/e/genome/index_salmon
fi

# Generate the sample list
ls data/*_1.fq.gz | sed 's:data/::; s/_1.fq.gz//' > sample_list.txt

# Perform Salmon quantification for each sample in the sample list
while IFS= read -r sample; do
    output_dir="results/${sample}"

    # Check if the output directory already exists
    if [ -d "$output_dir" ]; then
        echo "Skipping ${sample} as the output directory already exists."
    else
        salmon quant --threads 16 -i /mnt/e/genome/index_salmon -l A \
            -1 "data/${sample}_1.fq.gz" -2 "data/${sample}_2.fq.gz" -o "$output_dir"
    fi
done < sample_list.txt
