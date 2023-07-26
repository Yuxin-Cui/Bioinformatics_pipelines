#!/bin/bash
salmon --version

# Check if the file ref_index/pos.bin exists
if [ ! -e "ref_genome/ref_index/pos.bin" ]; then
  echo "Index files not found. Running run_salmon_step1.sh..."
  # Execute run_salmon_step1.sh
  bash run_salmon_step1.sh
else
  echo "Index files found. Proceeding with sequence alignment..."
fi


# Alignment
index=ref_genome/ref_index
for fn in data/Sample_{1..15}A;
do
  samp=`basename ${fn}`
  echo "Processing sample ${samp}"
  salmon quant -i $index -l A \
               -1 ${fn}/${samp}_1.fq.gz \
               -2 ${fn}/${samp}_2.fq.gz \
               -p 10 --validateMappings -o quants/${samp}_quant
done 
