#!/bin/bash

set -e
apt list --installed | grep -w -E 'samtools|bcftools|subread'
BASE_DIR=$(pwd)
echo working in $BASE_DIR
mkdir -p results/bam results/bcf results/vcf results/trimmed results/orphaned

# Trim the sequence reads
for file in data/*_1.fastq.gz; do
    SRR=$(basename $file _1.fastq.gz)
    echo working on $SRR
    
    # Check if trimmed files already exist
    if [ -e results/trimmed/${SRR}_1.trim.fastq.gz ] && [ -e results/trimmed/${SRR}_2.trim.fastq.gz ]; then
        echo "Trimmed files already exist for $SRR. Skipping trimming step."
        continue
    fi
    
    TrimmomaticPE $BASE_DIR/data/${SRR}_1.fastq.gz $BASE_DIR/data/${SRR}_2.fastq.gz \
                  $BASE_DIR/results/trimmed/${SRR}_1.trim.fastq.gz $BASE_DIR/results/orphaned/${SRR}_1.untrim.fastq.gz \
                  $BASE_DIR/results/trimmed/${SRR}_2.trim.fastq.gz $BASE_DIR/results/orphaned/${SRR}_2.untrim.fastq.gz \
                  SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 \
                  -threads 16
done

# Build index using Subread
genome=$BASE_DIR/data/genome/ecoli_rel606.fasta
subread-buildindex -o $genome.index $genome

# Align the sequence directly to BAM using Subread
for fq1 in results/trimmed/*_1.trim.fastq.gz; do
     SRR=$(basename $fq1 _1.trim.fastq.gz)
     echo "working with $SRR"
     fq1=results/trimmed/${SRR}_1.trim.fastq.gz
     fq2=results/trimmed/${SRR}_2.trim.fastq.gz
     sam=results/sam/${SRR}.sam
     bam=results/bam/${SRR}.bam
     raw_bcf=results/bcf/${SRR}_raw.bcf.gz
     variants=results/vcf/${SRR}_variants.vcf
     final_variants=results/vcf/${SRR}_final_variants.vcf
     
     # Use Subread for alignment
     subread-align -t 1 -i $genome.index -r $fq1 -R $fq2 -o $bam -M 4000 -T 8
     #samtools view -S -b -o $bam $sam
     samtools sort -o ${bam%.bam}-sorted.bam $bam  # Sort the BAM file
     samtools index ${bam%.bam}-sorted.bam  # Index the sorted BAM file
     bcftools mpileup -Oz -o $raw_bcf -f $genome ${bam%.bam}-sorted.bam --threads 16
     bcftools call --ploidy 1 -m -v -o $variants $raw_bcf --threads 16
     vcfutils.pl varFilter $variants > $final_variants
     ls -hl $final_variants
done
