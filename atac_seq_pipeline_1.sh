# Create working environment
# conda install -c conda-forge openjdk=20.0.0
# conda install -c bioconda -y fastqc
# conda install -c bioconda -y trimmomatic
# conda install -c bioconda -y subread
# conda install -c bioconda -y samtools
# conda install -c bioconda -y picard
# conda install -c bioconda -y macs2
# conda install -c bioconda -y bedtools
# conda install -c bioconda -y deeptools
# conda install -c bioconda -y homer

# mkdir
mkdir -p ./data/ ./genome/ data/trimming results/fastqc results/sam results/bam 

# Download the data
cd data
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_1.fastq.gz

# Create a sample id list
ls SRR*gz | cut -f 1 -d "_" | sort | uniq  > sample_ids.txt
cd ..

# Run fastqc or falco
## fastqc -t 8 *fastq.gz -o data/fastqc

for sample_id in $(cat data/sample_ids.txt); do
    mkdir -p results/fastqc/"$sample_id"_1
    mkdir -p results/fastqc/"$sample_id"_2
    falco "data/${sample_id}_1.fastq.gz" -o "results/fastqc/${sample_id}_1" -t 16
    falco "data/${sample_id}_2.fastq.gz" -o "results/fastqc/${sample_id}_2" -t 16
done

# Run multiqc 
cd results/fastqc
multiqc .
firefox multiqc_report.html

# Trimmingcd ..
cd ../../data/trimming

for sample_id in $(cat ../sample_ids.txt); do
    echo "Processing sample: $sample_id"
    trimmomatic PE -threads 2 \
    ../${sample_id}_1.fastq.gz ../${sample_id}_2.fastq.gz \
    ${sample_id}_tr_1P.fastq.gz ${sample_id}_tr_1U.fastq.gz \
    ${sample_id}_tr_2P.fastq.gz ${sample_id}_tr_2U.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:/usr/share/trimmomatic/NexteraPE-PE.fa:2:40:15
done

# Download reference genome sequence file and build index in R/Rstudio:
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
mainChromosomes <- paste0("chr",c(1:21,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Hsapiens.UCSC.hg38[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
writeXStringSet(mainChrSeqSet,
                "BSgenome.Hsapiens.UCSC.hg38.mainChrs.fa")

library(Rsubread)
buildindex("BSgenome.Hsapiens.UCSC.hg38.mainChrs",
           "BSgenome.Hsapiens.UCSC.hg38.mainChrs.fa",
           indexSplit = TRUE,
           memory = 4000)

# Run alignment in Ubuntu:
export SUBREAD_INDEX=$(pwd)/BSgenome.Hsapiens.UCSC.hg38.mainChrs
for sample_id in $(cat data/sample_ids.txt); do
    echo "running $sample_id"
    subread-align -i "$SUBREAD_INDEX" \
        -r data/"${sample_id}"_1.fastq.gz \
        -R data/"${sample_id}"_2.fastq.gz \
        -t 1 -o results/bam/"${sample_id}".bam \
        -T 8
done

# Sort and index
mkdir -p results/sorted_bam results/mappedrReads
for sample_id in $(cat data/sample_ids.txt); do
    echo "running $sample_id"
    samtools sort -@ 8 results/bam/"${sample_id}".bam -o - | \
    samtools index - results/sorted_bam/"${sample_id}".sorted.bam && \
    samtools idxstats results/sorted_bam/"${sample_id}".sorted.bam > results/mappedReads/"${sample_id}".mappedReads.txt
    echo "Done processing $sample_id"
done

# Remove duplicate
cd results/sorted_bam

for sample_id in $(cat ../../data/sample_ids.txt); do
    echo "running $sample_id"
    samtools rmdup "${sample_id}.sorted.bam" "${sample_id}.rmdup.bam"
    samtools index "${sample_id}.rmdup.bam"
    samtools flagstat "${sample_id}.rmdup.bam" > "./${sample_id}.rmdup.stat"
done

# Secondary clean up and QC
for sample_id in $(cat ../../data/sample_ids.txt); do
  samtools view -h -f 2 -q 30 ${sample_id}.rmdup.bam \
  | grep -v -e "mitochondria" -e "\\*" \
  | samtools sort -O bam -@ 10 -o ${sample_id}.last.sorted.bam - \
  && samtools index ${sample_id}.last.sorted.bam \
  && samtools flagstat ${sample_id}.last.sorted.bam > ./${sample_id}.last.stat
done

# Calculate insert size
for sample_id in $(cat ../../data/sample_ids.txt); do

    # Calculate insert size statistics using samtools and awk
    samtools view -@ 4 ./${sample_id}.last.sorted.bam | \
        awk '$2==99 || $2==147' | \
        awk -F '\t' '{print $9}' | \
        awk '{if($1 > 0) print $1;}' > ./${sample_id}.insert_sizes.txt

    # Calculate mean insert size and insert size standard deviation
    mean_insert_size=$(awk '{ total += $1; count++ } END { print total/count }' ./${sample_id}.insert_sizes.txt)
    insert_size_sd=$(awk -v mean="$mean_insert_size" '{ sum += ($1 - mean)^2 } END { print sqrt(sum / (NR-1)) }' ./${sample_id}.insert_sizes.txt)

    # Print insert size metrics
    echo "Filename: ${sample_id}"
    echo "Mean Insert Size: ${mean_insert_size}"
    echo "Insert Size Standard Deviation: ${insert_size_sd}"
done

# Save the codes below as "plot_histogram.py", and then run "python3 plot_histogram.py"

import matplotlib.pyplot as plt

filename = "SRR891269"  # Replace with the appropriate filename

insert_sizes = []
with open(f"./{filename}.insert_sizes.txt", "r") as f:
    for line in f:
        insert_sizes.append(int(line.strip()))

plt.hist(insert_sizes, bins=50, color='blue', edgecolor='black')
plt.title("Insert Size Distribution")
plt.xlabel("Insert Size")
plt.ylabel("Frequency")
plt.savefig(f"./{filename}.insert_size_histogram.pdf")
plt.show()

# peak calling
for sample_id in $(cat ../../data/sample_ids.txt); do
   effective_genome_size=119481543

   # Convert BAM to BED using bedtools
   bedtools bamtobed -i ${sample_id}.last.sorted.bam > ${sample_id}.last.bed

   # Call peaks using macs2
   macs2 callpeak \
      -t ${sample_id}.last.bed \
      -g ${effective_genome_size} \
      --nomodel --shift -100 --extsize 200 \
      -n ${sample_id}.last \
      -q 0.01 --outdir ./peaks
done

# Create .bw file for visualization using igv
for sample_id in $(cat ../../data/sample_ids.txt); do
    bamCoverage --bam ${sample_id}.last.sorted.bam -o ${sample_id}.last.bw 
        --binSize 10 
        --normalizeUsing RPGC 
        --effectiveGenomeSize ${effective_genome_size} 
        --ignoreForNormalization chrX 
        --extendReads
done

# igv

