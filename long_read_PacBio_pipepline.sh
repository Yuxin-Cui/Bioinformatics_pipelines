#!/bin/bash

# Step 1: Create necessary directories
mkdir -p input output human_ref

# Step 2: Set a memory limit (e.g., 8GB) using ulimit
echo "Setting memory limit..."
ulimit -v 6388608  # 8GB in KB

# Step 3: Check if the 'isoseq' conda environment exists, and create it if it doesn't
if conda env list | grep -q "isoseq"; then
  echo "The 'isoseq' environment already exists. Activating it..."
else
  echo "The 'isoseq' environment does not exist. Creating and activating it..."
  conda create --name isoseq -y
fi

# Initialize Conda in case it's not set up yet
if ! conda info --envs &> /dev/null; then
  echo "Initializing Conda for the shell..."
  conda init
  source ~/.bashrc  # or source ~/.zshrc for zsh
fi

# Activate the isoseq environment
conda activate isoseq

# Step 4: Check if the required tools (isoseq, lima, pbmm2) are available
check_tool_installed() {
  tool_name=$1
  if ! conda list | grep -q "$tool_name"; then
    echo "$tool_name is not installed. Installing..."
    conda install -c bioconda "$tool_name" -y
  else
    echo "$tool_name is already installed."
  fi
}

# Check if the required tools are installed
check_tool_installed isoseq
check_tool_installed lima
check_tool_installed pbmm2

# Step 5: Check if the hg38.fa file already exists in the input directory
if [ ! -f human_ref/hg38.fa ]; then
  echo "Reference genome (hg38.fa) not found in input/. Downloading it..."
  wget -P human_ref/ https://downloads.pacbcloud.com/public/dataset/ISMB_workshop/hg38.fa
else
  echo "Reference genome (hg38.fa) already exists. Skipping download."
fi

# Step 6: Download public example data if not already present in the input directory
if [ ! -d "input/isoseq3" ]; then
  echo "Downloading public example data..."
  wget -r -np -nH --cut-dirs=3 -R "index.html*" -P input/ https://downloads.pacbcloud.com/public/dataset/ISMB_workshop/isoseq3/
else
  echo "Public example data already exists in input/isoseq3/. Skipping download."
fi

# Step 7: Create human reference genome index (if not already indexed)
#echo "Creating human reference genome index..."
#pbmm2 index human_ref/hg38.fa human_ref/hg38.mmi -j 4

# Step 8: Process the .bam files in the input folder and its subdirectories
echo "Searching for .bam files in the 'input' directory and subdirectories..."
find input/ -type f -name "*.ccs.bam" | while read bam_file; do
  if [ -f "$bam_file" ]; then
    echo "Found BAM file: $bam_file"
    
    # Step 9: Demultiplexing
    echo "Running demultiplexing on $bam_file..."
    lima --isoseq --dump-clips --peek-guess -j 8 "$bam_file" input/isoseq3/isoseq_primers.fasta output/alz.demult.bam
    
    # Step 10: Refining the reads
    echo "Refining the reads..."
    isoseq3 refine --require-polya output/alz.demult.5p--3p.bam input/isoseq3/isoseq_primers.fasta output/alz.flnc.bam
    
    # Step 11: View the output files
    echo "Listing output files..."
    ls -ltrh output/

    # Step 12: Clustering the reads
    echo "Running clustering..."
    isoseq3 cluster output/alz.flnc.bam output/alz.polished.bam --verbose --use-qvs

    # Step 13: Aligning using pbmm2
    # Minimal RAM 64Gb required
    echo "Aligning reads with pbmm2..."
    pbmm2 align human_ref/hg38.fa output/alz.polished.hq.bam output/alz.aligned.bam -j 4 --preset ISOSEQ –-sort --log-level INFO
  
  else
    echo "No BAM files found in the 'input' directory or subdirectories."
  fi
done

echo "Pipeline finished!"
