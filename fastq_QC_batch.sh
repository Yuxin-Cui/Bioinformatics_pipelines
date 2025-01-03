# Make directories if needed
mkdir -p 00ref 01data 02QC 03Mapping 04Quantification

# Make a hyperlink
OUT=02QC

# Make a sample list
find $(pwd) -type f -name "*.fq.gz" > $(pwd)/$OUT/samplelist.txt

# Set up the right environment. 
# It appears python 3.12 is not compatible for the running of falco
conda activate py311

# Run falco
# conda install -c bioconda falco
for i in $(cat $OUT/samplelist.txt)
do
   base_name=$(basename "$i" .fq.gz)
   falco $i -o "$OUT/$base_name" -t 12
done

# Generate a merged QC summary using MultiQC
# MultiQC will automatically detect the falco output files
multiqc $OUT -o $OUT/multiqc_report

# Optional QC tool (slower)
#for i in $(cat $OUT/samplelist.txt)
#do
#   fastp -i "$i" -h "$OUT/${base_name}_fastp.html" -j "$OUT/${base_name}_fastp.json" -t 12 --disable_length_filtering
#done