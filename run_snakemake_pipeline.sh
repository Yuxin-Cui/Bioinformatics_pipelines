# Make sure the snakefile is validated using snakemake -n
# conda activate snakemake
snakemake --snakefile snakefile_bwa_2  -c 16
eog plots/quals.svg
snakemake --snakefile snakefile_bwa_2 --rulegraph | dot -Tsvg > rulegraph.svg
eog rulegraph.svg
 
