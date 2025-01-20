#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define the parameters for the input and output directories
params.input_dir = "$baseDir/data"
params.temp_dir = "$baseDir/quant"
params.fasta_file = "$baseDir/human_protein_and_universal_contaminants.fasta"
params.input_file = "" // Add a parameter for the input file
params.publish_dir = "$baseDir/result/nf_output"

process diann_analysis {

    // Define the input files for the process
    input:
    path input_files
    path fasta_file

    // Publish the output to a specified directory (e.g., published_results)
    publishDir(params.publish_dir, mode: 'copy') // Can use 'link' or 'move' as well

    // Define the container (optional: use Docker or Singularity)
    container 'diann-linux:latest'

    // Define the script that runs DiaNN
    script:
    """
    mkdir -p ${params.publish_dir} ${params.temp_dir}  // Ensure output and temp directories exist

    diann-linux \
        --f "${input_files}" \
        --lib "" \
        --threads 12 \
        --verbose 2 \
        --out "report.tsv" \
        --predictor \
        --qvalue 0.01 \
        --matrices \
        --temp "${params.temp_dir}" \
        --fasta "${fasta_file}" \
        --fasta-search \
        --met-excision \
        --cut "K*,R*" \
        --missed-cleavages 1 \
        --min-pep-len 5 \
        --max-pep-len 30 \
        --min-pr-mz 400 \
        --max-pr-mz 1200 \
        --min-pr-charge 1 \
        --max-pr-charge 4 \
        --unimod4 \
        --var-mods 1 \
        --var-mod "UniMod:35,15.994915,M" \
        --reanalyse \
        --smart-profiling
    """
}

// Define the channels (after process definitions)
Channel
    .fromPath(params.input_file ?: params.input_dir + "/20100611_Velos1_TaGe_SA_Hela_2.mzML") // Use the input file if provided
    .set { input_files }

workflow {
    // Call the diann_analysis process and explicitly pass the input channels
    diann_analysis(input_files, params.fasta_file)
}