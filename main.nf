#!/usr/bin/env nextflow

/*
========================================================================================
    BEST-IN-CLASS VIROME ANALYSIS PIPELINE
========================================================================================
    Author: AI Agent for AS3243
    Version: 2.0
    Description: Production-ready virome analysis pipeline with comprehensive tools
    
    Tools included:
    - GeNomad: Virus identification and annotation
    - CheckV: Virus quality assessment
    - iPHoP: Host prediction for viruses
    - BACPHLIP: Bacterial/phage lifestyle prediction
    - vConTACT2: Viral clustering and taxonomy
========================================================================================
*/

nextflow.enable.dsl = 2

// Pipeline header
def helpMessage() {
    log.info"""
    =========================================
     VIROME ANALYSIS PIPELINE v2.0
    =========================================
    Usage:
    nextflow run main.nf [options]
    
    Required arguments:
      --input_dir         Path to directory containing FASTA files
      --outdir            Output directory path
    
    Optional arguments:
      --genomad_db        GeNomad database path
      --iphop_db          iPHoP database path  
      --bacphlip_hmmsearch BACPHLIP hmmsearch binary path
      --cluster_one_jar   Cluster ONE JAR file path
      --help              Show this help message
    
    Profiles:
      -profile slurm      Run on SLURM cluster
      -profile test       Run with test data
      -profile conda      Use conda environments
      
    Example:
      nextflow run main.nf -profile slurm --input_dir data --outdir results
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Parameter validation
def checkParams() {
    def errors = []
    
    // Check required parameters
    if (!params.input_dir) errors.add("Missing required parameter: --input_dir")
    if (!params.outdir) errors.add("Missing required parameter: --outdir")
    
    // Check database paths
    if (!file(params.genomad_db).exists()) errors.add("GeNomad database not found: ${params.genomad_db}")
    if (!file(params.iphop_db).exists()) errors.add("iPHoP database not found: ${params.iphop_db}")
    if (!file(params.bacphlip_hmmsearch).exists()) errors.add("BACPHLIP hmmsearch not found: ${params.bacphlip_hmmsearch}")
    if (!file(params.cluster_one_jar).exists()) errors.add("ClusterONE JAR not found: ${params.cluster_one_jar}")
    
    // Check input directory
    if (!file(params.input_dir).exists()) errors.add("Input directory not found: ${params.input_dir}")
    
    if (errors.size() > 0) {
        log.error "Parameter validation failed:"
        errors.each { log.error "  - ${it}" }
        exit 1
    }
    
    log.info "✓ Parameter validation passed"
}

// Validate parameters
checkParams()

// Print pipeline information
log.info """
=========================================
 VIROME ANALYSIS PIPELINE v2.0
=========================================
Input directory     : ${params.input_dir}
Output directory    : ${params.outdir}
GeNomad database    : ${params.genomad_db}
iPHoP database      : ${params.iphop_db}
BACPHLIP hmmsearch  : ${params.bacphlip_hmmsearch}
ClusterONE JAR      : ${params.cluster_one_jar}
Max memory          : ${params.max_memory}
Max CPUs            : ${params.max_cpus}
Max time            : ${params.max_time}
=========================================
"""

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process RUN_GENOMAD {
    tag "$sample_id"
    label 'process_high'
    conda 'genomad'
    
    publishDir "${params.outdir}/${sample_id}/genomad", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(contigs)
    
    output:
    tuple val(sample_id), path("${sample_id}_genomad/${contigs.baseName}_summary/${contigs.baseName}_virus.fna"), emit: virus_fasta
    tuple val(sample_id), path("${sample_id}_genomad/**/*_summary.tsv"), emit: summary 
    tuple val(sample_id), path("${sample_id}_genomad"), emit: genomad_dir
    
    when:
    contigs.size() > 0
    
    script:
    """
    export CUDA_VISIBLE_DEVICES=""  

    # Create output directory
    mkdir -p ${sample_id}_genomad
    
    # Run GeNomad
    genomad end-to-end \\
        --cleanup \\
        --splits 8 \\
        ${contigs} \\
        ${sample_id}_genomad \\
        ${params.genomad_db}
    
    # Validate output
    if [ ! -f "${sample_id}_genomad/${contigs.baseName}_summary/${contigs.baseName}_virus.fna" ]; then
        echo "Warning: No viruses found by GeNomad for ${sample_id}"
        touch "${sample_id}_genomad/${contigs.baseName}_summary/${contigs.baseName}_virus.fna"
    fi
    """
}

process RUN_CHECKV {
    tag "$sample_id"
    label 'process_medium'
    conda 'checkv'
    
    publishDir "${params.outdir}/${sample_id}/checkv", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(virus_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_checkv/viruses.fna"), emit: checkv_fasta
    tuple val(sample_id), path("${sample_id}_checkv"), emit: checkv_dir
    
    when:
    virus_fasta.size() > 0
    
    script:
    """
    export CHECKVDB="${params.checkv_db}"
    # Create output directory
    mkdir -p ${sample_id}_checkv
    
    # Check if input file is not empty
    if [ -s "${virus_fasta}" ]; then
        # Run CheckV
        checkv end_to_end \\
            ${virus_fasta} \\
            ${sample_id}_checkv \\
            -t ${task.cpus}
        
        # Rename output for consistency
        if [ -f "${sample_id}_checkv/proviruses.fna" ] && [ -f "${sample_id}_checkv/viruses.fna" ]; then
            cat ${sample_id}_checkv/proviruses.fna ${sample_id}_checkv/viruses.fna > ${sample_id}_checkv/viruses_combined.fna
            mv ${sample_id}_checkv/viruses_combined.fna ${sample_id}_checkv/viruses.fna
        elif [ -f "${sample_id}_checkv/proviruses.fna" ]; then
            mv ${sample_id}_checkv/proviruses.fna ${sample_id}_checkv/viruses.fna
        elif [ ! -f "${sample_id}_checkv/viruses.fna" ]; then
            touch ${sample_id}_checkv/viruses.fna
        fi
    else
        echo "Warning: Empty virus FASTA file for ${sample_id}"
        mkdir -p ${sample_id}_checkv
        touch ${sample_id}_checkv/viruses.fna
    fi
    """
}

process RUN_IPHOP {
    tag "$sample_id"
    label 'process_high'
    conda 'iphop'
    
    publishDir "${params.outdir}/${sample_id}/iphop", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(checkv_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_iphop"), emit: iphop_dir
    
    when:
    checkv_fasta.size() > 0
    
    script:
    """
    # Check if input file has content
    if [ -s "${checkv_fasta}" ]; then
        # Run iPHoP
        iphop predict \\
            --fa_file ${checkv_fasta} \\
            --db_dir ${params.iphop_db} \\
            --out_dir ${sample_id}_iphop \\
            --num_threads ${task.cpus}
    else
        echo "Warning: Empty CheckV FASTA file for ${sample_id}"
        mkdir -p ${sample_id}_iphop
        echo "No predictions - empty input file" > ${sample_id}_iphop/Host_prediction_to_genus_m90.csv
    fi
    """
}

process RUN_BACPHLIP {
    tag "$sample_id"
    label 'process_medium'
    conda 'bacphlip'
    
    publishDir "${params.outdir}/${sample_id}/bacphlip", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(checkv_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_bacphlip_results.tsv"), emit: bacphlip_results
    
    when:
    checkv_fasta.size() > 0
    
    script:
    """
    # Check if input file has content
    if [ -s "${checkv_fasta}" ]; then
        # Run BACPHLIP
        bacphlip \\
            --multi_fasta ${checkv_fasta} \\
            --local_hmmsearch ${params.bacphlip_hmmsearch} \\
            --threads ${task.cpus} \\
            --output_dir bacphlip_temp
        
        # Move results to expected location
        if [ -f "bacphlip_temp/bacphlip_results.tsv" ]; then
            mv bacphlip_temp/bacphlip_results.tsv ${sample_id}_bacphlip_results.tsv
        else
            echo -e "sequence_id\\tprediction\\tconfidence" > ${sample_id}_bacphlip_results.tsv
            echo "Warning: BACPHLIP produced no results for ${sample_id}" >> ${sample_id}_bacphlip_results.tsv
        fi
    else
        echo "Warning: Empty CheckV FASTA file for ${sample_id}"
        echo -e "sequence_id\\tprediction\\tconfidence" > ${sample_id}_bacphlip_results.tsv
        echo "No predictions - empty input file" >> ${sample_id}_bacphlip_results.tsv
    fi
    """
}

process RUN_VCONTACT2 {
    tag "$sample_id"
    label 'process_high'
    conda 'vcontact2'
    module 'prodigal/2.6.3-GCCcore-11.3.0'
    
    publishDir "${params.outdir}/${sample_id}/vcontact2", mode: params.publish_dir_mode
    
    input:
    tuple val(sample_id), path(checkv_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_vcontact_out"), emit: vcontact_dir
    
    when:
    checkv_fasta.size() > 0
    
    script:
    """
    # Check if input file has content
    if [ -s "${checkv_fasta}" ]; then
        # Create working directory
        mkdir -p ${sample_id}_vcontact_temp
        
        # Run Prodigal to predict genes
        prodigal \\
            -i ${checkv_fasta} \\
            -a ${sample_id}_vcontact_temp/${sample_id}_proteins.faa \\
            -o ${sample_id}_vcontact_temp/${sample_id}_genes.gff \\
            -p meta
        
        # Create gene2genome mapping
        vcontact2_gene2genome \\
            -p ${sample_id}_vcontact_temp/${sample_id}_proteins.faa \\
            -o ${sample_id}_vcontact_temp/gene2genome.csv \\
            -s 'Prodigal-FAA'
        
        # Run vConTACT2
        vcontact2 \\
            --raw-proteins ${sample_id}_vcontact_temp/${sample_id}_proteins.faa \\
            --rel-mode 'Diamond' \\
            --proteins-fp ${sample_id}_vcontact_temp/gene2genome.csv \\
            --db 'ProkaryoticViralRefSeq211-Merged' \\
            --pcs-mode MCL \\
            --vcs-mode ClusterONE \\
            --c1-bin ${params.cluster_one_jar} \\
            --output-dir ${sample_id}_vcontact_out \\
            --threads ${task.cpus}
        
        # Clean up temporary files
        rm -rf ${sample_id}_vcontact_temp
    else
        echo "Warning: Empty CheckV FASTA file for ${sample_id}"
        mkdir -p ${sample_id}_vcontact_out
        echo "No analysis - empty input file" > ${sample_id}_vcontact_out/analysis_log.txt
    fi
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Create input channel from FASTA files
    input_ch = Channel
        .fromPath("${params.input_dir}/*.fasta")
        .map { file -> 
            def sample_id = file.baseName
            return tuple(sample_id, file)
        }
    
    // Check if any input files were found
    input_ch
        .ifEmpty { error "No FASTA files found in ${params.input_dir}/*.fasta" }
        .view { sample_id, file -> 
            log.info "Processing sample: ${sample_id} (${file})"
            return null
        }
    
    // Run GeNomad
    RUN_GENOMAD(input_ch)
    
    // Run CheckV on GeNomad results
    RUN_CHECKV(RUN_GENOMAD.out.virus_fasta)
    
    // Branch CheckV results for parallel processing
    checkv_results = RUN_CHECKV.out.checkv_fasta
    
    // Run parallel analyses
    RUN_IPHOP(checkv_results)
    RUN_BACPHLIP(checkv_results)
    RUN_VCONTACT2(checkv_results)
    
    // Log completion
    RUN_IPHOP.out.iphop_dir
        .concat(RUN_BACPHLIP.out.bacphlip_results)
        .concat(RUN_VCONTACT2.out.vcontact_dir)
        .collect()
        .view { 
            log.info "Pipeline completed successfully!"
            log.info "Results available in: ${params.outdir}"
            return null
        }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    =========================================
     VIROME PIPELINE COMPLETED
    =========================================
    Status:     ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:   ${workflow.duration}
    CPU hours:  ${workflow.stats.computeTimeFmt}
    Work dir:   ${workflow.workDir}
    Exit status: ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: 'None'}
    =========================================
    Results directory: ${params.outdir}
    
    Generated outputs:
    - GeNomad results: virus identification and annotation
    - CheckV results: virus quality assessment  
    - iPHoP results: host prediction
    - BACPHLIP results: lifestyle prediction
    - vConTACT2 results: viral clustering
    
    Reports available in: ${params.tracedir}/
    =========================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    =========================================
     PIPELINE EXECUTION FAILED
    =========================================
    Error message: ${workflow.errorMessage}
    Exit status: ${workflow.exitStatus}
    Work directory: ${workflow.workDir}
    
    Check the log files in the work directory for more details.
    =========================================
    """.stripIndent()
}
