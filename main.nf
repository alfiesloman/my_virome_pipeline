#!/usr/bin/env nextflow

/*
========================================================================================
    BEST-IN-CLASS VIROME ANALYSIS PIPELINE
========================================================================================
    Author: Alfie Sloman
    Tools included:
    - GeNomad: Virus identification and annotation
    - CheckV: Virus quality assessment
    - iPHoP: Host prediction for viruses
    - BACPHLIP: phage lifestyle prediction
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
    
    # Debug: Check input file
    echo "DEBUG: Input file ${contigs} size: \$(wc -c < ${contigs})"
    echo "DEBUG: Input file ${contigs} sequences: \$(grep -c '^>' ${contigs} || echo 0)"
    
    # Run GeNomad
    genomad end-to-end \\
        --cleanup \\
        --splits 8 \\
        ${contigs} \\
        ${sample_id}_genomad \\
        ${params.genomad_db}
    
    # Debug: Check GeNomad output
    virus_file="${sample_id}_genomad/${contigs.baseName}_summary/${contigs.baseName}_virus.fna"
    if [ -f "\$virus_file" ]; then
        echo "DEBUG: GeNomad virus file size: \$(wc -c < \$virus_file)"
        echo "DEBUG: GeNomad virus sequences: \$(grep -c '^>' \$virus_file || echo 0)"
    else
        echo "DEBUG: GeNomad virus file not found, creating empty file"
        mkdir -p "${sample_id}_genomad/${contigs.baseName}_summary"
        touch "\$virus_file"
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
    
    # Debug: Check input file
    echo "DEBUG: CheckV input file ${virus_fasta} size: \$(wc -c < ${virus_fasta})"
    echo "DEBUG: CheckV input sequences: \$(grep -c '^>' ${virus_fasta} || echo 0)"
    
    # Create output directory
    mkdir -p ${sample_id}_checkv
    
    # Check if input file is not empty and has sequences
    if [ -s "${virus_fasta}" ] && [ \$(grep -c '^>' ${virus_fasta} || echo 0) -gt 0 ]; then
        echo "DEBUG: Running CheckV on ${virus_fasta}"
        
        # Run CheckV
        checkv end_to_end \\
            ${virus_fasta} \\
            ${sample_id}_checkv \\
            -t ${task.cpus}
        
        # Debug: Check CheckV outputs
        echo "DEBUG: CheckV output files:"
        ls -la ${sample_id}_checkv/
        
        # Combine proviruses and viruses if both exist
        if [ -f "${sample_id}_checkv/proviruses.fna" ] && [ -f "${sample_id}_checkv/viruses.fna" ]; then
            echo "DEBUG: Combining proviruses and viruses"
            cat ${sample_id}_checkv/proviruses.fna ${sample_id}_checkv/viruses.fna > ${sample_id}_checkv/viruses_combined.fna
            mv ${sample_id}_checkv/viruses_combined.fna ${sample_id}_checkv/viruses.fna
        elif [ -f "${sample_id}_checkv/proviruses.fna" ]; then
            echo "DEBUG: Using proviruses as viruses"
            mv ${sample_id}_checkv/proviruses.fna ${sample_id}_checkv/viruses.fna
        elif [ ! -f "${sample_id}_checkv/viruses.fna" ]; then
            echo "DEBUG: No CheckV output, creating empty file"
            touch ${sample_id}_checkv/viruses.fna
        fi
    else
        echo "DEBUG: Empty or invalid virus FASTA file for ${sample_id}"
        mkdir -p ${sample_id}_checkv
        touch ${sample_id}_checkv/viruses.fna
    fi
    
    # Final debug: Check output file
    echo "DEBUG: Final CheckV output file size: \$(wc -c < ${sample_id}_checkv/viruses.fna)"
    echo "DEBUG: Final CheckV output sequences: \$(grep -c '^>' ${sample_id}_checkv/viruses.fna || echo 0)"
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
    # Debug: Check input file
    echo "DEBUG: iPHoP input file ${checkv_fasta} size: \$(wc -c < ${checkv_fasta})"
    echo "DEBUG: iPHoP input sequences: \$(grep -c '^>' ${checkv_fasta} || echo 0)"
    
    # Check if input file has content and sequences
    if [ -s "${checkv_fasta}" ] && [ \$(grep -c '^>' ${checkv_fasta} || echo 0) -gt 0 ]; then
        echo "DEBUG: Running iPHoP"
        # Run iPHoP
        iphop predict \\
            --fa_file ${checkv_fasta} \\
            --db_dir ${params.iphop_db} \\
            --out_dir ${sample_id}_iphop \\
            --num_threads ${task.cpus}
    else
        echo "DEBUG: Empty CheckV FASTA file for ${sample_id}"
        mkdir -p ${sample_id}_iphop
        echo "No predictions - empty input file" > ${sample_id}_iphop/Host_prediction_to_genus_m90.csv
    fi
    """
}

process RUN_BACPHLIP {
    tag "$sample_id"
    label 'process_medium'
    conda 'bioconda::bacphlip bioconda::hmmer python=3.8 "numpy<1.20" "pandas<1.3"'

    publishDir "${params.outdir}/${sample_id}/bacphlip", mode: params.publish_dir_mode

    input:
        tuple val(sample_id), path(checkv_fasta)

    output:
        tuple val(sample_id), path("${sample_id}_bacphlip_results.tsv"), emit: bacphlip_results
        path "${sample_id}_bacphlip.log"

    when:
        checkv_fasta.size() > 0

    script:
    """
    set -euo pipefail
    
    # Debug: Check input file
    echo "DEBUG: BACPHLIP input file ${checkv_fasta} size: \$(wc -c < ${checkv_fasta})" > ${sample_id}_bacphlip.log
    
    seq_count=\$(grep -c '^>' ${checkv_fasta} 2>/dev/null || echo 0)
    echo "DEBUG: BACPHLIP input sequences: \$seq_count" >> ${sample_id}_bacphlip.log
    
    if [ ! -s "${checkv_fasta}" ] || [ \$seq_count -eq 0 ]; then
        echo "DEBUG: Empty input file or no sequences" >> ${sample_id}_bacphlip.log
        printf "sequence_id\\tprediction\\tconfidence\\n" > ${sample_id}_bacphlip_results.tsv
        exit 0
    fi

    # Determine whether to use --multi_fasta flag based on sequence count
    if [ \$seq_count -gt 1 ]; then
        echo "DEBUG: Running BACPHLIP with --multi_fasta flag (\$seq_count sequences)" >> ${sample_id}_bacphlip.log
        bacphlip \
            -i ${checkv_fasta} \
            -f \
            --multi_fasta \
            2>> ${sample_id}_bacphlip.log
    else
        echo "DEBUG: Running BACPHLIP without --multi_fasta flag (single sequence)" >> ${sample_id}_bacphlip.log
        bacphlip \
            -i ${checkv_fasta} \
            -f \
            2>> ${sample_id}_bacphlip.log
    fi

    # BACPHLIP creates output as inputfile.bacphlip, so rename it
    if [ -f "${checkv_fasta}.bacphlip" ]; then
        mv ${checkv_fasta}.bacphlip ${sample_id}_bacphlip_results.tsv
        echo "DEBUG: Successfully renamed ${checkv_fasta}.bacphlip to ${sample_id}_bacphlip_results.tsv" >> ${sample_id}_bacphlip.log
    else
        echo "ERROR: BACPHLIP failed to create ${checkv_fasta}.bacphlip file" >> ${sample_id}_bacphlip.log
        printf "sequence_id\\tprediction\\tconfidence\\n" > ${sample_id}_bacphlip_results.tsv
    fi
    """
}
process RUN_VCONTACT_PRODIGAL {
    tag "$sample_id"
    label 'process_medium'
    conda 'bioconda::prodigal' // Simple environment!

    // Add the publish line because the results may be useful
    publishDir "${params.outdir}/${sample_id}/prodigal", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(checkv_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_proteins.faa"), path("${sample_id}_genes.gff"), path(checkv_fasta)

    script:
    """
    # Debug: Check input file
    echo "DEBUG: Prodigal input file ${checkv_fasta} size: \$(wc -c < ${checkv_fasta})"
    echo "DEBUG: Prodigal input sequences: \$(grep -c '^>' ${checkv_fasta} || echo 0)"
    
    # Check if input file has sequences
    if [ -s "${checkv_fasta}" ] && [ \$(grep -c '^>' ${checkv_fasta} || echo 0) -gt 0 ]; then
        # Run Prodigal to predict genes
        prodigal \\
            -i ${checkv_fasta} \\
            -a ${sample_id}_proteins.faa \\
            -o ${sample_id}_genes.gff \\
            -p meta
    else
        echo "DEBUG: Empty input, creating empty output files"
        touch ${sample_id}_proteins.faa
        touch ${sample_id}_genes.gff
    fi
    """
}

process RUN_VCONTACT_GENE2GENOME {
    tag "$sample_id"
    label 'process_high'
    // This env just needs the old numpy/pandas. NO diamond.
    conda 'bioconda::vcontact2 python=3.8 "numpy<1.20" "pandas<1.3"'
    
    input:
    tuple val(sample_id), path(proteins), path(genes), path(checkv_fasta)
    
    output:
    tuple val(sample_id), path(proteins), path("gene2genome.csv"), path(checkv_fasta)

    script:
    """
    # Debug: Check input files
    echo "DEBUG: Gene2genome proteins file ${proteins} size: \$(wc -c < ${proteins})"
    echo "DEBUG: Gene2genome protein sequences: \$(grep -c '^>' ${proteins} || echo 0)"
    
    # Check if proteins file has content
    if [ -s "${proteins}" ] && [ \$(grep -c '^>' ${proteins} || echo 0) -gt 0 ]; then
        # This process ONLY runs vcontact2_gene2genome
        vcontact2_gene2genome \\
            -p ${proteins} \\
            -o gene2genome.csv \\
            -s 'Prodigal-FAA'
    else
        echo "DEBUG: Empty proteins file, creating empty gene2genome.csv"
        echo "protein_id,contig_id,keywords" > gene2genome.csv
    fi
    """
}

process RUN_VCONTACT_MAIN {
    tag "$sample_id"
    label 'process_high'
    conda "bioconda::vcontact2 bioconda::diamond mcl=14.137 python=3.9 scipy=1.7.3 pandas=1.3.5 'numpy<1.23'"
    
    publishDir "${params.outdir}/${sample_id}/vcontact2", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(proteins), path(gene2genome), path(checkv_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_vcontact_out"), emit: vcontact_dir

    when:
    checkv_fasta.size() > 0
    
    script:
    """
    # Debug: Check input files
    echo "DEBUG: vConTACT2 input file ${checkv_fasta} size: \$(wc -c < ${checkv_fasta})"
    echo "DEBUG: vConTACT2 input sequences: \$(grep -c '^>' ${checkv_fasta} || echo 0)"
    echo "DEBUG: vConTACT2 proteins file ${proteins} size: \$(wc -c < ${proteins})"
    echo "DEBUG: vConTACT2 gene2genome file ${gene2genome} size: \$(wc -c < ${gene2genome})"
    
    # Check if input file has content and proteins exist
    if [ -s "${checkv_fasta}" ] && [ -s "${proteins}" ] && [ \$(grep -c '^>' ${checkv_fasta} || echo 0) -gt 0 ]; then
        # This process runs the MAIN vcontact2 command
        vcontact2 \\
            --raw-proteins ${proteins} \\
            --rel-mode 'Diamond' \\
            --proteins-fp ${gene2genome} \\
            --db 'ProkaryoticViralRefSeq211-Merged' \\
            --pcs-mode MCL \\
            --vcs-mode ClusterONE \\
            --c1-bin ${params.cluster_one_jar} \\
            --output-dir ${sample_id}_vcontact_out \\
            --threads ${task.cpus}
    else
        echo "DEBUG: Empty input files for ${sample_id}"
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
    
    // Run new THREE-step VCONTACT2 chain
    RUN_VCONTACT_PRODIGAL(checkv_results)
    RUN_VCONTACT_GENE2GENOME(RUN_VCONTACT_PRODIGAL.out)
    RUN_VCONTACT_MAIN(RUN_VCONTACT_GENE2GENOME.out)
    
    // Log completion
    RUN_IPHOP.out.iphop_dir
        .concat(RUN_BACPHLIP.out.bacphlip_results)
        .concat(RUN_VCONTACT_MAIN.out.vcontact_dir) // <-- Make sure this points to RUN_VCONTACT_MAIN
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