#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline          MGAP
 * Version           v3.2 (DSL2)
 * Description       Microbial Genome Assembly Pipeline
 */

// Default parameter to ensure logic works if flag is missing
params.no_trim = false

// Logging/Info block
log.info """
================================================================================
                                    NF-MGAP
                                     v3.2
================================================================================
fastq        : ${params.fastq}
ref          : ${params.ref}
spades       : ${params.spades}
executor     : ${params.executor}
trimming     : ${params.no_trim ? "SKIPPED" : "ENABLED (FastP)"}
================================================================================
"""

/*
======================================================================
    PROCESS DEFINITIONS
======================================================================
*/

process INDEX_REFERENCE {
    label "index"

    input:
    path ref

    output:
    path "ref.ABACAS", emit: abacas_ref

    script:
    """
    contig_count=\$(grep -c '>' ${ref})
    echo -e "Joining contigs for ABACAS\n"
    if [ \${contig_count} == 1 ]; then
        cp ${ref} ref.ABACAS
    else
        perl ${projectDir}/bin/joinMultifasta.pl ${ref} ref.ABACAS
    fi
    """
}

process FASTP {
    label "fastp"
    tag "$id"
    publishDir "./Outputs/QC", mode: 'copy', pattern: "*.html"

    input:
    tuple val(id), path(forward), path(reverse)

    output:
    tuple val(id), path("${id}_1.fq.gz"), path("${id}_2.fq.gz"), emit: trimmed_reads
    path "${id}.fastp.json", emit: json
    path "${id}.fastp.html", emit: html

    script:
    // Allows you to override the path in config (params.FASTP) or defaults to 'fastp'
    def fastp_cmd = params.FASTP ?: 'fastp'
    """
    ${fastp_cmd} \
      --in1 ${forward} \
      --in2 ${reverse} \
      --out1 ${id}_1.fq.gz \
      --out2 ${id}_2.fq.gz \
      --thread ${task.cpus} \
      --detect_adapter_for_pe \
      --length_required 36 \
      --json ${id}.fastp.json \
      --html ${id}.fastp.html
    """
}

process ASSEMBLY_WITH_REF {
    label "assembly"
    tag "$id"
    publishDir "./Outputs/", mode: 'copy', pattern: "*final.fasta", overwrite: true

    input:
    tuple val(id), path(r1), path(r2)
    path reference
    path abacas_ref

    output:
    tuple val(id), path("${id}_final.fasta"), emit: assembly

    script:
    """
    bash ${projectDir}/bin/assemble.sh ${id} ${projectDir} $task.cpus no ${params.ref} $params.spades ${task.memory.toGiga()} $params.image
    """
}

process ASSEMBLY_NO_REF {
    label "assembly"
    tag "$id"
    publishDir "./Outputs/", mode: 'copy', pattern: "*final.fasta", overwrite: true

    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_final.fasta"), emit: assembly

    script:
    // Note: I kept your path as 'assemble.sh' here (vs 'bin/assemble.sh' above)
    // to match your original script, but check if this needs to be standardized.
    """
    bash ${projectDir}/assemble.sh ${id} ${projectDir} $task.cpus no none $params.spades ${task.memory.toGiga()} $params.image
    """
}

/*
======================================================================
    WORKFLOW LOGIC
======================================================================
*/

workflow {
    // 1. Validate inputs
    reads_ch = Channel.fromFilePairs("${params.fastq}", flat: true)
                      .ifEmpty { error "Cannot find read files: ${params.fastq}" }

    // 2. Logic to Handle Trimming (FastP) vs Skipping

    if (params.no_trim) {
        // Skip Process: Pass raw reads directly to assembly
        assembly_input_ch = reads_ch
    } else {
        // Run Process: Run FastP first
        FASTP(reads_ch)
        assembly_input_ch = FASTP.out.trimmed_reads
    }

    // 3. Conditional Logic based on Reference existence
    if (params.ref) {

        ref_file = file(params.ref)
        if( !ref_file.exists() ) { error "Reference file not found: ${params.ref}" }

        // Run Indexing
        INDEX_REFERENCE(ref_file)

        // Run Assembly using the decided input channel (raw or trimmed)
        ASSEMBLY_WITH_REF(
            assembly_input_ch,
            ref_file,
            INDEX_REFERENCE.out.abacas_ref
        )

    } else {

        // Run Assembly without Reference
        ASSEMBLY_NO_REF(assembly_input_ch)

    }
}
