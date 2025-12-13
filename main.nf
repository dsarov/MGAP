#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline          MGAP
 * Version           v3.0 (DSL2)
 * Description       Microbial Genome Assembly Pipeline
 */

// Logging/Info block
log.info """
================================================================================
                                    NF-MGAP
                                     v3.0
================================================================================
fastq        : ${params.fastq}
ref          : ${params.ref}
spades       : ${params.spades}
executor     : ${params.executor}
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

process TRIMMOMATIC {
    label "trimmomatic"
    tag "$id"

    input:
    tuple val(id), path(forward), path(reverse)

    output:
    tuple val(id), path("${id}_1.fq.gz"), path("${id}_2.fq.gz"), emit: trimmed_reads

    script:
    // Note: Ensure params.TRIMMOMATIC is defined in config or defaults
    def trim_cmd = params.TRIMMOMATIC ?: 'trimmomatic'
    """
    ${trim_cmd} PE -threads ${task.cpus} ${forward} ${reverse} \
    ${id}_1.fq.gz ${id}_1_u.fq.gz ${id}_2.fq.gz ${id}_2_u.fq.gz \
    ILLUMINACLIP:${projectDir}/resources/trimmomatic/all_adapters.fa:2:30:10: \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

    rm ${id}_1_u.fq.gz ${id}_2_u.fq.gz
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

    // 2. Run Trimmomatic on all reads
    TRIMMOMATIC(reads_ch)

    // 3. Conditional Logic based on Reference existence
    if (params.ref) {

        ref_file = file(params.ref)
        if( !ref_file.exists() ) { error "Reference file not found: ${params.ref}" }

        // Run Indexing
        INDEX_REFERENCE(ref_file)

        // Run Assembly using the Reference and the Trimmed Reads
        ASSEMBLY_WITH_REF(
            TRIMMOMATIC.out.trimmed_reads,
            ref_file,
            INDEX_REFERENCE.out.abacas_ref
        )

    } else {

        // Run Assembly without Reference
        ASSEMBLY_NO_REF(TRIMMOMATIC.out.trimmed_reads)

    }
}
