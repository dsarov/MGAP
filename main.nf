#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline          MGAP
 * Version           v3.8 (DSL2)
 * Description       Microbial Genome Assembly Pipeline (with Kraken2 and CheckM2)
 */

// Default parameter values
params.kraken_db = null
params.checkm2_db = null

log.info """
================================================================================
                                    NF-MGAP
                                     v3.8
================================================================================
fastq             : ${params.fastq}
ref               : ${params.ref}
spades            : ${params.spades}
megahit           : ${!params.spades}
executor          : ${params.executor}
skip trimming     : ${params.notrim}
kraken_db         : ${params.kraken_db}
checkm2_db        : ${params.checkm2_db}
================================================================================
"""

/*
======================================================================
    PROCESS DEFINITIONS
======================================================================
*/

process FASTP {
    label "fastp"
    tag "$id"
    publishDir "./Outputs/QC/Fastp", mode: 'copy', pattern: "*.html"

    input:
    tuple val(id), path(forward), path(reverse)

    output:
    tuple val(id), path("${id}_1.fq.gz"), path("${id}_2.fq.gz"), emit: trimmed_reads
    path "${id}.fastp.json", emit: json
    path "${id}.fastp.html", emit: html

    script:
    """
    fastp \
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

process RENAME_READS {
    label "rename"
    tag "$id"

    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_1.fq.gz"), path("${id}_2.fq.gz"), emit: renamed_reads

    script:
    """
    ln -s ${r1} ${id}_1.fq.gz
    ln -s ${r2} ${id}_2.fq.gz
    """
}

process KRAKEN2 {
    label "kraken"
    tag "$id"
    publishDir "./Outputs/QC/Kraken2", mode: 'copy'

    input:
    tuple val(id), path(r1), path(r2)
    path k2_db

    output:
    path "${id}.kraken2_report.txt", emit: report
    path "${id}.kraken2.output", emit: raw_output

    script:
    """
    kraken2 \
        --db ${k2_db} \
        --threads ${task.cpus} \
        --report ${id}.kraken2_report.txt \
        --paired ${r1} ${r2} \
        --output ${id}.kraken2.output
    """
}

process ASSEMBLY_WITH_REF {
    label "assembly"
    tag "$id"
    publishDir "./Outputs/Assembly", mode: 'copy', pattern: "*final.fasta", overwrite: true

    input:
    tuple val(id), path(r1), path(r2)
    path reference

    output:
    tuple val(id), path("${id}_final.fasta"), emit: assembly

    script:
    """
    bash ${projectDir}/bin/assemble.sh ${id} ${projectDir} $task.cpus no ${reference} $params.spades ${task.memory.toGiga()}
    """
}

process ASSEMBLY_NO_REF {
    label "assembly"
    tag "$id"
    publishDir "./Outputs/Assembly", mode: 'copy', pattern: "*final.fasta", overwrite: true

    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_final.fasta"), emit: assembly

    script:
    """
    bash ${projectDir}/bin/assemble.sh ${id} ${projectDir} $task.cpus no none $params.spades ${task.memory.toGiga()}
    """
}

process CHECKM2 {
    label "checkm2"
    tag "$id"
    publishDir "./Outputs/QC/CheckM2", mode: 'copy'

    input:
    tuple val(id), path(assembly)
    path checkm2_db

    output:
    path "${id}_checkm2_report/*", emit: report

    script:
    """
    checkm2 predict \
        --threads ${task.cpus} \
        --input ${assembly} \
        --output-directory ${id}_checkm2_report \
        --db-path ${checkm2_db}
    """
}

process GENERATE_SUMMARY {
    label "summary"
    publishDir "./Outputs/", mode: 'copy'

    input:
    path assemblies

    output:
    path "assembly_stats.csv"

    script:
    """
    #!/usr/bin/env python3
    import os

    def calculate_n50(lengths):
        sorted_len = sorted(lengths, reverse=True)
        total_len = sum(lengths)
        csum = 0
        for l in sorted_len:
            csum += l
            if csum >= total_len / 2:
                return l
        return 0

    with open("assembly_stats.csv", "w") as out:
        out.write("Sample,Num_Contigs,Total_Size,Max_Contig,N50\\n")
        input_files = "${assemblies}".split()
        for fpath in input_files:
            sample = os.path.basename(fpath).replace('_final.fasta', '')
            lengths = []
            try:
                with open(fpath, 'r') as fasta:
                    seq_len = 0
                    for line in fasta:
                        line = line.strip()
                        if line.startswith(">"):
                            if seq_len > 0: lengths.append(seq_len)
                            seq_len = 0
                        else:
                            seq_len += len(line)
                    if seq_len > 0: lengths.append(seq_len)
            except IOError:
                continue
            if not lengths:
                out.write(f"{sample},0,0,0,0\\n")
                continue
            out.write(f"{sample},{len(lengths)},{sum(lengths)},{max(lengths)},{calculate_n50(lengths)}\\n")
    """
}

/*
======================================================================
    WORKFLOW LOGIC
======================================================================
*/

workflow {
    reads_ch = Channel.fromFilePairs("${params.fastq}", flat: true)
                      .ifEmpty { error "Cannot find read files: ${params.fastq}" }

    // 1. Pre-processing
    if (params.notrim) {
        RENAME_READS(reads_ch)
        assembly_input_ch = RENAME_READS.out.renamed_reads
    } else {
        FASTP(reads_ch)
        assembly_input_ch = FASTP.out.trimmed_reads
    }

    // 2. Kraken2 (Raw Read Screening)
    if (params.kraken_db) {
        k2_db = file(params.kraken_db)
        KRAKEN2(assembly_input_ch, k2_db)
    }

    // 3. Assembly
    if (params.ref) {
        ref_file = file(params.ref)
        ASSEMBLY_WITH_REF(assembly_input_ch, ref_file)
        final_assemblies_ch = ASSEMBLY_WITH_REF.out.assembly
    } else {
        ASSEMBLY_NO_REF(assembly_input_ch)
        final_assemblies_ch = ASSEMBLY_NO_REF.out.assembly
    }

    // 4. CheckM2 (Assembly Quality/Contamination)
    if (params.checkm2_db) {
        cm2_db = file(params.checkm2_db)
        CHECKM2(final_assemblies_ch, cm2_db)
    }

    // 5. Statistics
    GENERATE_SUMMARY(final_assemblies_ch.map{ it[1] }.collect())
}
