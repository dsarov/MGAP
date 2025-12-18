#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline          MGAP
 * Version           v3.10 (DSL2)
 * Description       Microbial Genome Assembly Pipeline (Kraken2 + CheckM1)
 */

// Default parameter values
params.kraken_db = null
params.checkm_data = null // Path to CheckM reference data folder

log.info """
================================================================================
                                    NF-MGAP
                                     v3.10
================================================================================
fastq             : ${params.fastq}
ref               : ${params.ref}
spades            : ${params.spades}
megahit           : ${!params.spades}
skip trimming     : ${params.notrim}
kraken_db         : ${params.kraken_db}
checkm_data       : ${params.checkm_data}
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
    path "${id}.fastp.html", emit: html

    script:
    """
    fastp --in1 ${forward} --in2 ${reverse} --out1 ${id}_1.fq.gz --out2 ${id}_2.fq.gz --thread ${task.cpus} --detect_adapter_for_pe --html ${id}.fastp.html
    """
}

process RENAME_READS {
    label "index"
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
    script:
    """
    kraken2 --db ${k2_db} --threads ${task.cpus} --report ${id}.kraken2_report.txt --paired ${r1} ${r2} --output /dev/null
    """
}

process ASSEMBLY_WITH_REF {
    label "assembly"
    tag "$id"
    publishDir "./Outputs/Assembly", mode: 'copy', pattern: "*final.fasta"
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
    publishDir "./Outputs/Assembly", mode: 'copy', pattern: "*final.fasta"
    input:
    tuple val(id), path(r1), path(r2)
    output:
    tuple val(id), path("${id}_final.fasta"), emit: assembly
    script:
    """
    bash ${projectDir}/bin/assemble.sh ${id} ${projectDir} $task.cpus no none $params.spades ${task.memory.toGiga()}
    """
}

process CHECKM {
    label "checkm"
    tag "$id"
    publishDir "./Outputs/QC/CheckM", mode: 'copy'

    input:
    tuple val(id), path(assembly)
    path checkm_data_dir

    output:
    path "${id}_checkm_results/", emit: folder
    path "${id}_checkm_stats.txt", emit: stats

    script:
    """
    # CheckM1 requires the data root to be set
    checkm data setRoot ${checkm_data_dir}

    # Run the lineage workflow
    checkm lineage_wf \
        -t ${task.cpus} \
        -x fasta \
        --file ${id}_checkm_stats.txt \
        . \
        ${id}_checkm_results
    """
}

process GENERATE_SUMMARY {
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

    if (params.notrim) {
        RENAME_READS(reads_ch)
        assembly_input_ch = RENAME_READS.out.renamed_reads
    } else {
        FASTP(reads_ch)
        assembly_input_ch = FASTP.out.trimmed_reads
    }

    if (params.kraken_db) {
        KRAKEN2(assembly_input_ch, file(params.kraken_db))
    }

    if (params.ref) {
        ASSEMBLY_WITH_REF(assembly_input_ch, file(params.ref))
        final_assemblies_ch = ASSEMBLY_WITH_REF.out.assembly
    } else {
        ASSEMBLY_NO_REF(assembly_input_ch)
        final_assemblies_ch = ASSEMBLY_NO_REF.out.assembly
    }

    if (params.checkm_data) {
        CHECKM(final_assemblies_ch, file(params.checkm_data))
    }

    GENERATE_SUMMARY(final_assemblies_ch.map{ it[1] }.collect())
}
