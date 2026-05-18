#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline          MGAP
 * Version           v3.5 (DSL2)
 * Description       Microbial Genome Assembly Pipeline
 */

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

// Renames raw reads if trimming is skipped so assembly script can find them
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
    bash ${projectDir}/bin/assemble.sh ${id} ${projectDir} $task.cpus no ${params.ref} $params.spades ${task.memory.toGiga()}
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
    bash ${projectDir}/bin/assemble.sh ${id} ${projectDir} $task.cpus no none $params.spades ${task.memory.toGiga()}
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
    import sys
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

    print("Writing assembly stats...")

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

            n_contigs = len(lengths)
            total_size = sum(lengths)
            max_contig = max(lengths)
            n50 = calculate_n50(lengths)

            out.write(f"{sample},{n_contigs},{total_size},{max_contig},{n50}\\n")
    """
}

/*
======================================================================
    WORKFLOW LOGIC
======================================================================
*/

workflow {


  log.info """
  ================================================================================
                                      NF-MGAP
                                       v3.6
  ================================================================================
  fastq             : ${params.fastq}
  ref               : ${params.ref}
  spades            : ${params.spades}
  megahit           : ${!params.spades}
  executor          : ${params.executor}
  skip trimming     : ${params.notrim}
  ================================================================================
  """

    reads_ch = Channel.fromFilePairs("${params.fastq}", flat: true)
                      .ifEmpty { error "Cannot find read files: ${params.fastq}" }

    // Logic to Handle Trimming vs Renaming
    if (params.notrim) {
        // Skip trimming, rename files via symlink
        RENAME_READS(reads_ch)
        assembly_input_ch = RENAME_READS.out.renamed_reads
    } else {
        // Run FastP
        FASTP(reads_ch)
        assembly_input_ch = FASTP.out.trimmed_reads
    }

    if (params.ref) {
        ref_file = file(params.ref)
        if( !ref_file.exists() ) { error "Reference file not found: ${params.ref}" }

        INDEX_REFERENCE(ref_file)

        ASSEMBLY_WITH_REF(
            assembly_input_ch,
            ref_file,
            INDEX_REFERENCE.out.abacas_ref
        )
        final_assemblies_ch = ASSEMBLY_WITH_REF.out.assembly

    } else {
        ASSEMBLY_NO_REF(assembly_input_ch)
        final_assemblies_ch = ASSEMBLY_NO_REF.out.assembly
    }

    GENERATE_SUMMARY(
        final_assemblies_ch.map{ it[1] }.collect()
    )
}
