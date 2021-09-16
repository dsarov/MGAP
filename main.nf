#!/usr/bin/env nextflow

/*
 *
 *  Pipeline            ARDaP
 *  Version             1.8.2
 *  Description         Antimicrobial resistance detection and prediction from WGS
 *  Authors             Derek Sarovich, Erin Price, Danielle Madden, Eike Steinig
 *
 */

log.info """
===============================================================================
                           NF-MGAP
                             v2.0
================================================================================

Optional Parameters:

    --fastq      Input PE read file wildcard (default: *_{1,2}.fastq.gz)

                 Currently this is set to $params.fastq
    --ref        ReferenceAlignment

    --executor   Change this flag for running in a HPC scheduler environment.
                 Default behavior is to run without a scheduler but a
                 wide range of schedulers are supported with nextflow.
                 Some of the supported schedulers include sge, pbs, pbspro,
                 slurm, lsf, moab, nqsii. For a full list please visit the
                 nextflow documentation

                 Currently executor is set to $params.executor



If you want to make changes to the default `nextflow.config` file
clone the workflow into a local directory and change parameters
in `nextflow.config`:

    nextflow clone dsarov/mgap outdir/

Update to the local cache of this workflow:

    nextflow pull dsarov/mgap

==================================================================
==================================================================
"""
//find ref and species specific databases

species=params.species
database_config_file="${baseDir}/Databases/Database.config"
ref_proc1="grep -w ${species} ${database_config_file}".execute()
ref_proc2="cut -f3".execute()
ref_proc1 | ref_proc2
ref_proc2.waitFor()
reftmp="${ref_proc2.text}"
ref=reftmp.trim()

database_proc1="grep -w ${species} ${database_config_file}".execute()
database_proc2="cut -f2".execute()
database_proc1 | database_proc2
database_proc2.waitFor()
databasetmp="${database_proc2.text}"
database=databasetmp.trim()
//println "${database}"
//println "${ref}"

params.reference="${baseDir}/Databases/${database}/${ref}"
params.resistance_db="${baseDir}/Databases/${database}/${database}.db"
params.card_db="${baseDir}/Databases/${database}/${database}_CARD.db"
params.snpeff="${database}"

fastq = Channel
  .fromFilePairs("${params.fastq}", flat: true)
	.ifEmpty { exit 1, """ Input read files could not be found.
Have you included the read files in the current directory and do they have the correct naming?
With the parameters specified, ARDaP is looking for reads named ${params.fastq}.
To fix this error either rename your reads to match this formatting or specify the desired format
when initializing ARDaP e.g. --fastq "*_{1,2}_sequence.fastq.gz"

"""
}

reference_file = file(params.reference)
if( !reference_file.exists() ) {
  exit 1, """
ARDaP can't find the reference file.
It is currently looking for this file --> ${params.reference}
Please check that this reference exists here --> ${baseDir}/Databases/${database}/${ref}
If this file doesn't exist either ARDaP is not configured to run with this reference/species
or there was an error during the installation process and ARDaP needs to be re-installed
"""
}

/*
======================================================================
      Part 1: create reference indices, dict files and bed files
======================================================================
*/

process IndexReference {

        label "index"

        input:
        file reference from reference_file

        output:
        file "ref.*" into ref_index_ch
        file "${reference}.fai" into ref_fai_ch1
        file "${reference.baseName}.dict" into ref_dict_ch1
        file "${reference}.bed" into refcov_ch

        script:
        if (ref!="none")
        """
        contig_count=`grep -c '>' ${ref}.fasta`
        echo -e "Joining contigs for ABACAS\n"
        if [ $contig_count == 1 ]; then
          mv ${ref}.fasta ${ref}ABACAS.fasta
        else
          perl $baseDir/bin/joinMultifasta.pl ${ref}.fasta ${ref}ABACAS.fasta"
        fi
        """
}



/*
=======================================================================
   Part 2A: Trim reads with light quality filter and remove adapters
=======================================================================
*/

process Trimmomatic {

    label "trimmomatic"
    tag {"$id"}

    input:
    set id, file(forward), file(reverse) from fastq

    output:
    set id, "${id}_1.fq.gz", "${id}_2.fq.gz" into downsample

    """
    trimmomatic PE -threads $task.cpus ${forward} ${reverse} \
    ${id}_1.fq.gz ${id}_1_u.fq.gz ${id}_2.fq.gz ${id}_2_u.fq.gz \
    ILLUMINACLIP:${baseDir}/resources/trimmomatic/all_adapters.fa:2:30:10: \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
    rm ${id}_1_u.fq.gz ${id}_2_u.fq.gz
    """
}

/*
=======================================================================
   Part 2A: Trim reads with light quality filter and remove adapters
=======================================================================
*/

process Assembly {

  label "Assembly"
  tag { "$id" }
  publishDir "./Outputs/", mode: 'copy', pattern: "*final.fasta", overwrite: true


      input:
      file reference from reference_file
      set id, "${id}_1.fq.gz", "${id}_2.fq.gz" from downsample

      output:
      set id, file("${id}_final.fasta")


      script:
      """
      bash assemble.sh ${id} ${reference} ${baseDir}

      """

}

workflow.onComplete {
	println ( workflow.success ? "\nDone! Result files are in --> ./Outputs\n \
  Antibiotic resistance reports are in --> ./Outputs/AbR_reports\n \
  If further analysis is required, bam alignments are in --> ./Outputs/bams\n \
  Phylogenetic tree and annotated merged variants are in --> ./Outputs/Phylogeny_and_annotation\n \
  Individual variant files are in --> ./Outputs/Variants/VCFs\n" \
  : "Oops .. something went wrong" )
}
