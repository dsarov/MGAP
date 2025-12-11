#!/usr/bin/env nextflow

/*
 *
 *  Pipeline            MGAP
 *  Version             v2.2
 *  Description         Microbial Genome Assembly Pipeline
 *  Authors             Derek Sarovich and Erin Price.
 *
 */

log.info """
================================================================================
                                    NF-MGAP
                                     v2.2
================================================================================

Optional Parameters:

    --fastq      Input PE read file wildcard (default: *_{1,2}.fastq.gz)

                 Currently this is set to $params.fastq

    --ref        Set to the name of the reference file.
                 Reference file used for reference assisted assembly using
                 ABACAS. For best results please set this to a closely related
                 reference (i.e. same species and sequence type is ideal)

                 Currently ref is set to $params.ref

    --spades     Set to true/false. Spades can alternatively be used for the initial assembly
                 instead of Velvet. This mode is suggested if the initial assembly is
                 very poor or fails to produce an assembly file of the correct size.
                 Most often this is caused by poor quality input data that needs
                 to be cleaned but it is worth a shot.

                Currently spades is set to $params.spades

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

fastq = Channel
  .fromFilePairs("${params.fastq}", flat: true)
	.ifEmpty { exit 1, """ Input read files could not be found.
Have you included the read files in the current directory and do they have the correct naming?
With the parameters specified, MGAP is looking for reads named ${params.fastq}.
To fix this error either rename your reads to match this formatting or specify the desired format
when initializing MGAP e.g. --fastq "*_{1,2}_sequence.fastq.gz"

"""
}

if (params.ref) {
  reference_file = file(params.ref)
  if( !reference_file.exists() ) {
    exit 1, """
  MGAP can't find the reference file.
  It is currently looking for this file --> ${params.ref}
  If this file doesn't exist, please download and copy to the analysis dirrectory
  """
  }
}

/*
======================================================================
      Part 1: join multifasta file into single contig
======================================================================
*/
if (params.ref) {
  process IndexReference {

        label "index"

        input:
        file ref from reference_file

        output:
        file "ref.ABACAS" into ref_index_ch

        """
        contig_count=\$(grep -c '>' ${ref})
        echo -e "Joining contigs for ABACAS\n"
        if [ \${contig_count} == 1 ]; then
          mv ${ref} ref.ABACAS
        else
          perl ${baseDir}/bin/joinMultifasta.pl ${ref} ref.ABACAS
        fi
        """
  }
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

      set id, "${id}_1.fq.gz", "${id}_2.fq.gz" into kraken, assemble

    script:
    """
    $params.TRIMMOMATIC PE -threads ${task.cpus} ${forward} ${reverse} \
    ${id}_1.fq.gz ${id}_1_u.fq.gz ${id}_2.fq.gz ${id}_2_u.fq.gz \
    ILLUMINACLIP:${baseDir}/resources/trimmomatic/all_adapters.fa:2:30:10: \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
    rm ${id}_1_u.fq.gz ${id}_2_u.fq.gz
    """
}


/*
=======================================================================
   Part 3: Run assembly script
=======================================================================
*/
if (params.ref) {
 process Assembly {

  label "assembly"
  tag { "$id" }
  publishDir "./Outputs/", mode: 'copy', pattern: "*final.fasta", overwrite: true

  input:
  file reference from reference_file
  set id, "${id}_1.fq.gz", "${id}_2.fq.gz" from assemble
  file "ref.ABACAS" from ref_index_ch

  output:
  set id, file("${id}_final.fasta")


  script:
  """
  bash assemble.sh ${id} ${baseDir} $task.cpus no ref $params.spades \"$task.memory\" $params.image
  """

}
} else {
  process Assembly_no_ref {

   label "assembly"
   tag { "$id" }
   publishDir "./Outputs/", mode: 'copy', pattern: "*final.fasta", overwrite: true


   input:
   set id, "${id}_1.fq.gz", "${id}_2.fq.gz" from assemble

   output:
   set id, file("${id}_final.fasta")


   script:
   """
   bash assemble.sh ${id} ${baseDir} $task.cpus no none $params.spades \"$task.memory\" $params.image
   """

 }
}
