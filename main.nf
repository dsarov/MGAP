#!/usr/bin/env nextflow

/*
 *
 *  Pipeline            MGAP
 *  Version             v2.0
 *  Description         Microbial Genome Assembly Pipeline
 *  Authors             Derek Sarovich, Erin Price,
 *
 */

log.info """
================================================================================
                                    NF-MGAP
                                     v2.1
================================================================================

Optional Parameters:

    --fastq      Input PE read file wildcard (default: *_{1,2}.fastq.gz)

                 Currently this is set to $params.fastq

    --ref        Reference file used for reference assisted assembly using
                 ABACAS. For best results please set this to a closely related
                 reference (i.e. same species and sequence type is ideal)

                 Currently ref is set to $params.ref

    --kraken     Kraken2 can be used to filter raw sequence reads if multiple
                 species contamination is suspected. To use this feature, please
                 set this flag to a specific genus and/or species. E.g. to filter
                 all reads that match Pseudomonas or "Pseudomonas aeruginosa". If
                 specifying species, please use quotes to correctly handle the space
                 between genus and species.

                 Currrently kraken is set to $params.kraken

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
  ARDaP can't find the reference file.
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
        file "ref.*" into ref_index_ch

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
   Part 2B: Optionally filter reads for specific species
=======================================================================
*/
if (params.kraken) {
 process Kraken {

    label "kraken"
    tag {"$id"}

    input:
    set id, file(forward), file(reverse) from kraken

    output:
    set id, "${id}_1.fq.gz", "${id}_2.fq.gz" into assemble

    script:
    """
    kraken2 --db ${KRAKEN_DB} --threads ${task.cpus} --use-names --output ${id}.output --paired ${id}_1.fq.gz ${id}_2.fq.gz
    while read line; do
      echo ${line} >> read.lst
    done < <(grep "${params.kraken}" ${id}.output | awk '{ print \$2 }')
    seqtk subseq ${id}_1.fq.gz read.lst > out_1.fq
    seqtk subseq ${id}_2.fq.gz read.lst > out_2.fq
    gzip out_1.fq
    gzip out_2.fq
    rm ${id}_1.fq.gz ${id}_2.fq.gz
    mv out_1.fq.gz ${id}_1.fq.gz
    mv out_2.fq.gz ${id}_2.fq.gz
    """
 }
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
      bash assemble.sh ${id} ${baseDir} $task.cpus no
      """

}
} else {
  process Assembly_no_ref {

   label "assembly"
   tag { "$id" }
   publishDir "./Outputs/", mode: 'copy', pattern: "*final.fasta", overwrite: true


       input:
       set id, "${id}_1.fq.gz", "${id}_2.fq.gz" from assemble
       file "ref.ABACAS" from ref_index_ch

       output:
       set id, file("${id}_final.fasta")


       script:
       """
       bash assemble.sh ${id} ${baseDir} $task.cpus no none
       """

 }
}
