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
                           NF-ARDaP
                             v1.8.2
================================================================================

Optional Parameters:

    --fastq      Input PE read file wildcard (default: *_{1,2}.fastq.gz)

                 Currently this is set to $params.fastq

    --species    Species specific database for resistance determination
                 (default: Burkholderia_pseudomallei)

                 Currently you are using $params.species

    --assemblies Optionally include a directory of assembled genomes in the
                 analysis. Set this parameter to 'true' if you wish to included
                 assembled genomes and place all assembled genomes in a
                 subdirectory called 'assemblies'. (default: false)

                 Currently mixtures is set to $params.assemblies

    --mixtures   Optionally perform within species mixtures analysis.
                 Set this parameter to 'true' if you are dealing with
                 multiple strains. (default: false)

                 Currently mixtures is set to $params.mixtures

    --size       ARDaP can optionally down-sample your read data to
                 run through the pipeline quicker. (default: 1000000)

                 Currently you are using $params.size

    --phylogeny  Please include if you would like a whole genome
                 phylogeny (FastTree2) and merged annotation files.
                 Note that this may take some time if you have a large
                 number of isolates (default: false)

                 Currently phylogeny is set to $params.phylogeny

    --executor   Change this flag for running in a HPC scheduler environment.
                 Default behavior is to run without a scheduler but a
                 wide range of schedulers are supported with nextflow.
                 Some of the supported schedulers include sge, pbs, pbspro,
                 slurm, lsf, moab, nqsii. For a full list please visit the
                 nextflow documentation

                 Currently executor is set to $params.executor

    --gwas       **Experimental**. If you have a database of GWAS co-ordinates
                 ARDaP can interrogate SNPs and indels across the entire genome
                 to identify novel mutations likely contributing to an antibiotic
                 resistance phenotype. For more information about this feature,
                 please contact the developers.

                 Currently gwas is set to $params.gwas

If you want to make changes to the default `nextflow.config` file
clone the workflow into a local directory and change parameters
in `nextflow.config`:

    nextflow clone dsarov/ardap outdir/

Update to the local cache of this workflow:

    nextflow pull dsarov/ardap

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

if (params.assemblies) {
  assembly_loc = Channel
    .fromPath("${params.assembly_loc}", checkIfExists: true)
    .ifEmpty {"No assembled genomes will be processed"}
    .map { file ->
      def id = file.name.toString().tokenize('.').get(0)
      return tuple(id, file)
    }
}

resistance_database_file = file(params.resistance_db)
if( !resistance_database_file.exists() ) {
  exit 1, "The resistance database file does no exist: ${params.resistance_db}"
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

card_db_file = file(params.card_db)

patient_meta_file = file(params.patientMetaData)
if( !patient_meta_file.exists() ) {
  exit 1, "The specified patient metadata file does not exist: ${params.patientMetaData}"
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

        """
        bwa index -a is -p ref $reference
        samtools faidx $reference
        picard CreateSequenceDictionary R=$reference O=${reference.baseName}.dict
        bedtools makewindows -g ${reference}.fai -w $params.window > ${reference}.bed
        """
}

/*
======================================================================
      Part 1B: create synthetic reads from reference files
======================================================================
*/

if (params.assemblies) {
  process Read_synthesis {

    label "art"
    tag {"$assembly.baseName"}

    input:
    set id, file(assembly) from assembly_loc

    output:
    set id, file("${assembly.baseName}_1_cov.fq.gz"), file("${assembly.baseName}_2_cov.fq.gz") into (alignment_assembly)

    """
    art_illumina -i ${assembly} -p -l 150 -f 30 -m 500 -s 10 -ss HS25 -na -o ${assembly.baseName}_out
    mv ${assembly.baseName}_out1.fq ${assembly.baseName}_1_cov.fq
    mv ${assembly.baseName}_out2.fq ${assembly.baseName}_2_cov.fq
    gzip ${assembly.baseName}_1_cov.fq
    gzip ${assembly.baseName}_2_cov.fq

    """
  }
}

/*
=======================================================================
Part 2: read processing, reference alignment and variant identification
=======================================================================
// Variant calling sub-workflow - basically SPANDx with a tonne of updates

*/

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
              Part 2B: Downsample reads to increase speed
=======================================================================
*/
process Downsample {

    label "ardap_default"
    tag { "$id" }
    // publishDir "./Clean_reads", mode: 'copy', overwrite: false

    input:
    set id, file(forward), file(reverse) from downsample

    output:
    set id, file("${id}_1_cov.fq.gz"), file("${id}_2_cov.fq.gz") into (alignment)

    script:
    if (params.size > 0) {
            """
            seqtk sample -s 11 ${forward} $params.size | gzip - > ${id}_1_cov.fq.gz
            seqtk sample -s 11 ${reverse} $params.size | gzip - > ${id}_2_cov.fq.gz
            """
     } else {
            // Rename files even if not downsampled to channel into Alignment
            """
            mv ${forward} ${id}_1_cov.fq.gz
            mv ${reverse} ${id}_2_cov.fq.gz
            """
      }
}
/*
=======================================================================
        Part 2C: Align reads against the reference with assemblies
=======================================================================
*/
if (params.assemblies) {
  process ReferenceAlignment_assembly {

    label "alignment"
    tag {"$id"}

    input:
    file ref_index from ref_index_ch
    set id, file(forward), file(reverse) from alignment.mix(alignment_assembly)
    file(card_ref) from Channel.fromPath("$baseDir/Databases/CARD/nucleotide_fasta_protein_homolog_model.fasta").collect()
    file card_db_ref from card_db_file

    output:
    set id, file("${id}.bam"), file("${id}.bam.bai") into dup
    set id, file("${id}.CARD_primary_output.txt") into abr_report_card_ch_1

    """
    bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a \
    -t $task.cpus ref ${forward} ${reverse} > ${id}.sam
    samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
    samtools sort -@ 1 -o ${id}.bam ${id}.bam_tmp
    samtools index ${id}.bam
    rm ${id}.sam ${id}.bam_tmp

    bwa index ${card_ref}
    samtools faidx ${card_ref}
    bedtools makewindows -g ${card_ref}.fai -w 90000 > card.coverage.bed
    bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a -t $task.cpus ${card_ref} ${forward} ${reverse} > ${id}.card.sam
    samtools view -h -b -@ 1 -q 1 -o bam_tmp ${id}.card.sam
    samtools sort -@ 1 -o ${id}.card.bam bam_tmp
    samtools index ${id}.card.bam
    bedtools coverage -a card.coverage.bed -b ${id}.card.bam > ${id}.card.bedcov
    bash SQL_queries_CARD.sh ${id} ${card_db_ref} ${baseDir}
    """

  }

} else {
/*
=======================================================================
               Part 2C: Align reads against the reference
=======================================================================
*/
  process ReferenceAlignment {

    label "alignment"
    tag {"$id"}
    publishDir "./Outputs/CARD", mode: 'copy', pattern: "*CARD_primary_output.txt", overwrite: false

    input:
    file ref_index from ref_index_ch
    set id, file(forward), file(reverse) from alignment
    file(card_ref) from Channel.fromPath("$baseDir/Databases/CARD/nucleotide_fasta_protein_homolog_model.fasta").collect()
    file card_db_ref from card_db_file

    output:
    set id, file("${id}.bam"), file("${id}.bam.bai") into dup
    set id, file("${id}.CARD_primary_output.txt") into abr_report_card_ch_1

    """
    bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a \
    -t $task.cpus ref ${forward} ${reverse} > ${id}.sam
    samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
    samtools sort -@ 1 -o ${id}.bam ${id}.bam_tmp
    samtools index ${id}.bam
    rm ${id}.sam ${id}.bam_tmp

    bwa index ${card_ref}
    samtools faidx ${card_ref}
    bedtools makewindows -g ${card_ref}.fai -w 90000 > card.coverage.bed
    bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a -t $task.cpus ${card_ref} ${forward} ${reverse} > ${id}.card.sam
    samtools view -h -b -@ 1 -q 1 -o bam_tmp ${id}.card.sam
    samtools sort -@ 1 -o ${id}.card.bam bam_tmp
    samtools index ${id}.card.bam
    bedtools coverage -a card.coverage.bed -b ${id}.card.bam > ${id}.card.bedcov
    bash SQL_queries_CARD.sh ${id} ${card_db_ref} ${baseDir}
    """

  }
}

/*
=======================================================================
                       Part 2D: De-duplicate bams
=======================================================================
*/
process Deduplicate {

    label "markduplicates"
    tag { "$id" }
    publishDir "./Outputs/bams", mode: 'copy', pattern: "*.bam*", overwrite: false

    input:
    set id, file(bam_alignment), file(bam_index) from dup
    file refcov from refcov_ch
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch_1

    output:
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") into (variantCalling, mixturePindel, variantcallingGVCF_ch)
    set id, file("output.per-base.bed.gz"), file("${id}.depth.txt") into coverageData
    set id, file("${id}.CARD_primary_output.txt") into abr_report_card_ch_2

    """
    gatk --java-options -Xmx${task.memory.toString().replaceAll(/[\sB]/,'')} MarkDuplicates -I "${id}.bam" -O ${id}.dedup.bam --REMOVE_DUPLICATES true \
    --METRICS_FILE ${id}.dedup.txt --VALIDATION_STRINGENCY LENIENT
    samtools index ${id}.dedup.bam

    mosdepth --by ${refcov} output ${id}.dedup.bam
    sum_depth=\$(zcat output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
    total_chromosomes=\$(zcat output.regions.bed.gz | awk '{print \$4}' | wc -l)
    echo "\$sum_depth/\$total_chromosomes" | bc > ${id}.depth.txt
    """
}

/*
=======================================================================
                        Part 2F: Variant identification
=======================================================================
*/
if (params.mixtures) {

  process VariantCallingMixture {

    label "gatk"
    tag { "$id" }
    publishDir "./Outputs/Variants/Annotated", mode: 'copy', pattern: "*.ALL.annotated.mixture.vcf", overwrite: false
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', pattern: "*.PASS.snps.indels.mixed.vcf", overwrite: false

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") from variantCalling
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch_2

    output:
    set id, file("${id}.ALL.annotated.mixture.vcf") into mixtureArdapProcessing
    file("pindel.out_D.vcf") into mixtureDeletionSummary
    file("pindel.out_TD.vcf") into mixtureDuplicationSummary
    set id, file("${id}.CARD_primary_output.txt") into abr_report_card_ch_3
    set id, file("${id}.PASS.snps.indels.mixed.vcf") into variants_publish_ch

    """
    gatk HaplotypeCaller -R ${reference} --I ${id}.dedup.bam -O ${id}.raw.snps.indels.mixed.vcf

    gatk VariantFiltration -R ${reference} -O ${id}.snps.indels.filtered.mixed.vcf -V ${id}.raw.snps.indels.mixed.vcf \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "StandardFilters"

    header=`grep -n "#CHROM" ${id}.snps.indels.filtered.mixed.vcf | cut -d':' -f 1`
    head -n "\$header" ${id}.snps.indels.filtered.mixed.vcf > snp_head
    cat ${id}.snps.indels.filtered.mixed.vcf | grep PASS | cat snp_head - > ${id}.PASS.snps.indels.mixed.vcf

    snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $params.snpeff ${id}.PASS.snps.indels.mixed.vcf > ${id}.ALL.annotated.mixture.vcf

    echo -e "${id}.dedup.bam\t250\tB" > pindel.bam.config
    pindel -f ${reference} -T $task.cpus -i pindel.bam.config -o pindel.out

    rm -f pindel.out_CloseEndMapped pindel.out_INT_final

    for f in pindel.out_*; do
      pindel2vcf -r ${reference} -R ${reference.baseName} -d ARDaP -p \$f -v \${f}.vcf -e 5 -is 15 -as 50000
      snpEff eff -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $params.snpeff \${f}.vcf > \${f}.vcf.annotated
    done

    """
  }

  process MixtureSummariesSQL {

    label "ardap_default"
    tag { "$id" }

    input:
    set id, file(variants) from mixtureArdapProcessing
    file(pindelD) from mixtureDeletionSummary
    file(pindelTD) from mixtureDuplicationSummary
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch_3

    output:
    set id, file("${id}.annotated.ALL.effects") into variants_all_ch
    set id, file("${id}.Function_lost_list.txt") into function_lost_ch1, function_lost_ch2
    set id, file("${id}.deletion_summary_mix.txt") into deletion_summary_mix_ch
    set id, file("${id}.duplication_summary_mix.txt") into duplication_summary_mix_ch
    set id, file("${id}.CARD_primary_output.txt") into abr_report_card_ch_4

    shell:

    '''
    echo 'Effects summary'

    awk '{if (match($0,"ANN=")){print substr($0,RSTART)}}' !{variants} > all.effects.tmp
    awk -F "|" '{ print $4,$10,$11,$15 }' all.effects.tmp | sed 's/c\\.//' | sed 's/p\\.//' | sed 's/n\\.//'> annotated.ALL.effects.tmp
    grep -E "#|ANN=" !{variants} > ALL.annotated.subset.vcf
    gatk VariantsToTable -V ALL.annotated.subset.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -GF AD -GF DP -O ALL.genotypes.subset.table
    tail -n +2 ALL.genotypes.subset.table | awk '{ print $5,$6,$7,$8 }' > ALL.genotypes.subset.table.headerless
    paste annotated.ALL.effects.tmp ALL.genotypes.subset.table.headerless > !{id}.annotated.ALL.effects

    echo 'Identification of high confidence mutations'

    grep '|HIGH|' !{variants} > ALL.func.lost
		awk '{if (match($0,"ANN=")){print substr($0,RSTART)}}' ALL.func.lost > ALL.func.lost.annotations
		awk -F "|" '{ print $4,$11,$15 }' ALL.func.lost.annotations | sed 's/c\\.//' | sed 's/p\\.//' | sed 's/n\\.//'> ALL.func.lost.annotations.tmp
		grep -E "#|\\|HIGH\\|" !{variants} > ALL.annotated.func.lost.vcf
    gatk VariantsToTable -V ALL.annotated.func.lost.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -GF AD -GF DP -O ALL.annotated.func.lost.table
		tail -n +2 ALL.annotated.func.lost.table | awk '{ print $5,$6,$7,$8 }' > ALL.annotated.func.lost.table.headerless
		paste ALL.func.lost.annotations.tmp ALL.annotated.func.lost.table.headerless > !{id}.Function_lost_list.txt

    echo 'Summary of deletions and duplications'

    grep -v '#' !{pindelD} | awk -v OFS="\t" '{ print $1,$2 }' > d.start.coords.list
		grep -v '#' !{pindelD} | gawk 'match($0, /END=([0-9]+);/,arr){ print arr[1]}' > d.end.coords.list
		grep -v '#' !{pindelD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $2 }' > mutant_depth.D
		grep -v '#' !{pindelD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $1+$2 }' > depth.D
		paste d.start.coords.list d.end.coords.list mutant_depth.D depth.D > !{id}.deletion_summary_mix.txt

    grep -v '#' !{pindelTD} | awk -v OFS="\t" '{ print $1,$2 }' > td.start.coords.list
		grep -v '#' !{pindelTD} | gawk 'match($0, /END=([0-9]+);/,arr){ print arr[1]}' > td.end.coords.list
		grep -v '#' !{pindelTD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $2 }' > mutant_depth.TD
		grep -v '#' !{pindelTD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $1+$2 }' > depth.TD
		paste td.start.coords.list td.end.coords.list mutant_depth.TD depth.TD > !{id}.duplication_summary_mix.txt

    '''

  }

} else {
    // Not a mixture
    //To do split GVCF calling when phylogeny isn't called

    process VariantCalling {

      label "gatk"
      tag { "$id" }
      publishDir "./Outputs/Variants/VCFs", mode: 'copy', pattern: "*FAIL*.vcf", overwrite: false
      publishDir "./Outputs/Variants/VCFs", mode: 'copy', pattern: "*PASS*.vcf", overwrite: false
      publishDir "./Outputs/Variants/Annotated", mode: 'copy', pattern: "*annotated*.vcf", overwrite: false

      input:
      file reference from reference_file
      file reference_fai from ref_fai_ch1
      file reference_dict from ref_dict_ch1
      set id, file(dedup_bam), file(dedup_index) from variantCalling
      set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch_2
      set id, file(perbase), file(depth) from coverageData

      output:
      set id, file("${id}.PASS.snps.vcf"), file("${id}.FAIL.snps.vcf") into filteredSNPs
      set id, file("${id}.PASS.indels.vcf"), file("${id}.FAIL.indels.vcf") into filteredindels
      set id, file("${id}.annotated.indel.effects") into annotated_indels_ch
      set id, file("${id}.annotated.snp.effects") into annotated_snps_ch
      set id, file("${id}.Function_lost_list.txt") into function_lost_ch1, function_lost_ch2
      set id, file("${id}.PASS.snps.annotated.vcf") into annotatedSNPs, annotated_snps_ch2
      set id, file("${id}.PASS.indels.annotated.vcf") into annotatedIndels, annotated_indels_ch2
      set id, file("${id}.deletion_summary.txt") into deletion_summary_ch
      set id, file("${id}.duplication_summary.txt") into duplication_summary_ch
      set id, file("${id}.CARD_primary_output.txt") into abr_report_card_ch_3

      script:
      """
      bash VariantCalling.sh ${id} ${reference} ${baseDir} ${params.snpeff}

      """

    }
}
/*
====================================================================
                              Part 3
  These processes will interrogate the SQL databases (except CARD)
  These have been split to run across different flavours of variants
                 so they can be run in parallel
=====================================================================
*/
 if (params.mixtures) {

  process SqlSnpsIndelsMix {

    label "queries"
    tag { "$id" }

    input:
    set id, file("${id}.annotated.ALL.effects") from variants_all_ch
    set id, file("${id}.Function_lost_list.txt") from function_lost_ch1
    file resistance_db from resistance_database_file
    file("patientMetaData.csv") from patient_meta_file
    set id, file("${id}.deletion_summary_mix.txt") from deletion_summary_mix_ch
    set id, file("${id}.duplication_summary_mix.txt") from duplication_summary_mix_ch
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch_4


    output:
    set id, file("${id}.AbR_output.final.txt") into r_report_ch
    file("patientMetaData.csv") into r_report_metadata_ch
    file("patientDrugSusceptibilityData.csv") into r_report_drug_data_ch

    script:
    """
    bash SQL_queries_SNP_indel_mix.sh ${id} ${resistance_db}
    bash SQL_queries_DelDupMix.sh ${id} ${resistance_db}
    bash AbR_reports_mix.sh ${id} ${resistance_db}
    """
  }
} else {
  process SqlSnpsIndels {

    label "genomic_queries"
    tag { "$id" }

    input:
    set id, file("${id}.annotated.indel.effects") from annotated_indels_ch
    set id, file("${id}.annotated.snp.effects") from annotated_snps_ch
    set id, file("${id}.Function_lost_list.txt") from function_lost_ch1
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch_3
    set id, file("${id}.duplication_summary.txt") from duplication_summary_ch
    set id, file("${id}.deletion_summary.txt") from deletion_summary_ch
    file("patientMetaData.csv") from patient_meta_file
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output.final.txt") into r_report_ch
    file("patientMetaData.csv") into r_report_metadata_ch
    file("patientDrugSusceptibilityData.csv") into r_report_drug_data_ch

    script:
    """
    bash SQL_queries_SNP_indel.sh ${id} ${resistance_db}
    bash SQL_queries_DelDup.sh ${id} ${resistance_db}
    bash AbR_reports.sh ${id} ${resistance_db}
    """
  }
}

if (params.gwas) {
  process GWASInterrogate {

    label "genomic_queries"
    tag { "$id" }
    publishDir "./Outputs/AbR_reports", mode: 'copy', overwrite: true

    input:
    set id, file("${id}.PASS.snps.annotated.vcf") from annotated_snps_ch2
    set id, file("${id}.PASS.indels.annotated.vcf") from annotated_indels_ch2
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output.GWAS_SNP.txt")
    set id, file("${id}.AbR_output.GWAS_indel.txt")

    script:
    """
    bash SQL_GWAS.sh ${id} ${resistance_db}
    """
  }
}

process R_report {

  label "report"
  tag { "$id" }
  publishDir "./Outputs/AbR_reports", mode: 'copy', overwrite: true

  input:
  set id, file("${id}.AbR_output.final.txt") from r_report_ch
  file("patientMetaData.csv") from r_report_metadata_ch
  file("patientDrugSusceptibilityData.csv") from r_report_drug_data_ch

  output:
  set id, file("${id}_report.html")
  set id, file("${id}.AbR_output.final.txt")

  script:
  """
  bash Report_html.sh ${species}
  """
}

/*
===========================================================================
= This process will combine all vcf files into a master VCF file
= Clean vcf files are concatenated and converted into a matrix for phylogeny programs
=
===========================================================================
*/

if (params.phylogeny) {
  process VariantCallingGVCF {

    label "gatk"
    tag { "$id" }
    publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: false, pattern: '*.gvcf'

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") from variantcallingGVCF_ch

    output:
    set id, file("${id}.raw.gvcf")
	  file("${id}.raw.gvcf") into gvcf_files

    """
    gatk HaplotypeCaller -R ${reference} -ERC GVCF --I ${id}.dedup.bam -O ${id}.raw.gvcf
    """
  }
  process Master_vcf {
    label "master_vcf"
    tag { "id" }
    publishDir "./Outputs/Master_vcf", mode: 'copy', overwrite: false

    input:
    file("*.raw.gvcf") from gvcf_files.collect()
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1

    output:
    set file("out.filtered.vcf"), file("out.vcf") into snp_matrix_ch

    script:
    """
    bash Master_vcf.sh ${reference.baseName}
    gatk VariantFiltration -R ${reference} -O out.filtered.vcf -V out.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "HaplotypeScoreFilter"
    """
  }
  process snp_matrix {
    label "snp_matrix"
    publishDir "./Outputs/Phylogeny_and_annotation", mode: 'copy', overwrite: false

    input:
    set file(filtered_vcf), file(out_vcf) from snp_matrix_ch

    output:
    file("Ortho_SNP_matrix.nex")
    file("MP_phylogeny.tre")
    file("ML_phylogeny.tre")
    file("All_SNPs_indels_annotated.txt")

    script:
    """
    bash SNP_matrix.sh $params.snpeff ${baseDir}
    """
  }
}

workflow.onComplete {
	println ( workflow.success ? "\nDone! Result files are in --> ./Outputs\n \
  Antibiotic resistance reports are in --> ./Outputs/AbR_reports\n \
  If further analysis is required, bam alignments are in --> ./Outputs/bams\n \
  Phylogenetic tree and annotated merged variants are in --> ./Outputs/Phylogeny_and_annotation\n \
  Individual variant files are in --> ./Outputs/Variants/VCFs\n" \
  : "Oops .. something went wrong" )
}
