#!/bin/bash


seq=$1
ref=$2
baseDir=$3

VelvOpt="${baseDir}/bin/velvet_1.2.10/contrrib/VelvetOptimiser-2.2.4/VelvetOptimiser.pl"
IMAGE="${baseDir}/bin/IMAGE_version2.4"

##starting and ending kmer for velvet optimiser
START_KMER=53	
END_KMER=75


##########################################################################
###                                                                    ###
###                          VELVET + OPTIMISER                        ###  
###                             WITH TRIMMED                           ###
###                                                                    ###
##########################################################################


echo -e "now running velvet optimiser with the following parameters\n"
echo -e "starting kmer = $START_KMER\n"
echo -e "ending kmer = $END_KMER\n"
${VelvOpt} -o -scaffolding yes -min_contig_lgth 1000 -s $START_KMER -e $END_KMER -f -shortPaired -fastq.gz ${seq}_merged.fastq.gz -t $NCPUS
mv ${seq}/velvTRIM/auto_data_*/contigs.fa ${seq}_velvet.scaff.fasta



##########################################################################
###                                                                    ###
###                            GAPFILLER                               ###
###                                                                    ###
##########################################################################

    echo -e "${seq}_Gapfiller\tbwa\t${seq}_1.fastq\t${seq}_2.fastq\t500\t0.25\tFR" > ${seq}/Gapfiller.txt
     ${seq}/ "perl $GAPFILLER/GapFiller.pl -l Gapfiller.txt -s ${seq}/${seq}_velvet.scaff.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b Velv_scaff"
    mv ${seq}/Velv_scaff/Velv_scaff.gapfilled.final.fa ${seq}/${seq}_velvet.fasta
	rm -rf ${seq}/Velv_scaff/


##########################################################################
###                                                                    ###
###                             ABACAS                                 ###
###                                                                    ###
##########################################################################

if [ "$contig_count" != 1 -a "$ref" != "none" ]; then	
    if [ ! -s ${seq}/${seq}mapped.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s ${seq}/${seq}_icorn.fasta -a ! -s ${seq}/${seq}_out.fasta ]; then
       perl $ABACAS -m -b -r ${ref}ABACAS.fasta -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped
      echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
      cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
	fi
fi
if [ "$contig_count" == 1 -a "$ref" != "none" ]; then
    if [ ! -s ${seq}/${seq}mapped.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s ${seq}/${seq}_icorn.fasta -a ! -s ${seq}/${seq}_out.fasta ]; then
       perl $ABACAS -m -b -r $PBS_O_WORKDIR/${ref}ABACAS.fasta -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped
      echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
      cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
	fi
fi
if [ "$ref" == "none" ]; then
  mv ${seq}/${seq}_velvet.fasta ${seq}/${seq}mapnunmap.fasta
fi


##########################################################################
###                                                                    ###
###                             IMAGE                                  ###
###                                                                    ###
##########################################################################

## include test for PAGIT assembly
if [ ! -s ${seq}/${seq}_IMAGE2_out.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s ${seq}/${seq}_icorn.fasta ]; then
   perl $IMAGE/image.pl -scaffolds ${seq}mapnunmap.fasta -prefix ${seq} -iteration 1 -all_iteration 3 -dir_prefix ite -kmer 81
   perl $IMAGE/restartIMAGE.pl ite3 71 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite6 61 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite9 51 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite12 41 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite15 31 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite18 21 3 partitioned
   $mv ite21/new.fa ${seq}_IMAGE2_out.fasta"
  # ${seq}/ite18 "perl $IMAGE/contigs2scaffolds.pl new.fa new.read.placed 300 500 scaffolds

  perl $IMAGE/image_run_summary.pl ite > ${seq}/IMAGE2.summary

 # mv ${seq}/ite18/scaffolds.fa ${seq}/${seq}_IMAGE2_out.fasta
  ## there are no quality control steps here to determine the best assembly. The program assumes that ite15 contains the best assembly.
  ## For Bp this is probably the case but may differ with other organisms  
 rm -rf ${seq}/ite*
 rm ${seq}/partitioned_1.fastq
 rm ${seq}/partitioned_2.fastq
 rm ${seq}/image.read.placed
 rm ${seq}/image.contigs.fa
fi  

##########################################################################
###                                                                    ###
###                             SSPACE                                 ###
###                                                                    ###
##########################################################################

  
if [ ! -s ${seq}/${seq}SSPACE.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then

  #TODO need to write SSPACE version check for different library file. The below is for SSPACEv3.0
  
  #echo -e "${seq}SSPACE\tbowtie\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" > ${seq}/library.txt
  
  #For SSPACE v2.0 basic
  echo -e "${seq}SSPACE\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" > ${seq}/library.txt
  
   ${seq} "perl $SSPACE -l ${seq}/library.txt -s ${seq}_IMAGE2_out.fasta"
  mv ${seq}/standard_output.final.scaffolds.fasta ${seq}/${seq}SSPACE.fasta
  rm -rf ${seq}/pairinfo
  rm -rf ${seq}/intermediate_results
  rm -rf ${seq}/bowtieoutput
  rm -rf ${seq}/reads
  rm ${seq}/standard_output.final.evidence
  rm ${seq}/standard_output.logfile.txt
  #mv ${seq}/standard_output.summary.txt ${seq}/SSPACE.summary.txt ##standard_output.summaryfile.txt is the correct name for this output
fi



### SSPACE test ############
## This will skip the next step if SSPACE doesn't find anything to scaffold in your assembly, which caused a gapfiller crash



if [ -s ${seq}/standard_output.summaryfile.txt ]; then
  SSPACE_test=`grep 'Total number of N' ${seq}/standard_output.summaryfile.txt |tail -n1 |awk '{print $6}'`
  if [ "$SSPACE_test" == 0 ]; then
   cp ${seq}/${seq}SSPACE.fasta ${seq}/${seq}_gap2.fasta
  fi
fi


##########################################################################
###                                                                    ###
###                            GAPFILLER 2                             ###
###  This step is skipped is SSPACE doesn't find anything to scaffold  ###
###                                                                    ###
##########################################################################
if [ ! -s ${seq}/${seq}_gap2.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
     ${seq}/ "perl $GAPFILLER/GapFiller.pl -l Gapfiller.txt -s ${seq}/${seq}SSPACE.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b SSPACE_scaff"
    mv ${seq}/SSPACE_scaff/SSPACE_scaff.gapfilled.final.fa ${seq}/${seq}_gap2.fasta
	rm -rf ${seq}/SSPACE_scaff/
fi

##########################################################################
###                                                                    ###
###                Remove contigs <1kb and image cleanup               ###
###                                                                    ###
##########################################################################
if [ ! -s ${seq}/${seq}_pilon.fasta -a "$long" == "no" ]; then
   $CONVERT_PROJECT -f fasta -t fasta -x 1000 -R Contig ${seq}/${seq}_gap2.fasta ${seq}/${seq}_pilon
  echo -e "Project has been filtered to remove contigs less than 1kb in size \n"
  
fi
if [ ! -s ${seq}/${seq}_pilon.fasta -a "$long" == "yes" ]; then
    mv ${seq}/${seq}_gap2.fasta ${seq}/${seq}_pilon.fasta
	echo -e "Project includes all contigs including <1kb in size\n"
fi	

if [ -d ${seq}/ite12 -a -s ${seq}/${seq}_pilon.fasta ]; then
  rm -rf ite*
fi


##########################################################################
###                                                                    ###
###                                PILON                               ###
###                                                                    ###
##########################################################################

#create bam file before running pilon
if [ ! -s ${seq}/${seq}_pilon.fasta.bwt ]; then
  bwa index ${seq}/${seq}_pilon.fasta
  else
  echo "Found ref index for Pilon"
fi
if [ ! -s ${seq}/${seq}.sam ]; then
    bwa mem -R '@RG\tID:Assembly\tSM:${seq}\tPL:ILLUMINA' -a -t 4 ${seq}/${seq}_pilon.fasta ${seq}/${seq}_1.fastq ${seq}/${seq}_2.fastq > ${seq}/${seq}.sam
  else
    echo "Found bam file for pilon"  
fi
if [ ! -s ${seq}/${seq}.bam ]; then
  	    samtools view -h -b -@ 1 -q 1 -o ${seq}.bam.tmp ${seq}.sam && samtools sort -@ 1 -o ${seq}.bam ${seq}.bam.tmp
		rm ${seq}.bam.tmp ${seq}.sam
fi
if [ ! -s ${seq}/${seq}.bam.bai ]; then
    samtools index ${seq}/${seq}.bam
fi

if [ ! -s ${seq}/pilon.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
   $JAVA -jar pilon --genome ${seq}/${seq}_pilon.fasta --frags ${seq}/${seq}.bam
fi 
if [ -s ${seq}/pilon.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
  mv ${seq}/pilon.fasta ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta
fi

## cleanup

if [ -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
rm -rf ${PBS_O_WORKDIR}/tmp/${seq}
exit 0
else 
exit 1
fi

