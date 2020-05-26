#!/bin/bash


##########################################################################
###                                                                    ###
###                           ABACAS REF                               ###
###                                                                    ###
##########################################################################
if [ "$ref" != "none" ]; then
    contig_count=`grep -c '>' ${ref}.fasta`

    if [ $contig_count -gt 1 ]; then
	echo -e "Joining contigs for ABACAS\n"
      perl $SCRIPTPATH/bin/joinMultifasta.pl ${ref}.fasta ${ref}ABACAS.fasta
    fi
    if [ ! -s $PBS_O_WORKDIR/${ref}ABACAS.fasta -a $contig_count == 1 ]; then
       ln -s $PBS_O_WORKDIR/${ref}.fasta $PBS_O_WORKDIR/${ref}ABACAS.fasta
    fi
fi

##########################################################################
###                                                                    ###
###                            READ HANDLING                           ###
###                                                                    ###
##########################################################################

if [ ! -s  ${seq}_merged.fastq.gz -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
   $JAVA -jar $TRIM PE -phred33 -threads $NCPUS ${seq}_1_sequence.fastq.gz ${seq}_2_sequence.fastq.gz  ${seq}_1.fastq.gz  ${seq}_1.tmp.fastq.gz  ${seq}_2.fastq.gz  ${seq}_2.tmp.fastq.gz ILLUMINACLIP:$TRIM_DIR/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
  
   gunzip  ${seq}_1.fastq.gz #TRIM specific
   gunzip  ${seq}_2.fastq.gz #TRIM specific
  
  rm ${seq}_2.tmp.fastq.gz  ${seq}_1.tmp.fastq.gz #TRIM specific

  $SHUFFLE ${seq}_1.fastq ${seq}_2.fastq ${seq}_merged.fastq #trim specific for velvet
  
  gzip ${seq}_merged.fastq #trim specific for velvet
  echo -e "Illumina sequences have been merged for Velvet assembly\n\n"
fi

##########################################################################
###                                                                    ###
###                          VELVET + OPTIMISER                        ###  
###                             WITH TRIMMED                           ###
###                                                                    ###
##########################################################################

if [ ! -s  ${seq}_velvet.scaff.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then #removed test for ! -s  ${seq}_velvet.fasta.bak
  echo -e "now running velvet optimiser with the following parameters\n"
  echo -e "starting kmer = $START_KMER\n"
  echo -e "ending kmer = $END_KMER\n"
  cd $PBS_O_WORKDIR/tmp/${seq}
  mkdir velvTRIM
  $VELVETOPT -o \"-scaffolding yes -min_contig_lgth 1000\" -s $START_KMER -e $END_KMER -f \"-shortPaired -fastq.gz  ${seq}_merged.fastq.gz\" -t $NCPUS
  mv  velvTRIM/auto_data_*/contigs.fa  ${seq}_velvet.scaff.fasta
  
fi

##########################################################################
###                                                                    ###
###                            GAPFILLER                               ###
###                                                                    ###
##########################################################################
if [ ! -s  ${seq}_velvet.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a -s  ${seq}_velvet.scaff.fasta ]; then
    echo -e "${seq}_Gapfiller\tbwa\t${seq}_1.fastq\t${seq}_2.fastq\t500\t0.25\tFR" >  Gapfiller.txt
    log_eval   "perl $GAPFILLER/GapFiller.pl -l Gapfiller.txt -s  ${seq}_velvet.scaff.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b Velv_scaff"
    mv  Velv_scaff/Velv_scaff.gapfilled.final.fa  ${seq}_velvet.fasta
	rm -rf  Velv_scaff/
fi

##########################################################################
###                                                                    ###
###                             ABACAS                                 ###
###                                                                    ###
##########################################################################

if [ "$contig_count" != 1 -a "$ref" != "none" ]; then	
    if [ ! -s  ${seq}mapped.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s  ${seq}_icorn.fasta -a ! -s  ${seq}_out.fasta ]; then
      perl $ABACAS -m -b -r $PBS_O_WORKDIR/${ref}ABACAS.fasta -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped
      echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
      cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
	fi
fi
if [ "$contig_count" == 1 -a "$ref" != "none" ]; then
    if [ ! -s  ${seq}mapped.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s  ${seq}_icorn.fasta -a ! -s  ${seq}_out.fasta ]; then
       perl $ABACAS -m -b -r $PBS_O_WORKDIR/${ref}ABACAS.fasta -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped
      echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
      cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
	fi
fi
if [ "$ref" == "none" ]; then
  mv  ${seq}_velvet.fasta  ${seq}mapnunmap.fasta
fi


##########################################################################
###                                                                    ###
###                             IMAGE                                  ###
###                                                                    ###
##########################################################################

## include test for PAGIT assembly
if [ ! -s  ${seq}_IMAGE2_out.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s  ${seq}_icorn.fasta ]; then
   perl $IMAGE/image.pl -scaffolds ${seq}mapnunmap.fasta -prefix ${seq} -iteration 1 -all_iteration 3 -dir_prefix ite -kmer 81
   perl $IMAGE/restartIMAGE.pl ite3 71 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite6 61 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite9 51 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite12 41 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite15 31 3 partitioned
   perl $IMAGE/restartIMAGE.pl ite18 21 3 partitioned
   mv ite21/new.fa ${seq}_IMAGE2_out.fasta
  cd $PBS_O_WORKDIR/tmp/${seq}
  perl $IMAGE/image_run_summary.pl ite >  IMAGE2.summary
  cd $PBS_O_WORKDIR
 # mv  ite18/scaffolds.fa  ${seq}_IMAGE2_out.fasta
  ## there are no quality control steps here to determine the best assembly. The program assumes that ite15 contains the best assembly.
  ## For Bp this is probably the case but may differ with other organisms  
 rm -rf  ite*
 rm  partitioned_1.fastq
 rm  partitioned_2.fastq
 rm  image.read.placed
 rm  image.contigs.fa
fi  

##########################################################################
###                                                                    ###
###                             SSPACE                                 ###
###                                                                    ###
##########################################################################

  
if [ ! -s  ${seq}SSPACE.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then

  #TODO need to write SSPACE version check for different library file. The below is for SSPACEv3.0
  
  #echo -e "${seq}SSPACE\tbowtie\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" >  library.txt
  
  #For SSPACE v2.0 basic
  echo -e "${seq}SSPACE\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" >  library.txt
  
   perl $SSPACE -l  library.txt -s ${seq}_IMAGE2_out.fasta
  mv  standard_output.final.scaffolds.fasta  ${seq}SSPACE.fasta
  rm -rf  pairinfo
  rm -rf  intermediate_results
  rm -rf  bowtieoutput
  rm -rf  reads
  rm  standard_output.final.evidence
  rm  standard_output.logfile.txt
fi



### SSPACE test ############
## This will skip the next step if SSPACE doesn't find anything to scaffold in your assembly, which caused a gapfiller crash



if [ -s  standard_output.summaryfile.txt ]; then
  SSPACE_test=`grep 'Total number of N'  standard_output.summaryfile.txt |tail -n1 |awk '{print $6}'`
  if [ "$SSPACE_test" == 0 ]; then
   cp  ${seq}SSPACE.fasta  ${seq}_gap2.fasta
  fi
fi


##########################################################################
###                                                                    ###
###                            GAPFILLER 2                             ###
###  This step is skipped is SSPACE doesn't find anything to scaffold  ###
###                                                                    ###
##########################################################################
if [ ! -s  ${seq}_gap2.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
   perl $GAPFILLER/GapFiller.pl -l Gapfiller.txt -s  ${seq}SSPACE.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b SSPACE_scaff
    mv  SSPACE_scaff/SSPACE_scaff.gapfilled.final.fa  ${seq}_gap2.fasta
	rm -rf  SSPACE_scaff/
fi

##########################################################################
###                                                                    ###
###                Remove contigs <1kb and image cleanup               ###
###                                                                    ###
##########################################################################
if [ ! -s  ${seq}_pilon.fasta -a "$long" == "no" ]; then
  $CONVERT_PROJECT -f fasta -t fasta -x 1000 -R Contig  ${seq}_gap2.fasta  ${seq}_pilon
  echo -e "Project has been filtered to remove contigs less than 1kb in size \n"
  
fi
if [ ! -s  ${seq}_pilon.fasta -a "$long" == "yes" ]; then
    mv  ${seq}_gap2.fasta  ${seq}_pilon.fasta
	echo -e "Project includes all contigs including <1kb in size\n"
fi	

if [ -d  ite12 -a -s  ${seq}_pilon.fasta ]; then
  rm -rf  ite*
fi


##########################################################################
###                                                                    ###
###                                PILON                               ###
###                                                                    ###
##########################################################################

#create bam file before running pilon
if [ ! -s  ${seq}_pilon.fasta.bwt ]; then
   $BWA index ${seq}_pilon.fasta
  else
  echo "Found ref index for Pilon"
fi
if [ ! -s  ${seq}.sam ]; then
     $BWA mem -R '@RG\tID:Assembly\tSM:${seq}\tPL:ILLUMINA' -a -t 4  ${seq}_pilon.fasta  ${seq}_1.fastq  ${seq}_2.fastq >  ${seq}.sam
  else
    echo "Found bam file for pilon"  
fi
if [ ! -s  ${seq}.bam ]; then
  	     $SAMTOOLS view -h -b -@ 1 -q 1 -o  ${seq}.bam.tmp  ${seq}.sam && $SAMTOOLS sort -@ 1 -o  ${seq}.bam  ${seq}.bam.tmp
		rm ${seq}.bam.tmp  ${seq}.sam
fi
if [ ! -s  ${seq}.bam.bai ]; then
     $SAMTOOLS index  ${seq}.bam
fi

if [ ! -s  pilon.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
   $JAVA -jar $PILON --genome  ${seq}_pilon.fasta --frags  ${seq}.bam
fi 
if [ -s  pilon.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
  mv  pilon.fasta ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta
fi

## cleanup

if [ -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
rm -rf ${PBS_O_WORKDIR}/tmp/${seq}
exit 0
else 
exit 1
fi

