#!/bin/bash


seq=$1
ref=$2
baseDir=$3


source "$SCRIPTPATH"/MGAP.config
source "$SCRIPTPATH"/velvet_optimiser.config
source "$SCRIPTPATH"/scheduler.config


##########################################################################
###                                                                    ###
###                            READ HANDLING                           ###
###                                                                    ###
##########################################################################

if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_merged.fastq.gz -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
  log_eval $PBS_O_WORKDIR "$JAVA -jar $TRIM PE -phred33 -threads $NCPUS ${seq}_1_sequence.fastq.gz ${seq}_2_sequence.fastq.gz $PBS_O_WORKDIR/tmp/${seq}/${seq}_1.fastq.gz $PBS_O_WORKDIR/tmp/${seq}/${seq}_1.tmp.fastq.gz $PBS_O_WORKDIR/tmp/${seq}/${seq}_2.fastq.gz $PBS_O_WORKDIR/tmp/${seq}/${seq}_2.tmp.fastq.gz ILLUMINACLIP:$TRIM_DIR/adapters/all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" #comment out to remove TRIM step
  log_eval $PBS_O_WORKDIR "gunzip $PBS_O_WORKDIR/tmp/${seq}/${seq}_1.fastq.gz" #TRIM specific
  log_eval $PBS_O_WORKDIR "gunzip $PBS_O_WORKDIR/tmp/${seq}/${seq}_2.fastq.gz" #TRIM specific
  rm $PBS_O_WORKDIR/tmp/${seq}/${seq}_2.tmp.fastq.gz $PBS_O_WORKDIR/tmp/${seq}/${seq}_1.tmp.fastq.gz #TRIM specific
  #log_eval $PBS_O_WORKDIR "gunzip -c ${seq}_1_sequence.fastq.gz > $PBS_O_WORKDIR/tmp/${seq}/${seq}_1.fastq" #noTRIM
  #log_eval $PBS_O_WORKDIR "gunzip -c ${seq}_2_sequence.fastq.gz > $PBS_O_WORKDIR/tmp/${seq}/${seq}_2.fastq" #noTRIM
  #log_eval $PBS_O_WORKDIR/tmp/${seq} "$SHUFFLE ${seq}_1.fastq ${seq}_2.fastq ${seq}_merged.fastq" #noTRIM
  #log_eval $PBS_O_WORKDIR/tmp/${seq} "gzip ${seq}_merged.fastq" #noTRIM
  log_eval $PBS_O_WORKDIR/tmp/${seq} "$SHUFFLE ${seq}_1.fastq ${seq}_2.fastq ${seq}_merged.fastq" #trim specific for velvet
  
  log_eval $PBS_O_WORKDIR/tmp/${seq} "gzip ${seq}_merged.fastq" #trim specific for velvet
  echo -e "Illumina sequences have been merged for Velvet assembly\n\n"
fi

##########################################################################
###                                                                    ###
###                          VELVET + OPTIMISER                        ###  
###                             WITH TRIMMED                           ###
###                                                                    ###
##########################################################################

if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.scaff.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then #removed test for ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.fasta.bak
  echo -e "now running velvet optimiser with the following parameters\n"
  echo -e "starting kmer = $START_KMER\n"
  echo -e "ending kmer = $END_KMER\n"
  cd $PBS_O_WORKDIR/tmp/${seq}
  mkdir velvTRIM
  log_eval $PBS_O_WORKDIR/tmp/${seq}/velvTRIM "$VELVETOPT -o \"-scaffolding yes -min_contig_lgth 1000\" -s $START_KMER -e $END_KMER -f \"-shortPaired -fastq.gz $PBS_O_WORKDIR/tmp/${seq}/${seq}_merged.fastq.gz\" -t $NCPUS"
  mv $PBS_O_WORKDIR/tmp/${seq}/velvTRIM/auto_data_*/contigs.fa $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.scaff.fasta
  cd $PBS_O_WORKDIR
fi

##########################################################################
###                                                                    ###
###                            GAPFILLER                               ###
###                                                                    ###
##########################################################################
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.scaff.fasta ]; then
    echo -e "${seq}_Gapfiller\tbwa\t${seq}_1.fastq\t${seq}_2.fastq\t500\t0.25\tFR" > $PBS_O_WORKDIR/tmp/${seq}/Gapfiller.txt
    log_eval $PBS_O_WORKDIR/tmp/${seq}/ "perl $GAPFILLER/GapFiller.pl -l Gapfiller.txt -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.scaff.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b Velv_scaff"
    mv $PBS_O_WORKDIR/tmp/${seq}/Velv_scaff/Velv_scaff.gapfilled.final.fa $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.fasta
	rm -rf $PBS_O_WORKDIR/tmp/${seq}/Velv_scaff/
fi

##########################################################################
###                                                                    ###
###                             ABACAS                                 ###
###                                                                    ###
##########################################################################

if [ "$contig_count" != 1 -a "$ref" != "none" ]; then	
    if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}mapped.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_icorn.fasta -a ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_out.fasta ]; then
      log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $ABACAS -m -b -r $PBS_O_WORKDIR/${ref}ABACAS.fasta -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped"
      echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
      cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
	fi
fi
if [ "$contig_count" == 1 -a "$ref" != "none" ]; then
    if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}mapped.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_icorn.fasta -a ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_out.fasta ]; then
      log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $ABACAS -m -b -r $PBS_O_WORKDIR/${ref}ABACAS.fasta -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped"
      echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
      cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
	fi
fi
if [ "$ref" == "none" ]; then
  mv $PBS_O_WORKDIR/tmp/${seq}/${seq}_velvet.fasta $PBS_O_WORKDIR/tmp/${seq}/${seq}mapnunmap.fasta
fi


##########################################################################
###                                                                    ###
###                             IMAGE                                  ###
###                                                                    ###
##########################################################################

## include test for PAGIT assembly
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_IMAGE2_out.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta -a ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_icorn.fasta ]; then
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/image.pl -scaffolds ${seq}mapnunmap.fasta -prefix ${seq} -iteration 1 -all_iteration 3 -dir_prefix ite -kmer 81"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/restartIMAGE.pl ite3 71 3 partitioned"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/restartIMAGE.pl ite6 61 3 partitioned"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/restartIMAGE.pl ite9 51 3 partitioned"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/restartIMAGE.pl ite12 41 3 partitioned"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/restartIMAGE.pl ite15 31 3 partitioned"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $IMAGE/restartIMAGE.pl ite18 21 3 partitioned"
  log_eval $PBS_O_WORKDIR/tmp/${seq} "mv ite21/new.fa ${seq}_IMAGE2_out.fasta"
  #log_eval $PBS_O_WORKDIR/tmp/${seq}/ite18 "perl $IMAGE/contigs2scaffolds.pl new.fa new.read.placed 300 500 scaffolds"
  cd $PBS_O_WORKDIR/tmp/${seq}
  perl $IMAGE/image_run_summary.pl ite > $PBS_O_WORKDIR/tmp/${seq}/IMAGE2.summary
  cd $PBS_O_WORKDIR
 # mv $PBS_O_WORKDIR/tmp/${seq}/ite18/scaffolds.fa $PBS_O_WORKDIR/tmp/${seq}/${seq}_IMAGE2_out.fasta
  ## there are no quality control steps here to determine the best assembly. The program assumes that ite15 contains the best assembly.
  ## For Bp this is probably the case but may differ with other organisms  
 rm -rf $PBS_O_WORKDIR/tmp/${seq}/ite*
 rm $PBS_O_WORKDIR/tmp/${seq}/partitioned_1.fastq
 rm $PBS_O_WORKDIR/tmp/${seq}/partitioned_2.fastq
 rm $PBS_O_WORKDIR/tmp/${seq}/image.read.placed
 rm $PBS_O_WORKDIR/tmp/${seq}/image.contigs.fa
fi  

##########################################################################
###                                                                    ###
###                             SSPACE                                 ###
###                                                                    ###
##########################################################################

  
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}SSPACE.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then

  #TODO need to write SSPACE version check for different library file. The below is for SSPACEv3.0
  
  #echo -e "${seq}SSPACE\tbowtie\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" > $PBS_O_WORKDIR/tmp/${seq}/library.txt
  
  #For SSPACE v2.0 basic
  echo -e "${seq}SSPACE\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" > $PBS_O_WORKDIR/tmp/${seq}/library.txt
  
  log_eval $PBS_O_WORKDIR/tmp/${seq} "perl $SSPACE -l $PBS_O_WORKDIR/tmp/${seq}/library.txt -s ${seq}_IMAGE2_out.fasta"
  mv $PBS_O_WORKDIR/tmp/${seq}/standard_output.final.scaffolds.fasta $PBS_O_WORKDIR/tmp/${seq}/${seq}SSPACE.fasta
  rm -rf $PBS_O_WORKDIR/tmp/${seq}/pairinfo
  rm -rf $PBS_O_WORKDIR/tmp/${seq}/intermediate_results
  rm -rf $PBS_O_WORKDIR/tmp/${seq}/bowtieoutput
  rm -rf $PBS_O_WORKDIR/tmp/${seq}/reads
  rm $PBS_O_WORKDIR/tmp/${seq}/standard_output.final.evidence
  rm $PBS_O_WORKDIR/tmp/${seq}/standard_output.logfile.txt
  #mv $PBS_O_WORKDIR/tmp/${seq}/standard_output.summary.txt $PBS_O_WORKDIR/tmp/${seq}/SSPACE.summary.txt ##standard_output.summaryfile.txt is the correct name for this output
fi



### SSPACE test ############
## This will skip the next step if SSPACE doesn't find anything to scaffold in your assembly, which caused a gapfiller crash



if [ -s $PBS_O_WORKDIR/tmp/${seq}/standard_output.summaryfile.txt ]; then
  SSPACE_test=`grep 'Total number of N' $PBS_O_WORKDIR/tmp/${seq}/standard_output.summaryfile.txt |tail -n1 |awk '{print $6}'`
  if [ "$SSPACE_test" == 0 ]; then
   cp $PBS_O_WORKDIR/tmp/${seq}/${seq}SSPACE.fasta $PBS_O_WORKDIR/tmp/${seq}/${seq}_gap2.fasta
  fi
fi


##########################################################################
###                                                                    ###
###                            GAPFILLER 2                             ###
###  This step is skipped is SSPACE doesn't find anything to scaffold  ###
###                                                                    ###
##########################################################################
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_gap2.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
    log_eval $PBS_O_WORKDIR/tmp/${seq}/ "perl $GAPFILLER/GapFiller.pl -l Gapfiller.txt -s $PBS_O_WORKDIR/tmp/${seq}/${seq}SSPACE.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b SSPACE_scaff"
    mv $PBS_O_WORKDIR/tmp/${seq}/SSPACE_scaff/SSPACE_scaff.gapfilled.final.fa $PBS_O_WORKDIR/tmp/${seq}/${seq}_gap2.fasta
	rm -rf $PBS_O_WORKDIR/tmp/${seq}/SSPACE_scaff/
fi

##########################################################################
###                                                                    ###
###                Remove contigs <1kb and image cleanup               ###
###                                                                    ###
##########################################################################
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta -a "$long" == "no" ]; then
  log_eval $PBS_O_WORKDIR/tmp/${seq}/ "$CONVERT_PROJECT -f fasta -t fasta -x 1000 -R Contig $PBS_O_WORKDIR/tmp/${seq}/${seq}_gap2.fasta $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon"
  echo -e "Project has been filtered to remove contigs less than 1kb in size \n"
  
fi
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta -a "$long" == "yes" ]; then
    mv $PBS_O_WORKDIR/tmp/${seq}/${seq}_gap2.fasta $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta
	echo -e "Project includes all contigs including <1kb in size\n"
fi	

if [ -d $PBS_O_WORKDIR/tmp/${seq}/ite12 -a -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta ]; then
  rm -rf $PBS_O_WORKDIR/tmp/${seq}/ite*
fi


##########################################################################
###                                                                    ###
###                                PILON                               ###
###                                                                    ###
##########################################################################

#create bam file before running pilon
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta.bwt ]; then
  log_eval $PBS_O_WORKDIR "$BWA index $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta"
  else
  echo "Found ref index for Pilon"
fi
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}.sam ]; then
    log_eval $PBS_O_WORKDIR "$BWA mem -R '@RG\tID:Assembly\tSM:${seq}\tPL:ILLUMINA' -a -t 4 $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta $PBS_O_WORKDIR/tmp/${seq}/${seq}_1.fastq $PBS_O_WORKDIR/tmp/${seq}/${seq}_2.fastq > $PBS_O_WORKDIR/tmp/${seq}/${seq}.sam"
  else
    echo "Found bam file for pilon"  
fi
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam ]; then
  	    log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -@ 1 -q 1 -o $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam.tmp $PBS_O_WORKDIR/tmp/${seq}/${seq}.sam && $SAMTOOLS sort -@ 1 -o $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam.tmp"
		rm "$PBS_O_WORKDIR/tmp/${seq}/${seq}.bam.tmp $PBS_O_WORKDIR/tmp/${seq}/${seq}.sam"
fi
if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam.bai ]; then
    log_eval $PBS_O_WORKDIR "$SAMTOOLS index $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam"
fi

if [ ! -s $PBS_O_WORKDIR/tmp/${seq}/pilon.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
  log_eval $PBS_O_WORKDIR/tmp/${seq} "$JAVA -jar $PILON --genome $PBS_O_WORKDIR/tmp/${seq}/${seq}_pilon.fasta --frags $PBS_O_WORKDIR/tmp/${seq}/${seq}.bam"
fi 
if [ -s $PBS_O_WORKDIR/tmp/${seq}/pilon.fasta -a ! -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
  mv $PBS_O_WORKDIR/tmp/${seq}/pilon.fasta ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta
fi

## cleanup

if [ -s ${PBS_O_WORKDIR}/Assemblies/${seq}_final.fasta ]; then
rm -rf ${PBS_O_WORKDIR}/tmp/${seq}
exit 0
else 
exit 1
fi

