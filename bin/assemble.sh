#!/bin/bash


seq=$1
baseDir=$2
NCPUS=$3
long=$4

source ${baseDir}/configs/dependencies.config
#testing
#baseDir=~/bin/mgap/
#NCPUS=23
#ref=Pa_PA01.fasta

export PERL5LIB=/home/dsarovich/.cpan/build/Perl4-CoreLibs-0.004-0/lib/
##starting and ending kmer for velvet optimiser
START_KMER=53	
END_KMER=75

#need to add chmod +x -r ./mgap/

##########################################################################
###                                                                    ###
###                          VELVET + OPTIMISER                        ###  
###                             WITH TRIMMED                           ###
###                                                                    ###
##########################################################################
echo "Unzipping reads"
gunzip -c ${seq}_1.fq.gz > ${seq}_1.fastq
gunzip -c ${seq}_2.fq.gz > ${seq}_2.fastq
echo "done"


echo "Shuffling sequences"
echo "command = perl ${SHUFFLE} ${seq}_1.fastq ${seq}_2.fastq ${seq}_merged.fastq\n"
perl ${SHUFFLE} ${seq}_1.fastq ${seq}_2.fastq ${seq}_merged.fastq

echo -e "now running velvet optimiser with the following parameters\n"
echo -e "starting kmer = $START_KMER\n"
echo -e "ending kmer = $END_KMER\n"

echo "running velvet optimiser"
echo "command = perl ${VelvOpt} -o \"-scaffolding yes -min_contig_lgth 1000\" -s ${START_KMER} -e ${END_KMER} -f \"-shortPaired -fastq.gz ${seq}_merged.fastq\" -t $NCPUS"
perl ${VelvOpt} -o "-scaffolding yes -min_contig_lgth 1000" -s ${START_KMER} -e ${END_KMER} -f "-shortPaired -fastq.gz ${seq}_merged.fastq" -t $NCPUS
mv auto_data_*/contigs.fa ${seq}_velvet.scaff.fasta

##########################################################################
###                                                                    ###
###                            GAPFILLER                               ###
###                                                                    ###
##########################################################################

echo "Running gapfiller"
echo "Creating library file"
echo -e "${seq}_Gapfiller\tbwa\t${seq}_1.fastq\t${seq}_2.fastq\t500\t0.25\tFR"
echo -e "${seq}_Gapfiller\tbwa\t${seq}_1.fastq\t${seq}_2.fastq\t500\t0.25\tFR" > Gapfiller.txt

echo "command=perl ${GAPFILLER} -l Gapfiller.txt -s ${seq}_velvet.scaff.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b Velv_scaff"

perl ${GAPFILLER} -l Gapfiller.txt -s ${seq}_velvet.scaff.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b Velv_scaff

if [ -s Velv_scaff/Velv_scaff.gapfilled.final.fa ]; then
  mv Velv_scaff/Velv_scaff.gapfilled.final.fa ${seq}_velvet.fasta
else
  mv ${seq}_velvet.scaff.fasta ${seq}_velvet.fasta
fi

rm -rf ./Velv_scaff/
#TO DO need filecheck here for gapfiller

##########################################################################
###                                                                    ###
###                             ABACAS                                 ###
###                                                                    ###
##########################################################################

if [ "$contig_count" != 1 -a "$ref" != "none" ]; then
  echo "command=perl $ABACAS -m -b -r ref.ABACAS -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped"
  perl ${ABACAS} -m -b -r ref.ABACAS -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped
  echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
  cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
fi

if [ "$contig_count" == 1 -a "$ref" != "none" ]; then
  echo "command=perl $ABACAS -m -b -r ref.ABACAS -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped"
  perl ${ABACAS} -m -b -r ref.ABACAS -q ${seq}_velvet.fasta -p nucmer -o ${seq}mapped
  echo -e "Velvet assembly has been mapped against the reference using ABACAS\n\n"
  cat ${seq}mapped.fasta ${seq}mapped.contigsInbin.fas > ${seq}mapnunmap.fasta
fi

if [ "$ref" == "none" ]; then
  mv ${seq}_velvet.fasta ${seq}mapnunmap.fasta
fi


##########################################################################
###                                                                    ###
###                             IMAGE                                  ###
###                                                                    ###
##########################################################################

## include test for PAGIT assembly
if [ ! -s ${seq}_IMAGE2_out.fasta -a ! -s ${seq}_final.fasta -a ! -s${seq}_icorn.fasta ]; then
  echo "command=perl $IMAGE/image.pl -scaffolds ${seq}mapnunmap.fasta -prefix ${seq} -iteration 1 -all_iteration 3 -dir_prefix ite -kmer 81"
   perl $IMAGE/image.pl -scaffolds ${seq}mapnunmap.fasta -prefix ${seq} -iteration 1 -all_iteration 3 -dir_prefix ite -kmer 81
   echo "command=perl $IMAGE/restartIMAGE.pl ite3 71 3 partitioned"
   perl $IMAGE/restartIMAGE.pl ite3 71 3 partitioned
   echo "command=perl $IMAGE/restartIMAGE.pl ite6 61 3 partitioned"
   perl $IMAGE/restartIMAGE.pl ite6 61 3 partitioned
   echo "command=perl $IMAGE/restartIMAGE.pl ite9 51 3 partitioned"
   perl $IMAGE/restartIMAGE.pl ite9 51 3 partitioned
   echo "command=perl $IMAGE/restartIMAGE.pl ite12 41 3 partitioned"
   perl $IMAGE/restartIMAGE.pl ite12 41 3 partitioned
   echo "command=perl $IMAGE/restartIMAGE.pl ite15 31 3 partitioned"
   perl $IMAGE/restartIMAGE.pl ite15 31 3 partitioned
   echo "command=perl $IMAGE/restartIMAGE.pl ite18 21 3 partitioned"
   perl $IMAGE/restartIMAGE.pl ite18 21 3 partitioned
   
   mv ite21/new.fa ${seq}_IMAGE2_out.fasta

  perl $IMAGE/image_run_summary.pl ite > IMAGE2.summary

 rm -r ite*
 rm partitioned_1.fastq
 rm partitioned_2.fastq
fi  

##########################################################################
###                                                                    ###
###                             SSPACE                                 ###
###                                                                    ###
##########################################################################

  #TODO need to write SSPACE version check for different library file. The below is for SSPACEv3.0
  
  #echo -e "${seq}SSPACE\tbowtie\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" > ${seq}/library.txt
  
  #For SSPACE v2.0 basic
  echo -e "${seq}SSPACE\t${seq}_1.fastq\t${seq}_2.fastq\t200\t0.25\tFR" > library.txt
  
  perl $SSPACE -l library.txt -s ${seq}_IMAGE2_out.fasta
  mv standard_output.final.scaffolds.fasta ${seq}SSPACE.fasta
  rm -r pairinfo
  rm -r intermediate_results
  rm -r bowtieoutput
  rm -r reads
  rm standard_output.final.evidence
  rm standard_output.logfile.txt
  #mv ${seq}/standard_output.summary.txt ${seq}/SSPACE.summary.txt ##standard_output.summaryfile.txt is the correct name for this output


### SSPACE test ############
## This will skip the next step if SSPACE doesn't find anything to scaffold in your assembly, which caused a gapfiller crash



if [ -s standard_output.summaryfile.txt ]; then
  SSPACE_test=`grep 'Total number of N' standard_output.summaryfile.txt |tail -n1 |awk '{print $6}'`
  if [ "$SSPACE_test" == 0 ]; then
   cp ${seq}SSPACE.fasta ${seq}_gap2.fasta
  fi
fi


##########################################################################
###                                                                    ###
###                            GAPFILLER 2                             ###
###  This step is skipped is SSPACE doesn't find anything to scaffold  ###
###                                                                    ###
##########################################################################
if [ ! -s ${seq}_gap2.fasta ]; then
  perl $GAPFILLER -l Gapfiller.txt -s ${seq}SSPACE.fasta -m 20 -o 2 -r 0.7 -n 10 -d 50 -t 10 -T ${NCPUS} -i 3 -b SSPACE_scaff
    mv SSPACE_scaff/SSPACE_scaff.gapfilled.final.fa ${seq}_gap2.fasta
	rm -r SSPACE_scaff/
fi

##########################################################################
###                                                                    ###
###                Remove contigs <1kb and image cleanup               ###
###                                                                    ###
##########################################################################
if [ "$long" == "no" ]; then
  $CONVERT_PROJECT -f fasta -t fasta -x 1000 -R Contig ${seq}_gap2.fasta ${seq}_pilon
  echo -e "Project has been filtered to remove contigs less than 1kb in size \n" 
else 
  mv ${seq}_gap2.fasta ${seq}_pilon.fasta
  echo -e "Project includes all contigs including <1kb in size\n" 
fi

##########################################################################
###                                                                    ###
###                                PILON                               ###
###                                                                    ###
##########################################################################

#create bam file before running pilon
if [ ! -s ${seq}_pilon.fasta.bwt ]; then
  bwa index ${seq}_pilon.fasta
  else
  echo "Found ref index for Pilon"
fi
if [ ! -s ${seq}.sam ]; then
    bwa mem -R '@RG\tID:Assembly\tSM:${seq}\tPL:ILLUMINA' -a -t $NCPUS ${seq}_pilon.fasta ${seq}_1.fastq ${seq}_2.fastq > ${seq}.sam
  else
    echo "Found bam file for pilon"  
fi
if [ ! -s ${seq}.bam ]; then
  	    samtools view -h -b -@ 1 -q 1 -o ${seq}.bam.tmp ${seq}.sam && samtools sort -@ 1 -o ${seq}.bam ${seq}.bam.tmp
		rm ${seq}.bam.tmp ${seq}.sam
fi
if [ ! -s ${seq}.bam.bai ]; then
    samtools index ${seq}.bam
fi

if [ ! -s pilon.fasta ]; then
   java -jar ${PILON} --genome ${seq}_pilon.fasta --frags ${seq}.bam
   mv pilon.fasta ${seq}_final.fasta
fi 

exit 0

