#!/bin/bash

### little script to align with maq a read pair Illumina/Solexa lane onto a genome to produce a bam file for Artemis.
# samtools and bwa must be in the PATH

genome=$1
kmer=$2
stepsize=$3
fastqF=$4
fastqR=$5
resultname=$6
insertSize=$7

#randome stamp for temporary file name
tmp=$$
if [ -z "$SMALT_PARAMETER" ] ; then 
SMALT_PARAMETER=" -y 0.5 -r 0 -x "; export SMALT_PARAMETER
fi

if [ -z "$insertSize" ] ;
		then
		echo "Usage: <genome> <k-mer> <stepsize> <Illumina forward reads> <Illumina reverse reads; 0 if single Reads> <resultname> <Insertsize Range:i.e. 350, 0 if single> optional: for more specific parameter, set the SMALT_PARAMETER enrovimental variable"
		echo
		exit;
	fi;

if [ -z "$ICORN2_THREADS" ]; then
  ICORN2_THREADS=1;
fi

export ICORN2_THREADS

$ICORN2_HOME/smalt index -k $kmer -s $stepsize Genome.index.$tmp $genome > out.$tmp.txt

echo "SMALT mapping parameters are set of following (SMALT_PARAMETER): $SMALT_PARAMETER";
if [ "$fastqR" = "0" ] ; then
	$ICORN2_HOME/smalt map $SMALT_PARAMETER -n 3  -f samsoft  -o $resultname.$tmp.sam Genome.index.$tmp $fastqF >> out.$tmp.txt
else 
	
# Check if the insert size is in the correct format:
# Generating sam with SSAHA2 (>= version 2.5, must be in path)
	if [ ! -z "$ICORN2_VERBOSE" ] ; then 
		echo "	smalt map $SMALT_PARAMETER  -f samsoft -i $insertSize -o $resultname.$tmp.sam Genome.index.$tmp $fastqF $fastqR >> out.$resultname.$tmp.txt"
	fi
	
	$ICORN2_HOME/smalt map $SMALT_PARAMETER -n $ICORN2_THREADS  -f samsoft -i $insertSize -o $resultname.$tmp.sam Genome.index.$tmp $fastqF $fastqR >> out.$resultname.$tmp.txt
	smaltOK=$?
	if [ $smaltOK -ne 0 ] ; then
		echo "Problems with SMALT: $smaltOK"
		exit $smaltOK
	fi
	
fi
rm  Genome.index.$tmp.*
#if [ -z "$tmp.sam" ] ; then
#	echo "Problem generating Sam file.";
#	exit;
#fi

# Convert reference to fai for bam file
samtools faidx $genome


# SAM to BAM

if [ ! -z "$ICORN2_VERBOSE" ] ; then 
	echo "Doing some samtools now, until I have a ok looking bam file"
fi
read=$(echo $fastqF | sed 's/_1.fastq//g');
echo -e "@RG\tID:$read\tSM:1" > line.$tmp.txt
awk -va=$read '{if ($1~/^@/) {print} else  {print $0"\tRG:Z:"a}}'  $resultname.$tmp.sam | samtools view -b -t $genome.fai - > $resultname.$tmp.bam
rm $resultname.$tmp.sam

samtools view -H  $resultname.$tmp.bam > header.$tmp.txt  2> output.samtools.txt
echo -e "@RG\tID:$read\tSM:1" >>  header.$tmp.txt
samtools reheader header.$tmp.txt  $resultname.$tmp.bam  > $resultname.$tmp.2.bam 2> output.samtools.txt
rm $resultname.$tmp.bam

# <(cat line.$tmp.txt $resultname.$tmp.sam)
mkdir $resultname
#samtools import $genome.fai $resultname.$tmp.sam $resultname.$tmp.bam &> output.samtools.txt
#samtools view -S -b -h $resultname.$tmp.sam > $resultname.test.$tmp.bam
#order the bam file
samtools sort $resultname.$tmp.2.bam $resultname/out 2> output.samtools.txt
rm  $resultname.$tmp.2.bam
if [ ! -z "$ICORN2_VERBOSE" ] ; then  
	echo "Marking duplicates, need some memory here, possible reason for crashes";
	echo "	java -XX:MaxPermSize=512m -Xmx2000m -XX:ParallelGCThreads=1 -XX:-UseParallelGC -XX:+UseSerialGC -jar $ICORN2_HOME/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=$resultname/out.bam O=$resultname/out.sorted.markdup.bam"
	
fi 
java -XX:MaxPermSize=512m -Xmx2000m -XX:ParallelGCThreads=1 -XX:-UseParallelGC -XX:+UseSerialGC -jar $ICORN2_HOME/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=$resultname/out.bam O=$resultname/out.sorted.markdup.bam &> out.markduplicates.txt
return=$?
if [ "$return" != "0" ] ;
	then 
	echo "Sorry, mark duplicates failed again..."
	echo "java -XX:MaxPermSize=512m -Xmx2000m -XX:ParallelGCThreads=1 -XX:-UseParallelGC -XX:+UseSerialGC -jar $ICORN2_HOME/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=$resultname/out.bam O=$resultname/out.sorted.markdup.bam"
	exit 1;
fi;
#index the bam file, to get the bai file.
rm $resultname/out.bam

samtools index  $resultname/out.sorted.markdup.bam

###delete files
rm  line.$tmp.txt out.$tmp.txt  header.$tmp.txt

