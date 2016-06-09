#!/bin/bash

readRoot=$1
fragmentSize=$2
referenceORG=$3
start=$4
end=$5


User=$(whoami)


if [ ! -z "$ICORN2_HOME" ] ; then

#ICORN2_HOME=/nfs/users/nfs_t/tdo/Bin; export ICORN2_HOME
	echo "Please set the ICORN2_HOME directory. It is the directory where this program is copied to."
ABSOLUTE_PATH2=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)/`basename "${BASH_SOURCE[0]}"`
ABSOLUTE_PATH=$(echo $ABSOLUTE_PATH2 | sed 's/icorn2.serial.sh//g');
echo "As this was not done, we are setting it automatically: ICORN2_HOME=$ABSOLUTE_PATH; export ICORN2_HOME"
ICORN2_HOME=$ABSOLUTE_PATH; export ICORN2_HOME
echo "waiting for 5 seconds"
sleep 5;
echo $ICORN2_HOME
fi

if [ -z "$end" ] ; then
	echo "Oh $User, you have to call me:"
	echo "<readroot> <fragementSize> <reference> <iterationstart> <iterationStop>"
	echo 
	echo "iCORN2 expected read pairs! Name convention is readroot_1.fastq and readroot_2.fastq are the two fastq file names."
	echo "To continue a correct, after three iteration, just put the iterationstart to 4."
	echo
	echo "Extras:"
	echo "If you have a multithread computer, please set following variable ICORN2_THREADS, like ICORN2_THREADS=8; export ICORN2_THREADS. Multithreading will be used in the mapping stage in the GATK variant calls."
	echo "If you like to change the SMALT mapping parameter, change the variable SMALT_PARAMETER"
	echo "For debugging information, set the BASH variable ICORN2_VERBOSE to an interger"
	exit 1;
	
fi
 
read1="$readRoot"_1.fastq
read2="$readRoot"_2.fastq

if [ ! -f "$referenceORG" ] ; then
	echo "$User, where is the correct reference?"
	exit 1; 
fi

if [ -z "$ICORN2_THREADS" ]; then
  ICORN2_THREADS=1;
fi

export ICORN2_THREADS

### getting a local copy of the reference, without | in the name
refRoot=$(echo $referenceORG | perl -nle '@ar=split(/\//); print $ar[($ar[scalar(@ar)]-1)]');

if [ ! -f "ICORN2.$refRoot.$start" ] ; then
cat $referenceORG | perl -e 'while(<STDIN>) {if (/>(\S+)/){ chomp; s/\|/_/g; s/\./_/g; print "$_\n";} else {print}}' > ICORN2.$refRoot.$start
fi

reference=ICORN2.$refRoot

tmp=$$

echo $reference 
echo
for ((i=$start;$i<=$end;i++)); do
	directory=ICORN2_$i

	### call the mapper
	echo "Iteration ++++ $i";

	$ICORN2_HOME/icorn2.smalt.sh  ICORN2.$refRoot.$i 13 3  $read1 $read2 ICORN2_$i 1200
	return=$?
	if [ "$return" != "0" ] ;
		then 
		echo "See Error in smalt script: 	$ICORN2_HOME/icorn2.smalt.sh  ICORN2.$refRoot.$i 13 3  $read1 $read2 ICORN2_$i 1200";
		exit 1;
	fi;	
	### call the SNP caller
    cd ICORN2_$i
	ln -s ../$reference.$i ref.fa
### call the SNP caller
	if [ ! -z "$ICORN2_VERBOSE" ] ; then
	echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam        -o gatk_variants.realign.intervals -R  ref.fa -nt $ICORN2_THREADS"
	echo "	java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner         -I out.sorted.markdup.bam        -o gatk.realigned.bam -targetIntervals gatk_variants.realign.intervals   -R  ref.fa"
	echo "java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper       -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa -nt $ICORN2_THREADS" 
	fi
	java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam        -o gatk_variants.realign.intervals -R  ref.fa -nt $ICORN2_THREADS &> out.gatk.1.txt
	java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T IndelRealigner         -I out.sorted.markdup.bam        -o gatk.realigned.bam -targetIntervals gatk_variants.realign.intervals   -R  ref.fa &> out.gatk.2.txt
	java -jar $ICORN2_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper       -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa -nt $ICORN2_THREADS &> out.gatk.3.txt
	
	cd ..

	
	### call icorn starter
	$ICORN2_HOME/icorn2.correct.sh $directory $readRoot $fragmentSize $reference.$(($i+1)) > out.ICORN2.CORRECT.$i.o
	return=$?
	if [ "$return" != "0" ] ;
		then 
		echo "See Error in icron correction script: $ICORN2_HOME/icorn2.correct.sh $directory $readRoot $fragmentSize $reference.$(($i+1))"
		exit 1;
	fi;	
	echo "Current corrections:"
	perl $ICORN2_HOME/icorn2.collectResults.pl .
done
	

echo "To look in into a correction, open the file ending with .1 in artemis and load the gff file onto it..."

