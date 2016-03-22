#!/bin/bash

### will call the correction bit of icorn2

dir=$1
readRoot=$2
fragmentSize=$3
resultname=$4

## first set the variables
ICORN_HOME=/nfs/pathogen003/tdo/Tools/ICORN.SNP; export ICORN_HOME    
SNPOMATIC_HOME=$ICORN2_HOME; export SNPOMATIC_HOME    
PERL5LIB=$ICORN2_HOME:$PERL5LIB; export PERL5LIB 

cd $dir

echo perl $ICORN2_HOME/icorn2.Correct.pl ref.fa gatk_variants.variants.vcf ../$readRoot $fragmentSize $resultname	
perl $ICORN2_HOME/icorn2.Correct.pl ref.fa gatk_variants.variants.vcf ../$readRoot $fragmentSize $resultname
	
tmp=$$
### get it into correct sequence length
perl $ICORN2_HOME/fasta2singleLine.pl $resultname  $resultname.$tmp
perl $ICORN2_HOME/fasta2multipleLine.pl $resultname.$tmp $resultname 80
#/nfs/users/nfs_t/tdo/Bin/bam.correctLineLength.sh $resultname

echo cp $resultname ../

cp $resultname ../
cd ..
