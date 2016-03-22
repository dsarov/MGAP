#!/bin/bash

### this script will change into the correct directory, and run the GATK SNP calling

dir=$1

cd $dir

### call the SNP caller
java -jar /nfs/users/nfs_m/mh12/bin/GenomeAnalysisTK-2.0-35-g2d70733/GenomeAnalysisTK.jar -T RealignerTargetCreator -I out.sorted.markdup.bam -o tdo.gatk_variants.realign.intervals -R ref.fa -nt 8; 
java -jar /nfs/users/nfs_m/mh12/bin/GenomeAnalysisTK-2.0-35-g2d70733/GenomeAnalysisTK.jar -T IndelRealigner  -I out.sorted.markdup.bam -targetIntervals  tdo.gatk_variants.realign.intervals   -o gatk.realigned.bam -R ref.fa;
java -jar /nfs/users/nfs_m/mh12/bin/GenomeAnalysisTK-2.0-35-g2d70733/GenomeAnalysisTK.jar -T UnifiedGenotyper -I gatk.realigned.bam -pnrm POOL -o gatk_variants.variants.vcf -ploidy 1 -glm POOLBOTH -R ref.fa  -nt 8;

cd ..
