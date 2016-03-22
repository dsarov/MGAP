#####################################
### debug version: please send an email to tdo@sanger.ac.uk if problems occur.
#####################################

Dear icorn user!

Finally we are updating iCORN to the next version. It took us a bit of time to release this version. The reason is that it was mostly an implementation improvement, close to our LSF system. But with the recent PacBio data and the fact that iCORN 1 does not scale very well, we decided to release iCORN2.

New:
-	iCORN2 is based on SMALT (mapper), samtools, GATK, snp-o-matic and some PERL scripts.
-	It is much faster and more robuster
-	It can be run in parallel through threading.
-	It will be easy to implement different mappers or change the GAKT setting
-	Does not compute the gff position to the first version
-	Good for pacbio correction (but that’s not new as icorn1 did this pretty good already


#####################################
Usage:
Unzip the icorn2.xxx.tgz with 
tar xvzf icorn2.v0.95.tgz

You can set the ICORN2_HOME to the location you extracted icorn2, but you don’t have to. Just remember the directory you unzipped it (pwd) and add ICORN2. 

Test the program:
In the ICORN2 directory there is the example directory. Change to it 
cd ICORN2/example

perl ../icorn2.sh PAGIT_test 350 Query.contigs.fa 1 3


#####################################
Application:
We applied iCORN2 to several PacBio assemblies of bacteria. In general it takes several hours to do three iterations. Further we tested iCORN2 on Plasmodium pacBio assembly, around 2 days of run time, and correct many errors, even quiver ran… Never the less, due to the high AT content of the Plasmosdium test case, we could not converge on correction, mostly due to homopolymer tracks.

ICORN2 is regularly run on many project, just to 350mb genomes. We have an LSF version, please contact me, or look into the LSF directory - good luck!

#####################################
Third party program:

We depending on third party programs…
GATK
SMALT
samtools
SNP-o-AMTIC

