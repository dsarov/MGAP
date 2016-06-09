# MGAP

<i>What is MGAP?</i>

MGAP (Microbial Genome Assembler Pipeline) is an automated pipeline that incorporates several genome assembly programs into one tool for high-quality, reference-assisted assembly of microbial genomes using paired-end Illumina data.

<i>What resource managers can MGAP use?</i>

MGAP currently only works with the PBS (Torque) resource manager.

<i>MGAP workflow</i>

MGAP performs reference assisted assemblies on paired-end Illumina FASTQ reads, either in phred +33 or phred +64 quality score format. MGAP cannot be used on single-end data, or for NGS data generated on platforms other than Illumina.

To achieve high-quality assemblies, MGAP incorporates the following programs into its workflow:
- <b>Trimmomatic</b> [(Bolger et al., 2014)](http://bioinformatics.oxfordjournals.org/content/30/15/2114)
- <b>Velvet</b> [(Zerbino & Birney, 2008)](http://genome.cshlp.org/content/18/5/821.full)
- <b>VelvetOptimiser</b> (https://github.com/tseemann/VelvetOptimiser)
- <b>GapFiller</b> [(Boetzer & Pirovano, 2012)](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-6-r56)
- <b>ABACAS</b> [(Assefa et al., 2009)](http://bioinformatics.oxfordjournals.org/content/25/15/1968.long)
- <b>IMAGE</b> [(Tsai et al., 2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-4-r41)
- <b>SSPACE</b> [(Boetzer et al., 2011)](http://bioinformatics.oxfordjournals.org/content/27/4/578.long)
- <b>ICORN2</b> [(Otto et al., 2010)](http://bioinformatics.oxfordjournals.org/content/26/14/1704.full)
- <b>MIRA (convert_project)</b> (Chevreux et al., 1999)

<i>How do I run MGAP?</i>

USAGE: MGAP.sh -r [reference, without .fasta extension] -s [specify single strain]

If you would prefer MGAP to perform assemblies without using a reference to assist or if no reference is available, set the -r flag to "none".

MGAP expects reads to be paired-end Illumina data in the following format: STRAIN_1_sequence.fastq.gz (first pair) and STRAIN_2_sequence.fastq.gz (second pair).

<i>What does MGAP do?</i>

Prior to assembly, FASTQ reads are conservatively trimmed and filtered using Trimmomatic v0.35 to remove low-quality bases and Illumina adapter contamination. The Trimmomatic parameters in MGAP are: LEADING=3, TRAILING=3, SLIDINGWINDOW=4:15, MINLEN=36, and ILLUMINACLIP (for TruSeq2 paired-end adapters). These parameters can be altered as desired. 

Next, a draft scaffold assembly is  created using Velvet 1.2.10, with parameters optimised using VelvetOptimiser v2.2.4 at a default kmer range of 53 to 75; these kmers can also be altered if required. 

GapFiller v2.1.1 is then used to attempt to fill in the scaffolds created by Velvet. Following creation of the draft Velvet assemblies, the best assembly is then improved upon using ABACAS v1.3.1 and IMAGE v2.4.1. If the user provides a reference genome, ABACAS scaffolds the contigs against this reference. If no user-specified reference is provided, this ABACAS step is skipped. Scaffolded contigs are then attempted to be joined using IMAGE, which will break contigs that have been incorrectly joined. A second attempt is then made to scaffold contigs using SSPACE v3.0. GapFiller is then run again to attempt to fill in the scaffolds created by SSPACE. 

Finally, ICORN2 v0.95 corrects any insertion-deletion (indel) and SNP errors in the final assembly. Contigs <1,000bp are excluded from the final .fasta assembly output using the miraconvert tool in MIRA v4. 

<i>Who created MGAP?</i>

MGAP was written by Derek Sarovich ([@DerekSarovich](https://twitter.com/DerekSarovich)) and Erin Price ([@Dr_ErinPrice](https://twitter.com/Dr_ErinPrice)) at Menzies School of Health Research, Darwin, NT, Australia.

<i>What to do if I run into issues with MGAP?</i>

Please send bug reports to mshr.bioinformatics@gmail.com or derek.sarovich@menzies.edu.au.

We are in the process of uploading our assembly pipeline to GitHub. If you are interested in our project, please contact us.

Please check back shortly for our updated project!
