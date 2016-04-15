# Microbial Genome Assembler Pipeline (MGAP)

MGAP is an automated pipeline that incorporates several existing genome assembly tools for high-quality, reference-assisted assembly of  microbial genomes using paired-end Illumina data.

<i>MGAP workflow</i>

MGAP performs assemblies on paired-end Illumina FASTQ reads, either in phred +33 or phred +64 quality score format. MGAP cannot be used for single-end data, or for platforms other than Illumina.

To achieve high-quality assemblies, MGAP incorporates several excellent and freely available programs into its workflow:
- Trimmomatic
- Velvet
- VelvetOptimiser
- Gapfiller
- ABACAS
- IMAGE
- SSPACE
- ICORN2
- MIRA

Prior to assembly, FASTQ reads are conservatively trimmed and filtered using Trimmomatic v0.35 (Bolger AM et al., 2014) to remove low-quality bases and Illumina adapter contamination. The Trimmomatic parameters in MGAP are: LEADING=3, TRAILING=3, SLIDINGWINDOW=4:15, MINLEN=36, and ILLUMINACLIP (for TruSeq2 paired-end adapters). These parameters can be altered as desired. 

Next, a draft scaffold assembly is  created using Velvet 1.2.10 (Zerbino DR & Birney E, 2008), with parameters optimised using VelvetOptimiser v2.2.4 (https://github.com/tseemann/VelvetOptimiser) at a default kmer range of 53 to 75; these kmers can also be altered if required. 

GapFiller v2.1.1 (Boetzer M & Pirovano W, 2012) is then used to attempt to fill in the scaffolds created by Velvet. Following creation of the draft Velvet assemblies, the best assembly is then improved upon using ABACAS v1.3.1 (Assefa S et al., 2009) and IMAGE v2.4.1 (Tsai IJ et al., 2010). If the user provides a reference genome, ABACAS scaffolds the contigs against this reference. If no user-specified reference is provided, this ABACAS step is skipped. Scaffolded contigs are then attempted to be joined using IMAGE, which will break contigs that have been incorrectly joined. A second attempt is then made to scaffold contigs using SSPACE v3.0 (Boetzer M et al., 2011). GapFiller is then run again to attempt to fill in the scaffolds created by SSPACE. 

Finally, ICORN2 v0.95 (Otto TD et al., 2010) corrects any insertion-deletion (indel) and SNP errors in the final assembly. Contigs <1,000bp are excluded from the final .fasta assembly output using the miraconvert tool in MIRA4 (Chevreux B et al., 1999). 

USAGE: MGAP.sh -r <reference, without .fasta extension> -s <specify single strain>

MGAP expects reads to be paired-end Illumina data in the following format: STRAIN_1_sequence.fastq.gz (first pair) and STRAIN_2_sequence.fastq.gz (second pair). 

We are in the process of uploading our assembly pipeline to GitHub. If you are interested in our project, please contact us at mshr.bioinformatics@gmail.com or derek.sarovich@menzies.edu.au.

Please check back shortly for our updated project!
