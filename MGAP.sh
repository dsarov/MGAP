#!/bin/bash



Thanks for using MGAP!!



usage()
{
echo -e  "USAGE: MGAP.sh -r <reference, without .fasta extension> -s <specify single strain>"
}
help()
{
cat << _EOF_
Thanks for using MGAP
The microbial genome assembler is an automated assmbly pipeline for paired end Illumina data
The MGAP pipeline workflow is as follows:
FastQ data is first conservatively filtered with Trimmomatic, next a draft assembly is constructed using Velvet with parameters optimised using velvet optimiser.
This draft assembly is then scaffolded and draft contigs are stiched together if possible using ABACAS and IMAGE.


_EOF_
}


#TODO
#include help


help()
{


}

#Define path to SPANDx install
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if  [ $# -lt 1 ]
    then
	    usage
		exit 1
fi

# source dependencies
source "$SCRIPTPATH"/SPANDx.config 
source "$SCRIPTPATH"/qsub.config

#declare variables
declare -rx SCRIPT=${0##*/}

echo "scriptpath = $SCRIPTPATH"

OPTSTRING="h:r:s:"

declare SWITCH
ref=no
org=haploid
seq_directory="$PWD"
assemble=yes

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do 
		case $SWITCH in
		
		r) ref="$OPTARG"
		   echo "Reference = $ref"
		   ;;
         
		s) strain="$OPTARG"
		   echo "Only strain $strain will be processed. By default all strains within current directory will be processed"
		   ;;
			
		\?) usage
		    exit 1
		    ;;
		
		h) usage
		   exit 1
		   ;;
		   
		*) echo "script error: unhandled argument"
           exit 1
		   usage
		   ;;
		   
		
	esac
done

 
echo -e "\nThe following parameters will be used\n"
echo -e "-------------------------------------\n"
echo "Output directory that will contain assemblies = $seq_directory/Assemblies"

ref_index=${seq_directory}/${ref}.bwt #index file created with BWA
REF_INDEX_FILE=${seq_directory}/${ref}.fasta.fai #index created with SAMtools
REF_BED=${seq_directory}/${ref}.bed
REF_DICT=${seq_directory}/${ref}.dict #Dictionary file created with Picard tools

if [ ! $PBS_O_WORKDIR ]; then
        PBS_O_WORKDIR="$seq_directory"
fi

cd $PBS_O_WORKDIR

## file checks and program checks

if [ ! -s "$ref.fasta" ]; then
        echo -e "Couldn't find reference file. Please make sure that reference file is in the analysis directory\n"
        exit 1
    else
	    echo -e "Found FASTA file\n"
fi

ref_blank_lines=`grep -c '^$' $ref.fasta`

if [ "$ref_blank_lines" != 0 ]; then
	    echo -e "ERROR: Reference FASTA file is formatted incorrectly and must contain 0 blank lines. Blank lines will cause BWA and GATK to fail."
        exit 1
    else
	    echo -e "FASTA file looks to contain no blank lines. Good.\n"
fi
#TODO

# include tests for assembly module #
## Test for dependencies required by MGAP

SHUFFLE="/usr/local/velvet_1.2.10/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl";

##TO do 
# needs to be included as an optional multiple kmer scenario 
ASSASS="/usr/local/velvet_1.2.10/contrib/AssemblyAssembler1.3/AssemblyAssembler1.3.py";

VELVET_PATH="/usr/local/velvet_1.2.10";
CONVERT_PROJECT="/home/dsarovich/bin/convert_project";
ABACAS="/home/dsarovich/bin/PAGIT/ABACAS/abacas.pl";
VELVETH="/usr/local/velvet_1.2.10/velveth";
VELVETG="/usr/local/velvet_1.2.10/velvetg";
TRIM="/usr/local/Trimmomatic-0.35/trimmomatic-0.35.jar";
VELVETOPT="/home/dsarovich/bin/velvet_1.2.10/contrib/VelvetOptimiser-2.2.4/VelvetOptimiser.pl"

IMAGE="/home/dsarovich/bin/IMAGE_version2.4/";
ICORN2_HOME="/home/dsarovich/bin/ICORN2/";
export ICORN2_HOME
SSPACE="/home/dsarovich/bin/SSPACE-BASIC-2.0_linux-x86_64/SSPACE_Basic_v2.0.pl";

PAGIT_HOME="/home/dsarovich/bin/PAGIT";


bwa_test=`command -v "$BWA"`

java_test=`command -v "$JAVA"`

if [ -z "$bwa_test" ]; then
	    echo "ERROR: SPANDx requires BWA to function. Please make sure the correct path is specified in SPANDx.config or the executable is in your PATH"
		exit 1
fi

if [ -z "$java_test" ]; then
	    echo "ERROR: MGAP requires java. Please make sure java is available on your system. The PATH to java can be modified in the MGAP.config file"
		exit 1
fi


### Handling and checks for read files
if [ "$strain" == all ]; then
    sequences_tmp=(`find $PBS_O_WORKDIR/*_1_sequence.fastq.gz -printf "%f "`)
    sequences=("${sequences_tmp[@]/_1_sequence.fastq.gz/}")
    n=${#sequences[@]}
    if [ $n == 0 ]; then
        echo -e "Program couldn't find any sequence files to process"
        echo -e "Sequences must be in the format STRAIN_1_sequence.fastq.gz STRAIN_2_sequence.fastq.gz for paired end reads"
    	echo -e "and STRAIN_1_sequence.fastq.gz for single end data\n"
	    exit 1
    else
        echo -e "Sequences have been loaded into SPANDx\n"
    fi
fi

## check for read pairing and correct notation #need to test
if [ "$pairing" == PE -a "$strain" == all ]; then
	sequences_tmp2=(`find $PBS_O_WORKDIR/*_2_sequence.fastq.gz -printf "%f "`)
    sequences2=("${sequences_tmp2[@]/_2_sequence.fastq.gz/}")
    n2=${#sequences2[@]}
    if [ $n != $n2 ]; then
	    echo "Number of forward reads don't match number of reverse reads. Please check that for running in PE mode all read files have correctly named pairs"
		exit 1
	fi
	for (( i=0; i<n; i++ )); do
	    if [ ${sequences[$i]} != ${sequences2[$i]} ]; then
            echo "Names of forward reads don't match names of reverse reads. Please check that for running in PE mode all read files have correctly named pairs"
			echo "Forward read names are ${sequences[@]}"
			echo "Reverse read names are ${sequences2[@]}"
            exit 1
        fi
    done;
fi

if [ "$strain" != all -a "$strain" != none ]; then
    sequences="$strain"
	echo -e "SPANDx will process the single strain $strain.\n"
	if [ "$matrix" == yes ]; then
	    matrix=no
		echo -e "SPANDx will not run SNP matrix with only one strain specified. If a SNP matrix is desired please run SPANDx with more strains or if strains have already been run link the SNP VCF and BAM files to the appropriate directories and set -s to none.\n"
	fi
	if [ "$pairing" == PE ]; then
	    if [ ! -s ${strain}_1_sequence.fastq.gz ]; then
		    echo -e "ERROR: SPANDx cannot find sequence file ${strain}_1_sequence.fastq.gz\n"
			echo -e "Please make sure the correct sequence reads are in the sequence directory\n"
			exit 1
        fi
		if [ ! -s ${strain}_2_sequence.fastq.gz ]; then
		    echo -e "ERROR: SPANDx cannot find sequence file ${strain}_2_sequence.fastq.gz\n"
			echo -e "Please make sure the correct sequence reads are in the sequence directory\n"
			exit 1
        fi
	fi
	if [ "$pairing" == SE ]; then
		if [ ! -s ${strain}_1_sequence.fastq.gz ]; then
		    echo -e "ERROR: SPANDx cannot find sequence file ${strain}_1_sequence.fastq.gz\n"
			echo -e "Please make sure the correct sequence reads are in the sequence directory\n"
			exit 1
        fi
	fi
fi	
	
	
if [ ! -d "tmp" ]; then
	mkdir $seq_directory/tmp
fi
if [ ! -d "Assemblies" ]; then
	mkdir $seq_directory/Assemblies
fi

if [ -s qsub_ids.txt ]; then
    rm qsub_ids.txt
fi


## Indexing of the reference with SAMTools and BWA
skip () {
if [ ! -s "$ref_index" ]; then
	echo -e "Submitting qsub job for BWA reference indexing\n"
    cmd="$BWA index -a is -p ${seq_directory}/${ref} ${seq_directory}/${ref}.fasta"
    qsub_id=`qsub -N index -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "index\t$qsub_id" >> qsub_ids.txt
fi
if [ ! -s "$REF_INDEX_FILE" ]; then
    echo -e "Submitting qsub job for SAMtools reference indexing\n"
    cmd="$SAMTOOLS faidx ${seq_directory}/${ref}.fasta"
    qsub_id=`qsub -N SAM_index -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "SAM_index\t$qsub_id" >> qsub_ids.txt
fi
if [ ! -s "$REF_DICT" ]; then
    echo -e "Submitting qsub job for ${ref}.dict creation\n"
    cmd="$JAVA $SET_VAR $CREATEDICT R=${seq_directory}/${ref}.fasta O=$REF_DICT"
	qsub_id=`qsub -N PICARD_dict -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
	echo -e "PICARD_dict\t$qsub_id" >> qsub_ids.txt
fi
if [ ! -s "$REF_BED" -a qsub_ids.txt ]; then
    echo -e "Submitting qsub job for BED file construction with BEDTools\n"
	qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
	depend="-W depend=afterok${qsub_cat_ids}"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    qsub_id=`qsub -N BED_window -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "BED_window\t$qsub_id" >> qsub_ids.txt	
fi
if [ ! -s "$REF_BED" -a ! qsub_ids.txt ]; then
    echo -e "Submitting qsub job for BED file construction with BEDTools\n"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    qsub_id=`qsub -N BED_window -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "BED_window\t$qsub_id" >> qsub_ids.txt	
fi
}

##TODO
#include resource management for SGE and SLURM

assemble_single ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
    if [ ! -s ${PBS_O_WORKDIR}/Assemblies/${sequences}_final.fasta ]; then
		echo -e "Submitting qsub job for de novo assembly of $sequences\n"
        var="seq=$sequences,ref=$ref,seq_path=$seq_directory,kmer=$kmer,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N Assemble_sequences -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=$NCPUS,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/assembler.sh`
        echo -e "aln_sequences\t$qsub_array_id" >> qsub_array_ids.txt
	fi
        
fi
if [ ! -s qsub_ids.txt ]; then
    if [ ! -s ${PBS_O_WORKDIR}/Assemblies/${sequences}_final.fasta ]; then
		echo -e "Submitting qsub job for de novo assembly of ${sequences[$i]}\n"
	    var="seq=$sequences,ref=$ref,seq_path=$seq_directory,kmer=$kmer,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N Assemble_$sequences -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=$NCPUS,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/assembler.sh`
		echo -e "aln_$sequences\t$qsub_array_id" >> qsub_array_ids.txt
	fi
fi

}

assemble ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Assemblies/${sequences[$i]}_final.fasta ]; then
		        echo -e "Submitting qsub job for de novo assembly of ${sequences[$i]}\n"
                var="seq=${sequences[$i]},ref=$ref,seq_path=$seq_directory,kmer=$kmer,SCRIPTPATH=$SCRIPTPATH"
		        qsub_array_id=`qsub -N Assemble_${sequences[$i]} -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=$NCPUS,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/assembler.sh`
                echo -e "aln_${sequences[$i]}\t$qsub_array_id" >> qsub_assemble_ids.txt
				sleep 0.25
	        fi
        done
fi
if [ ! -s qsub_ids.txt ]; then
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Assemblies/${sequences[$i]}_final.fasta ]; then
		        echo -e "Submitting qsub job for de novo assembly of ${sequences[$i]}\n"
	    	    var="seq=${sequences[$i]},ref=$ref,seq_path=$seq_directory,kmer=$kmer,SCRIPTPATH=$SCRIPTPATH"
		        qsub_array_id=`qsub -N Assemble_${sequences[$i]} -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=$NCPUS,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/assembler.sh`
				echo -e "aln_${sequences[$i]}\t$qsub_array_id" >> qsub_assemble_ids.txt
				sleep 0.25
	        fi
        done
fi
	}
	
if [ "$assemble" == yes -a "$strain" == all ]; then
    assemble
fi
if [ "$assemble" == yes -a "$strain" != all ]; then
    assemble_single
fi

exit 0