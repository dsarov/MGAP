#!/bin/bash

usage()
{
echo -e  "USAGE: MGAP.sh -r <reference, without .fasta extension> -s <specify single strain>\n\n"
}
help()
{
cat << _EOF_

Thanks for using MGAP

Microbial genome assembler is an automated assembly pipeline for paired end Illumina data

The MGAP pipeline workflow is as follows:

Data Filtering:
FastQ data is first conservatively filtered with Trimmomatic to remove adapter contamination and very low quality base calls.

Initial Assembly Draft:
A draft assembly is constructed using Velvet with parameters optimised using velvet optimiser.

Assembly Improvement and Error Correction
Scaffolds created with Velvet are first attempted to be filled using Gapfiller.
Next the draft contigs are ordered against a reference genome using ABACAS.
The scaffolded and ordered contigs are stiched together if possible using IMAGE.
SSPACE is then run over the assembly to determine if any further joins can be made between the contigs
Gapfiller is then used to fill in the joins created with SSPACE
Finally ICORN is run in an attempt to fix any indels or SNP errors introduced during the assembly process

Optionally contigs less than 1kb are removed from the final assembly using mira_convert 
To retain contigs <1kb in length please set the -l flag to yes 

_EOF_
}


if  [ $# -lt 1 ]
    then
	    usage
		exit 1
fi

#Define path to MGAP install
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

# source dependencies
source "$SCRIPTPATH"/MGAP.config 
source "$SCRIPTPATH"/scheduler.config

#declare variables
declare -rx SCRIPT=${0##*/}

OPTSTRING="hr:s:l:"




declare SWITCH

#default behaviour options
ref=none
org=haploid
seq_directory="$PWD"
assemble=yes
pairing=PE
strain=all
long=no
h=yes
# Examine individual options
while getopts "$OPTSTRING" SWITCH; do 
		case $SWITCH in
		
		h) usage
           help 
		   exit 1
		   ;;
		   
		r) ref="$OPTARG"
		   echo "Reference = $ref"
		   ;;
         
		s) strain="$OPTARG"
		   echo "Only strain $strain will be processed. By default all strains within current directory will be processed"
		   ;;
		   
		l) long="$OPTARG"
		   echo "Contigs <1kb will be retained from final assembly = $long"
			;;
			
		\?) usage
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
if [ "$ref" != "none" ]; then


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
fi

VELVETG_TEST=`command -v "$VELVETG"`
VELVETH_TEST=`command -v "$VELVETH"`
JAVA_TEST=`command -v "$JAVA"`
SHUFFLE_TEST=`command -v "$SHUFFLE"`
ABACAS_TEST=`command -v "$ABACAS"`
VELVETOPT_TEST=`command -v "$VELVETOPT"`
IMAGE_TEST=`command -v "$IMAGE/image.pl"`
ICORN2_HOME_TEST=`command -v "$ICORN2_HOME/icorn2.sh"`
SSPACE_TEST=`command -v "$SSPACE"`

if [ -z "$VELVETG_TEST" ]; then
	    echo -e "ERROR: MGAP requires velvet to function. Please make sure the correct path is specified in MGAP.config\n"
		echo "MGAP is attempting to find velvetg here: $VELVETG"
		exit 1
fi
if [ -z "$VELVETH_TEST" ]; then
	    echo -e "ERROR: MGAP requires velvet to function. Please make sure the correct path is specified in MGAP.config\n"
		echo "MGAP is attempting to find velvetg here: $VELVETH"
		exit 1
fi
if [ -z "$JAVA_TEST" ]; then
	    echo "ERROR: MGAP requires java. Please make sure java is available on your system. The PATH to java can be modified in the MGAP.config file"
		exit 1
fi
if [ -z "$SHUFFLE_TEST" ]; then
	    echo "ERROR: MGAP requires shuffle sequences from velvet. Please make sure shuffle sequences is available on your system. The PATH to shuffle sequences can be modified in the MGAP.config file"
		exit 1
fi
if [ -z "$ABACAS_TEST" ]; then
	    echo "ERROR: MGAP requires ABACAS. Please make sure ABACAS is available on your system. The PATH to ABACAS can be modified in the MGAP.config file"
		exit 1
fi
if [ ! -f "$TRIM" ]; then
	    echo "ERROR: MGAP requires Trimmomatic. Please make sure trimmomatic is available on your system. The PATH to trimmomatic can be modified in the MGAP.config file"
		exit 1
fi
if [ -z "$VELVETOPT_TEST" ]; then
	    echo "ERROR: MGAP requires velvet optimiser. Please make sure velvet optimiser is available on your system. The PATH to velvet optimiser can be modified in the MGAP.config file"
		exit 1
fi
if [ -z "$IMAGE_TEST" ]; then
	    echo "ERROR: MGAP requires Image. Please make sure Image is available on your system. The PATH to Image can be modified in the MGAP.config file"
	    echo "MGAP is attempting to find image here: $IMAGE"
		exit 1
fi
if [ -z "$ICORN2_HOME_TEST" ]; then
	    echo "ERROR: MGAP requires iCorn2. Please make sure iCorn2 is available on your system. The PATH to iCorn2 can be modified in the MGAP.config file"
		exit 1
fi
if [ -z "$SSPACE_TEST" ]; then
	    echo "ERROR: MGAP requires SSPACE. Please make sure SSPACE is available on your system. The PATH to SSPACE can be modified in the MGAP.config file"
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
        echo -e "Sequences have been loaded into MGAP\n"
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

#create directory structure

if [ ! -d "tmp" ]; then
	mkdir "$seq_directory"/tmp
fi
if [ ! -d "Assemblies" ]; then
	mkdir "$seq_directory"/Assemblies
fi

if [ -s qsub_ids.txt ]; then
    rm qsub_ids.txt
fi

##TODO
#include resource management for SGE and SLURM

assemble_single ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
    if [ ! -s ${PBS_O_WORKDIR}/Assemblies/${strain}_final.fasta ]; then
		echo -e "Submitting qsub job for de novo assembly of ${strain}\n"
        var="seq=${strain},ref=$ref,seq_path=$seq_directory,kmer=$kmer,long=$long,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N Assemble_${strain} -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=$NCPUS,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/assembler.sh`
        echo -e "aln_sequences\t$qsub_array_id" >> qsub_array_ids.txt
	fi
        
fi
if [ ! -s qsub_ids.txt ]; then
    if [ ! -s ${PBS_O_WORKDIR}/Assemblies/${sequences}_final.fasta ]; then
		echo -e "Submitting qsub job for de novo assembly of ${strain}\n"
	    var="seq=${strain},ref=$ref,seq_path=$seq_directory,kmer=$kmer,long=$long,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N Assemble_${strain} -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=$NCPUS,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/assembler.sh`
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
                var="seq=${sequences[$i]},ref=$ref,seq_path=$seq_directory,kmer=$kmer,long=$long,SCRIPTPATH=$SCRIPTPATH"
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
	    	    var="seq=${sequences[$i]},ref=$ref,seq_path=$seq_directory,kmer=$kmer,long=$long,SCRIPTPATH=$SCRIPTPATH"
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
