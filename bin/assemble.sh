#!/bin/bash

# Arguments passed from Nextflow
seq=$1
baseDir=$2
NCPUS=$3
long=$4
ref=$5
spades=$6
mem=$7


# Load dependencies if config exists, otherwise rely on Conda PATH
if [ -f "${baseDir}/configs/dependencies.config" ]; then
    source ${baseDir}/configs/dependencies.config
fi

# Ensure basic tools are set if variables aren't loaded from config
PILON=${PILON:-pilon}

echo "=========================================="
echo " Preparing Reads"
echo "=========================================="

##########################################################################
###                             ASSEMBLY                               ###
##########################################################################

if [ "$spades" == "true" ]; then
    echo "Running Spades..."
    # Using existing logic for Spades
    spades.py -o ${seq}_spades -1 ${seq}_1.fq.gz -2 ${seq}_2.fq.gz -t $NCPUS -m $mem
    mv ${seq}_spades/contigs.fasta ${seq}_assembly_raw.fasta
else
    echo "Running MEGAHIT (Replacing Velvet)..."
    
    # Clean up previous run if exists
    rm -rf ${seq}_megahit_out
    
    # Run MEGAHIT
    megahit -1 ${seq}_1.fq.gz -2 ${seq}_2.fq.gz -o ${seq}_megahit_out -t $NCPUS
    
    # Move output to standard name
    mv ${seq}_megahit_out/final.contigs.fa ${seq}_assembly_raw.fasta
fi

##########################################################################
###                      SCAFFOLDING (RagTag)                          ###
###          Replaces: GapFiller, ABACAS, IMAGE, SSPACE                ###
##########################################################################

# RagTag requires a reference. If 'ref' is 'none' or missing, we skip scaffolding.
if [ "$ref" != "none" ] && [ -s "$ref" ]; then
    echo "Reference found: $ref"
    echo "Running RagTag to scaffold contigs..."
    
    # ragtag.py scaffold <reference> <query_assembly>
    ragtag.py scaffold -t $NCPUS -o ${seq}_ragtag_out "$ref" ${seq}_assembly_raw.fasta
    
    # Rename RagTag output
    mv ${seq}_ragtag_out/ragtag.scaffold.fasta ${seq}_scaffolded.fasta
    
else
    echo "No reference provided (or ref=none). Skipping RagTag scaffolding."
    # Just pass the raw assembly to the next step
    mv ${seq}_assembly_raw.fasta ${seq}_scaffolded.fasta
fi

##########################################################################
###                  FILTERING (<1kb Removal)                          ###
##########################################################################

# We use 'seqtk' here instead of the old CONVERT_PROJECT legacy tool
# seqtk is standard in bioconda and likely already in your env

if [ "$long" == "no" ]; then
    echo "Filtering contigs < 1000bp..."
    seqtk seq -L 1000 ${seq}_scaffolded.fasta > ${seq}_pilon_input.fasta
    echo -e "Filtered removed contigs less than 1kb.\n" 
else 
    mv ${seq}_scaffolded.fasta ${seq}_pilon_input.fasta
    echo -e "Keeping all contigs.\n" 
fi

##########################################################################
###                           PILON (Polishing)                        ###
##########################################################################

echo "Running Pilon Polishing..."

# 1. Index the assembly
bwa index ${seq}_pilon_input.fasta

# 2. Map reads to the assembly (Piping for speed/disk space)
bwa mem -t $NCPUS ${seq}_pilon_input.fasta ${seq}_1.fq.gz ${seq}_2.fq.gz \
    | samtools sort -@ $NCPUS -o ${seq}.bam -

# 3. Index BAM
samtools index ${seq}.bam

# 4. Run Pilon
# Note: Pilon command varies slightly depending on installation.
# If installed via Conda, 'pilon' command works. 
# If $PILON var is set (from config), we use that.

echo "Executing Pilon..."
# We use a java max heap limit of $mem (passed from Nextflow)
# If PILON variable contains "jar", use java -jar, else assume it's the executable wrapper
if [[ "$PILON" == *".jar"* ]]; then
    java -Xmx${mem}G -jar ${PILON} --genome ${seq}_pilon_input.fasta --frags ${seq}.bam --output ${seq}_final
else
    pilon --genome ${seq}_pilon_input.fasta --frags ${seq}.bam --output ${seq}_final -Xmx${mem}G
fi

# Pilon outputs with a .fasta extension automatically, so the result is ${seq}_final.fasta

# Cleanup
rm ${seq}.bam ${seq}.bam.bai
rm -rf ${seq}_megahit_out ${seq}_ragtag_out

exit 0