#!/bin/bash

cat << _EOF_

This script will take a directory of assembles and summarise assemble quality for the following metrics: N50, genome size and contig number

The script can be run without any flags and assumes you are already in the correct directory

Please make sure you are in the directory containing the fasta assemblies and that all assemblies have the .fasta extension

Any genomes that are missing this extension will not be included

Type Y to continue or anything else to quit

_EOF_

read ref_test
if [ "$ref_test" == "Y" -o "$ref_test" == "y" -o "$ref_test" == "yes" ]; then
  echo -e "Continuing\n\n"
else
  exit 1
fi

module load perl/5.26.0

for f in *.fasta; do 
var=$(/home/dsarovich/bin/contigs_N50.pl $f) 
N50=$(echo $var | awk '{ print $2 }')
contigs=$(echo $var | awk '{ print $4 }')
length=$(echo $var | awk '{ print $7 }')

echo -e "$f\t$N50\t$contigs\t$length" > $f.tmp.summary

done

echo -e "\tN50\tcontigs\tlength" > header

cat header *.tmp.summary > final_summary.txt

clean_up () {
rm *.tmp.summary
rm header

}

clean_up

exit