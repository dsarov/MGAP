#!/bin/bash

readRoot=$1
fragmentSize=$2
referenceORG=$3
start=$4
end=$5
capillary=$6 


User=$(whoami)

if [ -z "$end" ] ; then
	echo "Oh $User, you have to call me:"
	echo "<readroot> <fragementSize> <reference> <iterationstart> <iterationStop>"
	exit 1;
	
fi
 
read1="$readRoot"_1.fastq
read2="$readRoot"_2.fastq

if [ ! -f "$referenceORG" ] ; then
	echo "$User, where is the correct reference?"
	exit 1; 
fi

### getting a local copy of the reference, without | in the name
refRoot=$(echo $referenceORG | perl -nle '@ar=split(/\//); print $ar[($ar[scalar(@ar)]-1)]');

if [ ! -f "ICORN2.$refRoot.$start" ] ; then
cat $referenceORG | perl -e 'while(<STDIN>) {if (/>(\S+)/){ s/\|/_/; print "$_\n";} else {print}}' > ICORN2.$refRoot.$start
fi

reference=ICORN2.$refRoot

tmp=$$

echo $reference 
echo
for ((i=$start;$i<=$end;i++)); do
	directory=ICORN2_$i

	### call the mapper
	 if [ "$i" = "$start" ] ; then
	    /nfs/users/nfs_m/mh12/git/python3/map_splitter.py -o " -r 0 -y 0.7 " --split 20000 -c ICORN2.MAP.$$.$i -k 20 -s 3 smalt $reference.$i $directory $capillary 
	else
	    /nfs/users/nfs_m/mh12/git/python3/map_splitter.py -o " -r 0 -y 0.7 " --split 20000  --depend ICORN.CORRECT.$$.$(($i-1)) -c ICORN2.MAP.$$.$i -k 20 -s 3 smalt $reference.$i $directory $capillary
	fi
	
	### call the SNP caller
	bsub -w "ICORN2.MAP.$$.$i"  -J"ICORN2.SNP.$$.$i" -o out.ICORN2.SNP.$$.$i.o -e out.ICORN2.SNP.$$.$i.e -R "select[type==X86_64 && mem >4000] rusage[mem=4000]" -M4000000   -n6 -R "span[hosts=1]" ~tdo/Bin/icorn2.startSNP.sh $directory
	
	### call icorn starter
	bsub -w "ICORN2.SNP.$$.$i"  -J"ICORN.CORRECT.$$.$i" -o out.ICORN2.CORRECT.$$.$i.o -e out.ICORN2.CORRECT.$$.$i.e -R "select[type==X86_64 && mem >8000] rusage[mem=8000]" -M8000000 ~tdo/Bin/icorn2.correct.sh $directory $readRoot $fragmentSize $reference.$(($i+1))
	
done
	
