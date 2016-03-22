#!/bin/bash

readRoot=$1
fragmentSize=$2
referenceORG=$3
start=$4
end=$5
insertsizeMax=$6
tmp=$7

if [ -z "$insertsizeMax" ] ; then
	insertsizeMax=1000
fi 


if [ -z "$tmp" ] ; then
   tmp=$$
fi 

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

tmp=$tmp

echo $reference 
echo
for ((i=$start;$i<=$end;i++)); do
	directory=ICORN2_$i

	### call the mapper
	 if [ "$i" = "$start" ] ; then
	    /nfs/users/nfs_m/mh12/git/python3/map_splitter.py -o " -i 10000 -r 0 -x -y 0.8 " -c ICORN2.MAP.$tmp.$i -k 13 -s 3 smalt $reference.$i $directory $read1 $read2
	else
	    /nfs/users/nfs_m/mh12/git/python3/map_splitter.py -o " -i 10000 -r 0 -x -y 0.8 " --depend ICORN.CORRECT.$tmp.$(($i-1)) -c ICORN2.MAP.$tmp.$i -k 13 -s 3 smalt $reference.$i $directory $read1 $read2
	fi
	
	### call the SNP caller
	bsub -w "ICORN2.MAP.$tmp.$i"  -J"ICORN2.SNP.$tmp.$i" -o out.ICORN2.SNP.$tmp.$i.o -e out.ICORN2.SNP.$tmp.$i.e -R "select[type==X86_64 && mem >4000] rusage[mem=4000]" -M4000000   -n6 -R "span[hosts=1]" ~tdo/Bin/icorn2.startSNP.sh $directory
	
	### call icorn starter
	bsub -w "ICORN2.SNP.$tmp.$i"  -J"ICORN.CORRECT.$tmp.$i" -o out.ICORN2.CORRECT.$tmp.$i.o -e out.ICORN2.CORRECT.$tmp.$i.e -R "select[type==X86_64 && mem >12000] rusage[mem=12000]" -M12000000 ~tdo/Bin/icorn2.correct.sh $directory $readRoot $fragmentSize $reference.$(($i+1))
	
done
	
