########################################################################################################################################
#Marten Boetzer BaseClear B.v. 14-07-2012                                                                                              #
#perl subscript estimate_insert_size_v2.pl                                                                                             #
#This script;                                                                                                                          #
#  -Estimates median insert size by mapping paired-reads on contigs                                                                    #
#  If both reads of a pair is mapped on the same contig, the insert size and orientation is estimated. If all pairs are collected,     #                            #
#  the median insert size for each orientation is estimated.                                                                           #
#  Users can speed up the process by including the number of pairs to use for insert size estimation using the -n option,              #
#  as well as filtering the contigs using the -m option.                                                                               #
#                                                                                                                                      #
#  To run this script;                                                                                                                 #
#  perl estimate_insert_size_v2.pl -c <contigfile> -1 <readfile1> -2 <readfile2> -n <number_of_pairs>                                  #
#                                  -m <min_contig_size> -g <mismatches> -p <prefix>                                                    #
#                                                                                                                                      #
#  This script finds pairs on same contig, after collecting the pairs the median insert size is estimated for each found orientation   #
#  Only contigs are used if larger than <min_contig_size> .                                                                            #
#                                                                                                                                      #
#  Output is the median insert size and a file with distribution of the insert size, per orientation and for all orientations together.#
#  Also, number of pairs for each found orientation (FR, RF, FF and RR) are given.                                                     #                          #
########################################################################################################################################

use FindBin qw($Bin);
use File::Path;
use strict;

my $inputparameters = "usage:\nperl estimate_insert_size_v2.pl -c <contigfile> -1 <readfile1> -2 <readfile2> -n <number_of_pairs> -m <min_contig_size> -p <prefix>\n";
$inputparameters.= "\n============ Parameters ============\n";
$inputparameters.= "-c  Fasta file containing sequences for alignment of the paired-reads\n";
$inputparameters.= "-1  First read-file of the pair in fastA or fastQ format, may be in .gz format but requires 'gunzip' to be installed\n";
$inputparameters.= "-2  Second read-file of the pair in fastA or fastQ format, may be in .gz format but requires 'gunzip' to be installed\n";
$inputparameters.= "-n  Number of read-pairs to align against the reference sequence of -c\n";
$inputparameters.= "-m  Minimum sequence size of the reference sequence of -c, sequences smaller than this size will be removed\n";
$inputparameters.= "-g  Number of allowed mismatches in Bowtie (default: 0 mismatches)\n";
$inputparameters.= "-p  Prefix of the ouput files\n";

require "getopts.pl";

use vars qw($opt_1 $opt_2 $opt_m $opt_o $opt_v $opt_p $opt_k $opt_a $opt_z $opt_s $opt_b $opt_n $opt_l $opt_x $opt_u $opt_t $opt_T $opt_g $opt_r $opt_d $opt_S $opt_c);
&Getopts('1:2:m:o:v:p:k:a:z:s:b:n:l:x:u:t:T:g:r:d:S:c:');

my ($min_contig_size, $numpairs) = (0,0);

my $contigfile = $opt_c;
my $fileA = $opt_1;
my $fileB = $opt_2;
my $mismatches = 0;
$numpairs = $opt_n if($opt_n);
$min_contig_size = $opt_m if($opt_m);
$mismatches = $opt_g if($opt_g);
my $prefix = $opt_p;


die "$inputparameters\nERROR: Can't find contig file: $contigfile -- fatal\n" if(! -e $contigfile);
die "$inputparameters\nERROR: Can't find read file 1: $fileA -- fatal\n" if(! -e $fileA);
die "$inputparameters\nERROR: Can't find read file 2: $fileB -- fatal\n" if(! -e $fileB);
die "$inputparameters\nERROR: Please specify a prefix for the output files -- fatal\n" if($prefix eq "");

die "$inputparameters\nERROR: You've inserted '$numpairs', which does not seem to be an valid number. Exiting.\n" if(!($numpairs>=0) || !($numpairs =~ /^\d+$/));
die "$inputparameters\nERROR: You've inserted '$min_contig_size', which does not seem to be an valid number. Exiting.\n" if(!($min_contig_size>=0) || !($min_contig_size =~ /^\d+$/));
die "$inputparameters\nERROR: Only up to 2 mismatches allowed -- fatal\n" if($mismatches > 2);

print "\ncontig = $contigfile\n";
print "file1 = $fileA\n";
print "file2 = $fileB\n";
print "pairs = $numpairs \n";
print "min_contig_size = $min_contig_size\n";
print "mismatches = $mismatches\n";
print "prefix of output = $prefix\n\n";

my ($direction, $insertsize);
mkpath('bowtieoutput');

my $readsfile = prepareReads($fileA,$fileB);
open (CONT, $contigfile) || die "Can't open contig file $contigfile\n";

my ($seq,$name)=("","");
my ($contignum, $pairstep) = (0,100);

print "Finding pairs...\n";
open (BOWCONT, ">bowtieoutput/bowtiecontig.fa");
while (<CONT>) {
  chomp;
  $seq.=$_ if(eof(CONT));
  if (/\>(\S+)/ || eof(CONT)){
    if($seq ne ""){
       $contignum++;
       if(length($seq) >= $min_contig_size){
         print BOWCONT ">contig$contignum\n$seq\n";
       }

       $name = "";
       $seq = "";
    }
    $name = $1;
  }
  else {
     $seq .= $_;
  }
}
close BOWCONT;

&mapWithBowtie($name, "bowtieoutput/bowtiecontig.fa",$readsfile);

sub prepareReads{
   my ($fileA, $fileB) = @_;
  my $fastq = 0;
  
  if ($fileA =~ /\.gz$/) {
    print "opening .gz file\n";
    open(TEST, "gunzip -c  $fileA |");
    $name = <TEST>;
    close TEST;
    $fastq = 1 if ($name =~ /^[@]/);
    open(FILEA, "gunzip -c $fileA |") || die "can't open .gz $fileA with gunzip";
  }else{
    open(TEST,"< $fileA");
    $name = <TEST>;
    close TEST;
    $fastq = 1 if ($name =~ /^[@]/);
    open(FILEA, "< $fileA");
  }
  if ($fileB =~ /\.gz$/) {
    open(FILEB, "gunzip -c  $fileB |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(FILEB, "< $fileB");
  }

  print "processing reads to proper format...\n";
  my ($count, $step) = (0, 1000000);
  my $readfile = "bowtieoutput/bowtiereads.fa";
  #return $readfile;
  open (BOWIN, ">$readfile") || die "Can't write to single file $readfile -- fatal\n";

  while(<FILEA>) {
    <FILEB>;
    last if(++$count >= $numpairs && $numpairs > 0);
    if($count == $step){
      CounterPrint("\tRead $count...");
      $step = $step + 1000000;
    }
    my $seq1 = <FILEA>;
    chomp $seq1;
    my $seq2 = <FILEB>;
    chomp $seq2;
    $seq1 =~ s/\r//g;
    $seq2 =~ s/\r//g;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);
    print BOWIN ">read$count\n$seq1\n>read$count\n$seq2\n";
  }
  close BOWIN;
  close FILEA;
  close FILEB;
  print "\n";
  return $readfile;
}

sub mapWithBowtie{
  my ($curcontig, $contig,$readinput) = @_;
  my $bowtieout = "contig.bowtieIndex";
  my $Bin2 = substr($Bin, 0, -5);
  my $bowtiepath = "$Bin2"."/bowtie/bowtie";
  $bowtiepath =~ s/ /\\ /g;
  my $bowbuildpath = $bowtiepath."-build";
  system("$bowbuildpath $contig bowtieoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?"; # returns exit status values
  open(IN, "$bowtiepath -p 1 -v $mismatches -m 1 bowtieoutput/$bowtieout --suppress 6,7 -f $readinput --quiet -S --sam-nohead|") || die "Can't open bowtie output -- fatal\n";
  my ($prevread, $prevline);
  my ($foundoriHash,$foundoriPairHash, $allinsertsizes);
  my $founddistanceHash;
  my $orihash;
  $orihash->{0} = "F";
  $orihash->{16} = "R";
  open (CSV, ">distribution.txt") || die "Can't open distribution.txt for writing -- fatal";

  my ($numpairs,$numcorrectpairs)=(0,0);
  while(my $line1 = <IN>){
      my $line2 = <IN>;
      my @t1 = split(/\t/,$line1);
      next if($t1[2] eq "*");
      my @t2 = split(/\t/,$line2);
      next if($t2[2] eq "*" || $t1[2] ne $t2[2]);

      my ($start1, $start2) = ($t1[3],$t2[3]);
      my ($dir1, $dir2) = ($orihash->{$t1[1]},$orihash->{$t2[1]});
      my ($len1, $len2) = ((length($t1[9])),(length($t2[9])));

      if($t1[3] > $t2[3]){
        ($start1, $start2) = ($t2[3],$t1[3]);
        ($dir1, $dir2) = ($orihash->{$t2[1]},$orihash->{$t1[1]});
        my ($len1, $len2) = ((length($t2[9])),(length($t1[9])));

      }
      my $finalstart1 = $start1;
      my $finalend2 = $start2+$len2;

      my $diff = $finalend2-$finalstart1;

      $foundoriHash->{"$dir1$dir2"}++;
      $foundoriPairHash->{"$dir1$dir2"}{$diff}++;
      $allinsertsizes->{$diff}++;
  }


  open (FULL, ">$prefix.distribution_ALL.txt") || die "Can't open $prefix.distribution_ALL.txt for writing -- fatal";
  print FULL "insertsize\ttimes\n";
  foreach my $insert (sort {$a<=>$b} keys %$allinsertsizes){
    print FULL "$insert\t$allinsertsizes->{$insert}\n";
  }
  close FULL;


  foreach my $ori (keys %$foundoriPairHash){
    my $totalpairs = $foundoriHash->{$ori};

    my ($median_ins,$record) = (0,0);
   my ($total_is, $overal_is,$median_ins, $stdev,$sumX,$sumX2) = (0,0,0,0,0,0);
    my $median_bin = int($totalpairs/2);
    my $insert_size =$foundoriPairHash->{"$ori"};
    open (CSV, ">$prefix.distribution_$ori.txt") || die "Can't open $prefix.distribution_$ori.txt for writing -- fatal";
    print CSV "insertsize\ttimes\n";
    foreach my $is (sort {$a<=>$b} keys %$insert_size){
      for(my $i=0;$i<$insert_size->{$is};$i++){
        $record++;

       $sumX += $is;
       $sumX2 += ($is * $is);
        $median_ins = $is if($record >= $median_bin && $median_ins == 0);
      }
     $overal_is += ($is * $totalpairs);
      print CSV "$is\t$insert_size->{$is}\n";
    }
    close CSV;
    

   my ($mean_ins,$sigma) = (0,0);
   if($sumX > 0 && $record > 0){
     $mean_ins = int($sumX/$record);
     $sigma = sprintf("%.2f",sqrt($sumX2/$record - $mean_ins*$mean_ins));
   }
   close CSV;

    print "Summary: $ori\n";
    print "\tpairs = $totalpairs\n";
    print "\tmedian insert size = $median_ins\n";
    print "stdev = $sigma\n";

  }


}
###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}
