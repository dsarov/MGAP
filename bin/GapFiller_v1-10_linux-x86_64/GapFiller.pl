#!/usr/bin/env perl
#GapFiller v1.10 Marten Boetzer - Walter Pirovano, July 2012

#AUTHORS
# Marten Boetzer and Walter Pirovano (c) 2011-2012

#CITATION
#If you use GapFiller in a scientific publication, please cite:
# Boetzer, M. and Pirovano, W., Toward almost closed genomes with GapFiller, Genome Biology, 13(6), 2012
# walter.pirovano@baseclear.com

#LICENSE
#   GapFiller Copyright (c) 2011-2012 BaseClear B.V. The Netherlands. All rights reserved.
#   GapFiller can be freely used by academic institutes or non-profit organizations.
#   Commercial parties need to acquire a license. For more information about commercial licenses look at our website or email info@baseclear.com.

#TERMS
#   This software comes without any warranty. BaseClear cannot guarantee that this software will fulfill any of your particular purposes or needs.
#   BaseClear will not be liable for data loss or any damages related to the software.
#   This software and eventual license herein granted shall not be copied, shared or offered for re-sale without the explicit approval of BaseClear.

#DOCUMENTATION
#   README, MANUAL and TUTORIAL distributed with this software @ www.baseclear.com
#   Boetzer, M. and Pirovano, W., Toward almost closed genomes with GapFiller, Genome Biology, 13(6), 2012
#   http://www.baseclear.com/sequencing/data-analysis/bioinformatics-tools/
#   We hope this code is useful to you -- Please send comments & suggestions to Walter.Pirovano@baseclear.com
#   If you use either the GapFiller code or ideas, please cite our work appropriately and accurately


###CHANGES in version 1.10:
#  Fixed a bug where internal sequences were reverse complemented

###CHANGES in version 1.9:
#  Included the option -S to skip reading of the files if it was already done in a previous analysis
#  Hashing of the reads is now done during filling the gaps, instead of during the mapping stage. This saves memory and time. 

  use strict;
  use File::Basename;
  use Storable;
  use File::Path;
  use File::Copy;
  use FindBin qw($Bin);
  use threads;
  use Text::Wrap;
  use Getopt::Std;
  use threads::shared;
  my $totalReadsProcessed :shared;
  my $totalReadFiles :shared;

  $Text::Wrap::columns = 61;

use vars qw($opt_m $opt_o $opt_v $opt_p $opt_k $opt_a $opt_z $opt_s $opt_b $opt_n $opt_l $opt_x $opt_u $opt_t $opt_T $opt_g $opt_r $opt_d $opt_S $opt_i);
getopts('m:o:v:p:k:a:z:s:b:n:l:x:u:t:T:g:r:d:S:i:');
my ($keep, $trim, $gaps, $threads, $difference,$base_overlap,$min_overlap,$min_base_ratio,$base_name, $min_tig_overlap, $numiteration)= (0, 10,1,1, 50, 2,29, 0.7, "standard_output", 10,10);

  my $seplines = ("-" x 60)."\n";
  my $version = "v1.11";
#-------------------------------------------------READ OPTIONS

use constant USAGE =><<EOH;

Usage: $0 $version


============ General Parameters ============
-l  Library file containing two paired-read files with insert size, error and orientation indication.
-s  Fasta file containing scaffold sequences used for extension.

============ Extension Parameters ============
-m  Minimum number of overlapping bases with the edge of the gap (default -m $min_overlap)
-o  Minimum number of reads needed to call a base during an extension (default -o $base_overlap)
-r  Percentage of reads that should have a single nucleotide extension in order to close a gap in a scaffold (Default: $min_base_ratio)
-d  Maximum difference between the gapsize and the number of gapclosed nucleotides. Extension is stopped if it matches this parameter + gap size (default -d $difference, optional).
-n  Minimum overlap required between contigs to merge adjacent sequences in a scaffold (default -n $min_tig_overlap, optional)
-t  Number of reads to trim off the start and begin of the sequence (usually missambled/low-coverage reads) (default -t $trim, optional)
-i  Number of iterations to fill the gaps (default -i $numiteration, optional)

============ Bowtie Parameters ============
-g  Maximum number of allowed gaps during mapping with Bowtie. Corresponds to the -v option in Bowtie. (default -g $gaps, optional)

============ Additional Parameters ============
-T  Number of threads to run (default -T $threads)
-S  Skip reading of the input files again
-b  Base name for your output files (optional)

EOH


if(!($opt_l) || !($opt_s)) {
	print "ERROR: Parameter -l is required. Please insert a library file\n" if(!$opt_l);
	print "ERROR: Parameter -s is required. Please insert a scaffold fastA file\n" if(!$opt_s);
	
	die USAGE;
}

my $scaffold = $opt_s if($opt_s);
my $libraryfile;
$libraryfile = $opt_l if ($opt_l);
$min_overlap = $opt_m if ($opt_m);
$base_overlap = $opt_o if ($opt_o);
$base_name = $opt_b if($opt_b);
$min_tig_overlap = $opt_n if($opt_n);
$min_base_ratio = $opt_r if ($opt_r);
$threads = $opt_T if ($opt_T);
$gaps = $opt_g if ($opt_g || $opt_g eq 0);
$difference = $opt_d if($opt_d);
$trim = $opt_t if($opt_t || $opt_t eq 0);
$keep = $opt_S if($opt_S || $opt_S eq 0);
$numiteration = $opt_i if($opt_i);
mkpath("$base_name");
mkpath("$base_name/reads");

  my $files = "";
  my $textline = "Your inserted inputs on $version at ".getDate().":\n\t-s $scaffold\n\t-l $libraryfile\n\t-b $base_name\n\t-o $base_overlap\n\t-m $min_overlap\n";
  $textline .= "\t-r $min_base_ratio\n\t-n $min_tig_overlap\n\t-T $threads\n\t-g $gaps\n\t-d $difference\n\t-t $trim\n\t-i $numiteration\n\n";



  my $summaryfile = "$base_name/$base_name.summaryfile.final.txt";

  open (SUMFILE, ">$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  print SUMFILE "$textline";
  print "$textline";
  &printMessage("\n=>".getDate().": Reading and processing paired-read files\n");

  #process the library and read the input files with multiple threads
  open(FILELIB, "< $libraryfile") || die "Can't open $libraryfile -- fatal\n";
  my ($prevlib,$sequenceReadsInGaps,$filehash, $maxGapScaf);
  my $ctlib = 0;
  my $totalfiles = 0;
  while(<FILELIB>){
    my ($lib, $aligner, $fileA, $fileB, $insert_size, $insert_stdev, $orient) = split(/\s+/, $_);
    next if($lib eq '');
    die "ERROR: Invalid aligner in library $lib: $aligner. Should be either 'bowtie', 'bwa' or 'bwasw' -- fatal\n" if($aligner ne "bwa" && $aligner ne "bwasw" && $aligner ne "bowtie");
    die "ERROR: Invalid file in library $lib: $fileA -- fatal\n" if(!(-e $fileA));
    die "ERROR: Invalid file in library $lib: $fileB -- fatal\n" if(!(-e $fileB));
    die "ERROR: Orientation must have length of 2 characters and should contain one of the following; FR, FF, FR or RF. Your library $lib has orientation of $orient ...Exiting.\n" if(!(length($orient) == 2) || !($orient =~ /[FR][FR]/));

    $ctlib++;
    $prevlib = $lib;
    if(!$keep){
      my $thr = threads->create(\&generateInputFiles, $ctlib, $fileA, $fileB, "$aligner.reads.lib$ctlib", $orient);
      if(!($ctlib % $threads)){
        foreach my $thr (threads->list()) {
          my @res = $thr->join();
          my ($libs,$reads) = split(":",$res[0]);
          my @readarray = split(",",$reads);
          foreach my $readf (@readarray){
            $filehash->{$libs}{$readf}++;
          }
        }
      }
    }else{
      my @readfiles = <$base_name/reads/$base_name.$aligner.reads.lib$ctlib*>;
      if($#readfiles < 0){
        die "ERROR: no readfiles found for library $lib line:\n$_\nYou have inserted option -S 1, change this to -S 0\n";
      }
      foreach my $readfile (@readfiles){
        $filehash->{$ctlib}{$readfile}++;
      }
    }
  }
  if(!$keep){
    foreach my $thr (threads->list()) {
      my @res = $thr->join();
      my ($libs,$reads) = split(":",$res[0]);
      my @readarray = split(",",$reads);
      foreach my $readf (@readarray){
        $filehash->{$libs}{$readf}++;
      }
    }
  }
  close FILELIB;

###get the total number of files of all libraries
  foreach my $libs(keys %$filehash){
    my $filelist = $filehash->{$libs};
    foreach my $readfiles(keys %$filelist){
      $totalReadFiles++;
    }
  }
  #exit;
  print "\n";
  my $prevscaffold = $scaffold;
  my $iteration = 0;
  my ($finalFile, $finalSummary, $finalClose, $finalFill);
  mkpath("$base_name/intermediate_results");
  while(++$iteration <= $numiteration){
    $totalReadsProcessed=0;
    $scaffold = $prevscaffold;
    print SUMFILE "After iteration $iteration:\n$seplines\n\n";
    print "ITERATION $iteration:\n";
    my $closefile = "$base_name/intermediate_results/$base_name.closed.evidence.iteration$iteration.txt";
    my $fillFile = "$base_name/intermediate_results/$base_name.filled.iteration$iteration.txt";

    open (CLOSEFILE, ">$closefile") || die "Can't open $closefile -- fatal\n";
    $ctlib=0;
    open(FILELIB, "< $libraryfile") || die "Can't open $libraryfile -- fatal\n";
    &printMessage("\n=>".getDate().": Mapping reads to scaffolds, reading bowtie output and storing unmapped reads\n");
    my $ctthreads=0;
  #read the libraries again and map the reads to the scaffolds
  ###TO DO: able to skip this if already mapped, make an option named 'steps'!!!

    system("rm -rf $base_name/alignoutput/*");  #remove previous folder, if exists
    mkpath("$base_name/alignoutput");
    $prevlib = "";
    while(<FILELIB>){
      my ($lib, $aligner, $fileA, $fileB, $insert_size, $insert_stdev, $reverse) = split(/\s+/, $_);
      next if($lib eq '' || $prevlib eq $lib || $fileA eq "TAB");
      $ctlib++;
      $prevlib = $lib;
      my $min_allowed = -1 * ($insert_stdev * $insert_size);
      my $maxlib = int($insert_size - $min_allowed);
      my $contigFile = "$base_name/alignoutput/$base_name.$ctlib.gapclosure.fa";

      processScaffold($scaffold, $contigFile, $maxlib, $trim);
      my $library = "$lib.$ctlib";
      die "Scaffold file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));

#BWA!
      if($aligner eq "bwa" || $aligner eq "bwasw"){
        my @readfiles= <$base_name/reads/$base_name.$aligner.reads.lib$ctlib*>;
        my $bwapath = "$Bin"."/bwa/bwa";
        $bwapath =~ s/ /\\ /g;
        my $bwaout = $base_name . ".$library.bwaIndex";
        &printMessage("\n=>".getDate().": Building BWA index for library $lib\n");
        my $filesize = -s "$contigFile";
        my $index = "bwtsw";
        $index = "is" if($filesize <= 10000000);
        open(STDERR, ">$base_name/alignoutput/tmpbwa_logfile");
        system("$bwapath index -a $index $contigFile -p $base_name/alignoutput/$bwaout") == 0 || die "\nBwa error; $?"; # returns exit status values
        my $filenum=1;
        foreach my $readfile (@readfiles){
          if($aligner eq "bwa"){
            my $bwaoutputaln = "$base_name/alignoutput/$base_name.$library.$ctlib.$filenum.bwa";
            my $procline = "$bwapath samse $base_name/alignoutput/$bwaout $bwaoutputaln $readfile |";
            my $samseline = "$bwapath aln -i 0 $base_name/alignoutput/$bwaout $readfile > $bwaoutputaln";
            my $thr = threads->create(\&getUnmappedReads, $maxlib, $procline, $samseline);
          }else{
            my $procline = "$bwapath bwasw $base_name/alignoutput/$bwaout $readfile |";
            my $thr = threads->create(\&getUnmappedReads, $maxlib, $procline);
          }
          $filenum++;
          getUnmappedThreadResult() if(!(++$ctthreads % $threads));
        }
      }
      if($aligner eq "bowtie"){
        my @readfiles = <$base_name/reads/$base_name.$aligner.reads.lib$ctlib*>;
        my $bowtieout = "$base_name.$library.bowtieIndex";
        my $bowtiepath = "$Bin"."/bowtie/bowtie";
        $bowtiepath =~ s/ /\\ /g;
        my $bowbuildpath = $bowtiepath."-build"; 
        system("$bowbuildpath $contigFile $base_name/alignoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?"; # returns exit status values
        foreach my $readfile (@readfiles){
          my $procline = "$bowtiepath -v $gaps -m 1 $base_name/alignoutput/$bowtieout -f $readfile -S --sam-nohead --quiet |";
          my $thr = threads->create(\&getUnmappedReads, $maxlib, $procline);
          getUnmappedThreadResult() if(!(++$ctthreads % $threads));
        }
      }
      $finalSummary = $summaryfile;
      $finalClose = $closefile;
      $finalFill = $fillFile;
    }
    close FILELIB;
    getUnmappedThreadResult();
    CounterPrint("                                                     ");
  #=cut till this part should be included with the 'steps'

    $finalFile = "$base_name/intermediate_results/$base_name.gapfilled.iteration$iteration.fa";
    my $preNcount = goThroughScaffold($scaffold, $finalFile, $fillFile);

    my ($sumorig) = writesummaryfiles($scaffold);
    my ($sumclosed,$afterNcount) = writesummaryfiles($finalFile);
    
    my $nucClosed = ($preNcount - $afterNcount);
    
    print SUMFILE "\tClosed $nucClosed out of $preNcount nucleotides\n\n$seplines\n";
    print "\tClosed $nucClosed out of $preNcount nucleotides\n\n";

    print SUMFILE "\n$seplines Scaffold statistics:\n$seplines";
    print SUMFILE "\n\tOriginal:\n$sumorig";
  
    print SUMFILE "\tAfter gapclosing:\n$sumclosed$seplines\n";
    $prevscaffold = $finalFile;
    close CLOSEFILE;
    $trim = 0;
    
    if($preNcount == $nucClosed){
      print "\tAll gaps are closed. Exiting...\n\n";
      $iteration++;
      last;
    }
    if($nucClosed == 0){
      print "\tCould not close any more gaps. Exiting...\n\n";
      $iteration++;
      last;
    }
  }
  --$iteration;
  #my ($finalFile, $finalSummary, $finalClose, $finalFill);
  copy("$base_name/intermediate_results/$base_name.closed.evidence.iteration$iteration.txt", "$base_name/$base_name.closed.evidence.final.txt");
  copy("$base_name/intermediate_results/$base_name.filled.iteration$iteration.txt", "$base_name/$base_name.filled.final.txt");
  copy("$base_name/intermediate_results/$base_name.gapfilled.iteration$iteration.fa", "$base_name/$base_name.gapfilled.final.fa");




  my $time = (time - $^T);
  my $minutes = int ($time / 60);
  $time = $time % 60;
  print "Process run succesfully on ".getDate()." in $minutes"." minutes and $time"." seconds\n\n\n";
  print SUMFILE "\nProcess run succesfully on ".getDate()." in $minutes"." minutes and $time"." seconds\n\n\n";
  close SUMFILE;
  close CLOSEFILE;

#####################################    END OF MAIN SCRIPT

#process the read files and cut it into subsets of 1 million pairs
sub generateInputFiles{
  my ($lib,$fileA, $fileB, $outfile, $orient, $fname) = @_;
  my ($name,$seq1,$seq2);
  my ($counterext, $Ncount, $countsinglet, $fastq, $step) = (0,0,0,0,1000000);

  my ($ori_1,$ori_2) = split(//, $orient);
  open (OUTSINGLEFILE, ">$base_name/reads/$base_name.$outfile.1") || die "Can't write to single file $base_name/reads/$base_name.$outfile.1 -- fatal\n";
  my $files="$lib:$base_name/reads/$base_name.$outfile.1";
  if ($fileA =~ /\.gz$/) {
    open(TEST, "gunzip -c  $fileA |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(TEST, "< $fileA");
  }
  $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);
  if ($fileA =~ /\.gz$/) {
    open(FILEA, "gunzip -c  $fileA |") || die "can't open .gz $fileA with gunzip";
  }else{
    open(FILEA, "< $fileA");
  }
  if ($fileB =~ /\.gz$/) {
    open(FILEB, "gunzip -c  $fileB |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(FILEB, "< $fileB");
  }

  my $ctstep = 1;
  while(<FILEA>) {
    <FILEB>;
    $seq1 = uc(<FILEA>); chomp $seq1;
    $seq2 = uc(<FILEB>); chomp $seq2;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);

    if(++$countsinglet == $step){
      CounterPrint("reading files of library line $lib @ $countsinglet");
      $step = $step + 1000000;
      close OUTSINGLEFILE;
      $ctstep++;
      open (OUTSINGLEFILE, ">$base_name/reads/$base_name.$outfile.$ctstep") || die "Can't write to single file $base_name/reads/$base_name.$outfile.$ctstep -- fatal\n";
      $files.=",$base_name/reads/$base_name.$outfile.$ctstep";
    }
    $seq1 = reverseComplement($seq1) if($ori_1 eq "R");
    $seq2 = reverseComplement($seq2) if($ori_2 eq "F");
    print OUTSINGLEFILE ">read$countsinglet.1\n$seq1\n>read$countsinglet.2\n$seq2\n";
  }
  CounterPrint("                                            ");
  close OUTSINGLEFILE;
  close FILEB;
  close FILEA;
  return $files;
}

#process the scaffold file, retreive only the edges of the sequences surrounding the gap.
sub processScaffold{
  my ($scaffile,$processedFile, $maxlib, $trim) = @_;
  my ($seq, $name, $scafct) = ("","",0,0);
  open(SCAF,"$scaffile") || die "Can't open $scaffile -- fatal\n";
  open (BOWIN, ">$processedFile") || die "Can't write to file BOWIN -- fatal\n";

  while(<SCAF>){
    chomp;
    my $line = $_;
    $seq.= uc($line) if(eof(SCAF));
    if (/\>(\S+)/ || eof(SCAF)){
      if($seq ne ''){
        $scafct++;
        my @seqgaps = split(/[N]{1,}/, $seq);
        if($#seqgaps > 0){
          my $ctgcount=0;
          my $numgap = $#seqgaps;
          $maxGapScaf->{"scaf$scafct"} = $numgap;
          foreach my $ctgseq (@seqgaps){
            $ctgseq = substr($ctgseq,$trim) if($ctgcount > 0);
            $ctgseq = substr($ctgseq,0,length($ctgseq)-$trim);
            $ctgcount++;
            my $len = length($ctgseq);
            if($len > (($maxlib * 2))){
              my $upper = (length($ctgseq) - ($maxlib));
              my $first = substr($ctgseq, 0, $maxlib);
              my $second = substr($ctgseq, $upper);
              my $newseq = $first."NNN".$second;
              $newseq = $first if($ctgcount > $numgap);
              print BOWIN ">scaf$scafct.$ctgcount.$len\n";
              print BOWIN "$newseq\n";
            }else{
              $ctgseq = substr($ctgseq, 0, $maxlib) if($ctgcount > $numgap);
              print BOWIN ">scaf$scafct.$ctgcount.$len\n";
              print BOWIN "$ctgseq\n";
            }
          }
        }
      }
      $seq='';
      $name = $1;
    }else{
      $seq.=uc($line);
    }
  }
  close SCAF;
  close BOWIN;
}

#get reads that are not mapped within a given distance (insertsize + (insertsize+deviation))
#store all reads per gap in a hash for later extension
sub getUnmappedReads{
  my ($maxlib, $procline, $aligninput) = @_;
  system("$aligninput") if($aligninput ne "");
  open(IN, "$procline") || die "Can't open baw output: Process: $procline -- fatal\n";
  my $usereadhash;
  my $sub = ($maxlib + 3);
  while(my $line1 = <IN>){
    next if($line1 =~ /^@/);
    my $line2 = <IN>;
    my @t1 = split(/\t/,$line1);
    my @t2 = split(/\t/,$line2);
    next if($t1[2] eq "*" && $t2[2] eq "*");
    my ($scaf1,$gap1, $len1) = split(/\./, $t1[2]);
    my ($scaf2,$gap2, $len2) = split(/\./, $t2[2]);
    if($t1[3] > $maxlib && $len1 > (($maxlib * 2))){
      my $minsub = $len1 - (($maxlib*2)+3);
      $t1[3] = ($t1[3] + $minsub);
    }
    if($t2[3] > $maxlib && $len2 > (($maxlib * 2))){
      my $minsub = $len2 - (($maxlib*2)+3);
      $t2[3] = ($t2[3] + $minsub);
    }
    if($t1[2] ne "*" && ($t1[1] == 0 || $t1[1] == 16)){  #if first read is mapped on a contig
      my ($strand, $tig, $start, $seq) = ($t1[1],$t1[2],$t1[3],$t2[9]);
      if ($strand == 0) { ###step 1: first read is on + strand
        if(($start > (($len1-$maxlib)-length($seq)))){
          if($gap1 <= $maxGapScaf->{$scaf1}){
            if($t2[1] == 0 || $t2[1] == 16){
              if(($gap1+1) == $gap2 && $t2[3] < $maxlib+length($seq)){   #is read on next contig and within unallowed region?
                next;
              }
              elsif($gap1 == $gap2 && $t2[3] > (($len2-$maxlib)-length($seq))){    #is read on same contig and within unallowed region?
                next;
              }
            }
            $usereadhash->{"$scaf1.$gap1"}{$seq}++;
          }
        }
      }else{ ###step 2: first read is on - strand
        if(($start < ($maxlib+length($seq)))){
          my $prevgap1 = $gap1-1;
          if($prevgap1 >= 1){ #read is not on first contig of the scaffold, since pair can't be on contig before the first contig
            if($t2[1] == 0 || $t2[1] == 16){
              if(($prevgap1+1) == $gap2 && $t2[3] < $maxlib+length($seq)){ #is read on next contig and within unallowed region?
                next;
              }
              elsif($prevgap1 == $gap2 && $t2[3] > (($len2-$maxlib)-length($seq))){ #is read on same contig and within unallowed region?
                next;
              }
            }
            $usereadhash->{"$scaf1.$prevgap1"}{$seq}++;
          }
        }
      }
    }      #if second read is mapped on a contig
    if($t2[2] ne "*" && ($t2[1] == 0 || $t2[1] == 16)){
      my ($strand, $tig, $start, $seq) = ($t2[1],$t2[2],$t2[3],$t1[9]);
      if ($strand == 0) { ###step 1: first read is on + strand
        if(($start > (($len2-$maxlib)-length($seq)))){
          if($gap2 <= $maxGapScaf->{$scaf2}){
            if($t1[1] == 0 || $t1[1] == 16){
              if(($gap2+1) == $gap1 && $t1[3] < $maxlib+length($seq)){  #is read on next contig and within unallowed region?
                next;
              }
              elsif($gap2 == $gap1 && $t1[3] > (($len1-$maxlib)-length($seq))){  #is read on same contig and within unallowed region?
                next;
              }
            }
            $usereadhash->{"$scaf2.$gap2"}{$seq}++;
          }
        }
      }else{ ###step 2: first read is on - strand
        if(($start < ($maxlib+length($seq)))){
          $gap2--;
          if($gap2 >= 1){ #read is not on first contig of the scaffold, since pair can't be on contig before the first contig
            if($t1[1] == 0 || $t1[1] == 16){
              if(($gap2+1) == $gap1 && $t1[3] < $maxlib+length($seq)){  #is read on next contig and within unallowed region?
                next;
              }
              elsif($gap2 == $gap1 && $t1[3] > (($len1-$maxlib)-length($seq))){ #is read on same contig and within unallowed region?
                next;
              }
            }
            $usereadhash->{"$scaf2.$gap2"}{$seq}++;
          }
        }
      }
    }
  }
  close IN;
  $totalReadsProcessed++;
  my $perc = sprintf("%.1f", ($totalReadsProcessed/$totalReadFiles)*100);
  CounterPrint("...Processed ".($totalReadsProcessed*1000000)." paired-reads (~$perc%)");
  return $usereadhash;
}

#get the unmapped reads of all the threads and store it in a new hash, per gap
sub getUnmappedThreadResult{
  foreach my $thr (threads->list()) {
    my @ret = $thr->join();
    my $unmaphash = $ret[0];
    foreach my $scaf (keys %$unmaphash){
      my $list=$unmaphash->{$scaf};
     # print "SCAF = $scaf\n";
      open (OUT, ">>$base_name/alignoutput/$scaf.txt");
      my $num = 0;
      foreach my $gapread (keys %$list){
        $num++;

        print OUT (">$scaf.read\n$gapread\n" x $unmaphash->{$scaf}{$gapread});
      }
      close OUT;
    }
  }
}

sub goThroughScaffold{
  my ($scaffile, $finalfile, $fillFile) = @_;
  my ($seq, $prevhead, $ct) = ("","",0);
  open(SCAF,"$scaffile") || die "Can't open scaffile -- fatal\n";
  open (CLOSED, ">$finalFile") || die "Can't write to file -- fatal\n";
  open (FILL, ">$fillFile") || die "Can't write to file $fillFile-- fatal\n";
  my ($totalgap, $totalclosed, $scafct, $totalgapnuc)= (0,0,0,0);
  &printMessage("\n=>".getDate().": Filling gaps\n");
  while(<SCAF>){
    chomp;
    $seq.= uc($_) if(eof(SCAF));
    if (/\>(\S+)/ || eof(SCAF)){
      my $head = $1;
      if($seq ne ""){
        $scafct++;
        my @gap = split(/N+/, $seq);
        my @countgapper;
        while ($seq =~ /(N+)/g) {
          my $start = (pos($seq)-length($1));
          my $end = pos($seq);
          my $lengap = length($1);
          my $gapstring = "$start|$end|$lengap";
          push @countgapper, $gapstring;
        }
        my ($finalseq, $closed) = ("", 0);
        if($#gap >0){
          print FILL ">$prevhead:\n";
          $totalgap = $totalgap + $#gap;
          my ($totgapIncl, $gapnum) = (0,0);
          for (my $x=0;$x < $#gap;$x++){
            $gapnum++;
            my $scafgap = "scaf$scafct.$gapnum";
            CounterPrint("Filling: scaffold$scafct"."|gap$gapnum      ");
            print CLOSEFILE ">$prevhead"."|gap$gapnum\n";
            my $gapinfo = $countgapper[$x];
            my ($start,$end,$lengap) = split(/\|/,$gapinfo);
            $lengap = ($lengap+($trim*2));
            $totalgapnuc+=$lengap;
            $totgapIncl+=$trim;
            $start = $start-$totgapIncl;
            $end = $start+$lengap;
            print FILL "\tgapnum=$gapnum\t";
            print FILL "start=$start\t";
            print FILL "end=$end\t";
            print FILL "gapsize=$lengap\t";

            my $gap1 = $gap[$x];
            $gap1 = substr($gap1,0,length($gap1)-$trim);
            $gap1 = $finalseq if($x >0);
            my $gap2 = $gap[($x+1)];
            my $subgap1 = substr($gap1,-100);
            $gap2 = substr($gap2,$trim);
            $gap2 = substr($gap2,0,length($gap2)-$trim) if($gapnum < $#gap);;

            my ($gapseq, $closed) = fillGap($scafgap,$ct, $subgap1, $gap2, $lengap);
            my $subgap3 = substr($gapseq,100);
            print FILL "closed=$closed\n";

            $finalseq = $gap1.$subgap3 if($x == 0);
            $finalseq = $finalseq.$subgap3 if($x >0);
            $closed++ if($closed ne "no");
            $totalclosed++ if($closed ne "no");

            $ct++;
          }
          print CLOSED ">$prevhead\n".wrap('', '', $finalseq)."\n";
        }else{
          print CLOSED ">$prevhead\n". wrap('', '', $seq)."\n";
        }
      }
      $prevhead = $head;
      $seq='';
    }else{
      $seq.= uc($_);
    }
  }
  CounterPrint("                             ");
  print SUMFILE "GAPCLOSING STATS:\n$seplines";
  print SUMFILE "\n\tClosed $totalclosed out of $totalgap gaps \n";
  print "\n\tClosed $totalclosed out of $totalgap gaps \n";
  close SCAF;
  close EVI;
  close FILL;
  close CLOSED;
  return $totalgapnuc;
}

sub fillGap{
  my ($gapnum, $ct1, $ctg1, $ctg2, $gapsize) = @_;
  my $readfile = "$base_name/alignoutput/$gapnum.txt";
  my $seqstig;
  my $numreads = 0;
  open(READS,"$readfile");
  while(<READS>){
    my $read = <READS>;
    $numreads++;
    chomp $read;
    my $subct = 0;
    while($subct < length($read)-$min_overlap){
      my $subseq = substr($read, $subct, $min_overlap+1);
      if(index($subseq,"N") == -1){
        if(defined $seqstig->{$subseq}){
          $seqstig->{$subseq}++;
        }else{
          $subseq = reverseComplement($subseq);
          $seqstig->{$subseq}++;
        }
      }
      $subct++;
    }
  }
  print FILL "reads=$numreads\t";
  my $seqstig_backup = $seqstig;
  my $prevlen1 = length($ctg1);
  my $prevlen2 = length($ctg2);

  my ($newseq, $exitstatus) = MergeTwoSequences($ctg1, $ctg2);
  my $gapclosed = 0;
  my ($extend3output, $extend5output) = ("","");
  if($exitstatus <=0 || !(($gapclosed) >= $gapsize-$difference && ($gapclosed) < $gapsize+$difference)){
    $exitstatus = 0;
    my ($extended3prime, $extended5prime, $iteration, $prevclosed) = (1,1,1,0);
    while(($extended3prime || $extended5prime) && $exitstatus <= 0){
      ($seqstig, $newseq, $ctg1, $exitstatus,$gapclosed, $extended3prime, $extend3output) = gapFillExtension("3", $seqstig, $ctg1, $ctg2,$gapclosed,$gapsize);
      if($exitstatus <=0){
        last if(($gapclosed-$min_tig_overlap)>($gapsize+$difference));
        $ctg2 = reverseComplement($ctg2);
        $ctg1 = reverseComplement($ctg1);
        ($seqstig, $newseq, $ctg2, $exitstatus,$gapclosed, $extended5prime, $extend5output) = gapFillExtension("5", $seqstig, $ctg2, $ctg1,$gapclosed,$gapsize);
        $ctg2 = reverseComplement($ctg2);
        $ctg1 = reverseComplement($ctg1);
        $newseq= reverseComplement($newseq);
      }
      last if(($gapclosed-$min_tig_overlap)>($gapsize+$difference));
    }
  }
  if($exitstatus == 0){
    print CLOSEFILE "$extend3output" if($extend3output ne "");
    print CLOSEFILE "$extend5output" if($extend5output ne "");
    ($newseq, $exitstatus) = MergeTwoSequencesFinal($ctg1, $ctg2) ;
  }
  if($exitstatus <=0 || !(($gapclosed) >= $gapsize-$difference)){
    $exitstatus = 0;
    $newseq = "";
  }
  my $extsize1 =  (length($ctg1)-$prevlen1);
  my $extsize2 =  (length($ctg2)-$prevlen2);
  $extsize1 = 0 if($extsize1 <0);
  $extsize2 = 0 if($extsize2 <0);

  my $filled = ($extsize1+$extsize2)-$exitstatus;
  $filled = 0 if($filled<0);
  print FILL "filled=$filled\text1=$extsize1\text2=$extsize2\tmerged=$exitstatus\t";
  my $Ngap = 1;
  $Ngap = ($gapsize-$filled) if($filled < $gapsize);
  my $finalseq = $ctg1.("N" x $Ngap).$ctg2;
  if($newseq ne ""){
    $finalseq = $newseq;
    print CLOSEFILE "\tsequences are merged\n\n";
    print FILL "remaining=0\t";
  }else{
    print FILL "remaining=$Ngap\t";
    print CLOSEFILE "\tsequences are NOT merged\n\n";
  }
  my $closed = "no";
  $closed = "yes" if($exitstatus>0);
  return $finalseq, $closed;
}

sub gapFillExtension{

   my ($direction, $bin, $seq, $seq2, $totalclosed, $gapsize) = @_;
   my $extended = 0;
   my $subseq = substr($seq, -$min_overlap);
   my $stopoutput = "";

   my $revseq = reverseComplement($subseq);
   my $overhang;
   #get number of occurences for extension of A,C,G,T
   $overhang->{'A'} = $bin->{$subseq."A"}+$bin->{"T$revseq"};
   $overhang->{'C'} = $bin->{$subseq."C"}+$bin->{"G$revseq"};
   $overhang->{'G'} = $bin->{$subseq."G"}+$bin->{"C$revseq"};
   $overhang->{'T'} = $bin->{$subseq."T"}+$bin->{"A$revseq"};

   #obtain total coverage
   my $coverage = $overhang->{'A'}+$overhang->{'C'}+$overhang->{'G'}+$overhang->{'T'};
   my $info = "\tclosed=$totalclosed/$gapsize\tdir$direction: Total:$coverage\tA:$overhang->{'A'}\tT:$overhang->{'T'}\tG:$overhang->{'G'}\tC:$overhang->{'C'}";
   if ($coverage < $base_overlap){
     $stopoutput =  "$info => coverage too low\n";
     return $bin, "", $seq, 0, $totalclosed, 0,$stopoutput;
   }
   my ($ct_dna, $previous_bz) = (0, "");
   my $extend_nuc = "";
   #determine most likely extension
   BASE:
   foreach my $bz (sort {$overhang->{$b}<=>$overhang->{$a}} keys %$overhang){
      if($ct_dna == 1){## the two most abundant bases at that position
        my $bestcoverage = $overhang->{$previous_bz} + $overhang->{$bz};
        if($overhang->{$previous_bz} < $base_overlap){
          $stopoutput = "$info => coverage too low\n";
          return $bin, "", $seq, 0, $totalclosed, 0,$stopoutput;
        }
        elsif(($overhang->{$previous_bz} / $bestcoverage) >= $min_base_ratio){### a simple consensus btw top 2
          $extend_nuc = "$previous_bz";
          print CLOSEFILE "$info => extending with $previous_bz\n";
          last BASE;
        }
        else{
           $stopoutput = "$info => below ratio\n";
           return $bin, "", $seq, 0, $totalclosed, 0,$stopoutput;
         }
      }
      $previous_bz = $bz;
      $ct_dna++;
   }
   my $checkseq = $seq . $extend_nuc;
   $totalclosed++;
   #check if two sequences can merge, if so, check if they correspond with the estimated gap distance
   my ($newseq, $size) = MergeTwoSequences($checkseq, $seq2);
   if($size > 0){
     my $gaprangemin = $gapsize-$difference;
     my $gaprangemax = $gapsize+$difference;
     if(($totalclosed) >= $gaprangemin){
       $extended = 0;
       return $bin, $newseq, $checkseq, $size,$totalclosed, 0,$stopoutput;
     }else{
       print CLOSEFILE "\t\tsequences can be merged, but difference between gap and total closed\n";
     }
   }
   deleteData($bin, $subseq."$extend_nuc");
   $seq = $checkseq;
   $extended = 1;
   if(($totalclosed-$min_tig_overlap)>($gapsize+$difference)){
     $extended = 0;
   }
   return $bin, "", $seq, 0, $totalclosed, $extended,$stopoutput;
}

#try to merge two sequences by at least -n basepairs
sub MergeTwoSequences{
  my ($ctg1, $ctg2) = @_;

  my ($max_overlap, $newseq) = (100, "");
  while($max_overlap >= $min_tig_overlap){
    my $seq2 = $ctg2;
    my $seq1 = $ctg1;
    my $subseq2 = substr($ctg2,0,$max_overlap);
    my $subseq1 = substr($ctg1,-$max_overlap);
    if($subseq1 eq $subseq2){
      my $newctg1 = substr($ctg1,0,-$max_overlap);
      my $newctg2 = substr($ctg2,$max_overlap);
      $newseq = $newctg1.lc($subseq1).$newctg2;
      return ($newseq, length($subseq1));
    }
    $max_overlap--;
  }
  return ("", 0);
}

#try to merge two sequences by at least -n basepairs

sub MergeTwoSequencesFinal{
  my ($cont1, $cont2) = @_;
  my $max_trim = ($difference-$min_tig_overlap);
  my $trim = 1;
  my $suboverlap1 = substr($cont1,-$min_tig_overlap);
  my $suboverlap2 = substr($cont2,0,$min_tig_overlap);

  my $subcont1 = $cont1;
  $subcont1 = substr($cont1,-100) if(length($subcont1) > 100);
  my $res1 = index(substr($cont2,0,($min_tig_overlap+$difference)),$suboverlap1);
  my $res2 = rindex($subcont1,$suboverlap2);
  if($res1 > 0 || $res2 > 0){
    my ($max_overlap_start, $newseq) = ($res1, "");
    $max_overlap_start = (length($subcont1)-$res2) if((length($subcont1)-$res2) > $max_overlap_start);
    $max_overlap_start+=$min_tig_overlap;
    while($trim < $max_trim){
      my $max_overlap = $max_overlap_start;
      while($max_overlap >= $min_tig_overlap){
        my $subseq2 = substr($cont2,0,$max_overlap);
        my $subseq1 = substr($cont1,-$max_overlap);
        if($subseq1 eq $subseq2){
          my $newctg1 = substr($cont1,0,-$max_overlap);
          my $newctg2 = substr($cont2,$max_overlap);
          $newseq = $newctg1.lc($subseq1).$newctg2;
          return ($newseq, length($subseq1));
        }

        $max_overlap--;
      }
      $cont1 = substr($cont1,0,-1) if($res1 < 0);
      $cont2 = substr($cont2,1) if($res2 < 0);
      $trim++;
    }
  }
  return ("", 0);
}
###DELETE READ DATA IF IT HAS BEEN USED FOR FILLING A GAP
sub deleteData {
   my ($bin, $sequence) = @_;
   my $comp_seq = reverseComplement($sequence);
   delete $bin->{$sequence};
   delete $bin->{$comp_seq};
   return $bin;
}

#WRITE SUMMARY STATISTICS
sub writesummaryfiles2{
  my ($input_file, $sumfile) = @_;
  open (INFILE, $input_file) || die "Can't open input file $input_file.\n";

  my ($counter, $sum, $seq, $name, $foundN50, $sumN50, $totalNcount) = (0,0, "","", 0, 0);
  my (@line, @lengths);
  while (<INFILE>) {
    s/\r\n/\n/;
    chomp;
    $seq.= $_ if(eof(INFILE));
    if ($_ =~ /^[>]/ || eof(INFILE)) {
      if($counter > 0){
         push(@lengths, length($seq));
         $sum+= length($seq);
         my $Ncount = () = $seq =~ /[Nn]/g;
         $totalNcount += $Ncount;
         ($seq) = "";
      }
      $counter++;
    }
    else {
       $seq .= $_;
    }
  }
  $counter--;
  my $half_length = $sum/2;
  
  my @lengths2 = reverse sort { $a <=> $b } @lengths;
  
  for(my $i = 0; $i <= $#lengths && $foundN50 == 0; $i++)
  {
    $sumN50 += @lengths2[$i];
    if($sumN50 >= $half_length){
      $foundN50 = @lengths2[$i] if($sumN50 >= $half_length);
      last;
    }
  }
  $sumfile .= "\t\tTotal number of scaffolds = $counter\n";
  $sumfile .= "\t\tSum (bp) = ". $sum. "\n";
  $sumfile .= "\t\t\tTotal number of N's = $totalNcount\n";
  $sumfile .= "\t\t\tSum (bp) no N's = ". ($sum-$totalNcount)."\n";
  $sumfile .= "\t\tMax scaffold size = ". @lengths2[0]."\n";
  $sumfile .= "\t\tMin scaffold size = ". @lengths2[$#lengths]."\n";
  $sumfile .= "\t\tAverage scaffold size = ".int($sum/$counter)."\n";
  $sumfile .= "\t\tN50 = ". $foundN50. "\n\n";
  
  close (INFILE);
  close OUTFILE;

  return $sumfile;
}

###WRITE SUMMARY STATISTICS FOR ALL CONTIGS OR SCAFFOLDS
sub writesummaryfiles{
  my ($input_file, $sumfile) = @_;

  open (INFILE, $input_file) || die "Can't open input file $input_file.\n";

  my ($seq, $name,$counter,$sum,$totalNcount,$totallen, $totalGC, $totalgap) = ("","",0,0,0,0,0,0);
  my (@line, @lengths);
  while (<INFILE>) {
    chomp;
    $seq.=$_ if(eof(INFILE));
    if ($_ =~ /^[>]/ || eof(INFILE)) { 
      if($seq ne ""){
        $counter++;
         push(@lengths, length($seq));
         my $len = length($seq);
         my $Gcount = () = $seq =~ /G/g;
         $totalGC = $totalGC + $Gcount;
         my $Ccount = () = $seq =~ /C/g;
         $totalGC = $totalGC + $Ccount;
         $sum+= $len;
         my @gap = split(/N+/, $seq);
         $totalgap = $totalgap + $#gap;
         my $Ncount = () = $seq =~ /[Nn]/g;
         $totalNcount += $Ncount;
         $name = "";
         $seq = "";
      }
  
      $name = $_;
    }
    else {
       $seq .= $_;
    }               
  }
  
  my $half_length = $sum/2;
  my $N25 = $half_length/2;
  my $N75 = $half_length/2+$half_length;
  
  my @lengths2 = reverse sort { $a <=> $b } @lengths;

  my ($sumN50, $foundN50, $foundN25, $foundN75) = (0,0,0,0);
  for(my $i = 0; $i <= $#lengths; $i++)
  {
    $sumN50 += @lengths2[$i];
    $foundN50 = @lengths2[$i] if($sumN50 >= $half_length && $foundN50 == 0);
    $foundN25 = @lengths2[$i] if($sumN50 >= $N25 && $foundN25 == 0);
    $foundN75 = @lengths2[$i] if($sumN50 >= $N75 && $foundN75 == 0);
  }
  my $GCcontent = sprintf("%.2f", (($totalGC/($sum-$totalNcount))*100));

  $sumfile .= "\t\tTotal number of scaffolds = $counter\n";
  $sumfile .= "\t\tSum (bp) = ". $sum. "\n";
  $sumfile .= "\t\t\tTotal number of N's = $totalNcount\n";
  $sumfile .= "\t\t\tSum (bp) no N's = ". ($sum-$totalNcount)."\n";
  $sumfile .= "\t\tGC Content = $GCcontent\%\n";
  $sumfile .= "\t\tMax scaffold size = ". @lengths2[0]."\n";
  $sumfile .= "\t\tMin scaffold size = ". @lengths2[$#lengths]."\n";
  $sumfile .= "\t\tAverage scaffold size = ".int($sum/$counter)."\n";
  $sumfile .= "\t\tN25 = ". $foundN25. "\n";
  $sumfile .= "\t\tN50 = ". $foundN50. "\n";
  $sumfile .= "\t\tN75 = ". $foundN75. "\n\n";
  
  close (INFILE);
  close OUTFILE;

  return $sumfile, $totalNcount;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}

###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
}

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(SUMFILE), $| = 1)[0]);
  select((select(LOG), $| = 1)[0]);
  $|++;
}
###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGCatgc/TACGtacg/;
   return (reverse());
}

sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}
