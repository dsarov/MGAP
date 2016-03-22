#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Inline 'C';

my ($quality_threshold, $min_length, $window, $offset, $paired, $fileA, $fileB, $basename) = (10, 20, 0, 64, 0, "","","paired");
GetOptions (
  "quality:i" => \$quality_threshold,
  "min:i" => \$min_length,
  "window:i" => \$window,
  "offset:i" => \$offset,
  "file1:s" => \$fileA,
  "file2:s" => \$fileB,
  "base:s" => \$basename
);

my $printusage = "USAGE:
  -o (ascii offset, default 64 for Illumina FastQ files, set to 33 for Sanger FastQ files),
  -q (quality threshold, default 10),
  -m (minimum contiguous positions above the quality, default 20),
  -w (window size - each position's average quality is calculated based on the qualities of this many bases on either side, default 0)
  -file1 First read-file of the paired-reads in fastQ format
  -file2 Second read-file of the paired-reads in fastQ format
  -b base name for the output (e.g. -b pair will generate two files: pair_1.fastq and pair_2.fastq\n";

if ($fileA eq "" || $fileB eq "" || !(-e $fileA) || !(-e $fileB)){
  die $printusage."\nError: No correct input files specified\n";
}
open FILEA, "< $fileA" || die "$printusage\n\ncan't open $fileA\n";
open FILEB, "< $fileB" || die "$printusage\n\ncan't open $fileB\n";

my $outputfile1 = "$basename"."_1.fastq";
my $outputfile2 = "$basename"."_2.fastq";
open OUT1, "> $outputfile1" || die "$printusage\nError: can't write to $outputfile1\n";
open OUT2, "> $outputfile2" || die "$printusage\nError: can't write to $outputfile1\n";


my $toprint;
my $numpairedreads=0;
my $filtpairedreads = 0;
my $numnucleotides1=0;
my $filtnucleotides1= 0;
my $numnucleotides2=0;
my $filtnucleotides2= 0;
while (<FILEA>) {
	$numpairedreads++;
        my $seqheader = $_;
	my $seqstring = <FILEA>;
	my $qualheader = <FILEA>;
	die "$printusage\n$fileA does not seem to be a fastq file = $qualheader\n" unless $qualheader =~ /^\+/;
	chomp(my $qualstring = <FILEA>);

	my ($hqrun_start, $hqrun_length) = (0,0);
	maskqual($qualstring, $offset, $window, $quality_threshold, $hqrun_start, $hqrun_length);
	if ($hqrun_length < $min_length) {
	  for (1..4) { <FILEB> }
	  $toprint = "";
	  $filtpairedreads++;
	} else {
		$toprint = "$seqheader" . substr($seqstring,$hqrun_start,$hqrun_length) .
			"\n$qualheader" . substr($qualstring,$hqrun_start,$hqrun_length) . "\n";
                my $lenorig1 = length($seqstring)-1;
	        my $lenfilt1 = $hqrun_length;

		my $seqheader = <FILEB>;
	        my $seqstring = <FILEB>;
        	my $qualheader = <FILEB>;
        	die "$printusage\n$fileB does not seem to be a fastq file\n" unless $qualheader =~ /^\+/;
        	chomp(my $qualstring = <FILEB>);
        
        	my ($hqrun_start, $hqrun_length) = (0,0);
        	maskqual($qualstring, $offset, $window, $quality_threshold, $hqrun_start, $hqrun_length);
                if ($hqrun_length >= $min_length) {
                  my $lenorig2 = length($seqstring)-1;
	          my $lenfilt2 = $hqrun_length;
                  $numnucleotides1+=$lenorig1;
                  $filtnucleotides1+=$lenfilt1;
                  $numnucleotides2+=$lenorig2;
                  $filtnucleotides2+=$lenfilt2;

                  print OUT1 $toprint;
                  print OUT2 "$seqheader" . substr($seqstring,$hqrun_start,$hqrun_length) .
			"\n$qualheader" . substr($qualstring,$hqrun_start,$hqrun_length) . "\n";
                  
                  $toprint = "";
                }else{
                  $filtpairedreads++;
                }

	}
}
close OUT1;
close OUT2;

my $remainingreads = $numpairedreads - $filtpairedreads;
my $percfiltreads = sprintf("%.2f",($filtpairedreads/$numpairedreads)*100);
print "\ntotal number of reads = $numpairedreads\n";
print "total number reads remaining = $remainingreads\n";
print "total number reads filtered = $filtpairedreads\n";
print "\tpercentage filtered = $percfiltreads%\n\n";

my $percfilt1 = sprintf("%.2f",100-($filtnucleotides1/$numnucleotides1)*100);
my $percfilt2 = sprintf("%.2f",100-($filtnucleotides2/$numnucleotides2)*100);
print "total number of nucleotides read 1 = $numnucleotides1\n";
print "total number of nucleotides read 1 after trimming = $filtnucleotides1\n\tpercentage trimmed= $percfilt1%\n\n";
print "total number of nucleotides read 2 = $numnucleotides2\n";
print "total number of nucleotides read 2 after trimming = $filtnucleotides2\n\tpercentage trimmed= $percfilt2%\n\n";





__END__
__C__
void maskqual(char* qstring, int o, int w, int q, SV* hqrun_start, SV* hqrun_length) {
	
	int i;
	char qplus[strlen(qstring) + 2*w];

	for (i = 0; i < w; i++) {
		qplus[i] = q + o; /* won't affect quality sum */
	}
	for (i = w; i < strlen(qstring) + w; i++) {
		qplus[i] = qstring[i-w]; /* copy quality values */
	}
	for (i = strlen(qstring) + w; i < strlen(qstring) + 2*w; i++) {
		qplus[i] = q + o; /* won't affect quality sum */
	}

	int curr_hqrun_start = 0;
	int longest_hqrun_start = 0;
	int curr_hqrun_length = 0;
	int longest_hqrun_length = 0;

	for (i = w; i < strlen(qstring) + w; i++) {
		int sum = 0;
		int j;
		for (j = -w; j <= w; j++) {
			sum = sum + qplus[i+j] - o;
		}
		/* printf("%d ",sum/(2*w+1)); */
		if (sum >= q * (2*w+1)) {
			curr_hqrun_length++;
			if (curr_hqrun_length > longest_hqrun_length) {
				longest_hqrun_length = curr_hqrun_length;
				longest_hqrun_start = curr_hqrun_start;
			} 
		} else {
			curr_hqrun_length = 0;
			curr_hqrun_start = i-w+1;
		}
	}
	sv_setiv(hqrun_start,longest_hqrun_start);
	sv_setiv(hqrun_length,longest_hqrun_length);
}
