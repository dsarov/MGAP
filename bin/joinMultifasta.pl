#! /usr/bin/perl -w
#
# File: joinMultifasta.pl
# Time-stamp: <26-Sep-2011 15:20:28 tdo>
# $Id: $
#
# Copyright (C) 2012 by Pathogene Group, Wellcome Trust Sanger Institute
#
# Author: Thomas Dan Otto
#
# Description: This program will join a multifasta file to a union file for ABACAS 1. It will put a spacer (1000n500g1000n), between each replicon. 
#              After the run of ABACAS, the program split2Multifasta.pl will split query relative to the multifasta reference
# important, repass the file $result.splitinfo for the 

use strict;


if (scalar(@ARGV)<1) {
  die "usage: <Genome (multi-fasta)> <resultName>\n\nThis script helps to run ABACAS with multifasta reference. After running ABACAS use split2Multifasta.pl to split the query into the correct amount of sequences.\nPlease the file <resultname>.splitinfo.txt as it is required to split the query later.\n";
}

my $spacer;
my $fastaFile = shift;
my $resultName = shift;

for (1..100){ $spacer.="NNNNNNNNNN" }
$spacer.="\n";
for (1..500) { $spacer.="GGGGGGGGGG" }
$spacer.="\n";
for (1..100){ $spacer.="NNNNNNNNNN" }
$spacer.="\n";

my ($res, $splitinfo) =joinFasta($fastaFile,$spacer);

open F, "> $resultName" or die "Problem to write to file $resultName: $!\n";
print F $res;
close(F);
open F, "> $resultName.splitinfo.txt" or die "Problem to write to file $resultName: $!\n";
print F $splitinfo;
close(F);

sub joinFasta{
  my $file   = shift;
  my $spacer = shift;
 
  open F2, $file or die "Sorry, I couldn't find $file: $!  \n";
  
  my $splitinfo;
  
  my($name) = $file =~ /(\w+)/;	
  my $res=">union.$name\n";
  my $n;
  ### ignore the first description line
  $_=<F2>;
  while (!(/^>/)){
  	 $_=<F2>;	
  }
  ### get the first name
  /^>(\S+)/;
  $splitinfo.="$1\t1\t";
  
  my $length=0;
  while (<F2>) {
        chomp;
        if (/^>(\S+)/){
    		$splitinfo.=$length."\n$1\t";
        	$res.=$spacer;
        	$length+=length($spacer);
        	$splitinfo.=$length."\t";
        }
        else {
	    $length+=length($_);
		$res.=$_
       }
  }
  $splitinfo.=$length."\n";
  close(F2);
  return ($res,$splitinfo);
}

