#! /usr/bin/perl -w
#
# File: snp.correctString.pl
# Time-stamp: <03-Jul-2009 13:30:59 tdo>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#icornLib::
# Author: Thomas Dan Otto
#
# Description:
#
# tdo 21.01: include read pairs of different length for SNPoMatic
#

use Data::Dumper;
use strict;
use warnings;
use lib "$ENV{ICORN_HOME}";
use icornLib;
use POSIX;



my $genomeName = shift;
my $vcfFile    = shift;
my $readRoot   = shift;
my $readRoot1FragmentSize = shift;
my $resultName = shift;

my $minQual    = shift;
if (!defined($resultName)) {
  print "Usage: <genome> <VCFfile - from GATK> <root fastq file> <expected insert size> <Resultname> optional: <minQual default 60>\n";
	exit (-1);
}

if ( ! (-f $genomeName)){
	die " Sorry, genome doesn't exist\n";	
}

if (-f $readRoot.".SNPoMATIC_1.fastq"){
$readRoot=$readRoot.".SNPoMATIC";
}
if ( ! (-f $readRoot."_1.fastq" ) || ! (-f $readRoot."_2.fastq" )){
	die " Sorry, problem accession the fastq files\n";	
}

if (!(isdigit($readRoot1FragmentSize)) || $readRoot1FragmentSize <0 || $readRoot1FragmentSize > 20000){
	die "Please set a better expected inset size.\n";	
}

if (!defined($minQual)) {
  $minQual=60;
}
## Variable for version 2, including the shift and control with SNPoMatic

my $pileupRef = $genomeName.".ref.fa";
my $pileupQry = $genomeName.".qry.fa";



my ($ref_snp,$ref_del,$ref_ins,$ref_coverage)   = icornLib::getVCF($vcfFile,$minQual);

#get reference fasta
my $ref_str                   = icornLib::getFasta($genomeName);


  # Do the irstcorrection
  my ($ref_stats,$ref_shiftCoordinates);

  ($ref_str,$ref_stats)
	= icornLib::correctSequenceV2(0,$ref_str,
									 $ref_snp,
									 $ref_del,
									 $ref_ins);
									 

  icornLib::writeString($ref_str,"$resultName.tmp");

  # Generate SNPoMatic pileups of of both sequences
  icornLib::callSubSNPoMaticPairedV2($genomeName,$pileupRef,$readRoot,$readRoot1FragmentSize);   # in itereation might
                                              # not be needed, but
                                              # more robuts
  icornLib::callSubSNPoMaticPairedV2("$resultName.tmp",$pileupQry,$readRoot,$readRoot1FragmentSize);
  
  #reload the reference string, as it was changed:
  my $ref_str2   = icornLib::getFasta($genomeName);
  
  my $ref_pileupRef = icornLib::getPileupSNPoMaticV2($pileupRef,$ref_coverage);
  
  my $ref_pileupQry = icornLib::getPileupSNPoMaticV2($pileupQry,$ref_coverage);
  
  ## later the shift, as just for stats
 # my $ref_sizeChr   = icornLib::getSizeChr($ref_pileupRef);
  
  ($ref_str2,$ref_stats)    = icornLib::correctSequenceV2(2,$ref_str2,
												 $ref_snp,
												 $ref_del,
												 $ref_ins,
												 $ref_pileupRef,
												 $ref_pileupQry
												 );

  icornLib::writeString($ref_str2,$resultName);

  icornLib::writeStatsV2($ref_stats,$resultName,$ref_coverage,$ref_str2,
							$ref_pileupRef,$ref_pileupQry);

exit 0;
