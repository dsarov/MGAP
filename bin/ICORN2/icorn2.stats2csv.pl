#! /usr/bin/perl -w
#
# File: icorn2.stats2csv.pl
# Time-stamp: <08-Oct-2012 14:47:19 tdo>
# $Id: $
#
# Copyright (C) 2012 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#

while (<STDIN>) {

for my $i (1..7) {
  $_=<STDIN>;
  @ar=split(/\s+/);
  print "\t$ar[1]";  
}
print "\n";  
}
