#! /bin/bash
#
# File: bam.correctLineLength.pl
# Time-stamp: <13-Mar-2012 15:37:48 tdo>
# $Id: $
#
# Copyright (C) 2010 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
# 



name=$1

if [ -z "$name" ] ; then
	echo "Please give the (multi-) fasta to be corrected...";
exit;
fi

tmp=$$


fasta2singleLine.pl $name $name.$tmp
fasta2multipleLine.pl  $name.$tmp $name 80
rm $name.$tmp
