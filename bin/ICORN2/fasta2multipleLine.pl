#!/usr/bin/perl

# Martin Hunt 11/03/09
# mh12@sanger.ac.uk

# takes a one-line-per-sequence fasta and converts to the evil multi-line format

use strict;
use warnings;

if ($#ARGV < 1){
    print "usage: fasta2multipleLine.pl <input fasta file> <output fasta file> <line width>\n";
    exit;
}

my $infile     = $ARGV[0];
my $outfile    = $ARGV[1];
my $line_width = $ARGV[2];

my $out_string = "";
my $tmp_string;

open F, "<", $infile or die "Error opening $infile";

while (my $line = <F>){
    if ($line =~ /^>/ || length($line) <= $line_width) {
        $out_string .= $line;
    }
    else {
        while (1){
            $out_string .= substr($line,0,$line_width) . "\n";
            $line = substr($line, $line_width);
            if (length($line) < $line_width){
                 $out_string .= "$line";
                 last;
            }
        }
    }
}

close F;

$out_string =~ s/\n+/\n/g;

open F, ">", $outfile or die "Error opening $outfile";
print F $out_string;
close F;

