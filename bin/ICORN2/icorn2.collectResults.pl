#!/usr/bin/env perl

use strict;
use warnings;


if ($#ARGV != 0) {
    print STDERR "Usage: $0 ICORN_output_directory

Gathers stats for each iteration of ICORN.
Output is to terminal in tab-separated format
";
    exit(1);
}

my $root_dir = $ARGV[0];

# Use the listing of icorn2 correct*.o files to get how many iterations
my @files = glob("$root_dir/out.ICORN2.CORRECT.*.o");
my $total_iterations;
my $o_files_prefix;

if ($files[-1] =~ /(.*)\.(\d+)\.o$/) {
    $total_iterations = $2;
    $o_files_prefix = $1;
}
else {
    print STDERR "Error getting iteration number from filename '$files[-1]'\n";
    exit(1);
}


print "#iter\tperfect1\tperfect1%\tperfect2\tperfect2%\tSNP\tINS\tDEL\tHETERO\tRej.SNP\tRej.INS\tRej.DEL\n";

# gather stats on each iteration
for my $i (1..$total_iterations) {

    # get perfect mapped reads stats
    my $fname = "$o_files_prefix.$i.o";
    open F, $fname or die "Error opening '$fname'";
    my $perfect_1 = -1;
    my $perfect_1_pc;
    my $perfect_2;
    my $perfect_2_pc;
    while (<F>) {
        # Want to match a line like this:
        #Reading solexa pair data from x_1.fastq ... scanned 33000000 solexa reads, 17256697 (52.29%) matched , 17710665 total matches, 0 skipped.
        if (/^Reading solexa pair data from .*solexa reads, (\d+) \((\d+\.\d+)\%\)/) {
            if ($perfect_1 == -1) {
                $perfect_1 = $1;
                $perfect_1_pc = $2;
            }
            else {
                $perfect_2 = $1;
                $perfect_2_pc = $2;
            }
        }
    }

    close F;

    # get SNP/indel stats
    my @files = glob("$root_dir/ICORN2_$i/*.General.stats");
    if ($#files != 0) {
        print STDERR "Too many files! $root_dir/ICORN2_$i/*.General.stats\n";
        exit(1);
    }
    $fname = $files[0];
    my $type;
    my %counts;
    open F, $fname or die "Error opening '$fname'\n";
    while (my $line = <F>) {
        chomp $line;
        if ($line !~/\t/){
            $type = $line;
            $counts{$type} = 0;
        }
        else {
            my (undef, $n) = split(/\t/, $line);
            $counts{$type} += $n;
        }
    }

    close F;

    # print the results for this iteration
    print join("\t", $i,
        $perfect_1, $perfect_1_pc, $perfect_2, $perfect_2_pc,
        $counts{SNP}, $counts{INS}, $counts{DEL}, $counts{HETERO},
        $counts{'Rej.SNP'}, $counts{'Rej.INS'}, $counts{'Rej.DEL'}) . "\n";
}

