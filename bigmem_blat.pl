#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $outfile = pop @ARGV;
my $qryfile = pop @ARGV;
my $reffile = pop @ARGV;

my $num_cores = 4;
my $num_lines = 10000;

GetOptions (
  "cpus|cores:i" => \$num_cores,
  "lines:i"      => \$num_lines,
);

my $tmpdir =  "$outfile\_" . int(rand(10000));
mkdir $tmpdir or die $!;
my $parallel_command = "";
system ("fastaqual_multiline_to_singleline.pl $qryfile | split -l $num_lines -a 4 - $tmpdir/split_") and die $!;
system ("parallel -j $num_cores blat @ARGV $reffile {} {}.psl ::: $tmpdir/split_*") and die $!;
system ("cat $tmpdir/*.psl >$outfile; rm -rf $tmpdir") and die $!; 