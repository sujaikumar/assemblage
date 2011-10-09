#!/usr/bin/env perl

=head1 NAME

fastx_separate.pl

=head1 SYNOPSIS

fastx_separate.pl -l [2|4] -p "pattern1" -p "pattern2" -p "pattern3" <STDIN>

=head1 DESCRIPTION

creates
  pattern1.fasta (or fastq if -l 4)
  pattern2.fasta
  pattern3.fasta
by pulling out all the reads that match a pattern and writing them to the output file

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2011.09.19

=cut

use strict;
use warnings;
use FileHandle;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $lines = 2; #default is fasta
my $ext = "fasta";
my @patterns = ();
GetOptions (
    "pattern=s" => \@patterns,
    "lines:i"   => \$lines,
);
if    ($lines==4) { $ext = "fastq" }

die 'Usage: fastx_separate.pl -l [2|4] -p "pattern1" -p "pattern2" -p "pattern3" <STDIN>
' unless @patterns and ($lines == 2 or $lines == 4);

my %pattern_fh;
for my $pattern (@patterns) {
    $pattern_fh{$pattern} = FileHandle->new(">$pattern.$ext");
}

my ($fh, $found);
while (<>) {
    $found = 0;
    for my $pattern (@patterns) {
        $fh = $pattern_fh {$pattern};
        if (/$pattern/) {
            $found = 1;
            print $fh $_;
            for my $i (2..$lines) {
                $_ = <>;
                print $fh $_;
            }
            last;
        }
    }
    next if $found;
    for my $i (2..$lines) {
        $_ = <>;
    }
}
