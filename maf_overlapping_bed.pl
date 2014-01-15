#!/usr/bin/env perl

# takes a pairwise maf file and a bed file (cat two bed files, one for each species)
# outputs only those maf pairs where both overlap bed file

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($maf_file, $bed_file);
GetOptions (
  "m|maf:s" => \$maf_file,
  "b|bed:s" => \$bed_file,
);

################################################################################################
# load bed file into memory
# each chrontig is a hash element, intervals are pushed on in an array
# (intervals will always be in increasing order of start coord)

my %bed_intervals;

my $bed_file_fh = &read_fh ($bed_file);

while (<$bed_file_fh>) {
    push @{$bed_intervals{$1}}, [$2,$3] if /^(\S+)\t(\d+)\t(\d+)/;
}

for my $chrontig (keys %bed_intervals) {
    @{$bed_intervals{$chrontig}} = sort { $a->[0] <=> $b->[0] } @{$bed_intervals{$chrontig}};
}

################################################################################################
# open maf file
# read each three lines as one block.
# Store block coords
# Rapidly (binary?) search each chrontig for overlapping elements

my $maf_file_fh = &read_fh ($maf_file);

my ($found, $block, $chrontig, $strand, $start, $end);

print "##maf version=1\n";

while (<$maf_file_fh>) {
    next unless /^a score/;
    $found = 1;
    $block = $_;
    while ($_ = <$maf_file_fh>) {
        last if /^\s*$/;
        if (/^s\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\d+)/) {
            $block   .= $_;
            $chrontig = $1;
            $strand   = $4;
            if ($strand eq "+") {
                $start = $2;
                $end   = $2 + $3;
            }
            else {
                $start = $5 - $2 - $3;
                $end   = $5 - $2;
            }
            $found = 0 unless exists $bed_intervals{$chrontig} and scalar(&binary_search_intervals ( $bed_intervals{$chrontig}, $start, $end ));
            last unless $found;
        }
    }
    print $block . "\n" if $found;
}

print "\n";

################################################################################################
# subroutines
################################################################################################

sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}

sub binary_search_intervals
{
    my ($array_ref, $range_st, $range_en) = @_;
    my @sorted = @$array_ref;
    
    my $low_index = 0; my $high_index = $#sorted;
    my $found = 0;
    while ( $low_index <= $high_index and not $found )
    {
        my $mid_index = int ( ( $low_index + $high_index ) / 2 );
        if ( $sorted[$mid_index]->[0] <= $range_en and $sorted[$mid_index]->[1] >= $range_st )
        {
            $found = 1;
            return $mid_index + 1;
        }
        if ( $sorted[$mid_index]->[0] < $range_en )
        {
            $low_index = $mid_index + 1;
        } else
        {
            $high_index = $mid_index - 1;
        }
    }
    return $found;
}
