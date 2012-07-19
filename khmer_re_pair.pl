#!/usr/bin/env perl

=head1 NAME

khmer_re_pair.pl

=head1 SYNOPSIS

khmer_re_pair.pl -o original_interleaved_fasta -k khmer_filtered_fasta

=head1 DESCRIPTION

- Takes -o, an original paired (interleaved) fasta read file (one line per sequence)
- Takes the khmer_filtered_fasta file
- For each read in the khmer_filtered_fasta file, it prints out both the /1 and the /2 (without duplicates)
- The goal is to not lose pairing info which khmer does lose
- Assumes the khmer_filtered_fasta file was created from the original file, and is in the same order

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2012.07.08

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $original_file;
my $khmer_file;

GetOptions (
    "original=s"=> \$original_file,
    "khmer=s"   => \$khmer_file,
);

if (not $original_file and not $khmer_file) {
    print STDERR "Usage: khmer_re_pair.pl -o original_interleaved_fasta -k khmer_filtered_fasta\n";
    exit;
}

my $original_fh = &read_fh ($original_file);
my $khmer_fh    = &read_fh ($khmer_file);

my %khmer_hash;
my $pair = "";

#-------------------------------------------
# load khmer filtered file headers in memory

while (<$khmer_fh>) {
    die "-khmer file does not seem to have /1 or /2 as the first read's suffix" unless /^(>\S+?)\/\d\b/;
    $khmer_hash{$1} = 1;
    <$khmer_fh>; #ignore the sequence 
}

#-------------------------------------------
# go through original file and print out 
# all pairs where header is in khmer_hash

while (<$original_fh>) {
    die "-original file does not seem to have /1 or /2 as the first read's suffix" unless /^(>\S+?)\/\d\b/;
    $pair = $_;
    for (1..3) {$pair .= <$original_fh>}
    print $pair if $khmer_hash{$1};
}

########################################

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
