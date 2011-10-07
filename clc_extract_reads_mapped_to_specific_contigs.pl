#!/usr/bin/env perl

=head1 NAME

clc_extract_reads_mapped_to_specific_contigs.pl

=head1 SYNOPSIS

clc_extract_reads_mapped_to_specific_contigs.pl -c inputcasfile -id contigids1.txt -id contigids2.txt [-u]

=head1 DESCRIPTION

- Takes a clc .cas file as input
  (the reference contigs file referred to inside .cas has to be available on disk for the prog to work)
- Takes multiple files with contigids
- Returns the readnums that mapped to each contigid file (one readnum file per contigidfile)
- Readnums are returned in pairs. If a /1 read maps to one file and /2 to another file, then put both in /1's file
- If the same contigid if present in two files, both readnum output files will contain that readnum
  (i.e. readnums will be duplicated. so it is the user's responsibility to ensure contig ids are not duplicated) 
- -u will also print unmapped.readnums (reads that were either not mapped to ANY contig, or reads that were not mapped to
  any of the contigid lists

Needs:

assembly_info
assembly_table (both are in clc ngs cell suite - do not need license to run)

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2011.05.22

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use IO::File;

my $casfile;
my @contigidfiles;
my $fastafile; # reference fasta file taken from cas file or specified on the command line
my $unmapped = "";

GetOptions (
    "casfile=s"    => \$casfile,
    "idfiles=s{,}" => \@contigidfiles,
    "fastafile=s"  => \$fastafile,
    "unmapped"     => \$unmapped,
);

#-------------------------------

my @readfiles;

open CAS, "assembly_info $casfile |" or die $!;
if (not defined $fastafile) {
    while (<CAS>) {
        /Contig files:/ and ($_ = <CAS>) and /(\S+)/ and $fastafile = $1 and last
    }
}

while (<CAS>) {
    next unless /Read files:/;
    while (<CAS>) {
        last if /^\s*$/;
        /^\s*(\S+)\s+\[\s+(\d+)/ and push (@readfiles, $1);
    }
    last if /^\s*$/;
}
close CAS;

#-------------------------------
# load fasta contig file into memory by id (and store its seq number starting with 0)

my %contigid_num;
my $num = -1;
my $fh = &read_fh($fastafile);
while (<$fh>) {
    if (/^>(\S+)/) {
        $num++;
        $contigid_num{$1} = $num;
    }
}
close $fh;

#-------------------------------
# load each contigidfile (that we want to extract reads for) into memory
# create a filehandle to write to for each contigidfile called contigidfile.readnums

my %include_ids; # two level hash 
my %contigidreadnums_fh;

for my $contigidfile (@contigidfiles) {
    $contigidreadnums_fh{$contigidfile} = IO::File->new;
    open $contigidreadnums_fh{$contigidfile}, ">$contigidfile.mapped.readnums" or die $!;
    open ID, "<$contigidfile" or die $!;
    while (<ID>) {
        /^>?(\S+)/ and $include_ids{$contigidfile}{$contigid_num{$1}} = 1;
    }
}
$contigidreadnums_fh{unmapped} = IO::File->new;
if ($unmapped) {
    open $contigidreadnums_fh{unmapped}, ">unmapped.readnums" or die $!;
}

#-------------------------------
# Now, go through each line in the output of assembly_table .cas
# and each time you see a contig id, check which of the contigidfiles it belongs in,
# and print the readnum to that file

open CAS, "assembly_table $casfile |" or die $!;
my ($fr, $fc, $rr, $rc);
while (<CAS>) {
    # get forward read, forward contig, reverse read, reverse contig;
    if (/^(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/) {
        $fr = $1 + 1;
        $fc = $5;
    }
    $_ = <CAS>;
    if (/^(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/) {
        $rr = $1 + 1 ;
        $rc = $5;
    }
    if ($fc == -1 and $rc == -1 and $unmapped) {
        print {$contigidreadnums_fh{unmapped}} "$fr\n$rr\n";
        next;
    }

    $fc = $rc if $fc == -1;
    my $found = 0;
    for my $contigidfile (@contigidfiles) {
        if ($include_ids{$contigidfile}{$fc}) {
            print {$contigidreadnums_fh{$contigidfile}} "$fr\n$rr\n";
            $found = 1;
        }
    }
   if (not $found and $unmapped) { 
       print {$contigidreadnums_fh{unmapped}} "$fr\n$rr\n";
   }
}
close CAS;

#############################################################################

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

