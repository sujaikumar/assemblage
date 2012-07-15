#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max sum);

=head1 NAME

clc_len_cov_gc_insert.pl

=head1 SYNOPSIS

1. Create a single end assembly (clc_novo_assemble is fast) - nucleotide fasta file

   clc_novo_assemble -o clcse.fna -q libA_1.txt.gz libA_2.txt.gz 

2. Use clc_ref_assemble to align reads to this assembly (interleave reads using -i)

   clc_ref_assemble  -d clcse.fna -q -i libA_1.txt.gz libA_2.txt.gz -o clcse.fna.cas
   
3. Use this script to generate len cov gc files and insert estimate files (that you can then plot using R ggplot2)

   clc_len_cov_gc_insert.pl -c clcse.fna.cas

Requirements:

1. clc assembly cell 4.03 beta or greater, in your path

=head1 OPTIONS

-c <clc_ref_assemble_outut_in_CAS_format>

Required: Note: The cas file stores relative paths to the reference seq and reads.
This script expects those files to be in the same place.
You can get around this by symlinking file locations or
changing the file locations stored in the CAS file using CLC's
change_assembly_files utility

-i
Optional: if you want to generate insert size estimates
If you use this option, then this script assumes that interleaved reads were aligned
i.e. you used
   clc_ref_assemble -q -i read1_file read2_file ...
or clc_ref_assemble -q interleaved_read_file    ...
If you haven't done this, don't use this option.

-r <integer> <integer> ...
-r 0 50 100 150 200
Specifies range breaks for insert size stats file.
Does nothing if you don't use the -i option.
Will summarise number of pairs with insert size between 0 and 50, 50 and 100, etc
If you expect 300 bp inserts, then you might choose:
-r 0 200 400 1000
to see how many reads were between 0 and 200 (i.e too short),
how many were between 200 and 400 (just right),
how many were between 400 and 1000 (too long)
The same information can be extracted from the freq.txt file.

=cut

my ($casfile, $fastafile, $statsfile, $contigfile, $outfile, $libname);
my $delimiter = "_";
my @range_breaks = (0,100,200,300,400,500,1000,2000,5000,10000);
my $insertestimate = 0;
GetOptions (
    "casfile=s"  => \$casfile,
    "delimiter=s"=> \$delimiter,
    "fastafile:s"=> \$fastafile,
    "statsfile:s"=> \$statsfile,
    "range=i{,}" => \@range_breaks,
    "insert" => \$insertestimate,
    "outfile=s" => \$outfile,
    "libname=s" => \$libname,
);

$outfile   = $casfile unless $outfile;
$libname   = $casfile unless $libname;
$fastafile = $outfile . ".lencovgc.fna" unless $fastafile;
$statsfile = $outfile . ".lencovgc.txt" unless $statsfile;

#############################################################################

open CAS, "assembly_info $casfile |" or die $!;
while (<CAS>) { /Contig files:/ and ($_ = <CAS>) and /(\S+)/ and $contigfile =$1 and last }

my $fastaidhash = &fastafile2hash_id ($contigfile);

open  FASTA, ">$fastafile" or die $!;
open  STATS, ">$statsfile" or di$!;
print STATS  "read_set\tcontig_id\tlen\tcov\tgc\n";

my ($header,$desc,$length,$gccount,$nonatgc,$seq);
while (<CAS>) {
    next unless /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*$/;
    my $header  = $$fastaidhash{$1}{header};
    my $desc    = $$fastaidhash{$1}{desc};
    my $length  = $$fastaidhash{$1}{len};
    my $gccount = $$fastaidhash{$1}{gc};
    my $nonatgc = $$fastaidhash{$1}{nonatgc};
    my $seq     = $$fastaidhash{$1}{seq};
    print FASTA ">".$header.$delimiter.$length.$delimiter.$4.$delimiter.sprintf("%.2f",$gccount/($length-$nonatgc))."$desc\n$seq\n";
    print STATS "$libname\t$header\t$length\t$4\t".($gccount/($length-$nonatgc))."\n";
}
close CAS;
close FASTA;
close STATS;

#############################################################################

exit 0 unless $insertestimate;

open CAS, "assembly_table $casfile |" or die $!;
my (%read1, %read2, @F, $numpairs_unmapped, $numpairs_one_mapped, $numpairs_diff_mapped, $numpairs_bothmapped, %FR, %RF, $FF, $totalpairs);
while (<CAS>) {
    chomp;
    @F=split/\s+/;
    $read1{rstart} =  $F[2];
    $read1{rend}   =  $F[3];
    $read1{contig} =  $F[4];
    $read1{tstart} =  $F[5];
    $read1{tend}   =  $F[6];
    $read1{strand} = ($F[7] ? 0 : 1);
    $read1{unique} = ($F[8] == 1 ? 1 : 0);
    $_ = <CAS>;
    chomp;
    @F=split/\s+/;
    $read2{rstart} =  $F[2];
    $read2{rend}   =  $F[3];
    $read2{contig} =  $F[4];
    $read2{tstart} =  $F[5];
    $read2{tend}   =  $F[6];
    $read2{strand} = ($F[7] ? 0 : 1);
    $read2{unique} = ($F[8] == 1 ? 1 : 0);
    $totalpairs++;
    
    if ($read1{contig} eq $read2{contig} and $read1{contig} != -1) {
        $numpairs_bothmapped++;
        if ($read1{strand} == $read2{strand}) {
            $FF++
        }
        elsif     ( $read1{strand} ) {
            if    ( $read2{tend}   >= $read1{tstart} ) { $FR{ $read2{tend} - $read1{tstart} }++ }
            elsif ( $read2{tstart} <= $read1{tend} )   { $RF{ $read1{tend} - $read2{tstart} }++ }
        }
        elsif     ( $read2{strand} ) {
            if    ( $read1{tend}   >= $read2{tstart} ) { $FR{ $read1{tend} - $read2{tstart} }++ }
            elsif ( $read1{tstart} <= $read2{tend} )   { $RF{ $read2{tend} - $read1{tstart} }++ }
        }
        next;
    }
    elsif ($read1{contig} != -1 and $read2{contig} != -1 and $read1{contig} ne $read2{contig}) {
        $numpairs_diff_mapped++
    }
    elsif ($read1{contig} == -1 and $read2{contig} == -1) {  $numpairs_unmapped++}
    elsif ($read1{contig} == -1 xor $read2{contig} == -1) {$numpairs_one_mapped++};
}
close CAS;

open FREQ, ">$outfile.insert.freq.txt" or die $!;
open STAT, ">$outfile.insert.stat.txt" or die $!;

$numpairs_unmapped=0    if not defined $numpairs_unmapped;
$numpairs_one_mapped=0  if not defined $numpairs_one_mapped;
$numpairs_diff_mapped=0 if not defined $numpairs_diff_mapped;
$numpairs_bothmapped=0  if not defined $numpairs_bothmapped;
$FR{0} = 0 if scalar keys %FR == 0;
$RF{0} = 0 if scalar keys %RF == 0;
$FF    = 0 if not defined $FF;

foreach (sort { $a <=> $b } keys %FR) {
	print FREQ "$libname\tFR\t$_\t$FR{$_}\n"
}
foreach (sort { $a <=> $b } keys %RF) {
	print FREQ "$libname\tRF\t$_\t$RF{$_}\n"
}
print STAT "Total pairs                                       : $totalpairs
Pairs with neither read mapped                    : $numpairs_unmapped (". sprintf("%.1f",100*$numpairs_unmapped/$totalpairs). " %)
Pairs with one read mapped                        : $numpairs_one_mapped (". sprintf("%.1f",100*$numpairs_one_mapped/$totalpairs). " %)
Pairs with both reads mapping different contigs   : $numpairs_diff_mapped (". sprintf("%.1f",100*$numpairs_diff_mapped/$totalpairs). " %)
Pairs with both reads best mapping to same contig : $numpairs_bothmapped (". sprintf("%.1f",100*$numpairs_bothmapped/$totalpairs). " %)\n";
print STAT "Num of FR    pairs: " . sum(values %FR) . " (".sprintf("%.1f",100 * sum(values %FR)/$numpairs_bothmapped)." % of $numpairs_bothmapped pairs mapping to same contig)\n";
my %range_counts;
for my $ins (keys %FR) {
	for my $range_index (1..$#range_breaks) {
		if ($ins >= $range_breaks[$range_index -1] and $ins < $range_breaks[$range_index]) {
			$range_counts{FR}[$range_index] += $FR{$ins};
			last;
		}
	}
}
for my $range_index (1..$#range_breaks) {
	print STAT "    FR Pairs with insert size range ".$range_breaks[$range_index -1]."-".$range_breaks[$range_index] ." : ". (defined $range_counts{FR}[$range_index] ? $range_counts{FR}[$range_index] : 0) . "\n";
}
print STAT "Num of RF    pairs: " . sum(values %RF) . " (".sprintf("%.1f",100 * sum(values %RF)/$numpairs_bothmapped)." % of $numpairs_bothmapped pairs mapping to same contig)\n";
for my $ins (keys %RF) {
	for my $range_index (1..$#range_breaks) {
		if ($ins >= $range_breaks[$range_index -1] and $ins < $range_breaks[$range_index]) {
			$range_counts{RF}[$range_index] += $RF{$ins};
			last;
		}
	}
}
for my $range_index (1..$#range_breaks) {
	print STAT "    RF Pairs with insert size range ".$range_breaks[$range_index -1]."-".$range_breaks[$range_index] ." : ". (defined $range_counts{RF}[$range_index] ? $range_counts{RF}[$range_index] : 0) . "\n";
}

print STAT "Num of FF/RR pairs: " . $FF . "\n";

close FREQ;
close STAT;

exit 0;

#############################################################################

sub fastafile2hash_id {
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $header;
    my $id = 0;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $id++;
            $sequences{$id}{header} = $1;
            $sequences{$id}{desc}   = $2;
        }
        else {
            chomp $line;
            $sequences{$id}{seq}     .= $line;
            $sequences{$id}{len}     += length $line;
            $sequences{$id}{gc}      += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences{$id}{nonatgc} += ($line =~ tr/N/N/);
        }
    }
    close $fh;
    return \%sequences;
}

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
