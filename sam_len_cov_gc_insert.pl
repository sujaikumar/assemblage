#!/usr/bin/env perl

=head1 NAME

sam_len_cov_gc_insert.pl

=head1 SYNOPSIS

1.  Create a single end assembly (soap de novo, velvet, abyss etc should all work) - nucleotide fasta file
    eg:

    velveth 51 k51 -short -fastq.gz libA_1.txt.gz libA_2.txt.gz
    velvetg 51

    Note: If you already have a paired end assembly, then that's fine too. The advantage of a single end
    assembly is that you minimise the chances of assembly artifacts at this stage
   
2.  Use any short read aligner that gives SAM or BAM output.
   
3.  Use this script to generate len cov gc files and insert estimate files (that you can then plot using R ggplot2)

    sam_len_cov_gc_insert.pl -f contigfastafile -s samfile

    The SAMfile can be created on the fly from a BAMfile like this:

    sam_len_cov_gc_insert.pl -f contigfastafile -s <(samtools view bamfile)
    or
    samtools view bamfile | sam_len_cov_gc_insert.pl -f contigfastafile -s -
    
    Note: If you want to also generate accurate insert size estimates, make sure the SAM/BAM file is created with
    interleaved reads and provide the -i flag.

    sam_len_cov_gc_insert.pl -f contigfastafile -s samfile -i

=head1 OPTIONS

-s <samfile>

-o outputprefix
Optional: By default, all output files are created using the input SAMfile name as prefix,
or "out" if the SAMfile is provided as an iostream. Using this option explicitly sets the
output prefix

-i
Optional: if you want to generate insert size estimates
If you use this option, then this script assumes that interleaved reads were aligned
If the reads weren't interleaved while mapping, this option will return meaningless results

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

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2012.03.30

=cut

use strict;
use warnings;
use List::Util qw(min max sum);
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($contigfastafile, $outfile, $samfile, $libname, $fastafile, $lencovgcfile);
my $delimiter = "_";
my @range_breaks;
my $insertestimate = 0;
GetOptions (
    "samfile=s"  => \$samfile,
    "delimiter=s"=> \$delimiter,
    "f|c|contigfastafile:s"=> \$contigfastafile,
    "range=i{,}" => \@range_breaks,
    "i|insert" => \$insertestimate,
    "o|outfile=s" => \$outfile,
    "l|libname=s" => \$libname,
);

die "Usage: sam_len_cov_gc_pairing.pl -s SAMfile -f contigfastafile [-o outputprefix] [-i]\n" unless $contigfastafile and $samfile;

if (not $outfile) {
    $outfile = (-f $samfile) ? $samfile : "out";
};
$libname   = $outfile   unless $libname;
$fastafile = $outfile . ".lencovgc.fna";
$lencovgcfile = $outfile . ".lencovgc.txt";
@range_breaks = (0,100,200,300,400,500,1000,2000,5000,10000) unless @range_breaks;

# load contig file into memory:
print STDERR scalar localtime() . " - Loading contig fasta file $contigfastafile into memory ...\n";
my $fastahash = &fastafile2hash ($contigfastafile);
print STDERR scalar localtime() . " - Loading contig fasta file $contigfastafile into memory ... DONE\n";

#---------------------------------------------------------------
# process each line of samfile to get pairing info and cov info
open SAM, "$samfile" or die $!;

my (%read1, %read2, @F, $numpairs_unmapped, $numpairs_one_mapped, $numpairs_diff_mapped, $numpairs_bothmapped, %FR, %RF, $FF, $totalpairs, $readid, $unpaired_warning_printed);

my $reads_processed=0;
print STDERR scalar localtime() . " - Reading SAM input from $samfile ...\n";
while (<SAM>) {
	$reads_processed+=2;
	print STDERR "Processed $reads_processed reads by " . localtime() . "...\n" if $reads_processed % 1000000 == 0;
    $read1{contig} = $read2{contig} = "*";
    # ignore comment lines/ SAM headers
    next if /^\s*$/ or /^@/ or /^#/;
    chomp;
    @F=split/\t/;
    $readid = $1 if $F[0] =~ m|^(.*?)(/1)?$|;
    if (($F[1] & 4) != 4) {
        $read1{strand} = ($F[1] & 16) ? 0 : 1;
        $read1{contig} =  $F[2];
        $read1{tstart} =  $F[3];
        $read1{tend}   =  $F[3] - 1;
        while ($F[5] =~ /(\d+)/g) { $read1{tend} += $1 }
        die "$F[2] in samfile $samfile, but not in contigfile $contigfastafile\n" if not exists $$fastahash{$F[2]};
        $$fastahash{$F[2]}{cov} += ($read1{tend} - $read1{tstart} + 1)/$$fastahash{$F[2]}{len};
    }
    next if not $_ = <SAM>;
    chomp;
    @F=split/\t/;
    print STDERR "Warning: Output not paired\nInsert length estimates will probably be meaningless\n" and $unpaired_warning_printed = 1 if $F[0] !~ /^$readid/ and not $unpaired_warning_printed and $insertestimate;
    if (($F[1] & 4) != 4) {
        $read2{strand} = ($F[1] & 16) ? 0 : 1;
        $read2{contig} =  $F[2];
        $read2{tstart} =  $F[3];
        $read2{tend}   =  $F[3] - 1;
        while ($F[5] =~ /(\d+)/g) { $read2{tend} += $1 }
        die "$F[2] in samfile $samfile, but not in contigfile $contigfastafile\n" if not exists $$fastahash{$F[2]};
        $$fastahash{$F[2]}{cov} += ($read2{tend} - $read2{tstart} + 1)/$$fastahash{$F[2]}{len};
    }

    $totalpairs++;
    
    if (($read1{contig} eq $read2{contig}) and ($read1{contig} ne "*")) {
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
    elsif ($read1{contig} ne "*" and $read2{contig} ne "*" and $read1{contig} ne $read2{contig}) {
        $numpairs_diff_mapped++
    }
    elsif ($read1{contig} eq "*" and $read2{contig} eq "*") {  $numpairs_unmapped++}
    elsif ($read1{contig} eq "*" xor $read2{contig} eq "*") {$numpairs_one_mapped++};
}
print STDERR scalar localtime() . " - Reading SAM input from $samfile ... DONE\n";

#-------------------------------------------------------
# First, print all the len, cov, gc info and make plots

print STDERR scalar localtime() . " - Making len cov gc data file for R ...\n";
print STDERR scalar localtime() . " - Making len cov gc fasta file ...\n";

open  FASTA,    ">$fastafile" or die $!;
open  LENCOVGC, ">$lencovgcfile" or die $!;
#removed header from output, R script will have to provide headers itself
#print LENCOVGC  "read_set\tcontig_id\tlen\tcov\tgc\n";

for my $header (keys %{$fastahash}) { 
	my $length  = $$fastahash{$header}{len};
	my $gccount = $$fastahash{$header}{gc};
	my $desc    = $$fastahash{$header}{desc};
	my $seq     = $$fastahash{$header}{seq};
	my $nonatgc = $$fastahash{$header}{nonatgc};
	my $cov     = 0;
	$cov        = $$fastahash{$header}{cov} if exists $$fastahash{$header}{cov};
	$$fastahash{$header}{ALL_READS}{cov} += $cov;
	print LENCOVGC "$libname\t$header\t$length\t$cov\t".($gccount/($length-$nonatgc))."\n";
    print FASTA    ">".$header.$delimiter.$length.$delimiter.sprintf("%.4f",$cov).$delimiter.sprintf("%.4f",$gccount/($length-$nonatgc))."$desc\n$seq\n";
}

print STDERR scalar localtime() . " - Making len cov gc fasta file ... DONE\n";
print STDERR scalar localtime() . " - Making len cov gc data file for R ... DONE\n";

exit 0 unless $insertestimate;
#-----------------------------------------------
# If -i option used, then print all the pairing stats and distribution plots

print STDERR scalar localtime() . " - Printing insert frequencies and stats ...\n";

open HIST, ">$outfile.pairing.hist.txt" or die $!;
open STAT, ">$outfile.pairing.stat.txt" or die $!;

$numpairs_unmapped=0    if not defined $numpairs_unmapped;
$numpairs_one_mapped=0  if not defined $numpairs_one_mapped;
$numpairs_diff_mapped=0 if not defined $numpairs_diff_mapped;
$numpairs_bothmapped=0  if not defined $numpairs_bothmapped;
$FR{0}=0 if not defined $FR{0};
$RF{0}=0 if not defined $RF{0};
$FF=0 if not defined $FF;

foreach (sort { $a <=> $b } keys %FR) {
	print HIST "$libname\tFR\t$_\t$FR{$_}\n"
}
foreach (sort { $a <=> $b } keys %RF) {
	print HIST "$libname\tRF\t$_\t$RF{$_}\n"
}
print STAT "
Total pairs                                       : $totalpairs
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
print STAT "Num of FF/RR pairs: " . $FF . "\n-------------\n";

print STDERR scalar localtime() . " - Printing pairing frequencies and stats ... DONE\n";
#-----------------------------------------------


#############################################################################

sub fastafile2hash {
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $header;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $header = $1;
            $sequences{$header}{desc} = $2;
        }
        else {
            chomp $line;
            $sequences{$header}{seq}     .= $line;
            $sequences{$header}{len}     += length $line;
            $sequences{$header}{gc}      += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences{$header}{nonatgc} += ($line =~ tr/N/N/);
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
