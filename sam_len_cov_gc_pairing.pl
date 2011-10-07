#!/usr/bin/env perl

=head1 NAME

sam_len_cov_gc_pairing.pl

=head1 SYNOPSIS

sam_len_cov_gc_pairing.pl -s samfile -c contigfile [-o outputprefix]

=head1 DESCRIPTION

creates
  outputprefix.pairing.hist.txt
  outputprefix.pairing.stat.txt
  outputprefix.pairing.R
  outputprefix.pairing.png

samfile must be created using bwa in single end mode with interleaved reads
outputprefix is optional, by default all output files will be in the same directory as samfile
  with samfile.pairing as prefix

assumes all reads are from same library and creates just one plot.

script creates assembly_table -n output and outputs R readable data file and then runs R script that uses ggplot2:
TYPE    INSSIZE FREQ
--------------------
FF      20      101
FF      21      103
FF      22      110
etc
Type can be:
FF (both in same dir, same as RR)
FR (reads pointing to each other)
RF (reads pointing away from each other)
All distances measured from start of one read to start of other read
STDERR gets:
Count of read pairs not mapped at all (both reads in pair not mapped)
Count of read pairs where one mapped but not other
Count of read pairs mapped to different contigs
Count of read pairs mapped to same contig

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2011.09.14

=cut

use strict;
use warnings;
use List::Util qw(min max sum);
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($samfile, $outfile, $contigfile, $fastafile, $lencovgcfile) = ("","","");
my $delimiter = "_";
my @range_breaks = (0,100,200,300,400,500,1000,2000,5000,10000);
GetOptions (
    "samfile=s"     => \$samfile,
    "outfile=s"     => \$outfile,
    "contigfile=s"  => \$contigfile,
    "range=i{,}"    => \@range_breaks,
    "fastafile:s"   => \$fastafile,
    "lencovgcfile:s"=> \$lencovgcfile,
    "delimiter=s"=> \$delimiter,
);

die "Usage: sam_len_cov_gc_pairing.pl -s samefile -c contigfile [-o outputprefix]\n" unless $contigfile and $samfile;

$outfile = "$samfile" unless $outfile;
$outfile = "samfile"  unless -f $samfile;

$fastafile    = $outfile . ".lencovgc.fna" unless $fastafile;
$lencovgcfile = $outfile . ".lencovgc.txt" unless $lencovgcfile;

# load contig file into memory:
print STDERR scalar localtime . " - Loading $contigfile into memory ...\n";
my $fastahash = &fastafile2hash ($contigfile);
print STDERR scalar localtime . " - Loading $contigfile into memory ... DONE\n";

#---------------------------------------------------------------
# process each line of samfile to get pairing info and cov info
open SAM, "$samfile" or die $!;

my (%read1, %read2, @F, %numpairs_unmapped, %numpairs_one_mapped, %numpairs_diff_mapped, %numpairs_bothmapped, %FR, %RF, $lib, %FF, %totalpairs, $readid);

$lib = "library";

print STDERR scalar localtime . " - Reading SAM input from $samfile ...\n";
while (<SAM>) {
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
        die "$F[2] in samfile $samfile, but not in contigfile $contigfile\n" if not exists $$fastahash{$F[2]};
        $$fastahash{$F[2]}{$lib}{cov} += ($read1{tend} - $read1{tstart} + 1)/$$fastahash{$F[2]}{len};
    }
    $_ = <SAM>;
    chomp;
    @F=split/\t/;
    die "Output not paired\n" if $F[0] !~ m|^$readid(/2)?$|;
    if (($F[1] & 4) != 4) {
        $read2{strand} = ($F[1] & 16) ? 0 : 1;
        $read2{contig} =  $F[2];
        $read2{tstart} =  $F[3];
        $read2{tend}   =  $F[3] - 1;
        while ($F[5] =~ /(\d+)/g) { $read2{tend} += $1 }
        die "$F[2] in samfile $samfile, but not in contigfile $contigfile\n" if not exists $$fastahash{$F[2]};
        $$fastahash{$F[2]}{$lib}{cov} += ($read2{tend} - $read2{tstart} + 1)/$$fastahash{$F[2]}{len};
    }
    
    $totalpairs{$lib}++;
    die "Reads not paired correctly\n" if scalar keys %totalpairs >100;
    
    if (($read1{contig} eq $read2{contig}) and ($read1{contig} ne "*")) {
        $numpairs_bothmapped{$lib}++;
        if ($read1{strand} == $read2{strand}) {
            $FF{$lib}++
        }
        elsif     ( $read1{strand} ) {
            if    ( $read2{tend}   >= $read1{tstart} ) { $FR{$lib}{ $read2{tend} - $read1{tstart} }++ }
            elsif ( $read2{tstart} <= $read1{tend} )   { $RF{$lib}{ $read1{tend} - $read2{tstart} }++ }
        }
        elsif     ( $read2{strand} ) {
            if    ( $read1{tend}   >= $read2{tstart} ) { $FR{$lib}{ $read1{tend} - $read2{tstart} }++ }
            elsif ( $read1{tstart} <= $read2{tend} )   { $RF{$lib}{ $read2{tend} - $read1{tstart} }++ }
        }
        next;
    }
    elsif ($read1{contig} ne "*" and $read2{contig} ne "*" and $read1{contig} ne $read2{contig}) {
        $numpairs_diff_mapped{$lib}++
    }
    elsif ($read1{contig} eq "*" and $read2{contig} eq "*") {  $numpairs_unmapped{$lib}++}
    elsif ($read1{contig} eq "*" xor $read2{contig} eq "*") {$numpairs_one_mapped{$lib}++};
}
print STDERR scalar localtime . " - Reading SAM input from $samfile ... DONE\n";

#-----------------------------------------------
# First, print all the pairing stats and distribution plots

print STDERR scalar localtime . " - Printing pairing frequencies and stats ...\n";

open HIST, ">$outfile.pairing.hist.txt" or die $!;
open STAT, ">$outfile.pairing.stat.txt" or die $!;

for my $lib (sort keys %totalpairs) {

    $numpairs_unmapped{$lib}=0    if not defined $numpairs_unmapped{$lib};
    $numpairs_one_mapped{$lib}=0  if not defined $numpairs_one_mapped{$lib};
    $numpairs_diff_mapped{$lib}=0 if not defined $numpairs_diff_mapped{$lib};
    $numpairs_bothmapped{$lib}=0  if not defined $numpairs_bothmapped{$lib};
    $FR{$lib}{0}=0 if not defined $FR{$lib};
    $RF{$lib}{0}=0 if not defined $RF{$lib};
    $FF{$lib}=0    if not defined $FF{$lib};

    foreach (sort { $a <=> $b } keys %{$FR{$lib}}) {
        print HIST "$lib\tFR\t$_\t$FR{$lib}{$_}\n"
    }
    foreach (sort { $a <=> $b } keys %{$RF{$lib}}) {
        print HIST "$lib\tRF\t$_\t$RF{$lib}{$_}\n"
    }
    print STAT "$lib
    Total pairs                                       : $totalpairs{$lib}
    Pairs with neither read mapped                    : $numpairs_unmapped{$lib} (". sprintf("%.1f",100*$numpairs_unmapped{$lib}/$totalpairs{$lib}). " %)
    Pairs with one read mapped                        : $numpairs_one_mapped{$lib} (". sprintf("%.1f",100*$numpairs_one_mapped{$lib}/$totalpairs{$lib}). " %)
    Pairs with both reads mapping different contigs   : $numpairs_diff_mapped{$lib} (". sprintf("%.1f",100*$numpairs_diff_mapped{$lib}/$totalpairs{$lib}). " %)
    Pairs with both reads best mapping to same contig : $numpairs_bothmapped{$lib} (". sprintf("%.1f",100*$numpairs_bothmapped{$lib}/$totalpairs{$lib}). " %)\n";
    print STAT "Num of FR    pairs: " . sum(values %{$FR{$lib}}) . " (".sprintf("%.1f",100 * sum(values %{$FR{$lib}})/$numpairs_bothmapped{$lib})." % of $numpairs_bothmapped{$lib} pairs mapping to same contig)\n";
    my %range_counts;
    for my $ins (keys %{$FR{$lib}}) {
        for my $range_index (1..$#range_breaks) {
            if ($ins >= $range_breaks[$range_index -1] and $ins < $range_breaks[$range_index]) {
                $range_counts{FR}[$range_index] += $FR{$lib}{$ins};
                last;
            }
        }
    }
    for my $range_index (1..$#range_breaks) {
        print STAT "    FR Pairs with insert size range ".$range_breaks[$range_index -1]."-".$range_breaks[$range_index] ." : ". (defined $range_counts{FR}[$range_index] ? $range_counts{FR}[$range_index] : 0) . "\n";
    }
    print STAT "Num of RF    pairs: " . sum(values %{$RF{$lib}}) . " (".sprintf("%.1f",100 * sum(values %{$RF{$lib}})/$numpairs_bothmapped{$lib})." % of $numpairs_bothmapped{$lib} pairs mapping to same contig)\n";
    for my $ins (keys %{$RF{$lib}}) {
        for my $range_index (1..$#range_breaks) {
            if ($ins >= $range_breaks[$range_index -1] and $ins < $range_breaks[$range_index]) {
                $range_counts{RF}[$range_index] += $RF{$lib}{$ins};
                last;
            }
        }
    }
    for my $range_index (1..$#range_breaks) {
        print STAT "    RF Pairs with insert size range ".$range_breaks[$range_index -1]."-".$range_breaks[$range_index] ." : ". (defined $range_counts{RF}[$range_index] ? $range_counts{RF}[$range_index] : 0) . "\n";
    }
    print STAT "Num of FF/RR pairs: " . $FF{$lib} . "\n-------------\n";
}

print STDERR scalar localtime . " - Printing pairing frequencies and stats ... DONE\n";

print STDERR scalar localtime . " - Making pairing plots in R ...\n";
open  R, ">$outfile.R" or die $!;
print R <<RCOMMAND;
d=read.table("$outfile.pairing.hist.txt",header=FALSE)
library(ggplot2)
theme_set(theme_bw())
png("$outfile.hist.png",3000,500 * length(levels(d\$V1)))
pairlin <- qplot(V3,V4,data=d,facets = V1 ~ V2, xlab='insert size',ylab='freq') + xlim(0,1000)
pairlog <- qplot(V3,V4,data=d,facets = V1 ~ V2, xlab='insert size',ylab='freq') + scale_x_log10() + scale_y_log10()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)
print(pairlin, vp = vplayout(1, 1))
print(pairlog, vp = vplayout(1, 2))
dev.off()
RCOMMAND

system "R -f $outfile.R";

print STDERR scalar localtime . " - Making pairing plots in R ... DONE\n";

#-----------------------------------------------
# Second, print all the len, cov, gc info and make plots

print STDERR scalar localtime . " - Making len cov gc data file for R ...\n";
print STDERR scalar localtime . " - Making len cov gc fasta file ...\n";

open  FASTA,    ">$fastafile" or die $!;
open  LENCOVGC, ">$lencovgcfile" or die $!;
print LENCOVGC  "read_set\tcontig_id\tlen\tcov\tgc\n";

# for each readset:
for my $lib (sort keys %totalpairs) {
    for my $header (keys %{$fastahash}) { 
        my $length  = $$fastahash{$header}{len};
        my $gccount = $$fastahash{$header}{gc};
        my $nonatgc = $$fastahash{$header}{nonatgc};
        my $cov     = 0;
        $cov        = $$fastahash{$header}{$lib}{cov} if exists $$fastahash{$header}{$lib}{cov};
        $$fastahash{$header}{ALL_READS}{cov} += $cov;
        print LENCOVGC "$lib\t$header\t$length\t$cov\t".($gccount/($length-$nonatgc))."\n";
    }
}
for my $header (keys %{$fastahash}) { 
    my $desc    = $$fastahash{$header}{desc};
    my $length  = $$fastahash{$header}{len};
    my $gccount = $$fastahash{$header}{gc};
    my $nonatgc = $$fastahash{$header}{nonatgc};
    my $seq     = $$fastahash{$header}{seq};
    my $cov     = $$fastahash{$header}{ALL_READS}{cov};
    print FASTA    ">".$header.$delimiter.$length.$delimiter.$cov.$delimiter.sprintf("%.2f",$gccount/($length-$nonatgc))."$desc\n$seq\n";
    print LENCOVGC "ALL_READS\t$header\t$length\t$cov\t".($gccount/($length-$nonatgc))."\n";
}

print STDERR scalar localtime . " - Making len cov gc fasta file ... DONE\n";
print STDERR scalar localtime . " - Making len cov gc data file for R ... DONE\n";

print STDERR scalar localtime . " - Plotting len cov gc data using ggplot2 in R ...\n";

open  R, ">$lencovgcfile.R" or die $!;
print R <<RCOMMAND;
library(ggplot2)
theme_set(theme_bw())
d <- read.table(\"$lencovgcfile\",header=TRUE)
png(\"$lencovgcfile.png\", 3600, 600 * length(levels(d\$read_set)))
lencovdot     <- qplot(len,cov,data=d,log="xy",facets= read_set ~ .,alpha=I(1/20))
gccovdot      <- qplot(gc ,cov,data=d,log="y", facets= read_set ~ .,alpha=I(1/20)) + xlim(0,1)
gclendot      <- qplot(gc ,len,data=d,log="y", facets= read_set ~ .,alpha=I(1/20)) + xlim(0,1)
gccovdotlin500<- qplot(gc ,cov,data=d,         facets= read_set ~ .,alpha=I(1/20)) + xlim(0,1) + ylim(0,500)
gccovdotlin100<- qplot(gc ,cov,data=d,         facets= read_set ~ .,alpha=I(1/20)) + xlim(0,1) + ylim(0,100)
d <- d[d\$cov>0,]
lencovden     <- qplot(len,cov,data=d,log="xy",     facets= read_set ~ .) +             stat_density2d(geom = "tile",  aes(fill = ..density..),contour = F) + scale_fill_gradientn(colours=c("white","black","red"))
gccovden      <- qplot(gc ,cov,data=d,log="y",      facets= read_set ~ .) + xlim(0,1) + stat_density2d(geom = "tile",  aes(fill = ..density..),contour = F) + scale_fill_gradientn(colours=c("white","black","red"))
gclenden      <- qplot(gc ,len,data=d,log="y",      facets= read_set ~ .) + xlim(0,1) + stat_density2d(geom = "tile",  aes(fill = ..density..),contour = F) + scale_fill_gradientn(colours=c("white","black","red"))
gccovdenlin   <- qplot(gc ,cov,data=d,ylim=c(0,300),facets= read_set ~ .) + xlim(0,1) + stat_density2d(geom = "tile",  aes(fill = ..density..),contour = F) + scale_fill_gradientn(colours=c("white","black","red"))
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 6)))
vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)
print(lencovden,      vp = vplayout(1, 1))
print(gccovden,       vp = vplayout(1, 2))
print(gccovdot,       vp = vplayout(1, 3))
print(gccovdotlin500, vp = vplayout(1, 4))
print(gccovdotlin100, vp = vplayout(1, 5))
print(gclenden,       vp = vplayout(1, 6))
dev.off()
RCOMMAND

system("R -f $lencovgcfile.R");

print STDERR scalar localtime . " - Plotting len cov gc data using ggplot2 in R ... DONE\n";

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
