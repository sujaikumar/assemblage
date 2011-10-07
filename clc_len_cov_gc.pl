#!/usr/bin/env perl

=head1 NAME

clc_len_cov_gc.pl

=head1 SYNOPSIS

clc_len_cov_gc.pl -c inputcasfile

=head1 DESCRIPTION

- Takes a clc .cas file as input (the reference contigs file referred to inside .cas has to be available on disk for the prog to work)
- creates gc vs len vs cov plots
- if -f fastafile specified, creates fastafile (puts len cov and gc info in fasta header, separated by delimiter, default "_")

Needs:

assembly_info (part of clc ngs cell suite - does not need license to run)
R (with library ggplot2)

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2011.05.22

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($casfile, $fastafile, $statsfile, $contigfile) = ("","","","");
my $delimiter = "_";
my @readfiles;

GetOptions (
    "casfile=s"  => \$casfile,
    "delimiter=s"=> \$delimiter,
    "fastafile:s"=> \$fastafile,
    "statsfile:s"=> \$statsfile,
);

$fastafile = $casfile . ".lencovgc.fna" unless $fastafile;
$statsfile = $casfile . ".lencovgc.txt" unless $statsfile;

open CAS, "assembly_info $casfile |" or die $!;
while (<CAS>) { /Contig files:/ and ($_ = <CAS>) and /(\S+)/ and $contigfile =$1 and last }

while (<CAS>) {
    next unless /Read files:/;
    while (<CAS>) {
        last if /^\s*$/;
        /^\s*(\S+)/ and push @readfiles, $1;
    }
    last if /^\s*$/;
}

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
    print STATS "ALL_READS\t$header\t$length\t$4\t".($gccount/($length-$nonatgc))."\n";
}
close CAS;
close FASTA;

for my $i (0..$#readfiles) {
    open CAS, "assembly_info -i ".($i+1)." $casfile |" or die $!; # assembly_info -i takes the 1st, 2nd... etc readfiles and gives length/gc according to that
    my $readfile = $readfiles[$i];
    while (<CAS>) {
        next unless /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*$/;
        my $header  = $$fastaidhash{$1}{header};
        my $length  = $$fastaidhash{$1}{len};
        my $gccount = $$fastaidhash{$1}{gc};
        my $nonatgc = $$fastaidhash{$1}{nonatgc};
        print STATS "$readfile\t$header\t$length\t$4\t".($gccount/($length-$nonatgc))."\n";
    }
    close CAS;
}
close STATS;

open  R, ">$statsfile.R" or die $!;
print R <<RCOMMAND;
library(ggplot2)
theme_set(theme_bw())
d <- read.table(\"$statsfile\",header=TRUE)
png(\"$statsfile.png\", 3600, 600 * length(levels(d\$read_set)))
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

system("R -f $statsfile.R");

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
