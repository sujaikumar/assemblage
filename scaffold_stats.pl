#!/usr/bin/env perl

=head1 NAME

scaffold_stats.pl - Takes one or more fasta files as input and calculates some assembly stats for it

=head1 SYNOPSIS

scaffold_stats.pl -fasta /path/to/file1.fasta -fasta /path/to/file2.fasta -threshold 100 -delim "," -h

=head1 DESCRIPTION

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2012.05.14

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use List::Util qw(sum);

my @fastafiles;
my @thresholds;
my $d                   = "\t";
my $output_prefix       = "scaffold_stats";
my $numNtosplitscaffold = 10;
my $minContigLength     = 100;
my $keepRfiles          = 0;

GetOptions (
    "fastafile=s{,}" => \@fastafiles,
    "threshold=i{,}" => \@thresholds,
    "delimiter=s"    => \$d,
    "output=s"       => \$output_prefix,
    "N=i"            => \$numNtosplitscaffold,
    "c=i"            => \$minContigLength,
    "k|keep"         => \$keepRfiles,
);
@thresholds = (200,1000) unless @thresholds;
@thresholds = sort   {$a<=>$b}  @thresholds;

#---------------------------------

die <<USAGE

Usage: scaffold_stats.pl -f contigs1.fa contigs2.fa -t 200 1000 -d ","
-t Threshold values are optional [-t 200 1000].
   Sequences shorter than the lowest threshold are not considered 
-N Minimum number of consecutive Ns for a scaffold to be split into contigs at that point [10] 
-o Output prefix. ["scaffold_stats"]
-d Delimiter for tabular output [TAB]
-c minimum contig length to report

USAGE
unless @fastafiles;

#### load fastafiles into hashes for scaffolds, contigs, and Ns

my %scaffold_hash;
my %contig_hash;
my %N_hash;

my %all_stats;

open RDATA, ">$output_prefix.Rdata" or die $!;

for my $fastafile (@fastafiles) {
    
    # Read in scaffold sequences from fasta files, only keep those above the lowest threshold length
    $scaffold_hash {$fastafile} = &fastafile2scaffoldhash  ($fastafile, $thresholds[0]);

    # Convert to contig hash and N hash
    $contig_hash   {$fastafile} = &scaffoldhash2contighash ($scaffold_hash{$fastafile});
    $N_hash        {$fastafile} = &scaffoldhash2Nhash      ($scaffold_hash{$fastafile});

    # Get stats for scaffold hash >lowest threshold and for contig hash and N hash
    # (threshold is set to 0 because contigs and runs of Ns have
    # no lower limit, length limits only apply to the original fasta/scaffold sequences)
    $all_stats     {$fastafile} {scaffold} { $thresholds[0] }  = &seqhash2stats ($scaffold_hash{$fastafile}, $thresholds[0]);
    $all_stats     {$fastafile} {contig}   { 0 }               = &seqhash2stats ($contig_hash  {$fastafile}, 0);
    $all_stats     {$fastafile} {N}        { 0 }               = &seqhash2stats ($N_hash       {$fastafile}, 0);

    # print lengths for scaffolds, contigs and Ns to RDATA file
    map { print RDATA "$fastafile\tScaffold\t$_\n" } @{$all_stats {$fastafile} {scaffold} {$thresholds[0]} {lengths}} ;
    map { print RDATA "$fastafile\tContig\t$_\n" }   @{$all_stats {$fastafile} {contig}   {0}              {lengths}} ;
    map { print RDATA "$fastafile\tN\t$_\n" }        @{$all_stats {$fastafile} {N}        {0}              {lengths}} ;

    # get scaffold stats at other thresholds (if any)
    for my $t (1..$#thresholds) {
        $all_stats {$fastafile} {scaffold} { $thresholds[$t] } = &seqhash2stats ($scaffold_hash{$fastafile}, $thresholds[$t]);
    }

}

# now print scaffold stats for multiple thresholds in human readable format
print "Filenames$d" . join($d,@fastafiles) . "\n";
print "LongestScaffold";
for my $fastafile (@fastafiles) {
    print $d . $all_stats {$fastafile} {scaffold} { $thresholds[0] } {max_length};
}
for my $threshold (@thresholds) {
    print "\nFor scaffolds longer than $threshold bp:\n";
    print "Num"        ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {num_seqs}   }; print "\n";
    print "Span"       ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {span}       }; print "\n";
    print "Min"        ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {min_length} }; print "\n";
    print "Mean"       ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {mean_length}}; print "\n";
    print "N50"        ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {n50}        }; print "\n";
    print "NumN50"     ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {num_n50}    }; print "\n";
    print "GC"         ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {scaffold} { $threshold } {gc}         }; print "\n";
}

print "\nFor contigs longer than $minContigLength bp (scaffolds split at >= $numNtosplitscaffold Ns):\n";
print "LongestContig"  ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {max_length} }; print "\n";
print "Num"            ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {num_seqs}   }; print "\n";
print "Span"           ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {span}       }; print "\n";
print "Min"            ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {min_length} }; print "\n";
print "Mean"           ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {mean_length}}; print "\n";
print "N50"            ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {n50}        }; print "\n";
print "NumN50"         ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {num_n50}    }; print "\n";
print "GC"             ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {contig} { 0 } {gc}         }; print "\n";

print "\nFor runs of Ns (>= $numNtosplitscaffold Ns):\n";
print "Num"        ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {N} { 0 } {num_seqs} }; print "\n";
print "Span"       ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {N} { 0 } {span}     }; print "\n";
print "N50"        ; for my $fastafile (@fastafiles) { print $d . $all_stats {$fastafile} {N} { 0 } {n50}      }; print "\n";

# now make pretty R plots using ggplot2

&make_cumulative_ggplot2("$output_prefix.Rdata"); 

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

#############################################################################

sub fastafile2scaffoldhash {
    my $fastafile     = shift @_;
    my $length_cutoff = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $seqid;
    while (<$fh>)
    {
        next if /^\s*$/ or /^#/;

        if (/^>(\S+)/) {
            $seqid = $1;
        }
        else {
            chomp($sequences{$seqid}{seq} .= $_  );
        }
    }
    foreach (keys %sequences) { delete $sequences{$_} if length($sequences{$_}{seq}) < $length_cutoff }
    return \%sequences;
}

#############################################################################

sub scaffoldhash2contighash {

    my $scaffoldhash = shift @_;
    my %contighash;
    my $contigid = 0;
    for my $seqid (keys %{$scaffoldhash}) {
        my $seq   = $$scaffoldhash{$seqid}{seq};
        foreach (split (/N{$numNtosplitscaffold,}/i,$seq)) {
            next if /^\s*$/; # skip blank entries which you get if seq begins with Ns
            next if length($_) <$minContigLength;
            $contighash{++$contigid}{seq} = $_;
        }
    }
    return \%contighash;
}

#############################################################################

sub scaffoldhash2Nhash {

    my $scaffoldhash = shift @_;
    my %Nhash;
    my $Nid = 0;
    for my $seqid (keys %{$scaffoldhash}) {
        my $seq   = $$scaffoldhash{$seqid}{seq};
        foreach (split (/[atgc]+/i,$seq)) {
            next if /^\s*$/; # skip blank entries which you get if seq begins with Ns
            next if length($_) <$numNtosplitscaffold;
            $Nhash{++$Nid}{seq} = $_;
        }
    }
    return \%Nhash;
}

#############################################################################

sub seqhash2stats {
    my $seqhash           = shift @_;
    my $length_threshold  = shift @_;
    my %stats;
    my $gc_sum            = 0;
    my $atgc_sum          = 0;
    my @sorted_lengths;
    foreach (keys %{$seqhash} ) {
        my $seq = $$seqhash{$_}{seq};
        next if length ($seq) < $length_threshold;
        push @sorted_lengths, length ($seq);
        $gc_sum               += ($seq =~ tr/gcGC/gcGC/);
        $atgc_sum             += ($seq =~ tr/atgcATGC/atgcATGC/);
    }
    $stats{gc}            = $atgc_sum        ? sprintf("%.3f",$gc_sum/$atgc_sum)   : 0;
    @sorted_lengths       = sort {$b <=> $a} @sorted_lengths;
    $stats{lengths}       = \@sorted_lengths;
    $stats{min_length}    =  @sorted_lengths ? $sorted_lengths[$#sorted_lengths]   : 0;
    $stats{max_length}    =  @sorted_lengths ? $sorted_lengths[0]                  : 0;
    $stats{num_seqs}      =  @sorted_lengths ? @sorted_lengths                     : 0;
    $stats{span}          =  @sorted_lengths ? sum(@sorted_lengths)                : 0;
    $stats{mean_length}   =  @sorted_lengths ? int($stats{span}/$stats{num_seqs})  : 0;
    $stats{median_length} =  @sorted_lengths ? $sorted_lengths[$stats{num_seqs}/2] : 0;
    my ($csum, $nlen, $n50, $num_n50) = (0,0,0,0);
    for $nlen (@sorted_lengths) { $csum += $nlen; $num_n50++; $n50 = $nlen; last if $csum >= ($stats{span}/2) }
    $stats{n50}           = $n50;
    $stats{num_n50}       = $num_n50;

    return \%stats;
}

#############################################################################

sub make_cumulative_ggplot2 {

my $rdatafile = shift @_;
open RSCRIPT, ">$rdatafile.R" or die $!;
    
print RSCRIPT '
library(ggplot2)
library(grid)
args <- commandArgs(trailingOnly = TRUE)

data=read.delim(args[1],header=F,col.names=c("file","type","len"))
for (f in levels(data$file)) {
    for (t in levels(data$type)) {
        data[data$file==f & data$type==t,"cumsum"] =cumsum   (data[data$file==f & data$type==t,"len"])
        data[data$file==f & data$type==t,"Rank"]   =seq_along(data[data$file==f & data$type==t,"len"])
    }
}

theme_set(theme_bw())
Scaffold_curve<-qplot(Rank,cumsum,data=data[data$type=="Scaffold",],geom="line",color=file) +  ylab ("Cumulative Length") + 
  opts(legend.justification=c(1,0), legend.position=c(1,0), legend.title=theme_blank(), title="Scaffolds") +
  xlab ("Scaffolds ranked by size")
Contig_curve  <-qplot(Rank,cumsum,data=data[data$type=="Contig",],  geom="line",color=file) +  ylab ("Cumulative Length") + 
  opts(legend.justification=c(1,0), legend.position=c(1,0), legend.title=theme_blank(), title="Contigs") +
  xlab ("Contigs ranked by size")
N_curve       <-qplot(Rank,cumsum,data=data[data$type=="N",],       geom="line",color=file) +  ylab ("Cumulative Length") + 
  opts(legend.justification=c(1,0), legend.position=c(1,0), legend.title=theme_blank(), title="Ns") +
  xlab ("Blocks of Ns ranked by size")

png(paste(args[1],".png",sep=""),2700,900,res=200)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)

print(Scaffold_curve, vp = vplayout(1, 1))
print(Contig_curve,   vp = vplayout(1, 2))
print(N_curve,        vp = vplayout(1, 3))
dev.off()';
close RSCRIPT;
system "Rscript $rdatafile.R $rdatafile &>/dev/null";
unlink  "$rdatafile.R", $rdatafile unless $keepRfiles;
}
