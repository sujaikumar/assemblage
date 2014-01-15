#!/usr/bin/env perl

# takes a pairwise or multi maf file and 
# outputs a GFF for each species (assuming species is listed in chrontig name: eg:
# ce.CHROMOSOME_I... )
# -b outputs bed
# -g outputs gff (default)

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($maf_file, $bed, $gff, $prefix) = ("-",0,1,"");

GetOptions (
  "m|maf:s"    => \$maf_file,
  "b|bed"      => \$bed,
  "g|gff"      => \$gff,
  "p|prefix:s" => \$prefix,
);
$gff = 0 if $bed;

################################################################################################
# open maf file
# read lines beginning with "s " as one block.

my $maf_file_fh = &read_fh ($maf_file);

my ($block_count, $block_id, $species, $chrontig, $strand, $start, $end);

my %overall_gff; #hash to store GFFs. key = species;

while (<$maf_file_fh>) {
    next unless /^a score/;
	my %current_gff; #hash to store GFFs. key = species;
	$block_id++;
    $block_count = 0;
    while ($_ = <$maf_file_fh>) {
        last if /^\s*$/;
        if (/^s\s+(([^.]+)\.\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\d+)/) {
#                 1=chrontig,2=sp  3=st   4=len   5=strand 6=total
            $block_count++;
            $chrontig  = $1;
            $species   = $2;
            $strand    = $5;
            if ($strand eq "+") {
                $start = $3;
                $end   = $3 + $4;
            }
            else {
                $start = $6 - $3 - $4;
                $end   = $6 - $3;
            }
            $current_gff {$species} = join("\t",$chrontig,"tba.maf","cross_genome_match",$start + 1,$end,".",$strand,".","block$block_id") . "\n" if $gff;
            $current_gff {$species} = join("\t",$chrontig,$start,$end,"block=$block_id",".",$strand) . "\n" if $bed;
        }
    }
    if ($block_count > 1) {
        foreach (keys %current_gff) {
            $overall_gff{$_} .= $current_gff{$_}
        }
    }
}

foreach (keys %overall_gff) {
    my $filename = "$prefix$_.maf.to." . ($gff ? "gff" : "bed");
    open  OUT, ">$filename" or die $_;
    print OUT  $overall_gff{$_};
}

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
