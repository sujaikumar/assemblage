#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

my ($blast2goannotfile, $gff3file) =("","");

GetOptions (
  "blast2goannotfile:s" => \$blast2goannotfile,
  "gff3file:s"    => \$gff3file,
);

################

if (not $blast2goannotfile or not $gff3file) {
  print STDERR "\nUsage: merge_blast2goannot_gff.pl -b <blast2goannotfile> -g <gff file>
gff file should have 'transcript' or 'mRNA' features corresponding to blast2go protein names.
These features must have an ID= field in col 9 corresponding to blast2go protein names\n\n";
  exit 1;
}

my $blast2goannotfh = &read_fh ($blast2goannotfile);

my %note;
my %go;
my %ec;

while (<$blast2goannotfh>) {
  chomp;
  my @F = split /\t/;
  $note{$F[0]}.=$F[2]   if exists $F[2];
  $go  {$F[0]}{$F[1]}=1 if $F[1]=~/^GO:/;
  $ec  {$F[0]}{$F[1]}=1 if $F[1]=~/^EC:/;
}

################

my $gff3fh = &read_fh ($gff3file);

while (<$gff3fh>) {
  chomp;
  my @F = split /\t/;
  if (exists $F[2] and ($F[2] eq "mRNA" or $F[2] eq "transcript")) {
    my $id=$1 if $F[8]=~/ID=([^\;]+)/;
    if (exists $note{$id}) {
      print "$_;Note=$note{$id};Ontology_term=" . join(",",keys %{$go{$id}});
      print ";EC=" . join(",",keys %{$ec{$id}}) if exists $ec{$id};
      print "\n";
      next;
    }
  }
  print "$_\n";
}

################################

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
