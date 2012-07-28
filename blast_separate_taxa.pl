#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($b1file, $b2file, $scorediff, $bothfile) = ("","",50, "hitting.both");
my ($evalue, $bitscore, $length, $perc_id);
GetOptions (
  "b1:s" => \$b1file,
  "b2:s" => \$b2file,
  "together:s" => \$bothfile,
  "diff:i" => \$scorediff,
  "evalue:f" => \$evalue,
  "bitscore:f" => \$bitscore,
  "length:i" => \$length,
  "perc_id:f" => \$perc_id,
);

die "Usage: blast_separate_taxa.pl -b1 <blast1file> -b2 <blast2file> [-d <min bitscore difference, default 50>]" unless $b1file and $b2file;

my %b1hash;
my $b1fh = &read_fh($b1file);
while (<$b1fh>)
{
	chomp;
	my @hsp = split /\t/;
	next unless scalar @hsp >= 12; #at least 12 fields in a blast tabular hit
	next if defined $evalue and $hsp[10] > $evalue;
	next if defined $bitscore and $hsp[11] < $bitscore;
	next if defined $length and $hsp[3] < $length;
	next if defined $perc_id and $hsp[2] < $perc_id;
	$b1hash{$hsp[0]} = $hsp[11] unless exists $b1hash{$hsp[0]} and $b1hash{$hsp[0]} > $hsp[11]; #store best bit score
}
close $b1fh;
my %b2hash;
my $b2fh = &read_fh($b2file);
while (<$b2fh>)
{
	chomp;
	my @hsp = split /\t/;
	next unless scalar @hsp >= 12; #at least 12 fields in a blast tabular hit
	next if defined $evalue and $hsp[10] > $evalue;
	next if defined $bitscore and $hsp[11] < $bitscore;
	next if defined $length and $hsp[3] < $length;
	next if defined $perc_id and $hsp[2] < $perc_id;
	$b2hash{$hsp[0]} = $hsp[11] unless exists $b2hash{$hsp[0]} and $b2hash{$hsp[0]} > $hsp[11]; #store best bit score
}
close $b2fh;

open B1ONLY, ">$b1file.only" or die $!;
open B2ONLY, ">$b2file.only" or die $!;
open BOTH, ">$bothfile" or die $!;
for my $b1key (keys %b1hash)
{
	if (exists $b2hash{$b1key})
	{
		if ($b1hash{$b1key} >= $b2hash{$b1key} + $scorediff)
		{
			print B1ONLY "$b1key\t$b1hash{$b1key}\n";
		}
		elsif ($b2hash{$b1key} >= $b1hash{$b1key} + $scorediff)
		{
			print B2ONLY "$b1key\t$b2hash{$b1key}\n";
		}
		else
		{
			print BOTH "$b1key\t$b1hash{$b1key}\t$b2hash{$b1key}\n";
		}
		delete $b2hash{$b1key}; # b2hash{b1key} removed so that it isnt printed when printing b2 only
	}
	else
	{
		print B1ONLY "$b1key\t$b1hash{$b1key}\n";
	}
}
close B1ONLY;
close BOTH;

for my $b2key (keys %b2hash)
{
	print B2ONLY "$b2key\t$b2hash{$b2key}\n";
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
