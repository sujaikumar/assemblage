#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

my ($blastfile,$dbfile, $queryfile, $excludefile) = ("","","","");
my ($evalue, $perc_identity, $length, $bitscore, $max_target_seqs, $hsp_per_target, $qregex, $dregex) = ("","","","","","","","");
my ($qcov, $dcov, $combinehsp) = ("","","");

GetOptions (
  "db:s" => \$dbfile,
  "query:s" => \$queryfile,
  "blastfile:s" => \$blastfile,
  "excludefile:s" => \$excludefile,
  "perc_identity:f" => \$perc_identity,
  "evalue:f" => \$evalue,
  "length:i" => \$length,
  "score:f" => \$bitscore,
  "max_target_seqs:i" => \$max_target_seqs,
  "hsp_per_target:i" => \$hsp_per_target,
  "qcov:s" => \$qcov,
  "dcov:s" => \$dcov,
  "qregex:s" => \$qregex,
  "dregex:s" => \$dregex,
  "combinehsp" => \$combinehsp,
);

my ($current_target_seqs, $current_hsps, $qid, $did) = (0,0,"","");

if (!$blastfile or ($qcov and not $queryfile) or ($dcov and not $dbfile) or ($combinehsp and not $dcov and not $qcov)
	or ($combinehsp and $hsp_per_target))
{
	die <<USAGE;
Usage:
blastm8_filter.pl
	-b, --blastfile : <tabular/m8/outfmt6 blast output> MANDATORY, all others optional
	-ex,--excludefile:<tabular file with col1=contig, col2=st, col3=en. If any blast hit
	                  overlaps these intervals, it will be excluded
	-e, --evalue    : evalue cutoff eg. -e 1e-5 filters hits <= 1e-5
	-l, --length    : length cutoff eg. -l 100 filters hits with alignment length >= 100
	-s, --bitscore  : bitscore cutoff (-s 200 filters hits with bitscore >= 100
	-q, --qregex    : regexp to filter only some query sequences,
	                  eg -q "^clc" only picks queries beginning with clc
	-d, --dregex    : regexp to filter only some database/target/ref sequences,
	                  eg -d "_c|_s" only picks those database sequences that match _c or _s
	-m, --max_target_seqs  : assuming blast file is sorted as output by blast/blast+,
	                  -m 1 will pick only 1 target seq per query
	-h, --hsp_per_target   : assuming blast file is sorted as output by blast/blast+,
	                  -h 1 will pick only 1 hsp per query_target pair.
	                  Using with -m 1 guarantees only one hsp per query
	    --qcov      : proportion of query covered by hit (0.0-1.0). Must specify query fasta file
	    --query     : query fasta file
	    --dcov      : proportion of database covered by hit (0.0-1.0). Must specify db fasta file
	    --db        : db fasta file
	-c, --combinehsp: Use this ONLY if you want to combine all hsps for a hit when calculating qcov/dcov. Doesn't make sense
	                  to use this with -h 1
USAGE
}

my %exclude_hash;
if (-r $excludefile)
{
    my $exclude_fh = &read_fh($excludefile);
    while (<$exclude_fh>)
    {
        next unless /^>?(\S+?)[\t_ ](\d+)[\t_ ](\d+)$/;
        push @{$exclude_hash{$1}}, [$2,$3];
    }
    for my $chrontig (keys %exclude_hash)
    {
        @{$exclude_hash{$chrontig}} = sort { $a->[0] <=> $b->[0] } @{$exclude_hash{$chrontig}};
    }
}

my $qsequences = &fastafile2hash($queryfile) if $queryfile;
my $dsequences = &fastafile2hash($dbfile) if $dbfile;

my ($current_hit_q_string, $current_hit_d_string, $current_hit_toprint) = ("","","");

my $blastfh = &read_fh($blastfile);
while (my $line = <$blastfh>)
{
	next if $line =~ /^#/ or $line =~ /^\s*$/;
	my @blastm8 = split(/\t/, $line);
	die "Input not in blast -m 8 or blast+ -outfmt 6 tabular format\n" if scalar @blastm8 < 12;
	next if $evalue and eval($blastm8[10]) > $evalue;
	next if $length and $blastm8[3] < $length;
	next if $perc_identity and $blastm8[2] < $perc_identity;
	next if $bitscore and $blastm8[11] < $bitscore;
	next if $qregex and $blastm8[0] !~ /$qregex/;
	next if $dregex and $blastm8[1] !~ /$dregex/;
	my ($qst, $qen, $dst, $den) = ($blastm8[6],$blastm8[7],$blastm8[8],$blastm8[9]);
	($qst,$qen) = ($qen,$qst) if $qst > $qen;
	($dst,$den) = ($den,$dst) if $dst > $den;
	my $qlen = $qen - $qst + 1;
	my $dlen = $den - $dst + 1;
	if ($blastm8[0] =~ /^(.+?)_(\d+)_(\d+)$/) { $qst = $qst + $2 - 1; $qen = $qen + $2 - 1; $blastm8[0] = $1 }
	if ($blastm8[1] =~ /^(.+?)_(\d+)_(\d+)$/) { $dst = $dst + $2 - 1; $den = $den + $2 - 1; $blastm8[1] = $1 }

    next if exists $exclude_hash{$blastm8[0]} and &binary_search_intervals ( $exclude_hash{$blastm8[0]}, $qst, $qen );
    next if exists $exclude_hash{$blastm8[1]} and &binary_search_intervals ( $exclude_hash{$blastm8[1]}, $dst, $den );

	next if $qcov and not defined $$qsequences{$blastm8[0]};
	next if $dcov and not defined $$dsequences{$blastm8[1]};
    
	if (not $combinehsp)
	{
		next if $qcov and $qlen/length($$qsequences{$blastm8[0]}{seq}) < $qcov;
		next if $dcov and $dlen/length($$dsequences{$blastm8[1]}{seq}) < $dcov;
	}
	if ($qid ne $blastm8[0])
	{
		$current_target_seqs = 1;
	}
	if (($qid eq $blastm8[0] or $qid eq "") and $did ne $blastm8[1])
	{
		$current_target_seqs++;
	}
	if("$qid$did" ne "$blastm8[0]$blastm8[1]")
	{
		# new hit
		$current_hsps = 1;
		# check if old hsp needs to be printed if $combinehsp is on
		if ($combinehsp and "$qid$did")
		{
			print $current_hit_toprint unless
				$qcov and ($current_hit_q_string=~tr/1/1/)/length($$qsequences{$qid}{seq}) <$qcov or
				$dcov and ($current_hit_d_string=~tr/1/1/)/length($$dsequences{$did}{seq}) <$dcov;
		}
		$current_hit_toprint = $line;
		$current_hit_q_string = "0" x length($$qsequences{$blastm8[0]}{seq}) if $queryfile;
		$current_hit_d_string = "0" x length($$dsequences{$blastm8[1]}{seq}) if $dbfile;
	} else
	{
		$current_hsps++
	}
	$qid = $blastm8[0];
	$did = $blastm8[1];
	next if $max_target_seqs and $current_target_seqs > $max_target_seqs;
	next if $hsp_per_target  and $current_hsps > $hsp_per_target;
	if ($combinehsp)
	{
		substr($current_hit_q_string, $qst -1, $qlen, "1" x $qlen) if $queryfile;
		substr($current_hit_d_string, $dst -1, $dlen, "1" x $dlen) if $dbfile;
		$current_hit_toprint .= $line if $current_hsps >1;
	} else
	{
		print $line;
	}
}

###################################################################################################
sub fastafile2hash
{
	my $fastafile = shift @_;
	my $changecase = "N"; $changecase = shift @_ if @_;
	my %sequences;
	my $seqid = "";
	my $fh = &read_fh($fastafile);
	while (<$fh>)
	{
		if (/^>(\S+)(.*)/) {
			$seqid = $1;
			$sequences{$seqid}{desc} = $2;
		}
		elsif ($changecase eq "L") { chomp($sequences{$seqid}{seq} .= lc($_)) }
		elsif ($changecase eq "U") { chomp($sequences{$seqid}{seq} .= uc($_)) }
		elsif ($changecase eq "N") { chomp($sequences{$seqid}{seq} .= $_) }
	}
	close $fh;
	return \%sequences;
}
###################################################################################################

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
###################################################################################################

sub binary_search_intervals
{
    my ($array_ref, $range_st, $range_en) = @_;
    my @sorted = @$array_ref;
    
    my $low_index = 0; my $high_index = $#sorted;
    my $found = 0;
    while ( $low_index <= $high_index and not $found )
    {
        my $mid_index = int ( ( $low_index + $high_index ) / 2 );
        if ( $sorted[$mid_index]->[0] <= $range_en and $sorted[$mid_index]->[1] >= $range_st )
        {
            $found = 1;
            return $mid_index + 1;
        }
        if ( $sorted[$mid_index]->[0] < $range_en )
        {
            $low_index = $mid_index + 1;
        } else
        {
            $high_index = $mid_index - 1;
        }
    }
    return $found;
}
