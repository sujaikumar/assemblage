#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($opt_evalue, $opt_length, $opt_perc_identity, $opt_score);
my ($Length_proportion,$both, $queryfile, $dbfile, $twoway) = (0.9,"F","","","");
my ($merged_overlaps_file) = ("");
GetOptions (
  "evalue:f" => \$opt_evalue,
  "length:i" => \$opt_length,
  "score:f" => \$opt_score,
  "S|Similarity:f" => \$opt_perc_identity,
  "L|proportion:f" => \$Length_proportion,
  "both:s" => \$both,
  "i|query:s" => \$queryfile,
  "db:s" => \$dbfile,
  "merged:s" => \$merged_overlaps_file,
  "twoway" => \$twoway,
);

my %maskstrings;
my $blastinfile = "";
my $doblastclust = 0;
my ($queryseqlengths, $dbseqlengths) = ("","");

if (-r $queryfile and -r $dbfile)
{
	$doblastclust = 1;
	$opt_perc_identity = 90  unless $opt_perc_identity;
}

################################################################
# read in querylengths and dblengths
if ($doblastclust)
{
	&multilinefastafile2singleline($queryfile);
	$queryseqlengths = &fastafile2sequencelengths($queryfile);
	&multilinefastafile2singleline($dbfile);
	$dbseqlengths = &fastafile2sequencelengths($dbfile);
}

################################################################
# read in blastfile ensuring thresholds are met, create mask strings

my $tmp_blastfile = "link_blast.pl.tmp." . int(1000 * rand);
open TMP,">$tmp_blastfile" or die $!;
while (<>)
{
	next if /^#/;
	my @fields = split /\t/;
	next unless scalar @fields >= 12;
	my ($qid, $tid, $perc_identity, $length, $gaps, $mismatches, $qst, $qen, $tst, $ten, $evalue, $score) = @fields;
	next if $qid eq $tid; # don't compare a sequence to itself
	next if $opt_evalue and $evalue > $opt_evalue;
	next if $opt_length and $length < $opt_length;
	next if $opt_perc_identity and $perc_identity < $opt_perc_identity;
	if ($doblastclust)
	{
		next if ($both eq "T") and ($length/$$queryseqlengths{$qid} < $Length_proportion or $length/$$queryseqlengths{$tid} < $Length_proportion);
		next if ($both eq "F") and ($length/$$queryseqlengths{$qid} < $Length_proportion and $length/$$queryseqlengths{$tid} < $Length_proportion);
	}
	
	print TMP;

	next if -f $merged_overlaps_file; # i.e. don't use masked strings
	
	($qst, $qen) = ($qen, $qst) if $qst > $qen;
	($tst, $ten) = ($ten, $tst) if $tst > $ten;
	if ($qid =~ /^(.+?)_(\d+)_(\d+)$/) { $qst = $qst + $2 - 1; $qen = $qen + $2 - 1; $qid = $1 }
	if ($tid =~ /^(.+?)_(\d+)_(\d+)$/) { $tst = $tst + $2 - 1; $ten = $ten + $2 - 1; $tid = $1 }

	if (exists $maskstrings{$qid}) { $maskstrings{$qid} .= "0" x ($qen - length($maskstrings{$qid})) }
	else                           { $maskstrings{$qid} .= "0" x  $qen }
	substr($maskstrings{$qid}, $qst - 1, $qen - $qst + 1, "1"  x ($qen - $qst + 1));

	if (exists $maskstrings{$tid}) { $maskstrings{$tid} .= "0" x ($ten - length($maskstrings{$tid})) }
	else                           { $maskstrings{$tid} .= "0" x  $ten }
	substr($maskstrings{$tid}, $tst - 1, $ten - $tst + 1, "1"  x ($ten - $tst + 1));
}
close TMP;

################################################################
# for each qid and tid create array of strings of masked bits (in order, so chrI can have 1_245, 250_260 etc)

my %intervals;
if (-r $merged_overlaps_file)
{
    my $merged_overlaps_fh = &read_fh ($merged_overlaps_file);
    while (<$merged_overlaps_fh>)
    {
        push @{$intervals{$1}}, [$2,$3] if /^(\S+)\t(\d+)\t(\d+)/;
    }
    for my $chrontig (keys %intervals)
    {
        @{$intervals{$chrontig}} = sort { $a->[0] <=> $b->[0] } @{$intervals{$chrontig}};
    }
}
else {
    for my $id (keys %maskstrings)
    {
        my $mask = $maskstrings{$id};
        while ($mask =~ /(1+)/g)
        {
            push @{$intervals{$id}}, [pos($mask) - length($1) + 1, pos($mask)];
        }
    }
}

################################################################
# go through BLAST file again, this time, compare each qid and tid segment
# with array of interval strings for that qid/tid
# if overlapping, attach the id_intervalst_intervalen to the hash created for each entry

my %onewaylinks;
my %twowaylinks;
my $count = 0;
open TMP, "<$tmp_blastfile" or die $!;
while (<TMP>)
{
    chomp;
	next if /^#/;
	my @fields = split /\t/;
	#next unless scalar @fields >= 12;
	my ($qid, $tid, $perc_identity, $length, $gaps, $mismatches, $qst, $qen, $tst, $ten, $evalue, $score) = @fields;
	($qst, $qen) = ($qen, $qst) if $qst > $qen;
	($tst, $ten) = ($ten, $tst) if $tst > $ten;
	if ($qid =~ /^(.+?)_(\d+)_(\d+)$/) { $qst = $qst + $2 - 1; $qen = $qen + $2 - 1; $qid = $1 }
	if ($tid =~ /^(.+?)_(\d+)_(\d+)$/) { $tst = $tst + $2 - 1; $ten = $ten + $2 - 1; $tid = $1 }

    ($qst, $qen) = @{$intervals{$qid}[ &binary_search_intervals ( $intervals{$qid}, $qst, $qen ) - 1 ]};
    ($tst, $ten) = @{$intervals{$tid}[ &binary_search_intervals ( $intervals{$tid}, $tst, $ten ) - 1 ]};
	
	# code for storing perc id and length for each pair removed
	$onewaylinks{"$qid\_$qst\_$qen"}{"$tid\_$tst\_$ten"} = 1;
	$onewaylinks{"$tid\_$tst\_$ten"}{"$qid\_$qst\_$qen"} = 1 unless $twoway;
	print STDERR "$qid\_$qst\_$qen\t$tid\_$tst\_$ten\t$perc_identity\t$length\t$evalue\t$score\n";
#	print STDERR ($count % 1000000 ? ".":$count) unless $count % 100000;$count++;
}
close TMP;
unlink $tmp_blastfile;

for my $node1 (keys %onewaylinks)
{
	for my $node2 (keys %{$onewaylinks{$node1}})
	{
		if (exists $onewaylinks{$node1}{$node2} and exists $onewaylinks{$node2}{$node1})
		{
				$twowaylinks{$node1}{$node2} = $twowaylinks{$node2}{$node1} = 1;
		}
	}
}

%onewaylinks = ();

################################################################
# make clusters based on %twowaylinks

my (%cluster);
my $clusterid = 1;
for my $node (keys %twowaylinks)
{
	next if exists $cluster{$node};
	&visitnode($node); # updates global variable %cluster
	$clusterid++;
}

# each node maps to a clusterid, now map each clusterid to set of nodes
my %cluster2node;
for my $node (keys %cluster) { $cluster2node{$cluster{$node}}{$node} = 1 }

# print contents of each cluster
foreach (sort {scalar keys %{$cluster2node{$b}} <=> scalar keys %{$cluster2node{$a}}} keys %cluster2node) { print join(" ", keys %{$cluster2node{$_}}) . "\n" }

################################################################
################################################################
################################################################

sub visitnode
{
	my $node = shift @_;
	return if exists $cluster{$node};
	$cluster{$node} = $clusterid;
	foreach (keys %{$twowaylinks{$node}})
	{
		next if exists $cluster{$_};
		&visitnode($_)
	}
}

###################################################################################################

sub overlapping
{
	die "overlapping needs 4 arguments\n" unless scalar @_ == 4;
	my ($ast, $aen, $bst, $ben) = @_;
	if ($_[0] > $_[1]) { $ast = $_[1]; $aen = $_[0] }
	if ($_[2] > $_[3]) { $bst = $_[3]; $ben = $_[2] }
	if (($ast > $ben) or ($bst > $aen)) { return 0 }
	my $st = $ast < $bst ? $ast : $bst;
	my $en = $aen > $ben ? $aen : $ben;
	return [$st, $en];
}

###################################################################################################
sub multilinefastafile2singleline
{
	my $fastafile = shift @_;
	unless (`wc -l $fastafile | cut -f1 -d' '` == (`grep -c '>' $fastafile` * 2))
	{
		my $tmpfile = rand();
		open FA,  "<$fastafile" or die "Can't open $fastafile\n";
		open TMP, ">$tmpfile" or die "Can't open $tmpfile for multiline fasta to singleline\n";
		my $i = 0;
		while (<FA>) {
			$i++;
			if (/^>/) { print TMP "\n" unless $i == 1; print TMP }
			elsif (/^\s*$/) { next }
			elsif (/\d/) { chomp; print TMP; print TMP " " }
			else { chomp; print TMP }
		}
		print TMP "\n";
		close FA;
		close TMP;
		rename $tmpfile,$fastafile or die $!
	}
}
###################################################################################################

sub fastafile2sequencelengths
{
	my $fastafile = shift @_;
	my %sequencelengths;
	open FA, "<$fastafile" or die $!;
	while (<FA>)
	{
		next unless /^>(\S+)(.*)/;
		$sequencelengths{$1} = length(<FA>) - 1;
	}
	return \%sequencelengths;
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
#############################################################################

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
