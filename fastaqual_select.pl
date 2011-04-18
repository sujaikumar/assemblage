#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($fastafile,$sort,$numfasta,$length,$prefix,$regexp,$headerfile,$includefile,$excludefile,$delimiter,$case) = ("-","S","","","","","","",""," ","");
GetOptions (
  "fastafile:s" => \$fastafile,
  "sort:s" => \$sort,
  "numfasta:i" => \$numfasta,
  "length:i" => \$length,
  "prefix:s" => \$prefix,
  "regexp:s" => \$regexp,
  "includefile:s" => \$includefile,
  "excludefile:s" => \$excludefile,
  "delimiter:s" => \$delimiter,
  "case:s" => \$case,
);

$case = uc(substr($case,0,1)) if $case;

my (%include_headers, %exclude_headers);
if ($includefile)
{
    open FILE,"<$includefile" or die "Couldn't open includes file $includefile\n";
    while (<FILE>)
    {
        chomp;
        if (/^>?(\S+)\s(\d+)\s(\d+)\b/ or /^>?(\S+)_(\d+)_(\d+)\b/)
        {
            push @{$include_headers{$1}}, ($2 <= $3) ? [$2,$3] : [$3,$2];
        } elsif (/^>?(\S+)/)
        {
            $include_headers{$1}=1
        }
    }
    close FILE;
}
if ($excludefile)
{
    open FILE,"<$excludefile" or die "Couldn't open excludes file $excludefile\n";
    while (<FILE>)
    {
        if (/^>?(\S+)/)
        {
            $exclude_headers{$1}=1
        }
    }
    close FILE;
}

my $seqs = &fastafile2hash($fastafile,"N",$sort);
my @sortkeys;

if (uc($sort) eq "R") {
    @sortkeys = sort {$$seqs{$a}{order} <=> $$seqs{$b}{order}} keys %{$seqs};
} elsif (uc($sort) eq "A") {
    @sortkeys = sort keys %{$seqs};
} else {
    @sortkeys = sort {length($$seqs{$b}{seq}) <=> length($$seqs{$a}{seq})} keys %{$seqs};
}

my $num_printed = 0;
for my $chrontig (@sortkeys)
{
    last if $numfasta and $numfasta == $num_printed;
    if ($$seqs{$chrontig}{seq}=~ /\d+/) {
        # it's a qual file
        my @qual = split(/\s+/,$$seqs{$chrontig}{seq});
        next if $length and scalar @qual < $length;
    }
    else {
        next if $length and length($$seqs{$chrontig}{seq}) < $length;
    }   
    next if $excludefile and exists $exclude_headers{$chrontig};
    next if $includefile and not exists $include_headers{$chrontig};
    my $toprint;
    if (%include_headers and $include_headers{$chrontig} =~ /array/i)
    {
        for my $interval (@{$include_headers{$chrontig}})
        {
            $toprint  = ">$prefix$chrontig";
            my ($start, $stop) = @{$interval};
            $toprint .= "$delimiter$start$delimiter$stop";
            $toprint .= " $$seqs{$chrontig}{desc}" if $$seqs{$chrontig}{desc};
            $toprint .= "\n" . substr($$seqs{$chrontig}{seq}, $start - 1, $stop - $start + 1) . "\n";
            next if $regexp and $toprint !~ /$regexp/;
            print $toprint;
            $num_printed++;
        }
    }
    else 
    {
        $toprint  = ">$prefix$chrontig";
        $toprint .= " $$seqs{$chrontig}{desc}" if $$seqs{$chrontig}{desc};
        $$seqs{$chrontig}{seq} = lc($$seqs{$chrontig}{seq}) if $case eq "L";
        $$seqs{$chrontig}{seq} = uc($$seqs{$chrontig}{seq}) if $case eq "U";
        $$seqs{$chrontig}{seq} =~ tr/[A-Z][a-z]/[a-z][A-Z]/ if $case eq "I";
        $toprint .= "\n$$seqs{$chrontig}{seq}\n";
        next if $regexp and $toprint !~ /$regexp/;
        print $toprint;
        $num_printed++;
    }
}

#############################################################################

sub fastafile2hash
{
	my $fastafile  = shift @_;
	my $changecase = "N";
	my $order      = "S"; # S = same as input, or R = random
	$changecase    = substr(uc(shift @_),0,1) if @_;
	$order         = substr(uc(shift @_),0,1) if @_;
	my %sequences;
	my $fh = &read_fh($fastafile);
	my $seqid;
	my $seq_counter;
	while (<$fh>)
	{
		if (/^>(\S+)(.*)/) {
		    $seqid = $1;
		    $sequences{$seqid}{desc} = $2;
		    $sequences{$seqid}{order} = $order eq "S" ? $seq_counter++ : rand;
		}
		else {
		    if (/\d/) {
		        chomp($sequences{$seqid}{seq} .= " $_"); # add space to sep qual values
		        $sequences{$seqid}{seq} =~ s/^\s+//;
		        $sequences{$seqid}{seq} =~ s/\s+$//;
		        next;
		    }
		    chomp($sequences{$seqid}{seq} .= lc($_)) if $changecase eq "L";
		    chomp($sequences{$seqid}{seq} .= uc($_)) if $changecase eq "U";
		    chomp($sequences{$seqid}{seq} .= $_    ) if $changecase eq "N";
		}
	}
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
