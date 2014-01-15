#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($fasta_file, $maf_file);
GetOptions (
  "f|fasta:s" => \$fasta_file,
  "m|maf:s"   => \$maf_file,
);

################################################################################################
# open fasta file and load into memory

my $sequences = &fastafile2hash ( $fasta_file) ;

################################################################################################
# open maf file
# modify lines beginning with "s "
# for each s line, replace dashed line with substr from fasta file
# print all other lines as is

my $maf_file_fh = &read_fh ($maf_file);

my ($meta_info, $species, $chrontig, $strand, $start, $end, $dashed_seq, $undashed_seq);

while (<$maf_file_fh>) {
    print and next unless /^(s\s+(([^.]+)\.\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\d+)\s+)(\S+)/;
                          # 1=meta 2=chrontig,3=sp 4=st   5=len  6=strand 7=total 8=dashed_seq
	$meta_info = $1;
	$chrontig  = $2;
	$species   = $3;
	$strand    = $6;
	$dashed_seq= $8;
	if ($strand eq "+") {
		$start = $4;
		$end   = $4 + $5;
		$undashed_seq = substr($$sequences{$chrontig}{seq},$start, $end - $start);
	}
	else {
		$start = $7 - $4 - $5;
		$end   = $7 - $4;
		$undashed_seq = revcomp(substr($$sequences{$chrontig}{seq},$start, $end - $start));
	}
	die "$chrontig not found in $fasta_file\n" unless exists $$sequences{$chrontig};
	$dashed_seq =  &replace_dashed_seq_with_masked_seq ($dashed_seq, $undashed_seq);
	$dashed_seq =~ tr/kmrswvyKMRSWVY/nnnnnnnNNNNNNN/; # to remove weird characters in di contigs
	
	print $meta_info . $dashed_seq . "\n";
}


#############################################################################################

sub replace_dashed_seq_with_masked_seq {

	my $dashed_seq = shift @_;
	my $masked_seq = shift @_;

	my ($dashpos, $dashlen, @alldashes);	
	while ($dashed_seq =~ /([^-]*)(\-+)/g) {
		$dashpos=length($1);
		$dashlen=length($2);
		push @alldashes,[$dashpos,$dashlen];
	}
	my $seqpos = 0;
	my $final_seq = "";
	foreach (@alldashes) {
		$dashpos = $_->[0];
		$dashlen = $_->[1];
		$final_seq .= substr($masked_seq,$seqpos,$dashpos) . ("-" x $dashlen);
		$seqpos = $seqpos + $dashpos;
	}
	$final_seq .= substr($masked_seq,$seqpos);
	return $final_seq;
}

#############################################################################################

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

#############################################################################################

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
		    chomp($sequences{$seqid}{seq} .= lc($_)) if $changecase eq "L";
		    chomp($sequences{$seqid}{seq} .= uc($_)) if $changecase eq "U";
		    chomp($sequences{$seqid}{seq} .= $_    ) if $changecase eq "N";
		}
	}
	return \%sequences;
}

#############################################################################################

sub revcomp {
	$_ = shift @_;
	tr/atgcATGC/tacgTACG/;
	return scalar reverse $_;
}
