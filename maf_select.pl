#!/usr/bin/env perl

# select those blocks (and maybe report some stats on stderr)
# that have sequences from the species listed in -s on the command line
# eg: maf_select.pl -m tba.20120324.20sp.lc.remove.maf.gz -s ce cbg cr

# Other options
# -i include only the specified species in the final MAF output
# -iff EXACTLY these species have to be present (none extra)
# -l length of alignment (including dashes)

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use List::Util qw(max);

my $maf_file = "-";
my @species_array;
my $min_length   = 30;
my $min_identity = 0; # perfectly identical cols must be at least this proportion
my $min_rel_id   = 0; # relative identity across all cols must be at least this proportion

GetOptions (
	"m|maf:s"	     => \$maf_file,
	"s|species=s{,}" => \@species_array,
	"l|length=i"     => \$min_length,
	"i|identity=f"   => \$min_identity,
	"r|rel_id=f"     => \$min_rel_id, 
);

################################################################################################
# open maf file
# read in a block at a time

my $maf_file_fh = &read_fh ($maf_file);
#my $maf_file_fh = \*DATA;

my $block;

while (my $maf_line = <$maf_file_fh>) {
	print $maf_line and next unless $maf_line =~ /^a score/;
	my $seq_count = 0;
	$block = $maf_line;
	while ($maf_line = <$maf_file_fh>) {
		last if $maf_line =~ /^\s*$/;
		if ($maf_line =~ /^s\s+/) {
			$seq_count++;
			$block .= $maf_line;
		}
	}

	next if $seq_count <2;
	
	# turn block into data struct
	my $block_struct = &block_to_struct ($block);
	
	# do something with the block
	# select only those species in block that are in @species_list
	&block_select_species ($block_struct, \@species_array, $min_length, $min_identity, $min_rel_id);
	
}


#############################################################################################

sub block_select_species {

	my %block_hash   = %{shift @_};
	my @species_array= @{shift @_};
	my $min_length   = shift @_;
	my $min_identity = shift @_;
	my $min_rel_id   = shift @_;

	my %species_list;
	foreach (@species_array) { $species_list{$_} = 1 }

	my $block = $block_hash{"a"} . "\n";
	my %species_found;
	my @seqid_list = sort { index(join(",",@species_array),&get_species_from_seqid($a)) <=> index(join(",",@species_array),&get_species_from_seqid($b)) } keys %{$block_hash{"s"}};

	for my $seqid (@seqid_list) {
		$species_found{$block_hash{"s"}{$seqid}{species}} = 1 if exists $species_list{$block_hash{"s"}{$seqid}{species}};
		my $seq_to_print_meta = "s $block_hash{s}{$seqid}{chrontig} $block_hash{s}{$seqid}{start} $block_hash{s}{$seqid}{len} $block_hash{s}{$seqid}{strand} $block_hash{s}{$seqid}{total_len} ";
		$block .= $seq_to_print_meta . (" "x (50-length($seq_to_print_meta))) . "$block_hash{s}{$seqid}{seq}\n";
	}

	return if join(".",sort(keys %species_list)) ne join(".",sort(keys %species_found)) or
	    length($block_hash{"m"}) < $min_length;
	
	# calculate length and absolute identity of each alignment
	# absolute identity is defined as ID cols/TOTAL cols
	# ID cols are those alignment columns where all sequences are identical

	my @array_of_seqs;
	for my $seqid (sort keys %{$block_hash{"s"}}) {
		push @array_of_seqs, $block_hash{"s"}{$seqid}{seq};
	}
	
	my @trans_array_of_seqs = &transpose_array_of_strings(\@array_of_seqs);
	
	my $identical_cols = 0;
	my $consensus_identity = 0;

	my ($a,$c,$g,$t,$d);
	foreach (@trans_array_of_seqs) {
		$identical_cols++ if $_ eq substr($_,0,1) x length($_) and substr($_,0,1) ne "-";
		$a = tr/Aa/Aa/;
		$c = tr/Cc/Cc/;
		$g = tr/Gg/Gg/;
		$t = tr/Tt/Tt/;
		$d = tr/-/-/;
		$consensus_identity += max($a,$c,$g,$t) if max($a,$c,$g,$t) >length($_)/2;
	}

	my $block_length = scalar @trans_array_of_seqs;
	return if $identical_cols/$block_length < $min_identity;
	return if $consensus_identity/($block_length * (scalar @array_of_seqs)) < $min_rel_id;

	print $block . "\n";	
	print STDERR "$block_length\t" . $identical_cols/$block_length . "\t" . $consensus_identity/($block_length * (scalar @array_of_seqs)) . "\n";
}

#############################################################################################

sub block_to_struct {

# read into hash struct
# convert all lower case into "x"

	my @block = split("\n",shift @_);
	my %blockhash;
	$blockhash{"a"} = $block[0];

	my $mask = "";

	for (my $i = 1; $i <= $#block; $i++) {
		my ($s, $chrontig, $start, $len, $strand, $total_len, $seq) = split(/\s+/, $block[$i]);
		$seq =~ tr/a-z/x/;
		my $seqid = "$chrontig\_$start\_$len\_$strand\_$total_len";
		$blockhash{"s"}{$seqid}{chrontig} = $chrontig;
		$blockhash{"s"}{$seqid}{species}  = $1 if $chrontig =~ /^(.+?)\./; 
		$blockhash{"s"}{$seqid}{start}	  = $start;
		$blockhash{"s"}{$seqid}{len}	  = $len;
		$blockhash{"s"}{$seqid}{strand}	  = $strand;
		$blockhash{"s"}{$seqid}{total_len}= $total_len;
		$blockhash{"s"}{$seqid}{seq}	  = $seq;
#		print "$seq\t$seqid\n";
		$mask = "1" x length($seq) unless $mask; #only needs to be initialised the first time
		my $mask_pos = 0;
		while ($seq =~ /([^x]*)(x+)/g) {
			$mask_pos += length($1);
			substr($mask, $mask_pos,  length($2), "x" x length($2));
			$mask_pos += length($2);
		}
	}
#	print $mask . "\n\n";
	$blockhash{"m"} = $mask;
	return \%blockhash;
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

sub revcomp {
	$_ = shift @_;
	tr/atgcATGC/tacgTACG/;
	return scalar reverse $_;
}

#############################################################################################

sub transpose_array_of_strings {
	my $in_aref = shift @_;
	my @new_array;

	for (my $i=0; $i < length($in_aref->[0]); $i++) {
		my $new_row = join("",map { substr($_,$i,1) } @$in_aref);
		push @new_array, $new_row;
	}
	return @new_array;	
}

#############################################################################################

sub get_species_from_seqid {
	my $seqid = shift @_;
	$seqid =~ /^([^\.]+)/ and return $1;
	return "";
}
