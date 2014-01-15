#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $maf_file = "-";
my $min_length = 10;
GetOptions (
  "m|maf:s"	  => \$maf_file,
  "l|length:i"=> \$min_length,
);

################################################################################################
# open maf file
# read in a block at a time

my $maf_file_fh = &read_fh ($maf_file);
#my $maf_file_fh = \*DATA;

my ($a_line, $block);

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
	# split block into multiple blocks in an array
	&block_split_masked ($block_struct, $min_length);
	
}


#############################################################################################

sub block_split_masked {

# take a block struct and convert into array of blocks (only keep those blocks with at least 
# min_cols

	my $block_struct = shift @_;
	my $min_cols     = shift @_;
	
	my $blockhash    = {%{$block_struct}};

	my $mask = $$blockhash{"m"};
	my @blocks_split = ();

	my $unmask_pos   = 0;

	while ($mask =~ /([^1]*)(1+)/g) {

		my $slice_start = $unmask_pos + length($1);
		my $slice_len   = length($2);
		$unmask_pos    += length($1)  + length($2);
		
		next if $slice_len < $min_cols;
		my $slicehash = {%{$blockhash}};

		# crudely recalc score as a fraction of score based on slice length;
		my $score;
		$$blockhash{"a"} =~ /=(\S+)/ and $score = $1;
		$$slicehash{"a"} = "a score=" . ($score * $slice_len / length($mask));

		my $seq_count = 0;

		$$slicehash{"s"} = {%{$$blockhash{"s"}}};
		
		for my $seqid (sort keys %{$$slicehash{"s"}}) {
		
			$$slicehash{"s"}{$seqid} = {%{$$blockhash{"s"}{$seqid}}};		
			$$slicehash{"s"}{$seqid}{seq} = substr($$slicehash{"s"}{$seqid}{seq}, $slice_start, $slice_len);

			my $pre_dashes   = (substr($$blockhash{"s"}{$seqid}{seq}, 0, $slice_start) =~  tr/\-/\-/);
			my $post_dashes  = (substr($$blockhash{"s"}{$seqid}{seq}, $slice_start + $slice_len) =~  tr/\-/\-/);
			my $slice_dashes = ($$slicehash{"s"}{$seqid}{seq} =~ tr/\-/\-/);

			$$slicehash{"s"}{$seqid}{len}    	= $slice_len - $slice_dashes;

			if ($$slicehash{"s"}{$seqid}{strand} eq "+") {
				$$slicehash{"s"}{$seqid}{start} += $slice_start - $pre_dashes;
			}
			else {
				$$slicehash{"s"}{$seqid}{start}  = $$blockhash{"s"}{$seqid}{start} + $slice_start - $pre_dashes;
			}
			
			# now that all the slicing is done - skip if the entire sequence is --s (delete it from hash).
			delete $$slicehash{"s"}{$seqid} if $$slicehash{"s"}{$seqid}{seq} eq "-" x length($$slicehash{"s"}{$seqid}{seq});
		}
		
		next if scalar keys %{$$slicehash{"s"}} < 2;
		
		# now, check all sequences at the same time.
		# If any cols are entirely dashes, remove them (won't change coords)
		my @array_of_seqs;
		for my $seqid (sort keys %{$$slicehash{"s"}}) {
			push @array_of_seqs, $$slicehash{"s"}{$seqid}{seq};
		}
		my @trans_array_of_seqs = &transpose_array_of_strings(\@array_of_seqs);
		my @deletedashtrans_array_of_seqs;
		foreach (@trans_array_of_seqs) {
			push @deletedashtrans_array_of_seqs, $_ unless /^-+$/;
		}
		my @deletedash_array_of_seqs = &transpose_array_of_strings(\@deletedashtrans_array_of_seqs);
		for my $seqid (sort keys %{$$slicehash{"s"}}) {
			$$slicehash{"s"}{$seqid}{seq} = shift @deletedash_array_of_seqs;
			delete $$slicehash{"s"}{$seqid} if $$slicehash{"s"}{$seqid}{seq} =~ tr/ATGCatgc/ATGCatgc/ <$min_cols;
		}

		next if scalar keys %{$$slicehash{"s"}} < 2;
		
		print $$slicehash{"a"} . "\n";
		for my $seqid (sort keys %{$$slicehash{"s"}}) {
			my $seq_to_print_meta = "s $$slicehash{s}{$seqid}{chrontig} $$slicehash{s}{$seqid}{start} $$slicehash{s}{$seqid}{len} $$slicehash{s}{$seqid}{strand} $$slicehash{s}{$seqid}{total_len} ";
			print $seq_to_print_meta . (" "x (50-length($seq_to_print_meta))) . "$$slicehash{s}{$seqid}{seq}\n";
		}
		print "\n";
	}
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

__DATA__
##maf version=12 scoring=tba.v12
a score=496656.0
s ce.CHROMOSOME_II			6710640	 99 + 15279345 CCGGTCAATCCGTTTCTTCCTGGTAATCCAGGCATTCCTG--------------GAGCTCCAATTGGTCCTGGTGGTCCTTCGTCTCCTGGTGGTCCCTCCTTCAAGTCACCT-------
s as.Scaffold508			  25895	 99 +	297890 CCTGTGAGGCCGTTCCTTCCTGGCAGACCAGGAAGACCAG--------------GTGCACCAATAGGTCCAATTGGCCCTGGGTCACCAGGAGGTCCTTCTTTGAGCTCTCCT-------
s bm.Bmal_supercontig14975	1370572	 99 +  1380737 CCTGTCATTCCATTTCGACCAGGTAGTCCTGGAAGACCAG--------------GAGCACCAACGGGTCCAACTGGCCCTGGGTCTCCAGGCGGGCCTTCCTTAAGTTCACTT-------
s sr.RATTI_contig_74967		 253818	 99 +	417258 CCTGTAAGACCATTTCTTCCTGGTAAACCTGGTACTCCTG--------------GAGCTCCAATTGGTCCATTTGGACCTGGATCTCCTGGTGGACCCTCACGTACTTCTCCT-------
s mh.MhA1_Contig2086		  16040	 99 -	 22022 CCAGTCATTCCATTTCTACCAGGCAATCCAGGCACTCCCG--------------AAGCTCCAATAGGACCAGGCGGTCCAAGATCTCCAGGTGGACCTTCTCTGACTTCTCCT-------
s mi.Minc_Contig4258			542	 99 +	  6189 CCAGTCATTCCATTTCGACCGGGTAATCCAGGCACTCCTG--------------AAGCCCCAATAGGACCAGGCGGTCCAAGATCTCCAGGCGGACCTTCCCTAACTTCTCCt-------
s bx.scaffold00364			 532272	 99 +  1353723 CCGGTCATTCCATTTCGTCCTGGTAGTCCGGGAACTCCTG--------------GTGCTCCAATGGGTCCAATTGGTCCGGCATCTCCAGGTGGGCCTTCTCTCACTTCTCCT-------
s cbg.chrII				   11344647	 99 - 16627154 CCAGTCAATCCATTTCttcctggcaatcctggcattcctg--------------gTGCTCCAATTGGTCCTGGTGGTCCTTCGTCTCCTGGTGGTCCCTCCTTAAGATCTCCT-------
s cbn.Cbre_Contig143		 122588	 99 -	334039 CCAGTTAATCCATTTCTTCCTGGTAATCCTGGCATTCCTG--------------GTGCTCCAATTGGTCCTGGTGGCCCTTCGTCTCCTGGTGGTCCCTCCTTCAAATCTCCT-------
s csp11.Scaffold629		   16929684	 81 - 33334924 cccggtaatcctggcattcctggtgctccaattggtcctg--------------gtggtccttcatctcctggtggtcc------------------ctctttcaaatctcct-------
s csp5.Csp5_scaffold_00043	  88759 113 +	136946 ccggtcaatccatttcttcctggcaatcctggcattcctgnnnnngcattcctggtgctccaattggtcctggtggtccctcgtctcccggtggtcctTCCTTCAGATCTCCT-------
s cr.Crem_Contig6			 692208	 99 -  2010644 CCTGTCAATCCATTTCttcctggcaatcctggcattcctg--------------gTGCTCCAATTGGTCCTGGTGGTCCCTCATCTCCTGGTGGTCCCTCTTTCAAGTCTCCT-------
s cj.Cjap.Contig17051		  30522	 99 -	200117 CCTGTGAGCCCGTTTCGTCCCGGCAATCCAGGCATACCAG--------------GAGCTCCAATCGGTCCTGGTGGTCCCTCGTCTCCAGGTGGTCCCTCCTTAAGATCTCCT-------
s ca.Can_chrRNAPATH197		  11288	 99 -	 13758 CCAGTCAATCCATTTCttcctggcaatcctggcattcctg--------------gTGCACCAACTGGTCCTGGAACTCCTTCATCTCCTGGTGGCCCTTCTTTCAAATCTCCT-------
s pp.Ppa_Contig5			 185981	 99 -  3196617 CCAGGAAGCCCATTGCgaccaggcagaccaggcagaccaT--------------TGGCACCGATGGGACCGGGCAGACCGGAGTCACCGGGTGGTCCCTCTCGCAGTTCTCCC-------
s ts.GL622789				2821812	 99 -  6373445 CCAACAGGACCGTCTCGCCCTGGTAAGCCAGGTGGTCCTG--------------GATAACCAACGGGTCCCGGTTGTCCAGGATCTCCAGGAGGACC-------AGGTTCCATTTTGGAA

a score=712889.0
s ce.CHROMOSOME_II			8284835 132 + 15279345 AGGTGACTGGTGTTCAAAATCCAGGA---CATGCTGTCTCTATTCTTCTTCGTGGATCAAACAAGTTGGTGCTTGAGGAGGCTGACAGATCTATCCACGATGCTCTTTGTGTGATTCGATGCCTGGTTAAGAAGA
s as.Scaffold709			 137964 132 +	224850 aggtaacaggtgtgcaaaacccgggc---caagcagtgtcgaTCCTGCTTCGCGGTTCGAACAAGCTTGTGCTAGATGAGGCGGAGCGTTCGCTGCACGATGCACTTTGTGTACTGAGGTGTTTGGTGAAGAAGA
s bm.Bmal_supercontig14324	  14681 132 -	 17576 AGGTGACTGGTGTGCAAAACCCGGGT---CAAGCAGTTTCAGTTTTGATTCGTGGTTCTAACAAATTAGTATTGGAAGAAGCAGAACGATCAATTCATGATGCGCTTTGCGTTATACGATGTTTGGTAAAGAAGA
s sr.RATTI_contig_74886		  57848 132 -  1087124 AAGTTACAGGAATTCAAAATTCATCT---AATGCTGTTTCTATTCTTCTTCGTGGTTCAAATAAATTGGTACTTGATGAAGCTGAACGTTCAATGCATGATGCTTTATGTGTAATTAGGTGTCTTGTCAAGAAGA
s mh.MhA1_Contig1191		  17757 132 +	 19827 aGGTAACTGGCCTTAAA---TCGAGTGGAAATGCTGTTTCTATTTTAATAAGAGGATCAAGTAAACTTGTTTTGGAAGAAGCGGAACGATCTTTACATGATGCTTTATGTGTTATTCGTTGTTTAGTAAAGAAAA
s mi.Minc_Contig77			  56944 132 +	 57816 AGGTCACTGGTCTTAAA---TCAAGTGGAAATGCTGTTTCAATTTTAATAAGAGGTTCGAGTAAACTTGTTTTGGAAGAAGCTGAACGATCTTTACACGATGCCTTGTGTGTTATACGTTGTTTAGTGAAGAAGA
s bx.scaffold01109			 953219 132 +  1905649 AGGTCACTGGACTACACGCCGCAGGA---AAAGCCGTCTCAATTTTGTTGCGAGGATCGAATAAGTTGGTGTTGGACGAAGCCGATCGTTCTCTCCATGATGCCCTTTGTGTTCTCAGATGCTTGGTCAAAAAGA
s cbg.chrII				   12896229 132 - 16627154 AGGTCACCGGAGTCCAGAATCCAGGT---CATGCCGTGTCAATCCTTCTCCGTGGATCCAACAAGCTGGTGCTCGAAGAGGCCGATAGATCTCTCCATGATGCCCTCTGCGTAATCAGATGTTTGGTAAAGAGAA
s cbn.Cbre_Contig90			 166245 132 -	263337 AAGTCACCGGAGTTCAGAATCCAGGT---CACGCCGTCTCGATCCTTCTCCGCGGATCGAACAAGTTGGTGCTCGAAGAGGCCGACAGATCCCTCCACGATGCTCTCTGTGTTATTCGGTGTCTCGTCAAGAGAA
s csp11.Scaffold629			2595965 132 + 33334924 AGGTCACAGGAGTTCAGAATCCAGGA---CACGCCGTCTCGATTCTCCTTCGCGGATCCAACAAGCTCGTCCTCGAGGAGGCCGATCGTTCCCTCCATGACGCTCTCTGCGTCATCAGATGCCTCGTCAAAAGGA
s csp5.Csp5_scaffold_00320	  47762 132 -	 57910 AGGTTACCGGAGTCCAGAACCCAGGC---CACGCCGTCTCGATTCTTCTCCGTGGATCCAACAAATTGGTGCTTGAAGAGGCCGACAGATCCCTTCACGATGCCCTCTGTGTCATTAGATGTCTGGTCAAGAGAA
s cr.Crem_Contig1			2608979 132 -  3217881 AGGTAACTGGTGTTCAAAACCCAGGA---CATGCCGTCTCAATTCTTCTCCGTGGATCAAACAAGTTGGTTCTCGAGGAGGCTGACAGATCTCTCCATGATGCTCTCTGTGTTATTCGATGCCTCGTCAAGAGAA
s cj.Cjap.Contig3265		   2423 132 +	 55799 AGGTTACTGGAGTACAGAATCCAGGT---CATGCTGTATCGATTCTTCTCCGCGGATCCAACAAACTTGTCCTGGAAGAGGCCGACAGATCTCTTCATGATGCCCTTTGCGTCATCAGATGCCTCGTCAAGAAGA
s ca.Can_chrRNAPATHr21710	  12220 132 -	 16237 AAGTGACCGGCGTGCAAAATCCCGGA---CAAGCTGTTTCGATTCTCCTTCGCGGTTCGAATAAGCTCGTGTTGGAAGAGGCCGATCGATCACTTCACGATGCTTTGTGTGTGATTCGTTGTCTTgtcaagagaa
s pp.Ppa_Contig36			1268317 132 +  1378626 AGGTTACTGATATTCAATCATCTGGT---CCTGCTGTCTCCATTCTTCTAAGAGGATCTAATAAATTGGTATTAGAAGAGGCAGATCGATCTCTACACGATGCACTATGTGTCATTAGATGTCTCGTCAAGAAGA
s ts.GL622788				1404094 132 -  5663831 AATTTACTGGCGTGCAGAATCCAGGT---AAAACAGCGTCAATAGTGCTGCGTGGTTCCAATCGACTTGTGCTGGACGAAGCAGAACGTTCGTTGCATGATGCGCTTGCTGTGCTTCGTTGTTTATACAAGAAGA
