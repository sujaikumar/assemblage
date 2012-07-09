#!/usr/bin/perl

die "Usage: unshuffleSequences_fastx.pl [2|4] <paired_fasta_or_fastq_file>\n" unless scalar @ARGV == 2;
$lines_per_seq = shift @ARGV;
die "Expecting first argument to be 2 or 4" unless $lines_per_seq =~ /^2|4$/;
$filenameIn = shift @ARGV;

open $INFILE, "<$filenameIn" or die "$filenameIn not readable\n";

open $OUT1, "> $filenameIn\.1" or die "Cannot create $filenameIn\.1\n";
open $OUT2, "> $filenameIn\.2" or die "Cannot create $filenameIn\.2\n";

while(<$INFILE>) {
	print $OUT1 $_;
	for (2..$lines_per_seq) { $_ = <$INFILE>; print $OUT1 $_ };
	for (1..$lines_per_seq) { $_ = <$INFILE>; print $OUT2 $_ };
}
