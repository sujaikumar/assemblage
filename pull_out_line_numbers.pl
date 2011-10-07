#!/usr/bin/perl

use strict;
use warnings;

die "\n Usage: pull_out_line_numbers.pl <filewithlinenumbers> <filewithlines>\n" .
    "\n <filewithlinenumbers> should be sorted: use <(sort -n file) if not sorted\n" . 
    " <filewithlinenumbers> assumes line numbering begins at 1\n\n" unless scalar @ARGV == 2;

my $line_numbers_file = shift @ARGV;
my $text_file = shift @ARGV;

open LINE, "<$line_numbers_file" or die $!;
open TEXT, "<$text_file"         or die $!;

my $current_text = 0;

while (my $current_line = <LINE>) {
    die "<filewithlines> finished before <filewithlinenumbers>\n" if eof TEXT;
    die "<filewithlinenumbers> can only have one positive integer per line\n" if $current_line !~ /^[1-9]\d*$/;
    while (my $line = <TEXT>) {
        $current_text++;
        if ($current_text == $current_line) {
            print $line;
            last;
        }
    }
}
