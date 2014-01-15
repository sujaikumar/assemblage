#!/usr/bin/perl

# inspired by map_data_ids in maker-2.25 by Carson Holt/ Yandell Lab
# sujai.kumar@ed.ac.uk 2012-07-21

# Usage: map_data_ids.pl -m twocolmap.file -i input.file >mapped.file
# Takes -m map file (two columns, tab delimited, from<tab>to). loads it into hash
# Takes -i inputfile
# - looks up hash; 
# - replaces every column of -i inputfile where the column is a value in the from column, with the to column
# - columns defined by \\t or -d 'anyregexpforsplitting'\n";

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($map_file, $input_file, $col, $delimiter) = ("","-","all","\t");
GetOptions (
  "map:s"      => \$map_file,
  "input:s"    => \$input_file,
  "col:s"      => \$col,
  "delimiter:s"=> \$delimiter,
);

if (not $map_file) {
    print STDERR "\nUsage: map_data_ids.pl -m twocolmap.file -i input.file [-c colnumber] [-d delimiter] > mapped.file\n
Replaces every column of -i inputfile where the column is a value in the from column, with the to column\n
-m map file (two columns, tab delimited, from<tab>to). loads it into hash
-i inputfile (in which cols are to be replaced)
-c colnumber, optional, default all: 
-d delimiter for splitting, default '\\t'. use -c <number> to specify which col to replace
   or, use -c all to replace in all columns\n";
    exit 1;
}

my @F;

################################################################################################
# open map file and load into memory

my $map_file_fh = &read_fh ($map_file);
my %map;

while (<$map_file_fh>) {
    chomp;
    @F = split /\t/;
    $map{$F[0]} = $F[1];
}

################################################################################################
# open input data file and parse each col

my $input_file_fh = &read_fh ($input_file);

my ($line, $pos, $last);

if ( uc(substr($col,0,1)) eq "A" ) {
    while ($line = <$input_file_fh>) {
    
        $pos = 0;
        chomp($line);
        while ($line =~ /(.*?)($delimiter)/g) {
            print (exists $map { $1 } ? $map { $1 } : $1 );
            print $2;
            $pos = pos($line);
        }
        $last = substr($line, $pos);
        print (exists $map { $last } ? $map { $last } : $last );
        print "\n";
    }
}
else {
    while ($line = <$input_file_fh>) {
        @F = split /$delimiter/, $line;
        $F[$col - 1] = $map { $F[$col - 1] } if exists $map { $F[$col - 1] };
        print join ($delimiter, @F);
    }
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
