#!/usr/bin/perl

use strict;
use warnings;

# takes as input file with:
# chrontigname_st_en or chrontigname\tst\ten, and gives back chrontig with merged intervals:

my %maskstrings;

while (<>) {
    next unless /^>?\s*(\S+)[\t_ ](\d+)[\t_ ](\d+)\s*$/;
    
    my ($chrontig, $st, $en) = ($1, $2, $3);
    ($st, $en) = ($en, $st) if $st > $en;
    
	if (exists $maskstrings{$chrontig}) {
	    $maskstrings{$chrontig} .= "0" x ($en - length($maskstrings{$chrontig}))
	}
	else {
	    $maskstrings{$chrontig} .= "0" x  $en
	}
	substr($maskstrings{$chrontig}, $st - 1, $en - $st + 1, "1"  x ($en - $st + 1));
}

for my $id (keys %maskstrings)
{
   my $mask = $maskstrings{$id};
   while ($mask =~ /(1+)/g)
   {
       print $id . "\t" . (pos($mask) - length($1) + 1) . "\t" . pos($mask) . "\n";
   }
}
