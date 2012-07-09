#!/usr/bin/perl
while (<>) {
	if (/^\s*$/) {next}
	elsif (/^>/) {print "\n" unless ++$i == 1; print}
	elsif (/\d/){ chomp; print; print " " }
	else { chomp; print }
}
print "\n"