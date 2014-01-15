#!/usr/bin/env perl

@a = @ARGV;
for $a (@a) {
	for $b (@a) {
		print "$a\t$b\n" unless exists $done{$a}{$b};
		$done{$a}{$b} = $done{$b}{$a} = 1;
	}	
}
