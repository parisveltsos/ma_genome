#!/usr/bin/perl
use warnings;
use strict;

my $len = 0;
my $in = $ARGV[0];

open(IN, "<", $in) or die $!;
while(<IN>){
    chomp;
    if(/^[\>\@]/){
	if($len > 0){
	    print $len."\n";
	}print $_."\t";
	$len=0;
    }
    else{
	s/\s+//g;
	$len+=length($_);
    }
}close(IN);

print $len."\n";
