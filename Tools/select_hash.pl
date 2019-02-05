#!/usr/bin/perl
open IN, shift;
my %list;
while (<IN>){
	if(/^(\S+)/){
		$list{$1}=1;
	}
}
close IN;

open IN, shift;
while (<IN>){
	if(/^(\S+)/){
		print if (exists $list{$1});
	}
}
close IN;
