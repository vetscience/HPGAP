#!/usr/bin/perl
use strict;
use warnings;

open IN, shift;

my $a;my $b;my $c;
while(<IN>){
	my @a = split /\s+/;
	if($a[0] =~ /#/){
		for (my $i = 0; $i<@a;$i++){
			$a = $i if ($a[$i] =~ /GeneName/);
			$b = $i if ($a[$i] =~ /variants_effect_missense_variant/);
			$c = $i if ($a[$i] =~ /variants_effect_synonymous_variant/);
		}
	}else {
		print "$a[$a]\t$a[$b]\t$a[$c]\n";
	}
}
close OTP;
close IN;
