#!/user/bin/perl

use strict;
use warnings;

my $fl1 = shift;
my @qual;my @qd;my @dp;
my @mrs;my @mq;

open (IN, "zcat $fl1|") or die $!;
while (<IN>){
	if (/#/){print "$_";}
	my @a = split /\t/;
	if (/QD=([\d\.]+)/){if($1<20){next;}push @qd, $1;}
	if (/DP=([\d\.]+)/){if($1<5){next;}push @dp, $1;}
	push @qual, $a[5];
	if (/MQRankSum=([\d\.]+)/){push @mrs, $1;}
	if (/MQ=([\d\.]+)/){push @mq, $1;}
	print "$_";
}
