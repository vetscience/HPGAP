#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;
my $window=shift;
my $length=shift;

open IN, $file;
my $st=1;
my $ed=$window;
my $n=0;
my $sum=0;
my $r2=0;
while(<IN>){
	if ($ed > $length){last;}
	if (/CHR/){next;}
	my @a=split /\s+/;
	if ($a[1]>$ed){
		if ($n <= 2){$r2="nan"}else{$r2=$sum/$n;}
		print "$a[0]\t$st\t$ed\t$r2\t$n\n" if ($ed <= $length);
		$sum=0;
		$n=0;
		$st+=$window;
		$ed+=$window;
		while(1){
			if ($a[1]>$ed){
				if ($n <= 2){$r2="nan"}else{$r2=$sum/$n;}
				print "$a[0]\t$st\t$ed\t$r2\t$n\n" if ($ed <= $length);
				$st+=$window;
				$ed+=$window;
			}else{last;}
		}
	}
	if(($a[1]<=$ed)&&($a[1]>=$st)){
		if (($a[2]<=$ed)&&($a[2]>=$st)){
			$n++;
			$sum+=$a[4];
		}
	}
}
