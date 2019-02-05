#!/root/miniconda3/bin/perl
use strict;
use warnings;
#Transcript:CBG10375	|ATG|	cb25.fpc2310b	380	0
open IN, shift;
my %wd;
my %syn;
my %nonsyn;
my $size=5000;
my $flag = 0;
while (<IN>){
	my @a = split /\t+/;
	my @b;
	my @c;
	$syn{$a[0]}=\@b;
	$nonsyn{$a[0]}=\@c;
	for(my $i=1;$i<$a[1];$i=$i+$size){
		my $end = $i+$size-1;
		my $key = "$a[0]:$i:$end";
		$wd{$key}="NULL";
	}
}
close IN;

open IN, shift;
while(<IN>){
	my @a = split /\t+/;
	$syn{$a[2]}->[$a[3]]=$a[4];
	$nonsyn{$a[2]}->[$a[3]]=1-$a[4];
}
close IN;

foreach (keys %wd){
	my @a = split /:/;
	my $s_syn =0;
	my $s_nonsyn =0;
	for(my $i=$a[1];$i<=$a[2];$i++){
		$s_syn = $s_syn + $syn{$a[0]}->[$i] if (exists $syn{$a[0]}->[$i]);
		$s_nonsyn =  $s_nonsyn + $nonsyn{$a[0]}->[$i] if (exists $nonsyn{$a[0]}->[$i]);
	}
	$s_syn = "NULL" if ($s_syn ==  0);
	$s_nonsyn = "NULL" if ($s_nonsyn ==  0);
	print "$a[0]\t$a[1]\t$a[2]\t$s_syn\t$s_nonsyn\n";
}

