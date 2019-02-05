#!/root/miniconda3/bin/perl
use strict;
use warnings;

open IN,shift;
my %silent;
while(<IN>){
    next if (/CHROM/);
    my @a = split /\t/;
    my $key = "$a[0]:$a[1]:$a[2]";
	$silent{$key}=$a[4];
}
close IN;

my %sites;
open IN, shift;
print "CHROM\tBINSTART\tBINEND\tPISILENT\n";
while(<IN>){
    next if (/CHROM/);
    my @a = split /\t/;
    my $key = "$a[0]:$a[1]:$a[2]";
	my $size = $a[2]-$a[1]+1;
	if(exists $silent{$key}){
		my $pi_silent=$silent{$key}*$size/($size-$a[4]);
		print "$a[0]\t$a[1]\t$a[2]\t$pi_silent\n";
	}
}
close IN;
