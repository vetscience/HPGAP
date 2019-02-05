#!/root/miniconda3/bin/perl
use strict;
use warnings;

open IN,shift;
my %syn;
while(<IN>){
    next if (/CHROM/);
    my @a = split /\t/;
    my $key = "$a[0]:$a[1]:$a[2]";
	$syn{$key}=$a[4];
}
close IN;

open IN,shift;
my %nonsyn;
while(<IN>){
    next if (/CHROM/);
    my @a = split /\t/;
    my $key = "$a[0]:$a[1]:$a[2]";
	$nonsyn{$key}=$a[4];
}
close IN;

my %sites;
open IN, shift;
print "CHROM\tBINSTART\tBINEND\tPIS\tPIN\tPINPIS\n";
while(<IN>){
    next if (/CHROM/);
    my @a = split /\t/;
    my $key = "$a[0]:$a[1]:$a[2]";
	my $size = $a[2]-$a[1]+1;
	if((exists $syn{$key}) && (exists $nonsyn{$key})){
		my $pi_syn = $syn{$key}*$size/$a[3];
		my $pi_nonsyn = $nonsyn{$key}*$size/$a[4];
		my $ratio = $pi_nonsyn/$pi_syn;
		print "$a[0]\t$a[1]\t$a[2]\t$pi_syn\t$pi_nonsyn\t$ratio\n";
	}
}
close IN;
