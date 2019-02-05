#!/usr/bin/perl
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use Bio::PopGen::Statistics;

my $polylist=$ARGV[0];
my $divlist=$ARGV[1];

open IN, $polylist;
my %poly;
while(<IN>){
	chomp;
	my @a = split /\t/;
	$poly{$a[0]}{Pn}=$a[1];
	$poly{$a[0]}{Ps}=$a[2];
}
close IN;

open IN, $divlist;
#print "Gene name\tDn\tPn\tDs\tPs\n";
my %div;
while(<IN>){
	chomp;
	my @a = split /\t/;
	
	#Gene name Dn Pn Ds Ps
	#print "$a[0]\t$a[2]\t$poly{$a[0]}{Pn}\t$a[1]\t$poly{$a[0]}{Ps}\n";	
	next unless (exists $poly{$a[0]}{Pn});
	next unless (exists $poly{$a[0]}{Ps});

	$div{$a[0]}{Dn}=$a[2];
	$div{$a[0]}{Ds}=$a[1];
	next unless ($div{$a[0]}{Dn} >0);
	next unless ($div{$a[0]}{Ds} >0);	
	
	my @stat;
	$stat[0]=$poly{$a[0]}{Pn};
	$stat[1]=$div{$a[0]}{Dn};
	$stat[2]=$poly{$a[0]}{Ps};
	$stat[3]=$div{$a[0]}{Ds};
	my $MK = Bio::PopGen::Statistics-> mcdonald_kreitman_counts(@stat);
	
	if ($MK >= 0.05) {print $a[0]."\t".$MK."\t"."Neutral"."\n";}
	else{
		my $percentage_div = $div{$a[0]}{Dn}/($div{$a[0]}{Dn}+$div{$a[0]}{Ds}+0.01);
		my $percentage_poly = $poly{$a[0]}{Pn}/($poly{$a[0]}{Pn}+$poly{$a[0]}{Ps}+0.01);
		my $flag = "Neutral";
		my $value = $percentage_div - $percentage_poly;
		if ($percentage_div > $percentage_poly) {$flag = "Positive"};
		if ($percentage_poly > $percentage_div) {$flag = "Negative"};
		print $a[0]."\t".$MK."\t".$flag."\t"."$value"."\t"."Dn:$div{$a[0]}{Dn}\tPn:$poly{$a[0]}{Pn}\tDs:$div{$a[0]}{Ds}\tPs:$poly{$a[0]}{Ps}"."\n";
	}
}
close IN;

