#!/usr/bin/perl

use strict;
use warnings;

my $fl1 = shift;
my @qual;my @qd;my @dp;
my @mrs;my @mq;
my @af;
open (IN, "zcat $fl1|") or die $!;
while (<IN>){
	next if (/#/);
	my @a = split /\t/;
	if (/QD=([\d\.]+)/){push @qd, $1;}
	if (/DP=([\d\.]+)/){push @dp, $1;}
	push @qual, $a[5];
	if (/MQRankSum=([\d\.]+)/){push @mrs, $1;}
	if (/MQ=([\d\.]+)/){push @mq, $1;}
	if (/AF=([\d\.]+)/){push @af, $1;}
}
$fl1 =~ /^(\w+)\./;
my $p01 = int(@qual * .001);
my $p05 = int(@qual * .005);
my $p1 = int(@qual * .01);
my $p50 = int(@qual * .5);
my $p25 = int(@qual * .25);
my $p75 = int(@qual * .75);
my $p99	= int(@qual * .99);
my $p995 = int(@qual * .995);
my $p999 = int(@qual * .999);

print "Sample\tFeature\tMin\t0.1%\t0.5%\t1%\t25%\t50%\t75%\t99%\t99.5%\t99.9%\tMax\n";
if ($p1 >= 1) {
   @qual = sort {$a <=> $b} @qual;
   print "$1\tqual\t$qual[0]\t$qual[$p01]\t$qual[$p05]\t$qual[$p1]\t$qual[$p25]\t$qual[$p50]\t$qual[$p75]\t$qual[$p99]\t$qual[$p995]\t$qual[$p999]\t$qual[-1]\n";
}

if ($p1 >= 1) {
	@qd = sort {$a <=> $b} @qd;
	print "$1\tQD\t$qd[0]\t$qd[$p01]\t$qd[$p05]\t$qd[$p1]\t$qd[$p25]\t$qd[$p50]\t$qd[$p75]\t$qd[$p99]\t$qd[$p995]\t$qd[$p999]\t$qd[-1]\n";
}

if ($p1 >= 1) {
    @dp = sort {$a <=> $b} @dp;
    print "$1\tDP\t$dp[0]\t$dp[$p01]\t$dp[$p05]\t$dp[$p1]\t$dp[$p25]\t$dp[$p50]\t$dp[$p75]\t$dp[$p99]\t$dp[$p995]\t$dp[$p999]\t$dp[-1]\n";
}

$p01 = int(@mq * .001);
$p05 = int(@mq * .005);
$p1 = int(@mq * .01);
$p50 = int(@mq * .5);
$p25 = int(@mq * .25);
$p75 = int(@mq * .75);
$p99 = int(@mq * .99);
$p995 = int(@mq * .995);
$p999 = int(@mq * .999);
if ($p1 >= 1) {
    @mq = sort {$a <=> $b} @mq;
    print "$1\tMQ\t$mq[0]\t$mq[$p01]\t$mq[$p05]\t$mq[$p1]\t$mq[$p25]\t$mq[$p50]\t$mq[$p75]\t$mq[$p99]\t$mq[$p995]\t$mq[$p999]\t$mq[-1]\n";
}

$p01 = int(@mrs * .001);
$p05 = int(@mrs * .005);
$p1 = int(@mrs * .01);
$p50 = int(@mrs * .5);
$p25 = int(@mrs * .25);
$p75 = int(@mrs * .75);
$p99 = int(@mrs * .99);
$p995 = int(@mrs * .995);
$p999 = int(@mrs * .999);
if ($p1 >= 1) {
    @mrs = sort {$a <=> $b} @mrs;
    print "$1\tMSR\t$mrs[0]\t$mrs[$p01]\t$mrs[$p05]\t$mrs[$p1]\t$mrs[$p25]\t$mrs[$p50]\t$mrs[$p75]\t$mrs[$p99]\t$mrs[$p995]\t$mrs[$p999]\t$mrs[-1]\n";
}

$p01 = int(@af * .001);
$p05 = int(@af * .005);
$p1 = int(@af * .01);
$p50 = int(@af * .5);
$p25 = int(@af * .25);
$p75 = int(@af * .75);
$p99 = int(@af * .99);
$p995 = int(@af * .995);
$p999 = int(@af * .999);
if ($p1 >= 1) {
    @af = sort {$a <=> $b} @af;
    print "$1\tAF\t$af[0]\t$af[$p01]\t$af[$p05]\t$af[$p1]\t$af[$p25]\t$af[$p50]\t$af[$p75]\t$af[$p99]\t$af[$p995]\t$af[$p999]\t$af[-1]\n";
}
