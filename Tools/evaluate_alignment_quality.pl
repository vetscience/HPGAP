#!/usr/bin/perl
use strict;
use warnings;
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use Bio::SeqIO;

my @files = glob ("*.nuc.fas");
foreach my $file (sort @files){
	my $gene = $1 if ($file =~ /^(\S+)\.multi/);
	my $seqio  = Bio::SeqIO->new(-format => 'fasta', -file   => $file);
	my $flag = 0;
	while (my $seqobj = $seqio->next_seq) {
		my $aln = $seqobj->seq();
		my $n_gaps = () = ($aln =~ /-/g);
		if (($n_gaps*2) > (length $aln)){$flag++;}
	}
	print "$gene\n" if ($flag == 0);
}
