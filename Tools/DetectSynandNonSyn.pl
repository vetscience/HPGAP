#!/root/miniconda3/bin/perl
use warnings;
use strict;
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use Bio::SeqIO;

my $color = $ARGV[0];
my $outfile = $color."_syn_nonsyn.txt";
open(OUTFILE, '>', $outfile) or die "Could not open $outfile \n";
my $filename = $ARGV[1];

my %codon_table = (
'ttt' => 'F', 'ttc' => 'F', 'tta' => 'L', 'ttg' => 'L',
'ctt' => 'L', 'ctc' => 'L', 'cta' => 'L', 'ctg' => 'L',
'att' => 'I', 'atc' => 'I', 'ata' => 'I', 'atg' => 'M',
'gtt' => 'V', 'gtc' => 'V', 'gta' => 'V', 'gtg' => 'V',
'tct' => 'S', 'tcc' => 'S', 'tca' => 'S', 'tcg' => 'S',
'cct' => 'P', 'ccc' => 'P', 'cca' => 'P', 'ccg' => 'P',
'act' => 'T', 'acc' => 'T', 'aca' => 'T', 'acg' => 'T',
'gct' => 'A', 'gcc' => 'A', 'gca' => 'A', 'gcg' => 'A',
'tat' => 'Y', 'tac' => 'Y', 'taa' => 'Stop', 'tag' => 'Stop',
'cat' => 'H', 'cac' => 'H', 'caa' => 'Q', 'cag' => 'Q',
'aat' => 'N', 'aac' => 'N', 'aaa' => 'K', 'aag' => 'K',
'gat' => 'D', 'gac' => 'D', 'gaa' => 'E', 'gag' => 'E',
'tgt' => 'C', 'tgc' => 'C', 'tga' => 'Stop', 'tgg' => 'W',
'cgt' => 'R', 'cgc' => 'R', 'cga' => 'R', 'cgg' => 'R',
'agt' => 'S', 'agc' => 'S', 'aga' => 'R', 'agg' => 'R',
'ggt' => 'G', 'ggc' => 'G', 'gga' => 'G', 'ggg' => 'G'
);


my @alignments = glob("*.multi.cds.best.nuc.fas");
foreach my $algn (sort @alignments) {
	$algn =~ /^(.+).multi.cds.best.nuc.fas/;
	my $pc = $1;
	print OUTFILE $pc."\t";
	my $seqio  = Bio::SeqIO->new(-format => 'fasta', -file   => $algn);
	my @sequences;
	while (my $seqobj = $seqio->next_seq) {
		my $nuc = $seqobj->seq();
		push(@sequences, $nuc);
	}
	my $length = length($sequences[0]);
	my $codon_no = $length/3;
	my $synonymous = 0;
	my $nonsynonymous = 0;
	my $consensus_amino_acid;
	foreach my $i (0 ... $codon_no) {
		my $k = $i*3;
		my $consensus = lc(substr($sequences[0],$k,3));
		next if ($consensus =~ /-/);
		chomp $consensus;
		$consensus_amino_acid = $codon_table{$consensus};
		foreach my $seq (@sequences){
			my $codon = lc(substr($seq,$k,3));
			chomp $codon;
			next if ($codon =~ /-/);
			my $amino_acid = $codon_table{$codon};
			if ($codon ne $consensus) {
				if ($amino_acid eq $consensus_amino_acid) {
				$synonymous++;
				}
				if ($amino_acid ne $consensus_amino_acid) {
				$nonsynonymous++;
				}
			}
		}
	}
	print OUTFILE $synonymous."\t".$nonsynonymous."\n";
}
close OUTFILE;
