package PopGenome_Phylogeny;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use PopGenome_Shared;
use YAML::Tiny;
use Bio::SeqIO;

#################################
#			   #
#   	Step3_phylogeny    		#
#			   #
#################################
sub PHYLOGENY{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/03.GeneticRelationships/Phylogeny/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/03.GeneticRelationships/Phylogeny/";
	my $gzvcf = $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf};
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/phylogeny.sh\n";
	print SH "cd $outpath\n";
	print SH "rm -rf $outpath/RAxML_*\n";
	print SH "zcat $gzvcf|vcf-to-tab >$outpath/plink_data.tab\n";
	print SH "vcf_tab_to_fasta_alignment.pl -i $outpath/plink_data.tab > $outpath/plink_data.fasta\n";
	print SH "seqmagick convert $outpath/plink_data.fasta $outpath/plink_data.phy\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $outpath/plink_data.phy -n trees -T 24 -# 20 -p 12345\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $outpath/plink_data.phy -n boots -T 24 -# 100 -p 23456 -b 23456\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -f b -t RAxML_bestTree.trees -T 2 -z RAxML_bootstrap.boots -n consensus\n";
	print SH "sumtrees.py --percentages --min-clade-freq=0.50 --target=RAxML_bestTree.trees --output=result2.tre RAxML_bootstrap.boots";
	close SH;

	`sh $shpath/phylogeny.sh 1>$shpath/phylogeny.sh.o 2>$shpath/phylogeny.sh.e` unless ($skipsh ==1);
}

1;