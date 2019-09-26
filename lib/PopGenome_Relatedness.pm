package PopGenome_Relatedness;
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

#--------------------------------------Step02.SampleFiltering--------------------------------------
sub RELATEDNESS{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/02.SampleFiltering/Relatedness";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/02.SampleFiltering/Relatedness";
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}
	
	open SH, ">$shpath/Relatedness.sh";
	print SH "cd $outpath;vcftools --gzvcf $cfg{step1}{variant_filtering}{vcf} --relatedness2 --out $outpath/RD\n";
	close SH;

	#switch on the bash running
	`sh $shpath/Relatedness.sh 1>$shpath/Relatedness.sh.o 2>$shpath/Relatedness.sh.e` unless ($skipsh ==1);

	#INDV1	INDV2	N_AaAa	N_AAaa	N1_Aa	N2_Aa	RELATEDNESS_PHI
	#1-1	1-1	3.52698e+06	0	3.52698e+06	3.52698e+06	0.5
	#1-1	1-2	1.55085e+06	156628	3.52698e+06	3.4146e+06	0.178286
	
	open IN, "$outpath/RD.relatedness2";
	open OT, ">$outpath/RD.relatedness2.filtered.txt";
	while (<IN>){
		my @a = split /\t/;
		if($a[6] =~ /RELATEDNESS_PHI/){print OT $_;next;}
		if($a[0] eq $a[1]){next;}
		if($a[6] >= 0.25){print OT $_;} 
	}
	close IN;
	close OT;
}

1;