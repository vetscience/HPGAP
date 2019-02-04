#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use FindBin '$Bin';
use lib "$Bin/lib";
use PopGenome;
use YAML::Tiny;

VARIANT_FILTERING(shift);
sub VARIANT_FILTERING{
	my $yml_file = shift;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	
	my $outpath = "$cfg{args}{outdir}/01.QualityControl/";
	my $shpath = "$cfg{args}{outdir}/NewShell/01.QualityControl/step1g";
	my $reference = $cfg{ref}{db}{$cfg{ref}{choose}}{path};

	#my ($samplelist,$reference,$outpath,$shpath,$main,$mountpath,$image)=@_;
	#my %samplelist = %{$samplelist};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/variant.filtering.sh";
	print SH "gatk SelectVariants \\\n";
	print SH "	-R $reference \\\n";
	print SH "	-V $outpath/Combined.HC.vcf.gz \\\n";
	print SH "	--select-type-to-include SNP \\\n";
	print SH "	-O $outpath/Combined_raw_snps1st.vcf && echo \"** GVCF Combined_raw_snps1st done\" && \\\n";
	print SH "gatk VariantFiltration \\\n";
	print SH "	-R $reference \\\n";
	print SH "	-V $outpath/Combined_raw_snps1st.vcf \\\n";
	print SH "	--filter-expression \"$cfg{step1}{variant_filtering}{snp}\" \\\n";
	print SH "	--filter-name \"my_snp_filter\" \\\n";
	print SH "	-O $outpath/Combined_filtered_snps1st.vcf && echo \"** GVCF Combined_raw_snps1st done\" \n";

	print SH "gatk SelectVariants \\\n";
	print SH "	-R $reference \\\n";
	print SH "	-V $outpath/Combined.HC.vcf.gz \\\n";
	print SH "	--select-type-to-include INDEL \\\n";
	print SH "	--maxIndelSize 60 \\\n";
	print SH "	-O $outpath/Combined_raw_indels1st.vcf && echo \"** GVCF Combined_raw_snps1st done\" && \\\n";
	print SH "gatk VariantFiltration \\\n";
	print SH "	-R $reference \\\n";
	print SH "	-V $outpath/Combined_raw_indels1st.vcf \\\n";
	print SH "	--filter-expression \"$cfg{step1}{variant_filtering}{indel}\" \\\n";
	print SH "	--filter-name \"my_indel_filter\" \\\n";
	print SH "	-O $outpath/Combined_filtered_indels1st.vcf && echo \"** Combined_raw_snps1st done\" \n";
	print SH "vcftools --vcf $outpath/Combined_filtered_snps1st.vcf --remove-filtered-all --recode --recode-INFO-all --stdout \| bgzip -c > $outpath/PASS.SNP.vcf.gz";
	#next step is to select the snp sites with high genotype likelihood in all samples.
	close SH;
	`echo hi 1>echo.o 2>echo.e`;
	my $mount="";
	for (my $i=0;$i<@{$cfg{args}{mount}};$i++){$mount .= "-v $cfg{args}{mount}->[$i]->{hostpath}:$cfg{args}{mount}->[$i]->{dockerpath} ";}
	print "udocker run $mount--env=\"$cfg{args}{env}\" $cfg{args}{container}\n";
#	open MH, ">>$main"; print MH "docker run -ti  -v $mountpath --rm $image sh -c 'sh $shpath/variant.filtering.sh 1>$shpath/variant.filtering.sh.o 2>$shpath/variant.filtering.sh.e'\n";close MH;
}
