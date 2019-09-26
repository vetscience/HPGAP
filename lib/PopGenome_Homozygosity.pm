package PopGenome_Homozygosity;
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
#   	step4_Homozygosity   	#
#			   #
#################################
sub HOMOZYGOSITY{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/Homozygosity/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/Homozygosity/";

	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{high_confidence_vcf};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}

	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\n";
	}

	open SH, ">$shpath/Homozygosity.sh";
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 4);
		open OT, ">$outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		print SH "cd $outpath;vcftools --gzvcf $gzvcf --keep $outpath/$pop_name.list --het --stdout ",'|perl -ne \'if(/INDV/){print "INDV\tO_HOM\tE_HOM\tN_SITES\tPERCENTAGE\tF\n";}else{@a=split /\t/;$p=$a[1]/$a[3]*100;$p=sprintf("%.2f", $p);print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$p\t$a[4]";}\'', ">$pop_name.het.list\n";
	}	
	close SH;

	`sh $shpath/Homozygosity.sh 1>$shpath/Homozygosity.sh.o 2>$shpath/Homozygosity.sh.e` unless ($skipsh ==1);

}

1;