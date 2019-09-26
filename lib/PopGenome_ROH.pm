package PopGenome_ROH;
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
#   	step4_ROH			   	#
#			   #
#################################
sub ROH{

	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/ROH/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/ROH/";
	my $gzvcf = $cfg{step1}{variant_filtering}{high_confidence_vcf};
	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}

	$cfg{step4}{ROH}{windowsize} = 100 unless (defined $cfg{step4}{ROH}{windowsize});
	$cfg{step4}{ROH}{scaffold_number_limit} = 95 unless (defined $cfg{step4}{ROH}{scaffold_number_limit});
	$cfg{step4}{ROH}{scaffold_length_cutoff} = 1000000 unless (defined $cfg{step4}{ROH}{scaffold_length_cutoff});
	my $window_size = $cfg{step4}{ROH}{windowsize};
	my $scaffold_number_limit = $cfg{step4}{ROH}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{ROH}{scaffold_length_cutoff};

	### load reference and map scaffold to number list 
	my $i=1;my $j=0;
	open OT, ">$outpath/chr_map.list";
	foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
		print OT "$id\t$i\n";
		$i++;
		if (($genome->{len}{$id} >= $cfg{step4}{ROH}{scaffold_length_cutoff})&&($j < $cfg{step4}{ROH}{scaffold_number_limit})){
			$j++;
		}
	}
	close OT;

	$cfg{step4}{ROH}{scaffold_number_limit} = $j;

	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\n";
	}

	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 4);
		open SH, ">$shpath/ROH.sh";
		open OT, ">$outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;
		print SH "cd $outpath\n";
		print SH 'bcftools annotate --threads 10 --rename-chrs', " $outpath/chr_map.list" ," $gzvcf", '|perl -ne \'if (/#\S+ID=(\d+)/){if($1<=',"$scaffold_number_limit",'){print;}}elsif(/^#/){print;}elsif(/^(\d+)\s+/){if($1<= ',"$scaffold_number_limit",'){print;}}\'', ">$outpath/$pop_name.pre.vcf\n"; 
		print SH "vcftools --vcf $outpath/$pop_name.pre.vcf --keep $outpath/$pop_name.list --plink --out $pop_name.plink \n";
		print SH "plink --file $pop_name.plink --make-bed --chr-set $scaffold_number_limit no-xy no-mt no-y --out $pop_name.plink\n";
		print SH "plink --bfile $outpath/$pop_name.plink --chr-set $scaffold_number_limit no-xy no-mt no-y --homozyg --homozyg-window-het 2 --homozyg-kb $window_size --homozyg-gap 50 --homozyg-window-snp 100 --homozyg-window-missing 5 --out $pop_name\n";
		close SH; 
	}
	`sh $shpath/ROH.sh 1>$shpath/ROH.sh.o 2>$shpath/ROH.sh.e` unless ($skipsh ==1);
}

1;