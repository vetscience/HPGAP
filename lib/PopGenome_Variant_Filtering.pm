package PopGenome_Variant_Filtering;
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
#    Step 1g Variant filtering  #
#			   #
#################################
sub VARIANT_FILTERING{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	
	my $outpath = "$cfg{args}{outdir}/01.QualityControl/Combined";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/step1_variant_filtering";
	my $reference = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	my $genome=LOADREF($reference);

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/variant.filtering.s1.sh";
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
	print SH "vcftools --vcf $outpath/Combined_filtered_snps1st.vcf --remove-filtered-all --recode --recode-INFO-all --stdout \| bgzip -c > $outpath/PASS.SNP.vcf.gz\n";

	close SH;

	#switch on the bash running
	`sh $shpath/variant.filtering.s1.sh 1>$shpath/variant.filtering.s1.sh.o 2>$shpath/variant.filtering.s1.sh.e` unless ($skipsh ==1);
	
	###filter SNP by depth if needed
	my $sample_size = keys %{$cfg{fqdata}};
	open (IN, "zcat $outpath/PASS.SNP.vcf.gz|") or die $!;
	my @dp;
	while (<IN>){
    	next if (/#/);
    	if (/DP=([\d\.]+)/){push @dp, $1;}
	}
	my $p50 = int(@dp * .5);
	my $avgdp50 = 3*$dp[$p50]/$sample_size;
	close IN; 
	
	my $chrcmd="";
	###filter by chromosomes if needed
	if (defined $cfg{step1}{variant_filtering}{chr}){
		my @a = split /,/, $cfg{step1}{variant_filtering}{chr};
		foreach (@a){
			$chrcmd .= " --chr $_";
		}
	}

	#default SNV sites 1: PASS + max-meanDP: 3*$avgdp50 + max-missing: 0.8 + biallelic + selected chr (optional)
	open SH, ">$shpath/variant.filtering.s2.sh";
	print SH "vcftools --gzvcf $outpath/PASS.SNP.vcf.gz $chrcmd --min-meanDP 1 --max-meanDP $avgdp50 --max-missing 0.8 --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --stdout \| bgzip -c > $outpath/PASS.SNP.DP.vcf.gz\n";
	close SH;
	`sh $shpath/variant.filtering.s2.sh 1>$shpath/variant.filtering.s2.sh.o 2>$shpath/variant.filtering.s2.sh.e` unless ($skipsh ==1);
	$cfg{step1}{variant_filtering}{vcf}="$outpath/PASS.SNP.DP.vcf.gz";

	### prepare the setting for low LD prunning
	$cfg{step1}{variant_filtering}{scaffold_number_limit} = 95 unless (defined $cfg{step1}{variant_filtering}{scaffold_number_limit});
	$cfg{step1}{variant_filtering}{scaffold_length_cutoff} = 0 unless (defined $cfg{step1}{variant_filtering}{scaffold_length_cutoff});
	#set ld prunning cut off
	$cfg{step1}{variant_filtering}{ldwindowsize} = 10 unless (defined $cfg{step1}{variant_filtering}{ldwindowsize});
	$cfg{step1}{variant_filtering}{ldwindowstep} = 5 unless (defined $cfg{step1}{variant_filtering}{ldwindowstep});
	$cfg{step1}{variant_filtering}{ldcutoff} = 0.3 unless (defined $cfg{step1}{variant_filtering}{ldcutoff});

	### load reference and map scaffold to number list 
	my $i=1;my $j=0;
	open OT, ">$outpath/chr_map.list";
	foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
		print OT "$id\t$i\n";
		$i++;
		if (($genome->{len}{$id} >= $cfg{step1}{variant_filtering}{scaffold_length_cutoff})&&($j < $cfg{step1}{variant_filtering}{scaffold_number_limit})){
			$j++;
		}
	}
	close OT;
	$cfg{step1}{variant_filtering}{scaffold_number_limit} = $j;
	my $scaffold_number_limit = $cfg{step1}{variant_filtering}{scaffold_number_limit};

	#Based on default SNV sites + no singletons + no missing data + minDP 2 + minQ 30  (high quality SNPs)
	open SH, ">$shpath/variant.filtering.s3.sh";
	print SH "cd $outpath\n";
	print SH "vcftools --gzvcf $cfg{step1}{variant_filtering}{vcf} --singletons --stdout >$outpath/singletons.list\n";
	print SH "vcftools --gzvcf $cfg{step1}{variant_filtering}{vcf} --exclude-positions $outpath/singletons.list --max-missing 1 --max-alleles 2 --minDP 2 --minQ 30 --recode --recode-INFO-all --stdout |bgzip -c >$outpath/high_confidence.vcf.gz\n";

	#Based on high quality SNV sites + low LD
	print SH 'bcftools annotate --threads 10 --rename-chrs', " $outpath/chr_map.list" ," $outpath/high_confidence.vcf.gz", '|perl -ne \'if (/#\S+ID=(\d+)/){if($1<=',"$scaffold_number_limit",'){print;}}elsif(/^#/){print;}elsif(/^(\d+)\s+/){if($1<= ',"$scaffold_number_limit",'){print;}}\'|vcftools --vcf - --plink',"\n";
	
	print SH <<EOF;
plink --file out --make-bed --chr-set $scaffold_number_limit no-xy no-mt no-y
plink --bfile plink --indep-pairwise $cfg{step1}{variant_filtering}{ldwindowsize} $cfg{step1}{variant_filtering}{ldwindowstep} $cfg{step1}{variant_filtering}{ldcutoff} --chr-set $scaffold_number_limit no-xy no-mt no-y
plink --bfile plink --extract plink.prune.in --make-bed --out high_confidence_prunned --chr-set $scaffold_number_limit no-xy no-mt no-y
plink --bfile high_confidence_prunned --chr-set $scaffold_number_limit no-xy no-mt no-y --recode vcf --out high_confidence_pre
cat high_confidence_pre.vcf|perl -ne 'print unless (/CHROM/);if(/CHROM/){s/_\\S+//g;print;}' | bgzip -c >$outpath/high_confidence_prunned.vcf.gz
EOF
	close SH;

	#switch on the bash running
	`sh $shpath/variant.filtering.s3.sh 1>$shpath/variant.filtering.s3.sh.o 2>$shpath/variant.filtering.s3.sh.e` unless ($skipsh ==1);
	
	$cfg{step1}{variant_filtering}{plink_data}="$outpath/high_confidence_prunned";
	$cfg{step1}{variant_filtering}{high_confidence_vcf}="$outpath/high_confidence.vcf.gz";
	$cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}="$outpath/high_confidence_prunned.vcf.gz";
	
	### format a report
	my %report;
	if ($cfg{step1}{variant_filtering}{vcf} =~ /.gz$/) { open(IN, "gunzip -c $cfg{step1}{variant_filtering}{vcf} |") || die "can’t open pipe to $cfg{step1}{variant_filtering}{vcf}";}
		else { open(IN, $cfg{step1}{variant_filtering}{vcf}) || die "can’t open $cfg{step1}{variant_filtering}{vcf}";}

	$report{snv1}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv1}{number} ++;
	}
	close IN;
	$report{snv1}{singletons} = `wc -l $outpath/singletons.list | cut -d" " -f1`;
	$report{snv1}{singletons} = $report{snv1}{singletons} - 1;

	if ($cfg{step1}{variant_filtering}{high_confidence_vcf} =~ /.gz$/) { open(IN, "gunzip -c $cfg{step1}{variant_filtering}{high_confidence_vcf} |") || die "can’t open pipe to $cfg{step1}{high_confidence_vcf}{vcf}";}
		else { open(IN, $cfg{step1}{variant_filtering}{high_confidence_vcf}) || die "can’t open $cfg{step1}{variant_filtering}{high_confidence_vcf}";}
	$report{snv2}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv2}{number} ++;
	}
	close IN;

	if ($cfg{step1}{variant_filtering}{lowld_high_confidence_vcf} =~ /.gz$/) { open(IN, "gunzip -c $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf} |") || die "can’t open pipe to $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}";}
		else { open(IN, $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}) || die "can’t open $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}";}
	$report{snv3}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv3}{number} ++;
	}
	close IN;
	####
	open OT, ">$outpath/snv.summary.txt";
	print OT "SNV set\tSNV size\tSingleton size\n";
	print OT "PASS.SNP.DP.vcf.gz\t$report{snv1}{number}\t$report{snv1}{singletons}\n";
	print OT "high_confidence.vcf.gz\t$report{snv2}{number}\t0\n";
	print OT "high_confidence_prunned.vcf.gz\t$report{snv3}{number}\t0\n";
	close OT; 
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( "$cfg{args}{outdir}/allcfg.yml" );
}

1;