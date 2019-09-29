package PopGenome_Calibration;
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
#    Step 1d Variant calling    #
#			   #
#################################
sub CALIBRATION{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};

	my $ploidy = $cfg{args}{ploidy};

	my %samplelist = %{$cfg{fqdata}};

	my $reference = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	my $outpath = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}"; 
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$cfg{ref}{choose}";
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	open CL, ">$shpath/cmd_step1d.list";

	foreach my $sample (keys %samplelist){
		my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}
		open SH, ">$shpath/$sample.step1d.sh";

		print SH "#!/bin/sh\ncd $sample_outpath\n";		
		# MarkDuplicates
		print SH "gatk MarkDuplicates \\\n";
	  	print SH "	--INPUT $sample.sorted.bam \\\n";
	  	print SH "	--OUTPUT $sample.sorted.markdup.bam \\\n";
	  	print SH "	--METRICS_FILE $sample.sorted.markdup_metrics.txt && \\\n";
	  	print SH "rm -f $sample.sorted.bam && \\\n";
	  	print SH "echo \"** $sample.sorted.markdup.bam done **\" && \\\n";
	  	print SH "samtools index $sample.sorted.markdup.bam && echo \"** $sample.sorted.markdup.bam index done **\" \n";
	  	unless ($cfg{step1}{variant_calling_mode} eq 'fast'){
	  	# HaplotypeCaller
		  	print SH "gatk HaplotypeCaller \\\n";
		  	print SH "	-R $reference \\\n";
		  	print SH "	-ploidy $ploidy \\\n";
		 	print SH "	-I $sample.sorted.markdup.bam \\\n";
		  	print SH "	-O $sample.HC.g1st.vcf.gz && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";
	  	
		  	###########SNP extraction and filtering#######
		  	# SelectVariants
		  	print SH "gatk SelectVariants \\\n";
		  	print SH "	-R $reference \\\n";
		  	print SH "	-V $sample.HC.g1st.vcf.gz \\\n";
		  	print SH "	--select-type-to-include SNP \\\n";
			print SH "	-O $sample\_raw_snps1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";
			# VariantFiltration
			print SH "gatk VariantFiltration \\\n";
		  	print SH "	-R $reference \\\n";
		  	print SH "	-V $sample\_raw_snps1st.vcf \\\n";
		  	print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
			print SH "	--filter-name \"my_snp_filter\" \\\n";
			print SH "	-O $sample\_filtered_snps1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			###########INDEL extraction and filtering#######
			print SH "gatk SelectVariants \\\n";
		  	print SH "	-R $reference \\\n";
		  	print SH "	-V $sample.HC.g1st.vcf.gz \\\n";
		  	print SH "	--select-type-to-include INDEL \\\n";
			print SH "	-O $sample\_raw_indels1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			print SH "gatk VariantFiltration \\\n";
		  	print SH "	-R $reference \\\n";
		  	print SH "	-V $sample\_raw_indels1st.vcf \\\n";
		  	print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
			print SH "	--filter-name \"my_indel_filter\" \\\n";
			print SH "	-O $sample\_filtered_indels1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			print SH "bgzip -f $sample\_filtered_snps1st.vcf\n";
			print SH "tabix -f $sample\_filtered_snps1st.vcf.gz \n";
			print SH "bgzip -f $sample\_filtered_indels1st.vcf\n";
			print SH "tabix -f $sample\_filtered_indels1st.vcf.gz \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $reference \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
		  	print SH "	-O $sample.sorted.markdup.recal_data.table \\\n";
			print SH "	--known-sites $sample\_filtered_snps1st.vcf.gz \\\n";
			print SH "	--known-sites $sample\_filtered_indels1st.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

			print SH "gatk ApplyBQSR \\\n";
			print SH "	-R $reference \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.BQSR.bam \\\n";
			print SH "	--bqsr-recal-file $sample.sorted.markdup.recal_data.table && \\\n";
			print SH "samtools index $sample.sorted.markdup.BQSR.bam && echo \"** $sample.sorted.markdup.BQSR.bam index done **\" \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $reference \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR.bam  \\\n";
		  	print SH "	-O $sample.sorted.markdup.recal_data1st_after.table \\\n";
			print SH "	--known-sites $sample\_filtered_snps1st.vcf.gz \\\n";
			print SH "	--known-sites $sample\_filtered_indels1st.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";
		}
		close SH;
		print CL "sh $shpath/$sample.step1d.sh 1>$shpath/$sample.step1d.sh.o 2>$shpath/$sample.step1d.sh.e\n";
	}
	close CL;

	`parallel -j $cfg{args}{threads} < $shpath/cmd_step1d.list`;
}

1;