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


	if ($cfg{step1}{variant_calling_mode} eq 'joint_recalibration') {

		my %var;

		$var{reference} = $reference;
		$var{outpath} = $outpath;
		$var{shpath} = $shpath;
		$var{ploidy} = $ploidy;

		if ( !-d "$var{outpath}/JointBQSR/" ) {make_path "$var{outpath}/JointBQSR/" or die "Failed to create path: $var{outpath}/JointBQSR/";}

		open CL, ">$var{shpath}/joint_bqsr_s1.list";
		my $sample_gvcfs = "";
		foreach my $sample (keys %samplelist){
			my $sample_outpath="$var{outpath}/$sample"; 
			if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

			open SH, ">$var{shpath}/joint_bqsr_s1_$sample.sh";
			print SH "#!/bin/sh\ncd $sample_outpath\n";

			print SH "gatk HaplotypeCaller \\\n";
			print SH "	--emit-ref-confidence GVCF \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-ploidy $var{ploidy} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.HC.1st.gvcf.gz && echo \"** GVCF ${sample}.HC.1st.gvcf.gz done\" && \\\n";
		  	print SH "rm -f $sample.sorted.markdup.BQSR2nd.bam && echo \"** variant calling done **\" > $var{shpath}/$sample.variant_calling.finished.txt\n";
		  	
		  	close SH;

		  	print CL "sh $var{shpath}/joint_bqsr_s1_$sample.sh 1>$var{shpath}/joint_bqsr_s1_$sample.sh.o 2>$var{shpath}/joint_bqsr_s1_$sample.sh.e\n";

		  	$sample_gvcfs .= "	-V $var{outpath}/$sample/$sample.HC.1st.gvcf.gz \\\n";

		}
		close CL;

		`parallel -j $cfg{args}{threads} < $var{shpath}/joint_bqsr_s1.list`;
		#`perl $Bin/lib/qsub.pl -r -d $var{shpath}/joint_bqsr_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=5G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/joint_bqsr_s1.list` unless (defined $opts{skipsh});

		#### First round of joint calling START #####
		if ( !-d "$var{outpath}/JointBQSR/" ) {
			make_path "$var{outpath}/JointBQSR/" or die "Failed to create path: $var{outpath}/JointBQSR/";
		}

		open CL, ">$var{shpath}/joint_bqsr_s2.list";
		open SH, ">$var{shpath}/joint_bqsr_s2.sh";

		print SH "gatk CombineGVCFs \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "$sample_gvcfs";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling.HC.g1st.vcf.gz && echo \"** JointCalling.HC.g1st.vcf.gz done ** \"\n";

		print SH "gatk GenotypeGVCFs \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-ploidy $var{ploidy} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC.g1st.vcf.gz \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling.HC1st.vcf.gz && echo \"** JointCalling.HC1st.vcf.gz done ** \"\n";

		print SH "gatk SelectVariants \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC1st.vcf.gz \\\n";
		print SH "	--select-type-to-include SNP \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_snps1st.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_snps1st done\" && \\\n";

		print SH "gatk VariantFiltration \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling_raw_snps1st.vcf \\\n";
		print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
		print SH "	--filter-name \"my_snp_filter\" \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_filtered_snps1st.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_snps1st done\" \n";

		print SH "gatk SelectVariants \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC1st.vcf.gz \\\n";
		print SH "	--select-type-to-include INDEL \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_indels1st.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_snps1st done\" && \\\n";

		print SH "gatk VariantFiltration \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling_raw_indels1st.vcf \\\n";
		print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
		print SH "	--filter-name \"my_indel_filter\" \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_filtered_indels1st.vcf && echo \"** JointBQSR/JointCalling_raw_snps1st done\" \n";

		print SH "bgzip -f $var{outpath}/JointBQSR/JointCalling_filtered_snps1st.vcf\n";
		print SH "tabix -f $var{outpath}/JointBQSR/JointCalling_filtered_snps1st.vcf.gz \n";
		print SH "bgzip -f $var{outpath}/JointBQSR/JointCalling_filtered_indels1st.vcf\n";
		print SH "tabix -f $var{outpath}/JointBQSR/JointCalling_filtered_indels1st.vcf.gz \n";

		close SH;
		print CL "sh $var{shpath}/joint_bqsr_s2.sh 1>$var{shpath}/joint_bqsr_s2.sh.o 2>$var{shpath}/joint_bqsr_s2.sh.e\n";
		close CL;

		`parallel -j $cfg{args}{threads} < $var{shpath}/joint_bqsr_s2.list`;
		#`perl $Bin/lib/qsub.pl -r -d $var{shpath}/joint_bqsr_s2_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=10G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/joint_bqsr_s2.list` unless (defined $opts{skipsh});
		#### First round of joint calling END ####

		#### First round of recalibration START ####
		open CL, ">$var{shpath}/joint_bqsr_s3.list";
		foreach my $sample (keys %samplelist){
			my $sample_outpath="$var{outpath}/$sample";

			open SH, ">$var{shpath}/joint_bqsr_s3_$sample.sh";
			print SH "#!/bin/sh\ncd $sample_outpath\n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.recal_data.table \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_snps1st.vcf.gz \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_indels1st.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

			print SH "gatk ApplyBQSR \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.BQSR.bam \\\n";
			print SH "	--bqsr-recal-file $sample.sorted.markdup.recal_data.table && \\\n";
			print SH "samtools index $sample.sorted.markdup.BQSR.bam && echo \"** $sample.sorted.markdup.BQSR.bam index done **\" \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR.bam  \\\n";
			print SH "	-O $sample.sorted.markdup.recal_data1st_after.table \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_snps1st.vcf.gz \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_indels1st.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

			# HaplotypeCaller
			print SH "gatk HaplotypeCaller \\\n";
			print SH "  --emit-ref-confidence GVCF\\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-ploidy $var{ploidy} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR.bam \\\n";
			print SH "	-O $sample.HC.2nd.gvcf.gz && echo \"** GVCF $sample.HC.2nd.gvcf.gz done\" && \\\n";
			print SH "rm -f $sample.sorted.markdup.BQSR.bam\n";
			
			close SH;
			print CL "sh $var{shpath}/joint_bqsr_s3_$sample.sh 1>$var{shpath}/joint_bqsr_s3_$sample.sh.o 2>$var{shpath}/joint_bqsr_s3_$sample..sh.e\n";

			$sample_gvcfs .= "	-V $var{outpath}/$sample/$sample.HC.2nd.gvcf.gz \\\n";
		}
		close CL;

		`parallel -j $cfg{args}{threads} < $var{shpath}/joint_bqsr_s3.list`;
		#`perl $Bin/lib/qsub.pl -r -d $var{shpath}/joint_bqsr_s3_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=10G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/joint_bqsr_s3.list` unless (defined $opts{skipsh});
		#### First round of recalibration END ####

		#### Second round of joint calling START ####
		open CL, ">$var{shpath}/joint_bqsr_s4.list";
		open SH, ">$var{shpath}/joint_bqsr_s4.sh";

		print SH "gatk CombineGVCFs \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "$sample_gvcfs";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling.HC.g2nd.vcf.gz && echo \"** JointCalling.HC.g2nd.vcf.gz done ** \"\n";

		print SH "gatk GenotypeGVCFs \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-ploidy $var{ploidy} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC.g2nd.vcf.gz \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling.HC2nd.vcf.gz && echo \"** JointCalling.HC2nd.vcf.gz done ** \"\n";

		print SH "gatk SelectVariants \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC2nd.vcf.gz \\\n";
		print SH "	--select-type-to-include SNP \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_snps2nd.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_snps2nd done\" && \\\n";

		print SH "gatk VariantFiltration \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling_raw_snps2nd.vcf \\\n";
		print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
		print SH "	--filter-name \"my_snp_filter\" \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_filtered_snps2nd.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_snps2nd done\" \n";

		print SH "gatk SelectVariants \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC.g2nd.vcf.gz \\\n";
		print SH "	--select-type-to-include INDEL \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_indels2nd.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_indels2nd done\" && \\\n";

		print SH "gatk VariantFiltration \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling_raw_indels2nd.vcf \\\n";
		print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
		print SH "	--filter-name \"my_indel_filter\" \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_filtered_indels2nd.vcf && echo \"** JointBQSR/JointCalling_filtered_indels2nd.vcf done\" \n";

		print SH "bgzip -f $var{outpath}/JointBQSR/JointCalling_filtered_snps2nd.vcf\n";
		print SH "tabix -f $var{outpath}/JointBQSR/JointCalling_filtered_snps2nd.vcf.gz \n";
		print SH "bgzip -f $var{outpath}/JointBQSR/JointCalling_filtered_indels2nd.vcf\n";
		print SH "tabix -f $var{outpath}/JointBQSR/JointCalling_filtered_indels2nd.vcf.gz \n";

		close SH;
		print CL "sh $var{shpath}/joint_bqsr_s4.sh 1>$var{shpath}/joint_bqsr_s4.sh.o 2>$var{shpath}/joint_bqsr_s4.sh.e\n";
		close CL;

		`parallel -j $cfg{args}{threads} < $var{shpath}/joint_bqsr_s4.list`;
		#`perl $Bin/lib/qsub.pl -r -d $var{shpath}/joint_bqsr_s4_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=10G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/joint_bqsr_s4.list` unless (defined $opts{skipsh});
		#### Second round of joint calling END ####

		#### Second round of recalibration START ####
		open CL, ">$var{shpath}/joint_bqsr_s5.list";
		foreach my $sample (keys %samplelist){
			my $sample_outpath="$var{outpath}/$sample";

			open SH, ">$var{shpath}/joint_bqsr_s5_$sample.sh";
			print SH "#!/bin/sh\ncd $sample_outpath\n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.recal_data2nd.table \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_snps2nd.vcf.gz\\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_indels2nd.vcf.gz \n";

			print SH "gatk ApplyBQSR \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.BQSR2nd.bam \\\n";
			print SH "	--bqsr-recal-file $sample.sorted.markdup.recal_data2nd.table && \\\n";
			print SH "samtools index $sample.sorted.markdup.BQSR2nd.bam && echo \"** $sample.sorted.markdup.BQSR2nd.bam index done **\" \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR2nd.bam \\\n";
			print SH "	-O $sample.sorted.markdup.recal_data2nd_after.table \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_snps2nd.vcf.gz \\\n";
			print SH "	--known-sites $var{outpath}/JointBQSR/JointCalling_filtered_indels2nd.vcf.gz  \n";

			# HaplotypeCaller
			print SH "gatk HaplotypeCaller \\\n";
			print SH "  --emit-ref-confidence GVCF\\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-ploidy $var{ploidy} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR2nd.bam \\\n";
			print SH "	-O $sample.HC.3rd.gvcf.gz && echo \"** GVCF $sample.HC.3rd.gvcf.gz done\" && \\\n";
			print SH "rm -f $sample.sorted.markdup.BQSR2nd.bam\n";
			
			print SH "gatk AnalyzeCovariates \\\n";
			print SH "	--before-report-file $sample.sorted.markdup.recal_data.table \\\n";
			print SH "	--after-report-file $sample.sorted.markdup.recal_data1st_after.table \\\n";
			print SH "	--plots-report-file $sample.recalQC.1st.pdf\n";

			print SH "gatk AnalyzeCovariates \\\n";
			print SH "	--before-report-file $sample.sorted.markdup.recal_data2nd.table \\\n";
			print SH "	--after-report-file $sample.sorted.markdup.recal_data2nd_after.table \\\n";
			print SH "	--plots-report-file $sample.recalQC.2nd.pdf && echo \"** variant calling done **\" > $var{shpath}/$sample.variant_calling.finished.txt\n";

			close SH;

			print CL "sh $var{shpath}/joint_bqsr_s5_$sample.sh 1>$var{shpath}/joint_bqsr_s5_$sample.sh.o 2>$var{shpath}/joint_bqsr_s5_$sample.sh.e\n";
			$sample_gvcfs .= "	-V $var{outpath}/$sample/$sample.HC.3rd.gvcf.gz \\\n";
		}
		close CL;

		`parallel -j $cfg{args}{threads} < $var{shpath}/joint_bqsr_s5.list`;
		#`perl $Bin/lib/qsub.pl -r -d $var{shpath}/joint_bqsr_s5_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=5G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/joint_bqsr_s5.list` unless (defined $opts{skipsh});
		#### Second round of recalibration END ####

		#### Final round of joint calling START ####
		open CL, ">$var{shpath}/joint_bqsr_s6.list";
		open SH, ">$var{shpath}/joint_bqsr_s6.sh";

		print SH "gatk CombineGVCFs \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "$sample_gvcfs";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling.HC.g3rd.vcf.gz && echo \"** JointCalling.HC.g3rd.vcf.gz done ** \"\n";

		print SH "gatk GenotypeGVCFs \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-ploidy $var{ploidy} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC.g3rd.vcf.gz \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling.HC3rd.vcf.gz && echo \"** JointCalling.HC3rd.vcf.gz done ** \"\n";

		print SH "gatk SelectVariants \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC3rd.vcf.gz \\\n";
		print SH "	--select-type-to-include SNP \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_snps3rd.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_snps3rd done\" && \\\n";

		print SH "gatk VariantFiltration \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling_raw_snps3rd.vcf \\\n";
		print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
		print SH "	--filter-name \"my_snp_filter\" \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_filtered_snps3rd.vcf && echo \"** GVCF JointBQSR/JointCalling_filtered_snps3rd done\" \n";

		print SH "gatk SelectVariants \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling.HC3rd.vcf.gz \\\n";
		print SH "	--select-type-to-include INDEL \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_indels3rd.vcf && echo \"** GVCF JointBQSR/JointCalling_raw_indels3rd done\" && \\\n";

		print SH "gatk VariantFiltration \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-V $var{outpath}/JointBQSR/JointCalling_raw_indels3rd.vcf \\\n";
		print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
		print SH "	--filter-name \"my_indel_filter\" \\\n";
		print SH "	-O $var{outpath}/JointBQSR/JointCalling_raw_indels3rd.vcf && echo \"** JointBQSR/JointCalling_raw_indels3rd done\" \n";

		print SH "bgzip -f $var{outpath}/JointBQSR/JointCalling_filtered_snps3rd.vcf\n";
		print SH "tabix -f $var{outpath}/JointBQSR/JointCalling_filtered_snps3rd.vcf.gz \n";
		print SH "bgzip -f $var{outpath}/JointBQSR/JointCalling_filtered_indels3rd.vcf\n";
		print SH "tabix -f $var{outpath}/JointBQSR/JointCalling_filtered_indels3rd.vcf.gz \n";

		close SH;
		print CL "sh $var{shpath}/joint_bqsr_s6.sh 1>$var{shpath}/joint_bqsr_s6.sh.o 2>$var{shpath}/joint_bqsr_s6.sh.e\n";
		close CL;

		`parallel -j $cfg{args}{threads} < $var{shpath}/joint_bqsr_s6.list`;
		#`perl $Bin/lib/qsub.pl -r -d $var{shpath}/joint_bqsr_s6_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=5G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/joint_bqsr_s6.list` unless (defined $opts{skipsh});
		#### Final round of joint calling END ####
	}else{

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
	}
	close CL;

	`parallel -j $cfg{args}{threads} < $shpath/cmd_step1d.list`;
}

1;