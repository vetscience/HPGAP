package PopGenome;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use YAML::Tiny;
use Bio::SeqIO;

#-------------------------------------- Prepare reference index --------------------------------------

############################
#			   #
#    Step 0 Indexing       #
#			   #
############################
sub INDEXING{

	#######
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	
	my $outpath = "$cfg{args}{outdir}/00.INDEXING/"; 
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/00.INDEXING";
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}
	
	open SH, ">$shpath/index.sh";
	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		my $reference = $cfg{ref}{db}{$temp_ref}{path};
		my $genome=LOADREF($reference);

		my $ref_base=basename($reference);
		my $ref_dir=dirname($reference);
		my $ref_name=$ref_base;
		$ref_name =~ s/\.fasta$//;
		$ref_name =~ s/\.fa$//;

		if (! -e "$outpath/$ref_base.bwt"){
			if ( !-d "$outpath" ) { make_path $shpath or die "Failed to create path: $outpath";}
			print SH "#!/bin/sh\ncd $outpath\n";
			#if ( -e $mtgenome){
			#	`ln -s $mtgenome $outpath/`;
			#	print ID "cat $reference $mtgenome >$ref_base\n";
			#	print ID "bwa index $ref_base $ref_name \n";
			#	print ID "picard CreateSequenceDictionary R=$ref_base O=$ref_name.dict \n";
			#	print ID "samtools faidx $ref_base \n";
			#}else{
			`cp -f $reference $outpath/` if ( !-e "$outpath/$ref_base" );
			print SH "bwa index $ref_base $ref_name \n";
			print SH "picard CreateSequenceDictionary R=$ref_base O=$ref_name.dict \n";
			print SH "samtools faidx $ref_base \n";
		#}
		}
		
		$cfg{ref}{db}{$temp_ref}{path} = "$outpath/$ref_base";
	}
	close SH;
	`sh $shpath/index.sh` unless ($skipsh ==1);

	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $yml_file );
}
#-------------------------------------- Analysis -----------------------------------------------

############################
#			   #
#   Step 1a Data filtering #
#			   #
############################
sub DATA_FILTERING{

	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	
	my $outpath = "$cfg{args}{outdir}/01.QualityControl/read_filtering/"; 
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_filtering/";
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	my %samplelist = %{$cfg{fqdata}};

	open CL, ">$shpath/cmd_step1a.list";
	foreach my $sample (keys %samplelist){
		my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}
		open SH, ">$shpath/$sample.step1a.sh";
		print SH "#!/bin/sh\ncd $sample_outpath\n";
		foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
			my $read;
			if ($samplelist{$sample}{rawdata}{$lib}{fq1} =~ /gz$/){
				$read = `gunzip -c $samplelist{$sample}{rawdata}{$lib}{fq1}|head -n 2|tail -n 1`;
			}else{
				$read = `cat $samplelist{$sample}{rawdata}{$lib}{fq1}|head -n 2|tail -n 1`;
			}
			$samplelist{$sample}{rawdata}{$lib}{Length} = split //, $read;
			$samplelist{$sample}{rawdata}{$lib}{Length}=int($samplelist{$sample}{rawdata}{$lib}{Length}*0.7);
			if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
				print SH "reformat.sh overwrite=true in1=$samplelist{$sample}{rawdata}{$lib}{fq1} in2=$samplelist{$sample}{rawdata}{$lib}{fq2} bhist=$lib\_bhist.txt qhist=$lib\_qhist.txt aqhist=$lib\_aqhist.txt lhist=$lib\_lhist.txt  gchist=$lib\_gchist.txt && \\\n";
				print SH "trimmomatic PE -threads $cfg{args}{threads} -phred$samplelist{$sample}{rawdata}{$lib}{Phred} $samplelist{$sample}{rawdata}{$lib}{fq1} $samplelist{$sample}{rawdata}{$lib}{fq2} $lib\_1.filt.fq.gz $lib\_1.filt.unpaired.fq.gz $lib\_2.filt.fq.gz $lib\_2.filt.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:$samplelist{$sample}{rawdata}{$lib}{Length} && \\\n";
				print SH "reformat.sh overwrite=true in1=$lib\_1.filt.fq.gz in2=$lib\_2.filt.fq.gz bhist=$lib\_bhist.filt.txt qhist=$lib\_qhist.filt.txt aqhist=$lib\_aqhist.filt.txt lhist=$lib\_lhist.filt.txt gchist=$lib\_gchist.filt.txt \n";
			}if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
				print SH "reformat.sh overwrite=true in1=$samplelist{$sample}{rawdata}{$lib}{fq1} bhist=$lib\_bhist.txt qhist=$lib\_qhist.txt aqhist=$lib\_aqhist.txt lhist=$lib\_lhist.txt  gchist=$lib\_gchist.txt && \\\n";
				print SH "trimmomatic SE -threads $cfg{args}{threads} -phred$samplelist{$sample}{rawdata}{$lib}{Phred} $samplelist{$sample}{rawdata}{$lib}{fq1} $lib\_1.filt.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:$samplelist{$sample}{rawdata}{$lib}{Length} && \\\n";
				print SH "reformat.sh overwrite=true in1=$lib\_1.filt.fq.gz bhist=$lib\_bhist.filt.txt qhist=$lib\_qhist.filt.txt aqhist=$lib\_aqhist.filt.txt lhist=$lib\_lhist.filt.txt gchist=$lib\_gchist.filt.txt \n";			
			}
		}
		close SH;
		print CL "sh $shpath/$sample.step1a.sh 1>$shpath/$sample.step1a.sh.o 2>$shpath/$sample.step1a.sh.e\n";
	}
	close CL;
	`parallel -j $cfg{args}{threads} < $shpath/cmd_step1a.list` unless ($skipsh ==1);
}

############################
#			   #
#    Step 1b Mapping       #
#			   #
############################
sub MAPPING{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};

	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		my $outpath = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
		if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
		my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
		if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

		my $reference = $cfg{ref}{db}{$temp_ref}{path};	
		open CL, ">$shpath/cmd_step1b.list";
		foreach my $sample (keys %samplelist){
			my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

			open SH, ">$shpath/$sample.step1b.sh";		
			print SH "#!/bin/sh\ncd $sample_outpath\n";
			foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
				print SH "bwa mem $reference ../../read_filtering/$sample/$lib\_1.filt.fq.gz ../../read_filtering/$sample/$lib\_2.filt.fq.gz -t $cfg{args}{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ $cfg{args}{threads} -F 4 - -o $lib\_filt.bam && \\\n" if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE");
				print SH "bwa mem $reference ../../read_filtering/$sample/$lib\_1.filt.fq.gz -t $cfg{args}{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ 10 -F 4 - -o $lib\_filt.bam && \\\n" if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE");
				print SH "samtools sort -@ $cfg{args}{threads} $lib\_filt.bam -o $lib\_filt.sort.bam --output-fmt BAM && \\\n";
				print SH "rm -f $lib\_filt.bam\n";
			}

			if (keys %{$samplelist{$sample}{rawdata}} == 1){
				foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){print SH "mv $lib\_filt.sort.bam $sample.sorted.bam\n";}
				print SH "samtools stats -@ $cfg{args}{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** bam.stats.txt done **\"\n";
			}

			if (keys %{$samplelist{$sample}{rawdata}} > 1){
				print SH "samtools merge -nr -@ $cfg{args}{threads} $sample.bam *_filt.sort.bam && echo \"** $sample.sort.bam done **\" && rm -f *_filt.sort.bam\n";
				print SH "samtools sort -@ $cfg{args}{threads} $sample.bam -o $sample.sorted.bam --output-fmt BAM && echo \"** $sample.sorted.bam done **\" && rm -f $sample.bam\n";
				print SH "samtools stats -@ $cfg{args}{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** bam.stats.txt done **\"\n";
			}
			close SH;
			print CL "sh $shpath/$sample.step1b.sh 1>$shpath/$sample.step1b.sh.o 2>$shpath/$sample.step1b.sh.e \n";
		}
		close CL;
		`parallel -j $cfg{args}{threads} < $shpath/cmd_step1b.list` unless ($skipsh ==1);
	}

	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $yml_file );
}
#################################
#			   #
#    Step 1c Report summary	    #
#			   #
#################################
sub READ_REPORT{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};

	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		my $outpath = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
		if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
		my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
		if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

		my $report_outpath="$outpath/Report"; if ( !-d $report_outpath ) {make_path $report_outpath or die "Failed to create path: $report_outpath";}

		my $report_sample_outpath="$outpath/Report/Samples"; if ( !-d $report_sample_outpath ) {make_path $report_sample_outpath or die "Failed to create path: $report_outpath";}

		foreach my $sample (keys %samplelist){
			
			my $sample_report_outpath="$outpath/Report/Samples/$sample"; if ( !-d $sample_report_outpath ) {make_path $sample_report_outpath or die "Failed to create path: $sample_report_outpath";}

			`cp $outpath/$sample/bam.stats.txt $sample_report_outpath`;
			`cp $cfg{args}{outdir}/01.QualityControl/read_filtering/$sample/*hist.txt $sample_report_outpath`;
			`cp $cfg{args}{outdir}/01.QualityControl/read_filtering/$sample/*hist.filt.txt $sample_report_outpath`;
			`cat $outpath/$sample/bam.stats.txt|grep COV|grep -v "#" >$sample_report_outpath/COV.stat.txt`;
			`cat $outpath/$sample/bam.stats.txt|grep IS|grep -v "#" >$sample_report_outpath/IS.stat.txt`;
			`cat $outpath/$sample/bam.stats.txt|grep ^SN | cut -f 2-|sed "s/ //g;" >$sample_report_outpath/SN.stat.txt`;
		}		

		open SH, ">$shpath/read_report.sh";
		print SH "#!/bin/sh\ncd $report_outpath\n";
		print SH "Rscript --vanilla $Bin/lib/ReadSummary.R $report_outpath $report_sample_outpath $temp_ref.Summary.xls\n";
		`sh $shpath/read_report.sh 1>$shpath/read_report.sh.o 2>$shpath/read_report.sh.e`;
	}
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $yml_file );
}

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

		close SH;
		print CL "sh $shpath/$sample.step1d.sh 1>$shpath/$sample.step1d.sh.o 2>$shpath/$sample.step1d.sh.e\n";
	}
	close CL;

	`parallel -j $cfg{args}{threads} < $shpath/cmd_step1d.list`;
}
#################################
#			   #
#    Step 1e Variant calling    #
#			   #
#################################
sub VARIANT_CALLING{
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

	open CL, ">$shpath/cmd_step1e.list";

	foreach my $sample (keys %samplelist){
		my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}
		open SH, ">$shpath/$sample.step1e.sh";

		print SH "#!/bin/sh\ncd $sample_outpath\n";

	  	# HaplotypeCaller
	  	print SH "gatk HaplotypeCaller \\\n";
	  	print SH "	-R $reference \\\n";
	  	print SH "	-ploidy $ploidy \\\n";
	  	print SH "	-I $sample.sorted.markdup.BQSR.bam \\\n";
	  	print SH "	-O $sample.HC.g2nd.vcf.gz && echo \"** GVCF $sample.HC.g.vcf.gz done\" && \\\n";
	  	print SH "rm -f $sample.sorted.markdup.BQSR.bam\n";
	  	
	  	###########SNP extraction and filtering#######
	  	# SelectVariants
	  	print SH "gatk SelectVariants \\\n";
	  	print SH "	-R $reference \\\n";
	  	print SH "	-V $sample.HC.g2nd.vcf.gz \\\n";
	  	print SH "	--select-type-to-include SNP \\\n";
		print SH "	-O $sample\_raw_snps2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";
		# VariantFiltration
		print SH "gatk VariantFiltration \\\n";
	  	print SH "	-R $reference \\\n";
	  	print SH "	-V $sample\_raw_snps2nd.vcf \\\n";
	  	print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
		print SH "	--filter-name \"my_snp_filter\" \\\n";
		print SH "	-O $sample\_filtered_snps2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

		###########INDEL extraction and filtering#######
		print SH "gatk SelectVariants \\\n";
	  	print SH "	-R $reference \\\n";
	  	print SH "	-V $sample.HC.g2nd.vcf.gz \\\n";
	  	print SH "	--select-type-to-include INDEL \\\n";
		print SH "	-O $sample\_raw_indels2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

		print SH "gatk VariantFiltration \\\n";
	  	print SH "	-R $reference \\\n";
	  	print SH "	-V $sample\_raw_indels2nd.vcf \\\n";
	  	print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
		print SH "	--filter-name \"my_indel_filter\" \\\n";
		print SH "	-O $sample\_filtered_indels2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

		print SH "bgzip -f $sample\_filtered_snps2nd.vcf\n";
		print SH "tabix -f $sample\_filtered_snps2nd.vcf.gz \n";
		print SH "bgzip -f $sample\_filtered_indels2nd.vcf\n";
		print SH "tabix -f $sample\_filtered_indels2nd.vcf.gz \n";

		print SH "gatk BaseRecalibrator \\\n";
		print SH "	-R $reference \\\n";
		print SH "	-I $sample.sorted.markdup.bam \\\n";
	  	print SH "	-O $sample.sorted.markdup.recal_data2nd.table \\\n";
		print SH "	--known-sites $sample\_filtered_snps2nd.vcf.gz \\\n";
		print SH "	--known-sites $sample\_filtered_indels2nd.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

		print SH "gatk ApplyBQSR \\\n";
		print SH "	-R $reference \\\n";
		print SH "	-I $sample.sorted.markdup.bam \\\n";
		print SH "	-O $sample.sorted.markdup.BQSR2nd.bam \\\n";
		print SH "	--bqsr-recal-file $sample.sorted.markdup.recal_data2nd.table && \\\n";
		print SH "samtools index $sample.sorted.markdup.BQSR2nd.bam && echo \"** $sample.sorted.markdup.BQSR2nd.bam index done **\" \n";

		print SH "gatk BaseRecalibrator \\\n";
		print SH "	-R $reference \\\n";
		print SH "	-I $sample.sorted.markdup.BQSR2nd.bam \\\n";
	  	print SH "	-O $sample.sorted.markdup.recal_data2nd_after.table \\\n";
		print SH "	--known-sites $sample\_filtered_snps2nd.vcf.gz \\\n";
		print SH "	--known-sites $sample\_filtered_indels2nd.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

		print SH "gatk HaplotypeCaller \\\n";
		print SH "	--emit-ref-confidence GVCF \\\n";
		print SH "	-R $reference \\\n";
		print SH "	-ploidy $ploidy \\\n";
		print SH "	-I $sample.sorted.markdup.BQSR2nd.bam \\\n";
		print SH "	-O $sample.HC.gvcf.gz && echo \"** GVCF ${sample}.HC.g.vcf.gz done\" && \\\n";
	  	print SH "rm -f $sample.sorted.markdup.BQSR2nd.bam\n";

	  	print SH "gatk AnalyzeCovariates \\\n";
	  	print SH "	--before-report-file $sample.sorted.markdup.recal_data.table \\\n";
	  	print SH "	--after-report-file $sample.sorted.markdup.recal_data1st_after.table \\\n";
	  	print SH "	--plots-report-file $sample.recalQC.1st.pdf\n";

	  	print SH "gatk AnalyzeCovariates \\\n";
	  	print SH "	--before-report-file $sample.sorted.markdup.recal_data2nd.table \\\n";
	  	print SH "	--after-report-file $sample.sorted.markdup.recal_data2nd_after.table \\\n";
	  	print SH "	--plots-report-file $sample.recalQC.2nd.pdf\n";

		close SH;
		print CL "sh $shpath/$sample.step1e.sh 1>$shpath/$sample.step1e.sh.o 2>$shpath/$sample.step1e.sh.e\n";
	}
	close CL;
	`parallel -j $cfg{args}{threads} < $shpath/cmd_step1e.list` unless ($skipsh ==1);
}
#################################
#			   #
#    Step 1f Variant calling    #
#			   #
#################################
sub COMBINE_CALLING{

	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};

	my $ploidy = $cfg{args}{ploidy};

	my %samplelist = %{$cfg{fqdata}};

	my $reference = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	my $outpath = "$cfg{args}{outdir}/01.QualityControl/"; 
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/Combined";
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/final.calling.sh";
	## Joint genotyping
	## First, merge all the gvcf results, then perform GenotypeGVCFs
	my $sample_gvcfs = "";
	foreach my $sample (keys %samplelist){
		$sample_gvcfs .= "	-V $outpath/read_mapping.$cfg{ref}{choose}/$sample/$sample.HC.gvcf.gz \\\n";
	}
	print SH "#!/bin/sh\ncd $outpath\n";
	# Consider replace CombineGCVFs with GenomicsDBimport
	# Through CombineGVCFs or GenomicsDBimport, we gather all the per-sample GVCFs and pass them to GenotypeGVCFs, which produces a set of joint-called SNP and indel calls ready for fitlering.

	print SH "gatk CombineGVCFs \\\n";
	print SH "	-R $reference \\\n";
	print SH "$sample_gvcfs";
	print SH "	-O $outpath/Combined/Combined.HC.g.vcf.gz && echo \"** Combined.HC.g.vcf.gz done ** \"\n";

	print SH "gatk GenotypeGVCFs \\\n";
	print SH "	-R $reference \\\n";
	print SH "	-ploidy $ploidy \\\n";
	print SH "	-V $outpath/Combined/Combined.HC.g.vcf.gz \\\n";
	print SH "	-O $outpath/Combined/Combined.HC.vcf.gz && echo \"** Combined.HC.vcf.gz done ** \"\n";

	# The established way to filter the raw variant callset is to use variant quality score recalibration (VQSR), which uses machine learning to identify annotation profiles of variants that are likely to be real, and assigns a VQSLOD score to each variant that is much more reliable than the QUAL score calculated by the caller. In the first step of this two-step process, the program builds a model based on training variants, then applies that model to the data to assign a well-calibrated probability to each variant call.
	# We can then use this variant quality score in the second step to filter the raw call set, thus producing a subset of calls with our desired level of quality, fine-tuned to balance specificity and sensitivity.
	## SNP mode
	close SH;
	`sh $shpath/final.calling.sh 1>$shpath/final.calling.sh.o 2>$shpath/final.calling.sh.e` unless ($skipsh ==1);
}
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
	if ($cfg{step1}{variant_filtering}{vcf} =~ /.gz$/) { open(IN, "gunzip -c $outpath/$cfg{step1}{variant_filtering}{vcf} |") || die "can’t open pipe to $cfg{step1}{variant_filtering}{vcf}";}
		else { open(IN, $cfg{step1}{variant_filtering}{vcf}) || die "can’t open $cfg{step1}{variant_filtering}{vcf}";}

	$report{snv1}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv1}{number} ++;
	}
	close IN;
	$report{snv1}{singletons} = `wc -l singletons.list`;
	$report{snv1}{singletons} = $report{snv1}{singletons} - 1;

	if ($cfg{step1}{variant_filtering}{high_confidence_vcf} =~ /.gz$/) { open(IN, "gunzip -c $outpath/$cfg{step1}{variant_filtering}{high_confidence_vcf} |") || die "can’t open pipe to $outpath/$cfg{step1}{high_confidence_vcf}{vcf}";}
		else { open(IN, $cfg{step1}{variant_filtering}{high_confidence_vcf}) || die "can’t open $outpath/$cfg{step1}{variant_filtering}{high_confidence_vcf}";}
	$report{snv2}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv2}{number} ++;
	}
	close IN;

	if ($cfg{step1}{variant_filtering}{lowld_high_confidence_vcf} =~ /.gz$/) { open(IN, "gunzip -c $outpath/$cfg{step1}{variant_filtering}{lowld_high_confidence_vcf} |") || die "can’t open pipe to $outpath/$cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}";}
		else { open(IN, $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}) || die "can’t open $outpath/$cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}";}
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
#################################
#			   #
#    Step 1g VARIANT_COMPARISON #
#			   #
#################################
sub VARIANT_COMPARISON{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};

	my $outpath = "$cfg{args}{outdir}/01.QualityControl/Comparison";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/Comparison";
	my $reference = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	my $gzvcf1 = $cfg{step1}{variant_filtering}{vcf};
	my $gzvcf2 = "$cfg{args}{outdir}/Input/c.briggsae_snps.vcf.gz";
	my $dustbed = "/home/darcy/PopGen_WorkFlow/Cb_Validation/Input/dust.bed";
	my $cdsbed = "/home/darcy/PopGen_WorkFlow/Cb_Validation/Input/cds.bed";

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	my $vcf1_base = basename($gzvcf1); $vcf1_base =~ s/.vcf.gz//g;
	my $vcf2_base = basename($gzvcf2); $vcf2_base =~ s/.vcf.gz//g;

	open SH, ">$shpath/comparison.sh";
	print SH "cd $outpath\n";
	print SH "#zcat $gzvcf1  | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | \"LC_ALL=C sort -k1,1 -k2,2n\"}'  |bgzip -c >$outpath/$vcf1_base.sorted.vcf.gz\n";
	print SH "#zcat $gzvcf2  | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | \"LC_ALL=C sort -k1,1 -k2,2n\"}'  |bgzip -c >$outpath/$vcf2_base.sorted.vcf.gz\n";

	print SH "#vcftools --gzvcf $outpath/$vcf1_base.sorted.vcf.gz --gzdiff $outpath/$vcf2_base.sorted.vcf.gz --diff-site --stdout >$outpath/compare.diff.sites_in_files\n";

	print SH "#cat $outpath/compare.diff.sites_in_files|grep -v CHROM|awk '{if(\$4 == 1){print;}}' |cut -f 1,2 >$outpath/unique.$vcf1_base.position.list\n";
	print SH "#cat $outpath/compare.diff.sites_in_files|grep -v CHROM|awk '{if(\$4 == 2){print;}}' |cut -f 1,3 >$outpath/unique.$vcf2_base.position.list\n";
	print SH "cat $outpath/compare.diff.sites_in_files|awk '{if(\$4 == \"B\"){print;}}'|grep -v CHROM|cut -f 1,3 >$outpath/common.position.list\n";

	open CL, ">$shpath/cmd1.list";
	print CL "vcftools --gzvcf $outpath/$vcf1_base.sorted.vcf.gz --positions $outpath/common.position.list --recode --recode-INFO-all --stdout \| bgzip -c >$outpath/$vcf1_base.sorted.common.vcf.gz\n";
	print CL "vcftools --gzvcf $outpath/$vcf2_base.sorted.vcf.gz --positions $outpath/common.position.list --recode --recode-INFO-all --stdout \| bgzip -c >$outpath/$vcf2_base.sorted.common.vcf.gz\n";
	print CL "vcftools --gzvcf $outpath/$vcf1_base.sorted.vcf.gz --positions $outpath/unique.$vcf1_base.position.list --recode --recode-INFO-all --stdout \| bgzip -c >$outpath/$vcf1_base.sorted.unique.vcf.gz\n";
	print CL "vcftools --gzvcf $outpath/$vcf2_base.sorted.vcf.gz --positions $outpath/unique.$vcf2_base.position.list --recode --recode-INFO-all --stdout \| bgzip -c >$outpath/$vcf2_base.sorted.unique.vcf.gz\n";
	close CL;
	print SH "parallel -j $cfg{args}{threads} < $shpath/cmd1.list\n";

	print SH "StatVCF.pl $outpath/$vcf1_base.sorted.common.vcf.gz >$outpath/$vcf1_base.common.stat\n";
	print SH "StatVCF.pl $outpath/$vcf1_base.sorted.unique.vcf.gz >$outpath/$vcf1_base.unique.stat\n";

	open CL, ">$shpath/cmd2.list";
	my @pre = ("$vcf1_base.sorted.common", "$vcf1_base.sorted.unique","$vcf2_base.sorted.common","$vcf2_base.sorted.unique");
	foreach my $pre (@pre){
		print CL "vcftools --gzvcf $outpath/$pre.vcf.gz --TsTv-summary --stdout >$outpath/$pre.tstv.summary\n";
		print CL "vcftools --gzvcf $outpath/$pre.vcf.gz --bed $dustbed --recode --recode-INFO-all --stdout |bgzip -c >$outpath/$pre.dust.vcf.gz\n";
		print CL "vcftools --gzvcf $outpath/$pre.vcf.gz --bed $cdsbed --recode --recode-INFO-all --stdout |bgzip -c >$outpath/$pre.cds.vcf.gz\n";
	}
	close CL;
	print SH "parallel -j $cfg{args}{threads} < $shpath/cmd2.list\n";
	close SH;

	`sh $shpath/comparison.sh 1>$shpath/comparison.sh.o 2>$shpath/comparison.sh.e` unless ($skipsh ==1);

	my %comp;
	foreach my $pre (@pre){
		open IN, "$outpath/$pre.tstv.summary";
		while (<IN>){
			if (/Ts\s+(\d+)/){$comp{$pre}{Ts}=$1;}
			if (/Tv\s+(\d+)/){$comp{$pre}{Tv}=$1;}
		}
		$comp{$pre}{TsTv}=$comp{$pre}{Ts}/$comp{$pre}{Tv};
		close IN;

		open (IN, "zcat $outpath/$pre.dust.vcf.gz|") or die $!;
		my $i=0;
		while (<IN>){
			next if (/#/);
			$i++;
		}
		$comp{$pre}{dust}=$i;
		close IN;

		open (IN, "zcat $outpath/$pre.cds.vcf.gz|") or die $!;
		$i=0;
		while (<IN>){
			next if (/#/);
			$i++;
		}
		$comp{$pre}{cds}=$i;
		close IN;
	}

	open STAT, ">$outpath/stat.txt";
	print STAT "Name\tTs/Tv\tDust\tCDS\n";
	foreach my $pre (keys %comp){
		print STAT "$pre\t$comp{$pre}{TsTv}\t$comp{$pre}{dust}\t$comp{$pre}{cds}\n";
	}
	close STAT;
	#To be done
}
#---------------------------------------------------------------------------------------------------


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
#--------------------------------------03.GeneticRelationships--------------------------------------

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
#################################
#			   #
#   	Step3_admixture    		#
#			   #
#################################
sub ADMIXTURE{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/03.GeneticRelationships/Admixture/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/03.GeneticRelationships/Admixture/";
	my $plink_data = $cfg{step1}{variant_filtering}{plink_data};
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/step3_admixture.sh";
	
	print SH "cd $outpath\n";
	print SH "cp -fr $plink_data* ./\n";
	open CL, ">$shpath/step3_admixture.cmd.list";
	for (my $k=0; $k<9;$k++){
		print CL "admixture --cv $plink_data.bed $k 1>$shpath/$k.admixture.o 2>$shpath/$k.admixture.e\n";
	}
	close CL;

	# generate the sample list file (for the sample order in the graph)
	`cp -f $Bin/lib/admixture.R $shpath/admixture.R`;
	if ((defined $cfg{step3}{admixture}{samplelist})&&(-e $cfg{step3}{admixture}{samplelist})){
		`cp -f $cfg{step3}{admixture}{samplelist} $outpath/sample.list`;
	}else{
		my @a = keys %{$cfg{population}};
		open SL, ">$outpath/sample.list";
		foreach (@a){
			print SL "$_\n";
		}
		close SL;
	}

	# generate the arguments for Rscript
	my @p;
	push @p, $outpath; 
	push @p, "$outpath/sample.list";
	push @p, "K";
	my $p = join (' ',@p);

	print SH "parallel -j $cfg{args}{threads} < $shpath/step3_admixture.cmd.list\n";
	print SH "cat $shpath/*admixture.o|grep \"CV error\" >$outpath/CV.error.txt\n";
	print SH "Rscript --vanilla $shpath/admixture.R $p\n";
	close SH;

	`sh $shpath/step3_admixture.sh 1>$shpath/step3_admixture.sh.o 2>$shpath/step3_admixture.sh.e` unless ($skipsh ==1);
}

sub PCA{
	#plink --bfile ../SNP_prunning/prunedData --pca --chr-set 95 no-xy no-mt no-y
}
#---------------------------------------04.InterPopulation------------------------------------------
#################################
#			   					#
#      step4 Divergence         #
#			   					#
#################################
sub DIVERGENCE{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/04.InterPopulation/Divergence/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/04.InterPopulation/Divergence/";

	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{vcf};

	$cfg{step4}{divergence}{windowsize} = 50000 unless (defined $cfg{step4}{divergence}{windowsize});
	$cfg{step4}{divergence}{scaffold_number_limit} = 6 unless (defined $cfg{step4}{divergence}{scaffold_number_limit});
	$cfg{step4}{divergence}{scaffold_length_cutoff} = 1000000 unless (defined $cfg{step4}{divergence}{scaffold_length_cutoff});
	
	my $window_size = $cfg{step4}{divergence}{windowsize};
	my $scaffold_number_limit = $cfg{step4}{divergence}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{divergence}{scaffold_length_cutoff};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}
	
	my $pop_list = "$outpath/pop.list";
	open OT, ">$pop_list";
	foreach my $sample (keys %{$cfg{population}}){
		print OT "$sample\t$cfg{population}{$sample}{'presumed population'}\n";
	}
	close OT;
	###
	`cp -f $Bin/lib/Divergence.R $shpath/Divergence.R`;

	open SH, ">$shpath/Divergence.sh";
	print SH "cd $outpath\n";
	print SH "source activate py27\n";
	print SH "#python /root/genomics_general/VCF_processing/parseVCF.py -i $gzvcf --skipIndels --minQual 0 --gtf flag=DP min=0 | gzip > $outpath/output.geno.gz\n";
	### generate population pairwise information
	open CL, ">$shpath/divergence.cmd.list";
	for (my $i = 0; $i < @{$cfg{step4}{divergence}{poppair}}; $i++){
		my $pop1 = $cfg{step4}{divergence}{poppair}->[$i]->{pop1};
		my $pop2 = $cfg{step4}{divergence}{poppair}->[$i]->{pop2};
		my $poppair = "$pop1-vs-$pop2";
		
		print SH "python /root/genomics_general/popgenWindows.py -w $window_size -m 5 -g output.geno.gz -f haplo -o $outpath/$poppair.csv.gz -T 5 -p $pop1 -p $pop2 --popsFile $pop_list\n" if ($cfg{args}{ploidy} == 1);
		print SH "python /root/genomics_general/popgenWindows.py -w $window_size -m 5 -g output.geno.gz -f phased -o $outpath/$poppair.csv.gz -T 5 -p $pop1 -p $pop2 --popsFile $pop_list\n" if ($cfg{args}{ploidy} >= 2);
		my $i=1;
		foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
			if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				my @p;
				push @p, $outpath; 
				push @p, "$id.$poppair.csv";
				push @p, "$id.$poppair.divergence.png";
				push @p, $window_size;
				push @p, $id;
				my $p = join (' ',@p);
			
				open IDSH, ">$shpath/$id.$poppair.divergence.sh\n";
				print IDSH "zcat $poppair.csv.gz\|grep -w \"$id\" >$outpath/$id.$poppair.csv\n";
				print IDSH "Rscript --vanilla $shpath/Divergence.R $p\n";
				close IDSH;

				print CL "sh $shpath/$id.$poppair.divergence.sh 1>$shpath/$id.$poppair.divergence.sh.o 2>$shpath/$id.$poppair.divergence.sh.e\n";
				$i++;
			}
		}
	}
	close CL;

	print SH "source deactivate\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/divergence.cmd.list\n";
	close SH;

	`sh $shpath/Divergence.sh 1>$shpath/Divergence.sh.o 2>$shpath/Divergence.sh.e` unless ($skipsh ==1);
}	
#---------------------------------------05.IntraPopulation------------------------------------------

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
#################################
#			   #
#   	step4_LD			   	#
#
#################################
sub LD{
	my ($yml_file,$skipsh) = @_;
	my $ori_gzvcf;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/LD/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/LD/";

	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{high_confidence_vcf};
	if ($cfg{args}{ploidy} == 1 ){
		$ori_gzvcf = $gzvcf;
		$gzvcf = "$outpath/diploid.high_confidence.vcf.gz";
	}

	$cfg{step4}{ld}{scaffold_number_limit} = 5000 unless (defined $cfg{step4}{ld}{scaffold_number_limit});
	$cfg{step4}{ld}{scaffold_length_cutoff} = 5000 unless (defined $cfg{step4}{ld}{scaffold_length_cutoff});
	my $scaffold_number_limit = $cfg{step4}{ld}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{ld}{scaffold_length_cutoff};

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

	`cp -f $Bin/lib/LD.R $shpath/LD.R`;
	open SH, ">$shpath/LD.sh";
	if ($cfg{args}{ploidy} == 1 ){
		print SH <<EOF;
zcat $ori_gzvcf|java -jar /root/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes(variant.getGenotypes().stream().map(G->!G.isCalled()?GenotypeBuilder.createMissing(G.getSampleName(),2):G).map(G->G.isCalled() && G.getPloidy()==1?new GenotypeBuilder(G).alleles(Arrays.asList(G.getAllele(0),G.getAllele(0))).make():G).collect(Collectors.toList())).attribute("AC",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 && !G.getAllele(0).isReference()?2:G.getAlleles().size()).sum()).attribute("AN",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 ?2:G.getAlleles().size()).sum()).make();'|bgzip -c >$gzvcf
EOF
	}
	open CL, ">$shpath/LD.cmd.list";
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		open OT, ">$outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		#write a script to generate a vcf file for each population
		open CLSH, ">$shpath/$pop_name.LD.sh";
		print CLSH "cd $outpath\n";
		print CLSH "vcftools --gzvcf $gzvcf --keep $outpath/$pop_name.list --recode --stdout |bgzip -c >$outpath/$pop_name.SNP.vcf.gz\n";
		print CLSH "vcftools --gzvcf $pop_name.SNP.vcf.gz --singletons --stdout >$outpath/$pop_name.out.singletons\n";
		print CLSH "vcftools --gzvcf $outpath/$pop_name.SNP.vcf.gz --exclude-positions $outpath/$pop_name.out.singletons --recode --stdout |bgzip -c  >$outpath/$pop_name.SNP.noSingle.vcf.gz\n";
		
		my $i=1;
		open PCL, ">$shpath/$pop_name.LD.cmd.list";

#		open RIN, ">$outpath/$pop_name.R.input";
		open GRIN, ">$outpath/$pop_name.GLD.R.input";
		foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
			if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){

				open IDSH, ">$shpath/$pop_name.$id.LD.sh";
				print IDSH "cd $outpath\n";
				print IDSH "vcftools --gzvcf $pop_name.SNP.noSingle.vcf.gz --chr $id --recode --stdout |bgzip -c  >$pop_name.$id.SNP.noSingle.vcf.gz\n";
				print IDSH "rm -rf $pop_name.$id.SNP.noSingle.beagl*\n";
				print IDSH "beagle gt=$pop_name.$id.SNP.noSingle.vcf.gz out=$pop_name.$id.SNP.noSingle.beagle\n";
#				print IDSH "vcftools --gzvcf $pop_name.$id.SNP.noSingle.beagle.vcf.gz --hap-r2 --ld-window-bp 1000000 --stdout |grep -v nan > $pop_name.$id.LD_window_1M.list\n";
				print IDSH "vcftools --gzvcf $pop_name.$id.SNP.noSingle.beagle.vcf.gz --geno-r2 --ld-window-bp 1000000 --stdout |grep -v nan | perl -ne '\@a=split /\\t+/;if (/CHR/){print;}elsif(((\$a[2]-\$a[1])>5000)&&(\$i!= 100)){\$i++;}elsif((\$a[2]-\$a[1])<=5000){print;}elsif(((\$a[2]-\$a[1])>5000 )&& (\$i ==100)){print;\$i=0;}'> $pop_name.$id.GLD_window_1M.list\n";
				print IDSH "window_LD.pl $pop_name.$id.GLD_window_1M.list $cfg{step4}{slidingwindow}{windowsize} ",$genome->{len}{$id}," > $pop_name.$id.GLD_window.stats\n";
#				print IDSH "cat $pop_name.$id.LD_window_1M.list",'|sed 1,1d | awk -F " " \'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($3-$2),$5}\''," >$pop_name.$id.LD_window_1M.summary\n";
				print IDSH "cat $pop_name.$id.GLD_window_1M.list",'|sed 1,1d | awk -F " " \'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($3-$2),$5}\''," >$pop_name.$id.GLD_window_1M.summary\n";
				close IDSH;

				print PCL "sh $shpath/$pop_name.$id.LD.sh 1>$shpath/$pop_name.$id.LD.sh.o 2>$shpath/$pop_name.$id.LD.sh.e\n";
#				print RIN "$pop_name.$id.LD_window_1M.summary\n";
				print GRIN "$pop_name.$id.GLD_window_1M.summary\n";
				$i++;
			}
		}
#		close RIN;
		close GRIN;

		close PCL;
		print CLSH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.LD.cmd.list\n"; # parallel run PCL 
		
#		my @p;
#		push @p, $outpath; 
#		push @p, "$pop_name.R.input";
#		push @p, "$pop_name.LD.png";
#		my $p = join (' ',@p);
#		print CLSH "Rscript --vanilla $shpath/LD.R $p\n";

		my @gp;
		push @gp, $outpath; 
		push @gp, "$pop_name.GLD.R.input";
		push @gp, "$pop_name.GLD.png";
		my $gp = join (' ',@gp);
		print CLSH "Rscript --vanilla $shpath/LD.R $gp\n";

		close CLSH;

		print CL "sh $shpath/$pop_name.LD.sh 1>$shpath/$pop_name.LD.sh.o 2>$shpath/$pop_name.LD.sh.e\n";		
	}
	close CL;
	print SH "parallel -j $cfg{args}{threads} < $shpath/LD.cmd.list\n";
	close SH;

	`sh $shpath/LD.sh 1>$shpath/LD.sh.o 2>$shpath/LD.sh.e` unless ($skipsh ==1);

}

#################################
#			   #
#   	step4_SlidingWindow 	#
#			   #
#################################
sub SLIDINGWINDOW{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/Slidingwindow/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/Slidingwindow/";

	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{vcf};

	$cfg{step4}{slidingwindow}{windowsize} = 5000 unless (defined $cfg{step4}{slidingwindow}{windowsize});
	$cfg{step4}{slidingwindow}{scaffold_number_limit} = 10000 unless (defined $cfg{step4}{slidingwindow}{scaffold_number_limit});
	$cfg{step4}{slidingwindow}{scaffold_length_cutoff} = 5000 unless (defined $cfg{step4}{slidingwindow}{scaffold_length_cutoff});
	my $window_size = $cfg{step4}{slidingwindow}{windowsize};
	my $scaffold_number_limit = $cfg{step4}{slidingwindow}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{slidingwindow}{scaffold_length_cutoff};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}

	### circle all the populations
	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\n";
	}

	`cp -f $Bin/lib/Diversity.R $shpath/Diversity.R`;
	open SH, ">$shpath/Slidingwindow.sh";
	print SH "cd $outpath\n";
	
	my $localgzvcf = $gzvcf;
	# prepare for a psudo-diploid vcf file for haploid 
#	my $localgzvcf = "$outpath/di_PASS.SNP.DP.vcf.gz";
#	if ($cfg{args}{ploidy} == 1){
#	print SH <<EOF;
#zcat $gzvcf|java -jar /root/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes(variant.getGenotypes().stream().map(G->!G.isCalled()?GenotypeBuilder.createMissing(G.getSampleName(),2):G).map(G->G.isCalled() && G.getPloidy()==1?new GenotypeBuilder(G).alleles(Arrays.asList(G.getAllele(0),G.getAllele(0))).make():G).collect(Collectors.toList())).attribute("AC",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 && !G.getAllele(0).isReference()?2:G.getAlleles().size()).sum()).attribute("AN",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 ?2:G.getAlleles().size()).sum()).make();'|bgzip -c >$localgzvcf
#EOF
#}else{
#	print SH "cp -f $gzvcf $localgzvcf\n";
#}	
	
	#generate a scaffold size list
	open SL, ">$outpath/chr.size.list";
	foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
		print SL "$id\t$genome->{len}{$id}\n";
	}close SL;

	`grep CDS $cfg{step4}{slidingwindow}{gff} >$outpath/CDS.gff`;
	#generate Bed files
	my $i = 1;
	foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
		if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
			open BED, ">$outpath/$id.bed";
			for (my $j = 0;($j + $window_size) <= $genome->{len}{$id};$j+=$window_size ){
				my $start = $j; my $end = $j + $window_size; 
				print BED "$id\t$start\t$end\n";
			}
			close BED;
			`bedtools intersect -a $outpath/$id.bed -b $outpath/CDS.gff -wo >$outpath/$id.intersected.bed`;
		}
	}

	open CL1, ">$shpath/Slidingwindow.cmd1.list";
	open CL2, ">$shpath/Slidingwindow.cmd2.list";
	open CL3, ">$shpath/Slidingwindow.cmd3.list";
	open CL4, ">$shpath/Slidingwindow.cmd4.list";
	open CL5, ">$shpath/Slidingwindow.cmd5.list";

	#loop for each available population
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		open OT, ">$outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;
		# generate a vcf file for each population
		print CL1 "vcftools --gzvcf $localgzvcf --keep $outpath/$pop_name.list --recode --stdout |bgzip -c >$outpath/$pop_name.SNP.vcf.gz\n";
		print CL1 "vcftools --gzvcf $outpath/snpEff.nonsyn.vcf.gz --keep $outpath/$pop_name.list --recode --stdout |bgzip -c >$outpath/$pop_name.snpEff.nonsyn.vcf.gz\n";
		print CL1 "vcftools --gzvcf $outpath/snpEff.syn.vcf.gz --keep $outpath/$pop_name.list --recode --stdout |bgzip -c >$outpath/$pop_name.snpEff.syn.vcf.gz\n";
		
		# calculate syn and nonsyn diversity values in each sliding window
		print CL2 "vcftools --gzvcf $outpath/$pop_name.snpEff.nonsyn.vcf.gz --window-pi $window_size --stdout >$outpath/$pop_name.snpEff.nonsyn.pi.list\n";
		print CL2 "vcftools --gzvcf $outpath/$pop_name.snpEff.syn.vcf.gz --window-pi $window_size --stdout >$outpath/$pop_name.snpEff.syn.pi.list\n";
		print CL3 "PiSynNonSyn.pl $outpath/$pop_name.snpEff.nonsyn.pi.list $outpath/$pop_name.snpEff.nonsyn.pi.list SynSitesBin.list >$outpath/$pop_name.PiSynNonSyn.list\n";
		
		# plot diversity
		my $i = 1;
		foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
			if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				
				my $len=$genome->{len}{$id};
				my $bigwin=11*$window_size;
				# get the diploid slidingwindow stats
				print CL4 "python /root/diploSHIC/diploSHIC.py fvecVcf diploid $outpath/$pop_name.SNP.vcf.gz $id $len $outpath/$id.$pop_name.SNP.fvec --winSize $bigwin --statFileName $outpath/$id.$pop_name.stat.result 1>$shpath/$id.$pop_name.stat.result.o 2>$shpath/$id.$pop_name.stat.result.e\n" if ($cfg{args}{ploidy} == 2);
				# get the haploid slidingwindow stats
				print CL4 "python /root/diploSHIC/diploSHIC.py fvecVcf haploid $outpath/$pop_name.SNP.vcf.gz $id $len $outpath/$id.$pop_name.SNP.fvec --winSize $bigwin --ancFileName $cfg{ref}{db}{$cfg{ref}{choose}}{path} --statFileName $outpath/$id.$pop_name.stat.result 1>$shpath/$id.$pop_name.stat.result.o 2>$shpath/$id.$pop_name.stat.result.e\n" if ($cfg{args}{ploidy} == 1);

				# calculate GC content, generate the table
				open GC, ">$outpath/$id.$pop_name.GCstat.result";
				for (my $j = 0;($j + $window_size) <= $genome->{len}{$id};$j+=$window_size ){
					my $GCcount = uc(substr($genome->{seq}{$id},$j,$window_size)) =~ tr/GC//;
					my $start = $j + 1; my $end = $j + $window_size; my $GCcontent = $GCcount/$window_size;
					print GC "$id\t$start\t$end\t$GCcontent\n";
				}
				close GC;

				my @p;
				push @p, $outpath; 
				push @p, "$id.$pop_name.Pi.list";
				push @p, "$id.$pop_name.SNPdensity.list";
				push @p, "$id.$pop_name.TajimaD.list";
				push @p, "$id.$pop_name.diversity.png";
				push @p, "$id.$pop_name.TajimaD.png";
				push @p, $pop{$pop_name}{count};
				push @p, $window_size;
				push @p, $id;
				my $p = join (' ',@p);
				
				open IDSH, ">$shpath/$id.$pop_name.Slidingwindow.sh";
				print IDSH "grep -w \"$id\" $outpath/$pop_name.Pi.list >$outpath/$id.$pop_name.Pi.list\n";
				print IDSH "grep -w \"$id\" $outpath/$pop_name.SNPdensity.list >$outpath/$id.$pop_name.SNPdensity.list\n";
				print IDSH "grep -w \"$id\" $outpath/$pop_name.TajimaD.list >$outpath/$id.$pop_name.TajimaD.list\n";
				print IDSH "Rscript --vanilla $shpath/Diversity.R $p\n";
				close IDSH;
				print CL5 "sh $shpath/$id.$pop_name.Slidingwindow.sh 1>$shpath/$id.$pop_name.Slidingwindow.sh.o 2>$shpath/$id.$pop_name.Slidingwindow.sh.e\n";

				$i++;
			}
		}
	}
	close CL1;
	close CL2;
	close CL3;
	close CL4;
	close CL5;
	###
	### annotate  vcf file
	print SH "snpEff $cfg{step4}{slidingwindow}{snpeff_species} $localgzvcf |bgzip -c >$outpath/snpEff.vcf.gz\n";
	### extract nonsyn snvs
	print SH "zcat $outpath/snpEff.vcf.gz |SnpSift filter \"(ANN[*].BIOTYPE has 'protein_coding') & ((ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT = 'start_lost')| (ANN[*].EFFECT = 'stop_gained')| (ANN[*].EFFECT = 'stop_lost'))\" |perl -ne 'if(/#/){print;}elsif(/ANN/){s/ANN\\S+/./g;print;}else{print;}'|bgzip -c >$outpath/snpEff.nonsyn.vcf.gz\n";
	### extract syn snvs
	print SH "zcat $outpath/snpEff.vcf.gz |SnpSift filter \"(ANN[*].BIOTYPE has 'protein_coding') & ((ANN[*].EFFECT = 'synonymous_variant') | (ANN[*].EFFECT = 'stop_retained_variant') | (ANN[*].EFFECT = 'start_retained'))\" |perl -ne 'if(/#/){print;}elsif(/ANN/){s/ANN\\S+/./g;print;}else{print;}'|bgzip -c >$outpath/snpEff.syn.vcf.gz\n";
	print SH "Number_syn_nonsyn_sites.pl $cfg{ref}{db}{$cfg{ref}{choose}}{path} $cfg{step4}{slidingwindow}{gff} test >$outpath/nonsyn.sites.list\n";
	print SH "SynSitesBin.pl chr.size.list $outpath/nonsyn.sites.list |grep -v NULL|sort -k1,1 -k2,2n >$outpath/SynSitesBin.list\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/Slidingwindow.cmd1.list\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/Slidingwindow.cmd2.list\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/Slidingwindow.cmd3.list\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/Slidingwindow.cmd4.list\n";
	print SH "#parallel -j $cfg{args}{threads} < $shpath/Slidingwindow.cmd5.list\n";
	close SH;

	`sh $shpath/Slidingwindow.sh 1>$shpath/Slidingwindow.sh.o 2>$shpath/Slidingwindow.sh.e` unless ($skipsh ==1);
	
	foreach my $pop_name (keys %pop){
		open ALLSTAT, ">$outpath/$pop_name.allstat.txt";
		next unless ($pop{$pop_name}{count} > 6);
		my $i = 1;
		foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
			if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				#$outpath/$id.$pop_name.stat.result
				#$outpath/../LD/$pop_name.$id.GLD_window.stats
				#$outpath/$id.$pop_name.GCstat.result
				open IDALLSTAT, ">$outpath/$id.$pop_name.allstat.txt";
				my %stats;
				open STAT, "$outpath/$id.$pop_name.stat.result";
				while(<STAT>){
					next if (/chrom/);
					chomp;
					my @a = split /\s+/;
					$stats{$a[1]}{diversity}="$a[3]\t$a[4]\t$a[5]";
					#print "$a[0]\t$a[1]\t$a[2]\t$stats{$a[1]}{diversity}\n";
				}
				close STAT;

				open STAT, "$outpath/$id.intersected.bed";
				while(<STAT>){
					next if (/chrom/);
					chomp;
					my @a = split /\s+/;
					my $start = $a[1]+1;
					my $p = $a[-1]/$window_size;
					if (exists $stats{$start}{p_coding}){
						$stats{$start}{p_coding}+= $p; 
					}else{
						$stats{$start}{p_coding}=$p;
					}
					#print "$a[0]\t$a[1]\t$a[2]\t$stats{$a[1]}{diversity}\n";
				}
				close STAT;

				open STAT, "$outpath/../LD/$pop_name.$id.GLD_window.stats";
				while(<STAT>){
					chomp;
					my @a = split /\t+/;
					$stats{$a[1]}{r2}="$a[3]";
				}
				close STAT;

				open STAT, "$outpath/$id.$pop_name.GCstat.result";
				while(<STAT>){
					chomp;
					my @a = split /\s+/;
					$stats{$a[1]}{r2} = "nan" unless (exists $stats{$a[1]}{r2});
					$stats{$a[1]}{p_coding} = 0 unless (exists $stats{$a[1]}{p_coding});
					if (exists $stats{$a[1]}{diversity}){
						print IDALLSTAT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$stats{$a[1]}{p_coding}\t$stats{$a[1]}{diversity}\t$stats{$a[1]}{r2}\n";
						print ALLSTAT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$stats{$a[1]}{p_coding}\t$stats{$a[1]}{diversity}\t$stats{$a[1]}{r2}\n";
					}
				}
				close STAT;
				close IDALLSTAT;
			}
			$i++;
		}
		close ALLSTAT;

		my %stats;
		open STAT, "$outpath/$pop_name.PiSynNonSyn.list";
		while(<STAT>){
			next if (/CHROM/);
			chomp;
			my @a = split /\t+/;
			$stats{$a[0]}{$a[1]}{pnps}="$a[3]\t$a[4]\t$a[5]";
		}
		close STAT;

		open ALLSTAT, ">$outpath/final.$pop_name.allstat.txt";
		print ALLSTAT "chrom\tstart\tend\tgc_content\tcoding_percentage\tpi\ttheta\ttajimsD\tr2\tps\tpn\tpnps\n";
		open STAT, "$outpath/$pop_name.allstat.txt";
		while(<STAT>){
			chomp;
			my @a = split /\s+/;
			$stats{$a[0]}{$a[1]}{pnps} = "nan\tnan\tnan" unless (exists $stats{$a[0]}{$a[1]}{pnps});
			print ALLSTAT "$_\t$stats{$a[0]}{$a[1]}{pnps}\n";
		}
		close STAT;
		close ALLSTAT;
		`Rscript --vanilla $Bin/lib/Correlation.R $outpath $outpath/final.$pop_name.allstat.txt`;
	}
	
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $yml_file );
}

#---------------------------------------05.DemographicHistory------------------------------------------
sub SFS{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/SFS/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/SFS/";
	my $slidingwindowpath = "$cfg{args}{outdir}/05.IntraPopulation/Slidingwindow/";

	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{vcf};
	my $localgzvcf = $gzvcf;
	#$localgzvcf = $gzvcf if ($cfg{args}{ploidy} == 2);

	$cfg{step4}{slidingwindow}{windowsize} = 5000 unless (defined $cfg{step4}{slidingwindow}{windowsize});
	$cfg{step4}{slidingwindow}{scaffold_number_limit} = 6 unless (defined $cfg{step4}{slidingwindow}{scaffold_number_limit});
	$cfg{step4}{slidingwindow}{scaffold_length_cutoff} = 1000000 unless (defined $cfg{step4}{slidingwindow}{scaffold_length_cutoff});
	$cfg{step4}{discoal}{hard_simulation_times} = 2000 unless (defined $cfg{step4}{discoal}{hard_simulation_times});
	$cfg{step4}{discoal}{soft_simulation_times} = 500 unless (defined $cfg{step4}{discoal}{soft_simulation_times});
	$cfg{step4}{discoal}{neut_simulation_times} = 2000 unless (defined $cfg{step4}{discoal}{neut_simulation_times});
	my $window_size = $cfg{step4}{slidingwindow}{windowsize};
	my $scaffold_number_limit = $cfg{step4}{slidingwindow}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{slidingwindow}{scaffold_length_cutoff};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/SFS.sh";

	open CL1, ">$shpath/SFS.cmd1.list";
	#loop for each available population
	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\tpop1\n";
	}
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		open OT, ">$outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		my $pop_size = 2* $pop{$pop_name}{count};
		if ($cfg{args}{ploidy} == 1 ){ $pop_size = $pop{$pop_name}{count};}
		open IDSH, ">$shpath/$pop_name.SFS.sh";
		print IDSH "mkdir -p $outpath/$pop_name && cd $outpath/$pop_name\n";
		### project SFS from the file
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.SNP.vcf.gz -p $outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $outpath/$pop_name/total_sfs && mv $outpath/$pop_name/output $outpath/$pop_name/total_sfs \n";
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.snpEff.nonsyn.vcf.gz -p $outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $outpath/$pop_name/nonsyn_sfs && mv $outpath/$pop_name/output $outpath/$pop_name/nonsyn_sfs \n";
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.snpEff.syn.vcf.gz -p $outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $outpath/$pop_name/syn_sfs  && mv $outpath/$pop_name/output $outpath/$pop_name/syn_sfs \n";

		print IDSH "cp -f $Bin/lib/1PopBot20Mb.est $outpath/$pop_name/ && cp -f $Bin/lib/1PopBot20Mb.tpl $outpath/$pop_name/ && cp -f $outpath/$pop_name/total_sfs/fastsimcoal2/pop1_MAFpop0.obs $outpath/$pop_name/1PopBot20Mb_MAFpop0.obs\n";
		print IDSH "sed -i `s/18/$pop_size/g` $outpath/$pop_name/1PopBot20Mb_MAFpop0.obs\n";
		for (my $i = 0; $i <10; $i++){
			print IDSH "fsc26 -t 1PopBot20Mb.tpl -e 1PopBot20Mb.est -n 10000 -d -M -L 40 -q -0 -c 40 -B 40 --foldedSFS\ncp -r 1PopBot20Mb 1PopBot20Mb.b$i\n";
		}
		print IDSH "cat $outpath/$pop_name/1PopBot20Mb.b*/1PopBot20Mb.bestlhoods | grep -v NCUR |sort -k1,1n > $outpath/$pop_name/EstimatedNe.list\n";
		close IDSH;
		print CL1 "sh $shpath/$pop_name.SFS.sh 1>$shpath/$pop_name.SFS.sh.o 2>$shpath/$pop_name.SFS.sh.e\n";
	}
	close CL1;
	print SH "parallel -j $cfg{args}{threads} < $shpath/SFS.cmd1.list\n";
	close SH;

	`sh $shpath/SFS.sh 1>$shpath/SFS.sh.o 2>$shpath/SFS.sh.e` unless ($skipsh ==1);
	
	foreach my $pop_name (keys %pop){
		open SH, ">$shpath/$pop_name.Simulation.sh";
		next unless ($pop{$pop_name}{count} > 6);
		open NELT, "$outpath/$pop_name/EstimatedNe.list";
		my @nelt = <NELT>;
		close NELT;

		my $length = @nelt;
		my @a = split /\s+/, $nelt[$length/2];
		print "$a[0]\n";

		my $m_rate=0.000000025;
		my $Bigwindowsize = 11*$cfg{step4}{slidingwindow}{windowsize};
		my $Ne = $a[0];
		my $Nan = $a[1];
		my $Nbot = $a[2];
		my $Tbot = $a[3];
		my $Tan = $a[6];

		my $Pt_min = 4 * $Ne * $m_rate * $Bigwindowsize/3.16227766;
		my $Pt_max = 4 * $Ne * $m_rate * $Bigwindowsize * 3.16227766;
		my $Pre_min = 4 * $Ne * $m_rate * $Bigwindowsize/3.16227766;
		my $Pre_max = 4 * $Ne * $m_rate * $Bigwindowsize * 3.16227766;
		my $Pa_min = 4 * $Ne * 0.01/3.16227766;
		my $Pa_max = 4 * $Ne * 1 * 3.16227766;
		my $T = 10000/$Ne;

		my $step = 1/11;
		my $init = 0.5/11;
		my $pop_size = 2*$pop{$pop_name}{count};

		open SIMUCL, ">$shpath/$pop_name.Simulation.cmd.list";
		open FVECCL, ">$shpath/$pop_name.SimFvec.cmd.list";
		for (my $i=0; $i<11;$i++){
			my $x = $init + $step*$i;

			# simulate hard sweeps
			print SIMUCL "/root/discoal/discoal $pop_size $cfg{step4}{discoal}{hard_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
			print SIMUCL " -x $x >$outpath/$pop_name/hard_$i.msOut\n";
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim diploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 2);
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim haploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 1);
			# simulate soft sweeps
			print SIMUCL "/root/discoal/discoal $pop_size $cfg{step4}{discoal}{soft_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 -Pf 0 0.1 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
			print SIMUCL " -x $x >$outpath/$pop_name/soft_$i.msOut\n";
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim diploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 2);
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim haploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 1);
		}
		print SIMUCL "/root/discoal/discoal $pop_size $cfg{step4}{discoal}{neut_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max;
		print SIMUCL " -en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
		print SIMUCL " >$outpath/$pop_name/neut.msOut\n";
		print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim diploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 2);
		print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim haploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 1);
		close SIMUCL;

		print SH "cd $outpath/$pop_name/\n";
		print SH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.Simulation.cmd.list\n";
		print SH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.SimFvec.cmd.list\n";
		print SH "mkdir -p $outpath/$pop_name/rawFVFiles && mv $outpath/$pop_name/*.fvec rawFVFiles/\n";
		print SH "mkdir -p $outpath/$pop_name/trainingSets\n";
		print SH "python /root/diploSHIC/diploSHIC.py makeTrainingSets rawFVFiles/neut.fvec rawFVFiles/soft rawFVFiles/hard 5 0,1,2,3,4,6,7,8,9,10 trainingSets/\n";
		print SH "mkdir -p $outpath/$pop_name/updatedSets\n";
		print SH "less trainingSets/linkedHard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/linkedHard.fvec\n";
		print SH "less trainingSets/hard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/hard.fvec\n";
		print SH "less trainingSets/linkedSoft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/linkedSoft.fvec\n";
		print SH "less trainingSets/neut.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/neut.fvec\n";
		print SH "less trainingSets/soft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/soft.fvec\n";
		print SH "python /root/diploSHIC/diploSHIC.py train updatedSets/ updatedSets/ bfsModel\n";
		print SH "mkdir -p $outpath/$pop_name/observedFVFiles && cp $slidingwindowpath/*.$pop_name.SNP.fvec $outpath/$pop_name/observedFVFiles/\n";

		my $i = 1;
		open PREDCL, ">$shpath/$pop_name.predict.cmd.list";
		foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
			if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				my $len=$genome->{len}{$id};
				print PREDCL "python /root/diploSHIC/diploSHIC.py predict bfsModel.json bfsModel.weights.hdf5 $outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.fvec $outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.preds\n";
				$i++;
			}
		}
		close PREDCL;
		print SH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.predict.cmd.list\n";
		close SH;
		`sh $shpath/$pop_name.Simulation.sh 1>$shpath/$pop_name.Simulation.sh.o 2>$shpath/$pop_name.Simulation.sh.e` unless ($skipsh ==1);
	}
}

#---------------------------------------06.Selection------------------------------------------------
sub MKTEST{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.Selection/MKtest/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.Selection/MKtest/";

	#my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{high_confidence_vcf};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	print "$cfg{step6}{mktest}{mclgroup}\n";
	open IN, $cfg{step6}{mktest}{mclgroup};
	my @label = ($cfg{step6}{mktest}{label1}, $cfg{step6}{mktest}{label2});
	print "$cfg{step6}{mktest}{label1}\t$cfg{step6}{mktest}{label2}\n";
	open OT, ">$outpath/sco_pair.list";
	while (<IN>){
		my $flag = 0;
		foreach my $label (@label){
			my $c = () = ($_ =~ /$label\|/g);
			$flag ++ if ($c == 1);
		}
		if ($flag == @label){
        	foreach my $label (@label){
				if ($_ =~ /$label\|(\S+)/){print OT "$1\t";}
			}
			print OT "\n";
		}
	}
	close OT;
	close IN;

	my $genepair = "$outpath/sco_pair.list";
	my $cds1 = $cfg{step6}{mktest}{cds1};
	my $cds2 = $cfg{step6}{mktest}{cds2};

	my $cdsio1  = Bio::SeqIO->new(-format => 'fasta', -file   => $cds1);
	my %cds1;
	while (my $cdsobj = $cdsio1->next_seq){ $cds1{$cdsobj->display_id}=$cdsobj->seq;}

	my $cdsio2  = Bio::SeqIO->new(-format => 'fasta', -file   => $cds2);
	my %cds2;
	while (my $cdsobj = $cdsio2->next_seq){ $cds2{$cdsobj->display_id}=$cdsobj->seq;}

	open CL, ">$shpath/aln.cmd.list";
	if ( !-d "$outpath/prank_alignment" ) {make_path "$outpath/prank_alignment" or die "Failed to create path: $outpath/prank_alignment";}
	if ( !-d "$shpath/prank_alignment" ) {make_path "$shpath/prank_alignment" or die "Failed to create path: $shpath/prank_alignment";}
	open IN, $genepair;
	while (<IN>){
    	my @a = split /\t/;
    	my $out = $a[0].".multi.cds";
    	
    	open OUT, ">$outpath/prank_alignment/$out";
    	print OUT ">$a[0]\n$cds1{$a[0]}\n";
    	print OUT ">$a[1]\n$cds2{$a[1]}\n";
    	close OUT;
   		
   		print CL "prank -d=$outpath/prank_alignment/$out -o=$outpath/prank_alignment/$out -translate -F 1>$shpath/prank_alignment/$out.prank.sh.o 2>$shpath/prank_alignment/$out.prank.sh.e &\n";    	
	}
	close IN;
	close CL;

	open SH, ">$shpath/MKtest.sh";
	print SH "cd $outpath/prank_alignment\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/aln.cmd.list\n";
	print SH "sed -i 's/NNN/---/g' $outpath/prank_alignment/*.nuc.fas\n";
	print SH "evaluate_alignment_quality.pl >$outpath/filtered.SCO.list\n";
	print SH "DetectSynandNonSyn.pl black $outpath/prank_alignment/multi.cds.best.nuc.fas\n";
	print SH "select_hash.pl $outpath/filtered.SCO.list black_syn_nonsyn.txt >$outpath/Divergence.txt\n";
	print SH "snpEff $cfg{step6}{mktest}{snpeff_species} $gzvcf |bgzip -c >$outpath/snpEff.vcf.gz\n";
	print SH "step6.polymorphism.pl $outpath/prank_alignment/snpEff_genes.txt >$outpath/polymorphism.txt\n";
	print SH "MK-test.pl $outpath/polymorphism.txt $outpath/Divergence.txt >$outpath/MK-test.SCO.result\n";
	close SH;

	`sh $shpath/MKtest.sh 1>$shpath/MKtest.sh.o 2>$shpath/MKtest.sh.e` unless ($skipsh ==1);

	open IN, "$outpath/polymorphism.txt";
	my $piN = 0; 
	my $piS = 0;
	while (<IN>){
		my @a = split /\s+/;
		#print "$a[1]\t$a[2]\n";
		$piN = $piN + $a[1];
		$piS = $piS + $a[2];
	}
	close IN;
	print "piN/PiS: ", $piN/$piS, "\n";
}

#################################
#			  #
#   	Step6_FST Outliers	   	#
#			  #
#################################
sub FSTOUTLIERS{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/06.Selection/Fst_outliers/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/06.Selection/Fst_outliers/";

	my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{vcf};


}

#----------------------------------- other subroutines ---------------------------------------------
sub LOADREF{
	my ($reference) = @_;
	my %genome;
	open REF, "$reference";
	my $id;
	while (<REF>){
		chomp;
		if (/\>(\S+)/){
			$id = $1;
			$genome{seq}{$id}="";
		}
		elsif(/(\w+)/){
			$genome{seq}{$id}.=$1;
		}
	}
	close REF;
	$genome{sumlen}=0;
	foreach my $id(keys %{$genome{seq}}){
		$genome{len}{$id}=length($genome{seq}{$id});
		$genome{sumlen} += $genome{len}{$id};
	}
	return \%genome;
}

1;