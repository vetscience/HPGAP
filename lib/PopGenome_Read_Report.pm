package PopGenome_Read_Report;
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
#    Step 1c Report summary	    #
#			   #
#################################
sub READ_REPORT{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};

	my %rs;
	my %rs_ref;
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
			`cat $outpath/$sample/bam.stats.txt|grep ^COV | cut -f 2- >$sample_report_outpath/COV.stat.txt`;
			`cat $outpath/$sample/bam.stats.txt|grep ^IS | cut -f 2- >$sample_report_outpath/IS.stat.txt`;
			`cat $outpath/$sample/bam.stats.txt|grep ^SN | cut -f 2-|sed "s/ //g;" >$sample_report_outpath/SN.stat.txt`;
		}		

		open SH, ">$shpath/read_report.sh";
		print SH "#!/bin/sh\ncd $report_outpath\n";
		print SH "Rscript --vanilla $Bin/lib/ReadSummary.R $report_outpath $report_sample_outpath $temp_ref.Summary.xls\n";
		`sh $shpath/read_report.sh 1>$shpath/read_report.sh.o 2>$shpath/read_report.sh.e`;

		open RS, "$report_outpath/$temp_ref.Summary.xls";
		<RS>;
		while(<RS>){
			my @a = split /\t/;
			my $sample_id = $a[0];
			my $mapped_rate = $a[7];
			my $match_rate = 1 - $a[8];
			my $paired_rate = $a[9];
			my $cover_area = $a[12];

			$rs{mapped_rate}{$sample_id}{$temp_ref} = $mapped_rate;
			$rs{match_rate}{$sample_id}{$temp_ref} = $match_rate;
			$rs{paired_rate}{$sample_id}{$temp_ref} = $paired_rate;
			$rs{cover_area}{$sample_id}{$temp_ref} = $cover_area;

			$rs_ref{$temp_ref} = 0

		}
		close RS;

	}
	
	foreach my $meric (keys %rs){
		foreach my $sample_id(keys %{$rs{$meric}}){
			my $temp_top = 0;
			my $temp_top_name = "";
			foreach my $temp_ref(keys %{$rs{$meric}{$sample_id}}){
				if ($rs{$meric}{$sample_id}{$temp_ref} >= $temp_top){
					$temp_top = $rs{$meric}{$sample_id}{$temp_ref};
					$temp_top_name = $temp_ref;
				}	
			}
			$rs_ref{$temp_top_name}++;
		}
	}

	my $temp_top = 0;
	my $temp_top_name = "";
	foreach my $temp_ref(keys %rs_ref){
		if ($rs_ref{$temp_ref}>$temp_top){
			$temp_top = $rs_ref{$temp_ref};
			$temp_top_name = $temp_ref;
		}
	}

	open RS_REF, ">$cfg{args}{outdir}/01.QualityControl/selected_genome.txt";
	print RS_REF "$temp_top_name\t$temp_top\n";
	close RS_REF;

	$cfg{ref}{choose} = $temp_top_name unless (exists $cfg{ref}{choose});
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $yml_file );
}

1;