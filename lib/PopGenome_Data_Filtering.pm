package PopGenome_Data_Filtering;
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

1;