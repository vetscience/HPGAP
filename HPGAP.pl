#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use FindBin '$Bin';
use lib "$Bin/lib";
use lib '/root/miniconda3/lib/site_perl/5.26.2';
#use PopGenome;
use YAML::Tiny;

#######
require PopGenome_Indexing;
require PopGenome_Data_Filtering;
require PopGenome_Mapping;
require PopGenome_Read_Report;
require PopGenome_Calibration;
require PopGenome_Variant_Calling;
require PopGenome_Combine_Calling;
require PopGenome_Variant_Filtering;
require PopGenome_Relatedness;
require PopGenome_Phylogeny;
require PopGenome_Admixture;
require PopGenome_Homozygosity;
require PopGenome_ROH;
require PopGenome_LD;
require PopGenome_Slidingwindow;
require PopGenome_SFS;
require PopGenome_MKTEST;
#######
my ($config, $step, $run, $skipsh, $help, $run_flag);

GetOptions (
	"config=s" => \$config,
	"step=s" => \$step, #string
	"run=s" => \$run,
	"skipsh" => \$skipsh,
	"help"  => \$help)   # flag
or die ("Error in command line arguments\n");
if (defined $run){ $run_flag = 1;}
#------------------------------------- default parameters ---------------------------------------
#$filter ||= "-l 15 -m 3 -p ATGC,10 -n 1 -z";
$run ||= '';
$skipsh ||= 0;
$step ||= '0:indexing;1:read_filtering,read_mapping,read_report,recalibration,variant_calling,combine_calling,variant_filtering;3:phylogeny,admixture;4:homozygosity,roh,ld,slidingwindow,sfs';
#--------------------------------------- get step options ----------------------------------------
#print "$verbose\n";
if ( defined ($help) ) {
	&usage;
	exit;
}	
unless ( defined ($config) ) {
	&usage;
	exit;
}	
#my $yml_file = shift;
my $yaml = YAML::Tiny->read( $config );
my %cfg = %{$yaml->[0]};
unless (exists $cfg{args}{threads}){$cfg{args}{threads}=40}

my $allcfg = "$cfg{args}{outdir}/allcfg.yml";

unless ( -e $allcfg){
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $allcfg );
}

my %step;
if(defined $step){
	my @step = split /;/, $step;
	foreach my $sopt(@step){
		chomp $sopt;
		my $sid = $1 if ($sopt =~ /^(\d+)\:/);
		$sopt =~ s/^\d+\://g;
		my @a = split /,/, $sopt;
		foreach (@a){
			my $opt = $1 if(/(\w+)/);
			$step{$sid}{$opt}=1 if (defined $opt);
		}
	}
}

#-------------------------------------- Main steps  --------------------------------------------------
unless (defined $run_flag){
	my $main = "$cfg{args}{outdir}/HPGAP.main.sh";
	open MH, ">$main"; print MH "#!/bin/sh\ncd $cfg{args}{outdir}\n";
	my $udocker_cmd="udocker run ";
		for (my $i=0;$i<@{$cfg{args}{mount}};$i++){
			if (exists $cfg{args}{mount}->[$i]->{host_tmp}){
				$udocker_cmd .= "-v $cfg{args}{mount}->[$i]->{host_tmp}:/tmp ";
			}elsif (exists $cfg{args}{mount}->[$i]->{host_path}){
				$udocker_cmd .= "-v $cfg{args}{mount}->[$i]->{host_path}:$cfg{args}{mount}->[$i]->{host_path} ";
			}
		}
		$udocker_cmd .= "--env=\"$cfg{args}{env}\" $cfg{args}{container} /bin/bash -c ";
		
	#01.indexing
	print MH "time $udocker_cmd 'HPGAP.pl --run step0_indexing --config $allcfg' && echo 'indexing done'\n" if (exists $step{0}{indexing});

	print MH "time $udocker_cmd 'HPGAP.pl --run step1_read_filtering --config $allcfg'  && echo 'read_filtering done'\n" if (exists $step{1}{read_filtering});

	print MH "time $udocker_cmd 'HPGAP.pl --run step1_read_mapping --config $allcfg' && echo 'read_mapping done'\n" if (exists $step{1}{read_mapping});
	
	print MH "time $udocker_cmd 'HPGAP.pl --run step1_read_report --config $allcfg' && echo 'read_report done'\n" if (exists $step{1}{read_report});

	print MH "time $udocker_cmd 'HPGAP.pl --run step1_recalibration --config $allcfg' && echo 'recalibration done'\n" if (exists $step{1}{recalibration});

	print MH "time $udocker_cmd 'HPGAP.pl --run step1_variant_calling --config $allcfg' && echo 'variant_calling done'\n" if (exists $step{1}{variant_calling});

	print MH "time $udocker_cmd 'HPGAP.pl --run step1_combine_calling --config $allcfg' && echo 'combine_calling done'\n" if (exists $step{1}{combine_calling});

	print MH "time $udocker_cmd 'HPGAP.pl --run step1_variant_filtering --config $allcfg' && echo 'variant_filtering done'\n" if (exists $step{1}{variant_filtering});


	#print "$udocker_cmd 'HPGAP.pl --run step1_comparison --config $allcfg'\n" if (exists $step{1}{Comparison});
	#PopGenome::VARIANT_COMPARISON($allcfg,$skipsh) if ($run eq 'step1_comparison');
	#############################################
	#											#
	#	Require statistics of variants	
	#											#
	#############################################

	#02.SampleFiltering
	print MH "time $udocker_cmd 'HPGAP.pl --run step2_relatedness --config $allcfg' && echo 'relatedness done'\n" if (exists $step{2}{relatedness});

	print MH "time $udocker_cmd 'HPGAP.pl --run step3_phylogeny --config $allcfg' && echo 'phylogeny done'\n" if (exists $step{3}{phylogeny});

	print MH "time $udocker_cmd 'HPGAP.pl --run step3_admixture --config $allcfg' && echo 'admixture done'\n" if (exists $step{3}{admixture});

	print MH "time $udocker_cmd 'HPGAP.pl --run step4_homozygosity --config $allcfg' && echo 'homozygosity done'\n" if (exists $step{5}{homozygosity});

	print MH "time $udocker_cmd 'HPGAP.pl --run step4_roh --config $allcfg' && echo 'roh done'\n" if (exists $step{4}{roh});

	print MH "time $udocker_cmd 'HPGAP.pl --run step4_ld --config $allcfg' && echo 'ld done'\n" if (exists $step{4}{ld});

	#05.IntraPopulation
	print MH "time $udocker_cmd 'HPGAP.pl --run step4_slidingwindow --config $allcfg' && echo 'slidingwindow done'\n" if (exists $step{4}{slidingwindow});

	print MH "time $udocker_cmd 'HPGAP.pl --run step4_sfs --config $allcfg' && echo 'sfs done'\n" if (exists $step{4}{sfs});

	#06.Selection
	print MH "time $udocker_cmd 'HPGAP.pl --run step6_mktest --config $allcfg' && echo 'MKtest done'\n" if (exists $step{6}{MKtest});
	close MH;
}

	PopGenome_Indexing::INDEXING($allcfg,$skipsh) if ($run eq 'step0_indexing');

	PopGenome_Data_Filtering::DATA_FILTERING($allcfg,$skipsh) if ($run eq 'step1_read_filtering');

	PopGenome_Mapping::MAPPING($allcfg,$skipsh) if ($run eq 'step1_read_mapping');

	PopGenome_Read_Report::READ_REPORT($allcfg,$skipsh) if ($run eq 'step1_read_report');

	PopGenome_Calibration::CALIBRATION($allcfg,$skipsh) if ($run eq 'step1_recalibration');

	PopGenome_Variant_Calling::VARIANT_CALLING($allcfg,$skipsh) if ($run eq 'step1_variant_calling');

	PopGenome_Combine_Calling::COMBINE_CALLING($allcfg,$skipsh) if ($run eq 'step1_combine_calling');

	PopGenome_Variant_Filtering::VARIANT_FILTERING($allcfg,$skipsh) if ($run eq 'step1_variant_filtering');


	#############################################
	#											#
	#	Require statistics of variants	
	#											#
	#############################################

	#02.SampleFiltering
	PopGenome_Relatedness::RELATEDNESS($allcfg,$skipsh) if ($run eq 'step2_relatedness');

	PopGenome_Phylogeny::PHYLOGENY($allcfg,$skipsh) if ($run eq 'step3_phylogeny');
	
	PopGenome_Admixture::ADMIXTURE($allcfg,$skipsh) if ($run eq 'step3_admixture');

	PopGenome_Homozygosity::HOMOZYGOSITY($allcfg,$skipsh) if ($run eq 'step4_homozygosity');

	PopGenome_ROH::ROH($allcfg,$skipsh) if ($run eq 'step4_roh');

	PopGenome_LD::LD($allcfg,$skipsh) if ($run eq 'step4_ld');

	#05.IntraPopulation
	PopGenome_Slidingwindow::SLIDINGWINDOW($allcfg,$skipsh) if ($run eq 'step4_slidingwindow');

	PopGenome_SFS::SFS($allcfg,$skipsh) if ($run eq 'step4_sfs');

	#06.Selection
	PopGenome_MKTEST::MKTEST($allcfg,$skipsh) if ($run eq 'step6_mktest');

#----------------------------------- usage sub progamm ------------------------------------------

sub usage{
	print STDERR <<USAGE;

Description
	For Popolation genetic analysis in helminths.

Version 
	06 Feb 2019: Version v1.0.0

Author
	Daxi Wang
	email: darcywdx\@gmail.com
	please contact me if you find any bug.
							   	
Usage
	HPGAP.pl --config <path to the .yml config file> --run <String> [-options]
	
	--run <String> use this option choose one of steps below (the option should not be used with --step at the same time)
		step0_indexing
		step1_read_filtering
		step1_read_mapping
		step1_read_report
		step1_recalibration
		step1_variant_calling
		step1_combine_calling
		step1_variant_filtering
		step2_relatedness
		step3_phylogeny
		step3_admixture
		step4_homozygosity
		step4_roh
		step4_ld
		step4_slidingwindow
		step4_sfs

	--config path to the .yml config file
	
	--step <String>	specified steps separated by semicolon(;). The names of analyses in each step are separated by comma (,);
		(e.g. "0:indexing;1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;3:phylogeny,admixture;4:homozygosity,roh,ld,slidingwindow,sfs").

		All the avaliable analyses in each step: 
			0:indexing;
			1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;
			3:phylogeny,admixture;
			4:homozygosity,roh,ld,slidingwindow,sfs
			6:mktest

	--skipsh use this to skip running bash inside the pipeline
	
	--help

Note 
	1. No symbolic link to the files outside the mounted volume, which means all the data files themselves should be located within the mounted volume.
	2. For each pair of fastq data, the second colomn (library or flowcell code) should always be unique.
	3. All the paths need to be absolute path

Example 
	To be done
USAGE
}
