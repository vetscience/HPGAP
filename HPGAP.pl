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
use PopGenome;
use YAML::Tiny;

my ($config, $step, $run, $skipsh, $verbose);

GetOptions (
	"config=s" => \$config,
	"step=s" => \$step, #string
	"run=s" => \$run,
	"skipsh" => \$skipsh,
	"verbose"  => \$verbose)   # flag
or die ("Error in command line arguments\n");

#------------------------------------- default parameters ---------------------------------------
#$filter ||= "-l 15 -m 3 -p ATGC,10 -n 1 -z";
$run ||= '';
$skipsh ||= 0;
#--------------------------------------- get step options ----------------------------------------
#print "$verbose\n";
unless ( defined ($config) ) {
	&usage;
	exit;
}	

#my $yml_file = shift;
my $yaml = YAML::Tiny->read( $config );
my %cfg = %{$yaml->[0]};

`cp -f $config $cfg{args}{outdir}/allcfg.yml` unless ( -e "$cfg{args}{outdir}/allcfg.yml");
my $allcfg = "$cfg{args}{outdir}/allcfg.yml";

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
my $main = "$cfg{args}{outdir}/HPGAP.main.sh";
open MH, ">$main"; print MH "#!/bin/sh\ncd $cfg{args}{outdir}\n";
my $udocker_cmd="udocker run ";
	for (my $i=0;$i<@{$cfg{args}{mount}};$i++){
		$udocker_cmd .= "-v $cfg{args}{mount}->[$i]->{hostpath}:$cfg{args}{mount}->[$i]->{dockerpath} ";
	}
	$udocker_cmd .= "--env=\"$cfg{args}{env}\" $cfg{args}{container} /bin/bash -c ";
	
#$genome=PopGenome::LOADREF($reference);

#Quality control and variant calling
#if (exists $step{0}{A});
#if (exists $step{1}{A});
#if (exists $step{1}{B});
#if (exists $step{1}{C});
#if (exists $step{1}{D});
#if (exists $step{1}{E});
#if (exists $step{1}{F});
#print "$udocker_cmd 'HPGAP.pl --run step1_variant_filtering --config $allcfg'\n" if (exists $step{1}{variant_filtering});

#01.indexing
PopGenome::INDEXING($allcfg,$skipsh) if ($run eq 'step0_indexing');

PopGenome::DATA_FILTERING($allcfg,$skipsh) if ($run eq 'step1_read_filtering');

PopGenome::MAPPING($allcfg,$skipsh) if ($run eq 'step1_read_mapping');

PopGenome::CALIBRATION($allcfg,$skipsh) if ($run eq 'step1_recalibration');

PopGenome::VARIANT_CALLING($allcfg,$skipsh) if ($run eq 'step1_variant_calling');

PopGenome::COMBINE_CALLING($allcfg,$skipsh) if ($run eq 'step1_combine_calling');

print "$udocker_cmd 'HPGAP.pl --run step1_variant_filtering --config $allcfg'\n" if (exists $step{1}{variant_filtering});
PopGenome::VARIANT_FILTERING($allcfg,$skipsh) if ($run eq 'step1_variant_filtering');


print "$udocker_cmd 'HPGAP.pl --run step1_comparison --config $allcfg'\n" if (exists $step{1}{Comparison});
PopGenome::VARIANT_COMPARISON($allcfg,$skipsh) if ($run eq 'step1_comparison');
#############################################
#											#
#	Require statistics of variants	
#											#
#############################################
#udocker run -v /home/darcy/tmp:/tmp -v /mnt/NAS/Clonorchis_sinensis/:/mnt/NAS/Clonorchis_sinensis/ -v /home/darcy/PopGen_WorkFlow/:/home/darcy/PopGen_WorkFlow/ --env="PATH=/root/admixture_linux-1.3.0/:/root/gatk:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin/:/usr/sbin:/usr/bin:/sbin:/bin:/home/darcy/PopGen_WorkFlow/Pipeline/:/home/darcy/PopGen_WorkFlow/Pipeline/lib:/home/darcy/PopGen_WorkFlow/Pipeline/Tools" HPGAP_c1 /bin/bash -c 'HPGAP.pl --config /home/darcy/PopGen_WorkFlow/Cb_Validation/allcfg.yml --run step5_ld'

#02.SampleFiltering
print "$udocker_cmd 'HPGAP.pl --run step2_relatedness --config $allcfg'\n" if (exists $step{2}{Relatedness});
PopGenome::RELATEDNESS($allcfg,$skipsh) if ($run eq 'step2_relatedness');
#03.GeneticRelationships
#print "$udocker_cmd 'HPGAP.pl --run step3_snp_prunning --config $allcfg'\n" if (exists $step{3}{SNP_prunning});
#PopGenome::SNP_prunning($allcfg,$skipsh) if ($run eq 'step3_snp_prunning');

print "$udocker_cmd 'HPGAP.pl --run step3_phylogeny --config $allcfg'\n" if (exists $step{3}{Phylogeny});
PopGenome::PHYLOGENY($allcfg,$skipsh) if ($run eq 'step3_phylogeny');

print "$udocker_cmd 'HPGAP.pl --run step3_admixture --config $allcfg'\n" if (exists $step{3}{Admixture});
PopGenome::ADMIXTURE($allcfg,$skipsh) if ($run eq 'step3_admixture');

#04.InterPopulation
print "$udocker_cmd 'HPGAP.pl --run step4_divergence --config $allcfg'\n" if (exists $step{4}{Divergence});
PopGenome::DIVERGENCE($allcfg,$skipsh) if ($run eq 'step4_divergence');

#05.IntraPopulation
print "$udocker_cmd 'HPGAP.pl --run step5_diversity --config $allcfg'\n" if (exists $step{5}{slidingwindow});
PopGenome::SLIDINGWINDOW($allcfg,$skipsh) if ($run eq 'step5_slidingwindow');
#LD: Only use diploid individuals (vcftools).
PopGenome::SFS($allcfg,$skipsh) if ($run eq 'step5_sfs');


print "$udocker_cmd 'HPGAP.pl --run step5_homozygosity --config $allcfg'\n" if (exists $step{5}{Homozygosity});
PopGenome::HOMOZYGOSITY($allcfg,$skipsh) if ($run eq 'step5_homozygosity');
#if (exists $step{5}{Homozygosity});
print "$udocker_cmd 'HPGAP.pl --run step5_ld --config $allcfg'\n" if (exists $step{5}{LD});
PopGenome::LD($allcfg,$skipsh) if ($run eq 'step5_ld');
#LD: Only use diploid individuals (vcftools).
print "$udocker_cmd 'HPGAP.pl --run step5_roh --config $allcfg'\n" if (exists $step{5}{ROH});
PopGenome::ROH($allcfg,$skipsh) if ($run eq 'step5_roh');
#06.DemographicHistory
#if (exists $step{6}{PSMC});

#06.Selection
print "$udocker_cmd 'HPGAP.pl --run step6_mktest --config $allcfg'\n" if (exists $step{6}{MKtest});
PopGenome::MKTEST($allcfg,$skipsh) if ($run eq 'step6_mktest');

print "$udocker_cmd 'HPGAP.pl --run step6_fstoutliers --config $allcfg'\n" if (exists $step{6}{FSToutliers});
PopGenome::FSTOUTLIERS($allcfg,$skipsh) if ($run eq 'step6_fstoutliers');


close MH;
#----------------------------------- usage sub progamm ------------------------------------------
sub usage{
	print STDERR <<USAGE;

Description
	For Popolation genetic analysis in helminths.

Version 
	13 Dec 2018: Version v1.0.0

Author
	Daxi Wang
	email: darcywdx\@gmail.com
	please contact me if you find any bug.
							   	
Usage
	HPGAP.pl --samplelist sample.list --reference reference.fa [-options]
	
	--variant_filtering
	--step <String>			specified steps , separated by commas (e.g., "0:A;1:A,B,C,E,F,G").
		0.	Indexing
		1a. Read filtering
		1b.	Read mapping
		1c.	Variant calling
		2.	Sample filtering
		3. 	Inferring genetic relationships
		4. 	Defining inter-population structures
		5.	Intra-population statistics
		6.	Demographic inference
		7.	Identification of loci under natural selection
		8.	Generating reports

Note 
	1. No symbolic link to the files outside the mounted volume, which means all the data files themselves should be located within the mounted volume.
	2. For each pair of fastq data, the second colomn (library or flowcell code) should always be unique.
	3. All the paths need to be absolute path

Example 
	To be done
USAGE
}
