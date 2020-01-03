#!/usr/bin/perl -w 

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use YAML::Tiny;

my $usage="
    Usage:
    Options:
        -i <input>   input list_file contain 2 col (or 3 col for PE),name and raw_path(must be fastq format)
        -s <seqtype> input sequence type [SE]|PE
        -o <opath>   output file [./result.yml]
        -j <qsub>    qsub the job [Y] or N run on the nodes directly
        -m <mem>     memory needed[2G]
        -t <thread>  thread needed[4]
        -q <queue>   queue [st.q]
        -P <prj>     Project Number [P18Z10200N0119]
        -h|?Help!
    Example:perl $0 -i name_path.list -j Y 
";
my ($input,$seqtype,$outpath,$job,$mem,$thread,$queue,$prj,$help);
GetOptions(
    '-i=s' => \$input,
    '-s=s' => \$seqtype,
    '-o=s' => \$outpath,
    '-j=s' => \$job,
    '-m=s' => \$mem,
    '-t=i' => \$thread,
    '-q=s' => \$queue,
    '-P=s' => \$prj,
    'h|?'  => \$help,
);


if($help or !$input){die "$usage\n";}
$outpath ||= "./result.yml";
$job ||="Y";
$queue ||= "st.q";
$prj ||= "P18Z10200N0119";
$mem ||= "2G";
$thread ||= 4;
$seqtype ||= "SE";

my $yaml = YAML::Tiny->read( $config );
my %cfg = %{$yaml->[0]};
unless (exists $cfg{args}{threads}){$cfg{args}{threads}=40}

open IN,"$input" or die $!; 
while(<Raw>){
    chomp;
    if($seqtype eq "SE"){
        my($name,$path) = (split /\s+/,$_)[0,1];
        $raw{$name} = $path;
    }elsif($seqtype eq "PE"){
        my($name,$path1,$path2) = (split /\s+/,$_)[0,1,2];
  
    }
}
close IN;

my $allcfg = "$cfg{args}{outdir}/allcfg.yml";

unless ( -e $allcfg){
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $allcfg );
}