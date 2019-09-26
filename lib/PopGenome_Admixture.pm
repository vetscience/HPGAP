package PopGenome_Admixture;
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

1;