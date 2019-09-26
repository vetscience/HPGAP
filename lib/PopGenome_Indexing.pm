package PopGenome_Indexing;
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

#------------------- Prepare reference index --------------------------------------

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

1;