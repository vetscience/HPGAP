package PopGenome_Mapping;
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

1;