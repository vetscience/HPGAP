package PopGenome_LD;
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

1;