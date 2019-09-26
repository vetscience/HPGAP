package PopGenome_Slidingwindow;
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
#   	step4_SlidingWindow 	#
#			   #
#################################
sub SLIDINGWINDOW{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/Slidingwindow/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/Slidingwindow/";

	my $genome = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
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

1;