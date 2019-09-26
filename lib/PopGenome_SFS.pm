package PopGenome_SFS;
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

#---------------------------------------05.DemographicHistory------------------------------------------
sub SFS{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.IntraPopulation/SFS/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.IntraPopulation/SFS/";
	my $slidingwindowpath = "$cfg{args}{outdir}/05.IntraPopulation/Slidingwindow/";

	my $genome = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{vcf};
	my $localgzvcf = $gzvcf;
	#$localgzvcf = $gzvcf if ($cfg{args}{ploidy} == 2);

	$cfg{step4}{slidingwindow}{windowsize} = 5000 unless (defined $cfg{step4}{slidingwindow}{windowsize});
	$cfg{step4}{slidingwindow}{scaffold_number_limit} = 6 unless (defined $cfg{step4}{slidingwindow}{scaffold_number_limit});
	$cfg{step4}{slidingwindow}{scaffold_length_cutoff} = 1000000 unless (defined $cfg{step4}{slidingwindow}{scaffold_length_cutoff});
	$cfg{step4}{discoal}{hard_simulation_times} = 2000 unless (defined $cfg{step4}{discoal}{hard_simulation_times});
	$cfg{step4}{discoal}{soft_simulation_times} = 500 unless (defined $cfg{step4}{discoal}{soft_simulation_times});
	$cfg{step4}{discoal}{neut_simulation_times} = 2000 unless (defined $cfg{step4}{discoal}{neut_simulation_times});
	my $window_size = $cfg{step4}{slidingwindow}{windowsize};
	my $scaffold_number_limit = $cfg{step4}{slidingwindow}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{slidingwindow}{scaffold_length_cutoff};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	open SH, ">$shpath/SFS.sh";

	open CL1, ">$shpath/SFS.cmd1.list";
	#loop for each available population
	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\tpop1\n";
	}
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		open OT, ">$outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		my $pop_size = 2* $pop{$pop_name}{count};
		if ($cfg{args}{ploidy} == 1 ){ $pop_size = $pop{$pop_name}{count};}
		open IDSH, ">$shpath/$pop_name.SFS.sh";
		print IDSH "mkdir -p $outpath/$pop_name && cd $outpath/$pop_name\n";
		### project SFS from the file
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.SNP.vcf.gz -p $outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $outpath/$pop_name/total_sfs && mv $outpath/$pop_name/output $outpath/$pop_name/total_sfs \n";
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.snpEff.nonsyn.vcf.gz -p $outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $outpath/$pop_name/nonsyn_sfs && mv $outpath/$pop_name/output $outpath/$pop_name/nonsyn_sfs \n";
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.snpEff.syn.vcf.gz -p $outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $outpath/$pop_name/syn_sfs  && mv $outpath/$pop_name/output $outpath/$pop_name/syn_sfs \n";

		print IDSH "cp -f $Bin/lib/1PopBot20Mb.est $outpath/$pop_name/ && cp -f $Bin/lib/1PopBot20Mb.tpl $outpath/$pop_name/ && cp -f $outpath/$pop_name/total_sfs/fastsimcoal2/pop1_MAFpop0.obs $outpath/$pop_name/1PopBot20Mb_MAFpop0.obs\n";
		print IDSH "sed -i `s/18/$pop_size/g` $outpath/$pop_name/1PopBot20Mb_MAFpop0.obs\n";
		for (my $i = 0; $i <10; $i++){
			print IDSH "fsc26 -t 1PopBot20Mb.tpl -e 1PopBot20Mb.est -n 10000 -d -M -L 40 -q -0 -c 40 -B 40 --foldedSFS\ncp -r 1PopBot20Mb 1PopBot20Mb.b$i\n";
		}
		print IDSH "cat $outpath/$pop_name/1PopBot20Mb.b*/1PopBot20Mb.bestlhoods | grep -v NCUR |sort -k1,1n > $outpath/$pop_name/EstimatedNe.list\n";
		close IDSH;
		print CL1 "sh $shpath/$pop_name.SFS.sh 1>$shpath/$pop_name.SFS.sh.o 2>$shpath/$pop_name.SFS.sh.e\n";
	}
	close CL1;
	print SH "parallel -j $cfg{args}{threads} < $shpath/SFS.cmd1.list\n";
	close SH;

	`sh $shpath/SFS.sh 1>$shpath/SFS.sh.o 2>$shpath/SFS.sh.e` unless ($skipsh ==1);
	
	foreach my $pop_name (keys %pop){
		open SH, ">$shpath/$pop_name.Simulation.sh";
		next unless ($pop{$pop_name}{count} > 6);
		open NELT, "$outpath/$pop_name/EstimatedNe.list";
		my @nelt = <NELT>;
		close NELT;

		my $length = @nelt;
		my @a = split /\s+/, $nelt[$length/2];
		print "$a[0]\n";

		my $m_rate=0.000000025;
		my $Bigwindowsize = 11*$cfg{step4}{slidingwindow}{windowsize};
		my $Ne = $a[0];
		my $Nan = $a[1];
		my $Nbot = $a[2];
		my $Tbot = $a[3];
		my $Tan = $a[6];

		my $Pt_min = 4 * $Ne * $m_rate * $Bigwindowsize/3.16227766;
		my $Pt_max = 4 * $Ne * $m_rate * $Bigwindowsize * 3.16227766;
		my $Pre_min = 4 * $Ne * $m_rate * $Bigwindowsize/3.16227766;
		my $Pre_max = 4 * $Ne * $m_rate * $Bigwindowsize * 3.16227766;
		my $Pa_min = 4 * $Ne * 0.01/3.16227766;
		my $Pa_max = 4 * $Ne * 1 * 3.16227766;
		my $T = 10000/$Ne;

		my $step = 1/11;
		my $init = 0.5/11;
		my $pop_size = 2*$pop{$pop_name}{count};

		open SIMUCL, ">$shpath/$pop_name.Simulation.cmd.list";
		open FVECCL, ">$shpath/$pop_name.SimFvec.cmd.list";
		for (my $i=0; $i<11;$i++){
			my $x = $init + $step*$i;

			# simulate hard sweeps
			print SIMUCL "/root/discoal/discoal $pop_size $cfg{step4}{discoal}{hard_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
			print SIMUCL " -x $x >$outpath/$pop_name/hard_$i.msOut\n";
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim diploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 2);
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim haploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 1);
			# simulate soft sweeps
			print SIMUCL "/root/discoal/discoal $pop_size $cfg{step4}{discoal}{soft_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 -Pf 0 0.1 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
			print SIMUCL " -x $x >$outpath/$pop_name/soft_$i.msOut\n";
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim diploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 2);
			print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim haploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 1);
		}
		print SIMUCL "/root/discoal/discoal $pop_size $cfg{step4}{discoal}{neut_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max;
		print SIMUCL " -en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
		print SIMUCL " >$outpath/$pop_name/neut.msOut\n";
		print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim diploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 2);
		print FVECCL "python /root/diploSHIC/diploSHIC.py fvecSim haploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($cfg{args}{ploidy} == 1);
		close SIMUCL;

		print SH "cd $outpath/$pop_name/\n";
		print SH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.Simulation.cmd.list\n";
		print SH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.SimFvec.cmd.list\n";
		print SH "mkdir -p $outpath/$pop_name/rawFVFiles && mv $outpath/$pop_name/*.fvec rawFVFiles/\n";
		print SH "mkdir -p $outpath/$pop_name/trainingSets\n";
		print SH "python /root/diploSHIC/diploSHIC.py makeTrainingSets rawFVFiles/neut.fvec rawFVFiles/soft rawFVFiles/hard 5 0,1,2,3,4,6,7,8,9,10 trainingSets/\n";
		print SH "mkdir -p $outpath/$pop_name/updatedSets\n";
		print SH "less trainingSets/linkedHard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/linkedHard.fvec\n";
		print SH "less trainingSets/hard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/hard.fvec\n";
		print SH "less trainingSets/linkedSoft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/linkedSoft.fvec\n";
		print SH "less trainingSets/neut.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/neut.fvec\n";
		print SH "less trainingSets/soft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$outpath/$pop_name/updatedSets/soft.fvec\n";
		print SH "python /root/diploSHIC/diploSHIC.py train updatedSets/ updatedSets/ bfsModel\n";
		print SH "mkdir -p $outpath/$pop_name/observedFVFiles && cp $slidingwindowpath/*.$pop_name.SNP.fvec $outpath/$pop_name/observedFVFiles/\n";

		my $i = 1;
		open PREDCL, ">$shpath/$pop_name.predict.cmd.list";
		foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
			if (($genome->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				my $len=$genome->{len}{$id};
				print PREDCL "python /root/diploSHIC/diploSHIC.py predict bfsModel.json bfsModel.weights.hdf5 $outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.fvec $outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.preds\n";
				$i++;
			}
		}
		close PREDCL;
		print SH "parallel -j $cfg{args}{threads} < $shpath/$pop_name.predict.cmd.list\n";
		close SH;
		`sh $shpath/$pop_name.Simulation.sh 1>$shpath/$pop_name.Simulation.sh.o 2>$shpath/$pop_name.Simulation.sh.e` unless ($skipsh ==1);
	}
}

1;