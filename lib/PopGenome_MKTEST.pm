package PopGenome_MKTEST;
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

sub MKTEST{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my $outpath = "$cfg{args}{outdir}/05.Selection/MKtest/";
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/05.Selection/MKtest/";

	#my $genome = LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{step1}{variant_filtering}{high_confidence_vcf};

	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

	print "$cfg{step6}{mktest}{mclgroup}\n";
	open IN, $cfg{step6}{mktest}{mclgroup};
	my @label = ($cfg{step6}{mktest}{label1}, $cfg{step6}{mktest}{label2});
	print "$cfg{step6}{mktest}{label1}\t$cfg{step6}{mktest}{label2}\n";
	open OT, ">$outpath/sco_pair.list";
	while (<IN>){
		my $flag = 0;
		foreach my $label (@label){
			my $c = () = ($_ =~ /$label\|/g);
			$flag ++ if ($c == 1);
		}
		if ($flag == @label){
        	foreach my $label (@label){
				if ($_ =~ /$label\|(\S+)/){print OT "$1\t";}
			}
			print OT "\n";
		}
	}
	close OT;
	close IN;

	my $genepair = "$outpath/sco_pair.list";
	my $cds1 = $cfg{step6}{mktest}{cds1};
	my $cds2 = $cfg{step6}{mktest}{cds2};

	my $cdsio1  = Bio::SeqIO->new(-format => 'fasta', -file   => $cds1);
	my %cds1;
	while (my $cdsobj = $cdsio1->next_seq){ $cds1{$cdsobj->display_id}=$cdsobj->seq;}

	my $cdsio2  = Bio::SeqIO->new(-format => 'fasta', -file   => $cds2);
	my %cds2;
	while (my $cdsobj = $cdsio2->next_seq){ $cds2{$cdsobj->display_id}=$cdsobj->seq;}

	open CL, ">$shpath/aln.cmd.list";
	if ( !-d "$outpath/prank_alignment" ) {make_path "$outpath/prank_alignment" or die "Failed to create path: $outpath/prank_alignment";}
	if ( !-d "$shpath/prank_alignment" ) {make_path "$shpath/prank_alignment" or die "Failed to create path: $shpath/prank_alignment";}
	open IN, $genepair;
	while (<IN>){
    	my @a = split /\t/;
    	my $out = $a[0].".multi.cds";
    	
    	open OUT, ">$outpath/prank_alignment/$out";
    	print OUT ">$a[0]\n$cds1{$a[0]}\n";
    	print OUT ">$a[1]\n$cds2{$a[1]}\n";
    	close OUT;
   		
   		print CL "prank -d=$outpath/prank_alignment/$out -o=$outpath/prank_alignment/$out -translate -F 1>$shpath/prank_alignment/$out.prank.sh.o 2>$shpath/prank_alignment/$out.prank.sh.e &\n";    	
	}
	close IN;
	close CL;

	open SH, ">$shpath/MKtest.sh";
	print SH "cd $outpath/prank_alignment\n";
	print SH "parallel -j $cfg{args}{threads} < $shpath/aln.cmd.list\n";
	print SH "sed -i 's/NNN/---/g' $outpath/prank_alignment/*.nuc.fas\n";
	print SH "evaluate_alignment_quality.pl >$outpath/filtered.SCO.list\n";
	print SH "DetectSynandNonSyn.pl black $outpath/prank_alignment/multi.cds.best.nuc.fas\n";
	print SH "select_hash.pl $outpath/filtered.SCO.list black_syn_nonsyn.txt >$outpath/Divergence.txt\n";
	print SH "snpEff $cfg{step6}{mktest}{snpeff_species} $gzvcf |bgzip -c >$outpath/snpEff.vcf.gz\n";
	print SH "step6.polymorphism.pl $outpath/prank_alignment/snpEff_genes.txt >$outpath/polymorphism.txt\n";
	print SH "MK-test.pl $outpath/polymorphism.txt $outpath/Divergence.txt >$outpath/MK-test.SCO.result\n";
	close SH;

	`sh $shpath/MKtest.sh 1>$shpath/MKtest.sh.o 2>$shpath/MKtest.sh.e` unless ($skipsh ==1);

	open IN, "$outpath/polymorphism.txt";
	my $piN = 0; 
	my $piS = 0;
	while (<IN>){
		my @a = split /\s+/;
		#print "$a[1]\t$a[2]\n";
		$piN = $piN + $a[1];
		$piS = $piS + $a[2];
	}
	close IN;
	print "piN/PiS: ", $piN/$piS, "\n";
}

1;