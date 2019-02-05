#!/usr/bin/perl

use strict;
use warnings;

my $fl1 = shift;
my $fl2 = shift;
my $outmat=shift;
my $vcf1 = read_vcf($fl1);
my $vcf2 = read_vcf($fl2);

open OM, ">$outmat";
print OM "CHROM\tPOSITION\tSTATUS\tREF\tALT1\tALT2";
foreach my $i (keys $vcf1->{sample}){print OM "\t$vcf1->{sample}{$i}{id}";}
print OM "\n";

foreach my $i (keys $vcf1->{sample}){
	$vcf1->{sample}{$i}{gt1}=0;
	$vcf1->{sample}{$i}{gt2}=0;
	$vcf1->{sample}{$i}{same}=0;
	$vcf1->{sample}{$i}{alt1_m}=0;
	$vcf1->{sample}{$i}{alt2_m}=0;
	$vcf1->{sample}{$i}{unqalt1_m}=0;
	$vcf1->{sample}{$i}{unqalt2_m}=0;
	$vcf1->{sample}{$i}{unqalt1}=0;
	$vcf1->{sample}{$i}{unqalt2}=0;
	$vcf1->{sample}{$i}{alt1}=0;
	$vcf1->{sample}{$i}{alt2}=0;
	$vcf1->{sample}{$i}{alt_same}=0;
	$vcf1->{sample}{$i}{novcf1}=0;
	$vcf1->{sample}{$i}{novcf2}=0;
}
foreach my $chr (keys %{$vcf1->{loc}}){
	foreach my $pos (keys %{$vcf1->{loc}{$chr}}){
		if (exists $vcf2->{loc}{$chr}{$pos}){
			print OM "$chr\t$pos\tB\t$vcf1->{loc}{$chr}{$pos}{ref}\t$vcf1->{loc}{$chr}{$pos}{alt}\t$vcf2->{loc}{$chr}{$pos}{alt}";
			if($vcf1->{loc}{$chr}{$pos}{alt} eq $vcf2->{loc}{$chr}{$pos}{alt}){
				print OM "\tSAME";
			}elsif(($vcf1->{loc}{$chr}{$pos}{alt} =~ /,/) && ($vcf2->{loc}{$chr}{$pos}{alt} =~ /,/)){
				my @al1 = split /,/, $vcf1->{loc}{$chr}{$pos}{alt};
				my @al2 = split /,/, $vcf2->{loc}{$chr}{$pos}{alt};
				
				my %al1;
				$al1{$_}++ for (@al1);
				my %al2;
				$al2{$_}++ for (@al2);
				
				my $flag = 0; 
				foreach my $al1 (@al1){
					$flag ++ unless (exists $al2{$al1});
				}
				foreach my $al2 (@al2){
					$flag ++ unless (exists $al1{$al2});
				}
				print OM "\tSAME" if ($flag == 0);
				print OM "\tDIFF" if ($flag != 0);
			}else{
				print OM "\tDIFF";
			}
			########calculate sample SNP#################
			foreach my $i (keys $vcf1->{sample}){
				
				print OM "\t",$vcf1->{loc}{$chr}{$pos}{$i};
				$vcf1->{sample}{$i}{total}++;
				if(($vcf1->{loc}{$chr}{$pos}{$i} =~ /\w+/) && ($vcf2->{loc}{$chr}{$pos}{$i} =~ /\w+/)){
					#count genotyped number
					$vcf1->{sample}{$i}{gt1}++;
					$vcf1->{sample}{$i}{gt2}++;
					
					my @gt1 = split /\//, $vcf1->{loc}{$chr}{$pos}{$i};
					my @gt2 = split /\//, $vcf2->{loc}{$chr}{$pos}{$i};
					my %gt1;$gt1{$_}++ for (@gt1);
					my %gt2;$gt2{$_}++ for (@gt2);
					
					my $f_gt = 0;
					my $f_alt1 = 0;
					my $f_alt2 = 0;
					foreach my $gt1 (@gt1){
						$f_gt ++ unless (exists $gt2{$gt1});
						$f_alt1 ++ if ($gt1 ne $vcf1->{loc}{$chr}{$pos}{ref});
					}
					foreach my $gt2 (@gt2){
						$f_gt ++ unless (exists $gt1{$gt2});
						$f_alt2 ++ if ($gt2 ne $vcf2->{loc}{$chr}{$pos}{ref});
					}
					#if the called genotypes, including ref and alt, are the same
					$vcf1->{sample}{$i}{same}++ if ($f_gt == 0);
					
					#if the called genotype in vcf 1 is different from the ref (which means alternative allele)
					if ($f_alt1 > 0){
						$vcf1->{sample}{$i}{alt1}++;
						#the number of sites with alternative alleles in vcf1, while vcf2 remains reference aty the same site
						$vcf1->{sample}{$i}{unqalt1}++ if ($f_alt2 == 0);
					}

					#if the called genotype in vcf 2 is different from the ref (which means alternative allele)
					if ($f_alt2 > 0){
						$vcf1->{sample}{$i}{alt2}++;
						$vcf1->{sample}{$i}{unqalt2}++ if ($f_alt1 == 0);
					}
					
					if (($f_alt1 > 0)&&($f_alt2 > 0)&&($f_gt == 0)){
						$vcf1->{sample}{$i}{alt_same}++;
					}
				#if the vcf1 is genotyped, while vcf2 is not 
				}if(($vcf1->{loc}{$chr}{$pos}{$i} =~ /\w+/) && ($vcf2->{loc}{$chr}{$pos}{$i} !~ /\w+/)){
					
					#count genotyped sites
					$vcf1->{sample}{$i}{gt1}++;
					my @gt1 = split /\//, $vcf1->{loc}{$chr}{$pos}{$i};
					my %gt1;$gt1{$_}++ for (@gt1);
					my $f_gt  = 0 ;
					foreach my $gt1 (@gt1){
						$f_gt ++ if ($gt1 ne $vcf1->{loc}{$chr}{$pos}{ref});
					}
					#if the called genotype in vcf 1 is different from the ref
					if ($f_gt > 0){
						$vcf1->{sample}{$i}{alt1}++;
						$vcf1->{sample}{$i}{unqalt1}++;
						#count if missing in the other vcf
						$vcf1->{sample}{$i}{alt1_m}++;
						$vcf1->{sample}{$i}{unqalt1_m}++;
					}
				
				#if the vcf2 is genotyped, while vcf1 is not
				}if(($vcf1->{loc}{$chr}{$pos}{$i} !~ /\w+/) && ($vcf2->{loc}{$chr}{$pos}{$i} =~ /\w+/)){
					$vcf1->{sample}{$i}{gt2}++;
					my @gt1 = split /\//, $vcf2->{loc}{$chr}{$pos}{$i};
					my %gt1;$gt1{$_}++ for (@gt1);
					my $f_gt  = 0 ;
					foreach my $gt1 (@gt1){
						$f_gt ++ if ($gt1 ne $vcf2->{loc}{$chr}{$pos}{ref});
					}
					#if the called genotype in vcf 2 is different from the ref
					if ($f_gt > 0){
						$vcf1->{sample}{$i}{alt2}++;
						$vcf1->{sample}{$i}{unqalt2}++;
						#count if missing in the other vcf
						$vcf1->{sample}{$i}{alt2_m}++;
						$vcf1->{sample}{$i}{unqalt2_m}++;			
					}
				}
			}
			print OM "\n";
			##########END################################
		# 
		}else{
			print OM "$chr\t$pos\t1\t$vcf1->{loc}{$chr}{$pos}{ref}\t$vcf1->{loc}{$chr}{$pos}{alt}\t-";
			print OM "\tNA";
			foreach my $i (keys $vcf1->{sample}){
				print OM "\t",$vcf1->{loc}{$chr}{$pos}{$i};
				$vcf1->{sample}{$i}{total}++;
				if($vcf1->{loc}{$chr}{$pos}{$i} =~ /\w+/){
					$vcf1->{sample}{$i}{gt1}++;
					my @gt1 = split /\//, $vcf1->{loc}{$chr}{$pos}{$i};
					my %gt1;$gt1{$_}++ for (@gt1);
					my $f_gt  = 0 ;
					foreach my $gt1 (@gt1){
						$f_gt ++ if ($gt1 ne $vcf1->{loc}{$chr}{$pos}{ref});
                    }
					if ($f_gt > 0){
						$vcf1->{sample}{$i}{alt1}++;
						$vcf1->{sample}{$i}{unqalt1}++;
						$vcf1->{sample}{$i}{novcf2}++;
					}
				}
			}
			print OM "\n";
		}
	}
}

# loop for the vcf2 sites that do not exist in vcf1
foreach my $chr (keys %{$vcf2->{loc}}){
	foreach my $pos (keys %{$vcf2->{loc}{$chr}}){
		unless (exists $vcf1->{loc}{$chr}{$pos}){
			print OM "$chr\t$pos\t2\t$vcf2->{loc}{$chr}{$pos}{ref}\t-\t$vcf2->{loc}{$chr}{$pos}{alt}";
			print OM "\tNA";
			foreach my $i (keys $vcf2->{sample}){
				print OM "\t",$vcf2->{loc}{$chr}{$pos}{$i};
				$vcf1->{sample}{$i}{total}++;
				if($vcf2->{loc}{$chr}{$pos}{$i} =~ /\w+/){
					$vcf1->{sample}{$i}{gt2}++;
					my @gt1 = split /\//, $vcf2->{loc}{$chr}{$pos}{$i};
					my %gt1;$gt1{$_}++ for (@gt1);
					my $f_gt  = 0 ;
					foreach my $gt1 (@gt1){
						$f_gt ++ if ($gt1 ne $vcf2->{loc}{$chr}{$pos}{ref});
					}
					if ($f_gt > 0){
						$vcf1->{sample}{$i}{alt2}++;
						$vcf1->{sample}{$i}{unqalt2}++;
						$vcf1->{sample}{$i}{novcf1}++;
					}
				}
			}
			print OM "\n";
		}
	}
}

#Output the sample statistics
print "SAMPLE\tTOTAL\tGT1\tGT2\tSAME\t";
print "ALT1(ALT1_MISSING)\tALT2(ALT2_MISSING)\tALT_SAME\t";
print "UNQALT1(UNQALT1_MISSING)\tUNQALT2(UNQALT2_MISSING)\t";
print "NOVCF1\tNOVCF2\n";
foreach my $i (keys $vcf2->{sample}){
	print "$vcf1->{sample}{$i}{id}\t$vcf1->{sample}{$i}{total}\t$vcf1->{sample}{$i}{gt1}\t$vcf1->{sample}{$i}{gt2}\t$vcf1->{sample}{$i}{same}\t";
	print "$vcf1->{sample}{$i}{alt1}($vcf1->{sample}{$i}{alt1_m})\t$vcf1->{sample}{$i}{alt2}($vcf1->{sample}{$i}{alt2_m})\t$vcf1->{sample}{$i}{alt_same}\t";
	print "$vcf1->{sample}{$i}{unqalt1}($vcf1->{sample}{$i}{unqalt1_m})\t$vcf1->{sample}{$i}{unqalt2}($vcf1->{sample}{$i}{unqalt2_m})\t";
	print "$vcf1->{sample}{$i}{novcf1}\t$vcf1->{sample}{$i}{novcf2}\n";
}

sub read_vcf {
	my ($fl1)=@_;
	open (IN, "zcat $fl1|") or die $!;
	my %vcf;
	while (<IN>){
		if (/CHROM/){
			if (/FORMAT\s+(\S.*)/){
				my @a= split /\s+/, $1;
				my $index = 9;
				foreach (@a){
					$vcf{sample}{$index}{id}=$_;		
					$index ++;
				}
			}
		}elsif( $_ !~ /#/){
			#I       151     cbwivar00000001 T       A       13804   PASS    .       GT      .       0       0       .       .       .       .       0 
			my @a= split /\s+/;	
			my $chr = $a[0];
			my $pos = $a[1];
			my $ref = $a[3];
			my $alt = $a[4];
			#from $a[9]
			$vcf{loc}{$chr}{$pos}{ref}=$ref;
			$vcf{loc}{$chr}{$pos}{alt}=$alt;
			for (my $i = 9; $i <  @a; $i++){
				$a[$i]=~s/\s+//g;
				if($a[$i]=~ /^([\.\d]*\/*[\.\d]*)/) {$vcf{loc}{$chr}{$pos}{$i}=$1;}
			}
			my @allele;
			@allele = split /,/, $vcf{loc}{$chr}{$pos}{alt};
			unshift @allele, $vcf{loc}{$chr}{$pos}{ref};
			foreach my $i (keys $vcf{sample}){
				for(my $n=0; $n<4;$n++){
					$vcf{loc}{$chr}{$pos}{$i}=~ s/$n/$allele[$n]/g;
				}
			}
		}
	}
	close IN;
	return \%vcf;
}
